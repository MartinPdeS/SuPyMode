#!/usr/bin/env python
# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import Any, Optional, Sequence, Tuple
from MPSPlots.styles import mps
from MPSPlots import colormaps
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from SuPyMode.binary.interface_supermode import SUPERMODE
from SuPyMode.utils import interpret_mode_of_interest
from SuPyMode.binary.interface_taper import AlphaProfile


@dataclass
class _ProfileEntry:
    profile: Any
    label: Optional[str] = None
    color: Optional[str] = None
    linestyle: str = "--"
    linewidth: float = 2.0


@dataclass(frozen=True)
class PlotterOptions:
    """
    Configuration controlling how objects are selected and combined for plotting.
    """

    mode_of_interest: list[SUPERMODE] | str = "all"
    combination: list | str = "pairs"
    show_crossings: bool = False


class Plotter:
    """
    Plotter configured by a single plot kind, fed by explicit add methods.

    The contract is simple:
    - Plot kind decides what is rendered.
    - add_profile, add_superset, add_supermode decide what data participates.
    - show renders.

    Important behavior
    ------------------
    If at least one profile is added, and the plot kind supports profile overlay
    (currently adiabatic), the profile curves are overlaid automatically.
    There is no add_profiles flag.

    Supported kinds
    ---------------
    - "adiabatic": plots adiabatic coupling curves from supersets and/or supermodes,
      and overlays profile adiabatic curves if profiles were added.
    - "index": plots effective index vs ITR from supersets and/or supermodes.
    - "beta": plots beta vs ITR from supersets and/or supermodes.
    - "eigenvalue": plots eigenvalue vs ITR from supersets and/or supermodes.
    """

    def __init__(
        self,
        *,
        kind: str,
        use_mps_style: bool = True,
        options: Optional[PlotterOptions] = None,
    ):
        self.kind = self._normalize_kind(kind)
        self.use_mps_style = use_mps_style
        self.options = PlotterOptions() if options is None else options

        self._supersets: list[Any] = []
        self._extra_supermodes: list[SUPERMODE] = []
        self._profiles: list[_ProfileEntry] = []

        self._supermodes_from_supersets: dict[tuple[Any, Any], SUPERMODE] = {}

    @staticmethod
    def _normalize_kind(kind: str) -> str:
        return kind.lower().replace(" ", "").replace("_", "")

    def add_profile(
        self,
        *profiles,
        label: str | None = None,
        color: str | None = None,
        linestyle: str = "--",
        linewidth: float = 2.0,
    ):
        """
        Add one or more taper profiles.

        Styling arguments apply to all profiles passed in this call.
        """
        for profile in profiles:
            self._profiles.append(
                _ProfileEntry(
                    profile=profile,
                    label=label,
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                )
            )
        return self

    def add_superset(self, *supersets: Any, **override_options) -> "Plotter":
        """
        Add one or more supersets.

        Superset is expected to expose:
        - supermodes
        - interpret_combination
        - add_crossings_to_ax (optional, used for index and beta when show_crossings is True)

        Any override_options update PlotterOptions for this Plotter instance.
        """
        if override_options:
            self.options = PlotterOptions(
                **{**self.options.__dict__, **override_options}
            )

        for superset in supersets:
            self._supersets.append(superset)
            for mode in getattr(superset, "supermodes", []):
                self._supermodes_from_supersets[self._mode_key(mode)] = mode

        return self

    def add_supermode(self, *supermodes: SUPERMODE) -> "Plotter":
        """
        Add one or more standalone supermodes.

        Dedup rules:
        - If the supermode already exists inside an added superset, it is ignored.
        - If it was already added explicitly, it is ignored.
        """
        for mode in supermodes:
            key = self._mode_key(mode)

            if key in self._supermodes_from_supersets:
                continue

            already_added = any(
                self._mode_key(existing) == key for existing in self._extra_supermodes
            )
            if already_added:
                continue

            self._extra_supermodes.append(mode)

        return self

    def show(
        self,
        *,
        axes: Optional[plt.Axes] = None,
        show: bool = True,
        **override_options,
    ):
        """
        Render the configured plot.

        Parameters
        ----------
        axes : matplotlib.axes.Axes, optional
            Axis to draw on. If None, a figure and axis are created.
        show : bool
            If True, calls plt.show().
        **override_options
            Any PlotterOptions fields to override for this call.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the plot.
        """
        if override_options:
            self.options = PlotterOptions(
                **{**self.options.__dict__, **override_options}
            )

        context = plt.style.context(mps) if self.use_mps_style else _NullContext()

        with context:
            figure, axes = _ensure_axes(axes, figsize=(10, 4))

            if self.kind == "adiabatic":
                self._plot_adiabatic(axes)
            elif self.kind == "index":
                self._plot_scalar_representation(
                    axes,
                    representation_name="index",
                    y_label="Effective index",
                    crossings_key="index",
                )
            elif self.kind == "beta":
                self._plot_scalar_representation(
                    axes,
                    representation_name="beta",
                    y_label="Propagation constant",
                    crossings_key="beta",
                )
            elif self.kind == "eigenvalue":
                self._plot_scalar_representation(
                    axes,
                    representation_name="eigenvalue",
                    y_label="Eigenvalue",
                    crossings_key="eigen_value",
                )
            else:
                raise ValueError(f"Unknown plot kind: {self.kind}")

            if show:
                plt.show()

            return figure

    def _plot_adiabatic(self, axes: plt.Axes) -> None:
        """
        Plot adiabatic coupling curves.

        Sources:
        - For each added superset: plots pairwise adiabatic curves based on interpret_combination.
        - For explicitly added modes: plots all pairs among them.
        - For each added profile: overlays profile.adiabatic versus profile.itr_list automatically.
        """
        any_curve_plotted = False

        for superset in self._supersets:
            mode_list = interpret_mode_of_interest(
                superset=superset,
                mode_of_interest=self.options.mode_of_interest,
            )
            combination = superset.interpret_combination(
                mode_of_interest=mode_list,
                combination=self.options.combination,
            )

            for mode_0, mode_1 in combination:
                mode_0.adiabatic.plot(ax=axes, other_supermode=mode_1, show=False)
                any_curve_plotted = True

        if len(self._extra_supermodes) >= 2:
            for i in range(len(self._extra_supermodes)):
                for j in range(i + 1, len(self._extra_supermodes)):
                    mode_0 = self._extra_supermodes[i]
                    mode_1 = self._extra_supermodes[j]
                    mode_0.adiabatic.plot(ax=axes, other_supermode=mode_1, show=False)
                    any_curve_plotted = True

        for entry in self._profiles:
            profile = entry.profile

            axes.plot(
                profile.itr_list,
                profile.adiabatic,
                color=entry.color,
                linestyle=entry.linestyle,
                linewidth=entry.linewidth,
                label=(
                    entry.label
                    if entry.label is not None
                    else getattr(profile, "label", "profile")
                ),
            )
            any_curve_plotted = True

        if not any_curve_plotted:
            raise RuntimeError(
                "Nothing to plot for kind='adiabatic'. Add a superset, at least two supermodes, or a profile."
            )

        axes.set(xlabel="ITR", ylabel="Adiabatic criterion")
        axes.legend()

    def _plot_scalar_representation(
        self,
        axes: plt.Axes,
        *,
        representation_name: str,
        y_label: str,
        crossings_key: str,
    ) -> None:
        """
        Plot a scalar representation (index, beta, eigenvalue) from supersets and/or modes.

        This expects the representation object to implement:
        - plot(ax=..., show=False)
        """
        any_curve_plotted = False

        for superset in self._supersets:
            mode_list = interpret_mode_of_interest(
                superset=superset,
                mode_of_interest=self.options.mode_of_interest,
            )

            for mode in mode_list:
                getattr(mode, representation_name).plot(ax=axes, show=False)
                any_curve_plotted = True

            if self.options.show_crossings and hasattr(superset, "add_crossings_to_ax"):
                superset.add_crossings_to_ax(
                    ax=axes,
                    mode_of_interest=mode_list,
                    data_type=crossings_key,
                )

        for mode in self._extra_supermodes:
            getattr(mode, representation_name).plot(ax=axes, show=False)
            any_curve_plotted = True

        if not any_curve_plotted:
            raise RuntimeError(
                f"Nothing to plot for kind='{self.kind}'. Add a superset or at least one supermode."
            )

        axes.set(xlabel="ITR", ylabel=y_label)
        axes.legend()

    @staticmethod
    def _mode_key(mode: Any) -> tuple[Any, Any]:
        """
        Stable identifier for mode deduplication.

        Prefer mode.ID if available, else fallback to hash(mode).
        """
        if hasattr(mode, "ID"):
            try:
                return tuple(mode.ID)  # (mode_number, solver_number)
            except TypeError:
                pass
        return ("hash", hash(mode))


def _ensure_axes(axes: Optional[plt.Axes], *, figsize: tuple[float, float]):
    if axes is not None:
        return axes.figure, axes
    figure, axes = plt.subplots(1, 1, figsize=figsize)
    return figure, axes


class _NullContext:
    def __enter__(self):
        return None

    def __exit__(self, exc_type, exc, tb):
        return False


def generate_propagation_gif(
    profile: AlphaProfile,
    output_directory: str = "./new_gif.gif",
    dpi: int = 100,
    fps: int = 20,
    number_of_frames: int = 200,
    dark_background: bool = True,
) -> None:
    """
    Generates an animated GIF of light propagation in a taper structure.

    Parameters:
        output_directory (str): Path where the GIF will be saved.
        dpi (int): Dots per inch for the output GIF.
        fps (int): Frames per second for the animation.
        number_of_frames (int): Total number of frames in the animation.
        dark_background (bool): If True, use a dark background for the GIF.

    Returns:
        None
    """
    figure, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.set_xlabel(
        "Propagation axis [mm]", color="white" if dark_background else "black"
    )
    style_context = "dark_background" if dark_background else "default"

    sub_sampling_factor = int(profile.distance.size / number_of_frames)

    sub_distance = profile.distance[::sub_sampling_factor] * 1e3
    sub_radius = profile.radius[::sub_sampling_factor]
    sub_itr_list = profile.itr_list[::sub_sampling_factor]
    with plt.style.context(style_context):

        def init_func() -> tuple:
            line_0 = ax.plot(sub_distance, sub_radius, color="black")
            line_1 = ax.plot(sub_distance, -sub_radius, color="black")

            line_2 = ax.fill_between(
                sub_distance, +sub_radius, -sub_radius, color="lightblue", alpha=0.8
            )

            return [*line_0, *line_1, line_2]

        def animate(slice_number: int) -> tuple:
            position = sub_distance[slice_number]
            itr = sub_itr_list[slice_number]
            title = f"[slice: {slice_number} - ITR: {itr:.3f}]"

            if slice_number > 0:
                ax.lines[-1].remove()

            line_0 = ax.set_title(title, color="white")

            line_1 = ax.axvline(position, linestyle="--", linewidth=2, color="red")

            return line_0, line_1

    animation = FuncAnimation(
        fig=figure,
        func=animate,
        init_func=init_func,
        blit=True,
        repeat=True,
        frames=sub_itr_list.size,
    )

    animation.save(output_directory, dpi=dpi, writer=PillowWriter(fps=fps))


def get_3_figures(
    *,
    figure_size: Tuple[float, float] = (8, 6),
    left: float = 0.1,
    right: float = 0.9,
    bottom: float = 0.1,
    top: float = 0.9,
    wspace: float = 1.0,
    hspace: float = 0.15,
):
    """
    Create a 3-axes figure layout:
    - top row: one axis spanning both columns
    - bottom row: two square axes (images)

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax0 : matplotlib.axes.Axes
        Top axis spanning full width.
    ax1 : matplotlib.axes.Axes
        Bottom-left image axis.
    ax2 : matplotlib.axes.Axes
        Bottom-right image axis.
    """
    fig = plt.figure(figsize=figure_size)

    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(1, 1),
        height_ratios=(1, 1),
        left=left,
        right=right,
        bottom=bottom,
        top=top,
        wspace=wspace,
        hspace=hspace,
    )

    ax0 = fig.add_subplot(gs[0, :])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])

    for ax in (ax1, ax2):
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(
            axis="both",
            which="both",
            top=False,
            bottom=False,
            right=False,
            left=False,
            labelleft=False,
            labelbottom=False,
            grid_alpha=0,
        )

    return fig, ax0, ax1, ax2


@dataclass(frozen=True)
class _AnimationSettings:
    n_step: int = 20
    dpi: float = 100.0
    fps: int = 50
    interval_ms: int = 40
    outfile: str = "propagation.gif"


class ModePropagationGifCreator:
    """
    Create an animated GIF of mode fields along a taper profile.

    The API is intentionally close to your original implementation:
    - constructor takes (superset, profile, max_number_of_mode, dark_background)
    - make_animation(n_step, dpi, fps) writes "propagation.gif" by default

    Improvements
    ------------
    - Stable artists: uses persistent QuadMesh objects and updates their data,
      avoiding ax.clear() which is slow and breaks blitting.
    - Robust figure layout: fields are created from a single method and stored.
    - Profile marker: uses one Line2D updated each frame.
    - Optional progress printing without spamming by default.

    Notes
    -----
    This expects:
    - profile.distance, profile.radius, profile.length
    - profile.master_interpolation_z_to_itr(z) -> float
    - superset.get_slice_structure(itr=..., add_symmetries=True) -> object with .fields
    - superset[mode].stylized_name
    """

    def __init__(
        self,
        superset,
        profile,
        max_number_of_mode: Optional[int] = None,
        dark_background: bool = True,
    ):
        self.superset = superset
        self.profile = profile
        self.dark_background = dark_background

        self.number_of_mode = (
            len(superset.supermodes)
            if max_number_of_mode is None
            else int(max_number_of_mode)
        )

        self.figure = None
        self.grid_spec = None
        self.ax_profile = None
        self.field_axes: list = []
        self._field_meshes: list = []
        self._profile_line = None

        self.generate_figure()

    # -------------------------
    # Figure construction
    # -------------------------
    def generate_figure(self, unit_size: Tuple[float, float] = (3, 6)) -> None:
        figure_size = (unit_size[0] * self.number_of_mode, unit_size[1])
        self.figure = plt.figure(figsize=figure_size)

        self.grid_spec = self.figure.add_gridspec(
            2,
            self.number_of_mode,
            left=0.1,
            right=0.95,
            bottom=0.1,
            top=0.9,
            wspace=0.5,
            hspace=0.35,
        )

        self._generate_profile_ax()
        self._generate_field_axes()

    def _generate_profile_ax(self) -> None:
        self.ax_profile = self.figure.add_subplot(self.grid_spec[0, :])
        self.ax_profile.set_xlabel("Propagation distance z")
        self.ax_profile.set_ylabel("Coupler profile")

        top = self.profile.radius
        bottom = -self.profile.radius

        self.ax_profile.plot(self.profile.distance, top, color="black")
        self.ax_profile.plot(self.profile.distance, bottom, color="black")
        self.ax_profile.fill_between(
            self.profile.distance, top, bottom, color="lightblue", alpha=0.8
        )

        self._profile_line = self.ax_profile.axvline(0.0, linestyle="--", color="red")

    def _generate_field_axes(self) -> None:
        self.field_axes = []
        self._field_meshes = []

        for mode_index in range(self.number_of_mode):
            ax = self.figure.add_subplot(self.grid_spec[1, mode_index])
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(self.superset[mode_index].stylized_name)

            # Placeholder mesh. We will replace data on first populate.
            mesh = ax.pcolormesh(
                [[0.0, 0.0], [0.0, 0.0]], cmap=colormaps.blue_black_red, shading="auto"
            )

            self.field_axes.append(ax)
            self._field_meshes.append(mesh)

    # -------------------------
    # Frame evaluation
    # -------------------------
    def _get_fields_at_z(self, z: float):
        itr = self.profile.master_interpolation_z_to_itr(z)
        slice_structure = self.superset.get_slice_structure(
            itr=itr, add_symmetries=True
        )
        return itr, slice_structure.fields

    def populate_axes(self, z: float) -> None:
        itr, fields = self._get_fields_at_z(z)

        self._profile_line.set_xdata([z, z])
        self.ax_profile.set_title(f"Z-distance: {z:>5.3e}    ITR: {itr:>5.3f}")

        for mesh, field in zip(self._field_meshes, fields):
            # QuadMesh stores a flattened array; set_array expects 1D.
            mesh.set_array(field.ravel())

        self.figure.canvas.draw_idle()

    def update_axes(self, z: float) -> Sequence:
        itr, fields = self._get_fields_at_z(z)

        self._profile_line.set_xdata([z, z])
        self.ax_profile.set_title(f"Z-distance: {z:>5.3e}    ITR: {itr:>5.3f}")

        updated_artists = [self._profile_line]

        for mesh, field in zip(self._field_meshes, fields):
            mesh.set_array(field.ravel())
            updated_artists.append(mesh)

        return updated_artists

    # -------------------------
    # Animation
    # -------------------------
    def make_animation(self, n_step: int = 20, dpi: float = 100, fps: int = 50):
        settings = _AnimationSettings(n_step=n_step, dpi=float(dpi), fps=int(fps))

        self.populate_axes(z=0.0)

        z_max = float(self.profile.length)
        step = z_max / float(settings.n_step)

        def animate(frame_index: int):
            z = frame_index * step
            return self.update_axes(z)

        animation = FuncAnimation(
            self.figure,
            animate,
            interval=settings.interval_ms,
            blit=True,
            repeat=True,
            frames=settings.n_step,
        )

        animation.save(
            settings.outfile,
            dpi=settings.dpi,
            writer=PillowWriter(fps=settings.fps),
        )
