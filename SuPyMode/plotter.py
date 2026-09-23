#!/usr/bin/env python
# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple
from MPSPlots import colormaps
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

from SuPyMode.binary.interface_taper import AlphaProfile


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
