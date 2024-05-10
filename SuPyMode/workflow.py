#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import List, Union, Optional
from dataclasses import dataclass
from pathlib import Path

from FiberFusing import Geometry, BackGround
from SuPyMode.solver import SuPySolver
from FiberFusing.fiber import catalogue as fiber_catalogue
from SuPyMode.profiles import AlphaProfile  # noqa: F401
from FiberFusing import configuration  # noqa: F401

from PyFinitDiff.finite_difference_2D import Boundaries
from pathvalidate import sanitize_filepath


def prepare_simulation_geometry(
        wavelength: float,
        clad_structure: object,
        fiber_list: List[object],
        capillary_tube: Optional[object] = None,
        fusion_degree: Union[float, str] = 'auto',
        fiber_radius: Optional[float] = None,
        x_bounds: Union[str, List[float]] = '',
        y_bounds: Union[str, List[float]] = '',
        clad_index: Union[float, str] = 'silica',
        core_position_scrambling: float = 0,
        index_scrambling: float = 0,
        resolution: int = 150,
        rotation: float = 0,
        boundary_pad_factor: float = 1.2,
        gaussian_filter: float = 0,
        background_index: float = 1) -> Geometry:
    """
    Prepares and returns the simulation geometry for optical fiber configurations,
    incorporating fused structures, optional capillary tubes, and specific boundary conditions.

    Args:
        wavelength (float): Wavelength for refractive index calculation.
        clad_structure (object): Cladding or structural template for fibers.
        fiber_list (List[object]): List of fibers to be included in the simulation.
        capillary_tube (Optional[object]): Optional capillary structure to add.
        fusion_degree (Union[float, str]): Degree of fusion, specifying overlap between fibers.
        fiber_radius (Optional[float]): Uniform radius for all fibers, if specified.
        x_bounds (Union[str, List[float]]): X-axis boundary conditions or limits.
        y_bounds (Union[str, List[float]]): Y-axis boundary conditions or limits.
        clad_index (Union[float, str]): Refractive index for the cladding material, can be a known string identifier.
        core_position_scrambling (float): Random displacement added to fiber core positions to simulate imperfections.
        index_scrambling (float): Noise level to simulate index inhomogeneity.
        resolution (int): Resolution of the geometry mesh grid.
        rotation (float): Rotation angle for the structure (degrees).
        boundary_pad_factor (float): Padding factor for boundary adjustments.
        gaussian_filter (float): Standard deviation for Gaussian blur to smooth sharp transitions.
        background_index (float): Background refractive index for the simulation area.

    Returns:
        Geometry: Configured geometry object ready for simulation.
    """
    if isinstance(clad_index, str) and clad_index.lower() == 'silica':
        index = fiber_catalogue.get_silica_index(wavelength=wavelength)
    elif isinstance(clad_index, (float, int)):
        index = clad_index
    else:
        raise ValueError("Invalid clad_index: must be either 'silica' or a numeric index value.")

    background = BackGround(index=background_index)

    clad_instance = prepare_fused_structure(
        clad_class=clad_structure,
        fiber_radius=fiber_radius,
        fusion_degree=fusion_degree,
        index=index,
        core_position_scrambling=core_position_scrambling,
        rotation=rotation
    )

    geometry = Geometry(
        background=background,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        resolution=resolution,
        index_scrambling=index_scrambling,
        boundary_pad_factor=boundary_pad_factor,
        gaussian_filter=gaussian_filter
    )

    if capillary_tube is not None:
        geometry.add_structure(capillary_tube)

    if clad_instance is not None:
        geometry.add_structure(clad_instance)

    if clad_instance is not None:
        for fiber, core in zip(fiber_list, clad_instance.cores):
            fiber.set_position(core)

    geometry.add_fiber(*fiber_list)

    return geometry


def prepare_fused_structure(
        clad_class: type,
        fiber_radius: float,
        fusion_degree: float,
        index: float,
        core_position_scrambling: float,
        rotation: float) -> object:
    """
    Prepare and returns a clad instance according to the clad class given as input.

    :param      clad_class:                The clad class
    :type       clad_class:                type
    :param      fiber_radius:              The fiber radius
    :type       fiber_radius:              float
    :param      fusion_degree:             The fusion degree
    :type       fusion_degree:             float
    :param      index:                     The index
    :type       index:                     float
    :param      core_position_scrambling:  The core position scrambling
    :type       core_position_scrambling:  float
    :param      rotation:                  The rotation
    :type       rotation:                  float

    :returns:   The clad instance
    :rtype:     object
    """
    if clad_class is None:
        return None

    clad_instance = clad_class(
        fiber_radius=fiber_radius,
        fusion_degree=fusion_degree,
        index=index,
        core_position_scrambling=core_position_scrambling
    )

    if rotation != 0:
        clad_instance.rotate(rotation)

    return clad_instance


@dataclass
class Workflow():
    """
    A class to configure and execute optical simulations using finite difference methods on specified fiber geometries.

    Attributes:
        wavelength (float): Wavelength at which the simulation is evaluated.
        resolution (int): Resolution of the simulation mesh.
        fiber_radius (float): Radius for the fused clad structure.
        n_sorted_mode (int): Number of modes that are computed and sorted.
        n_added_mode (int): Additional modes computed beyond the sorted modes for increased accuracy.
        itr_final (float): Final Inverse Taper Ratio (ITR) for mode evaluation.
        itr_initial (float): Initial ITR for mode evaluation.
        n_step (int): Number of steps to iterate through the ITR section.
        fusion_degree (Union[float, str]): Fusion degree for the clad fused structure; 'auto' for automatic calculation.
        clad_rotation (float): Rotation of the clad structure in degrees.
        accuracy (int): Accuracy level of the finite difference method.
        debug_mode (int): Debug mode level for verbose output during computations.
        auto_label (bool): If True, automatically labels supermodes.
        generate_report (bool): If True, generates a PDF report of the simulation results.
        save_superset (bool): If True, saves the computed superset instance for later use.
        fiber_list (List[object]): List of fibers included in the optical structure.
        boundaries (List[Boundaries]): Boundary conditions applied to the simulation.
        capillary_tube (Optional[object]): Additional capillary structure to include in the simulation.
        clad_structure (Optional[object]): Initial optical structure used for the simulation.
        x_bounds (str): Boundary conditions along the x-axis.
        y_bounds (str): Boundary conditions along the y-axis.
        air_padding_factor (float): Factor for padding the structure with air to prevent boundary effects.
        gaussian_filter_factor (Optional[float]): Gaussian blurring factor applied to the structure for smoothing.

    Plotting Flags:
        Various flags that determine which aspects of the simulation are plotted.

    Methods:
        __post_init__: Initializes the simulation geometry and solver upon object creation.
        plot: Plots various simulation outputs based on the provided plot type.
        save_superset_instance: Saves the computed superset to a file for later use.
        generate_pdf_report: Generates a comprehensive PDF report of all relevant simulation data and results.
        _get_auto_generated_filename: Generates a filename based on the simulation parameters.
    """

    #  Geometry arguments --------------------------
    wavelength: float
    clad_rotation: float = 0
    capillary_tube: object = None
    resolution: int = 100
    clad_structure: object = None
    fiber_list: list = tuple()
    fiber_radius: float = 62.5e-6
    fusion_degree: float = 'auto'
    x_bounds: str = ''
    y_bounds: str = ''
    air_padding_factor: float = 1.2
    gaussian_filter_factor: float = None

    #  Solver arguments --------------------------
    n_sorted_mode: int = 4
    n_added_mode: int = 4
    itr_final: float = 0.05
    itr_initial: float = 1.0
    n_step: int = 500
    extrapolation_order: int = 2
    core_position_scrambling: float = 0
    index_scrambling: float = 0
    boundaries: list = (Boundaries(),)
    accuracy: int = 2

    #  Plot arguments --------------------------
    plot_geometry: bool = False
    plot_cladding: bool = False
    plot_field: bool = False
    plot_adiabatic: bool = False
    plot_coupling: bool = False
    plot_beating_length: bool = False
    plot_eigen_values: bool = False
    plot_index: bool = False
    plot_beta: bool = False

    #  Extra arguments --------------------------
    debug_mode: int = 1
    auto_label: bool = False
    generate_report: bool = False
    save_superset: bool = False

    def __post_init__(self):
        """Initializes the simulation geometry and solver, and optionally plots the initial setup if enabled."""
        self.geometry = prepare_simulation_geometry(
            wavelength=self.wavelength,
            clad_structure=self.clad_structure,
            fiber_list=self.fiber_list,
            capillary_tube=self.capillary_tube,
            fusion_degree=self.fusion_degree,
            fiber_radius=self.fiber_radius,
            resolution=self.resolution,
            y_bounds=self.y_bounds,
            x_bounds=self.x_bounds,
            rotation=self.clad_rotation,
            gaussian_filter=self.gaussian_filter_factor
        )

        if self.plot_cladding:
            self.clad_structure.plot().show()

        if self.plot_geometry:
            self.geometry.plot().show()

        self._initialize_solver_()

        if self.plot_field:
            self.plot(plot_type='field').show()

        if self.plot_adiabatic:
            self.plot(plot_type='adiabatic').show()

        if self.plot_coupling:
            self.plot(plot_type='normalized-coupling').show()

        if self.plot_beating_length:
            self.plot(plot_type='beating-length').show()

        if self.plot_eigen_values:
            self.plot(plot_type='eigen-value').show()

        if self.plot_index:
            self.plot(plot_type='index').show()

        if self.plot_beta:
            self.plot(plot_type='beta').show()

        if self.generate_report:
            self.generate_pdf_report()

        if self.save_superset:
            self.save_superset_instance()

    @property
    def superset(self):
        return self.solver.superset

    def _initialize_solver_(self) -> None:
        """Initializes the solver with the set geometry and starts the mode computation process."""

        self.solver = SuPySolver(
            geometry=self.geometry,
            tolerance=1e-20,
            max_iter=5000,
            accuracy=self.accuracy,
            debug_mode=self.debug_mode,
            extrapolation_order=self.extrapolation_order
        )

        self.solver.init_superset(
            wavelength=self.wavelength,
            n_step=self.n_step,
            itr_initial=self.itr_initial,
            itr_final=self.itr_final
        )

        for boundary in self.boundaries:
            self.solver.add_modes(
                n_added_mode=self.n_added_mode,
                n_sorted_mode=self.n_sorted_mode,
                boundaries=boundary,
                auto_label=self.auto_label
            )

        self.solver.superset.sort_modes(sorting_method='beta')

    def get_superset(self):
        return self.solver.superset

    def _get_auto_generated_filename_(self) -> str:
        """
        Returns an auton-generated filename taking account for:
        |  structure name
        |  fiber names
        |  resolution of the simulations
        |  wavelength

        :returns:   The automatic generated filename.
        :rtype:     str
        """
        fiber_name = "".join(fiber.__class__.__name__ for fiber in self.fiber_list)

        filename = (
            f"structure={self.clad_structure.__class__.__name__}_"
            f"{fiber_name}_"
            f"resolution={self.resolution}_"
            f'wavelength={self.wavelength}'
        )

        return filename.replace('.', '_')

    def save_superset_instance(self, filename: str = 'auto', directory: str = 'auto') -> Path:
        """Saves the superset instance to a file, defaulting to an auto-generated filename if not specified."""
        if filename == 'auto':
            filename = self._get_auto_generated_filename_()

        filename = Path(filename + '.pdf')

        filename = sanitize_filepath(filename)

        self.solver.superset.save_instance(
            filename=filename,
            directory=directory
        )

        return filename

    def generate_pdf_report(self, filename: str = 'auto', **kwargs) -> Path:
        """
        Generate a pdf file compiling the essential computed components.

        :param      filename:  The filename
        :type       filename:  str
        :param      kwargs:    The keywords arguments
        :type       kwargs:    dictionary

        :returns:   The path directory of the report
        :rtype:     Path
        """
        if filename == 'auto':
            filename = self._get_auto_generated_filename_()

        filename = Path(filename).with_suffix('.pdf')

        filename = sanitize_filepath(filename)

        self.solver.superset.generate_pdf_report(filename=filename, **kwargs)

        return filename

    def plot(self, *args, **kwargs):
        """Plots various types of data from the simulation based on the specified plot type."""

        return self.solver.superset.plot(*args, **kwargs)


# -
