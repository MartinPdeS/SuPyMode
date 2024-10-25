#!/usr/bin/env python
# -*- coding: utf-8 -*-
from FiberFusing.fiber import catalogue as fiber_catalogue  # noqa:
from SuPyMode.profiles import AlphaProfile  # noqa: F401
from FiberFusing import configuration  # noqa: F401

from typing import List, Union, Optional, Tuple
from pathlib import Path
from FiberFusing import Geometry, BackGround
from FiberFusing.fiber.generic_fiber import GenericFiber
from SuPyMode.solver import SuPySolver
from PyOptik import MaterialBank
from PyFinitDiff.finite_difference_2D import Boundaries
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict

config_dict = ConfigDict(
    extra='forbid',
    strict=True,
    arbitrary_types_allowed=True,
    kw_only=True,
)


def prepare_simulation_geometry(
        wavelength: float,
        clad_structure: object,
        fiber_list: List[object],
        capillary_tube: Optional[object] = None,
        fusion_degree: Union[float, str] = 'auto',
        fiber_radius: Optional[float] = None,
        x_bounds: Union[str, Tuple[float, float]] = '',
        y_bounds: Union[str, Tuple[float, float]] = '',
        clad_index: Union[float, str] = 'silica',
        core_position_scrambling: float = 0,
        index_scrambling: float = 0,
        resolution: int = 150,
        rotation: float = 0,
        air_padding_factor: float = 1.2,
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

    Raises:
        ValueError: If clad_index is not 'silica' or a numeric value.
    """

    def get_clad_index(clad_index: float, wavelength: float):
        """Retrieve the cladding index based on the input type."""
        if isinstance(clad_index, str) and clad_index.lower() == 'silica':
            return MaterialBank.fused_silica.compute_refractive_index(wavelength)
        elif isinstance(clad_index, (float, int)):
            return clad_index
        else:
            raise ValueError("Invalid clad_index: must be either 'silica' or a numeric index value.")

    index = get_clad_index(clad_index, wavelength)
    background = BackGround(index=background_index)

    clad_instance = prepare_fused_structure(
        clad_class=clad_structure,
        fiber_radius=fiber_radius,
        fusion_degree=fusion_degree,
        index=index,
        core_position_scrambling=core_position_scrambling,
        rotation=rotation
    )

    structures = []
    if capillary_tube is not None:
        structures.append(capillary_tube)
    if clad_instance is not None:
        structures.append(clad_instance)

    if clad_instance is not None:
        for fiber, core in zip(fiber_list, clad_instance.cores):
            fiber.set_position((core.x, core.y))

    structures.extend(fiber_list)

    geometry = Geometry(
        background=background,
        x_bounds=x_bounds,
        additional_structure_list=structures,
        y_bounds=y_bounds,
        resolution=resolution,
        index_scrambling=index_scrambling,
        boundary_pad_factor=air_padding_factor,
        gaussian_filter=gaussian_filter
    )

    return geometry


def prepare_fused_structure(
        clad_class: type,
        fiber_radius: float,
        fusion_degree: Union[float | str],
        index: float,
        core_position_scrambling: float,
        rotation: float) -> object:
    """
    Prepares and returns a clad instance according to the provided clad class and configuration parameters.

    Args:
        clad_class (type): The class representing the cladding structure.
        fiber_radius (float): The radius of the fiber.
        fusion_degree (float): The degree of fusion for the fibers.
        index (float): The refractive index of the cladding material.
        core_position_scrambling (float): Random displacement added to fiber core positions to simulate imperfections.
        rotation (float): Rotation angle for the structure (in degrees).

    Returns:
        Optional[object]: An instance of the clad structure configured with the provided parameters, or None if no clad class is provided.

    Raises:
        ValueError: If any of the required parameters are not provided or invalid.
    """
    if clad_class is None:
        return None

    clad_instance = clad_class(
        fiber_radius=fiber_radius,
        fusion_degree=fusion_degree,
        index=index
    )

    clad_instance.randomize_core_position(random_factor=core_position_scrambling)

    if rotation != 0:
        clad_instance.rotate(rotation)

    return clad_instance


@dataclass(config=config_dict)
class Workflow():
    """
    Configures and executes optical simulations using finite difference methods on specified fiber geometries.

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
        fiber_list (List[GenericFiber]): List of fibers included in the optical structure.
        boundaries (List[Boundaries]): Boundary conditions applied to the simulation.
        capillary_tube (Optional[object]): Additional capillary structure to include in the simulation.
        clad_structure (Optional[object]): Initial optical structure used for the simulation.
        x_bounds (Union[str, List[float]]): Boundary conditions along the x-axis.
        y_bounds (Union[str, List[float]]): Boundary conditions along the y-axis.
        air_padding_factor (float): Factor for padding the structure with air to prevent boundary effects.
        gaussian_filter_factor (Optional[float]): Gaussian blurring factor applied to the structure for smoothing.

    Plotting Flags:
        plot_geometry (bool): If True, plots the simulation geometry.
        plot_cladding (bool): If True, plots the cladding structure.
        plot_field (bool): If True, plots the field distribution.
        plot_adiabatic (bool): If True, plots adiabatic transitions.
        plot_coupling (bool): If True, plots coupling information.
        plot_beating_length (bool): If True, plots the beating length.
        plot_eigen_values (bool): If True, plots eigenvalues.
        plot_index (bool): If True, plots the refractive index distribution.
        plot_beta (bool): If True, plots propagation constants.

    Methods:
        __post_init__: Initializes the simulation geometry and solver upon object creation.
        plot: Plots various simulation outputs based on the provided plot type.
        save_superset_instance: Saves the computed superset to a file for later use.
        generate_pdf_report: Generates a comprehensive PDF report of all relevant simulation data and results.
    """

    # Geometry attributes
    wavelength: float
    boundaries: List[Boundaries]
    resolution: int
    clad_rotation: Optional[float] = 0
    capillary_tube: Optional[object] = None
    clad_structure: Optional[object] = None
    fiber_list: Optional[List[GenericFiber]] = ()
    fiber_radius: Optional[float] = 62.5e-6
    fusion_degree: Optional[Union[float, str]] = 'auto'
    x_bounds: Union[str, Tuple[float, float]] = 'centering'
    y_bounds: Union[str, Tuple[float, float]] = 'centering'
    air_padding_factor: Optional[float] = 1.2
    gaussian_filter_factor: Optional[float] = None

    # Solver attributes
    n_sorted_mode: int = 4
    n_step: int = 500
    itr_final: float = 0.05
    n_added_mode: Optional[int] = 4
    itr_initial: Optional[float] = 1.0
    extrapolation_order: Optional[int] = 2
    core_position_scrambling: Optional[float] = 0
    index_scrambling: Optional[float] = 0
    accuracy: Optional[int] = 2

    # Plotting flags
    plot_geometry: Optional[bool] = False
    plot_cladding: Optional[bool] = False
    plot_field: Optional[bool] = False
    plot_adiabatic: Optional[bool] = False
    plot_coupling: Optional[bool] = False
    plot_beating_length: Optional[bool] = False
    plot_eigen_values: Optional[bool] = False
    plot_index: Optional[bool] = False
    plot_beta: Optional[bool] = False

    # Extra attributes
    debug_mode: Optional[int] = 1
    auto_label: Optional[bool] = False
    generate_report: Optional[bool] = False
    save_superset: Optional[bool] = False

    def __post_init__(self):
        """
        Initializes the simulation geometry and solver, and optionally plots the initial setup if enabled.
        """
        self.geometry = self.prepare_simulation_geometry()
        self._plot_initial_setup()
        self._initialize_solver()
        self._plot_simulation_outputs()
        self._finalize_workflow()

    def prepare_simulation_geometry(self):
        """Prepares the simulation geometry."""
        return prepare_simulation_geometry(
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
            air_padding_factor=self.air_padding_factor,
            gaussian_filter=self.gaussian_filter_factor
        )

    def _plot_initial_setup(self):
        """Plots the initial simulation setup if the respective flags are enabled."""
        if self.plot_cladding and self.clad_structure:
            self.clad_structure.plot()
        if self.plot_geometry:
            self.geometry.plot()

    @property
    def superset(self):
        return self.solver.superset

    def _initialize_solver(self):
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

    def get_superset(self):
        return self.solver.superset

    def _plot_simulation_outputs(self):
        """Plots the simulation outputs based on the respective flags."""
        if self.plot_field:
            self.plot('field')
        if self.plot_adiabatic:
            self.plot('adiabatic')
        if self.plot_coupling:
            self.plot('normalized-coupling')
        if self.plot_beating_length:
            self.plot('beating-length')
        if self.plot_eigen_values:
            self.plot('eigen-value')
        if self.plot_index:
            self.plot('index')
        if self.plot_beta:
            self.plot('beta')

    def _finalize_workflow(self):
        """Finalizes the workflow by generating reports and saving superset if enabled."""
        if self.generate_report:
            self.generate_pdf_report()
        if self.save_superset:
            self.save_superset_instance()

    def plot(self, *args, **kwargs):
        """Plots various types of data from the simulation based on the specified plot type."""
        return self.solver.superset.plot(*args, **kwargs)

    def save_superset_instance(self, filename: str = 'auto', directory: str = 'auto') -> Path:
        """
        Saves the superset instance to a file, defaulting to an auto-generated filename if not specified.

        Args:
            filename (str): Filename for the saved instance.
            directory (str): Directory for saving the instance.

        Returns:
            Path: Path to the saved file.
        """
        self.solver.superset.save_instance(
            filename=filename,
            directory=directory
        )

        return filename

    def generate_pdf_report(self, filename: str = 'auto', directory: str = 'auto', **kwargs) -> Path:
        """
        Generates a PDF report of all relevant simulation data and results.

        Args:
            filename (str): Filename for the report.
            **kwargs: Additional arguments for the report generation.

        Returns:
            Path: Path to the generated PDF report.
        """
        self.solver.superset.generate_pdf_report(
            filename=filename,
            directory=directory,
            **kwargs
        )

        return filename

# -
