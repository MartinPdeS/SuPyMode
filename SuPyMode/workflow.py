#!/usr/bin/env python
# -*- coding: utf-8 -*-

from FiberFusing.fiber import FiberLoader  # noqa:
from SuPyMode.profiles import AlphaProfile  # noqa: F401
from FiberFusing.profile import Profile, StructureType # noqa: F401
from FiberFusing.graded_index import GradedIndex  # noqa: F401
from PyFinitDiff import BoundaryValue  # noqa: F401

from FiberFusing import DomainAlignment
from typing import List, Union, Optional, Tuple
from pathlib import Path
from FiberFusing import Geometry, BackGround
from FiberFusing.fiber.generic_fiber import GenericFiber
from SuPyMode.solver import SuPySolver
from PyFinitDiff.finite_difference_2D import Boundaries
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict

fiber_loader = FiberLoader()

config_dict = ConfigDict(
    extra='forbid',
    strict=True,
    arbitrary_types_allowed=True,
    kw_only=True,
)


@dataclass(config=config_dict)
class Workflow():
    """
    Configures and executes optical simulations using finite difference methods on specified fiber geometries.

    Parameters
    ----------
    wavelength : float
        Wavelength at which the simulation is evaluated.
    resolution : int
        Resolution of the simulation mesh.
    n_sorted_mode : int
        Number of modes that are computed and sorted.
    n_added_mode : int
        Additional modes computed beyond the sorted modes for increased accuracy.
    itr_final : float
        Final Inverse Taper Ratio (ITR) for mode evaluation.
    itr_initial : float
        Initial ITR for mode evaluation.
    n_step : int
        Number of steps to iterate through the ITR section.
    clad_rotation : float
        Rotation of the clad structure in degrees.
    accuracy : int
        Accuracy level of the finite difference method.
    debug_mode : int
        Debug mode level for verbose output during computations.
    auto_label : bool
        If True, automatically labels supermodes.
    generate_report : bool
        If True, generates a PDF report of the simulation results.
    save_superset : bool
        If True, saves the computed superset instance for later use.
    fiber_list : List[GenericFiber])
        List of fibers included in the optical structure.
    boundaries : List[Boundaries]
        Boundary conditions applied to the simulation.
    capillary_tube : Optional[object]
        Additional capillary structure to include in the simulation.
    clad_structure : Optional[object]
        Initial optical structure used for the simulation.
    x_bounds : Union[DomainAlignment, List[float]])
        Boundary conditions along the x-axis.
    y_bounds : Union[DomainAlignment, List[float]])
        Boundary conditions along the y-axis.
    air_padding_factor : float
        Factor for padding the structure with air to prevent boundary effects.
    gaussian_filter_factor : Optional[float]
        Gaussian blurring factor applied to the structure for smoothing.

    """

    # Geometry attributes
    wavelength: float
    boundaries: List[Boundaries]
    resolution: int
    clad_rotation: Optional[float] = 0
    capillary_tube: Optional[object] = None
    clad_structure: Optional[object] = None
    background_index: Optional[float] = 1.0
    fiber_list: Optional[List[GenericFiber]] = ()
    x_bounds: Union[DomainAlignment, Tuple[float, float]] = DomainAlignment.CENTERING
    y_bounds: Union[DomainAlignment, Tuple[float, float]] = DomainAlignment.CENTERING
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
    index_scrambling: Optional[float] = 0.
    accuracy: Optional[int] = 2

    # Extra attributes
    debug_mode: Optional[int] = 1
    auto_label: Optional[bool] = False

    @property
    def superset(self):
        return self.solver.superset

    def run_solver(self):
        """Initializes the solver with the set geometry and starts the mode computation process."""
        self.solver = SuPySolver(
            geometry=self.geometry,
            tolerance=1e-20,
            max_iteration=5000,
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

    def plot(self, *args, **kwargs):
        """Plots various types of data from the simulation based on the specified plot type."""
        return self.solver.superset.plot(*args, **kwargs)

    def save_superset_instance(self, filename: str = 'auto', directory: str = 'auto') -> Path:
        """
        Saves the superset instance to a file, defaulting to an auto-generated filename if not specified.

        Parameters
        ----------
        filename : str
            The name of the file to save the report as. Defaults to 'auto', which generates a timestamped filename.
        directory : str
            The directory where the report will be saved. Defaults to 'auto', which uses the current

        Returns
        -------
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

        Parameters
        ----------
        filename : str
            The name of the file to save the report as. Defaults to 'auto', which generates a timestamped filename.
        directory : str
            The directory where the report will be saved. Defaults to 'auto', which uses the current

        Returns
        -------
        Path: Path to the generated PDF report.
        """
        self.solver.superset.generate_pdf_report(
            filename=filename,
            directory=directory,
            **kwargs
        )

        return filename

    def initialize_geometry(self, plot: bool = False) -> Geometry:
        """
        Prepares and returns the simulation geometry for optical fiber configurations,
        incorporating fused structures, optional capillary tubes, and specific boundary conditions.

        Parameters
        ----------
        plot_geometry : bool
            Flag indicating whether to plot the geometry after initialization.

        Returns
        -------
        Geometry: Configured geometry object ready for simulation.
        """
        background = BackGround(refractive_index=self.background_index)

        self.prepare_fused_structure()

        structures = []
        if self.capillary_tube is not None:
            structures.append(self.capillary_tube)
        if self.clad_structure is not None:
            structures.append(self.clad_structure)

        if self.index_scrambling != 0.:
            for fiber in self.fiber_list:
                fiber.randomize_refractive_index(factor=self.index_scrambling)

        structures.extend(self.fiber_list)

        self.geometry = Geometry(
            x_bounds=self.x_bounds,
            y_bounds=self.y_bounds,
            resolution=self.resolution,
            index_scrambling=self.index_scrambling,
            boundary_pad_factor=self.air_padding_factor,
            gaussian_filter=self.gaussian_filter_factor
        )

        self.geometry.add_structure(background, *structures)

        self.geometry.initialize()

        if plot:
            self.geometry.plot()

    def prepare_fused_structure(self) -> None:
        """
        Prepares and returns a clad instance according to the provided clad class and configuration parameters.
        """
        if self.clad_structure is None:
            return

        self.clad_structure.randomize_core_positions(random_factor=self.core_position_scrambling)

        if self.clad_rotation != 0:
            self.clad_structure.rotate(self.clad_rotation)

# -
