#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from pathlib import Path

from FiberFusing import Geometry, BackGround
from SuPyMode.solver import SuPySolver
from FiberFusing.fiber import catalogue as fiber_catalogue
from SuPyMode.profiles import AlphaProfile
from FiberFusing import configuration

from PyFinitDiff.boundaries import Boundaries2D
from pathvalidate import sanitize_filepath


def prepare_simulation_geometry(
        wavelength: float,
        clad_structure: object,
        fiber_list: list,
        capillary_tube: object = None,
        fusion_degree: float = None,
        fiber_radius: float = None,
        x_bounds: str | list = '',
        y_bounds: str | list = '',
        clad_index: float | str = 'silica',
        core_position_scrambling: float = 0,
        index_scrambling: float = 0,
        resolution: int = 150,
        rotation: float = 0,
        boundary_pad_factor: float = 1.2,
        gaussian_filter: float = 0,
        background_index: float = 1) -> Geometry:
    """
    Prepare and returns the processed geometry for simulation using SuPyMode.

    :param      wavelength:                Wavelength at which evaluate the computation
    :type       wavelength:                float
    :param      clad_structure:            Initial optical structure
    :type       clad_structure:            object
    :param      fiber_list:                The fiber list
    :type       fiber_list:                list
    :param      capillary_tube:            Additional optical structure such as clad to add
    :type       capillary_tube:            object
    :param      fusion_degree:             Fusion degree for the clad fused structure
    :type       fusion_degree:             float
    :param      fiber_radius:              The fiber radius for the fused clad structure, all radii are assumed the same here.
    :type       fiber_radius:              float
    :param      x_bounds:                  The x-axis boundaries.
    :type       x_bounds:                  str
    :param      y_bounds:                  The y-axis boundaries.
    :type       y_bounds:                  str
    :param      clad_index:                The fused clad refractive index.
    :type       clad_index:                str
    :param      core_position_scrambling:  The core position scrambling.
    :type       core_position_scrambling:  float
    :param      resolution:                The rasterisation resolution for the geometry.
    :type       resolution:                int
    :param      background_index:          The background refractive index.
    :type       background_index:          float

    :returns:   The simulation geometry
    :rtype:     Geometry
    """

    # assert (x_bounds in ['left', 'right', '']) or is, f"Invalid 'x_bounds' input: {x_bounds}, value has to be in ['left', 'rigth']."
    # assert y_bounds in ['top', 'bottom', ''], f"Invalid 'y_bounds' input: {y_bounds}, value has to be in ['top', 'bottom']."

    if clad_index.lower() == 'silica':
        index = fiber_catalogue.get_silica_index(wavelength=wavelength)

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
    #  Geometry arguments --------------------------
    wavelength: float
    """ Wavelenght at which evaluate the computation """
    clad_rotation: float = 0
    """ Rotation of the clad structure [degree] """
    capillary_tube: object = None
    """ Additional optical structure such as clad to add """
    resolution: int = 100
    """ Discretization of the mesh [resolution x resolution] """
    clad_structure: object = None
    """ Initial optical structure """
    fiber_list: list = tuple()
    """ List of the fiber to add to the optical structure """
    fiber_radius: float = 62.5e-6
    """ Fiber radius for the clad fused structure """
    fusion_degree: float = None
    """ Fusion degree for the clad fused structure """
    x_bounds: str = ''
    """ X-boundaries """
    y_bounds: str = ''
    """ Y-boundaries """
    air_padding_factor: float = 1.2
    """ Padding factor for air around the optica structure, preferable over 1.2 """
    gaussian_filter_factor: float = None
    """ Gaussian blurring of the optical structure """

    #  Solver arguments --------------------------
    n_sorted_mode: int = 4
    """ Number of mode that are computed """
    n_added_mode: int = 4
    """ Number of mode that are computed additionally to the sorted modes """
    itr_final: float = 0.05
    """ Final ITR at which evaluate the modes """
    itr_initial: float = 1.0
    """ Start ITR at which evaluate the modes """
    n_step: int = 500
    """ Discretization of the z-profile """
    extrapolation_order: int = 2
    """ Eigen_value extrapolation for slice solving """
    core_position_scrambling: float = 0
    """ Scrambling of the clad core position """
    index_scrambling: float = 0
    """ Scrambling of the structure refractive index """
    boundaries: list = (Boundaries2D(),)
    """ List of boundaries cndition to which evaluate to modes """
    accuracy: int = 2
    """ Accuracy of the finit-difference set of value """

    #  Plot arguments --------------------------
    plot_geometry: bool = False
    """ Plot the computed geometry mesh prior computation """
    plot_cladding: bool = False
    """ Plot the cladding structure prior computation """
    plot_field: bool = False
    """ Plot the mode field after computation """
    plot_adiabatic: bool = False
    """ Plot the adiabatic criterion after computation """
    plot_coupling: bool = False
    """ Plot the mode coupling after computation """
    plot_beating_length: bool = False
    """ Plot the mode beating length after computation """
    plot_eigen_values: bool = False
    """ Plot the computed eigen_values after computation """
    plot_index: bool = False
    """ Plot the computed effective index after computation """
    plot_beta: bool = False
    """ Plot the computed propagation constant after computation """

    #  Extra arguments --------------------------
    debug_mode: int = 1
    """ Level of debug mode printing [0, 1, 2, 3]"""
    auto_label: bool = False
    """ Enable auto labeling of the supermodes """
    generate_report: bool = False
    """ Generate final pdf reports containing geometry, fields, coupling, adiabatic criterions """
    save_superset: bool = False
    """ Save the created superset instance into a pickle file for further use """

    def __post_init__(self):
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

        self.solver.superset.sorting_modes(sorting_method='beta')

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
        """
        Saves a superset instance in the form of a serialized files using the picles library.

        :param      filename:   The filename
        :type       filename:   str
        :param      directory:  The directory
        :type       directory:  str

        :returns:   The path directory of the saved instance
        :rtype:     Path
        """
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

        self.solver.superset.generate_pdf_report(
            filename=filename,
            **kwargs
        )

        return filename

    def plot(self, *args, **kwargs):
        return self.solver.superset.plot(*args, **kwargs)


# -
