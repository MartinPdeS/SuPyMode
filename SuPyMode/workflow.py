
from dataclasses import dataclass

from FiberFusing import Geometry, BackGround
from SuPyMode.solver import SuPySolver
from FiberFusing.fiber import catalogue as fiber_catalogue
from SuPyMode.profiles import AlphaProfile
from FiberFusing import configuration

from PyFinitDiff.boundaries import Boundaries2D


@dataclass
class Workflow():
    fiber_list: list
    """ List of the fiber to add to the optical structure """
    wavelength: float
    """ Wavelenght at which evaluate the computation """
    clad_structure: object = None
    """ Initial optical structure """
    clad_rotation: float = 0
    """ Rotation of the clad structure [degree] """
    capillary_tube: object = None
    """ Additional optical structure such as clad to add """
    resolution: int = 100
    """ Discretization of the mesh [resolution x resolution] """
    n_sorted_mode: int = 4
    """ Number of mode that are computed """
    n_added_mode: int = 4
    """ Number of mode that are computed additionally to the sorted modes """
    gaussian_filter_factor: float = None
    """ Gaussian blurring of the optical structure """
    fiber_radius: float = 62.5e-6
    """ Fiber radius for the clad fused structure """
    fusion_degree: float = None
    """ Fusion degree for the clad fused structure """
    air_padding_factor: float = 1.2
    """ Padding factor for air around the optica structure, preferable over 1.2 """
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
    debug_mode: bool = False
    """ Enable debug mode printing """
    auto_label: bool = False
    """ Enable auto labeling of the supermodes """
    boundaries: list = (Boundaries2D(),)
    """ List of boundaries cndition to which evaluate to modes """
    accuracy: int = 2
    """ Accuracy of the finit-difference set of value """
    x_bounds: str = 'centering'
    """ X-boundaries """
    y_bounds: str = 'centering'
    """ Y-boundaries """
    generate_report: bool = False
    """ Generate final pdf reports containing geometry, fields, coupling, adiabatic criterions """
    save_superset: bool = False
    """ Save the created superset instance into a pickle file for further use """

    def __post_init__(self):
        core_positions = self._initialize_cladding_structure_()

        self._initialize_fiber_structures_(core_positions=core_positions)

        self._initialize_geometry_()

        self._initialize_solver_()

        if self.plot_field:
            self.plot(plot_type='field')

        if self.plot_adiabatic:
            self.plot(plot_type='adiabatic')

        if self.plot_coupling:
            self.plot(plot_type='normalized-coupling')

        if self.generate_report:
            self.generate_pdf_report()

        if self.save_superset:
            self.save_superset_instance()

    @property
    def superset(self):
        return self.solver.superset

    def _initialize_cladding_structure_(self) -> list:
        """
        Initializes the cladding structure to be added to mesh geometry.
        Returns a list of the core positions

        :returns:   list of core positions
        :rtype:     list
        """
        if self.clad_structure is None:
            return [(0, 0) for _ in self.fiber_list]

        silica_index = fiber_catalogue.get_silica_index(wavelength=self.wavelength)

        kwargs = dict(
            fiber_radius=self.fiber_radius,
            index=silica_index,
            core_position_scrambling=self.core_position_scrambling
        )

        if self.fusion_degree is None:
            self.clad_structure = self.clad_structure(**kwargs)
        else:
            self.clad_structure = self.clad_structure(**kwargs, fusion_degree=self.fusion_degree)

        self.clad_structure.rotate(angle=self.clad_rotation)

        core_positions = self.clad_structure.cores

        return core_positions

    def _initialize_fiber_structures_(self, core_positions: list = None) -> None:
        """
        Initializes the fiber structures to be added to the mesh geometry.

        :returns:   No returns
        :rtype:     None
        """
        new_fiber_list = []
        for n, (position, fiber) in enumerate(zip(core_positions, self.fiber_list)):
            fiber.set_position(position=position)
            new_fiber_list.append(fiber)

        self.fiber_list = new_fiber_list

        if self.plot_cladding:
            self.clad_structure.plot().show()

    def _initialize_geometry_(self) -> None:
        """
        Initializes the mesh geometry.

        :returns:   No returns
        :rtype:     None
        """
        self.geometry = Geometry(
            background=BackGround(index=1),
            x_bounds=self.x_bounds,
            y_bounds=self.y_bounds,
            resolution=self.resolution,
            index_scrambling=self.index_scrambling,
            gaussian_filter=self.gaussian_filter_factor,
            boundary_pad_factor=self.air_padding_factor
        )

        if self.capillary_tube is not None:
            self.geometry.add_structure(self.capillary_tube)

        if self.clad_structure is not None:
            self.geometry.add_structure(self.clad_structure)

        self.geometry.add_fiber(*self.fiber_list)

        if self.plot_geometry:
            figure = self.geometry.plot()
            figure.show_colorbar = True
            figure.show()

    def _initialize_solver_(self) -> None:

        self.solver = SuPySolver(
            geometry=self.geometry,
            tolerance=1e-20,
            max_iter=5000,
            show_iteration=self.debug_mode,
            accuracy=self.accuracy,
            show_eigenvalues=False,
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

    def save_superset_instance(self, filename: str = 'auto', directory: str = 'auto') -> None:
        """
        Saves a superset instance in the form of a serialized files using the picles library.

        :param      filename:   The filename
        :type       filename:   str
        :param      directory:  The directory
        :type       directory:  str

        :returns:   No returns
        :rtype:     None
        """
        if filename == 'auto':
            filename = self._get_auto_generated_filename_()

        self.solver.superset.save_instance(
            filename=filename,
            directory=directory
        )

    def generate_pdf_report(self, filename: str = 'auto', **kwargs) -> None:
        """
        Generate a pdf file compiling the essential computed components.

        :param      filename:  The filename
        :type       filename:  str
        :param      kwargs:    The keywords arguments
        :type       kwargs:    dictionary

        :returns:   No returns
        :rtype:     None
        """
        if filename == 'auto':
            filename = self._get_auto_generated_filename_()

        self.solver.superset.generate_pdf_report(
            filename=filename + ".pdf",
            **kwargs
        )

    def plot(self, *args, **kwargs):
        return self.solver.superset.plot(*args, **kwargs)


# -
