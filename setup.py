# https://stackoverflow.com/questions/42585210/extending-setuptools-extension-to-use-cmake-in-setup-py

import os
import sys
import pathlib

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig
from setuptools.command.build_py import build_py as build_py_orig


major, minor, micro = sys.version_info[:3]


class build_py(build_py_orig):
    def run(self):
        self.run_command("build_ext")
        return super().run()


class CMakeExtension(Extension):

    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):
    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        root_directory = pathlib.Path().absolute()

        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))

        extdir.mkdir(parents=True, exist_ok=True)

        cmake_args = [f'-DPYBIND11_PYTHON_VERSION={major}.{minor}']

        os.chdir(str(build_temp))

        command = ['cmake', str(root_directory)] + cmake_args

        self.spawn(command)

        if not self.dry_run:
            build_command = ['cmake', '--build', '.', '-j4']
            self.spawn(build_command)

        os.chdir(str(root_directory))


setup(
    ext_modules=[CMakeExtension(name='SuPyMode/binary')],
    cmdclass={"build_ext": build_ext, 'build_py': build_py},
)

# -
