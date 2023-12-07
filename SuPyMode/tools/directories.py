#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import SuPyMode


__all__ = [
    'root_path',
    'project_path',
    'test_path',
    'instance_directory',
    'version_path',
    'validation_data_path',
    'doc_path',
    'static_doc_path',
    'logo_path',
    'doc_css_path',
]

root_path = Path(SuPyMode.__path__[0])

project_path = root_path.parents[0]

test_path = project_path.joinpath('tests')

instance_directory = root_path.joinpath('superset_instances')

reports_path = root_path.joinpath('reports')

version_path = root_path.joinpath('VERSION')

validation_data_path = root_path.joinpath('validation_data')

doc_path = root_path.parents[0].joinpath('docs')

static_doc_path = doc_path.joinpath('images')

logo_path = static_doc_path.joinpath('logo.png')

doc_css_path = doc_path.joinpath('source/_static/default.css')

rtd_example = 'https://pymiesim.readthedocs.io/en/latest/examples.html'

rtd_material = 'https://pymiesim.readthedocs.io/en/latest/material.html'

rtd_lp_mode = 'https://pymiesim.readthedocs.io/en/latest/LPModes.html'


if __name__ == '__main__':
    for path_name in __all__:
        path = locals()[path_name]
        print(path)
        assert path.exists(), f"Path {path_name} do not exists"

# -
