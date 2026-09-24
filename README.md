# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| SuPyMode/plotter.py                             |      121 |       73 |       10 |        0 |     38.17% |107-142, 186-203, 209-224, 227-240, 243-259, 265-269, 272-281, 284-295, 301-321 |
| SuPyMode/propagation.py                         |       64 |       55 |       12 |        0 |     11.84% |28-34, 58-101, 107-114, 146-199 |
| SuPyMode/representation/adiabatic.py            |       29 |        2 |        6 |        2 |     88.57% |    54, 96 |
| SuPyMode/representation/beating\_length.py      |       23 |        2 |        4 |        2 |     85.19% |    74, 80 |
| SuPyMode/representation/field.py                |      133 |       42 |       50 |       19 |     61.20% |55, 75, 77-\>84, 88-89, 95-96, 102-103, 145, 169, 203-\>206, 228-261, 282-286, 317, 320, 326-327, 332, 335, 339, 342, 391-\>exit, 444-448, 467-\>470, 470-\>477, 477-\>480, 480-\>483 |
| SuPyMode/representation/normalized\_coupling.py |       26 |        2 |        6 |        2 |     87.50% |    52, 85 |
| SuPyMode/superset.py                            |      201 |       73 |       50 |        3 |     59.36% |78, 90-93, 109-111, 154, 174-179, 190-199, 226-235, 257-276, 292-296, 316-332, 367-410, 439-452, 507-508, 593-594, 617-620 |
| SuPyMode/superset\_plots.py                     |      127 |       23 |       56 |       11 |     77.05% |29, 60, 63-64, 69, 76, 202-213, 289-290, 417, 425, 427, 429, 432-437 |
| SuPyMode/utils.py                               |       67 |       18 |       22 |        6 |     66.29% |40, 56-\>61, 91-99, 116-122, 145, 176, 181, 191-193 |
| SuPyMode/workflow.py                            |       77 |        8 |       16 |        4 |     84.95% |146, 165-167, 212, 217-218, 236, 250 |
| **TOTAL**                                       |  **947** |  **298** |  **240** |   **49** | **64.70%** |           |

4 files skipped due to complete coverage.


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/MartinPdeS/SuPyMode/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/MartinPdeS/SuPyMode/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2FMartinPdeS%2FSuPyMode%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.