# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| SuPyMode/profiles.py                            |       79 |        5 |        6 |        1 |     92.94% |128-138, 190->189 |
| SuPyMode/propagation.py                         |       64 |       55 |       12 |        0 |     11.84% |28-34, 58-101, 107-114, 146-199 |
| SuPyMode/representation/adiabatic.py            |       29 |        2 |        6 |        2 |     88.57% |    54, 96 |
| SuPyMode/representation/beating\_length.py      |       23 |        2 |        4 |        2 |     85.19% |    74, 80 |
| SuPyMode/representation/beta.py                 |       17 |        1 |        2 |        1 |     89.47% |        59 |
| SuPyMode/representation/eigen\_value.py         |       16 |        1 |        2 |        1 |     88.89% |        58 |
| SuPyMode/representation/field.py                |      134 |       50 |       50 |       18 |     55.43% |56, 76, 78->85, 89-90, 96-97, 103-104, 146, 170, 204->207, 229-262, 283-287, 318, 321, 327-328, 333, 336, 340, 343, 375-393, 445-449, 468->471, 471->478, 478->481, 481->484 |
| SuPyMode/representation/index.py                |       16 |        1 |        2 |        1 |     88.89% |        58 |
| SuPyMode/representation/normalized\_coupling.py |       26 |        2 |        6 |        2 |     87.50% |    52, 85 |
| SuPyMode/special.py                             |       95 |       95 |        8 |        0 |      0.00% |     3-251 |
| SuPyMode/superset.py                            |      203 |       74 |       50 |        3 |     59.29% |60, 84, 96-99, 115-117, 160, 180-185, 196-205, 232-241, 263-282, 298-302, 322-338, 376-418, 447-460, 515-516, 601-602, 625-628 |
| SuPyMode/superset\_plots.py                     |      104 |       24 |       50 |        8 |     71.43% |53, 97, 141, 171-182, 258-259, 366, 368, 370, 373-378, 450-463 |
| SuPyMode/tools/special.py                       |       71 |       71 |        8 |        0 |      0.00% |     4-158 |
| SuPyMode/tools/utils.py                         |       96 |       96 |       42 |        0 |      0.00% |     4-271 |
| SuPyMode/utils.py                               |      112 |       43 |       44 |       11 |     56.41% |22-29, 63-71, 77-82, 86-91, 98-110, 115, 134-140, 163, 192, 208->213, 242, 245-246, 251, 258, 291, 296, 306-308 |
| SuPyMode/workflow.py                            |       77 |        8 |       16 |        4 |     84.95% |143, 160-165, 210, 215-216, 234, 246 |
| **TOTAL**                                       | **1193** |  **530** |  **310** |   **54** | **51.96%** |           |

1 file skipped due to complete coverage.


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