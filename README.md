# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| SuPyMode/plotter.py                             |      121 |       73 |       10 |        0 |     38.17% |107-142, 186-203, 209-224, 227-240, 243-259, 265-269, 272-281, 284-295, 301-321 |
| SuPyMode/propagation.py                         |       64 |       55 |       12 |        0 |     11.84% |28-34, 58-101, 107-114, 146-199 |
| SuPyMode/representation/adiabatic.py            |       29 |        2 |        6 |        2 |     88.57% |    54, 96 |
| SuPyMode/representation/beating\_length.py      |       23 |        2 |        4 |        2 |     85.19% |    74, 80 |
| SuPyMode/representation/beta.py                 |       17 |        1 |        2 |        1 |     89.47% |        59 |
| SuPyMode/representation/eigen\_value.py         |       16 |        1 |        2 |        1 |     88.89% |        58 |
| SuPyMode/representation/field.py                |      134 |       50 |       50 |       18 |     55.43% |56, 76, 78->85, 89-90, 96-97, 103-104, 146, 170, 204->207, 229-262, 283-287, 318, 321, 327-328, 333, 336, 340, 343, 375-393, 445-449, 468->471, 471->478, 478->481, 481->484 |
| SuPyMode/representation/index.py                |       16 |        1 |        2 |        1 |     88.89% |        58 |
| SuPyMode/representation/normalized\_coupling.py |       26 |        2 |        6 |        2 |     87.50% |    52, 85 |
| SuPyMode/superset.py                            |      201 |       73 |       50 |        3 |     59.36% |78, 90-93, 109-111, 154, 174-179, 190-199, 226-235, 257-276, 292-296, 316-332, 367-410, 439-452, 507-508, 593-594, 617-620 |
| SuPyMode/superset\_plots.py                     |      126 |       23 |       56 |       11 |     76.92% |24, 55, 58-59, 64, 71, 197-208, 284-285, 412, 420, 422, 424, 427-432 |
| SuPyMode/utils.py                               |       67 |       18 |       22 |        6 |     66.29% |40, 56->61, 91-99, 116-122, 145, 176, 181, 191-193 |
| SuPyMode/workflow.py                            |       77 |        8 |       16 |        4 |     84.95% |146, 165-167, 212, 217-218, 236, 250 |
| **TOTAL**                                       |  **948** |  **309** |  **240** |   **51** | **63.30%** |           |

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