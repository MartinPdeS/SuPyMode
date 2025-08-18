# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ------: | --------: |
| SuPyMode/helper.py                              |       64 |       19 |       24 |        7 |     61% |28-86, 150, 159->162, 163, 233->237, 246->249, 250, 253 |
| SuPyMode/mode\_label.py                         |       42 |        4 |       16 |        3 |     88% |23->30, 103, 139, 152, 163 |
| SuPyMode/profiles.py                            |      256 |       37 |       26 |        8 |     84% |53, 68, 136, 159-166, 180, 232-239, 251-257, 301, 322, 325, 328, 456, 479->482, 489-493, 505-506, 512-513, 519-520, 526-527, 533-534, 594, 695-707, 755->754 |
| SuPyMode/propagation.py                         |       64 |       55 |       12 |        0 |     12% |21-27, 45-68, 72-79, 109-154 |
| SuPyMode/representation/adiabatic.py            |       19 |        1 |        2 |        1 |     90% |        57 |
| SuPyMode/representation/base.py                 |       74 |       17 |       12 |        3 |     74% |48, 66, 94-95, 106, 123-124, 140, 152, 164, 176, 212, 224, 252, 268, 279, 290 |
| SuPyMode/representation/beating\_length.py      |       13 |        0 |        0 |        0 |    100% |           |
| SuPyMode/representation/beta.py                 |       11 |        0 |        0 |        0 |    100% |           |
| SuPyMode/representation/eigen\_value.py         |       11 |        0 |        0 |        0 |    100% |           |
| SuPyMode/representation/field.py                |      111 |       42 |       40 |       13 |     54% |59, 83, 117->120, 140-157, 178-182, 211, 214, 220-221, 226, 229, 233, 236, 270-290, 337-341, 358->361, 361->365, 365->368, 368->371 |
| SuPyMode/representation/index.py                |       10 |        0 |        0 |        0 |    100% |           |
| SuPyMode/representation/normalized\_coupling.py |       16 |        1 |        2 |        1 |     89% |        55 |
| SuPyMode/solver.py                              |       58 |        1 |        6 |        1 |     97% |        62 |
| SuPyMode/supermode.py                           |      109 |       16 |       18 |        7 |     82% |215, 228-231, 245, 280, 305, 346, 348->352, 353-354, 357-358, 361-362, 395 |
| SuPyMode/superset.py                            |      217 |       71 |       62 |        9 |     63% |58, 94, 106-109, 170, 188-193, 204-215, 242-251, 269-287, 303-307, 325-341, 377-409, 436-447, 498-499, 578->581, 602-605, 673->675, 675->677, 678, 681->683, 683->685, 685->exit |
| SuPyMode/superset\_plots.py                     |       93 |       19 |       42 |        6 |     74% |41, 152-155, 221, 325, 327, 329, 332-337, 408-422 |
| SuPyMode/utils.py                               |      113 |       43 |       44 |       11 |     57% |27-34, 65-73, 77-82, 86-91, 96-108, 113, 130-136, 159, 187, 199->204, 230, 233-234, 237, 242, 271, 274, 284-286 |
| SuPyMode/workflow.py                            |      122 |       18 |       36 |       14 |     78% |157, 159, 197, 199, 201, 203, 205, 207, 209, 214, 216, 220, 237-242, 282, 288-289, 318 |
|                                       **TOTAL** | **1403** |  **344** |  **342** |   **84** | **71%** |           |


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