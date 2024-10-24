# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ------: | --------: |
| SuPyMode/helper.py                              |       48 |        3 |       16 |        5 |     88% |68, 77->80, 81, 151->155, 168 |
| SuPyMode/mode\_label.py                         |       46 |        4 |       16 |        3 |     89% |30->37, 110, 146, 159, 170 |
| SuPyMode/profiles.py                            |      251 |       37 |       28 |       10 |     83% |53, 68, 136, 159-166, 180, 232-239, 251-257, 301, 322, 325, 328, 456, 479->482, 489-493, 505-506, 512-513, 519-520, 526-527, 533-534, 594, 695-707, 736->739, 739->742, 742->745 |
| SuPyMode/propagation.py                         |       65 |       55 |       12 |        0 |     13% |20-26, 39-62, 66-73, 96-141 |
| SuPyMode/representation/adiabatic.py            |       17 |        2 |        4 |        2 |     81% |     7, 56 |
| SuPyMode/representation/base.py                 |       76 |       15 |       14 |        2 |     79% |4, 48, 66, 123-124, 140, 152, 164, 176, 212, 224, 252, 268, 279, 290 |
| SuPyMode/representation/beating\_length.py      |       15 |        3 |        2 |        1 |     76% |  7, 62-63 |
| SuPyMode/representation/beta.py                 |       13 |        1 |        2 |        1 |     87% |         7 |
| SuPyMode/representation/eigen\_value.py         |       13 |        1 |        2 |        1 |     87% |         7 |
| SuPyMode/representation/field.py                |      112 |       35 |       42 |       15 |     61% |7, 59, 83, 117->120, 140-157, 178-182, 211, 214, 220-221, 226, 229, 233, 236, 289->exit, 337-341, 358->361, 361->365, 365->368, 368->371 |
| SuPyMode/representation/index.py                |       12 |        1 |        2 |        1 |     86% |         7 |
| SuPyMode/representation/normalized\_coupling.py |       18 |        2 |        4 |        2 |     82% |     7, 55 |
| SuPyMode/solver.py                              |       63 |        2 |        6 |        1 |     96% |     57-58 |
| SuPyMode/supermode.py                           |      119 |       19 |       32 |        9 |     81% |164, 177-180, 194, 229, 254, 295, 297->301, 302-303, 306-307, 310-311, 344, 378, 383-384 |
| SuPyMode/superset.py                            |      217 |       71 |       62 |        9 |     63% |58, 94, 106-109, 170, 188-193, 204-215, 242-251, 269-287, 303-307, 325-341, 377-409, 436-447, 498-499, 578->581, 602-605, 673->675, 675->677, 678, 681->683, 683->685, 685->exit |
| SuPyMode/superset\_plots.py                     |       94 |       17 |       42 |        4 |     77% |41, 152-155, 221, 332, 337-342, 413-427 |
| SuPyMode/utils.py                               |      116 |       45 |       46 |       12 |     56% |8-9, 27-34, 65-73, 77-82, 86-91, 96-108, 113, 130-136, 159, 187, 199->204, 230, 233-234, 237, 242, 271, 274, 284-286 |
| SuPyMode/workflow.py                            |      128 |       19 |       40 |       14 |     79% |75-78, 103, 152, 278, 280, 318, 320, 322, 324, 326, 328, 330, 335, 337, 341, 354-359 |
|                                       **TOTAL** | **1423** |  **332** |  **372** |   **92** | **72%** |           |


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