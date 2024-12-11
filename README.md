# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ------: | --------: |
| SuPyMode/helper.py                              |       64 |        5 |       24 |        9 |     84% |68, 75->78, 79, 81->84, 150, 159->162, 163, 233->237, 250 |
| SuPyMode/mode\_label.py                         |       46 |        4 |       16 |        3 |     89% |30->37, 110, 146, 159, 170 |
| SuPyMode/profiles.py                            |      251 |       37 |       28 |       10 |     83% |54, 69, 137, 160-167, 181, 233-240, 252-258, 302, 323, 326, 329, 457, 480->483, 490-494, 506-507, 513-514, 520-521, 527-528, 534-535, 595, 696-708, 738->741, 741->744, 744->exit |
| SuPyMode/propagation.py                         |       65 |       55 |       12 |        0 |     13% |20-26, 39-62, 66-73, 96-141 |
| SuPyMode/representation/adiabatic.py            |       21 |        2 |        4 |        2 |     84% |     7, 57 |
| SuPyMode/representation/base.py                 |       76 |       15 |       14 |        2 |     79% |4, 48, 66, 123-124, 140, 152, 164, 176, 212, 224, 252, 268, 279, 290 |
| SuPyMode/representation/beating\_length.py      |       15 |        3 |        2 |        1 |     76% |  7, 62-63 |
| SuPyMode/representation/beta.py                 |       13 |        1 |        2 |        1 |     87% |         7 |
| SuPyMode/representation/eigen\_value.py         |       13 |        1 |        2 |        1 |     87% |         7 |
| SuPyMode/representation/field.py                |      112 |       35 |       42 |       15 |     61% |7, 59, 83, 117->120, 140-157, 178-182, 211, 214, 220-221, 226, 229, 233, 236, 289->exit, 337-341, 358->361, 361->365, 365->368, 368->371 |
| SuPyMode/representation/index.py                |       12 |        1 |        2 |        1 |     86% |         7 |
| SuPyMode/representation/normalized\_coupling.py |       18 |        2 |        4 |        2 |     82% |     7, 55 |
| SuPyMode/solver.py                              |       59 |        0 |        4 |        0 |    100% |           |
| SuPyMode/supermode.py                           |      119 |       19 |       32 |        9 |     81% |164, 177-180, 194, 229, 254, 295, 297->301, 302-303, 306-307, 310-311, 344, 378, 383-384 |
| SuPyMode/superset.py                            |      217 |       71 |       62 |        9 |     63% |58, 94, 106-109, 170, 188-193, 204-215, 242-251, 269-287, 303-307, 325-341, 377-409, 436-447, 498-499, 578->581, 602-605, 673->675, 675->677, 678, 681->683, 683->685, 685->exit |
| SuPyMode/superset\_plots.py                     |       93 |       17 |       42 |        4 |     77% |41, 152-155, 221, 327, 332-337, 408-422 |
| SuPyMode/utils.py                               |      116 |       45 |       46 |       12 |     56% |8-9, 27-34, 65-73, 77-82, 86-91, 96-108, 113, 130-136, 159, 187, 199->204, 230, 233-234, 237, 242, 271, 274, 284-286 |
| SuPyMode/workflow.py                            |      130 |       19 |       40 |       14 |     79% |75-78, 94, 154, 280, 282, 320, 322, 324, 326, 328, 330, 332, 337, 339, 343, 356-361 |
|                                       **TOTAL** | **1440** |  **332** |  **378** |   **95** | **72%** |           |


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