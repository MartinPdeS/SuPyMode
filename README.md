# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| SuPyMode/mode\_label.py                         |       42 |        4 |       16 |        3 |     87.93% |23->30, 103, 139, 152, 163 |
| SuPyMode/profiles.py                            |      249 |       36 |       26 |        8 |     84.00% |52, 67, 135, 158-165, 179, 231-238, 250-256, 300, 321, 324, 327, 455, 478->481, 488-492, 504-505, 511-512, 518-519, 525-526, 532-533, 593, 690-700, 747->746 |
| SuPyMode/propagation.py                         |       64 |       55 |       12 |        0 |     11.84% |21-27, 45-68, 72-79, 109-154 |
| SuPyMode/representation/adiabatic.py            |       19 |        1 |        2 |        1 |     90.48% |        57 |
| SuPyMode/representation/base.py                 |       74 |       17 |       12 |        3 |     74.42% |48, 66, 94-95, 106, 123-124, 140, 152, 164, 176, 212, 224, 252, 268, 279, 290 |
| SuPyMode/representation/field.py                |      111 |       42 |       40 |       13 |     54.30% |59, 83, 117->120, 140-157, 178-182, 211, 214, 220-221, 226, 229, 233, 236, 270-290, 337-341, 358->361, 361->365, 365->368, 368->371 |
| SuPyMode/representation/normalized\_coupling.py |       16 |        1 |        2 |        1 |     88.89% |        55 |
| SuPyMode/solver.py                              |       58 |        1 |        6 |        1 |     96.88% |        62 |
| SuPyMode/special.py                             |       71 |       71 |        8 |        0 |      0.00% |     4-158 |
| SuPyMode/supermode.py                           |      104 |       16 |       18 |        7 |     81.15% |215, 228-231, 245, 280, 305, 346, 348->352, 353-354, 357-358, 361-362, 395 |
| SuPyMode/superset.py                            |      205 |       72 |       50 |        3 |     60.39% |63, 99, 111-114, 175, 195-200, 211-220, 247-256, 278-297, 313-317, 337-353, 391-433, 462-475, 528-529, 614-615, 638-641 |
| SuPyMode/superset\_plots.py                     |      104 |       24 |       50 |        8 |     71.43% |59, 103, 147, 177-188, 264-265, 372, 374, 376, 379-384, 456-469 |
| SuPyMode/tools/special.py                       |       71 |       71 |        8 |        0 |      0.00% |     4-158 |
| SuPyMode/tools/utils.py                         |       96 |       96 |       42 |        0 |      0.00% |     4-253 |
| SuPyMode/utils.py                               |      113 |       43 |       44 |       11 |     56.69% |27-34, 65-73, 77-82, 86-91, 96-108, 113, 130-136, 159, 187, 199->204, 230, 233-234, 237, 242, 271, 274, 284-286 |
| SuPyMode/workflow.py                            |       77 |        8 |       16 |        4 |     84.95% |143, 160-165, 210, 215-216, 234, 246 |
| **TOTAL**                                       | **1519** |  **558** |  **352** |   **63** | **59.33%** |           |

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