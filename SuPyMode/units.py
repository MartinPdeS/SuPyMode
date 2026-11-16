"""Shared TypedUnit registry and validators.

Legacy bare dimensional values are interpreted in SI units (workflow rotation
angles use degrees). Explicit quantities are validated before numerical
conversion. Refractive indices and taper ratios require plain numbers.
"""

from TypedUnit import (
    Angle,
    Finite,
    Length,
    Matrix,
    NoUnit,
    Positive,
    Quantity,
    Scalar,
    Vector,
    unit_validator,
    ureg,
    validate_annotation,
    validate_quantity,
)

__all__ = [
    "Angle",
    "Finite",
    "Length",
    "Matrix",
    "NoUnit",
    "Positive",
    "Quantity",
    "Scalar",
    "Vector",
    "unit_validator",
    "ureg",
    "validate_annotation",
    "validate_quantity",
]
