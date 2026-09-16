"""Clean-room exact checker for the public radial-stress tensor/force seam.

This independently reimplements only the short R31--R34 formulas exposed by the
public corrected Navier--Stokes reader.  It does not import source-bundle code.
"""

from __future__ import annotations
from fractions import Fraction
from typing import Sequence, Tuple

Vec3 = Tuple[Fraction, Fraction, Fraction]
MatDerivs = Tuple[Vec3, Vec3, Vec3]

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_REPO_COMMIT = "fab69fdc4ac197159b8e6ae8d73a82bde2b20d55"
RESEARCH_STATE_SOURCE_COMMIT = "e0a04c078eaf0209e40b3d47d6f6dffb3f2a3e7f"
SOURCE_COMPONENT = (
    "proof_sources/lanes/web_bundle_audit/tex_handoff/portable_upstream/"
    "radial_stress_reconstruction_audit.tex"
)
SOURCE_COMPONENT_SHA256 = (
    "dd6c119670316f1804226723caaefe062a580e40d2101efab86d823b5858ed8a"
)
BUNDLE_SHA256 = "43b128e24f395327b2dd0f9874ca1ff120a52d625ab454f7df97d51f10c328c5"
DOI = "10.5281/zenodo.22678406"

PAPER_EXACT_VELOCITY_AVAILABLE = False
FULL_RECONSTRUCTION = False


def _q(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, int)):
        raise TypeError(f"{name} must be an exact int or Fraction")
    return Fraction(value)


def _vadd(*vectors: Sequence[Fraction]) -> Vec3:
    return tuple(sum((v[i] for v in vectors), Fraction(0)) for i in range(3))  # type: ignore[return-value]


def _vscale(c: Fraction, v: Sequence[Fraction]) -> Vec3:
    return tuple(c * v[i] for i in range(3))  # type: ignore[return-value]


def cylindrical_basis(
    x: Fraction | int, y: Fraction | int, r: Fraction | int
) -> tuple[Vec3, Vec3, Vec3, MatDerivs, MatDerivs]:
    """Return (e_r,e_theta,e_z,d e_r,d e_theta) at an exact point."""
    x = _q(x, "x")
    y = _q(y, "y")
    r = _q(r, "r")
    if r <= 0:
        raise ValueError("r must be positive")
    if x * x + y * y != r * r:
        raise ValueError("x^2+y^2 must equal r^2 exactly")
    er = (x / r, y / r, Fraction(0))
    et = (-y / r, x / r, Fraction(0))
    ez = (Fraction(0), Fraction(0), Fraction(1))
    r3 = r * r * r
    der = (
        (Fraction(1, 1) / r - x * x / r3, -x * y / r3, Fraction(0)),
        (-x * y / r3, Fraction(1, 1) / r - y * y / r3, Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    det = (
        (x * y / r3, Fraction(1, 1) / r - x * x / r3, Fraction(0)),
        (-Fraction(1, 1) / r + y * y / r3, -x * y / r3, Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    return er, et, ez, der, det


def divergence_scalar_outer(
    sigma: Fraction | int,
    grad_sigma: Sequence[Fraction | int],
    a: Sequence[Fraction | int],
    da: Sequence[Sequence[Fraction | int]],
    b: Sequence[Fraction | int],
    db: Sequence[Sequence[Fraction | int]],
) -> Vec3:
    """Compute div(sigma * a tensor b) directly by Cartesian product rule."""
    sigma = _q(sigma, "sigma")
    gs = tuple(_q(v, "grad_sigma") for v in grad_sigma)
    av = tuple(_q(v, "a") for v in a)
    bv = tuple(_q(v, "b") for v in b)
    dav = tuple(tuple(_q(v, "da") for v in row) for row in da)
    dbv = tuple(tuple(_q(v, "db") for v in row) for row in db)
    if (
        len(gs) != 3
        or len(av) != 3
        or len(bv) != 3
        or len(dav) != 3
        or len(dbv) != 3
    ):
        raise ValueError("all vectors/derivative tables must be three-dimensional")
    out = [Fraction(0), Fraction(0), Fraction(0)]
    for i in range(3):
        total = Fraction(0)
        for j in range(3):
            total += gs[j] * av[i] * bv[j]
            total += sigma * dav[j][i] * bv[j]
            total += sigma * av[i] * dbv[j][j]
        out[i] = total
    return tuple(out)  # type: ignore[return-value]


def cylindrical_force_components(
    r: Fraction | int,
    sigma1: Fraction | int,
    dr_sigma1: Fraction | int,
    dz_sigma1: Fraction | int,
    sigma2: Fraction | int,
    dr_sigma2: Fraction | int,
) -> Vec3:
    """R34 cylindrical components (radial, azimuthal, axial)."""
    r = _q(r, "r")
    if r <= 0:
        raise ValueError("r must be positive")
    sigma1 = _q(sigma1, "sigma1")
    dr_sigma1 = _q(dr_sigma1, "dr_sigma1")
    dz_sigma1 = _q(dz_sigma1, "dz_sigma1")
    sigma2 = _q(sigma2, "sigma2")
    dr_sigma2 = _q(dr_sigma2, "dr_sigma2")
    return (
        dz_sigma1,
        dr_sigma2 + 2 * sigma2 / r,
        dr_sigma1 + sigma1 / r,
    )


def radial_inverse_rhs(
    source: Fraction | int, moment_term: Fraction | int
) -> Fraction:
    """Public R9 right-hand side -F + b_e M_e F."""
    return -_q(source, "source") + _q(moment_term, "moment_term")


def cartesian_tensor_divergence(
    x: Fraction | int,
    y: Fraction | int,
    r: Fraction | int,
    sigma1: Fraction | int,
    dr_sigma1: Fraction | int,
    dz_sigma1: Fraction | int,
    sigma2: Fraction | int,
    dr_sigma2: Fraction | int,
) -> Vec3:
    """Direct Cartesian divergence of the public R32 symmetric tensor."""
    x = _q(x, "x")
    y = _q(y, "y")
    r = _q(r, "r")
    s1 = _q(sigma1, "sigma1")
    s1r = _q(dr_sigma1, "dr_sigma1")
    s1z = _q(dz_sigma1, "dz_sigma1")
    s2 = _q(sigma2, "sigma2")
    s2r = _q(dr_sigma2, "dr_sigma2")
    er, et, ez, der, det = cylindrical_basis(x, y, r)
    zeros: MatDerivs = (
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    grad1 = (s1r * x / r, s1r * y / r, s1z)
    grad2 = (s2r * x / r, s2r * y / r, Fraction(0))
    terms = (
        divergence_scalar_outer(s2, grad2, er, der, et, det),
        divergence_scalar_outer(s2, grad2, et, det, er, der),
        divergence_scalar_outer(s1, grad1, er, der, ez, zeros),
        divergence_scalar_outer(s1, grad1, ez, zeros, er, der),
    )
    return _vadd(*terms)


def cylindrical_to_cartesian(
    force: Sequence[Fraction | int],
    x: Fraction | int,
    y: Fraction | int,
    r: Fraction | int,
) -> Vec3:
    """Map exact cylindrical (r,theta,z) components to Cartesian components."""
    if len(force) != 3:
        raise ValueError("force must have three components")
    er, et, ez, _, _ = cylindrical_basis(x, y, r)
    fr, ft, fz = (_q(v, "force") for v in force)
    return _vadd(_vscale(fr, er), _vscale(ft, et), _vscale(fz, ez))
