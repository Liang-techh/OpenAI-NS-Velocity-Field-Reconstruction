"""Conservative reconstruction ledger, not a theorem prover.

The audit surface is derived from ``status.construction_status`` so ``audit``
and ``status`` cannot silently disagree about which construction boundaries have
landed. A green test run, caller flag, or residual-defined force still cannot
upgrade the truth status.
"""
from __future__ import annotations
import copy

from .status import SOURCE_PINS, construction_status

SOURCES = {
    "paper_title": "Finite Time Blowup for Navier-Stokes",
    "paper_url": SOURCE_PINS["paper_url"],
    "paper_sha256": None,
    "paper_hash_status": "not computed; equations inspected in rendered PDF pages",
    "upstream_repository": "https://github.com/" + SOURCE_PINS["upstream_lean_repository"],
    "upstream_commit": SOURCE_PINS["upstream_commit_metadata_observed"],
    "upstream_observed_date": "2026-09-10",
    "lean_build_verified": SOURCE_PINS["lean_compiled"],
    "lean_crosscheck_status": "commit pinned only; theorem mapping and local Lean build pending",
}

_STAGE_SOURCES = {
    0: "Eq. (4.1)",
    1: "Thm. 4.6; Appendices A/B/C; heat component A.32-A.38",
    2: "Section 5; Eqs. (5.1)-(5.8); Eq. (5.27)",
    3: "Sections 6-7; Eqs. (6.1)-(6.11), (7.1)-(7.11)",
    4: "Section 7; Eqs. (7.24)-(7.30)",
    5: "Section 8",
    6: "Section 9",
    7: "Section 10; Eqs. (10.4), (10.20)-(10.21); endpoint Taylor-Borel extension",
    8: "tests and diagnostic reports",
}

_PAPER_EXACT_STATUSES = {"paper-exact", "paper-exact-formula"}


class IncompleteReconstructionError(RuntimeError):
    pass


def status_report() -> dict:
    """Return the audit view of the same fail-closed truth surface as ``status``."""
    runtime = construction_status()
    stages = []
    for stage in runtime["stages"]:
        entry = {
            "stage": stage["id"],
            "name": stage["name"],
            "status": stage["status"],
            "complete": stage["status"] in _PAPER_EXACT_STATUSES,
            "source": _STAGE_SOURCES[stage["id"]],
        }
        if "implemented" in stage:
            entry["implemented"] = stage["implemented"]
        if "remaining" in stage:
            entry["remaining"] = stage["remaining"]
        stages.append(entry)

    return {
        "schema_version": 1,
        "full_reconstruction": runtime["paper_exact_velocity_available"],
        "paper_exact_velocity_available": runtime["paper_exact_velocity_available"],
        "numerical_checks_are_proofs": False,
        "sources": copy.deepcopy(SOURCES),
        "stages": copy.deepcopy(stages),
        "blockers": [stage["name"] for stage in stages if not stage["complete"]],
    }


def require_complete_reconstruction() -> None:
    report = status_report()
    if not report["full_reconstruction"]:
        raise IncompleteReconstructionError(
            "Full paper reconstruction is unavailable: " + "; ".join(report["blockers"])
        )
