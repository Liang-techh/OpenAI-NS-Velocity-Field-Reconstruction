"""Conservative reconstruction ledger, not a theorem prover.

A green test run, a caller's paper_exact flag, or a residual-defined force must
never upgrade this ledger. Missing constructors need actual source-mapped
implementations and independent evidence before their entries can be changed.
"""
from __future__ import annotations
import copy

SOURCES = {
    "paper_title": "Finite Time Blowup for Navier-Stokes",
    "paper_url": "https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf",
    "paper_sha256": None,
    "paper_hash_status": "not computed; equations inspected in rendered PDF pages",
    "upstream_repository": "https://github.com/openai/NavierStokesAndEuler",
    "upstream_commit": "f9e8bc5b38b6e212696e8a30e3e91517af887bbd",
    "upstream_observed_date": "2026-09-10",
    "lean_build_verified": False,
    "lean_crosscheck_status": "commit pinned only; theorem mapping and local Lean build pending",
}

_STAGES = (
    (0, "similarity geometry", "implemented-formula", True, "Eq. (4.1)"),
    (1, "complete leading profile", "partial: kinematics, natural-axis rescaling and heat exterior", False,
     "Thm. 4.6; Appendices A/B/C; heat component A.32-A.38"),
    (2, "recursive background coefficients", "assembly only; recursive solver missing", False,
     "Section 5"),
    (3, "charts, support separation and transported phases", "partial: fixed dyadic geometry only", False,
     "Eqs. (6.1)-(6.6); remaining Section 6 missing"),
    (4, "oscillatory stress realization", "missing", False, "Section 7"),
    (5, "compact mean corrections", "missing", False, "Section 8"),
    (6, "residual-improvement iteration", "missing", False, "Section 9"),
    (7, "final compact field and smooth force", "composition only; paper inputs missing", False,
     "Section 10"),
    (8, "independent verification", "partial: numerical regression, no full-field certificate", False,
     "tests and diagnostic reports"),
)


class IncompleteReconstructionError(RuntimeError):
    pass


def status_report() -> dict:
    stages = [{"stage": i, "name": name, "status": status, "complete": complete,
               "source": source} for i, name, status, complete, source in _STAGES]
    return {"schema_version": 1, "full_reconstruction": all(s["complete"] for s in stages),
            "numerical_checks_are_proofs": False, "sources": copy.deepcopy(SOURCES),
            "stages": stages,
            "blockers": [s["name"] for s in stages if not s["complete"]]}


def require_complete_reconstruction() -> None:
    report = status_report()
    if not report["full_reconstruction"]:
        raise IncompleteReconstructionError("Full paper reconstruction is unavailable: "
                                             + "; ".join(report["blockers"]))
