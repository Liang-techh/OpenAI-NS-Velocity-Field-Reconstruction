"""CLI entry point combining the current audit/demo flow with v0.2 status/verify."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # PR #5's diagnostic CLI remains available without replacing the current
    # provenance-aware audit/demo semantics.
    if argv and argv[0] in {"status", "verify"}:
        from .cli import main as diagnostic_main
        return diagnostic_main(argv)

    from .provenance import (
        status_report,
        require_complete_reconstruction,
        IncompleteReconstructionError,
    )
    parser = argparse.ArgumentParser(prog="ns-reconstruct")
    sub = parser.add_subparsers(dest="command", required=True)
    audit = sub.add_parser("audit", help="inspect the conservative completion/provenance ledger")
    audit.add_argument("--output", type=Path, help="optional JSON report file")
    audit.add_argument("--require-paper-exact", action="store_true",
                       help="exit 2 if any required stage is incomplete")
    demo = sub.add_parser("demo", help="export explicitly labelled toy diagnostics")
    demo.add_argument("--output", type=Path, default=Path("results"))
    demo.add_argument("--require-paper-exact", action="store_true",
                      help="refuse incomplete reconstruction before writing data")
    args = parser.parse_args(argv)
    try:
        if args.command == "audit":
            report = status_report()
            if args.output is not None:
                from .demo import write_json
                write_json(args.output, report)
            print(json.dumps(report, indent=2, allow_nan=False))
            if args.require_paper_exact:
                require_complete_reconstruction()
            return 0

        if args.require_paper_exact:
            require_complete_reconstruction()
        from .demo import generate_demo
        report = generate_demo(args.output)
        print(json.dumps({
            "report": str(args.output / "report.json"),
            "full_reconstruction": report["full_reconstruction"],
            "all_diagnostic_checks_passed": report["all_diagnostic_checks_passed"],
        }, indent=2))
        return 0 if report["all_diagnostic_checks_passed"] else 1
    except IncompleteReconstructionError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except (ValueError, ArithmeticError, OSError) as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
