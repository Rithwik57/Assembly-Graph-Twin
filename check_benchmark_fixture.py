"""Run a tiny positive-control benchmark and validate its expected evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from benchmark import benchmark_evidence_fingerprint, run_benchmark


FIXTURE_DIR = Path(__file__).with_name("fixtures") / "known_bubble"


def run_fixture_check(output: Path | None = None) -> dict[str, Any]:
    reads = FIXTURE_DIR / "reads.fasta"
    reference = FIXTURE_DIR / "reference.fasta"
    options = {
        "read_paths": [reads],
        "reference_path": reference,
        "k": 3,
        "seed": 7,
        "max_reads": None,
    }
    first = run_benchmark(**options)
    second = run_benchmark(**options)
    first_fingerprint = benchmark_evidence_fingerprint(first)
    second_fingerprint = benchmark_evidence_fingerprint(second)

    expected = {
        "bubble_candidate_count_at_least": 1,
        "uniquely_labeled_candidate_count_at_least": 1,
        "label_semantics": "reference_compatibility_proxy",
        "deterministic_evidence": True,
    }
    failures: list[str] = []
    if first["metrics"]["candidate_count"] < expected["bubble_candidate_count_at_least"]:
        failures.append("positive control did not produce a bubble candidate")
    if first["metrics"]["uniquely_labeled_candidate_count"] < expected["uniquely_labeled_candidate_count_at_least"]:
        failures.append("no candidate received a unique reference-compatibility proxy label")
    if first["protocol"].get("label_semantics") != expected["label_semantics"]:
        failures.append("report does not identify reference labels as compatibility proxies")
    if first_fingerprint != second_fingerprint:
        failures.append("candidate evidence changed between identical seeded runs")

    result = {
        "check": "known_bubble_positive_control",
        "status": "passed" if not failures else "failed",
        "expected": expected,
        "observed": {
            "candidate_count": first["metrics"]["candidate_count"],
            "uniquely_labeled_candidate_count": first["metrics"]["uniquely_labeled_candidate_count"],
            "accuracy_on_resolved_labeled": first["metrics"]["accuracy_on_resolved_labeled"],
            "false_resolutions": first["metrics"]["false_resolutions"],
            "seed": first["protocol"]["seed"],
            "k": first["protocol"]["k"],
            "read_fingerprint_sha256": first["inputs"]["read_fingerprint_sha256"],
            "evidence_fingerprint_sha256": first_fingerprint,
        },
        "failures": failures,
        "benchmark_report": first,
    }
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).parent / "reports" / "checks" / "known_bubble_fixture_report.json",
        help="machine-readable check and benchmark report (default: %(default)s)",
    )
    args = parser.parse_args()
    try:
        result = run_fixture_check(args.output)
    except (OSError, ValueError, RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({key: value for key, value in result.items() if key != "benchmark_report"}, sort_keys=True))
    print(f"Wrote fixture benchmark report to {args.output}")
    if result["status"] != "passed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
