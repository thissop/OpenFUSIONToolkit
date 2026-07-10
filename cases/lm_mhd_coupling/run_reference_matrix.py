#!/usr/bin/env python3
"""Run the reduced LM-MHD reference artifact matrix.

This is a reproducibility driver for diagnostic scripts that are already backed
by unit tests. It records command status and elapsed time; it does not turn the
diagnostics into validation claims.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class ReferenceCommand:
    name: str
    path: Path
    args: tuple[str, ...] = ()


PLASMA_INTERFACE_COMMANDS = (
    ReferenceCommand("plasma_interface_sheath_scan", Path("cases/lm_mhd_coupling/plasma_interface_sheath_scan.py")),
    ReferenceCommand(
        "plasma_interface_thermal_response_scan",
        Path("cases/lm_mhd_coupling/plasma_interface_thermal_response_scan.py"),
    ),
    ReferenceCommand(
        "plasma_interface_current_feedback_scan",
        Path("cases/lm_mhd_coupling/plasma_interface_current_feedback_scan.py"),
    ),
    ReferenceCommand(
        "plasma_interface_marangoni_scan",
        Path("cases/lm_mhd_coupling/plasma_interface_marangoni_scan.py"),
    ),
    ReferenceCommand(
        "plasma_interface_free_surface_gate_scan",
        Path("cases/lm_mhd_coupling/plasma_interface_free_surface_gate_scan.py"),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=("plasma-interface",), default="plasma-interface")
    parser.add_argument("--results", default="cases/lm_mhd_coupling/results", help="results directory")
    parser.add_argument("--keep-going", action="store_true", help="run remaining commands after a failure")
    return parser.parse_args()


def commands_for_suite(suite: str) -> tuple[ReferenceCommand, ...]:
    if suite == "plasma-interface":
        return PLASMA_INTERFACE_COMMANDS
    raise ValueError(f"unknown suite {suite!r}")


def run_command(command: ReferenceCommand, results: Path) -> dict[str, str | int | float]:
    log_path = results / f"reference_matrix_{command.name}.log"
    argv = [sys.executable, str(ROOT / command.path), *command.args]
    start = time.monotonic()
    with log_path.open("w", newline="\n") as log:
        log.write(f"$ {' '.join(argv)}\n\n")
        completed = subprocess.run(
            argv,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    elapsed = time.monotonic() - start
    return {
        "name": command.name,
        "command": " ".join(argv),
        "returncode": completed.returncode,
        "elapsed_s": elapsed,
        "log": str(log_path.relative_to(ROOT)),
    }


def write_run_csv(results: Path, rows: list[dict[str, str | int | float]]) -> None:
    with (results / "reference_matrix_runs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "command", "returncode", "elapsed_s", "log"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(results: Path, suite: str, rows: list[dict[str, str | int | float]]) -> None:
    failed = [row for row in rows if row["returncode"] != 0]
    with (results / "reference_matrix_summary.md").open("w", newline="\n") as handle:
        handle.write("# LM-MHD reference matrix run\n\n")
        handle.write(f"Suite: `{suite}`\n\n")
        handle.write("This driver reruns reduced diagnostics and records command status. ")
        handle.write("It is a reproducibility artifact, not a physics validation result.\n\n")
        handle.write("| command | return code | elapsed [s] | log |\n")
        handle.write("|---|---:|---:|---|\n")
        for row in rows:
            handle.write(
                f"| `{row['name']}` | `{row['returncode']}` | "
                f"`{float(row['elapsed_s']):.3f}` | `{row['log']}` |\n"
            )
        handle.write("\n")
        if failed:
            handle.write(f"Result: FAIL ({len(failed)} command(s) failed).\n")
        else:
            handle.write("Result: PASS (all commands returned zero).\n")
        handle.write("\nRun `cases/lm_mhd_coupling/lm_mhd_artifact_manifest.py` after this driver ")
        handle.write("to refresh artifact hashes and PNG sanity metadata.\n")


def main() -> None:
    args = parse_args()
    results = ROOT / args.results
    results.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str | int | float]] = []
    for command in commands_for_suite(args.suite):
        row = run_command(command, results)
        rows.append(row)
        print(f"{command.name}: returncode={row['returncode']} elapsed={float(row['elapsed_s']):.3f}s")
        if row["returncode"] != 0 and not args.keep_going:
            break
    write_run_csv(results, rows)
    write_summary(results, args.suite, rows)
    if any(row["returncode"] != 0 for row in rows):
        raise SystemExit(1)
    print(f"Wrote reference matrix summary to {results}")


if __name__ == "__main__":
    main()
