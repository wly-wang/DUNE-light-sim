#!/usr/bin/env python3

"""Prepare, submit, and inspect a scratch-based radial G4 campaign on Eddie."""

import argparse
import ast
import csv
import itertools
import json
import math
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


LOW_Z_MIN_CM = 0.0
LOW_Z_JOIN_CM = 300.0
HIGH_Z_JOIN_CM = 5500.0
HIGH_Z_MAX_CM = 5800.0

# X=0 and X=10 are omitted because fit_validate_eddie rejects X <= 10 cm.
X_VALUES_CM = (
    20.0, 30.0, 40.0, 50.0, 100.0, 150.0, 200.0, 250.0,
    280.0, 300.0, 310.0, 320.0, 330.0, 340.0, 350.0,
)

Y_VALUES_CM = (
    31.135, 46.703, 62.270, 77.838, 93.406, 108.973, 124.541,
    140.108, 155.676, 171.243, 186.811, 202.379, 217.946, 233.514,
    249.082, 264.649, 280.217, 295.784, 311.352, 326.920, 342.487,
    358.055, 373.623, 389.190, 404.758, 420.325, 435.893, 451.461,
    467.028, 482.596, 498.163, 513.731, 529.298, 544.866, 560.433,
    576.001, 591.569,
)

DEFAULT_SPACING_CM = 25.0
DEFAULT_JOB_LIMIT = 1000
PHOTONS_PER_EVENT = 100000
EVENTS_PER_ORIGIN = 100
PHOTON_ENERGY_EV = 9.69


def utc_now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def radial_z_values(spacing):
    width = LOW_Z_JOIN_CM - LOW_Z_MIN_CM
    n_intervals = round(width / spacing)
    if spacing <= 0.0 or not math.isclose(
        n_intervals * spacing, width, rel_tol=0.0, abs_tol=1.0e-9
    ):
        raise ValueError(
            f"Z spacing must divide the 300 cm radial depth exactly; got {spacing}"
        )

    offsets = [(index + 0.5) * spacing for index in range(n_intervals)]
    low_z = sorted(LOW_Z_JOIN_CM - offset for offset in offsets)
    high_z = sorted(HIGH_Z_JOIN_CM + offset for offset in offsets)
    return low_z + high_z


def read_existing_z(path):
    text = path.read_text()
    match = re.search(r"^Z\s*=\s*(\[[^\n]+\])", text, flags=re.MULTILINE)
    if not match:
        raise RuntimeError(f"Could not find a Z list in {path}")
    return [float(value) for value in ast.literal_eval(match.group(1))]


def git_commit(script_path):
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=script_path.parent.parent,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def ensure_empty_campaign(path):
    if path.exists() and any(path.iterdir()):
        raise RuntimeError(
            f"Campaign directory is not empty: {path}\n"
            "Choose a new campaign name; existing jobs will not be overwritten."
        )
    path.mkdir(parents=True, exist_ok=True)


def write_steering(path, positions):
    with path.open("w") as stream:
        stream.write(
            "  x      y     z      t      dx      dy      dz      dt      "
            "p         dp        n\n"
        )
        for x, y, z in positions:
            stream.write(
                f"  {x:.3f}    {y:.3f}    {z:.3f}    0.0    0.0    0.0    "
                f"0.0    0.0    {PHOTON_ENERGY_EV:.2f}    0.25    "
                f"{PHOTONS_PER_EVENT}\n"
            )


def shell_command(parts):
    return " ".join(shlex.quote(str(part)) for part in parts)


def write_job_script(
    path,
    job_index,
    job_directory,
    fcl_name,
    n_events,
    setup_script,
    container_image,
    apptainer,
):
    bind_paths = (
        "/cvmfs",
        str(Path.home()),
        str(job_directory.parent),
        "/opt",
        "/run/user",
        "/etc/hostname",
        "/etc/hosts",
        "/etc/krb5.conf",
    )
    inner = (
        f"source {shlex.quote(str(setup_script))}; "
        f"cd {shlex.quote(str(job_directory))}; "
        "echo Current directory: $(pwd); "
        "echo Starting full optical simulation; "
        f"lar -c ../{shlex.quote(fcl_name)} -n {n_events}"
    )
    command = shell_command(
        (
            apptainer,
            "exec",
            "--cleanenv",
            "-B",
            ",".join(bind_paths),
            container_image,
            "/bin/bash",
            "-lc",
            inner,
        )
    )

    with path.open("w") as stream:
        stream.write("#!/bin/bash\n")
        stream.write(f"#$ -N rd25_{job_index:04d}\n")
        stream.write("#$ -cwd\n")
        stream.write("#$ -l h_rt=10:00:00\n")
        stream.write("#$ -l h_vmem=8G\n")
        stream.write(f"#$ -o {job_directory}/sge.stdout\n")
        stream.write(f"#$ -e {job_directory}/sge.stderr\n")
        stream.write("set -u\n")
        stream.write("SECONDS=0\n")
        stream.write("rm -f SUCCESS FAILED\n")
        stream.write(f"{command} > log 2>&1\n")
        stream.write("exit_code=$?\n")
        stream.write("duration=$SECONDS\n")
        stream.write(
            'echo "Job finished with status ${exit_code} after '
            '${duration} seconds" >> log\n'
        )
        stream.write(
            'if [ "$exit_code" -eq 0 ] && [ -s semi_analytical_model.root ]; then\n'
        )
        stream.write("  touch SUCCESS\n")
        stream.write("else\n")
        stream.write("  touch FAILED\n")
        stream.write("fi\n")
        stream.write("exit $exit_code\n")
    path.chmod(0o755)


def prepare(args):
    script_path = Path(__file__).resolve()
    campaign = args.campaign_dir.expanduser().resolve()
    source_fcl = args.fcl.expanduser().resolve()
    existing_grid = args.existing_grid.expanduser().resolve()
    # Preserve /home/$USER in the generated command. On Eddie it may resolve to
    # /exports/eddie3_homes_local, which is not one of the explicit container
    # bind paths used by these jobs.
    setup_script = args.setup_script.expanduser().absolute()

    if not source_fcl.is_file():
        raise RuntimeError(f"FHiCL file not found: {source_fcl}")
    if not existing_grid.is_file():
        raise RuntimeError(f"Existing grid file not found: {existing_grid}")
    if not setup_script.is_file():
        raise RuntimeError(f"DUNE setup script not found: {setup_script}")

    candidate_z_values = radial_z_values(args.z_spacing)
    existing_z = read_existing_z(existing_grid)
    skipped_existing_z = []
    z_values = []
    for z in candidate_z_values:
        nearest = min(existing_z, key=lambda old_z: abs(old_z - z))
        if abs(z - nearest) <= args.duplicate_tolerance:
            skipped_existing_z.append(
                {
                    "candidate_z_cm": z,
                    "existing_z_cm": nearest,
                    "separation_cm": abs(z - nearest),
                }
            )
        else:
            z_values.append(z)
    if not z_values:
        raise RuntimeError("Every candidate Z plane overlaps the existing sample")

    ensure_empty_campaign(campaign)
    campaign_fcl = campaign / source_fcl.name
    shutil.copy2(source_fcl, campaign_fcl)

    positions = list(itertools.product(X_VALUES_CM, Y_VALUES_CM, z_values))
    random.Random(args.seed).shuffle(positions)
    chunk_size = math.ceil(len(positions) / args.job_limit)
    chunks = [
        positions[start:start + chunk_size]
        for start in range(0, len(positions), chunk_size)
    ]

    with (campaign / "positions.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("X_cm", "Y_cm", "Z_cm"))
        writer.writerows(positions)

    with (campaign / "z_values.txt").open("w") as stream:
        for z in sorted(z_values):
            nearest = min(existing_z, key=lambda old_z: abs(old_z - z))
            stream.write(
                f"{z:.3f}  nearest_existing={nearest:.3f}  "
                f"separation={abs(z - nearest):.3f}\n"
            )

    jobs = []
    for job_index, chunk in enumerate(chunks):
        job_directory = campaign / f"job{job_index:04d}"
        job_directory.mkdir()
        steering = job_directory / "myLightSourceSteering.txt"
        job_script = job_directory / "jobscript.sh"
        write_steering(steering, chunk)
        write_job_script(
            job_script,
            job_index,
            job_directory,
            campaign_fcl.name,
            len(chunk) * EVENTS_PER_ORIGIN,
            setup_script,
            args.container_image,
            args.apptainer,
        )
        jobs.append(
            {
                "index": job_index,
                "directory": str(job_directory),
                "positions": len(chunk),
                "events": len(chunk) * EVENTS_PER_ORIGIN,
            }
        )

    manifest = {
        "created_utc": utc_now(),
        "repository_commit": git_commit(script_path),
        "campaign_directory": str(campaign),
        "source_fcl": str(source_fcl),
        "campaign_fcl": str(campaign_fcl),
        "z_spacing_cm": args.z_spacing,
        "candidate_z_values_cm": sorted(candidate_z_values),
        "z_values_cm": sorted(z_values),
        "skipped_existing_z": skipped_existing_z,
        "x_values_cm": list(X_VALUES_CM),
        "y_values_cm": list(Y_VALUES_CM),
        "photons_per_event": PHOTONS_PER_EVENT,
        "events_per_origin": EVENTS_PER_ORIGIN,
        "photons_per_origin": PHOTONS_PER_EVENT * EVENTS_PER_ORIGIN,
        "random_seed": args.seed,
        "requested_job_limit": args.job_limit,
        "actual_jobs": len(jobs),
        "positions_per_full_job": chunk_size,
        "total_positions": len(positions),
        "total_events": len(positions) * EVENTS_PER_ORIGIN,
        "jobs": jobs,
    }
    with (campaign / "campaign_manifest.json").open("w") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")

    print("Campaign prepared; no jobs were submitted.")
    print(f"Campaign: {campaign}")
    low_z_count = sum(z < LOW_Z_JOIN_CM for z in z_values)
    high_z_count = sum(z > HIGH_Z_JOIN_CM for z in z_values)
    print(
        f"New Z planes: {len(z_values)} "
        f"({low_z_count} low-Z, {high_z_count} high-Z)"
    )
    if skipped_existing_z:
        for skipped in skipped_existing_z:
            print(
                "Skipped existing plane: "
                f"candidate {skipped['candidate_z_cm']:.3f} cm, "
                f"existing {skipped['existing_z_cm']:.3f} cm"
            )
    print(f"X values: {len(X_VALUES_CM)}")
    print(f"Y values: {len(Y_VALUES_CM)}")
    print(f"Total origins: {len(positions)}")
    print(f"Jobs: {len(jobs)} (up to {chunk_size} origins per job)")
    print(f"Photons per origin: {PHOTONS_PER_EVENT * EVENTS_PER_ORIGIN}")
    print("Inspect campaign_manifest.json, z_values.txt, and several job scripts.")
    print("When satisfied, run the submit subcommand first without --confirm.")


def load_manifest(campaign):
    manifest_path = campaign / "campaign_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"Campaign manifest not found: {manifest_path}")
    with manifest_path.open() as stream:
        return json.load(stream)


def submit(args):
    campaign = args.campaign_dir.expanduser().resolve()
    manifest = load_manifest(campaign)
    submission_path = campaign / (
        "retry_submissions.tsv" if args.retry_failed else "submissions.tsv"
    )
    already_submitted = set()
    initial_submission_path = campaign / "submissions.tsv"
    if initial_submission_path.exists():
        with initial_submission_path.open() as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if fields:
                    already_submitted.add(int(fields[0]))

    if args.retry_failed:
        pending = [
            job for job in manifest["jobs"]
            if (Path(job["directory"]) / "FAILED").exists()
        ]
    else:
        pending = [
            job for job in manifest["jobs"]
            if int(job["index"]) not in already_submitted
            and not (Path(job["directory"]) / "SUCCESS").exists()
        ]
    selected = pending
    if args.limit is not None:
        if args.limit <= 0:
            raise ValueError("--limit must be a positive integer")
        selected = pending[:args.limit]

    print(f"Prepared jobs: {len(manifest['jobs'])}")
    print(f"Already submitted: {len(already_submitted)}")
    if args.retry_failed:
        print(f"Failed jobs eligible for retry: {len(pending)}")
    else:
        print(f"Pending submission: {len(pending)}")
    print(f"Selected for this submission: {len(selected)}")
    if not args.confirm:
        print("Dry run only. Add --confirm to submit these jobs with qsub.")
        return

    with submission_path.open("a") as submission_stream:
        for job in selected:
            job_index = int(job["index"])
            job_directory = Path(job["directory"])
            result = subprocess.run(
                ["qsub", str(job_directory / "jobscript.sh")],
                cwd=job_directory,
                check=True,
                capture_output=True,
                text=True,
            )
            response = result.stdout.strip().replace("\t", " ")
            submission_stream.write(
                f"{job_index}\t{utc_now()}\t{response}\n"
            )
            submission_stream.flush()
            if args.retry_failed:
                (job_directory / "FAILED").unlink(missing_ok=True)
            print(response)


def status(args):
    campaign = args.campaign_dir.expanduser().resolve()
    manifest = load_manifest(campaign)
    counts = {
        "success": 0,
        "failed": 0,
        "root_without_marker": 0,
        "incomplete": 0,
    }
    completed_roots = []
    for job in manifest["jobs"]:
        job_directory = Path(job["directory"])
        output = job_directory / "semi_analytical_model.root"
        if (job_directory / "SUCCESS").exists() and output.is_file() and output.stat().st_size > 0:
            counts["success"] += 1
            completed_roots.append(output)
        elif (job_directory / "FAILED").exists():
            counts["failed"] += 1
        elif output.is_file() and output.stat().st_size > 0:
            counts["root_without_marker"] += 1
        else:
            counts["incomplete"] += 1

    file_list = campaign / "completed_root_files.txt"
    with file_list.open("w") as stream:
        for path in completed_roots:
            stream.write(f"{path}\n")

    print(f"Campaign: {campaign}")
    print(f"Total jobs: {len(manifest['jobs'])}")
    for label, count in counts.items():
        print(f"{label}: {count}")
    print(f"Merge-ready file list: {file_list}")


def build_parser():
    script_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="Create jobs without submitting")
    prepare_parser.add_argument("--campaign-dir", type=Path, required=True)
    prepare_parser.add_argument("--z-spacing", type=float, default=DEFAULT_SPACING_CM)
    prepare_parser.add_argument("--job-limit", type=int, default=DEFAULT_JOB_LIMIT)
    prepare_parser.add_argument("--seed", type=int, default=20261118)
    prepare_parser.add_argument("--duplicate-tolerance", type=float, default=0.5)
    prepare_parser.add_argument(
        "--fcl", type=Path,
        default=script_directory / "protoduneHD_v6_refactored.fcl",
    )
    prepare_parser.add_argument(
        "--existing-grid", type=Path,
        default=script_directory / "origin_grid_lists.txt",
    )
    prepare_parser.add_argument(
        "--setup-script", type=Path,
        default=Path.home() / "setup/dune_light_sim_v10_16_00d00_setup.sh",
    )
    prepare_parser.add_argument(
        "--apptainer",
        default="/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer",
    )
    prepare_parser.add_argument(
        "--container-image",
        default="/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest",
    )
    prepare_parser.set_defaults(func=prepare)

    submit_parser = subparsers.add_parser("submit", help="Submit prepared jobs")
    submit_parser.add_argument("--campaign-dir", type=Path, required=True)
    submit_parser.add_argument(
        "--limit", type=int,
        help="Submit at most this many pending jobs (use 1 for a pilot)",
    )
    submit_parser.add_argument(
        "--retry-failed", action="store_true",
        help="Submit only jobs currently carrying a FAILED marker",
    )
    submit_parser.add_argument("--confirm", action="store_true")
    submit_parser.set_defaults(func=submit)

    status_parser = subparsers.add_parser("status", help="Summarize campaign outputs")
    status_parser.add_argument("--campaign-dir", type=Path, required=True)
    status_parser.set_defaults(func=status)
    return parser


def main():
    args = build_parser().parse_args()
    try:
        args.func(args)
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
