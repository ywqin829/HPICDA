"""HPICDA pipeline runner — index-driven (2026-09-05 normalization).

Replaces the Stage-0-era skeleton. The authoritative execution record is
runner/PIPELINE_INDEX.yaml + 18_checkpoints/*.json + the execution notebook.

Honesty note: the full pipeline was NEVER rerun unattended end-to-end.
Stages marked rerunnable=yes can be re-executed standalone from packaged
inputs; rerunnable=fetch requires 03_functional_screens/raw/fetch_raw_data.sh
first; rerunnable=manual involves human gates / interactive revision and
will NOT be run by this driver.

Usage:
  python runner/run_pipeline.py --list                 # show index + status
  python runner/run_pipeline.py --verify               # re-hash outputs vs checkpoints
  python runner/run_pipeline.py --run <stage_id>       # rerun one rerunnable=yes stage
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from checkpoint import read_checkpoint, sha256_file  # noqa: E402


def load_index():
    import yaml
    with open(os.path.join(HERE, "PIPELINE_INDEX.yaml"), encoding="utf-8") as f:
        return yaml.safe_load(f)


def cmd_list(idx):
    for s in idx["stages"]:
        ck = s.get("checkpoint")
        status = "no-checkpoint"
        if ck:
            rec = read_checkpoint(str(ck))
            status = rec.get("status") if rec else "checkpoint-missing"
        print(f"{s['stage_id']:32s} [{s.get('rerunnable','?'):6s}] "
              f"gate={s.get('gate') or '-':14s} {status}  {s.get('script')}")


def cmd_verify(idx):
    """Re-hash every output recorded in checkpoints; report drift."""
    n_ok = n_bad = n_missing = 0
    for s in idx["stages"]:
        ck = s.get("checkpoint")
        if not ck:
            continue
        rec = read_checkpoint(str(ck))
        if not rec:
            continue
        for path, expect in (rec.get("output_sha256") or {}).items():
            p = path if os.path.isabs(path) else os.path.join(ROOT, path)
            if not os.path.exists(p):
                print(f"[missing] {s['stage_id']}: {path}")
                n_missing += 1
            elif sha256_file(p) != expect:
                print(f"[DRIFT]   {s['stage_id']}: {path}")
                n_bad += 1
            else:
                n_ok += 1
    print(f"\nverify: ok={n_ok} drift={n_bad} missing={n_missing}")
    return 0 if (n_bad == 0 and n_missing == 0) else 1


def cmd_run(idx, stage_id):
    stage = next((s for s in idx["stages"] if s["stage_id"] == stage_id), None)
    if stage is None:
        print(f"unknown stage_id {stage_id}; see --list")
        return 2
    if stage.get("rerunnable") != "yes":
        print(f"[halt] {stage_id} is rerunnable={stage.get('rerunnable')}: "
              f"{stage.get('note','requires manual steps / external data')}")
        return 2
    gate = stage.get("gate")
    if gate:
        g = idx["gates"].get(gate, {})
        gfile = os.path.join(ROOT, g.get("file", ""))
        if not os.path.exists(gfile):
            print(f"[halt] {stage_id} requires {gate} approval file: {g.get('file')}")
            return 2
    script = stage.get("script")
    if not script or "+" in script:
        print(f"[halt] {stage_id} has no single runnable script entry")
        return 2
    from run_stage import run_stage
    return run_stage(stage_id, script.strip())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--run", metavar="STAGE_ID")
    args = ap.parse_args()
    idx = load_index()
    if args.list:
        cmd_list(idx)
        return 0
    if args.verify:
        return cmd_verify(idx)
    if args.run:
        return cmd_run(idx, args.run)
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
