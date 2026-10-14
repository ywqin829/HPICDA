"""HPICDA stage runner (v1.1 blueprint sections 57, 70, 71).

Runs a single stage script, enforces human gates, writes checkpoint, and prints
the mandatory per-stage report. Never auto-passes a human gate.
"""
import os
import sys
import subprocess
import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkpoint import write_checkpoint, is_stale, read_checkpoint  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG = os.path.join(ROOT, "00_protocol", "project_config.yaml")
PROTOCOL = os.path.join(ROOT, "00_protocol", "MASTER_EXECUTION_BLUEPRINT_v1.1.md")
LOGS = os.path.join(ROOT, "19_logs")
os.makedirs(LOGS, exist_ok=True)

# stage_id -> (script, human_gate_required or None)
HUMAN_GATES = {
    "stage0_decision": "HG0_stage0_feasibility",
    "stage1_ontology": "HG1_endpoint_direction_semantics",
    "stage2_ingest": "HG2_independence_shared_control",
    "stage12_freeze": "HG3_biological_candidate_freeze",
    "stage13_holdout": "HG4_holdout_unlock",
}


def gate_approved(gate_name):
    p = os.path.join(ROOT, "18_human_reviews", f"{gate_name}.approved")
    return os.path.exists(p)


def run_stage(stage_id, script_rel, inputs=None, upstream=None, extra_args=None):
    inputs = inputs or []
    upstream = upstream or []
    script = os.path.join(ROOT, script_rel)
    gate = HUMAN_GATES.get(stage_id)
    gate_status = "NONE"
    if gate:
        gate_status = "APPROVED" if gate_approved(gate) else "PENDING"
        if gate_status != "APPROVED":
            print(f"STAGE: {stage_id}\nSTATUS: WAITING_HUMAN_REVIEW\nHUMAN_GATE_STATUS: {gate}={gate_status}")
            return 2

    log = os.path.join(LOGS, f"{stage_id}.log")
    cmd = [sys.executable, script] + (extra_args or [])
    with open(log, "w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=ROOT)
    status = "done" if proc.returncode == 0 else "failed"
    outputs = []
    rec = write_checkpoint(stage_id, script, inputs, CONFIG, PROTOCOL, outputs, upstream, status=status)
    print("\n".join([
        f"STAGE: {stage_id}",
        f"STATUS: {status.upper()}",
        f"INPUTS: {len(inputs)}",
        f"OUTPUTS: {len(outputs)}",
        f"GATES: {'gate ' + gate + '=' + gate_status if gate else 'none'}",
        f"HUMAN_GATE_STATUS: {gate_status}",
        f"NEXT_STAGE: {'proceed' if status == 'done' else 'halt'}",
        f"LOG: {log}",
    ]))
    return proc.returncode


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: run_stage.py <stage_id> <script_rel> [args...]")
        sys.exit(64)
    sys.exit(run_stage(sys.argv[1], sys.argv[2], extra_args=sys.argv[3:]))
