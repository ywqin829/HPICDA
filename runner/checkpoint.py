"""HPICDA checkpoint framework (v1.1 blueprint section 59).

Each stage writes a checkpoint file under 18_checkpoints/ recording SHA256 of
script, inputs, config, protocol, and outputs. A stage is stale if any upstream
hash changed; downstream stages are then invalidated (no false-complete).
"""
import hashlib
import json
import os
import datetime

CKPT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "18_checkpoints")
os.makedirs(CKPT_DIR, exist_ok=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _ckpt_path(stage_id):
    return os.path.join(CKPT_DIR, f"{stage_id}.json")


def write_checkpoint(stage_id, script_path, input_paths, config_path, protocol_path,
                     output_paths, upstream_stage_ids, status="done", extra=None):
    rec = {
        "stage_id": stage_id,
        "status": status,  # done / failed / skipped
        "script_sha256": sha256_file(script_path) if script_path and os.path.exists(script_path) else None,
        "input_sha256": {p: sha256_file(p) for p in input_paths if os.path.exists(p)},
        "config_sha256": sha256_file(config_path) if config_path and os.path.exists(config_path) else None,
        "protocol_sha256": sha256_file(protocol_path) if protocol_path and os.path.exists(protocol_path) else None,
        "output_sha256": {p: sha256_file(p) for p in output_paths if os.path.exists(p)},
        "upstream": {s: read_checkpoint(s).get("output_sha256") for s in upstream_stage_ids if read_checkpoint(s)},
        "timestamp_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "protocol_version": "1.1",
    }
    if extra:
        rec.update(extra)
    with open(_ckpt_path(stage_id), "w") as f:
        json.dump(rec, f, indent=2)
    # marker files per blueprint
    for ext in (".done", ".failed", ".skipped"):
        p = os.path.join(CKPT_DIR, f"{stage_id}{ext}")
        if os.path.exists(p):
            os.remove(p)
    open(os.path.join(CKPT_DIR, f"{stage_id}.{status}"), "w").write(rec["timestamp_utc"])
    return rec


def read_checkpoint(stage_id):
    p = _ckpt_path(stage_id)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return None


def is_stale(stage_id):
    """A stage is stale if missing, failed, or any upstream output hash changed."""
    rec = read_checkpoint(stage_id)
    if rec is None or rec.get("status") != "done":
        return True
    for up_id, up_hash in (rec.get("upstream") or {}).items():
        up_rec = read_checkpoint(up_id)
        if up_rec is None or up_rec.get("status") != "done":
            return True
        if up_hash is not None and up_rec.get("output_sha256") != up_hash:
            return True
    return False
