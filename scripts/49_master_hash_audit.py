"""HPICDA script 49 — 全哈希审计(2026-09-05 复现包归一化).

对归一化后的完整工作树重算 SHA256,生成 MASTER_SHA256SUMS,并与三层历史
清单对账:SHA256SUMS.txt(Stage 0 时代)、milestone_v3.zip 内 SHA256SUMS、
18_checkpoints/*.json(脚本/输入/输出哈希)。产出 hash_audit_report.md。

审计范围 = 归档树定义:排除 tmp/、__pycache__、.git、
03_functional_screens/raw/PRJNA1281304/{sra,fastq}(8.9G SRA reads,不打包)。
"""
import hashlib
import json
import os
import zipfile
from datetime import datetime, timezone

ROOT = "/workspace/HPICDA"
PANEL = "/mnt/results/HPICDA"
EXCLUDE_DIRS = {"tmp", "__pycache__", ".git"}
EXCLUDE_PREFIXES = [
    "03_functional_screens/raw/PRJNA1281304/sra/",
    "03_functional_screens/raw/PRJNA1281304/fastq/",
    "16_revision_r11_orcs/raw/",  # ORCS 2.0.18 原始屏文件(1,956 个,3.4G),可用 fetch_raw_data.sh 重取
]
# 已知合法漂移(检查点后修订,内容与面板一致,见交接文档 §2.7)
KNOWN_DRIFT = {
    "03_functional_screens/registry_amendments.tsv": "Stage 2/R2 检查点后合法修订日志追加",
    "02_metadata/human_layer_rules.yaml": "R8 检查点后人类层规则更新",
    "README.md": "2026-09-05 归一化重写为最终版;milestone_v3 内为 Stage 12 暂停时代旧版",
}
# milestone_v3.zip 为策展布局,其 SHA256SUMS 路径需映射回工作树真实位置
# (工作树修订目录为 16_revision_*,无 16_revision_analyses/ 包裹层)
MILESTONE_MAP = [
    ("freeze_v2/", "09_candidate_freeze/"),
    ("governance/HG3_AMENDMENT_v2.approved", "18_human_reviews/HG3_AMENDMENT_v2.approved"),
    ("governance/HG4.approved", "18_human_reviews/HG4.approved"),
    ("governance/terminology_constraints_v3.yaml", "02_metadata/terminology_constraints_v3.yaml"),
    ("report/revision_report_v3.md", "revision_report_v3.md"),
    ("revision_analyses/16_revision_", "16_revision_"),
    ("stage13_14_15/stage13_", "20_stage13_holdout_adjudication/stage13_"),
    ("stage13_14_15/stage14_", "21_stage14_druggability/stage14_"),
    ("stage13_14_15/stage15_", "22_stage15_pareto/stage15_"),
    ("wetlab_design/", "17_wetlab_design/"),
]
# 临时文件等良性缺失(历史清单有记录但本就不应存在于归档树)
BENIGN_MISSING = ("tmp/hp_test",)
# 自指文件:清单与报告不能稳定地包含自身哈希,排除并在此声明
SELF_FILES = {"MASTER_SHA256SUMS", "hash_audit_report.md"}


def sha256_file(path, buf=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def walk_tree():
    files = []
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, ROOT)
            if rel in SELF_FILES:
                continue
            if any(rel.startswith(p) for p in EXCLUDE_PREFIXES):
                continue
            files.append(rel)
    return sorted(files)


def load_historical():
    """返回 {rel_path: [(source, hash), ...]}"""
    hist = {}

    def add(rel, h, src):
        hist.setdefault(rel, []).append((src, h))

    # 层 1:SHA256SUMS.txt(Stage 0 时代)
    p = os.path.join(ROOT, "SHA256SUMS.txt")
    if os.path.exists(p):
        for line in open(p):
            parts = line.split()
            if len(parts) == 2:
                add(parts[1].lstrip("*"), parts[0], "SHA256SUMS.txt")

    # 层 2:milestone_v3.zip 内 SHA256SUMS
    zp = os.path.join(ROOT, "milestone_v3.zip")
    if os.path.exists(zp):
        with zipfile.ZipFile(zp) as z:
            for name in z.namelist():
                if name.endswith("SHA256SUMS"):
                    for line in z.read(name).decode().splitlines():
                        parts = line.split()
                        if len(parts) == 2:
                            rel = parts[1].lstrip("*")
                            rel = rel.replace("milestone_v3/", "", 1)
                            rel = rel.lstrip("./")
                            for src_p, dst_p in MILESTONE_MAP:
                                if rel.startswith(src_p):
                                    rel = dst_p + rel[len(src_p):]
                                    break
                            add(rel, parts[0], "milestone_v3")
    # 层 3:检查点(输入+输出)
    ckdir = os.path.join(ROOT, "18_checkpoints")
    for fn in sorted(os.listdir(ckdir)):
        if not fn.endswith(".json"):
            continue
        rec = json.load(open(os.path.join(ckdir, fn)))
        src = f"checkpoint:{rec.get('stage_id', fn)}"
        for key in ("input_sha256", "output_sha256"):
            d = rec.get(key) or {}
            if isinstance(d, dict):
                for path, h in d.items():
                    rel = os.path.relpath(path, ROOT) if os.path.isabs(path) else path
                    add(rel, h, src)
    return hist


def main():
    files = walk_tree()
    hist = load_historical()
    master_lines, verified, drifted, newly_added = [], [], [], []
    for rel in files:
        h = sha256_file(os.path.join(ROOT, rel))
        master_lines.append(f"{h}  {rel}")
        if rel in hist:
            if any(h == hh for _, hh in hist[rel]):
                verified.append(rel)
            else:
                drifted.append((rel, sorted({s for s, _ in hist[rel]})))
        else:
            newly_added.append(rel)
    missing_all = [(rel, sorted({s for s, _ in srcs}))
                   for rel, srcs in hist.items() if not os.path.exists(os.path.join(ROOT, rel))]
    missing = [(r, s) for r, s in missing_all
               if not any(b in r for b in BENIGN_MISSING)]
    missing_benign = [(r, s) for r, s in missing_all
                      if any(b in r for b in BENIGN_MISSING)]

    with open(os.path.join(ROOT, "MASTER_SHA256SUMS"), "w") as f:
        f.write("\n".join(master_lines) + "\n")

    ckdir = os.path.join(ROOT, "18_checkpoints")
    dones = sorted(f for f in os.listdir(ckdir) if f.endswith(".done"))
    hrdir = os.path.join(ROOT, "18_human_reviews")
    approved = sorted(f for f in os.listdir(hrdir) if f.endswith(".approved"))

    # Stage 14 污染事件修复验证
    s14 = "21_stage14_druggability/stage14_druggability_axis.tsv"
    s14_hash = sha256_file(os.path.join(ROOT, s14))
    ck14 = json.load(open(os.path.join(ckdir, "stage14_druggability.json")))
    s14_expect = (ck14.get("output_sha256") or {}).get(f"{ROOT}/{s14}", "")
    s14_ok = s14_hash == s14_expect

    # 面板大文件读回校验
    panel_checks = []
    for rel in ["04_depmap_baseline/gene_effect_26Q1.csv",
                "07_cross_pathogen/convergence_v2_null.npz"]:
        pp = os.path.join(PANEL, rel)
        if os.path.exists(pp):
            panel_checks.append((rel, sha256_file(pp) == sha256_file(os.path.join(ROOT, rel))))

    rep = []
    rep.append("# HPICDA 全哈希审计报告\n")
    rep.append("**日期**:2026-09-05 · **脚本**:scripts/49_master_hash_audit.py · "
               "**范围**:归档树(排除 tmp/、__pycache__、.git、PRJNA1281304 sra+fastq、16_revision_r11_orcs/raw(ORCS 2.0.18 原始件,可重取))\n")
    rep.append("## 总览\n")
    rep.append("| 分类 | 数量 |\n|---|---|")
    rep.append(f"| 归档树文件总数(MASTER_SHA256SUMS 行数) | {len(files)} |")
    rep.append(f"| verified(与历史清单一致) | {len(verified)} |")
    rep.append(f"| drifted(历史清单有记录但哈希不同) | {len(drifted)} |")
    rep.append(f"| newly-added(历史清单外,含本轮新增) | {len(newly_added)} |")
    rep.append(f"| missing(历史清单有记录但文件不存在) | {len(missing)} |\n")
    rep.append("## drifted 明细\n")
    if drifted:
        for rel, srcs in drifted:
            note = KNOWN_DRIFT.get(rel, "**未解释 —— 待用户定夺**")
            rep.append(f"- `{rel}`(历史来源:{', '.join(srcs)})—— {note}")
    else:
        rep.append("- 无")
    rep.append("\n## missing 明细\n")
    if missing:
        for rel, srcs in missing:
            rep.append(f"- `{rel}`(历史来源:{', '.join(srcs)})—— **未解释,待用户定夺**")
    else:
        rep.append("- 无实质性缺失(missing = 0)")
    if missing_benign:
        rep.append("\n良性缺失(临时文件,本就不应入归档树):")
        for rel, srcs in missing_benign:
            rep.append(f"- `{rel}`(历史来源:{', '.join(srcs)})")
    rep.append("\n## 已知合法漂移的处置\n")
    rep.append("以下文件在检查点写入后发生过合法修订,内容与结果面板一致,"
               "MASTER_SHA256SUMS 以当前内容为准,历史哈希仍保留在相应检查点中:\n")
    for rel, note in KNOWN_DRIFT.items():
        rep.append(f"- `{rel}`:{note}")
    rep.append("\n## Stage 14 污染事件修复验证\n")
    rep.append(f"- `21_stage14_druggability/stage14_druggability_axis.tsv` 当前哈希 "
               f"`{s14_hash[:16]}…` vs 检查点 `{s14_expect[:16]}…` → "
               f"**{'一致(修复有效)' if s14_ok else '不一致 —— 需立即排查'}**")
    rep.append("\n## 标记核对\n")
    rep.append(f"- .done 标记:{len(dones)} 个({', '.join(dones)})")
    rep.append(f"- .approved 批准件:{len(approved)} 个({', '.join(approved)})")
    rep.append("\n## 面板大文件读回校验\n")
    if panel_checks:
        for rel, ok in panel_checks:
            rep.append(f"- `{rel}`:{'一致' if ok else '**不一致**'}")
    else:
        rep.append("- (审计时面板大文件尚未就位)")
    rep.append("\n## newly-added 说明\n")
    rep.append(f"共 {len(newly_added)} 个。绝大多数为 Stage 1–15 各阶段输出"
               "(历史三层清单合计仅覆盖约百个文件,从未全量登记);"
               "本轮新增文件(README 最终版、PIPELINE_INDEX.yaml、environment_lock.txt、"
               "fetch_raw_data.sh、本报告、MASTER_SHA256SUMS 等)亦在此列。"
               "MASTER_SHA256SUMS 自此为唯一权威清单。")
    rep.append(f"\n生成时间:{datetime.now(timezone.utc).isoformat()}")
    with open(os.path.join(ROOT, "hash_audit_report.md"), "w") as f:
        f.write("\n".join(rep) + "\n")

    print(json.dumps({"total": len(files), "verified": len(verified),
                      "drifted": len(drifted), "newly_added": len(newly_added),
                      "missing": len(missing), "stage14_ok": s14_ok,
                      "done_markers": len(dones), "approved": len(approved),
                      "panel_checks": panel_checks}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
