# R13 — CRISPRcleanR CNV 校正(评审 §6 技术对照)
# 设计(预注册于脚本头):
#   GSE178699 (TKOv3 文库, hg19 坐标来自 Target.Site):
#     每细胞系: 感染(3 重复) vs mock(3 重复均值) → Δ;mock(3 重复) vs T0 → β_mock
#   GSE210689 (Brunello 文库, 100% 序列重叠验证, GRCh38 坐标来自 Addgene 注释):
#     Lat vs UI;Pro vs UI(UI 均值作参照)
#   校正: ccr.NormfoldChanges(min_reads=30, ScalingByTotalReads) → ccr.GWclean(nperm=1000)
#   基因级: ccr.geneMeanFCs
#   输出: 校正前后 gene LFC,供 39_r13_overlay.py 做经验零假设与候选叠加
# 注意: NormfoldChanges 以第 1 个样本列为对照(ncontrols=1);mock/T0/UI 用重复均值(取整)
suppressMessages(library(CRISPRcleanR))
ROOT <- "/workspace/HPICDA"
OUT <- file.path(ROOT, "16_revision_r13_crisprcleanr")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()

# ---------- 1. GSE178699 注释(TKOv3, Target.Site 解析) ----------
cat("== GSE178699 ==\n")
raw <- read.delim(gzfile(file.path(ROOT, "03_functional_screens/raw/GSE178699/GSE178699_covid_screens_readCounts.txt.gz")),
                  check.names = FALSE)
# Target.Site: chr19:58863885-58863904_A1BG_+
ts <- regmatches(raw$Target.Site, regexec("^chr([^:]+):(\\d+)-(\\d+)_(.+)_(.)$", raw$Target.Site))
ts <- do.call(rbind, lapply(ts, function(x) if (length(x) == 6) x[2:6] else rep(NA, 5)))
ann178 <- data.frame(CODE = raw$GENE_CLONE, GENES = raw$GENE,
                     EXONE = NA, CHRM = ts[, 1], STRAND = ts[, 5],
                     STARTpos = as.numeric(ts[, 2]), ENDpos = as.numeric(ts[, 3]),
                     seq = sub(".*_", "", raw$GENE_CLONE), stringsAsFactors = FALSE)
ann178 <- ann178[!is.na(ann178$STARTpos) & !is.na(ann178$CHRM), ]
# 排除多基因映射 guide(同一 seq 对应多个基因)
dup_seq <- ann178$seq[duplicated(ann178$seq) & !is.na(ann178$seq)]
ann178 <- ann178[!(ann178$seq %in% dup_seq & ann178$seq != ""), ]
cat("GSE178699 注释:", nrow(ann178), "guides\n")

cell_lines <- list(
  VERO     = list(t0 = "VERO.T0",     mock = c("VERO.T7A", "VERO.T7B", "VERO.T7C"),
                  inf = c("VERO.T16A.SARS.IC90", "VERO.T16B.SARS.IC90", "VERO.T16C.SARS.IC90")),
  HuH7     = list(t0 = "HuH7.T0",     mock = c("HuH7.T14A", "HuH7.T14B", "HuH7.T14C"),
                  inf = c("HuH7.T15A.SARS2", "HuH7.T15B.SARS2", "HuH7.T15C.SARS2")),
  CALU3    = list(t0 = "CALU3.T0",    mock = c("CALU3.T25A", "CALU3.T25B", "CALU3.T25C"),
                  inf = c("CALU3.T43A.lowMOI.SARS", "CALU3.T43B.lowMOI.SARS", "CALU3.T43C.lowMOI.SARS")),
  CACO2    = list(t0 = "CACO2.T0",    mock = c("CACO2.T39A", "CACO2.T39B", "CACO2.T39C"),
                  inf = c("CACO2.T37A.MOI10", "CACO2.T37B.MOI10", "CACO2.T37C.MOI10")),
  HEK293A_T = list(t0 = "HEK293A_T.T0", mock = c("HEK293A_T.T12A", "HEK293A_T.T12B", "HEK293A_T.T12C"),
                  inf = c("HEK293A_T.T12A.MOI2", "HEK293A_T.T12B.MOI2", "HEK293A_T.T12C.MOI2")))

run_contrast <- function(counts, ann, ctrl_cols, trt_cols, tag, outdir, min_reads = 30) {
  # 对照 = 重复均值(取整);处理 = 各重复
  ctrl <- round(rowMeans(counts[, ctrl_cols, drop = FALSE]))
  df <- data.frame(sgRNA = ann$CODE, gene = ann$GENES, ctrl = ctrl,
                   counts[, trt_cols, drop = FALSE], check.names = FALSE)
  colnames(df)[4:ncol(df)] <- paste0("trt", seq_along(trt_cols))
  ann_idx <- ann; rownames(ann_idx) <- ann_idx$CODE  # logFCs2chromPos 需 CODE 行名
  fnc <- ccr.NormfoldChanges(Dframe = df, display = FALSE, outdir = outdir,
                             min_reads = min_reads, EXPname = tag, libraryAnnotation = ann_idx)
  fcs <- fnc$logFCs
  sorted <- ccr.logFCs2chromPos(fcs, ann_idx)
  sorted <- sorted[!is.na(sorted$CHR) & is.finite(sorted$avgFC), ]
  corrected <- ccr.GWclean(sorted, label = tag, display = FALSE, saveTO = NULL,
                           nperm = 1000, verbose = 1)
  corr <- corrected$corrected_logFCs   # 列: CHR,startp,endp,genes,avgFC(原始),BP,correction,correctedFC(校正后),guideIdx
  # 原始 guide 级 avgFC(重复均值)
  raw_avg <- rowMeans(fcs[, grep("^trt", colnames(fcs)), drop = FALSE], na.rm = TRUE)
  guide_out <- data.frame(sgrna = fcs$sgRNA, gene = fcs$gene,
                          raw = as.numeric(raw_avg), stringsAsFactors = FALSE)
  guide_out <- guide_out[match(rownames(corr), guide_out$sgrna), ]
  guide_out$corrected <- as.numeric(corr$correctedFC)
  guide_out <- guide_out[!is.na(guide_out$sgrna), ]
  write.table(guide_out, file.path(outdir, paste0("guide_lfc_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # 基因级均值(校正前后)
  gm <- aggregate(cbind(raw, corrected) ~ gene, data = guide_out, FUN = mean, na.rm = TRUE)
  colnames(gm) <- c("gene", "raw", "corrected")
  write.table(gm, file.path(outdir, paste0("gene_lfc_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(tag, "完成;基因数:", nrow(gm), "\n")
  invisible(gm)
}

dir178 <- file.path(OUT, "gse178699"); dir.create(dir178, showWarnings = FALSE)
counts178 <- raw[match(ann178$CODE, raw$GENE_CLONE), ]
for (cl in names(cell_lines)) {
  cfg <- cell_lines[[cl]]
  run_contrast(counts178, ann178, cfg$mock, cfg$inf, paste0(cl, "_inf_vs_mock"), dir178)
  run_contrast(counts178, ann178, cfg$t0, cfg$mock, paste0(cl, "_mock_vs_t0"), dir178)
}
cat("GSE178699 完成 (", round(difftime(Sys.time(), t0, units = "mins"), 1), " min)\n")

# ---------- 2. GSE210689(Brunello,坐标映射) ----------
cat("== GSE210689 ==\n")
mtb <- read.delim(file.path(ROOT, "03_functional_screens/raw/GSE210689/counts_mageck.txt"), check.names = FALSE)
br <- read.delim(file.path(OUT, "raw/brunello_library.txt"), check.names = TRUE, quote = "")
br <- br[!is.na(br$Position.of.Base.After.Cut..1.based.) & br$Position.of.Base.After.Cut..1.based. != "", ]
# NC_ 访问号 → 染色体名(去版本号映射;Addgene 注释为 GRCh37 版本号)
nc_base <- sub("\\..*$", "", br$Genomic.Sequence)
nc2chr <- c("NC_000001"="1","NC_000002"="2","NC_000003"="3","NC_000004"="4",
            "NC_000005"="5","NC_000006"="6","NC_000007"="7","NC_000008"="8",
            "NC_000009"="9","NC_000010"="10","NC_000011"="11","NC_000012"="12",
            "NC_000013"="13","NC_000014"="14","NC_000015"="15","NC_000016"="16",
            "NC_000017"="17","NC_000018"="18","NC_000019"="19","NC_000020"="20",
            "NC_000021"="21","NC_000022"="22","NC_000023"="X","NC_000024"="Y")
br$chr <- nc2chr[nc_base]
br <- br[!is.na(br$chr), ]
cutpos <- as.numeric(br$Position.of.Base.After.Cut..1.based.)
ann210 <- data.frame(CODE = br$sgRNA.Target.Sequence, GENES = br$Target.Gene.Symbol,
                     EXONE = NA, CHRM = br$chr,
                     STRAND = ifelse(br$Strand == "sense", "+", "-"),
                     STARTpos = cutpos - 10, ENDpos = cutpos + 10,
                     seq = br$sgRNA.Target.Sequence, stringsAsFactors = FALSE)
dup2 <- ann210$CODE[duplicated(ann210$CODE)]
ann210 <- ann210[!(ann210$CODE %in% dup2), ]
mtb2 <- mtb[match(ann210$CODE, mtb$sgRNA), ]
mtb2 <- mtb2[!is.na(mtb2$sgRNA), ]
ann210 <- ann210[ann210$CODE %in% mtb2$sgRNA, ]
mtb2 <- mtb2[match(ann210$CODE, mtb2$sgRNA), ]
cat("GSE210689 注释:", nrow(ann210), "guides\n")

dir210 <- file.path(OUT, "gse210689"); dir.create(dir210, showWarnings = FALSE)
# GSE210689 测序深度极低(UI 臂中位 1 read/guide),min_reads 降为 1(探索性,低置信)
run_contrast(mtb2, ann210, c("UI1", "UI2", "UI3"), c("Lat1", "Lat2", "Lat3"), "Lat_vs_UI", dir210, min_reads = 1)
run_contrast(mtb2, ann210, c("UI1", "UI2", "UI3"), c("Pro1", "Pro2", "Pro3"), "Pro_vs_UI", dir210, min_reads = 1)

cat("全部完成,总耗时", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
