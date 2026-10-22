#!/usr/bin/env python3
"""Build one review PDF containing cover letter, manuscript and figures.

The supplementary workbook is embedded as a PDF attachment. The full S5 call
table remains on Zenodo because it exceeds the initial-submission file limit.
"""

from __future__ import annotations

import html
import json
import os
import re
import tempfile
from pathlib import Path

import fitz
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
)


ROOT = Path(
    os.environ.get("HPICDA_ROOT", Path(__file__).resolve().parents[1])
).resolve()
REV = ROOT / "25_manuscript_package/revision_v2_2026-10-21"
UPLOAD = REV / "upload_set_2026-10-21"
COVER = REV / "cover_letter_MicrobiologySpectrum_v2.md"
MANUSCRIPT = REV / "HPICDA_manuscript_v2.md"
FIGURES = UPLOAD
TABLES = ROOT / "25_manuscript_package/tables/TABLES_INDEX.json"
WORKBOOK = UPLOAD / "Supplementary_Tables_S1-S21b.xlsx"
OUTPUT = UPLOAD / "HPICDA_Spectrum_combined_review.pdf"

pdfmetrics.registerFont(TTFont("Arial", r"C:/Windows/Fonts/arial.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Bold", r"C:/Windows/Fonts/arialbd.ttf"))
pdfmetrics.registerFont(TTFont("Arial-Italic", r"C:/Windows/Fonts/ariali.ttf"))


STYLES = getSampleStyleSheet()
STYLES.add(
    ParagraphStyle(
        name="ArialBody",
        fontName="Arial",
        fontSize=9.2,
        leading=13,
        alignment=TA_JUSTIFY,
        spaceAfter=6,
    )
)
STYLES.add(
    ParagraphStyle(
        name="ArialTitle",
        fontName="Arial-Bold",
        fontSize=16,
        leading=20,
        alignment=TA_CENTER,
        spaceAfter=12,
    )
)
STYLES.add(
    ParagraphStyle(
        name="ArialH1",
        fontName="Arial-Bold",
        fontSize=13,
        leading=16,
        spaceBefore=12,
        spaceAfter=7,
    )
)
STYLES.add(
    ParagraphStyle(
        name="ArialH2",
        fontName="Arial-Bold",
        fontSize=11.2,
        leading=14,
        spaceBefore=10,
        spaceAfter=5,
    )
)
STYLES.add(
    ParagraphStyle(
        name="ArialH3",
        fontName="Arial-Bold",
        fontSize=10,
        leading=13,
        spaceBefore=8,
        spaceAfter=4,
    )
)
STYLES.add(
    ParagraphStyle(
        name="ArialBullet",
        fontName="Arial",
        fontSize=9.1,
        leading=12.5,
        leftIndent=14,
        firstLineIndent=-9,
        spaceAfter=3,
    )
)


def inline_markup(text: str) -> str:
    text = html.escape(text, quote=False)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", text)
    text = text.replace("`", "")
    return text


def markdown_story(path: Path, add_page_break: bool) -> list:
    story: list = []
    paragraph_lines: list[str] = []

    def flush() -> None:
        if paragraph_lines:
            story.append(Paragraph(inline_markup(" ".join(paragraph_lines)), STYLES["ArialBody"]))
            paragraph_lines.clear()

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.rstrip()
        stripped = line.strip()
        if not stripped:
            flush()
            continue
        if stripped in {"---", "***", "___"}:
            flush()
            story.append(Spacer(1, 6))
            continue
        heading = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if heading:
            flush()
            level = len(heading.group(1))
            style = "ArialTitle" if level == 1 else f"ArialH{min(level - 1, 3)}"
            story.append(Paragraph(inline_markup(heading.group(2)), STYLES[style]))
            continue
        bullet = re.match(r"^[-*]\s+(.*)$", raw)
        if bullet:
            flush()
            story.append(Paragraph("• " + inline_markup(bullet.group(1)), STYLES["ArialBullet"]))
            continue
        numbered = re.match(r"^\s*(\d+)\.\s+(.*)$", raw)
        if numbered:
            flush()
            story.append(
                Paragraph(
                    f"{numbered.group(1)}. " + inline_markup(numbered.group(2)),
                    STYLES["ArialBullet"],
                )
            )
            continue
        paragraph_lines.append(stripped)
    flush()
    if add_page_break:
        story.append(PageBreak())
    return story


def table_index_story() -> list:
    entries = json.loads(TABLES.read_text(encoding="utf-8"))
    story = [
        PageBreak(),
        Paragraph("Supplementary table index", STYLES["ArialH1"]),
        Paragraph(
            "The complete manageable table set is attached to this PDF as "
            "Supplementary_Tables_S1-S21b.xlsx. The full S5 call table is archived at "
            "Zenodo DOI 10.5281/zenodo.23550165 because its compressed file exceeds the "
            "15 MB initial-submission limit. S5b is the per-screen summary.",
            STYLES["ArialBody"],
        ),
        Spacer(1, 5),
    ]
    for entry in entries:
        story.append(
            Paragraph(
                f"<b>{html.escape(str(entry['id']))}</b>. "
                f"{inline_markup(str(entry['caption']))} "
                f"({entry['rows']} rows × {entry['cols']} columns)",
                STYLES["ArialBody"],
            )
        )
    return story


def build_text_pdf(path: Path) -> None:
    document = SimpleDocTemplate(
        str(path),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=16 * mm,
        bottomMargin=16 * mm,
        title="HPICDA Microbiology Spectrum submission",
    )
    story = [
        Paragraph("HPICDA — Microbiology Spectrum submission", STYLES["ArialTitle"]),
        Paragraph("Cover letter", STYLES["ArialH1"]),
    ]
    story += markdown_story(COVER, add_page_break=True)
    story += [
        Paragraph("Manuscript", STYLES["ArialH1"]),
    ]
    story += markdown_story(MANUSCRIPT, add_page_break=False)
    story += table_index_story()
    document.build(story)


def combine() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        text_pdf = Path(tmp) / "text.pdf"
        build_text_pdf(text_pdf)
        output = fitz.open()
        text_doc = fitz.open(text_pdf)
        output.insert_pdf(text_doc)

        figure_files = [
            UPLOAD / "fig1_study_architecture.pdf",
            UPLOAD / "fig2_baseline_vs_infection_landscape.pdf",
            UPLOAD / "fig3_cross_pathogen_map.pdf",
            UPLOAD / "fig4_human_context.pdf",
            UPLOAD / "fig5_modules_mechanism.pdf",
            UPLOAD / "fig6_tractability.pdf",
            UPLOAD / "fig7_holdout.pdf",
        ]
        for index, figure in enumerate(figure_files, start=1):
            figure_doc = fitz.open(figure)
            page = output.new_page(width=landscape(A4)[0], height=landscape(A4)[1])
            page.insert_text(
                (28, 28),
                f"Figure {index}",
                fontsize=12,
                fontname="helv",
                color=(0.12, 0.10, 0.08),
            )
            target = fitz.Rect(24, 40, page.rect.width - 24, page.rect.height - 24)
            page.show_pdf_page(target, figure_doc, 0, keep_proportion=True)

        output.set_metadata(
            {
                "title": "HPICDA Microbiology Spectrum submission",
                "subject": "Combined cover letter, manuscript, figures and supplementary table index",
            }
        )
        if WORKBOOK.exists():
            output.embfile_add(
                WORKBOOK.name,
                WORKBOOK.read_bytes(),
                filename=WORKBOOK.name,
            )
        output.save(OUTPUT, garbage=4, deflate=True)
        output.close()
        text_doc.close()
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    combine()
