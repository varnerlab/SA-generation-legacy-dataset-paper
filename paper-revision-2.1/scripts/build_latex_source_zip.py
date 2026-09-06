#!/usr/bin/env python3
"""Build the clean, compiler-ready LaTeX source archive for SNAPP."""

from pathlib import Path
import re
import shutil
import subprocess
import zipfile


REVISION_DIR = Path(__file__).resolve().parent.parent
REPOSITORY_DIR = REVISION_DIR.parent
ARCHIVE = REPOSITORY_DIR / "Varner-LaTeX-Source-Clean-v2.1.zip"

SECTION_NAMES = (
    "abstract",
    "introduction",
    "results",
    "discussion",
    "methods",
    "acknowledgments",
    "floats",
)

FIGURE_NAMES = (
    "validate_bio_correlations.pdf",
    "validate_pregnancy_progression.pdf",
    "validate_cross_visit_corr.pdf",
    "sim_phase_diagram.pdf",
    "validate_conditioned_features_v2.pdf",
    "validate_mechanistic_combined_tf_only_v2.pdf",
    "downstream_utility_scatter_v2.pdf",
)

SUPPLEMENTARY_LABEL_PATTERN = re.compile(
    r"^\\newlabel\{(?P<label>[^}]+)\}\{\{(?P<number>[^{}]*)\}\{",
    re.MULTILINE,
)

SUPPLEMENTARY_REFERENCE_PATTERNS = (
    re.compile(r"\\supptable\{([^}]+)\}"),
    re.compile(r"\\suppfigure\{([^}]+)\}"),
    re.compile(r"\\suppfigures\{([^}]+)\}\{([^}]+)\}"),
)


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required source file is missing: {path}")
    return path


def build_supplementary_label_source(
    aux_text: str, manuscript_text: str
) -> bytes:
    """Convert LaTeX-generated external labels into SNAPP-safe TeX source."""
    label_numbers = {
        match.group("label"): match.group("number")
        for match in SUPPLEMENTARY_LABEL_PATTERN.finditer(aux_text)
    }

    referenced_labels: set[str] = set()
    for pattern in SUPPLEMENTARY_REFERENCE_PATTERNS:
        for match in pattern.finditer(manuscript_text):
            referenced_labels.update(match.groups())

    missing_labels = referenced_labels - label_numbers.keys()
    if missing_labels:
        missing = ", ".join(sorted(missing_labels))
        raise RuntimeError(f"Supplementary labels missing from compiled aux: {missing}")

    lines = [
        "% Generated automatically from supplementary-information.aux.",
        "% Do not enter or edit supplementary numbers manually.",
    ]
    for label in sorted(referenced_labels):
        number = label_numbers[label]
        lines.append("\\newlabel{" + label + "}{{" + number + "}{}{}{}{}}")
    lines.append("")
    return "\n".join(lines).encode()


def main() -> None:
    kpsewhich = shutil.which("kpsewhich")
    if kpsewhich is None:
        raise RuntimeError("kpsewhich is required to locate naturemag.bst")

    bst_output = subprocess.check_output(
        [kpsewhich, "naturemag.bst"], text=True
    ).strip()
    nature_bst = require_file(Path(bst_output))

    manuscript = require_file(REVISION_DIR / "main.tex").read_text()
    if not manuscript.startswith("\\documentclass"):
        raise RuntimeError("Expected main.tex to begin with \\documentclass")

    external_document = "\\externaldocument{supplementary-information}"
    if manuscript.count(external_document) != 1:
        raise RuntimeError(
            "Expected exactly one supplementary external-document declaration"
        )
    manuscript = manuscript.replace(
        external_document, "\\input{supplementary-labels.tex}"
    )

    supplementary_reference_macros = (
        (
            "\\newcommand{\\supptable}[1]{Supplementary Table~\\ref{#1}}",
            "\\newcommand{\\supptable}[1]{Supplementary Table~\\ref*{#1}}",
        ),
        (
            "\\newcommand{\\suppfigure}[1]{Supplementary Figure~\\ref{#1}}",
            "\\newcommand{\\suppfigure}[1]{Supplementary Figure~\\ref*{#1}}",
        ),
        (
            "\\newcommand{\\suppfigures}[2]{Supplementary Figures~\\ref{#1}--\\ref{#2}}",
            "\\newcommand{\\suppfigures}[2]{Supplementary Figures~\\ref*{#1}--\\ref*{#2}}",
        ),
    )
    for linked_macro, unlinked_macro in supplementary_reference_macros:
        if manuscript.count(linked_macro) != 1:
            raise RuntimeError(
                "Expected exactly one supplementary reference macro definition"
            )
        manuscript = manuscript.replace(linked_macro, unlinked_macro)

    section_sources = {
        f"sections/{name}.tex": require_file(
            REVISION_DIR / "sections" / f"{name}.tex"
        ).read_bytes()
        for name in SECTION_NAMES
    }
    manuscript_sections = "\n".join(
        source.decode() for source in section_sources.values()
    )
    supplementary_labels = build_supplementary_label_source(
        require_file(
            REVISION_DIR / "supplementary-information.aux"
        ).read_text(),
        manuscript_sections,
    )

    clean_main = (
        "% !TeX program = pdflatex\n"
        "% !TeX root = main.tex\n"
        "% This is the single main LaTeX file for the SNAPP source archive.\n"
        "\\def\\cleanversion{1}\n"
        + manuscript
    )

    readme = """SPRINGER NATURE SNAPP LATEX SOURCE PACKAGE

MAIN LATEX FILE: main.tex
COMPILER: pdflatex

The archive contains exactly one LaTeX file with a document class and document
environment: main.tex. It builds the clean, unmarked manuscript.

For a manual build, run pdflatex twice:

  pdflatex -interaction=nonstopmode main.tex
  pdflatex -interaction=nonstopmode main.tex

The compiled bibliography (main.bbl) is supplied. The supplementary label map
(supplementary-labels.tex) is generated automatically from the current
standalone supplement. The manuscript continues to cite symbolic label names;
supplementary numbers are not entered manually. The generated TeX map avoids
depending on an auxiliary .aux file that an online compiler may discard.
"""

    members: dict[str, bytes] = {
        "main.tex": clean_main.encode(),
        "main.bbl": require_file(
            REVISION_DIR / "manuscript-clean.bbl"
        ).read_bytes(),
        "references.bib": require_file(
            REVISION_DIR / "references.bib"
        ).read_bytes(),
        "naturemag.bst": nature_bst.read_bytes(),
        "supplementary-labels.tex": supplementary_labels,
        "README.txt": readme.encode(),
    }

    members.update(section_sources)

    for name in FIGURE_NAMES:
        source = require_file(REVISION_DIR / "sections" / "figures" / name)
        members[f"sections/figures/{name}"] = source.read_bytes()

    temporary_archive = ARCHIVE.with_suffix(".zip.tmp")
    with zipfile.ZipFile(
        temporary_archive, "w", compression=zipfile.ZIP_DEFLATED
    ) as archive:
        for archive_name in sorted(members):
            archive.writestr(archive_name, members[archive_name])

    temporary_archive.replace(ARCHIVE)
    print(f"Created {ARCHIVE} with {len(members)} files")


if __name__ == "__main__":
    main()
