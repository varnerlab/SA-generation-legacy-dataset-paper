#!/usr/bin/env python3
"""Build the clean, compiler-ready LaTeX source archive for SNAPP."""

from pathlib import Path
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


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required source file is missing: {path}")
    return path


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

The compiled bibliography (main.bbl) and the LaTeX-generated supplementary
label data (supplementary-information.aux) are supplied. Supplementary citation
numbers remain generated from the labels in the standalone supplement; they
are not entered manually in the manuscript source.
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
        "supplementary-information.aux": require_file(
            REVISION_DIR / "supplementary-information.aux"
        ).read_bytes(),
        "README.txt": readme.encode(),
    }

    for name in SECTION_NAMES:
        source = require_file(REVISION_DIR / "sections" / f"{name}.tex")
        members[f"sections/{name}.tex"] = source.read_bytes()

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
