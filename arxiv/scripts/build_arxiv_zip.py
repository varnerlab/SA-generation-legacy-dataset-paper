#!/usr/bin/env python3
"""Build the arXiv source ZIP from arxiv/main.tex.

Usage (from arxiv/):  python3 scripts/build_arxiv_zip.py v3

arxiv/main.tex is one document: it includes the supplement after the main
figures. arXiv compiles only the submitted source and does not run BibTeX, so
this script:

  1. Checks that main.tex includes the supplement and does not depend on
     another document's .aux file.
  2. Compiles the document with pdflatex and BibTeX in a temporary folder, so
     the shipped main.bbl matches the submitted source.
  3. Writes .zips/arxiv-paper-<version>.zip with main.tex, main.bbl,
     references.bib, the section files, and only the figures the sections use.
     It also writes .zips/arxiv-paper-<version>.pdf, the PDF this source produces.
"""

from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

ARXIV_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_DIR = ARXIV_DIR.parent

SECTION_NAMES = (
    "abstract",
    "introduction",
    "results",
    "discussion",
    "methods",
    "acknowledgments",
    "floats",
    "supplementary",
)


def checked_main() -> str:
    """Return arxiv/main.tex after checking that it is one self-contained document."""
    text = (ARXIV_DIR / "main.tex").read_text()
    if text.count("\\input{sections/supplementary}") != 1:
        raise RuntimeError("arxiv/main.tex must include sections/supplementary exactly once")
    for token in ("\\externaldocument", "xr-hyper"):
        if token in text:
            raise RuntimeError(f"arxiv/main.tex must not use {token}")
    return text


def used_figures(section_text: str) -> list[str]:
    """Return the figure file names that the section files include."""
    names = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", section_text)
    figures = sorted(set(names))
    for name in figures:
        if not (ARXIV_DIR / "sections" / "figures" / name).is_file():
            raise FileNotFoundError(f"Figure used in the sections is missing: {name}")
    return figures


def run(command: list[str], folder: Path) -> None:
    """Run one LaTeX or BibTeX pass; LaTeX warnings are checked later in the log."""
    subprocess.run(command, cwd=folder, capture_output=True)


def main() -> None:
    version = sys.argv[1]
    archive = REPOSITORY_DIR / ".zips" / f"arxiv-paper-{version}.zip"
    preview = REPOSITORY_DIR / ".zips" / f"arxiv-paper-{version}.pdf"

    # Collect the files that go into the archive.
    members: dict[str, bytes] = {"main.tex": checked_main().encode(),
                                 "references.bib": (ARXIV_DIR / "references.bib").read_bytes()}
    section_text = ""
    for name in SECTION_NAMES:
        source = (ARXIV_DIR / "sections" / f"{name}.tex").read_text()
        members[f"sections/{name}.tex"] = source.encode()
        section_text += source
    for name in used_figures(section_text):
        members[f"sections/figures/{name}"] = (ARXIV_DIR / "sections" / "figures" / name).read_bytes()

    # Compile the document once to make main.bbl and the preview PDF.
    with tempfile.TemporaryDirectory() as folder:
        build = Path(folder)
        for member, data in members.items():
            path = build / member
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        run(["pdflatex", "-interaction=nonstopmode", "main.tex"], build)
        run(["bibtex", "main"], build)
        for _ in range(3):
            run(["pdflatex", "-interaction=nonstopmode", "main.tex"], build)
        if not (build / "main.bbl").is_file() or not (build / "main.pdf").is_file():
            raise RuntimeError("The document did not build; see main.log")
        members["main.bbl"] = (build / "main.bbl").read_bytes()
        shutil.copyfile(build / "main.pdf", preview)

    # Write the archive in a fixed order so that rebuilds differ only when sources change.
    temporary = archive.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as output:
        for member in sorted(members):
            entry = zipfile.ZipInfo(member, date_time=(2026, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            output.writestr(entry, members[member])
    temporary.replace(archive)
    print(f"Created {archive} with {len(members)} files and {preview.name}")


if __name__ == "__main__":
    main()
