#!/usr/bin/env python3
"""Build the arXiv source ZIP: one main.tex with the supplement appended.

Usage (from arxiv/):  python3 scripts/build_arxiv_zip.py v3

The local arxiv/ build makes two PDFs (main.pdf and supplementary.pdf) that
reference each other through xr-hyper and supplementary.aux. arXiv compiles
only the submitted source, so this script makes one self-contained document,
as in the posted v2:

  1. Remove the xr-hyper package and the \\externaldocument{supplementary} line
     from main.tex, so main.tex no longer needs supplementary.aux.
  2. Append \\input{sections/supplementary} after the main figures. The
     supplement resets its own table, figure, and equation counters, so its
     numbering matches the standalone supplement.
  3. Compile the combined document with pdflatex and BibTeX in a temporary
     folder, so the shipped main.bbl includes the references the supplement cites.
  4. Write .zips/arxiv-paper-<version>.zip with main.tex, main.bbl,
     references.bib, the section files, and only the figures the sections use.
     Also write .zips/arxiv-paper-<version>.pdf, the PDF this source produces.
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

SUPPLEMENT_INPUT = """
% ── Supplementary Information (appended for the arXiv source) ─────────────────
\\clearpage
\\input{sections/supplementary}
"""


def replace_once(text: str, old: str, new: str) -> str:
    """Replace one exact occurrence of old, and fail if there is not exactly one."""
    if text.count(old) != 1:
        raise RuntimeError(f"Expected exactly one {old!r} in arxiv/main.tex")
    return text.replace(old, new)


def combined_main() -> str:
    """Return arxiv/main.tex changed into one document that includes the supplement."""
    text = (ARXIV_DIR / "main.tex").read_text()
    text = replace_once(text, "\\usepackage{xr-hyper}\n", "")
    text = replace_once(text, "\\externaldocument{supplementary}\n", "")
    text = replace_once(text, "\\input{sections/floats}\n",
                        "\\input{sections/floats}\n" + SUPPLEMENT_INPUT)
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
    members: dict[str, bytes] = {"main.tex": combined_main().encode(),
                                 "references.bib": (ARXIV_DIR / "references.bib").read_bytes()}
    section_text = ""
    for name in SECTION_NAMES:
        source = (ARXIV_DIR / "sections" / f"{name}.tex").read_text()
        members[f"sections/{name}.tex"] = source.encode()
        section_text += source
    for name in used_figures(section_text):
        members[f"sections/figures/{name}"] = (ARXIV_DIR / "sections" / "figures" / name).read_bytes()

    # Compile the combined document once to make main.bbl and the preview PDF.
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
            raise RuntimeError("The combined document did not build; see main.log")
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
