#!/usr/bin/env python3
"""Check an arXiv source ZIP by compiling it alone, the way arXiv does.

Usage (from arxiv/):  python3 scripts/check_arxiv_zip.py v3

The ZIP is unpacked into an empty temporary folder and built with pdflatex
three times using the shipped main.bbl (arXiv does not run BibTeX). The check
fails if the build needs a file outside the ZIP, leaves a reference undefined,
or prints "??", or if the supplement is missing from the PDF.
"""

from pathlib import Path
import re
import subprocess
import sys
import tempfile
import zipfile

REPOSITORY_DIR = Path(__file__).resolve().parents[2]


def main() -> None:
    version = sys.argv[1]
    archive = REPOSITORY_DIR / ".zips" / f"arxiv-paper-{version}.zip"
    if not archive.is_file():
        sys.exit(f"FAIL: {archive} does not exist")

    failures = []
    with tempfile.TemporaryDirectory() as folder:
        build = Path(folder)
        with zipfile.ZipFile(archive) as source:
            source.extractall(build)

        # Exactly one top-level LaTeX file, so arXiv compiles one document.
        top_level = [p.name for p in build.glob("*.tex") if "\\documentclass" in p.read_text()]
        if top_level != ["main.tex"]:
            failures.append(f"top-level LaTeX files are {top_level}, expected ['main.tex']")
        if not (build / "main.bbl").is_file():
            failures.append("main.bbl is missing")

        # The sources must not depend on another document or carry reviewer markup.
        sources = "\n".join(p.read_text() for p in build.rglob("*.tex"))
        for token in ("\\externaldocument", "\\rone{", "\\rtwo{", "\\rboth{", "\\rthree{"):
            if token in sources:
                failures.append(f"sources contain {token}")

        for _ in range(3):
            subprocess.run(["pdflatex", "-interaction=nonstopmode", "main.tex"],
                           cwd=build, capture_output=True)
        log = (build / "main.log").read_text(errors="replace")
        pdf = build / "main.pdf"
        if not pdf.is_file():
            failures.append("pdflatex did not produce main.pdf")
        else:
            text = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True,
                                  text=True).stdout
            if "??" in text:
                failures.append("the PDF contains ??")
            for needed in ("Supplementary Information",
                           "Supplementary Bandwidth and Plausibility Results",
                           "Supplementary Table 18: Two ways to keep hormone values positive",
                           "de-identified real patient data"):
                if needed not in re.sub(r"\s+", " ", text):
                    failures.append(f"the PDF does not contain {needed!r}")
        for pattern in (r"undefined", r"multiply defined", r"File `[^']+' not found",
                        r"^! "):
            if re.search(pattern, log, re.MULTILINE):
                failures.append(f"main.log matches {pattern!r}")

    if failures:
        for failure in failures:
            print("FAIL:", failure)
        sys.exit(1)
    print(f"PASS: {archive.name} compiles alone with the supplement and no unresolved references")


if __name__ == "__main__":
    main()
