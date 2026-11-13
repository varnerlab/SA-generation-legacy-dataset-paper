#!/usr/bin/env python3
"""Flag bad page breaks in a PDF: a small piece of a paragraph, or part of a word,
left alone at the top or bottom of a page.

Usage (from arxiv/):  python3 scripts/check_page_breaks.py main.pdf

The check reads each page with pdftotext in layout mode. In that text, the first
line of each paragraph is indented by a few spaces, and the last line of a
paragraph is short and ends with ".", ":", or ")". A page is reported when:
  - it starts with the last one or two lines of a paragraph from the previous
    page (for example, one word or half a sentence);
  - it ends with the first line of a paragraph that continues on the next page;
  - it ends with a hyphenated word part, so one word is split across two pages;
  - it ends inside a web address that continues on the next page.
Caption and table lines are not paragraph text, so lines that start with
"Figure", "Table", "Supplementary", or "Algorithm" are skipped. Reference
pages use a hanging indent (an entry's later lines are indented), so on those
pages the check reports an entry's last line alone at the top of a page, or an
entry's first line alone at the bottom.
"""

from pathlib import Path
import statistics
import subprocess
import sys

CAPTION_STARTS = ("Figure", "Table", "Supplementary", "Algorithm")
PARAGRAPH_INDENT = range(2, 9)   # leading spaces that mark a paragraph's first line


def page_lines(pdf: Path, page: int) -> list[str]:
    """Return the non-empty text lines of one page, with leading spaces kept and
    the page number removed."""
    text = subprocess.run(["pdftotext", "-layout", "-f", str(page), "-l", str(page), str(pdf), "-"],
                          capture_output=True, text=True).stdout
    lines = [line.rstrip() for line in text.splitlines() if line.strip()]
    if lines and lines[-1].strip().isdigit():
        lines = lines[:-1]
    return lines


def indent(line: str) -> int:
    """Return the number of leading spaces in a line."""
    return len(line) - len(line.lstrip(" "))


def is_paragraph_start(line: str) -> bool:
    """Return True when the line is the indented first line of a paragraph."""
    return indent(line) in PARAGRAPH_INDENT


def is_paragraph_end(line: str, typical: float) -> bool:
    """Return True when the line is the short last line of a paragraph."""
    text = line.strip()
    return text.endswith((".", ":", ")")) and len(text) < 0.80 * typical


def is_reference_line(line: str) -> bool:
    """Return True when the line is a later line of a reference entry (hanging indent).
    pdftotext indents these lines by one or two spaces."""
    return 1 <= indent(line) <= 8


def is_reference_page(lines: list[str]) -> bool:
    """Return True when many lines are indented, as in a hanging-indent reference list."""
    indented = sum(1 for line in lines if is_reference_line(line))
    return indented > 0.4 * len(lines)


def is_body_text(line: str) -> bool:
    """Return True for paragraph text, and False for captions and table lines."""
    return not line.strip().startswith(CAPTION_STARTS)


def bad_breaks(pdf: Path) -> list[str]:
    """Return one message for each bad page break in the PDF."""
    info = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    pages = int(next(line.split()[-1] for line in info.splitlines() if line.startswith("Pages:")))
    messages = []
    previous = page_lines(pdf, 1)
    for page in range(2, pages + 1):
        lines = page_lines(pdf, page)
        if not lines or not previous:
            previous = lines
            continue
        typical = statistics.median(len(line.strip()) for line in lines)
        first = lines[0]
        last = previous[-1]

        if is_reference_page(lines) and is_reference_page(previous):
            # A reference entry's last line alone at the top, or its first line alone at the bottom.
            if is_reference_line(first) and (len(lines) == 1 or indent(lines[1]) == 0):
                messages.append(f"page {page} starts with the last line of a reference: {first.strip()[:40]!r}")
            if indent(last) == 0 and is_reference_line(first):
                messages.append(f"page {page - 1} ends with the first line of a reference: {last.strip()[:40]!r}")
            previous = lines
            continue

        # The page starts with the last one or two lines of a paragraph.
        if is_body_text(first) and indent(first) == 0:
            if is_paragraph_end(first, typical):
                messages.append(f"page {page} starts with the last line of a paragraph: {first.strip()!r}")
            elif len(lines) > 1 and indent(lines[1]) == 0 and is_paragraph_end(lines[1], typical):
                messages.append(f"page {page} starts with the last two lines of a paragraph: "
                                f"{lines[1].strip()!r}")

        # The previous page ends with the first line of a paragraph.
        if is_body_text(last) and is_paragraph_start(last) and indent(first) == 0 \
                and is_body_text(first):
            messages.append(f"page {page - 1} ends with the first line of a paragraph: {last.strip()[:40]!r}")

        # A word or a web address is split across the two pages.
        last_text = last.strip()
        if len(last_text) > 1 and last_text.endswith("-") and last_text[-2].isalpha():
            messages.append(f"page {page - 1} ends with a split word: {last_text[-30:]!r}")
        last_token = last_text.split()[-1]
        if ("://" in last_token or "github.com" in last_token) and first.strip()[0].islower():
            messages.append(f"page {page - 1} ends inside a web address: {last_token!r}")
        previous = lines
    return messages


def main() -> None:
    messages = bad_breaks(Path(sys.argv[1]))
    for message in messages:
        print("FAIL:", message)
    if messages:
        sys.exit(1)
    print(f"PASS: no paragraph fragments or split words across pages in {sys.argv[1]}")


if __name__ == "__main__":
    main()
