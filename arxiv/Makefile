TEX = pdflatex
BIB = bibtex
MAIN = main
ARXIV_VERSION ?= v3

.PHONY: all clean arxiv-package

# The preprint is one document: main.tex includes the supplement after the main
# figures. Three LaTeX passes after BibTeX resolve every cross-reference.
all:
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(BIB) $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN)

# Build the arXiv source ZIP in ../.zips/, then compile that ZIP alone the way
# arXiv does and check the result.
arxiv-package:
	python3 scripts/build_arxiv_zip.py $(ARXIV_VERSION)
	python3 scripts/check_arxiv_zip.py $(ARXIV_VERSION)

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log $(MAIN).out \
	      $(MAIN).toc $(MAIN).lof $(MAIN).lot $(MAIN).fls $(MAIN).fdb_latexmk \
	      $(MAIN).synctex.gz $(MAIN).nav $(MAIN).snm $(MAIN).vrb \
	      $(MAIN).run.xml $(MAIN)-blx.bib $(MAIN).pdf
