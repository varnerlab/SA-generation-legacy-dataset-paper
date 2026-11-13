TEX = pdflatex
BIB = bibtex
MAIN = main
SUPP = supplementary
ARXIV_VERSION ?= v3

.PHONY: all clean arxiv-package

# The manuscript and supplement cross-reference one another. Build the manuscript
# once to create main.aux, build the supplement, then rebuild the manuscript so
# every Supplementary Table and Figure reference resolves.
all:
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(BIB) $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(TEX) -interaction=nonstopmode $(SUPP)
	$(TEX) -interaction=nonstopmode $(SUPP)
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN)

# Build the arXiv source ZIP (one main.tex with the supplement appended) in
# ../.zips/, then compile that ZIP alone the way arXiv does and check the result.
arxiv-package:
	python3 scripts/build_arxiv_zip.py $(ARXIV_VERSION)
	python3 scripts/check_arxiv_zip.py $(ARXIV_VERSION)

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log $(MAIN).out \
	      $(MAIN).toc $(MAIN).lof $(MAIN).lot $(MAIN).fls $(MAIN).fdb_latexmk \
	      $(MAIN).synctex.gz $(MAIN).nav $(MAIN).snm $(MAIN).vrb \
	      $(MAIN).run.xml $(MAIN)-blx.bib $(MAIN).pdf \
	      $(SUPP).aux $(SUPP).log $(SUPP).out $(SUPP).pdf
