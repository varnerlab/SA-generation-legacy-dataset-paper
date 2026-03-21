TEX = pdflatex
BIB = bibtex
MAIN = main

.PHONY: all clean

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex references.bib sections/*.tex
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(BIB) $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN)
	$(TEX) -interaction=nonstopmode $(MAIN)

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log $(MAIN).out \
	      $(MAIN).toc $(MAIN).lof $(MAIN).lot $(MAIN).fls $(MAIN).fdb_latexmk \
	      $(MAIN).synctex.gz $(MAIN).nav $(MAIN).snm $(MAIN).vrb \
	      $(MAIN).run.xml $(MAIN)-blx.bib $(MAIN).pdf
