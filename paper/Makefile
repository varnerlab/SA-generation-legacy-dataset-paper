TEX = pdflatex
BIB = bibtex
MAIN = main
SUPP = supplementary

.PHONY: all clean

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

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log $(MAIN).out \
	      $(MAIN).toc $(MAIN).lof $(MAIN).lot $(MAIN).fls $(MAIN).fdb_latexmk \
	      $(MAIN).synctex.gz $(MAIN).nav $(MAIN).snm $(MAIN).vrb \
	      $(MAIN).run.xml $(MAIN)-blx.bib $(MAIN).pdf \
	      $(SUPP).aux $(SUPP).log $(SUPP).out $(SUPP).pdf
