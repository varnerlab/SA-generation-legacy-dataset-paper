TEX = pdflatex
BIB = bibtex
MAIN = main
MARKED = manuscript-marked
CLEAN = manuscript-clean
SUPP = supplementary-information
RESPONSE_DIR = ../peer-review-feedback
RELATED_DIR = related-files
RESPONSE_TEX = $(RELATED_DIR)/response-to-reviewers.tex
RESPONSE_PDF = $(RELATED_DIR)/response-to-reviewers.pdf

.PHONY: all clean figures-upload response

all:
	$(TEX) -jobname=$(MARKED) -interaction=nonstopmode '\input{$(MAIN).tex}'
	$(BIB) $(MARKED)
	$(TEX) -jobname=$(MARKED) -interaction=nonstopmode '\input{$(MAIN).tex}'
	$(TEX) -jobname=$(CLEAN) -interaction=nonstopmode '\def\cleanversion{1}\input{$(MAIN).tex}'
	$(BIB) $(CLEAN)
	$(TEX) -jobname=$(CLEAN) -interaction=nonstopmode '\def\cleanversion{1}\input{$(MAIN).tex}'
	$(TEX) -jobname=$(SUPP) -interaction=nonstopmode '\def\cleanversion{1}\input{supplementary.tex}'
	$(TEX) -jobname=$(SUPP) -interaction=nonstopmode '\def\cleanversion{1}\input{supplementary.tex}'
	$(TEX) -jobname=$(MARKED) -interaction=nonstopmode '\input{$(MAIN).tex}'
	$(TEX) -jobname=$(MARKED) -interaction=nonstopmode '\input{$(MAIN).tex}'
	$(TEX) -jobname=$(CLEAN) -interaction=nonstopmode '\def\cleanversion{1}\input{$(MAIN).tex}'
	$(TEX) -jobname=$(CLEAN) -interaction=nonstopmode '\def\cleanversion{1}\input{$(MAIN).tex}'
	$(MAKE) figures-upload response

figures-upload:
	mkdir -p figures-upload
	cp sections/figures/validate_bio_correlations.pdf figures-upload/Figure1.pdf
	cp sections/figures/validate_pregnancy_progression.pdf figures-upload/Figure2.pdf
	cp sections/figures/validate_cross_visit_corr.pdf figures-upload/Figure3.pdf
	cp sections/figures/sim_phase_diagram.pdf figures-upload/Figure4.pdf
	cp sections/figures/validate_conditioned_features_v2.pdf figures-upload/Figure5.pdf
	cp sections/figures/validate_mechanistic_combined_tf_only_v2.pdf figures-upload/Figure6.pdf
	cp sections/figures/downstream_utility_scatter_v2.pdf figures-upload/Figure7.pdf

response:
	$(MAKE) -C $(RESPONSE_DIR) all
	mkdir -p $(RELATED_DIR)
	cp $(RESPONSE_DIR)/response-to-reviewers.tex $(RESPONSE_TEX)
	cp $(RESPONSE_DIR)/response-to-reviewers.pdf $(RESPONSE_PDF)

clean:
	rm -f $(MARKED).aux $(MARKED).bbl $(MARKED).blg $(MARKED).log $(MARKED).out \
	      $(MARKED).pdf $(CLEAN).aux $(CLEAN).bbl $(CLEAN).blg $(CLEAN).log \
	      $(CLEAN).out $(CLEAN).pdf $(SUPP).aux $(SUPP).log $(SUPP).out $(SUPP).pdf
	$(MAKE) -C $(RESPONSE_DIR) clean
