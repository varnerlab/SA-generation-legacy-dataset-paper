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
SUBMISSION_DIR = submission-files
ARCHIVE = ../Varner-paper-technical-correction-v2.1.zip

.PHONY: all clean figures-upload response submission-files package

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
	$(MAKE) figures-upload response submission-files

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

submission-files:
	mkdir -p $(SUBMISSION_DIR)/manuscript-file \
	         $(SUBMISSION_DIR)/related-files \
	         $(SUBMISSION_DIR)/supplementary-information \
	         $(SUBMISSION_DIR)/figures \
	         $(SUBMISSION_DIR)/rebuttal-letter
	cp $(CLEAN).pdf $(SUBMISSION_DIR)/manuscript-file/
	cp $(MARKED).pdf $(SUBMISSION_DIR)/related-files/
	cp $(SUPP).pdf $(SUBMISSION_DIR)/supplementary-information/
	cp figures-upload/Figure*.pdf $(SUBMISSION_DIR)/figures/
	cp $(RESPONSE_PDF) $(SUBMISSION_DIR)/rebuttal-letter/
	cp README.txt $(SUBMISSION_DIR)/related-files/
	cp upload-checklist.md $(SUBMISSION_DIR)/

package: all
	zip -rq -FS $(ARCHIVE) $(SUBMISSION_DIR)

clean:
	rm -f $(MARKED).aux $(MARKED).bbl $(MARKED).blg $(MARKED).log $(MARKED).out \
	      $(MARKED).pdf $(CLEAN).aux $(CLEAN).bbl $(CLEAN).blg $(CLEAN).log \
	      $(CLEAN).out $(CLEAN).pdf $(SUPP).aux $(SUPP).log $(SUPP).out $(SUPP).pdf
	$(MAKE) -C $(RESPONSE_DIR) clean
