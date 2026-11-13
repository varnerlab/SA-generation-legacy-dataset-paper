LATEXMK ?= latexmk

SOURCE := response-to-reviewers.tex
PDF := response-to-reviewers.pdf
SOURCE2 := response-to-reviewers-round2.tex
PDF2 := response-to-reviewers-round2.pdf

LATEXMK_FLAGS := \
	-xelatex \
	-interaction=nonstopmode \
	-halt-on-error \
	-file-line-error

.PHONY: all pdf round2 clean rebuild

all: $(PDF) $(PDF2)

pdf: $(PDF)

round2: $(PDF2)

$(PDF): $(SOURCE)
	$(LATEXMK) $(LATEXMK_FLAGS) $(SOURCE)

$(PDF2): $(SOURCE2)
	$(LATEXMK) $(LATEXMK_FLAGS) $(SOURCE2)

rebuild:
	$(MAKE) clean
	$(MAKE) all

clean:
	$(LATEXMK) -c $(SOURCE) $(SOURCE2)
