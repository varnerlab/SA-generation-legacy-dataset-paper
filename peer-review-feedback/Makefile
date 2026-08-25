LATEXMK ?= latexmk

SOURCE := response-to-reviewers.tex
PDF := response-to-reviewers.pdf

LATEXMK_FLAGS := \
	-xelatex \
	-interaction=nonstopmode \
	-halt-on-error \
	-file-line-error

.PHONY: all pdf clean rebuild

all: $(PDF)

pdf: $(PDF)

$(PDF): $(SOURCE)
	$(LATEXMK) $(LATEXMK_FLAGS) $(SOURCE)

rebuild:
	$(MAKE) clean
	$(MAKE) all

clean:
	$(LATEXMK) -c $(SOURCE)
