PANDOC ?= pandoc
LATEXMK ?= latexmk

SOURCE := response-to-reviewers.md
TEX := response-to-reviewers.tex
PDF := response-to-reviewers.pdf

PANDOC_FLAGS := \
	--from=gfm+hard_line_breaks+tex_math_dollars \
	--to=latex \
	--standalone \
	--wrap=preserve \
	--variable=documentclass=article \
	--variable=fontsize=11pt \
	--variable=papersize=letter \
	--variable=geometry=margin=0.5in \
	--variable=mainfont=Arial \
	--variable=colorlinks=true \
	--variable=linkcolor=black \
	--variable=filecolor=black \
	--variable=citecolor=black \
	--variable=urlcolor=black \
	--metadata=pagetitle:'Response to Reviewers'

LATEXMK_FLAGS := \
	-xelatex \
	-interaction=nonstopmode \
	-halt-on-error \
	-file-line-error

.PHONY: all tex pdf clean rebuild

all: $(PDF)

tex: $(TEX)

pdf: $(PDF)

$(TEX): $(SOURCE) Makefile
	$(PANDOC) $(SOURCE) $(PANDOC_FLAGS) --output=$(TEX)

$(PDF): $(TEX)
	$(LATEXMK) $(LATEXMK_FLAGS) $(TEX)

rebuild:
	$(MAKE) clean
	$(MAKE) all

clean:
	$(LATEXMK) -c $(TEX)
