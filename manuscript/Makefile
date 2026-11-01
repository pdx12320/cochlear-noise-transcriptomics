TECTONIC ?= tectonic
PYTHON ?= python3

all: main.pdf main_scientific_reports.pdf supplementary.pdf cover_letter.pdf

main.pdf: main.tex preamble.tex preamble_srep.tex scientific_reports.bst references.bib $(wildcard sections/*.tex) figures_main.tex $(wildcard figures/Figure[1-6]_*.pdf)
	$(TECTONIC) main.tex

main_scientific_reports.pdf: main.pdf
	cp main.pdf main_scientific_reports.pdf

supplementary.pdf: supplementary.tex preamble.tex $(wildcard figures/FigureS*.pdf)
	$(TECTONIC) supplementary.tex

cover_letter.pdf: cover_letter.tex preamble.tex
	$(TECTONIC) cover_letter.tex

figures:
	$(PYTHON) code/redraw_figures.py
	$(PYTHON) code/plot_cell_atlases.py
	$(PYTHON) code/plot_sample_robustness.py
	$(PYTHON) code/plot_annotation_sensitivity.py

.PHONY: all figures
