SCIENTIFIC REPORTS MANUSCRIPT AND ANALYSIS PACKAGE — v1.2.2

AUTHOR
Dengxiang Yang
Southern University of Science and Technology, Shenzhen, China
Correspondence: 12412602@mail.sustech.edu.cn

INITIAL SUBMISSION FILES
main.pdf: upload as the single Manuscript File; six main figures and legends are included.
main_scientific_reports.pdf is an identical, journal-labelled copy; upload only one.
supplementary.pdf: upload separately as Supplementary Information (Figures S1–S8 and Table S1).
cover_letter.pdf: upload as Cover Letter.
The submission-system title, abstract, authors and declarations should match these files.
Originality and exclusive-submission declarations are completed by the author in the submission system.

LATEX SOURCE
main.tex is the manuscript entry point. main_scientific_reports.tex is an alias.
Compile main.tex, supplementary.tex and cover_letter.tex separately with an existing
Tectonic, latexmk or Overleaf installation. The delivered files use Tectonic 0.17.0.
scientific_reports.bst supplies the journal reference-list style; references.bib is the source bibliography.
The main manuscript includes continuous review line numbers and page numbers.

JOURNAL FORMAT
Introduction, Results, Discussion, Methods; Data Availability and Code Availability
before References; Acknowledgements, Funding, Author contributions and Additional
Information after References; main figures and legends at the end. The AI-use
statement appears in Methods. Supplementary captions are in the separate supplement.
Official guide: https://www.nature.com/srep/author-instructions/submission-guidelines
Formatting was checked against the supplied initial-submission checklist.

FIGURES AND SOURCE DATA
figures/: fourteen vector PDF figures, matching SVGs and high-resolution PNGs.
Panel letters use bold lowercase. Scientific labels, numerical values and plot geometry
are preserved. figures_main.tex and supplementary.tex contain the authoritative captions.
SUPPLEMENTARY_FIGURE_MAP.tsv links published figure numbers to preserved source-table prefixes.
source_data/: exact plotted-panel tables and the source for Supplementary Table S1.
source_inputs/ and revision_inputs/: frozen figure-reproduction inputs.
Run code/redraw_figures.py, code/plot_cell_atlases.py, code/plot_sample_robustness.py
and code/plot_annotation_sensitivity.py using code/requirements.txt.

ANALYSIS
../analysis/GSE160639/: raw-count DESeq2 analysis and 27 targeted pathway tests.
../analysis/GSE168041/sensitivity/: 48 animal-omission fits and 56 focused pathway estimates.
../analysis/GSE252285/annotation_sensitivity/: derived annotation and matched-null sensitivity.
../analysis/: adopted primary single-cell workflows, inputs and results.
Scientific data, thresholds, model settings and analysis outputs are unchanged.

PUBLIC ARCHIVE
https://github.com/pdx12320/cochlear-noise-transcriptomics/releases/tag/v1.2.2
Analysis objects and runtime sources are in data collection v1.0.0.
Complete sensitivity records are in data collection v1.1.0.
Both collections provide the data used by the current paper.

VALIDATION
BUILD_VALIDATION.json records document structure, counts, compilation and visual checks.
STATISTICAL_VALIDATION.json records scientific checks.
../manifests/REPOSITORY_FILES.tsv covers the repository payload.
