# Third-party notices

## Public source studies

Raw sequencing and deposited matrices originate from GEO GSE168041, GSE160639 and GSE252285. Cite their source records and linked publications; use the public accessions listed in `data/README.md`. This archive contains derived analysis results and selected processed objects, not newly generated animal data.

## MSigDB

Frozen gene sets use MSigDB 2026.1.Mm. MSigDB contents are copyright Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California, distributed under CC BY 4.0 with collection-specific additional terms: https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp . The archived definitions used here are mouse Hallmark, Reactome and GO Biological Process sets. They retain exact membership or are filtered to the declared analysis targets; no KEGG/BioCarta/STKE sets are used in the manuscript.

Attribute Castanza et al. (2023), doi:10.1038/s41592-023-02014-7; Liberzon et al. (2015), doi:10.1016/j.cels.2015.12.004; and the relevant Reactome/GO resources. Membership filtering and rank-universe intersections are recorded in the analysis outputs.

## Reactome

Reactome data are available under CC0 (https://reactome.org/license). Software and illustrations may have different licenses. The manuscript cites the Reactome knowledgebase paper; this archive does not redistribute Reactome illustrations.

## Software

R, Python, Seurat, UCell, Scanpy, PyDESeq2, DESeq2, fgsea, scDblFinder and other dependencies retain their own licenses. Their installed environments are documented in version and source manifests. The release asset `Bulk_R_Source_Packages_20261101.tar` preserves original R package source distributions and their included licenses. Those licenses continue to govern the third-party package sources; the repository MIT license covers author-owned code. Conda runtime binaries are obtained from the exact upstream URLs in the environment lock.

## Manuscript

Manuscript text: copyright 2026 Dengxiang Yang, all rights reserved except permissions applicable to individual cited or externally supplied materials. The repository code license does not apply to the manuscript text.

## Bibliography style

`manuscript/scientific_reports.bst` is a distinctly named derivative of Patrick W. Daly's `unsrtnat.bst`, distributed under the LaTeX Project Public License. The original copyright and license notice remain in the style file.
