# Sample-robustness sensitivity analysis planned before execution

## Scope and inferential status

Retain the complete eight-animal main analysis for every tissue–cell-type model. Evaluate six models: fibrocyte/mesenchymal, endothelial and immune/myeloid pseudobulks in the spiral ganglion (SGN) and lateral wall. Fit eight leave-one-animal-out (LOO) models per group, each retaining three animals from one condition and four from the other. These 48 overlapping fits describe sensitivity to individual animals; they are not independent replications, new confirmatory tests, or estimates of statistical power. Never exclude an animal from the main analysis based on these results.

Primary reporting covers the TNF–NF-kB Hallmark pathway in all six models and the Reactome MAPK-family pathway in lateral-wall endothelial cells. Report every omitted sample_id/GSM/condition, retained group sizes, NES, nominal P, FDR, direction, leading edge, and estimation warnings. Show original full-sample estimates alongside all LOO estimates. Report direction retention and the range of NES/FDR without defining a new pass/fail robustness threshold. A nonsignificant LOO estimate alone is not evidence of an opposite effect.

## Locked inputs and fitting

Read existing frozen pseudobulk integer counts, sample metadata, gene annotation, full-analysis DE tables and fitted-object normalization records. Preserve each full model's fitted gene IDs (the original total-count-at-least-10 mask applied across all eight samples). Do not reselect genes after deleting an animal, change cell annotations, recompute pseudobulks, or reapply cell/sample coverage gates. Any gene becoming all zero or acquiring an undefined statistic is retained in the complete DE output and explicitly recorded.

Use the existing Python 3.12.13/PyDESeq2 0.5.4 environment and original `~condition` design, `Noise` versus `Naive` contrast, alpha 0.05, Cook's filtering/refitting and independent filtering. Re-estimate normalization and dispersion on each seven-animal dataset. Retain the original PyDESeq2 constructor/default inference settings and seed 20260825. Use at most eight inference workers as originally specified, with BLAS/OpenMP threads limited to one. Gene-level FDR remains model-specific.

Use existing R 4.5.3/fgsea 1.36.2, frozen native mouse MSigDB 2026.1.Mm Hallmark plus Reactome memberships, signed Wald ranking, and the original duplicate rule (largest absolute statistic, then signed statistic, then lexical gene ID). Preserve minSize 10, maxSize 500, weighting exponent 1, eps 0, standard two-sided score, multilevel sample size 101, 1,000 preliminary permutations, one fgsea process, and seed 20260825 reset per fit. The original full-sample eligible testing family is locked separately for each model: SGN endothelial 926, fibrocyte/mesenchymal 984, immune/myeloid 954; lateral-wall endothelial 929, fibrocyte/mesenchymal 954, immune/myeloid 898. All family members are evaluated, not just the displayed pathways. If an original member becomes untestable under the unchanged overlap limits, report it as unavailable and retain the original family size in BH adjustment; do not substitute a newly eligible pathway. No FDR family pooling across LOO iterations is used or claimed.

## Per-animal descriptive displays

Use the original full-fit normalized pseudobulk counts for animal-level displays, checking raw counts divided by archived size factors against the stored layer. Export representative genes already discussed in the manuscript (`Ccl5`, `Ifih1`, `Mcl1`, `Ifit2`, `Map2k2`, `Mapk3`, `Dusp2`) for every relevant sample. Also export all fixed pathway-member normalized counts.

For descriptive pathway profiles, log1p-transform normalized counts, standardize each gene across the eight original animals using sample SD (ddof = 1), and average the per-gene standardized expression within the fixed pathway. Resolve duplicate symbols by the full-model ranking rule once and keep that choice fixed. Omit zero-variance genes from the standardized mean and record their number. The average is unsigned and gives each eligible member equal weight; it is not weighted by the gene's condition-effect direction. These profiles are not GSEA NES and receive no cell-level significance test. Display animal identities, original condition, pseudobulk cell count and size factor.

## Bounded pilot and resource plan

Available server resources: 64 logical CPUs, approximately 156 GiB available memory; the data volume has only about 23 GiB free. Create a new analysis directory and keep new outputs below 0.5 GiB by storing compressed tables and avoiding 48 large fitted objects. Original data and results remain read-only.

Run one pilot first: lateral-wall fibrocyte/mesenchymal, omitting `Lateral_wall_naive_1` (GSM5124299), the sample previously identified in the descriptive PCA. The pilot is included in the final 48 fits if it completes and passes validation. Bound the fit and enrichment separately by 10-minute timeouts, with an initial 16 GiB virtual-memory ceiling and a single active fit. Record elapsed time, peak memory, warnings and result dimensions. Extrapolate the observed cost before launching the remaining fits; report a material runtime/resource increase to the coordinating agent. No dependency installation is authorized by this plan or needed if the recorded environments verify.

## Outputs and validation

All research outputs are in English. Supply complete compressed DE/GSEA tables for each omitted animal; input and output SHA-256 manifests; software/parameter records; a sample-identity table; primary seven pathway/model summaries with all eight omissions; per-animal normalized-expression tables and figures; and a concise limitations-aware report. Independently check contrast/sample identities, fixed masks/families, ranking construction and BH arithmetic. The main manuscript and original results are outside this task's write scope.

## Execution amendment after the pilot

The pilot completed with 23.9 seconds total DE wall time (12.9 seconds inside fitting/output) and 16.1 seconds total enrichment wall time. Independent input, size-factor, Wald-P, gene-BH, pathway-BH, ranking and membership checks passed. The fixed 954-pathway family remained fully testable. Two implementation-only compatibility issues were corrected during the pilot: gzip input now uses base R instead of requiring optional R.utils, and a data.table validation expression was made explicit. Failed attempts are retained separately and are not scientific outputs.

After the coordinating agent authorized up to two independent concurrent tasks, a resource check found load average 10.62 on 64 CPUs, approximately 153 GiB available memory, and no sbatch/qsub commands. The main queue handles the first 28 manifest entries; a second queue handles the disjoint last 20 entries. The main controller skips already completed tail jobs. At most two fits run simultaneously, each with the unchanged eight inference workers and one-thread BLAS/OpenMP limit (16 inference workers total); fgsea remains one process per task. Every job has its own output directory and seed reset. The expected wall time is approximately 18–22 minutes, within the one-hour controller budget. No scientific parameter changed.
