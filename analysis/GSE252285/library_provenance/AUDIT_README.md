# GSE252285 raw reconstruction audit

## Scope and invariants

- The project entry point is `/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815`.
- Because `/data` had only about 59 GB free at initialization, that path is a symlink to the dedicated new directory on `/path/to/analysis_storage` (about 6.1 TB free at initialization).
- Existing stage1 Seurat/Top20 analyses are out of scope and are not read, changed, or rerun.
- Cell calling and chemistry detection remain automatic. No `--expect-cells` or `--force-cells` parameter is used.
- The output-only flag `--create-bam=false` avoids large BAM files and does not alter barcode/UMI counting or cell calling.

## Provenance evidence

- GEO SOFT/MINiML: GSE252285 and GSM7998628/29/30 sample identities, treatment groups, 10x Chromium 3-prime library description, and GEO assembly annotation.
- NCBI SRA run XML: SRR/SRX/SRS/BioSample/GSM relationships, submitted-original FASTQ filenames/sizes/MD5 values, 2 x 150 bp read statistics, and library descriptors.
- ENA PRJNA1059254 run report: public canonical full-quality paired FASTQ URLs, sizes, and MD5 values.
- Paper: Jiang et al., *Cell Biology and Toxicology* (2024), DOI `10.1007/s10565-024-09912-2`, PMCID `PMC11408584`.

The SRA source objects are marked `Use Cloud Data Delivery` and were not anonymously accessible (HTTP 403 in the access test). The workflow therefore downloads the ENA canonical full-quality FASTQs and validates every compressed file against the ENA MD5. The NCBI submitted-original metadata are retained separately in `metadata/sra_submitted_original_fastq_manifest.tsv` so this archive transformation is explicit.

## Version and reference audit

- The paper says Gel Bead Kit v3, paired-end 150 bp, and Cell Ranger 2.1.0.
- The official GEO Matrix Market header directly records `software_version: cellranger-7.1.0`.
- Cell Ranger 2.1.0 predates support for the v3 chemistry, so the paper's software version is treated as a reporting conflict; reconstruction uses 7.1.0, matching the official matrix header.
- GEO labels the assembly GRCm39, but the official feature list is audited independently against both `refdata-gex-GRCm39-2024-A` and `refdata-gex-mm10-2020-A`. The pipeline proceeds only when exactly one candidate has an exact `(gene_id, gene_name)` set match to all 32,285 official features.

## Library and mapping gates

- Archive identity alone is not used to merge the four HHL SRRs. Before Cell Ranger, the workflow compares the top exact-v3-whitelist R1 barcodes across all six SRRs. B1 SRRs are accepted as technical partitions only if all within-HHL pairwise overlaps are strongly separated from between-SRX overlaps. The evidence is saved in `qc/run_partition_barcode_overlap.tsv` and `qc/run_partition_status.json`.
- Suffix mapping is based only on the required 3 x 3 barcode intersections after removing GEM suffixes. A mapping is resolved only for a one-to-one mutual best match with large intersection, adequate recalls/Jaccard, and at least 20-fold row/column separation from the next-best overlap. Otherwise the workflow stops with `CONDITION_MAPPING_RESOLVED=NO` and does not guess.

## Main outputs

- `metadata/`: GEO, SRA, ENA, paper, provenance, original-submission, and downloaded FASTQ manifests.
- `official/`: official barcodes, features, and matrix.
- `qc/`: FastQC, MultiQC, read-structure, and technical-partition audits.
- `cellranger/<library>/`: complete pipestances and requested web summary, metrics, raw matrix, filtered matrix, commands, and logs.
- `comparison/`: official matrix/reference audit, reconstructed metrics, required 3 x 3 overlap, mapping decision, and cell-count differences.
- `RESULT_SUMMARY.txt`: concise final validated result; it is created only after all required outputs and checksum/mapping gates pass.
