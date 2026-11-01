#!/usr/bin/env python3
"""Download and checksum all ENA FASTQs for GSE252285.

The script deliberately treats the four B1 SRRs as technical partitions of
one biological 10x library because they share SRX23063846, SAMN39202397, and
the GSM7998629 library name in the archive metadata.
"""

import csv
import concurrent.futures
import datetime
import hashlib
import os
import subprocess
import sys
import glob
import xml.etree.ElementTree as ET


PROJECT = os.environ.get(
    "PROJECT", "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
ENA = os.path.join(PROJECT, "metadata", "ENA_PRJNA1059254_read_run.tsv")
ARCHIVE = os.path.join(PROJECT, "fastq", "archive")
CR_INPUT = os.path.join(PROJECT, "fastq", "cellranger_input")
MANIFEST = os.path.join(PROJECT, "metadata", "fastq_manifest.tsv")
PROVENANCE = os.path.join(PROJECT, "metadata", "gsm_srr_library_provenance.tsv")
ORIGINAL_MANIFEST = os.path.join(
    PROJECT, "metadata", "sra_submitted_original_fastq_manifest.tsv"
)
MAX_DOWNLOADS = int(os.environ.get("FASTQ_DOWNLOAD_WORKERS", "4"))

LIBRARIES = {
    "GSM7998628": ("A1_Control", "Control"),
    "GSM7998629": ("B1_HHL", "HHL"),
    "GSM7998630": ("C1_NIHL", "NIHL"),
}


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(16 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def preserve_invalid(path, label):
    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    renamed = path + ".{}.{}".format(label, stamp)
    os.rename(path, renamed)
    print("preserved invalid file as {}".format(renamed), flush=True)


def download_one(url, final_path, expected_bytes, expected_md5):
    os.makedirs(os.path.dirname(final_path), exist_ok=True)
    if os.path.isfile(final_path):
        actual_bytes = os.path.getsize(final_path)
        actual_md5 = md5sum(final_path)
        if actual_bytes == expected_bytes and actual_md5 == expected_md5:
            print("verified existing {}".format(final_path), flush=True)
            return actual_bytes, actual_md5
        preserve_invalid(final_path, "invalid")

    part_path = final_path + ".part"
    for attempt in range(1, 4):
        print(
            "downloading attempt {}/3: {} -> {}".format(attempt, url, final_path),
            flush=True,
        )
        command = [
            "/usr/bin/curl",
            "-fsSL",
            "--retry",
            "50",
            "--retry-delay",
            "15",
            "--connect-timeout",
            "30",
            "--continue-at",
            "-",
            "--output",
            part_path,
            url,
        ]
        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError:
            if attempt == 3:
                raise
            continue

        actual_bytes = os.path.getsize(part_path)
        actual_md5 = md5sum(part_path)
        if actual_bytes == expected_bytes and actual_md5 == expected_md5:
            os.rename(part_path, final_path)
            print(
                "checksum OK: {} bytes md5={} {}".format(
                    actual_bytes, actual_md5, final_path
                ),
                flush=True,
            )
            return actual_bytes, actual_md5

        print(
            "checksum mismatch: expected {} / {}; got {} / {}".format(
                expected_bytes, expected_md5, actual_bytes, actual_md5
            ),
            file=sys.stderr,
            flush=True,
        )
        preserve_invalid(part_path, "checksum_mismatch")

    raise RuntimeError("failed checksum after retries: {}".format(final_path))


def read_ena():
    with open(ENA, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 6:
        raise RuntimeError("expected 6 SRR rows, observed {}".format(len(rows)))
    if set(row["sample_alias"] for row in rows) != set(LIBRARIES):
        raise RuntimeError("unexpected GSM set in ENA metadata")
    return rows


def write_provenance(rows):
    fields = [
        "gsm",
        "srr",
        "srx",
        "biosample",
        "sra_sample",
        "library",
        "condition",
        "technical_partition",
        "library_layout",
        "instrument_model",
        "read_count",
        "base_count",
        "evidence",
    ]
    with open(PROVENANCE, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        by_gsm = {}
        for row in rows:
            by_gsm.setdefault(row["sample_alias"], []).append(row)
        for gsm in sorted(by_gsm):
            library, condition = LIBRARIES[gsm]
            gsm_rows = sorted(by_gsm[gsm], key=lambda item: item["run_accession"])
            experiments = {item["experiment_accession"] for item in gsm_rows}
            biosamples = {item["sample_accession"] for item in gsm_rows}
            names = {item["library_name"] for item in gsm_rows}
            if len(experiments) != 1 or len(biosamples) != 1 or len(names) != 1:
                raise RuntimeError("SRRs for {} do not form one archive library".format(gsm))
            for partition, row in enumerate(gsm_rows, 1):
                writer.writerow(
                    {
                        "gsm": gsm,
                        "srr": row["run_accession"],
                        "srx": row["experiment_accession"],
                        "biosample": row["sample_accession"],
                        "sra_sample": row["secondary_sample_accession"],
                        "library": library,
                        "condition": condition,
                        "technical_partition": partition,
                        "library_layout": row["library_layout"],
                        "instrument_model": row["instrument_model"],
                        "read_count": row["read_count"],
                        "base_count": row["base_count"],
                        "evidence": "same SRX+BioSample+GSM library name => one 10x library",
                    }
                )


def write_sra_original_manifest():
    fields = [
        "gsm",
        "srr",
        "srx",
        "original_filename",
        "original_bytes",
        "original_md5",
        "cloud_uri",
        "archive_access_type",
        "download_decision",
    ]
    original_rows = []
    pattern = os.path.join(PROJECT, "metadata", "NCBI_SRR*_runinfo.xml")
    for path in sorted(glob.glob(pattern)):
        bundle = ET.parse(path).getroot()
        run = bundle.find("RUN")
        srr = run.get("accession")
        srx = run.find("EXPERIMENT_REF").get("accession")
        gsm_element = run.find(
            './Pool/Member/IDENTIFIERS/EXTERNAL_ID[@namespace="GEO"]'
        )
        for sra_file in run.findall("./SRAFiles/SRAFile"):
            if sra_file.get("supertype") != "Original":
                continue
            alternative = sra_file.find("Alternatives")
            original_rows.append(
                {
                    "gsm": gsm_element.text,
                    "srr": srr,
                    "srx": srx,
                    "original_filename": sra_file.get("filename"),
                    "original_bytes": sra_file.get("size"),
                    "original_md5": sra_file.get("md5"),
                    "cloud_uri": alternative.get("url"),
                    "archive_access_type": alternative.get("access_type"),
                    "download_decision": (
                        "NCBI source object requires Cloud Data Delivery; "
                        "download ENA canonical full-quality FASTQ and verify ENA MD5"
                    ),
                }
            )
    if len(original_rows) != 12:
        raise RuntimeError(
            "expected 12 submitted original FASTQ metadata rows; observed {}".format(
                len(original_rows)
            )
        )
    with open(ORIGINAL_MANIFEST, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(original_rows)


def fetch_task(task):
    actual_bytes, actual_md5 = download_one(
        task["url"],
        task["archive_path"],
        task["expected_bytes"],
        task["expected_md5"],
    )
    cr_path = task["cellranger_path"]
    final_path = task["archive_path"]
    library_dir = os.path.dirname(cr_path)
    relative_target = os.path.relpath(final_path, library_dir)
    if os.path.lexists(cr_path):
        if (
            not os.path.islink(cr_path)
            or os.path.realpath(cr_path) != os.path.realpath(final_path)
        ):
            raise RuntimeError("unexpected existing Cell Ranger input: {}".format(cr_path))
    else:
        os.symlink(relative_target, cr_path)

    result = dict(task)
    result.pop("sort_key")
    result.update(
        {
            "actual_bytes": actual_bytes,
            "actual_md5": actual_md5,
            "checksum_status": "PASS",
        }
    )
    return result


def main():
    os.makedirs(ARCHIVE, exist_ok=True)
    os.makedirs(CR_INPUT, exist_ok=True)
    rows = read_ena()
    write_provenance(rows)
    write_sra_original_manifest()

    by_gsm = {}
    for row in rows:
        by_gsm.setdefault(row["sample_alias"], []).append(row)

    download_tasks = []
    for gsm in sorted(by_gsm):
        library, condition = LIBRARIES[gsm]
        library_dir = os.path.join(CR_INPUT, library)
        os.makedirs(library_dir, exist_ok=True)
        for lane, row in enumerate(
            sorted(by_gsm[gsm], key=lambda item: item["run_accession"]), 1
        ):
            urls = row["fastq_ftp"].split(";")
            md5s = row["fastq_md5"].split(";")
            sizes = [int(value) for value in row["fastq_bytes"].split(";")]
            if not (len(urls) == len(md5s) == len(sizes) == 2):
                raise RuntimeError("expected paired FASTQs for {}".format(row["run_accession"]))
            for read_number, (url, expected_md5, expected_bytes) in enumerate(
                zip(urls, md5s, sizes), 1
            ):
                if not url.startswith("http"):
                    # The ENA endpoint is substantially faster over HTTP on this
                    # host. Content integrity is still enforced by exact byte
                    # count and the archive-provided MD5 before a file is accepted.
                    url = "http://" + url
                basename = os.path.basename(url)
                final_path = os.path.join(ARCHIVE, row["run_accession"], basename)
                cr_name = "{}_S1_L{:03d}_R{}_001.fastq.gz".format(
                    library, lane, read_number
                )
                cr_path = os.path.join(library_dir, cr_name)
                download_tasks.append(
                    {
                        "sort_key": (gsm, row["run_accession"], read_number),
                        "gsm": gsm,
                        "library": library,
                        "condition": condition,
                        "srr": row["run_accession"],
                        "srx": row["experiment_accession"],
                        "lane": "L{:03d}".format(lane),
                        "read": "R{}".format(read_number),
                        "url": url,
                        "expected_bytes": expected_bytes,
                        "expected_md5": expected_md5,
                        "archive_path": final_path,
                        "cellranger_path": cr_path,
                    }
                )

    if MAX_DOWNLOADS < 1 or MAX_DOWNLOADS > 4:
        raise RuntimeError("FASTQ_DOWNLOAD_WORKERS must be between 1 and 4")
    download_tasks.sort(key=lambda item: item["sort_key"])
    print(
        "starting {} checksum-gated FASTQ downloads with {} worker(s)".format(
            len(download_tasks), MAX_DOWNLOADS
        ),
        flush=True,
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_DOWNLOADS) as executor:
        manifest_rows = list(executor.map(fetch_task, download_tasks))

    fields = list(manifest_rows[0].keys())
    with open(MANIFEST, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(manifest_rows)

    total_bytes = sum(row["actual_bytes"] for row in manifest_rows)
    print(
        "FASTQ_DOWNLOAD_COMPLETE files={} bytes={} gib={:.3f}".format(
            len(manifest_rows), total_bytes, total_bytes / (1024.0 ** 3)
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
