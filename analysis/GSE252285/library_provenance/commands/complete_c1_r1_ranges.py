#!/usr/bin/env python3
"""Complete the remaining C1 R1 suffix with four checksum-gated HTTP ranges."""

import concurrent.futures
import datetime
import hashlib
import os
import shutil
import subprocess


PROJECT = os.environ.get(
    "PROJECT", "/path/to/cochlear_analysis/GSE252285_raw_reconstruction_20260815"
)
FINAL = os.path.join(
    PROJECT, "fastq", "archive", "SRR27387600", "SRR27387600_1.fastq.gz"
)
PART = FINAL + ".part"
URL = (
    "http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR273/000/"
    "SRR27387600/SRR27387600_1.fastq.gz"
)
EXPECTED_BYTES = 27305385846
EXPECTED_MD5 = "37836dc437f44c32453e82595709b3e7"
RANGE_ROOT = os.path.join(PROJECT, "tmp", "c1_r1_http_ranges")


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(16 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def preserve_invalid(path):
    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    renamed = path + ".invalid." + stamp
    os.rename(path, renamed)
    print("preserved invalid range as {}".format(renamed), flush=True)


def download_range(item):
    index, start, end = item
    expected = end - start + 1
    path = os.path.join(
        RANGE_ROOT, "range_{:02d}_{}_{}.bin".format(index, start, end)
    )
    if os.path.isfile(path) and os.path.getsize(path) == expected:
        print("verified existing range {} bytes {}".format(expected, path), flush=True)
        return path
    if os.path.exists(path):
        preserve_invalid(path)

    for attempt in range(1, 4):
        print(
            "range attempt {}/3: bytes {}-{} -> {}".format(
                attempt, start, end, path
            ),
            flush=True,
        )
        try:
            subprocess.check_call(
                [
                    "/usr/bin/curl",
                    "-fsSL",
                    "--retry",
                    "50",
                    "--retry-delay",
                    "15",
                    "--connect-timeout",
                    "30",
                    "--range",
                    "{}-{}".format(start, end),
                    "--output",
                    path,
                    URL,
                ]
            )
        except subprocess.CalledProcessError:
            if os.path.exists(path):
                preserve_invalid(path)
            if attempt == 3:
                raise
            continue
        actual = os.path.getsize(path)
        if actual == expected:
            print("range size PASS: {} bytes {}".format(actual, path), flush=True)
            return path
        print(
            "range size mismatch: expected {} observed {}".format(expected, actual),
            flush=True,
        )
        preserve_invalid(path)
    raise RuntimeError("range failed after retries: {}-{}".format(start, end))


def main():
    os.makedirs(RANGE_ROOT, exist_ok=True)
    if os.path.isfile(FINAL):
        actual_size = os.path.getsize(FINAL)
        actual_md5 = md5sum(FINAL)
        if actual_size == EXPECTED_BYTES and actual_md5 == EXPECTED_MD5:
            print("C1_R1_ALREADY_COMPLETE")
            return
        raise RuntimeError("unexpected existing final C1 R1")
    if not os.path.isfile(PART):
        raise RuntimeError("C1 R1 prefix part is missing")
    prefix_bytes = os.path.getsize(PART)
    if not 0 < prefix_bytes < EXPECTED_BYTES:
        raise RuntimeError("unexpected C1 R1 prefix size: {}".format(prefix_bytes))

    remaining = EXPECTED_BYTES - prefix_bytes
    ranges = []
    cursor = prefix_bytes
    for index in range(4):
        slots = 4 - index
        length = (EXPECTED_BYTES - cursor + slots - 1) // slots
        end = min(EXPECTED_BYTES - 1, cursor + length - 1)
        ranges.append((index + 1, cursor, end))
        cursor = end + 1
    if cursor != EXPECTED_BYTES:
        raise RuntimeError("range planning did not cover the expected file")

    print(
        "C1 R1 prefix={} remaining={} ranges={}".format(
            prefix_bytes, remaining, ranges
        ),
        flush=True,
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        paths = list(executor.map(download_range, ranges))

    if os.path.getsize(PART) != prefix_bytes:
        raise RuntimeError("C1 R1 prefix changed during range downloads")
    with open(PART, "ab") as output:
        for path, (_, start, end) in zip(paths, ranges):
            expected = end - start + 1
            if os.path.getsize(path) != expected:
                raise RuntimeError("range changed before assembly: {}".format(path))
            with open(path, "rb") as source:
                shutil.copyfileobj(source, output, length=16 * 1024 * 1024)
            print("assembled range {}-{}".format(start, end), flush=True)
        output.flush()
        os.fsync(output.fileno())

    actual_size = os.path.getsize(PART)
    actual_md5 = md5sum(PART)
    if actual_size != EXPECTED_BYTES or actual_md5 != EXPECTED_MD5:
        raise RuntimeError(
            "assembled C1 R1 failed validation: size={} md5={}".format(
                actual_size, actual_md5
            )
        )
    os.rename(PART, FINAL)
    print(
        "C1_R1_RANGE_COMPLETION_PASS bytes={} md5={}".format(
            actual_size, actual_md5
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
