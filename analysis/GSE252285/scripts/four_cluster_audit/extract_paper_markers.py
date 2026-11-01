import os
from pathlib import Path
import re
import pdfplumber

pdf_path = Path(os.environ.get("PAPER_PDF", "references/s10565-024-09912-2.pdf"))
terms = re.compile(
    r"SGN|spiral ganglion|neuron|Schwann|satellite|supporting|hair cell|"
    r"fibrocyte|chondro|endothelial|pericyte|smooth muscle|Car1|Ces2g|"
    r"marker|cell type|cluster",
    re.IGNORECASE,
)

with pdfplumber.open(pdf_path) as pdf:
    for page_number, page in enumerate(pdf.pages, start=1):
        text = page.extract_text() or ""
        lines = text.splitlines()
        hits = [i for i, line in enumerate(lines) if terms.search(line)]
        if not hits:
            continue
        selected = set()
        for i in hits:
            selected.update(range(max(0, i - 4), min(len(lines), i + 5)))
        print(f"\n===== PDF PAGE {page_number} =====")
        for i in sorted(selected):
            print(lines[i])
