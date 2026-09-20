#!/usr/bin/env python3
"""
setup_datasets.py - Download and assemble annotation databases for ProxiMate.

Usage:
    python setup_datasets.py                     # Download all to ./Datasets/
    python setup_datasets.py --output-dir /path  # Custom output directory
    python setup_datasets.py --force             # Re-download even if files exist
    python setup_datasets.py --skip biogrid hpa  # Skip specific datasets
"""

import argparse
from datetime import datetime, timezone
import gzip
import io
import os
import subprocess
import sys
import zipfile
import requests

# ---------------------------------------------------------------------------
# Organism Configuration
# ---------------------------------------------------------------------------

# Supported organisms for ProxiMate annotation.
# To add a new organism:
#   1. Add an entry here with its NCBI Taxonomy ID and feature flags
#   2. Sync this config in both setup_datasets.py and Scripts/annotator.py
#   3. Add the organism to the GUI dropdown in GUI/app.py (scoring panel)
#   4. If the organism has a species-specific database (like HPA for human),
#      add a download function in setup_datasets.py and conditional logic in annotator.py
ORGANISMS = {
    "human": {"organism_id": 9606, "has_hpa": True, "has_corum": True, "has_hcm": True},
    "mouse": {"organism_id": 10090, "has_hpa": False, "has_corum": False, "has_hcm": False},
    "yeast": {"organism_id": 559292, "has_hpa": False, "has_corum": False, "has_hcm": False},
}

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# UniProt REST API URL template — organism_id is filled per organism
UNIPROT_STREAM_URL_TEMPLATE = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=(organism_id:{organism_id})+AND+(reviewed:true)"
    "&format=tsv"
    "&fields=accession,gene_names,cc_subcellular_location,go_c,go_p,go_f,"
    "ft_motif,ft_region,ft_repeat,ft_compbias,ft_domain"
    "&compressed=true"
)
UNIPROT_FILENAME = "uniprot_anns.tsv"
UNIPROT_REQUIRED_COLUMNS = {
    "Entry", "Gene Names", "Subcellular location [CC]",
    "Gene Ontology (cellular component)",
    "Gene Ontology (biological process)",
    "Gene Ontology (molecular function)",
    "Motif", "Region",
    "Repeat", "Compositional bias", "Domain [FT]",
}

BIOGRID_ALL_URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/"
    "Latest-Release/BIOGRID-ALL-LATEST.tab3.zip"
)
BIOGRID_MV_URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/"
    "Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
)
BIOGRID_ALL_FILENAME = "BIOGRID-ALL.tab3.txt"
BIOGRID_MV_FILENAME = "BIOGRID-MV-Physical.tab3.txt"
BIOGRID_SUMMARY_FILENAME = "biogrid_summary.csv"
# The same summary with Human Cell Map (Go et al. 2021) evidence removed, built only for
# organisms flagged has_hcm above.  HCM is itself a BioID screen, so proximity-labeling
# runs scored against it look more "known" than they are.
BIOGRID_NO_HCM_SUMMARY_FILENAME = "biogrid_summary_no_hcm.csv"
HCM_PUBLICATION = "PUBMED:34079125"
BIOGRID_REQUIRED_COLUMNS = {
    "Organism ID Interactor A", "Organism ID Interactor B",
    "Experimental System Type", "Publication Source", "SWISS-PROT Accessions Interactor A",
    "SWISS-PROT Accessions Interactor B",
}

HPA_URL = "https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip"
HPA_FILENAME = "subcellular_location.tsv"
HPA_REQUIRED_COLUMNS = {"Gene name", "Main location"}

CORUM_URL = "https://mips.helmholtz-muenchen.de/corum/download/humanComplexes.txt.zip"
CORUM_FILENAME = "corum_humanComplexes.txt"
CORUM_REQUIRED_COLUMNS = {"complex_name", "subunits_uniprot_id"}

BUILD_INFO_FILENAME = "build_info.txt"

REQUEST_TIMEOUT = 300
CHUNK_SIZE = 8192
USER_AGENT = "ProxiMate-setup/1.0 (https://github.com/plutzer/ProxiMate)"

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def log(message: str) -> None:
    print(f"[setup_datasets] {message}", file=sys.stderr, flush=True)


def file_exists_and_nonempty(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def verify_tsv_columns(filepath: str, required_columns: set, encoding: str = "utf-8") -> bool:
    """Whether the header line of a tab-separated file names every required column."""
    try:
        with open(filepath, "r", encoding=encoding) as f:
            header_line = f.readline().strip()
        if not header_line:
            log(f"  Verification: {filepath} is empty")
            return False
        actual_columns = set(header_line.split("\t"))
        missing = required_columns - actual_columns
        if missing:
            log(f"  Verification FAILED for {os.path.basename(filepath)}")
            log(f"  Missing columns: {missing}")
            return False
        return True
    except UnicodeDecodeError:
        if encoding != "latin-1":
            return verify_tsv_columns(filepath, required_columns, "latin-1")
        log(f"  Verification: Could not decode {filepath}")
        return False
    except Exception as e:
        log(f"  Verification error: {e}")
        return False


def download(url: str, description: str) -> bytes:
    response = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    log(f"  {description}: Downloaded {len(response.content) / 1024 / 1024:.1f} MB")
    return response.content


def extract_from_zip(content: bytes, suffix: str, target_path: str, label: str) -> bool:
    """Extract a file matching *suffix* from a ZIP archive to target_path."""
    # Detect HTML response (e.g. BioGRID license wall)
    if content[:5] in (b"<html", b"<!DOC", b"<HTML", b"<!doc"):
        log(f"  {label}: Received HTML instead of ZIP — site may require manual download")
        return False
    try:
        with zipfile.ZipFile(io.BytesIO(content)) as zf:
            matches = [n for n in zf.namelist() if n.endswith(suffix)]
            if not matches:
                log(f"  {label}: No *{suffix} file in archive. Contents: {zf.namelist()}")
                return False
            source_name = matches[0]
            log(f"  {label}: Extracting {source_name}")
            os.makedirs(os.path.dirname(target_path) or ".", exist_ok=True)
            with zf.open(source_name) as src, open(target_path, "wb") as dst:
                while True:
                    chunk = src.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    dst.write(chunk)
        return True
    except zipfile.BadZipFile:
        log(f"  {label}: Downloaded file is not a valid ZIP")
        return False

# ---------------------------------------------------------------------------
# Per-dataset downloads
# ---------------------------------------------------------------------------

def download_uniprot(organism_dir: str, force: bool, organism_name: str, organism_id: int) -> bool:
    """Download UniProt annotations for a specific organism."""
    os.makedirs(organism_dir, exist_ok=True)
    target = os.path.join(organism_dir, UNIPROT_FILENAME)
    label = f"UniProt ({organism_name})"

    if file_exists_and_nonempty(target) and not force:
        log(f"{label}: {UNIPROT_FILENAME} already exists, skipping (use --force to re-download)")
        return verify_tsv_columns(target, UNIPROT_REQUIRED_COLUMNS)

    url = UNIPROT_STREAM_URL_TEMPLATE.format(organism_id=organism_id)
    log(f"{label}: Downloading reviewed proteome annotations (organism_id={organism_id})...")
    try:
        content = download(url, label)
        try:
            decompressed = gzip.decompress(content)
        except gzip.BadGzipFile:
            log(f"  {label}: Response was not gzip — using raw content")
            decompressed = content

        with open(target, "wb") as f:
            f.write(decompressed)
        log(f"{label}: Written {os.path.getsize(target):,} bytes to {UNIPROT_FILENAME}")

        if not verify_tsv_columns(target, UNIPROT_REQUIRED_COLUMNS):
            log(f"{label}: Column verification FAILED")
            return False
        log(f"{label}: Column verification passed")
        return True
    except requests.RequestException as e:
        log(f"{label}: Download failed — {e}")
        return False


def download_biogrid(output_dir: str, force: bool) -> bool:
    """Download BioGRID bulk files (shared across all organisms)."""
    downloads = [
        (BIOGRID_ALL_URL, BIOGRID_ALL_FILENAME, "BioGRID-ALL"),
        (BIOGRID_MV_URL, BIOGRID_MV_FILENAME, "BioGRID-MV"),
    ]
    all_ok = True
    for url, target_name, label in downloads:
        target = os.path.join(output_dir, target_name)
        if file_exists_and_nonempty(target) and not force:
            log(f"{label}: {target_name} already exists, skipping")
            continue
        log(f"{label}: Downloading...")
        try:
            content = download(url, label)
            if not extract_from_zip(content, ".tab3.txt", target, label):
                log(f"{label}: Try downloading manually from https://downloads.thebiogrid.org/BioGRID")
                all_ok = False
                continue
            log(f"{label}: Written {os.path.getsize(target):,} bytes")
            if not verify_tsv_columns(target, BIOGRID_REQUIRED_COLUMNS):
                log(f"{label}: Column verification FAILED")
                all_ok = False
        except requests.RequestException as e:
            log(f"{label}: Download failed — {e}")
            all_ok = False
    return all_ok


def download_hpa(organism_dir: str, force: bool) -> bool:
    """Download Human Protein Atlas subcellular location data (human only)."""
    os.makedirs(organism_dir, exist_ok=True)
    target = os.path.join(organism_dir, HPA_FILENAME)
    if file_exists_and_nonempty(target) and not force:
        log(f"HPA: {HPA_FILENAME} already exists, skipping")
        return verify_tsv_columns(target, HPA_REQUIRED_COLUMNS)

    log("HPA: Downloading subcellular location data...")
    try:
        content = download(HPA_URL, "HPA")
        if not extract_from_zip(content, ".tsv", target, "HPA"):
            return False
        log(f"HPA: Written {os.path.getsize(target):,} bytes")
        if not verify_tsv_columns(target, HPA_REQUIRED_COLUMNS):
            log("HPA: Column verification FAILED — the 'Main location' column may have been renamed in newer HPA versions")
            return False
        log("HPA: Column verification passed")
        return True
    except requests.RequestException as e:
        log(f"HPA: Download failed — {e}")
        return False


def download_corum(output_dir: str, force: bool) -> bool:
    """Download CORUM protein complexes (human only, stored at output_dir root)."""
    target = os.path.join(output_dir, CORUM_FILENAME)
    if file_exists_and_nonempty(target) and not force:
        log(f"CORUM: {CORUM_FILENAME} already exists, skipping")
        return verify_tsv_columns(target, CORUM_REQUIRED_COLUMNS, encoding="latin-1")

    log("CORUM: Downloading human protein complexes...")
    try:
        content = download(CORUM_URL, "CORUM")
    except requests.RequestException as e:
        log(f"CORUM: Download failed — {e}")
        log("CORUM: Please download manually from https://mips.helmholtz-muenchen.de/corum/")
        log(f"CORUM: Save the 'Human Complexes' file as {target}")
        return False

    if not extract_from_zip(content, ".txt", target, "CORUM"):
        return False
    if not verify_tsv_columns(target, CORUM_REQUIRED_COLUMNS, encoding="latin-1"):
        log("CORUM: Column verification FAILED")
        return False
    log(f"CORUM: Written {os.path.getsize(target):,} bytes, verification passed")
    return True

# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def run_preprocess_biogrid(output_dir: str, organism_name: str, organism_id: int,
                           exclude_hcm: bool = False) -> bool:
    """Run BioGRID preprocessing for a specific organism.

    With ``exclude_hcm`` the summary is written under the no-HCM filename with every
    Human Cell Map evidence row removed first.
    """
    biogrid_all = os.path.join(output_dir, BIOGRID_ALL_FILENAME)
    biogrid_mv = os.path.join(output_dir, BIOGRID_MV_FILENAME)

    if not file_exists_and_nonempty(biogrid_all) or not file_exists_and_nonempty(biogrid_mv):
        log("BioGRID preprocessing: Input files missing, skipping")
        return False

    script_dir = os.path.dirname(os.path.abspath(__file__))
    preprocess_script = os.path.join(script_dir, "preprocess_biogrid.py")
    if not os.path.isfile(preprocess_script):
        log(f"BioGRID preprocessing: Cannot find {preprocess_script}")
        return False

    organism_dir = os.path.join(output_dir, organism_name)
    os.makedirs(organism_dir, exist_ok=True)
    label = f"BioGRID preprocessing ({organism_name})"
    summary_filename = BIOGRID_SUMMARY_FILENAME
    cmd = [
        sys.executable, preprocess_script,
        "--biogrid_all", biogrid_all,
        "--biogrid_mv", biogrid_mv,
        "--output_dir", organism_dir,
        "--organism_id", str(organism_id),
    ]
    if exclude_hcm:
        label = f"BioGRID preprocessing ({organism_name}, no HCM)"
        summary_filename = BIOGRID_NO_HCM_SUMMARY_FILENAME
        cmd += ["--exclude_publication", HCM_PUBLICATION,
                "--output_filename", summary_filename]

    log(f"{label}: Generating {summary_filename} (organism_id={organism_id})...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            log(f"{label}: Failed (exit code {result.returncode})")
            if result.stderr:
                log(f"  stderr: {result.stderr[:500]}")
            return False

        summary = os.path.join(organism_dir, summary_filename)
        if not file_exists_and_nonempty(summary):
            log(f"{label}: {summary_filename} was not created")
            return False
        log(f"{label}: Generated {summary_filename} ({os.path.getsize(summary):,} bytes)")
        return True
    except subprocess.TimeoutExpired:
        log(f"{label}: Timed out after 300 seconds")
        return False

# ---------------------------------------------------------------------------
# Build info
# ---------------------------------------------------------------------------

def write_build_info(output_dir: str, results: dict) -> None:
    """Write build_info.txt: the build date and each download's outcome."""
    target = os.path.join(output_dir, BUILD_INFO_FILENAME)
    lines = [f"Build date: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}"]
    for name, success in results.items():
        lines.append(f"  {name}: {'OK' if success else 'FAILED'}")
    with open(target, "w") as f:
        f.write("\n".join(lines) + "\n")
    log(f"Build info written to {target}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download and assemble annotation databases for ProxiMate.",
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory for dataset files (default: ./Datasets/)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-download files even if they already exist",
    )
    parser.add_argument(
        "--skip", nargs="*", choices=["uniprot", "biogrid", "hpa", "corum"],
        default=[],
        help="Skip specific dataset downloads",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = (
        os.path.abspath(args.output_dir)
        if args.output_dir
        else os.path.join(script_dir, "..", "Datasets")
    )
    os.makedirs(output_dir, exist_ok=True)
    log(f"Output directory: {output_dir}")

    results = {}

    # UniProt: download per organism into {output_dir}/{organism}/
    if "uniprot" not in args.skip:
        for org_name, org_config in ORGANISMS.items():
            organism_dir = os.path.join(output_dir, org_name)
            results[f"UniProt ({org_name})"] = download_uniprot(
                organism_dir, args.force, org_name, org_config["organism_id"]
            )

    # BioGRID: download raw bulk files once (shared), then preprocess per organism
    if "biogrid" not in args.skip:
        results["BioGRID"] = download_biogrid(output_dir, args.force)
        if results["BioGRID"]:
            for org_name, org_config in ORGANISMS.items():
                results[f"BioGRID preprocessing ({org_name})"] = run_preprocess_biogrid(
                    output_dir, org_name, org_config["organism_id"]
                )
                if org_config["has_hcm"]:
                    results[f"BioGRID preprocessing ({org_name}, no HCM)"] = run_preprocess_biogrid(
                        output_dir, org_name, org_config["organism_id"], exclude_hcm=True
                    )

    # HPA: human only
    if "hpa" not in args.skip:
        human_dir = os.path.join(output_dir, "human")
        results["Human Protein Atlas"] = download_hpa(human_dir, args.force)

    # CORUM: human only, stored at output_dir root (git-tracked)
    if "corum" not in args.skip:
        results["CORUM"] = download_corum(output_dir, args.force)

    # Write build info log
    write_build_info(output_dir, results)

    log("=" * 50)
    log("SUMMARY")
    log("=" * 50)
    for name, success in results.items():
        status = "OK" if success else "FAILED"
        log(f"  {name}: {status}")

    failed = [n for n, s in results.items() if not s]
    if failed:
        log(f"\n{len(failed)} dataset(s) failed. See messages above.")
        return 1
    log(f"\nAll datasets downloaded successfully to {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
