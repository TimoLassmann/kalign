"""Dataset management for kalign benchmarks.

Supports BAliBASE, BRAliBASE, and BaliFam100 benchmark datasets.
"""

import shutil
import subprocess
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

BENCHMARKS_DIR = Path(__file__).parent
DATA_DIR = BENCHMARKS_DIR / "data"
DOWNLOADS_DIR = DATA_DIR / "downloads"


@dataclass
class BenchmarkCase:
    """A single benchmark test case (one MSA family)."""

    family: str
    dataset: str
    unaligned: Path
    reference: Path
    seq_type: str  # "protein" | "rna" | "dna"


# ---------------------------------------------------------------------------
# BAliBASE
# ---------------------------------------------------------------------------

BALIBASE_URL = "http://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz"
BALIBASE_DIR = DOWNLOADS_DIR / "bb3_release"


def balibase_download() -> None:
    """Download and extract BAliBASE R1-5."""
    DOWNLOADS_DIR.mkdir(parents=True, exist_ok=True)
    tarball = DOWNLOADS_DIR / "BAliBASE_R1-5.tar.gz"

    if BALIBASE_DIR.exists():
        return

    print(f"Downloading BAliBASE from {BALIBASE_URL} ...")
    urllib.request.urlretrieve(BALIBASE_URL, tarball)
    print("Extracting ...")
    with tarfile.open(tarball) as tf:
        tf.extractall(DOWNLOADS_DIR)
    tarball.unlink()


def balibase_is_available() -> bool:
    return BALIBASE_DIR.exists()


def balibase_cases() -> List[BenchmarkCase]:
    """Discover BAliBASE test cases (*.tfa paired with *.msf)."""
    cases = []
    for tfa in sorted(BALIBASE_DIR.rglob("*.tfa")):
        # Skip BBS (truncated) cases — only use full-length BB alignments
        if tfa.stem.startswith("BBS"):
            continue
        msf = tfa.with_suffix(".msf")
        if not msf.exists():
            continue
        # Determine RV category from path
        parts = tfa.parts
        rv = ""
        for p in parts:
            if p.startswith("RV"):
                rv = p
                break
        family = tfa.stem
        cases.append(
            BenchmarkCase(
                family=family,
                dataset=f"balibase_{rv}" if rv else "balibase",
                unaligned=tfa,
                reference=msf,
                seq_type="protein",
            )
        )
    return cases


# ---------------------------------------------------------------------------
# BRAliBASE
# ---------------------------------------------------------------------------

BRALIBASE_URLS = [
    ("data-set1",
     "https://web.archive.org/web/20160408045442id_/http://projects.binf.ku.dk/pgardner/bralibase/data-set1.tar.gz"),
    ("data-set2",
     "https://web.archive.org/web/20160408045453id_/http://projects.binf.ku.dk/pgardner/bralibase/data-set2.tar.gz"),
]


def bralibase_download() -> None:
    """Download and extract BRAliBASE set1 + set2."""
    DOWNLOADS_DIR.mkdir(parents=True, exist_ok=True)

    for name, url in BRALIBASE_URLS:
        target = DOWNLOADS_DIR / name
        if target.exists():
            continue
        tarball = DOWNLOADS_DIR / f"{name}.tar.gz"
        print(f"Downloading BRAliBASE {name} from {url} ...")
        urllib.request.urlretrieve(url, tarball)
        print("Extracting ...")
        with tarfile.open(tarball) as tf:
            tf.extractall(DOWNLOADS_DIR)
        tarball.unlink()


def bralibase_is_available() -> bool:
    return any((DOWNLOADS_DIR / name).exists() for name, _ in BRALIBASE_URLS)


def bralibase_cases() -> List[BenchmarkCase]:
    """Discover BRAliBASE test cases.

    data-set1 structure: {family}/structural/*.fa paired with {family}/unaligned/*.fa
      families: g2intron, rRNA, SRP, tRNA, U5
    data-set2 structure: structural/*.fasta paired with unaligned/*.fasta
    """
    cases = []
    for name, _ in BRALIBASE_URLS:
        base = DOWNLOADS_DIR / name
        if not base.exists():
            continue
        # Glob for both .fa and .fasta reference files
        refs = sorted(list(base.rglob("structural/*.fa")) +
                      list(base.rglob("structural/*.fasta")))
        for ref in refs:
            # Corresponding unaligned file — try same extension first, then alternative
            unaligned_dir = ref.parent.parent / "unaligned"
            unaligned = unaligned_dir / ref.name
            if not unaligned.exists():
                alt_ext = ".fasta" if ref.suffix == ".fa" else ".fa"
                unaligned = unaligned_dir / (ref.stem + alt_ext)
                if not unaligned.exists():
                    continue
            family = ref.stem
            # For data-set1, extract RNA family category from path
            # e.g. data-set1/tRNA/structural/aln1.fa -> category "tRNA"
            structural_parent = ref.parent.parent.name
            if structural_parent != name:
                # This is a subcategory (g2intron, rRNA, SRP, tRNA, U5)
                dataset = f"bralibase_{structural_parent}"
            else:
                dataset = f"bralibase_{name}"
            cases.append(
                BenchmarkCase(
                    family=family,
                    dataset=dataset,
                    unaligned=unaligned,
                    reference=ref,
                    seq_type="rna",
                )
            )
    return cases


# ---------------------------------------------------------------------------
# BaliFam100 (CC0 - committed to repo)
# ---------------------------------------------------------------------------

BALIFAM_DIR = DATA_DIR / "balifam100"
BALIFAM_REPO = "https://github.com/rcedgar/balifam.git"


def balifam_download() -> None:
    """Clone BaliFam100 repo (CC0 license) into data/balifam100/."""
    if BALIFAM_DIR.exists() and any(BALIFAM_DIR.iterdir()):
        return

    BALIFAM_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Cloning BaliFam100 from {BALIFAM_REPO} ...")
    subprocess.run(
        ["git", "clone", "--depth", "1", BALIFAM_REPO, str(BALIFAM_DIR)],
        check=True,
    )
    # Remove .git to keep it clean for committing
    git_dir = BALIFAM_DIR / ".git"
    if git_dir.exists():
        shutil.rmtree(git_dir)


def balifam_is_available() -> bool:
    return BALIFAM_DIR.exists() and (
        any(BALIFAM_DIR.rglob("*.fa")) or any(BALIFAM_DIR.rglob("*.fasta"))
    )


def balifam_cases() -> List[BenchmarkCase]:
    """Discover BaliFam100 test cases."""
    cases = []
    # BaliFam100 has in/ (unaligned) and ref/ (reference) subdirectories
    in_dir = BALIFAM_DIR / "in"
    ref_dir = BALIFAM_DIR / "ref"

    if not in_dir.exists() or not ref_dir.exists():
        # Try alternative structure
        for ref in sorted(BALIFAM_DIR.rglob("ref/*.fa")):
            in_file = ref.parent.parent / "in" / ref.name
            if not in_file.exists():
                in_file = ref.parent.parent / "in" / ref.with_suffix(".fasta").name
                if not in_file.exists():
                    continue
            cases.append(
                BenchmarkCase(
                    family=ref.stem,
                    dataset="balifam100",
                    unaligned=in_file,
                    reference=ref,
                    seq_type="protein",
                )
            )
        return cases

    for ref in sorted(ref_dir.glob("*")):
        if not ref.is_file():
            continue
        in_file = in_dir / ref.name
        if not in_file.exists():
            continue
        cases.append(
            BenchmarkCase(
                family=ref.stem,
                dataset="balifam100",
                unaligned=in_file,
                reference=ref,
                seq_type="protein",
            )
        )
    return cases


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

DATASETS = {
    "balibase": {
        "download": balibase_download,
        "is_available": balibase_is_available,
        "cases": balibase_cases,
    },
    "bralibase": {
        "download": bralibase_download,
        "is_available": bralibase_is_available,
        "cases": bralibase_cases,
    },
    "balifam100": {
        "download": balifam_download,
        "is_available": balifam_is_available,
        "cases": balifam_cases,
    },
}


def download_dataset(name: str) -> None:
    """Download a specific dataset by name, or 'all' for everything."""
    if name == "all":
        for ds in DATASETS.values():
            ds["download"]()
    elif name in DATASETS:
        DATASETS[name]["download"]()
    else:
        raise ValueError(f"Unknown dataset: {name}. Choose from: {list(DATASETS.keys())} or 'all'")


def get_cases(name: str, max_cases: Optional[int] = None) -> List[BenchmarkCase]:
    """Get benchmark cases for a dataset. Downloads if not available."""
    if name == "all":
        all_cases = []
        for ds in DATASETS.values():
            if ds["is_available"]():
                all_cases.extend(ds["cases"]())
        if max_cases:
            all_cases = all_cases[:max_cases]
        return all_cases

    if name not in DATASETS:
        raise ValueError(f"Unknown dataset: {name}. Choose from: {list(DATASETS.keys())}")

    ds = DATASETS[name]
    if not ds["is_available"]():
        ds["download"]()
    cases = ds["cases"]()
    if max_cases:
        cases = cases[:max_cases]
    return cases
