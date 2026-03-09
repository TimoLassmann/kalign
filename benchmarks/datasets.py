"""Dataset management for kalign benchmarks.

Supports BAliBASE, BRAliBASE, and BaliFam100 benchmark datasets.
"""

import shutil
import subprocess
import tarfile
import urllib.error
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
BALIBASE_URL_FALLBACK = "https://web.archive.org/web/20231117003408/https://www.lbgi.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz"
BALIBASE_DIR = DOWNLOADS_DIR / "bb3_release"


def _validate_tarball(path: Path) -> bool:
    """Check that a file is a valid gzip tarball (not an HTML error page)."""
    try:
        with tarfile.open(path) as tf:
            tf.getnames()
        return True
    except (tarfile.ReadError, tarfile.CompressionError, EOFError):
        return False


def balibase_download() -> None:
    """Download and extract BAliBASE R1-5."""
    DOWNLOADS_DIR.mkdir(parents=True, exist_ok=True)
    tarball = DOWNLOADS_DIR / "BAliBASE_R1-5.tar.gz"

    if BALIBASE_DIR.exists():
        return

    if tarball.exists() and not _validate_tarball(tarball):
        print(f"  Removing corrupt tarball: {tarball}")
        tarball.unlink()

    if not tarball.exists():
        for url in [BALIBASE_URL, BALIBASE_URL_FALLBACK]:
            try:
                print(f"Downloading BAliBASE from {url} ...")
                urllib.request.urlretrieve(url, tarball)
                if not _validate_tarball(tarball):
                    print(f"  Downloaded file is not a valid tarball (got HTML error page?)")
                    tarball.unlink()
                    continue
                break
            except (urllib.error.HTTPError, urllib.error.URLError) as e:
                print(f"  Failed: {e}")
                if tarball.exists():
                    tarball.unlink()
                continue
        else:
            raise RuntimeError(
                "Could not download BAliBASE from any URL.\n"
                "Please download manually and place at:\n"
                f"  {tarball}\n"
                "Or extract bb3_release/ into:\n"
                f"  {DOWNLOADS_DIR}/"
            )

    print("Extracting ...")
    with tarfile.open(tarball) as tf:
        tf.extractall(DOWNLOADS_DIR)


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
        if tarball.exists() and not _validate_tarball(tarball):
            print(f"  Removing corrupt tarball: {tarball}")
            tarball.unlink()
        if not tarball.exists():
            print(f"Downloading BRAliBASE {name} from {url} ...")
            try:
                urllib.request.urlretrieve(url, tarball)
            except (urllib.error.HTTPError, urllib.error.URLError) as e:
                raise RuntimeError(f"Could not download BRAliBASE {name}: {e}")
            if not _validate_tarball(tarball):
                tarball.unlink()
                raise RuntimeError(
                    f"Downloaded BRAliBASE {name} is not a valid tarball.\n"
                    f"Please download manually from:\n  {url}\n"
                    f"and place at:\n  {tarball}"
                )
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
# MDSA (Multiple DNA Sequence Alignments — Carroll et al. 2007)
# ---------------------------------------------------------------------------

MDSA_BASE_URL = "https://dna.cs.byu.edu/mdsas/tarred"
MDSA_DIR = DOWNLOADS_DIR / "mdsa"
# Skip PREFAB (all pairwise, 2 seqs — not useful for MSA benchmarking)
MDSA_DATABASES = ["balibase", "oxbench", "smart"]
MDSA_VERSION = "all"  # "100s" is too small (400 cases); "all" has 1,869 usable cases
MDSA_MAX_SEQS = 500   # Skip very large families (SMART has up to 1,769 seqs)


def mdsa_download() -> None:
    """Download and extract MDSA DNA alignment benchmark datasets."""
    MDSA_DIR.mkdir(parents=True, exist_ok=True)

    for db in MDSA_DATABASES:
        db_dir = MDSA_DIR / f"{db}_mdsa_{MDSA_VERSION}"
        if db_dir.exists() and any(db_dir.glob("*.afa")):
            continue

        tarball = MDSA_DIR / f"{db}_mdsa_{MDSA_VERSION}.tar.gz"
        if tarball.exists() and not _validate_tarball(tarball):
            print(f"  Removing corrupt tarball: {tarball}")
            tarball.unlink()
        if not tarball.exists():
            url = f"{MDSA_BASE_URL}/{db}_mdsa_{MDSA_VERSION}.tar.gz"
            print(f"Downloading MDSA {db} from {url} ...")
            urllib.request.urlretrieve(url, tarball)
            if not _validate_tarball(tarball):
                tarball.unlink()
                raise RuntimeError(
                    f"Downloaded MDSA {db} is not a valid tarball.\n"
                    f"Please download manually from:\n  {url}\n"
                    f"and place at:\n  {tarball}"
                )

        print(f"Extracting {db} ...")
        with tarfile.open(tarball) as tf:
            tf.extractall(MDSA_DIR)

    # Generate unaligned files by stripping gaps from reference .afa files
    unaligned_dir = MDSA_DIR / "unaligned"
    unaligned_dir.mkdir(exist_ok=True)

    for db in MDSA_DATABASES:
        db_dir = MDSA_DIR / f"{db}_mdsa_{MDSA_VERSION}"
        if not db_dir.exists():
            continue
        db_unaligned = unaligned_dir / db
        db_unaligned.mkdir(exist_ok=True)
        for afa in sorted(db_dir.glob("*.afa")):
            out = db_unaligned / afa.name
            if out.exists():
                continue
            _strip_gaps_fasta(afa, out)

    print(f"MDSA ready: {sum(1 for _ in unaligned_dir.rglob('*.afa'))} unaligned files")


def _strip_gaps_fasta(ref_path: Path, out_path: Path) -> None:
    """Read aligned FASTA, strip gap characters, write unaligned FASTA."""
    sequences = []
    current_header = None
    current_seq = []

    with open(ref_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    seq = "".join(current_seq).replace("-", "").replace(".", "")
                    sequences.append((current_header, seq))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            seq = "".join(current_seq).replace("-", "").replace(".", "")
            sequences.append((current_header, seq))

    with open(out_path, "w") as f:
        for header, seq in sequences:
            f.write(f"{header}\n{seq}\n")


def mdsa_is_available() -> bool:
    unaligned_dir = MDSA_DIR / "unaligned"
    return unaligned_dir.exists() and any(unaligned_dir.rglob("*.afa"))


def mdsa_cases() -> List[BenchmarkCase]:
    """Discover MDSA test cases.

    Reference = original .afa (aligned FASTA)
    Unaligned = gap-stripped version in unaligned/{db}/*.afa
    Skips PREFAB (pairwise) and very large families (>MDSA_MAX_SEQS).
    """
    cases = []
    unaligned_dir = MDSA_DIR / "unaligned"

    for db in MDSA_DATABASES:
        ref_dir = MDSA_DIR / f"{db}_mdsa_{MDSA_VERSION}"
        ua_dir = unaligned_dir / db

        if not ref_dir.exists() or not ua_dir.exists():
            continue

        for ref in sorted(ref_dir.glob("*.afa")):
            unaligned = ua_dir / ref.name
            if not unaligned.exists():
                continue

            # Count sequences and skip very large families
            n_seqs = sum(1 for line in open(ref) if line.startswith(">"))
            if n_seqs > MDSA_MAX_SEQS:
                continue
            # Skip trivial pairwise cases
            if n_seqs < 3:
                continue

            cases.append(
                BenchmarkCase(
                    family=ref.stem,
                    dataset=f"mdsa_{db}",
                    unaligned=unaligned,
                    reference=ref,
                    seq_type="dna",
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
    "mdsa": {
        "download": mdsa_download,
        "is_available": mdsa_is_available,
        "cases": mdsa_cases,
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
        for ds_name, ds in DATASETS.items():
            if not ds["is_available"]():
                try:
                    ds["download"]()
                except (urllib.error.URLError, OSError) as e:
                    print(f"WARNING: Could not download {ds_name} dataset: {e}")
                    continue
            if ds["is_available"]():
                all_cases.extend(ds["cases"]())
        if max_cases:
            all_cases = all_cases[:max_cases]
        return all_cases

    if name not in DATASETS:
        raise ValueError(f"Unknown dataset: {name}. Choose from: {list(DATASETS.keys())}")

    ds = DATASETS[name]
    if not ds["is_available"]():
        try:
            ds["download"]()
        except (urllib.error.URLError, OSError) as e:
            print(f"WARNING: Could not download {name} dataset: {e}")
            print("Skipping benchmark — dataset unavailable.")
            return []
    cases = ds["cases"]()
    if max_cases:
        cases = cases[:max_cases]
    return cases
