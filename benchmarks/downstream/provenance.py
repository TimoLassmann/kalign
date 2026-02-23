"""Version collection, hardware info, and result file naming for reproducible benchmarks.

Captures full provenance metadata (tool versions, hardware specs, git state)
and provides deterministic file naming for benchmark result files.
"""

import json
import os
import platform
import re
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional


@dataclass
class Provenance:
    """Full provenance record for a benchmark run."""

    timestamp: str
    kalign_version: str
    kalign_commit: str
    container_image: Optional[str]
    hostname: str
    cpu_model: str
    cpu_cores: int
    ram_gb: float
    os_version: str
    python_version: str
    tool_versions: Dict[str, str]
    parameters: Dict[str, object]


# ---------------------------------------------------------------------------
# Tool version collection
# ---------------------------------------------------------------------------

# Each entry: (binary, args, regex to extract version string)
_TOOL_VERSION_COMMANDS = {
    "kalign": ("kalign", ["-v"], r"(\d+\.\d+\.\d+\S*)"),
    "mafft": ("mafft", ["--version"], r"v?([\d.]+\S*)"),
    "muscle": ("muscle", ["--version"], r"([\d.]+\S*)"),
    "clustalo": ("clustalo", ["--version"], r"([\d.]+\S*)"),
    "hmmer": ("hmmbuild", ["-h"], r"HMMER\s+([\d.]+\S*)"),
    "iqtree": ("iqtree2", ["--version"], r"version\s+([\d.]+\S*)"),
    "hyphy": ("hyphy", ["--version"], r"([\d.]+\S*)"),
    "indelible": ("indelible", ["--help"], r"INDELible\s+[Vv]?([\d.]+\S*)"),
    "guidance2": ("guidance.pl", ["--version"], r"([\d.]+\S*)"),
}


def _run_version_command(binary: str, args: list, pattern: str) -> str:
    """Run a single version command, return parsed version or 'not found'."""
    try:
        result = subprocess.run(
            [binary] + args,
            capture_output=True,
            text=True,
            timeout=5,
        )
        # Some tools print version on stderr (e.g. mafft)
        combined = result.stdout + "\n" + result.stderr
        match = re.search(pattern, combined)
        if match:
            return match.group(1)
        # Command ran but regex didn't match -- return raw first line
        first_line = combined.strip().split("\n")[0].strip()
        return first_line if first_line else "unknown"
    except FileNotFoundError:
        return "not found"
    except subprocess.TimeoutExpired:
        return "timeout"
    except OSError:
        return "not found"


def collect_tool_versions() -> Dict[str, str]:
    """Collect version strings for all known external tools.

    Tools that are not installed are reported as ``"not found"`` rather
    than raising an error, so callers can always serialise the result.
    """
    versions: Dict[str, str] = {}
    for name, (binary, args, pattern) in _TOOL_VERSION_COMMANDS.items():
        versions[name] = _run_version_command(binary, args, pattern)
    return versions


# ---------------------------------------------------------------------------
# Hardware / environment helpers
# ---------------------------------------------------------------------------

def _get_cpu_model() -> str:
    """Return a human-readable CPU model string (macOS + Linux)."""
    system = platform.system()
    if system == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True, text=True, timeout=5,
            )
            model = result.stdout.strip()
            if model:
                return model
        except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
            pass
    elif system == "Linux":
        try:
            with open("/proc/cpuinfo") as fh:
                for line in fh:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
        except OSError:
            pass
    return platform.processor() or "unknown"


def _get_ram_gb() -> float:
    """Return total physical RAM in GiB, rounded to one decimal."""
    system = platform.system()
    if system == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True, text=True, timeout=5,
            )
            return round(int(result.stdout.strip()) / (1024 ** 3), 1)
        except (FileNotFoundError, subprocess.TimeoutExpired, OSError, ValueError):
            pass
    elif system == "Linux":
        try:
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemTotal"):
                        kb = int(line.split()[1])
                        return round(kb / (1024 ** 2), 1)
        except (OSError, ValueError):
            pass
    return 0.0


def _get_kalign_version() -> str:
    """Return the installed kalign Python package version."""
    try:
        import kalign as _kalign
        return getattr(_kalign, "__version__", "unknown")
    except ImportError:
        return "not installed"


def _get_kalign_commit() -> str:
    """Return the short git commit hash of the kalign repo."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
            cwd=Path(__file__).resolve().parents[2],  # repo root
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
        pass
    return "unknown"


def _get_container_image() -> Optional[str]:
    """Detect container image from environment if running inside one."""
    # Common env vars set by container runtimes / Dockerfiles
    for var in ("CONTAINER_IMAGE", "IMAGE_NAME"):
        value = os.environ.get(var)
        if value:
            return value
    # Check cgroup for docker/podman hints
    try:
        with open("/proc/1/cgroup") as fh:
            text = fh.read()
            if "docker" in text or "podman" in text or "containerd" in text:
                return "unknown-container"
    except OSError:
        pass
    return None


# ---------------------------------------------------------------------------
# Main provenance collector
# ---------------------------------------------------------------------------

def collect_provenance(parameters: dict) -> Provenance:
    """Collect all provenance information for the current environment.

    Parameters
    ----------
    parameters : dict
        Benchmark-specific parameters to record (e.g. dataset name,
        number of threads, alignment options).
    """
    now = datetime.now(timezone.utc)
    return Provenance(
        timestamp=now.strftime("%Y-%m-%dT%H:%M:%S%z"),
        kalign_version=_get_kalign_version(),
        kalign_commit=_get_kalign_commit(),
        container_image=_get_container_image(),
        hostname=platform.node(),
        cpu_model=_get_cpu_model(),
        cpu_cores=os.cpu_count() or 0,
        ram_gb=_get_ram_gb(),
        os_version=f"{platform.system()} {platform.release()}",
        python_version=platform.python_version(),
        tool_versions=collect_tool_versions(),
        parameters=parameters,
    )


# ---------------------------------------------------------------------------
# Result file naming and symlink management
# ---------------------------------------------------------------------------

def result_path(results_dir: Path, pipeline: str) -> Path:
    """Generate a unique result file path.

    Returns a path of the form::

        results_dir/<pipeline>/run_YYYYMMDD_HHMMSS_<commit>.json

    The directory is created if it does not exist.
    """
    now = datetime.now(timezone.utc)
    stamp = now.strftime("%Y%m%d_%H%M%S")
    commit = _get_kalign_commit()
    name = f"run_{stamp}_{commit}.json"

    out_dir = Path(results_dir) / pipeline
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / name


def update_latest_symlink(result_file: Path) -> None:
    """Create or update a ``latest.json`` symlink next to *result_file*.

    The symlink points to *result_file* using a relative target so that
    the whole results tree can be moved without breaking the link.
    """
    result_file = Path(result_file).resolve()
    link_path = result_file.parent / "latest.json"

    # Relative target so the symlink is portable
    target = result_file.name

    # Atomic replace: unlink then symlink (no rename for symlinks)
    try:
        link_path.unlink()
    except FileNotFoundError:
        pass
    link_path.symlink_to(target)


def load_latest_results(results_dir: Path) -> dict:
    """Load the JSON data from the ``latest.json`` symlink.

    Parameters
    ----------
    results_dir : Path
        Directory that contains the ``latest.json`` symlink (typically
        a pipeline sub-directory under the main results directory).

    Returns
    -------
    dict
        Parsed JSON contents of the latest result file.

    Raises
    ------
    FileNotFoundError
        If no ``latest.json`` symlink exists in *results_dir*.
    """
    link_path = Path(results_dir) / "latest.json"
    with open(link_path) as fh:
        return json.load(fh)
