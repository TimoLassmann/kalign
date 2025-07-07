"""
Command-line interface wrapper for Kalign.

This module provides a CLI entry point that wraps the native Kalign binary,
ensuring the command-line interface remains available after pip installation.
"""

import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


def find_kalign_binary() -> Optional[str]:
    """
    Find the Kalign binary in the system PATH or package directory.

    Returns
    -------
    str or None
        Path to kalign binary if found, None otherwise
    """
    # First, try to find in system PATH
    try:
        result = subprocess.run(["which", "kalign"], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except (subprocess.SubprocessError, FileNotFoundError):
        pass

    # Try common installation locations
    common_paths = [
        "/usr/local/bin/kalign",
        "/usr/bin/kalign",
        "/opt/homebrew/bin/kalign",
        "/opt/local/bin/kalign",
    ]

    for path in common_paths:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path

    return None


def main(args: Optional[List[str]] = None) -> int:
    """
    Main entry point for the kalign CLI wrapper.

    Parameters
    ----------
    args : list of str, optional
        Command-line arguments. If None, uses sys.argv[1:]

    Returns
    -------
    int
        Exit code from kalign binary
    """
    if args is None:
        args = sys.argv[1:]

    # Find the kalign binary
    kalign_bin = find_kalign_binary()

    if kalign_bin is None:
        print("Error: kalign binary not found in PATH", file=sys.stderr)
        print(
            "Please ensure Kalign is properly installed on your system", file=sys.stderr
        )
        print(
            "You can download it from: https://github.com/TimoLassmann/kalign",
            file=sys.stderr,
        )
        return 1

    # Execute kalign with the provided arguments
    try:
        # Use subprocess.run to forward all arguments and maintain interactivity
        result = subprocess.run([kalign_bin] + args)
        return result.returncode

    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        return 130

    except Exception as e:
        print(f"Error executing kalign: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
