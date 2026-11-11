from pathlib import Path

import pytest

import kalign


def test_read_sequences_swiss_extension(tmp_path: Path) -> None:
    """The .swiss extension maps to Biopython's "swiss" format."""
    pytest.importorskip("Bio")
    path = tmp_path / "entry.swiss"
    path.write_text(
        "ID   TEST_HUMAN              Reviewed;           8 AA.\n"
        "AC   P00001;\n"
        "DE   RecName: Full=Test protein;\n"
        "SQ   SEQUENCE   8 AA;  900 MW;  0000000000000000 CRC64;\n"
        "     MKTAYIAK\n"
        "//\n"
    )
    sequences, ids = kalign.io.read_sequences(path)
    assert sequences == ["MKTAYIAK"]
    assert len(ids) == 1
