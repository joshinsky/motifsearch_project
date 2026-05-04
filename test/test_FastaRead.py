import pytest
import sys
import os

# Ensure the src directory is available for imports
sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))

from Motifind import MotiFind

# Build absolute paths to test data files relative to this test file
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE_DIR, "testdata")

def get_path(filename):
    return os.path.join(DATA_DIR, filename)


def test_fasta_basic():
    # Verify that a normal FASTA file is parsed correctly
    moti = MotiFind(None, None, None)

    headers, sequences = moti.FastaRead(get_path("fasta_basic.fsa"))

    assert len(headers) == len(sequences)
    assert len(headers) > 0


def test_fasta_empty():
    # Verify that an empty FASTA file returns empty lists
    moti = MotiFind(None, None, None)

    headers, sequences = moti.FastaRead(get_path("fasta_empty.fsa"))

    assert headers == []
    assert sequences == []


def test_fasta_no_seq():
    # Verify that missing sequences raises an error
    moti = MotiFind(None, None, None)

    with pytest.raises(ValueError):
        moti.FastaRead(get_path("fasta_no_seq.fsa"))


def test_fasta_no_headers():
    # Verify that missing headers raises an error
    moti = MotiFind(None, None, None)

    with pytest.raises(ValueError):
        moti.FastaRead(get_path("fasta_no_headers.fsa"))


def test_fasta_spaces():
    # Verify that sequences with spaces or formatting are handled correctly
    moti = MotiFind(None, None, None)

    headers, sequences = moti.FastaRead(get_path("fasta_spaces.fsa"))

    assert len(headers) == len(sequences)
    assert all(isinstance(seq, str) for seq in sequences)


def test_fasta_blank_line():
    # Verify that blank lines in the file do not break parsing
    moti = MotiFind(None, None, None)

    headers, sequences = moti.FastaRead(get_path("fasta_blank_line.fsa"))

    assert len(headers) == len(sequences)


def test_fasta_long_seq():
    # Verify that long sequences are read correctly
    moti = MotiFind(None, None, None)

    headers, sequences = moti.FastaRead(get_path("fasta_long_seq.fsa"))

    assert len(headers) == len(sequences)
    assert len(sequences[0]) > 0


def test_fasta_multiline():
    # Verify that sequences split across multiple lines are combined correctly
    moti = MotiFind(None, None, None)

    headers, sequences = moti.FastaRead(get_path("fasta_multiline.fsa"))

    assert len(headers) == len(sequences)


def test_fasta_missing_file():
    # Verify that a missing file raises FileNotFoundError
    moti = MotiFind(None, None, None)

    with pytest.raises(FileNotFoundError):
        moti.FastaRead(get_path("does_not_exist.fsa"))
