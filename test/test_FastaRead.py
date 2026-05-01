import pytest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))

from Motifind import Motifind


def test_fasta_basic():
    moti = Motifind()
    path = "testdata/fasta_basic.fsa"

    headers, sequences = moti.FastaRead(path)

    assert len(headers) == len(sequences)
    assert len(headers) > 0


def test_fasta_empty():
    moti = Motifind()
    path = "testdata/fasta_empty.fsa"

    headers, sequences = moti.FastaRead(path)

    assert headers == []
    assert sequences == []


def test_fasta_no_seq():
    moti = Motifind()
    path = "testdata/fasta_no_seq.fsa"

    with pytest.raises(ValueError):
        moti.FastaRead(path)


def test_fasta_no_headers():
    moti = Motifind()
    path = "testdata/fasta_no_headers.fsa"

    with pytest.raises(ValueError):
        moti.FastaRead(path)


def test_fasta_spaces():
    moti = Motifind()
    path = "testdata/fasta_spaces.fsa"

    headers, sequences = moti.FastaRead(path)

    assert len(headers) == len(sequences)
    assert all(isinstance(seq, str) for seq in sequences)


def test_fasta_blank_line():
    moti = Motifind()
    path = "testdata/fasta_blank_line.fsa"

    headers, sequences = moti.FastaRead(path)

    assert len(headers) == len(sequences)


def test_fasta_long_seq():
    moti = Motifind()
    path = "testdata/fasta_long_seq.fsa"

    headers, sequences = moti.FastaRead(path)

    assert len(headers) == len(sequences)
    assert len(sequences[0]) > 0


def test_fasta_multiline():
    moti = Motifind()
    path = "testdata/fasta_multiline.fsa"

    headers, sequences = moti.FastaRead(path)

    assert len(headers) == len(sequences)


def test_fasta_missing_file():
    moti = Motifind()

    with pytest.raises(FileNotFoundError):
        moti.FastaRead("testdata/does_not_exist.fsa")
