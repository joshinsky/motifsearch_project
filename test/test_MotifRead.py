import pytest
import sys
import os
from pathlib import Path

#sys.path.append(os.path.join(os.path.dirname(__file__), "../src"))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))


from Motifind import MotiFind

TESTDATA = Path(__file__).resolve().parents[1] / "testdata"


# define test file names
@pytest.mark.parametrize("filename, expectation",[
	(TESTDATA / "motif_basic.txt", 					[('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'T', 'A'}, 5.0), ('char', {'T'}, 8.0)]),
	(TESTDATA / "motif_empty.txt", 					[]),
	(TESTDATA / "motif_extremely_long_motif.txt", 	[('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0), ('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0)]),
	(TESTDATA / "motif_gap_len_zero.txt", 			[('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 0, 0), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0)]),
	(TESTDATA / "motif_missing_bases.txt", 			IndexError),
	(TESTDATA / "motif_missing_file.txt", 			FileNotFoundError),
	(TESTDATA / "motif_missing_penalty.txt", 		IndexError),
	(TESTDATA / "motif_mixed_upper_lower.txt", 		[('char', {'g'}, 3.0), ('char', {'c'}, 3.0), ('char', {'c'}, 3.0), ('char', {'g'}, 6.0), ('char', {'c'}, 6.0), ('char', {'c'}, 6.0), ('char', {'g', 'A'}, 7.0), ('char', {'c'}, 6.0), ('char', {'c'}, 6.0), ('char', {'A'}, 9.0), ('char', {'U'}, 9.0), ('char', {'G'}, 9.0), ('char', {'G'}, 9.0)]),
	(TESTDATA / "motif_negative_gap_size.txt", 		ValueError),
	(TESTDATA / "motif_negative_penalties.txt", 	ValueError),
	(TESTDATA / "motif_one_line_with_many_tabs.txt", [('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'T', 'A'}, 5.0), ('char', {'T'}, 8.0)]),
	(TESTDATA / "motif_only_comments.txt", 			[]),
	(TESTDATA / "motif_only_gap.txt", 				[('gap', 5, 10)]),
	(TESTDATA / "motif_penalties_zero.txt", 		[('char', {'T'}, 0.0), ('char', {'T'}, 0.0), ('char', {'G'}, 0.0), ('char', {'A'}, 0.0), ('char', {'C'}, 0.0), ('char', {'A'}, 0.0), ('gap', 2, 2), ('char', {'T'}, 0.0), ('char', {'A'}, 0.0), ('char', {'T'}, 0.0), ('char', {'A'}, 0.0), ('char', {'T', 'A'}, 0.0), ('char', {'T'}, 0.0)]),
	(TESTDATA / "motif_random_linebreaks.txt", 		IndexError),
	(TESTDATA / "motif_space_separated.txt", 		IndexError),
	(TESTDATA / "motif_space_then_tab.txt", 		[('char', {'G'}, 10.0), ('char', {'*'}, 4.0), ('char', {'G'}, 10.0), ('char', {'K'}, 10.0), ('char', {'T', 'S'}, 10.0)]),
	(TESTDATA / "motif_starts_with_gap.txt", 		[('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'T', 'A'}, 5.0), ('char', {'T'}, 8.0)]),
	(TESTDATA / "motif_strange_motif.txt", 			[('char', {'c', 'h', 'I'}, 4.0), ('char', {'n', 'B', 'i'}, 3.0), ('char', {'n', 'E', 'i'}, 4.0), ('char', {'M', 'i', 'o', 'v', 't'}, 2.0), ('char', {'F', 'r', 'ü'}, 2.0), ('char', {'e', 'n', 'E', 'i'}, 2.0), ('char', {'a', 'g', 'M', 'n', 'i', 'o', 's', 'r', 'k', 'm', 'u'}, 1.0), ('char', {'f', 'H', 'i', '!', 'l', 'e'}, 3.0)]),
	(TESTDATA / "motif_two_following_gaps.txt", 	[('char', {'G'}, 10.0), ('gap', 4, 4), ('gap', 5, 5), ('char', {'G'}, 10.0), ('char', {'K'}, 10.0), ('char', {'T', 'S'}, 10.0)]),
	(TESTDATA / "motif_with_comma.txt", 			ValueError),
	(TESTDATA / "motif_with_float_penalties.txt", 	[('char', {'T'}, 7.4), ('char', {'T'}, 8.1), ('char', {'G'}, 6.2), ('char', {'A'}, 5.0), ('char', {'C'}, 5.6), ('char', {'A'}, 5.7), ('gap', 15, 21), ('char', {'T'}, 8.2), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0)]),
	(TESTDATA / "motif_with_space_instead_minus.txt", ValueError),
	(TESTDATA / "motif_with_underscores.txt", 		[('char', {'T'}, 7.0), ('char', {'T'}, 8.0), ('char', {'G'}, 6.0), ('char', {'A'}, 5.0), ('char', {'C'}, 5.0), ('char', {'A'}, 5.0), ('gap', 15, 21), ('char', {'T'}, 8.0), ('char', {'A'}, 8.0), ('char', {'T'}, 6.0), ('char', {'A'}, 6.0), ('char', {'A', 'T'}, 5.0), ('char', {'T'}, 8.0)]),
	(TESTDATA / "motif_with_written_penalties.txt", ValueError)
	])

def test_MotifRead(filename, expectation):

	mf = MotiFind(fasta_name="", motif_name="", max_penalty=0)

	# case 1: we expect an error
	if isinstance(expectation, type) and issubclass(expectation, Exception):
		with pytest.raises(expectation):
			motif = mf.MotifRead(filename)

	# case 2: we expect a valid output
	else:
		motif = mf.MotifRead(filename)
		assert motif == expectation



