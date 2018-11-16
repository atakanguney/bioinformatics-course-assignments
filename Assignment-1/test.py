from Bio import pairwise2
from PairwiseAlignment import read_scoring_matrix
from PairwiseAlignment import global_sequence_alignment
from PairwiseAlignment import local_sequence_alignment
from PairwiseAlignment import write_aligned_pairs
from PairwiseAlignment import read_pairs

def test_score(seq1, seq2, score, gap_opening, gap_extention, scoring_function, local=False):
    if local:
        true_alignments = pairwise2.align.localds(
            seq1, seq2, scoring_function, -1 * gap_opening, -1 * gap_extention)
    else:
        true_alignments = pairwise2.align.globalds(
            seq1, seq2, scoring_function, -1 * gap_opening, -1 * gap_extention)

    for _, _, true_score, _, _ in true_alignments:
        if true_score == score:
            return True

    return False

SCORING_FILE = "BLOSUM62.txt"
scoring_function = read_scoring_matrix(SCORING_FILE)

gap_opening = 11
gap_extention = 1

tests = [
    "test-pair1.txt",
    "test-pair2.txt",
    "test-pair3.txt",
]

seqs = [read_pairs(test) for test in tests]
i = 0
for seq1, seq2 in seqs:
    score, aligneds = global_sequence_alignment(seq1, seq2, gap_opening, gap_extention, scoring_function)

    print("Global Alignment")
    print("Score: {}".format(score))
    print("Sequence 1: {}".format(aligneds[0]))
    print("Sequence 2: {}".format(aligneds[1]))

    # Test score with BioPythons score
    assert test_score(seq1, seq2, score, gap_opening, gap_extention, scoring_function), "Wrong Global Score !"

    for k, a in enumerate(
            pairwise2.align.globalds(seq1, seq2, scoring_function, -1 * gap_opening, -1 * gap_extention)):
        write_aligned_pairs(
            tests[i][:-4] + "-global-result-ground-truth-{}".format(k + 1), a[0], a[1], a[2],
            foldername="test-ground-truth-results")

    write_aligned_pairs(
        tests[i][:-4] + "-global-result", aligneds[0], aligneds[1], score, foldername="test-results")

    score, aligneds = local_sequence_alignment(seq1, seq2, gap_opening, gap_extention, scoring_function)

    print("Local Alignment")
    print("Score: {}".format(score))
    print("Sequence 1: {}".format(aligneds[0]))
    print("Sequence 2: {}".format(aligneds[1]))

    assert test_score(seq1, seq2, score, gap_opening, gap_extention, scoring_function,
                        local=True), "Wrong Local Score !"

    for k, a in enumerate(
            pairwise2.align.localds(seq1, seq2, scoring_function, -1 * gap_opening, -1 * gap_extention)):
        write_aligned_pairs(
            tests[i][:-4] + "-local-result-ground-truth-{}".format(k + 1), a[0], a[1], a[2],
            foldername="test-ground-truth-results")

    write_aligned_pairs(tests[i][:-4] + "-local-result",
                        aligneds[0], aligneds[1], score, foldername="test-results")
    i += 1
