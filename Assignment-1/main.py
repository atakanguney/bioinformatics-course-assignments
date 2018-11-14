from copy import deepcopy


def read_scoring_matrix(filename):
    """Reads file of scoring matrix

    Parameters
    ----------
    filename: str
        File to be read

    Returns
    -------
        dict
            Scoring function as dict
    """
    with open(filename, "r") as file_:
        lines = []

        # Take uncommented lines
        for line in file_:
            if line[0] == "#":
                continue
            else:
                lines.append(line[1:].strip())

    # Make 2D array
    lines = [line.split() for line in lines]

    # Get amino acids names
    aminoacids = lines[0]

    # Get scoring matrix
    scores = lines[1:]
    scores = [[int(score) for score in score_line] for score_line in scores]

    # Construct scoring map
    scoring_function = {}
    for i, amino1 in enumerate(aminoacids):
        for j, amino2 in enumerate(aminoacids):
            scoring_function[(amino1, amino2)] = scores[i][j]

    return scoring_function


def global_sequence_alignment(sequence1, sequence2, gap_openning, gap_extention, scoring_function):
    """Needleman-Wunsh algorithm implementation for global alignment problem

    Make pairwise alignment according to Needleman-Wunsh algorithm

    Paramaters
    ----------
    sequence1: list
        List of aminoacids of first protein sequence

    sequence2: list
        List of aminoacids of second protein sequence

    gap_openning: int | float
        Gap openning penalty

    gep_extention: int | float
        Gap extention penalty

    scoring_function: dict
        Scoring function

    Returns
    -------
        score as int | float, aligned sequences as tuple, matrix constructed while running algorithm
    """

    length_1 = len(sequence1)
    length_2 = len(sequence2)

    # Initialize main-level matrix
    S = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]

    # First column initialization
    for i in range(1, length_1 + 1):
        S[i][0] = -1 * (gap_openning + (i - 1) * gap_extention)

    # First row initialization
    for j in range(1, length_2 + 1):
        S[0][j] = -1 * (gap_openning + (j - 1) * gap_extention)

    # Lower level
    S_lower = deepcopy(S)
    # Upper level
    S_upper = deepcopy(S)

    i = 1
    j = 1

    back_tracking = [[None] * length_2 for _ in range(length_1)]

    for i in range(1, length_1 + 1):
        for j in range(1, length_2 + 1):
            # Update S_upper
            S_upper[i][j] = max(S_upper[i - 1][j] -
                                gap_extention, S[i - 1][j] - gap_openning)

            # Update S_lower
            S_lower[i][j] = max(S_lower[i][j - 1] -
                                gap_extention, S[i][j - 1] - gap_openning)

            # Update S main
            S[i][j] = max(S[i - 1][j - 1] + scoring_function[(sequence1[i - 1],
                                                              sequence2[j - 1])], max(S_lower[i][j], S_upper[i][j]))

            # Save where S[i][j] is coming
            if S[i][j] == S[i - 1][j - 1] + scoring_function[(sequence1[i - 1], sequence2[j - 1])]:
                back_tracking[i - 1][j - 1] = "MAIN"
            elif S[i][j] == S_lower[i][j]:
                back_tracking[i - 1][j - 1] = "LOWER"
            elif S[i][j] == S_upper[i][j]:
                back_tracking[i - 1][j - 1] = "UPPER"

    # Back Tracking method
    def back_track(seq1, seq2, back_tracking):
        aligned_seq1 = ""
        aligned_seq2 = ""

        i = len(seq1) - 1
        j = len(seq2) - 1

        while i > -1 and j > -1:
            move = back_tracking[i][j]

            if move == "MAIN":
                aligned_seq1 = seq1.pop() + aligned_seq1
                aligned_seq2 = seq2.pop() + aligned_seq2
                i -= 1
                j -= 1
            elif move == "LOWER":
                aligned_seq1 = "_" + aligned_seq1
                aligned_seq2 = seq2.pop() + aligned_seq2
                j -= 1
            else:
                aligned_seq1 = seq1.pop() + aligned_seq1
                aligned_seq2 = "_" + aligned_seq2
                i -= 1

        if j > -1:
            while j > -1:
                aligned_seq1 = "_" + aligned_seq1
                aligned_seq2 = seq2.pop() + aligned_seq2
                j -= 1
        elif i > -1:
            while i > -1:
                aligned_seq1 = seq1.pop() + aligned_seq1
                aligned_seq2 = "_" + aligned_seq2
                i -= 1

        return aligned_seq1, aligned_seq2

    return S[-1][-1], back_track(sequence1, sequence2, back_tracking), S


def local_sequence_alignment(sequence1, sequence2, gap_openning, gap_extention, scoring_function):
    """Smith-Waterman  algorithm implementation for local alignment problem

    Make pairwise local alignment according to Smith-Waterman algorithm

    Paramaters
    ----------
    sequence1: list
        List of aminoacids of first protein sequence

    sequence2: list
        List of aminoacids of second protein sequence

    gap_openning: int | float
        Gap openning penalty

    gep_extention: int | float
        Gap extention penalty

    scoring_function: dict
        Scoring function

    Returns
    -------
        score as int | float, aligned sequences as tuple, matrix constructed while running algorithm
    """
    length_1 = len(sequence1)
    length_2 = len(sequence2)

    # Initialize main-level matrix
    S = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]

    # First column initialization
    for i in range(1, length_1 + 1):
        S[i][0] = -1 * (gap_openning + (i - 1) * gap_extention)

    # First row initialization
    for j in range(1, length_2 + 1):
        S[0][j] = -1 * (gap_openning + (j - 1) * gap_extention)

    # Lower level
    S_lower = deepcopy(S)
    # Upper level
    S_upper = deepcopy(S)

    i = 1
    j = 1

    back_tracking = [[None] * length_2 for _ in range(length_1)]

    for i in range(1, length_1 + 1):
        for j in range(1, length_2 + 1):
            # Update S_upper
            S_upper[i][j] = max(S_upper[i - 1][j] -
                                gap_extention, S[i - 1][j] - gap_openning)

            # Update S_lower
            S_lower[i][j] = max(S_lower[i][j - 1] -
                                gap_extention, S[i][j - 1] - gap_openning)

            # Update S main
            S[i][j] = max(0, max(S[i - 1][j - 1] + scoring_function[(sequence1[i - 1],
                                                                     sequence2[j - 1])], max(S_lower[i][j], S_upper[i][j])))

            # Save where S[i][j] is coming
            if S[i][j] == S[i - 1][j - 1] + scoring_function[(sequence1[i - 1], sequence2[j - 1])]:
                back_tracking[i - 1][j - 1] = "MAIN"
            elif S[i][j] == S_lower[i][j]:
                back_tracking[i - 1][j - 1] = "LOWER"
            elif S[i][j] == S_upper[i][j]:
                back_tracking[i - 1][j - 1] = "UPPER"

    # Back tracking function
    def back_track(seq1, seq2, back_tracking, maxi, maxj):
        aligned_seq1 = ""
        aligned_seq2 = ""

        i = maxi
        j = maxj

        while i > -1 and j > -1:
            move = back_tracking[i][j]

            if move == "MAIN":
                aligned_seq1 = seq1[i] + aligned_seq1
                aligned_seq2 = seq2[j] + aligned_seq2
                i -= 1
                j -= 1
            elif move == "LOWER":
                aligned_seq1 = "_" + aligned_seq1
                aligned_seq2 = seq2[j] + aligned_seq2
                j -= 1
            elif move == "UPPER":
                aligned_seq1 = seq1[i] + aligned_seq1
                aligned_seq2 = "_" + aligned_seq2
                i -= 1
            else:
                break

        return aligned_seq1, aligned_seq2

    def find_maxij(matrix):
        max_ = matrix[0][0]
        maxi = 0
        maxj = 0

        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if matrix[i][j] > max_:
                    max_ = matrix[i][j]
                    maxi = i
                    maxj = j

        return maxi - 1, maxj - 1

    maxi, maxj = find_maxij(S)

    return S[-1][-1], back_track(sequence1, sequence2, back_tracking, maxi, maxj), S


def read_pairs(filename):
    """Read pairs

    Parameters
    ----------
    filename: str
        Name of file of protein sequences

    Returns
    -------
        Protein sequences as tuple
    """
    with open(filename, "r") as file_:
        seq1 = file_.readline().strip()
        seq2 = file_.readline().strip()

    return seq1, seq2


def write_aligned_pairs(filename, sequence1, sequence2, score):
    """Writes resulst to file `filename`.txt

    Parameters
    ----------
    filename: str
        File to be written

    sequence1: list
        List of aligned sequence 1

    sequence2: list
        List of aligned sequence 2

    score: int | float
        Max aligned score
    """

    with open(filename + ".txt", "w") as file_:
        file_.write("Score: {}\n".format(score))
        file_.write("".join(sequence1) + "\n")
        file_.write("".join(sequence2) + "\n")


if __name__ == "__main__":
    SCORING_FILE = "BLOSUM62.txt"
    scoring_function = read_scoring_matrix(SCORING_FILE)

    # for amino1, amino2 in scoring_function:
    #     print("{} , {}: {}".format(amino1, amino2, scoring_function[(amino1, amino2)]))

    gap_openning = 4
    gap_extention = -1

    tests = [
        "test-pair1.txt",
        "test-pair2.txt",
        "test-pair3.txt",
    ]

    seqs = [read_pairs(test) for test in tests]
    i = 0
    for seq1, seq2 in seqs:
        score, aligneds, S = global_sequence_alignment(
            list(seq1), list(seq2), gap_openning, gap_extention, scoring_function)
        print("Global Alignment")
        print("Score: {}".format(score))
        print("Sequence 1: {}".format(aligneds[0]))
        print("Sequence 2: {}".format(aligneds[1]))
        # for line in S:
        #    print(line)

        write_aligned_pairs(
            tests[i][:-4] + "-global-result", aligneds[0], aligneds[1], score)
        score, aligneds, S = local_sequence_alignment(list(seq1), list(
            seq2), gap_openning, gap_extention, scoring_function)
        print("Local Alignment")
        print("Score: {}".format(score))
        print("Sequence 1: {}".format(aligneds[0]))
        print("Sequence 2: {}".format(aligneds[1]))
        write_aligned_pairs(tests[i][:-4] + "-local-result",
                            aligneds[0], aligneds[1], score)
        # for line in S:
        #    print(line)
        i += 1
