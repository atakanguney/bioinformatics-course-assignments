import os
import sys

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


# Back Tracking method
def _back_track(seq1, seq2, back_tracking, S, gap_opening, gap_extention, start=None, local=False):
    aligned_seq1 = ""
    aligned_seq2 = ""

    i = len(seq1)
    j = len(seq2)

    if start:
        i = start[0]
        j = start[1]

    while True:
        move = back_tracking[i][j]
        if move == "MAIN":
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif move == "LOWER":
            aligned_seq1 = "_" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1
        elif move == "UPPER":
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "_" + aligned_seq2
            i -= 1
        elif move == "LOWER_EXT":
            target_score = S[i][j]
            target = j
            # Find GAP opening
            for n in range(target):
                aligned_seq1 = "_" + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2
                j -= 1

                actual_score = S[i][j] - gap_opening - n * gap_extention

                if actual_score == target_score:
                    break

        elif move == "UPPER_EXT":
            target_score = S[i][j]
            target = i
            # Find GAP opening
            for n in range(target):
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = "_" + aligned_seq2
                i -= 1

                actual_score = S[i][j] - gap_opening - n * gap_extention

                if actual_score == target_score:
                    break
        else:
            if local:
                break
            # Reached to None
            while i and j:
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2
                i -= 1
                j -= 1
            while i:
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = "_" + aligned_seq2
                i -= 1
            while j:
                aligned_seq1 = "_" + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2
                j -= 1
            break

    return aligned_seq1, aligned_seq2


def global_sequence_alignment(sequence1, sequence2, gap_opening, gap_extention, scoring_function):
    """Needleman-Wunsh algorithm implementation for global alignment problem

    Make pairwise alignment according to Needleman-Wunsh algorithm

    Paramaters
    ----------
    sequence1: str
        List of aminoacids of first protein sequence

    sequence2: str
        List of aminoacids of second protein sequence

    gap_opening: int | float
        Gap openning penalty

    gep_extention: int | float
        Gap extention penalty

    scoring_function: dict
        Scoring function

    Returns
    -------
        score as int | float, aligned sequences as tuple, matrix constructed while running algorithm
    """

    # Convert strings to list
    sequence1 = list(sequence1)
    sequence2 = list(sequence2)

    # Get length of each protein
    length_1 = len(sequence1)
    length_2 = len(sequence2)

    # Initialize main-level matrix
    S = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]
    # First column initialization
    for i in range(1, length_1 + 1):
        S[i][0] = -1 * (gap_opening + (i - 1) * gap_extention)

    # First row initialization
    for j in range(1, length_2 + 1):
        S[0][j] = -1 * (gap_opening + (j - 1) * gap_extention)

    # Lower level - horizontal matrix
    S_lower = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]
    # Upper level - vertical matrix
    S_upper = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]

    # First Row initialization of horizontal matrix and vertical matrix
    for i in range(length_1 + 1):
        S_lower[i][0] = -float("inf")
        S_upper[i][0] = -float("inf")

    # First Column initialization of horizontal matrix and vertical matrix
    for j in range(length_2 + 1):
        S_lower[0][j] = -float("inf")
        S_upper[0][j] = -float("inf")

    back_tracking = [[None] * (length_2 + 1) for _ in range(length_1 + 1)]

    for i in range(1, length_1 + 1):
        for j in range(1, length_2 + 1):

            vertical_open = S[i - 1][j] - gap_opening
            vertical_extend = S_upper[i - 1][j] - gap_extention

            # Update S_upper
            S_upper[i][j] = max(vertical_open, vertical_extend)

            horizontal_open = S[i][j - 1] - gap_opening
            horizontal_extend = S_lower[i][j - 1] - gap_extention
            # Update S_lower
            S_lower[i][j] = max(horizontal_open, horizontal_extend)

            no_gap = S[i - 1][j - 1] + scoring_function[(sequence1[i - 1], sequence2[j - 1])]
            # Update S main
            S[i][j] = max(no_gap, vertical_open, vertical_extend,
                          horizontal_open, horizontal_extend)

            # Save where S[i][j] is coming
            if S[i][j] == no_gap:
                back_tracking[i][j] = "MAIN"
            if S[i][j] == vertical_open:
                back_tracking[i][j] = "UPPER"
            if S[i][j] == horizontal_open:
                back_tracking[i][j] = "LOWER"
            if S[i][j] == vertical_extend:
                back_tracking[i][j] = "UPPER_EXT"
            if S[i][j] == horizontal_extend:
                back_tracking[i][j] = "LOWER_EXT"

    return S[-1][-1], _back_track(sequence1, sequence2, back_tracking, S, gap_opening, gap_extention)


def local_sequence_alignment(sequence1, sequence2, gap_opening, gap_extention, scoring_function):
    """Smith-Waterman  algorithm implementation for local alignment problem

    Make pairwise local alignment according to Smith-Waterman algorithm

    Paramaters
    ----------
    sequence1: str
        List of aminoacids of first protein sequence

    sequence2: str
        List of aminoacids of second protein sequence

    gap_opening: int | float
        Gap openning penalty

    gep_extention: int | float
        Gap extention penalty

    scoring_function: dict
        Scoring function

    Returns
    -------
        score as int | float, aligned sequences as tuple, matrix constructed while running algorithm
    """

    # Convert strings to list
    sequence1 = list(sequence1)
    sequence2 = list(sequence2)

    # Get length of each protein sequence
    length_1 = len(sequence1)
    length_2 = len(sequence2)

    # Initialize main-level matrix
    S = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]

    # Lower level - horizontal matrix
    S_lower = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]
    S_upper = [[0] * (length_2 + 1) for _ in range(length_1 + 1)]

    # First Row initialization of horizontal matrix and vertical matrix
    for i in range(length_1 + 1):
        S_lower[i][0] = -float("inf")
        S_upper[i][0] = -float("inf")

    # Initialize back tracking matrix
    back_tracking = [[None] * (length_2 + 1) for _ in range(length_1 + 1)]

    for i in range(1, length_1 + 1):
        for j in range(1, length_2 + 1):
            vertical_open = S[i - 1][j] - gap_opening
            vertical_extend = S_upper[i - 1][j] - gap_extention

            # Update S_upper
            S_upper[i][j] = max(vertical_open, vertical_extend)

            horizontal_open = S[i][j - 1] - gap_opening
            horizontal_extend = S_lower[i][j - 1] - gap_extention

            # Update S_lower
            S_lower[i][j] = max(horizontal_open, horizontal_extend)

            no_gap = S[i - 1][j - 1] + scoring_function[(sequence1[i - 1], sequence2[j - 1])]

            # Update S main
            S[i][j] = max(0, no_gap,
                          vertical_open, vertical_extend,
                          horizontal_open, horizontal_extend)

            # Save where S[i][j] is coming
            if S[i][j] == 0:
                back_tracking[i][j] = None
            if S[i][j] == no_gap:
                back_tracking[i][j] = "MAIN"
            if S[i][j] == horizontal_open:
                back_tracking[i][j] = "LOWER"
            if S[i][j] == vertical_open:
                back_tracking[i][j] = "UPPER"
            if S[i][j] == vertical_extend:
                back_tracking[i][j] = "UPPER_EXT"
            if S[i][j] == horizontal_extend:
                back_tracking[i][j] = "LOWER_EXT"

    # Find maximum scor in the matrix
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

        return maxi, maxj, max_

    maxi, maxj, score = find_maxij(S)
    start = (maxi, maxj)

    return score, _back_track(sequence1, sequence2, back_tracking, S, gap_opening, gap_extention, start, True)


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


def write_aligned_pairs(filename, sequence1, sequence2, score, foldername=""):
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

    if foldername:
        if not os.path.exists(foldername):
            os.makedirs(foldername)
        foldername = foldername + "/"

    with open(foldername + filename + ".txt", "w") as file_:
        file_.write("Score: {}\n".format(score))
        file_.write("".join(sequence1) + "\n")
        file_.write("".join(sequence2) + "\n")


if __name__ == "__main__":
    args = sys.argv

    sequences_file = args[1]

    try:
        gap_opening = int(args[2])
        gap_extention = int(args[3])
    except ValueError:
        raise ValueError("gap opening and gap extention must be numeric")

    seq1,seq2 = read_pairs(sequences_file)
    method = args[4]

    SCORING_FILE = "BLOSUM62.txt"
    scoring_function = read_scoring_matrix(SCORING_FILE)

    if method == "global":
        score, alignments = global_sequence_alignment(seq1, seq2, gap_opening, gap_extention, scoring_function)
        print("Method: Global Alignment")
        print("Score: {}".format(score))
        print("Alignments")
        print("----------")
        print(alignments[0])
        print(alignments[1])
        write_aligned_pairs(sequences_file[:-4] + "-global-result", alignments[0], alignments[1], score, foldername="results")

    elif method == "local":
        score, alignments = local_sequence_alignment(seq1, seq2, gap_opening, gap_extention, scoring_function)
        print("Method: Local Alignment")
        print("Score: {}".format(score))
        print("Alignments")
        print("----------")
        print(alignments[0])
        print(alignments[1])
        write_aligned_pairs(sequences_file[:-4] + "-local-result", alignments[0], alignments[1], score, foldername="results")
