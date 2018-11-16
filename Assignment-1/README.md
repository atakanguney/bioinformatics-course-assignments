# Assignmet - 1 - Pairwise Sequence Alignment with Affine Gap Penalty

## Requirements

-   python==3.6.6
-   biopython==1.72 to run `test.py`

For detailed list of my conda environment packages is [here](requirements.txt).

## How to run program

-   Usage of alignments
    - Template
        - `python PairwiseAlignment.py "sequence-filename" "gap-opening" "gap-extention" "method"`
    -   Using script to get output, on the command line, run the command
        - Example usage
        - `python PairwiseAlignment.py test-pair1.txt 11 1 global`
        - `python PairwiseAlignment.py test-pair1.txt 11 1 local`

    -   Example output
        - 
        ```
            Method: Global Alignment
            Score: -1
            Alignments
            ----------
            PLEASANTLY
            _MEAN___LY
        ```
    -   It also creates a folder named "results" and saves the results as txt:
        -   For example, assuming your input sequence file named as "test-pair1.txt", program saves results as "test-pair1-global-result.txt" or "test-pair1-local-result.txt" with respect to the method that have been used.

## Using methods in your python scripy
In `PairwiseAlignment.py`, there are 2 methods for alignments:

-   For Global Alignment with affinity gap penalties:
    -   `global_sequence_alignment(sequence1, sequence2, gap_opening, gap_extention, scoring_function)`
    This function is an implementation of Needleman-Wunsh algorithm with affine gap penalties
    - `scoring_function` is a dictionary that contains aminoacid pairs as key and their similarities as values, in main program it is read from "BLOSUM62.txt" file, it can be changed with another similarity dictionary.

- For Local Alignment with affinity gap penalties:
    - `local_sequence_alignment(sequence1, sequence2, gap_opening, gap_extention, scoring_function)`
    This function is an implementation of Smith-Waterman algorithm with affine gap penalties.
    - `scoring_function` is a dictionary that contains aminoacid pairs as key and their similarities as values, in main program it is read from "BLOSUM62.txt" file, it can be changed with another similarity dictionary.

    - To use functions directly in your python script, assuming your working directory is same as the file `PairwiseAlignmeny.py` file:
        - `>>> from PairwiseAlignment import global_sequence_alignment, read_scoring_matrix`
        - `>>> scoring_function = read_scoring_matrix("BLOSUM62.txt")`
        - `>>> score, alignments = global_sequence_alignment("PLEASANTLY", "MEANLY", 11, 1, scoring_function)`

    - Returns
        - `>>> score`
        - `-1`
        - `>>> alignments`
        - `('PLEASANTLY', '_MEAN___LY')`

-   For detailed usage [see](test.py).
-   To get results of given pairs "test-pair1.txt", "test-pair2.txt", "test-pair3.txt", run test.py
-   Note: You need to install `biopython` package to run `test.py`
