from Bio import SeqIO
import numpy as np
import pandas as pd


def read_score_matrix(file_path):
    # Read matrix as a DataFrame, ignoring lines starting with "#"
    df = pd.read_csv(file_path, sep="\s+", comment="#", index_col=0)
    return df


def parse_fasta(file_path):
    # Parse sequences and identifiers from a FASTA file
    records = list(SeqIO.parse(file_path, "fasta"))
    return records[0].id, str(records[0].seq), records[1].id, str(records[1].seq)


def global_alignment(seq1, seq2, score_matrix, gap_penalty):
    # Initialize scoring and traceback matrices
    n, m = len(seq1), len(seq2)
    scores = np.zeros((n + 1, m + 1), dtype=int)
    traceback = np.zeros((n + 1, m + 1), dtype=int)

    # Fill the first row and column with gap penalties
    for i in range(1, n + 1):
        scores[i, 0] = i * gap_penalty
        traceback[i, 0] = 1  # Indicates a gap in seq2
    for j in range(1, m + 1):
        scores[0, j] = j * gap_penalty
        traceback[0, j] = 2  # Indicates a gap in seq1

    # Populate the scoring and traceback matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = (
                scores[i - 1, j - 1] + score_matrix.loc[seq1[i - 1], seq2[j - 1]]
            )
            delete_score = scores[i - 1, j] + gap_penalty
            insert_score = scores[i, j - 1] + gap_penalty

            # Choose the maximum score for the current cell
            scores[i, j] = max(match_score, delete_score, insert_score)

            # Set traceback direction based on the chosen score
            if scores[i, j] == match_score:
                traceback[i, j] = 0
            elif scores[i, j] == delete_score:
                traceback[i, j] = 1
            else:
                traceback[i, j] = 2

    # Traceback to reconstruct the alignment
    align1, align2 = "", ""
    i, j = n, m
    while i > 0 or j > 0:
        if traceback[i, j] == 0:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif traceback[i, j] == 1:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return [(align1, align2)], scores[n, m]


def local_alignment(seq1, seq2, score_matrix, gap_penalty):
    # Initialize scoring and traceback matrices
    n, m = len(seq1), len(seq2)
    scores = np.zeros((n + 1, m + 1), dtype=int)
    traceback = np.zeros((n + 1, m + 1), dtype=int)
    max_score = 0
    max_positions = []

    # Fill the scoring and traceback matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = (
                scores[i - 1, j - 1] + score_matrix.loc[seq1[i - 1], seq2[j - 1]]
            )
            delete_score = scores[i - 1, j] + gap_penalty
            insert_score = scores[i, j - 1] + gap_penalty

            # Use max(0, ...) to allow for local alignment
            scores[i, j] = max(0, match_score, delete_score, insert_score)

            # Update the max score and its positions
            if scores[i, j] > max_score:
                max_score = scores[i, j]
                max_positions = [(i, j)]
            elif scores[i, j] == max_score:
                max_positions.append((i, j))

            # Set traceback direction based on the chosen score
            if scores[i, j] == match_score:
                traceback[i, j] = 0
            elif scores[i, j] == delete_score:
                traceback[i, j] = 1
            elif scores[i, j] == insert_score:
                traceback[i, j] = 2

    # Traceback to reconstruct all alignments starting from max positions
    alignments = []
    for start_i, start_j in max_positions:
        align1, align2 = "", ""
        i, j = start_i, start_j
        while i > 0 and j > 0 and scores[i, j] >= 0:
            if traceback[i, j] == 0:
                align1 = seq1[i - 1] + align1
                align2 = seq2[j - 1] + align2
                i -= 1
                j -= 1
            elif traceback[i, j] == 1:
                align1 = seq1[i - 1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j - 1] + align2
                j -= 1

        alignments.append((align1, align2))

    # Sort alignments by sequence order
    alignments = sorted(alignments, key=lambda x: (x[0], x[1]))

    return alignments, max_score


def alignment(input_path, score_path, output_path, aln, gap):
    # Read the scoring matrix and sequences from input files
    score_matrix = read_score_matrix(score_path)
    id1, seq1, id2, seq2 = parse_fasta(input_path)

    # Choose alignment type and calculate alignment
    if aln == "global":
        alignments, score = global_alignment(seq1, seq2, score_matrix, gap)
    elif aln == "local":
        alignments, score = local_alignment(seq1, seq2, score_matrix, gap)

    # Write the results to the output file in FASTA format
    with open(output_path, "w") as f:
        for align1, align2 in alignments:
            f.write(f">{id1}\n{align1}\n")
            f.write(f">{id2}\n{align2}\n")
