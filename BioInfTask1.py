import re
import numpy as np
import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

def parse_sequences():
    """
    Parses the input from the user: either reads two sequences from a FASTA file,
    or takes two sequences as command-line arguments. Validates DNA sequences.

    Returns:
        tuple: Two uppercased DNA sequences (seq_a, seq_b).

    Raises:
        ValueError: If any of the arguments are invalid.
    """
    parser = argparse.ArgumentParser(description="Compare two sequences from input or a file")
    parser.add_argument("--file", type=str, help="Path to file containing sequences")
    parser.add_argument("-a", type=str, help="First sequence")
    parser.add_argument("-b", type=str, help="Second sequence")
    args = parser.parse_args()

    if args.file:
        sequences = []
        try:
            for record in SeqIO.parse(args.file, "fasta"):
                sequences.append(str(record.seq))
                if len(sequences) >= 2:
                    break
            if len(sequences) < 2:
                raise ValueError("FASTA file must contain at least two sequences.")
            a, b = sequences[:2]
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
    else:
        if not args.a or not args.b:
            raise ValueError("Both -a and -b must be provided when --file is not used.")
        a, b = args.a, args.b

        if not is_valid_dna(a) or not is_valid_dna(b):
            raise ValueError("Sequences must contain only A, C, G, or T characters.")

    return a.upper(), b.upper()


def is_valid_dna(seq):
    """Checks if the sequence contains only valid DNA characters: A, C, G, T.
    Arguments:
        seq (str): The DNA sequence to validate.

    Returns:
        bool: True if the sequence contains only A, C, G, or T (case-insensitive), False otherwise.
    "
    """
    return bool(re.fullmatch(r"[ACGTacgt]+", seq))

def generate_alignment_line(a, b):
    """
    Generates a string representing the visual alignment between two sequences.
    Matches are marked with '|', mismatches with '-', and gaps with ';'.

    Arguments:
        a (str): First aligned sequence.
        b (str): Second aligned sequence.

    Returns:
        str: Alignment visualization line.
    """
    alignment_line = []
    for x, y in zip(a, b):
        if x == y:
            alignment_line.append('|')  # Match
        elif x == '-' or y == '-':
            alignment_line.append(';')  # Gap
        else:
            alignment_line.append('-')  # Mismatch
    return ''.join(alignment_line)


def initialize_matrix(len_a, len_b, gap_penalty):
    """
    Initializes the score matrix with gap penalties applied to the first row and column.

    Arguments:
        len_a (int): Length of the first sequence.
        len_b (int): Length of the second sequence.
        gap_penalty (int): Penalty score for introducing gaps.

    Returns:
        numpy array: Initialized scoring matrix.
    """
    matrix = np.zeros((len_a + 1, len_b + 1))

    for i in range(1, len_a + 1):
        matrix[i][0] = matrix[i - 1][0] + gap_penalty
    for j in range(1, len_b + 1):
        matrix[0][j] = matrix[0][j - 1] + gap_penalty

    return matrix


def compute_scores(a, b, matrix, match_score, mismatch_score, gap_penalty):
    """
    Fills the scoring matrix based on matches, mismatches, and gaps according
    to the Needleman-Wunsch algorithm.

    This part is highly inspired by the example available on Wikipedia.
    https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm

    Arguments:
        a (str): First sequence.
        b (str): Second sequence.
        matrix (numpy array): Scoring matrix to update.
        match_score (int): Score for a match.
        mismatch_score (int): Score for a mismatch.
        gap_penalty (int): Penalty for a gap.
    """
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                score = match_score  # Match: score = +1
            else:
                score = mismatch_score  # Mismatch: score = 0

            matrix[i][j] = max(
                matrix[i - 1][j] + gap_penalty,  # Vertical gap
                matrix[i][j - 1] + gap_penalty,  # Horizontal gap
                matrix[i - 1][j - 1] + score  # Diagonal (match/mismatch)
            )


def traceback(a, b, matrix, match_score, mismatch_score, gap_penalty):
    """
    Traces back through the scoring matrix to reconstruct the best alignment path between the sequences.
    This one is also inspired by the example available on Wikipedia.

    Arguments:
        a (str): First sequence.
        b (str): Second sequence.
        matrix (numpy array): Completed scoring matrix.
        match_score (int): Score for a match.
        mismatch_score (int): Score for a mismatch.
        gap_penalty (int): Penalty for a gap.

    Returns:
        tuple: Aligned sequences (aligned_a, aligned_b) and list of path coordinates.
    """
    aligned_a, aligned_b = "", ""
    path_coords = []
    i, j = len(a), len(b)

    while i > 0 or j > 0:
        path_coords.append((i, j))
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + (
        match_score if a[i - 1] == b[j - 1] else mismatch_score):
            aligned_a = a[i - 1] + aligned_a
            aligned_b = b[j - 1] + aligned_b
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            aligned_a = a[i - 1] + aligned_a
            aligned_b = "-" + aligned_b
            i -= 1
        else:
            aligned_a = "-" + aligned_a
            aligned_b = b[j - 1] + aligned_b
            j -= 1

    path_coords.append((0, 0))
    path_coords.reverse()

    return aligned_a, aligned_b, path_coords


def needleman_wunsch(a, b, match_score=1, mismatch_score=0, gap_penalty=-1):
    """
    Coordinates the overall Needleman-Wunsch alignment process by calling matrix initialization,
    score computation, and traceback.

    Arguments:
        a (str): First sequence.
        b (str): Second sequence.
        match_score (int, optional): Score for a match. Defaults to 1.
        mismatch_score (int, optional): Score for a mismatch. Defaults to 0.
        gap_penalty (int, optional): Penalty for a gap. Defaults to -1.

    Returns:
        tuple: Aligned sequences, match count, gap count, identity percentage,
        score matrix, and path coordinates.
    """
    # Initialize scoring matrix
    matrix = initialize_matrix(len(a), len(b), gap_penalty)

    # Compute the matrix with scores
    compute_scores(a, b, matrix, match_score, mismatch_score, gap_penalty)

    # Perform traceback to get the optimal alignment
    aligned_a, aligned_b, path_coords = traceback(a, b, matrix, match_score, mismatch_score, gap_penalty)

    # Calculate match count, gap count, and identity
    match_count = sum(1 for x, y in zip(aligned_a, aligned_b) if x == y)
    gap_count = aligned_a.count('-') + aligned_b.count('-')
    identity = (match_count / len(aligned_a)) * 100

    return aligned_a, aligned_b, match_count, gap_count, identity, matrix, path_coords


def plot_alignment_matrix(matrix, coordinates, seq_a, seq_b, image_filename):
    """
    Generates a visual heatmap of the alignment score matrix with matplotlib, and plots the optimal alignment path.

    Arguments:
        matrix (numpy array): Scoring matrix.
        coordinates (list): List of path coordinates.
        seq_a (str): First sequence.
        seq_b (str): Second sequence.
        image_filename (str): Filename to save the generated plot.

    """
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.imshow(matrix, cmap='Blues')

    # Create a grid and label it with the scores
    ax.set_xticks(np.arange(len(seq_b) + 1))
    ax.set_yticks(np.arange(len(seq_a) + 1))

    ax.set_xticklabels([' '] + list(seq_b))
    ax.set_yticklabels([' '] + list(seq_a))

    # Annotate the matrix with the cell values
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            ax.text(j, i, str(int(matrix[i, j])), ha='center', va='center', fontsize=8)

    # Highlight the optimal alignment path
    path_y, path_x = zip(*coordinates)
    ax.plot(path_x, path_y, marker='x', color='mediumblue', linestyle='-', linewidth=2.5, markersize=5, alpha=0.7)

    # Set labels and title
    ax.set_xlabel("Sequence B")
    ax.set_ylabel("Sequence A")
    ax.set_title("Needleman-Wunsch Alignment Matrix with Path")

    # Tight layout to avoid overlap
    plt.tight_layout()

    # Save the plot
    try:
        plt.savefig(image_filename, dpi=300)
        plt.close()
    except Exception as e:
        print(f"Error saving image: {e}")


def save_to_pdf(seq_a, seq_b, aligned_a, aligned_b, match, gap, identity, match_score, mismatch_score, gap_penalty,
                graph_path):
    """
    Creates a PDF report summarizing the alignment parameters, sequences, alignment results, and an embedded matrix plot.

    Arguments:
        seq_a (str): First input sequence.
        seq_b (str): Second input sequence.
        aligned_a (str): First aligned sequence.
        aligned_b (str): Second aligned sequence.
        match (int): Total number of matches.
        gap (int): Total number of gaps.
        identity (float): Identity percentage.
        match_score (int): Match score.
        mismatch_score (int): Mismatch score.
        gap_penalty (int): Gap penalty.
        graph_path (str): Path to the alignment matrix image.

    """
    pdf_filename = "seq_comparison_report.pdf"
    c = canvas.Canvas(pdf_filename, pagesize=letter)
    width, height = letter
    margin = 50
    y = height - margin

    # Title of the PDF
    c.setFont("Helvetica-Bold", 16)
    c.drawString(100, y, "Needleman-Wunsch Sequence Comparison Report")
    y -= 50

    # Alignment Parameters Section
    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Alignment Parameters")
    y -= 20
    c.setFont("Helvetica", 12)
    c.drawString(100, y, f"• Match score:        {match_score}")
    y -= 20
    c.drawString(100, y, f"• Mismatch penalty:   {mismatch_score}")
    y -= 20
    c.drawString(100, y, f"• Gap penalty:        {gap_penalty}")
    y -= 30

    # Sequences Section
    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Sequences")
    y -= 20
    c.setFont("Helvetica", 12)
    c.drawString(100, y, f"> Sequence A: {seq_a}")
    y -= 20
    c.drawString(100, y, f"> Sequence B: {seq_b}")
    y -= 30

    # Alignment Summary Section
    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Alignment Summary")
    y -= 20
    c.setFont("Helvetica", 12)
    c.drawString(100, y, f"• Aligned length:     {len(aligned_a)}")
    y -= 20
    c.drawString(100, y, f"• Matches:            {match}")
    y -= 20
    c.drawString(100, y, f"• Gaps:               {gap}")
    y -= 20
    c.drawString(100, y, f"• Identity:           {identity:.2f}%")
    y -= 30

    # Alignment Representation Section
    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Alignment Representation")
    y -= 20
    max_line_width = 80
    alignment_line = generate_alignment_line(aligned_a, aligned_b)

    for i in range(0, len(aligned_a), max_line_width):
        c.drawString(100, y, aligned_a[i:i + max_line_width])
        y -= 20
        c.drawString(100, y, alignment_line[i:i + max_line_width])
        y -= 20
        c.drawString(100, y, aligned_b[i:i + max_line_width])
        y -= 30

    # Add graph to the PDF if the path is provided
    if graph_path:
        try:
            c.drawImage(graph_path, 100, y - 300, width=400, height=300)  # Adjust position and size
        except Exception as e:
            print(f"Error adding image: {e}")

    # Save the PDF
    c.save()
    print(f"PDF saved as {pdf_filename}")


if __name__ == "__main__":
    """
    Connects all the pieces: reads input, performs alignment, generates plot, and saves the final report as a PDF.
    """
    seq_a, seq_b = parse_sequences()
    if seq_a and seq_b:

        aligned_a, aligned_b, match, gap, identity, score_matrix, path_coords = needleman_wunsch(
            seq_a, seq_b, match_score=1, mismatch_score=0, gap_penalty=-1
        )

        graph_file = "alignment_matrix.png"
        plot_alignment_matrix(score_matrix, path_coords, seq_a, seq_b, graph_file)

        save_to_pdf(seq_a, seq_b, aligned_a, aligned_b, match, gap, identity, 1, 0, -1, graph_file)
