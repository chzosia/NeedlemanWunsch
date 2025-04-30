# DNA Sequence Alignment with Needleman-Wunsch

This project implements the **Needleman–Wunsch algorithm** for global alignment of DNA sequences. It can compare two sequences provided via command-line arguments or from a FASTA file, and outputs a visual matrix and a PDF report.

## Features

- Input via FASTA file or direct command-line sequences
- Validation of DNA sequences (A, C, G, T only)
- Global alignment using Needleman–Wunsch
- Matrix visualization with alignment path
- Auto-generated PDF report with summary and alignment plot

---

## How to Run

### Option 1: From FASTA file

```bash
python BioInfTask1.py --file your_sequences.fasta
```

### Option 2: From command-line

```bash
python BioInfTask1.py -a ACTGACTG -b ACTACTG
```
