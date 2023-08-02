# Smith-Waterman Sequence Alignment Tool README

The Smith-Waterman Sequence Alignment Tool is a Python-based utility that implements the Smith-Waterman local sequence alignment algorithm. This algorithm is widely used in bioinformatics for identifying similar regions between two nucleotide or protein sequences. In particular, it is used for comparing relatively short sequences against a longer sequence or genome.

## How to Use

To use the tool, you need to provide a query sequence, a subject sequence, a similarity matrix (also known as a substitution matrix), and a gap penalty.

```python
# Example query and subject sequences
query_sequence = "cacagtgactatgaggcctccttaactgtgccaaaattc"
subject_sequence = "tattgctgctgctgttactgcaaccatttatttcagtctaagaaattctcccatcaatggcagttcttttgtgaccacatggaagcatcatttaaaaaattattccaatagtttttggaggaaacatcatttttaataatgatggggcttctgggggtgctgccctagtaacaatcatgtatcttgtcataggcactgcaggacaaagccttgagcagccctctgaagtgacagctgtggaaggagccattgtccagataaactgcacgtaccagacatctgggttttatgggctgtcctggtaccagcaacatgatggcggagcacccacatttctttcttacaatgctctggatggtttggaggagacaggtcgtttttcttcattccttagtcgctctgatagttatggttacctccttctacaggagctccagatgaaagactctgcctcttacttctgcgctgtgagagacgcacagtgactatgaggcctccttaactgtgccaaaattcaaaagacaatcagtggagtacaggtgggcttgagaagttctagaacttcctgagtgtatctttgcttaccgtctaattttaaacttgtttagaatatgtagcatagggctctgagacggaaccccattattttctgaaccattgattgctaaagagaagga"
# Gap penalty
gap_penalty = -2

# Similarity matrix (substitution matrix) for DNA sequences
similarity_matrix = {
    'A': {'A': 5, 'C': -1, 'G': -2, 'T': -1, '-': gap_penalty},
    'C': {'A': -1, 'C': 5, 'G': -1, 'T': -2, '-': gap_penalty},
    'G': {'A': -2, 'C': -1, 'G': 5, 'T': -1, '-': gap_penalty},
    'T': {'A': -1, 'C': -2, 'G': -1, 'T': 5, '-': gap_penalty},
    '-': {'A': gap_penalty, 'C': gap_penalty, 'G': gap_penalty, 'T': gap_penalty, '-': gap_penalty}
}

# Call the Smith-Waterman algorithm
similarity_score_percent, aligned_query_length, indel_aligned_subject, snp_aligned_subject, first_match_index = smith_waterman(query_sequence, subject_sequence, similarity_matrix, gap_penalty)

# Display the aligned sequences in the desired format
display_alignment(query_sequence, subject_sequence, indel_aligned_subject, first_match_index)
```

## Function Definitions

### smith_waterman()

The main function is `smith_waterman(query_sequence, subject_sequence, similarity_matrix, gap_penalty)`. This function performs the Smith-Waterman local sequence alignment algorithm and returns:

- `similarity_score_percent`: The similarity score in percent. It is calculated as the number of matched characters minus the number of mismatched characters, divided by the length of the aligned query, all multiplied by 100.
- `aligned_query_length`: The length of the aligned query sequence.
- `indel_aligned_subject`: The aligned subject sequence considering both insertions/deletions (indels) and substitutions.
- `snp_aligned_subject`: The aligned subject sequence considering only substitutions.
- `first_match_index`: The index of the first match in the aligned sequences.

### display_alignment()

The `display_alignment(query_sequence, subject_sequence, aligned_subject, first_match_index)` function displays the query and subject sequence alignment. It takes the same inputs as `smith_waterman()`. It prints the alignment in a legible, structured format.

## Output

The script outputs the similarity score as a percentage and, if the similarity score is not 100%, it outputs the reference (subject) and query sequences for both substitutions (SNPs) and indels.
