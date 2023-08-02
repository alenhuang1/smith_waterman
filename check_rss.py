def smith_waterman(query_sequence, subject_sequence, similarity_matrix, gap_penalty):
    query_sequence = query_sequence.upper()
    subject_sequence = subject_sequence.upper()

    rows = len(query_sequence) + 1
    cols = len(subject_sequence) + 1
    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    traceback_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    max_score = 0
    max_i, max_j = 0, 0
    for i in range(1, rows):
        for j in range(1, cols):
            similarity_score = similarity_matrix[query_sequence[i - 1]][subject_sequence[j - 1]]
            diag_score = score_matrix[i - 1][j - 1] + similarity_score
            up_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(0, diag_score, up_score, left_score)

            if score_matrix[i][j] == 0:
                traceback_matrix[i][j] = 0
            elif score_matrix[i][j] == diag_score:
                traceback_matrix[i][j] = 1
            elif score_matrix[i][j] == up_score:
                traceback_matrix[i][j] = 2
            else:
                traceback_matrix[i][j] = 3

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_i, max_j = i, j

    indel_aligned_subject = ""
    snp_aligned_subject = ""
    i, j = max_i, max_j
    first_match_index = 0
    started_alignment = False  # Add this line before the while loop
    while traceback_matrix[i][j] != 0 and score_matrix[i][j] > 0:
        if traceback_matrix[i][j] == 1:  # Diagonal
            indel_aligned_subject = subject_sequence[j - 1] + indel_aligned_subject
            if query_sequence[i - 1] == subject_sequence[j - 1]:
                snp_aligned_subject = subject_sequence[j - 1] + snp_aligned_subject
                print(snp_aligned_subject)
                print(i,j)
                started_alignment = True  # Set this to True when the first match is encountered
                if(first_match_index == 0):
                    first_match_index = j
            else:
                if started_alignment:
                    snp_aligned_subject = '-' + snp_aligned_subject
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 2:  # Up
            indel_aligned_subject = '-' + indel_aligned_subject
            if started_alignment:
                snp_aligned_subject = '-' + snp_aligned_subject
            i -= 1
        else:  # Left
            if score_matrix[i][j] == score_matrix[i][j-1] + gap_penalty:
                indel_aligned_subject = '-' + indel_aligned_subject
                if started_alignment:
                    snp_aligned_subject = '-' + snp_aligned_subject
            else:
                indel_aligned_subject = subject_sequence[j - 1] + indel_aligned_subject
                if started_alignment:
                    snp_aligned_subject = 'subject_sequence[j - 1]' + snp_aligned_subject
            j -= 1
            first_match_index -= 1
        if traceback_matrix[i][j] == 0:
            started_alignment = False  # Set this to False when the alignment ends
    first_match_index = max_j - len(query_sequence) + 1

    while len(snp_aligned_subject) > len(query_sequence):
        snp_aligned_subject = snp_aligned_subject[1:]

    aligned_query_length = len(query_sequence)

    if len(indel_aligned_subject) != len(query_sequence):
        match_count = sum(1 for q, s in zip(query_sequence, snp_aligned_subject) if q == s)
        mismatch_count = sum(1 for q, s in zip(query_sequence, indel_aligned_subject) if s == '-')
    else:
        match_count = sum(1 for q, s in zip(query_sequence, indel_aligned_subject) if q == s)
        mismatch_count = sum(1 for q, s in zip(query_sequence, snp_aligned_subject) if q != s and q != '-')

    similarity_score_percent = (match_count - mismatch_count) / (aligned_query_length) * 100

    return similarity_score_percent, aligned_query_length, indel_aligned_subject, snp_aligned_subject, first_match_index



def display_alignment(query_sequence, subject_sequence, aligned_subject, first_match_index):
    i = 1
    while i < len(subject_sequence) + 1:
        print("Reference: " + str(i) + "      " + subject_sequence[i - 1:i + 79])
        if first_match_index is not None and i <= first_match_index <= i + 80:
            if first_match_index + len(query_sequence) > i + 80:
                print("Query    : " + str(i) + "      " + "-" * (first_match_index - i) + query_sequence[:(i + 80) - first_match_index])
                print()
                print("Reference: " + str(i + 80) + "      " + subject_sequence[i + 79:i + 159])
                print("Query    : " + str(i + 80) + "      " + query_sequence[(80 + i) - first_match_index:] + "-" * (80 - len(query_sequence[(80 + i) - first_match_index:])))
                print()
                i += 160
                continue
            else:
                print("Query    : " + str(i) + "      " + "-" * (first_match_index - i) + query_sequence + "-" * ((i + 80) - len(query_sequence) - first_match_index))
                i += 80
                continue
        else:
            if i <= len(subject_sequence) + 1 <= i + 80:
                print("Query    : " + str(i) + "     ", "-" * (len(subject_sequence) + 1 - i))
            else:
                print("Query    : " + str(i) + "     ", "-" * 80)
        print()
        i += 80


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
print()
print(f"Similarity Score: {similarity_score_percent}%")
if(similarity_score_percent != 100):
    print()
    print("SNP:")
    print(f"Ref  : {first_match_index, snp_aligned_subject}")
    print(f"Query: {first_match_index, query_sequence.upper()}")
    print()
    print("Indel:")
    print(f"Ref  : {first_match_index, indel_aligned_subject}")
    print(f"Query: {first_match_index, query_sequence.upper()}")
    print()