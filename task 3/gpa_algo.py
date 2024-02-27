import sys

def read_sequence_from_file(filename):
    with open(filename, 'r') as file:
        sequence = file.read().strip()
    return sequence

def global_alignment(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    m, n = len(seq1), len(seq2)
    # Initialize the DP matrix
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize the traceback matrix to keep track of alignment path
    traceback = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize the first row and column with gap penalties
    for i in range(1, m + 1):
        dp[i][0] = dp[i-1][0] + gap_penalty
        traceback[i][0] = 1
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j-1] + gap_penalty
        traceback[0][j] = 2
    
    # Fill the DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                match = dp[i-1][j-1] + match_score
            else:
                match = dp[i-1][j-1] + mismatch_penalty
            gap1 = dp[i-1][j] + gap_penalty
            gap2 = dp[i][j-1] + gap_penalty
            dp[i][j] = max(match, gap1, gap2)
            
            # Update the traceback matrix
            if dp[i][j] == match:
                traceback[i][j] = 3  # diagonal
            elif dp[i][j] == gap1:
                traceback[i][j] = 1  # up
            else:
                traceback[i][j] = 2  # left
                
    # Traceback to construct the aligned sequences
    aligned_seq1 = ''
    aligned_seq2 = ''
    i, j = m, n
    while i > 0 or j > 0:
        if traceback[i][j] == 3:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    return aligned_seq1, aligned_seq2, dp[m][n]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python alignment.py <input_file1> <input_file2>")
        sys.exit(1)

    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]

    sequence1 = read_sequence_from_file(input_file1)
    sequence2 = read_sequence_from_file(input_file2)

    alignment1, alignment2, score = global_alignment(sequence1, sequence2)

    print("Alignment Score:", score)
    print("Aligned Sequence 1:", alignment1)
    print("Aligned Sequence 2:", alignment2)
