from collections import defaultdict

def read_fasta(filepath):
    """Read multi-FASTA file, return dict of {identifier: sequence}."""
    records = {}
    with open(filepath, 'r') as file:
        seq_id = None
        seq_list = []
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_id:
                    records[seq_id] = ''.join(seq_list)
                seq_id = line[1:].split()[0]
                seq_list = []
            else:
                seq_list.append(line)
        if seq_id:
            records[seq_id] = ''.join(seq_list)
    return records

def count_records(records):
    """Count number of sequences in the FASTA records."""
    return len(records)

def length_stats(records):
    """Return lengths, longest and shortest sequence info."""
    lengths = {k: len(v) for k, v in records.items()}
    max_len = max(lengths.values())
    min_len = min(lengths.values())
    longest_seqs = [k for k, v in lengths.items() if v == max_len]
    shortest_seqs = [k for k, v in lengths.items() if v == min_len]
    return lengths, max_len, longest_seqs, min_len, shortest_seqs

def find_orfs(sequence, frame):
    """
    Find all ORFs in a given sequence for a specific forward reading frame (1, 2, or 3).
    Returns list of tuples (start_position_1_based, length).
    """
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    seq_len = len(sequence)
    frame_idx = frame - 1
    start_pos = None

    for i in range(frame_idx, seq_len - 2, 3):
        codon = sequence[i:i+3]
        if start_pos is None:
            if codon == start_codon:
                start_pos = i + 1  # 1-based numbering
        else:
            if codon in stop_codons:
                orf_len = i + 3 - start_pos + 1
                orfs.append((start_pos, orf_len))
                start_pos = None
    return orfs

def analyze_orfs(records, frame):
    """
    Analyze all sequences for ORFs in the given reading frame.
    Returns all ORFs per sequence, longest ORF length and sequence containing it.
    """
    all_orfs = {}
    longest_orf_length = 0
    longest_orf_id = None
    for seq_id, seq in records.items():
        orfs = find_orfs(seq, frame)
        all_orfs[seq_id] = orfs
        for (start_pos, length) in orfs:
            if length > longest_orf_length:
                longest_orf_length = length
                longest_orf_id = seq_id
    return all_orfs, longest_orf_length, longest_orf_id

def longest_orf_in_sequence(orfs_for_seq):
    """Return longest ORF (start_pos, length) in list of ORFs."""
    if not orfs_for_seq:
        return None
    return max(orfs_for_seq, key=lambda x: x[1])

def find_repeats(records, n):
    """
    Find all repeats of length n and count their occurrences across all sequences.
    Returns repeat counts, their positions, the most frequent repeat(s) and max count.
    """
    repeat_counts = defaultdict(int)
    positions = defaultdict(list)
    for seq_id, seq in records.items():
        seq_len = len(seq)
        for i in range(seq_len - n + 1):
            substr = seq[i:i+n]
            repeat_counts[substr] += 1
            positions[substr].append((seq_id, i))
    max_count = max(repeat_counts.values()) if repeat_counts else 0
    most_frequent_repeats = [r for r, count in repeat_counts.items() if count == max_count]
    return repeat_counts, positions, most_frequent_repeats, max_count

def main():
    import argparse

    parser = argparse.ArgumentParser(description="DNA multi-FASTA analyser for exam questions.")
    parser.add_argument('fasta_file', help="Path to the input multi-FASTA file.")
    args = parser.parse_args()

    # Read sequences
    records = read_fasta(args.fasta_file)

    # Question 1: Number of records
    num_records = count_records(records)
    print(f"Number of records in the file: {num_records}")

    # Question 2: Lengths and shortest/longest info
    lengths, max_len, longest_seqs, min_len, shortest_seqs = length_stats(records)
    print(f"Lengths of sequences:")
    for seq_id, length in lengths.items():
        print(f"  {seq_id}: {length} bp")

    print(f"Longest sequence length: {max_len} bp")
    print(f"Longest sequence(s): {', '.join(longest_seqs)}")

    print(f"Shortest sequence length: {min_len} bp")
    print(f"Shortest sequence(s): {', '.join(shortest_seqs)}")

    # Interactive part: ORF analysis query
    print("\nOpen Reading Frames (ORFs) analysis on forward strands.")
    while True:
        try:
            frame_input = input("\nEnter forward reading frame (1, 2, or 3) for ORF analysis or 'q' to quit: ").strip()
            if frame_input.lower() == 'q':
                break
            frame = int(frame_input)
            if frame not in [1, 2, 3]:
                print("Invalid frame. Please enter 1, 2, or 3.")
                continue

            all_orfs, longest_orf_len, longest_orf_seq_id = analyze_orfs(records, frame)
            print(f"\nLongest ORF length in frame {frame}: {longest_orf_len}")
            print(f"Sequence with longest ORF: {longest_orf_seq_id}")

            seq_id_input = input("Enter sequence identifier to find its longest ORF or press enter to skip: ").strip()
            if seq_id_input:
                if seq_id_input in all_orfs:
                    longest_orf = longest_orf_in_sequence(all_orfs[seq_id_input])
                    if longest_orf:
                        start_pos, length = longest_orf
                        print(f"Longest ORF in sequence '{seq_id_input}': Start position {start_pos}, Length {length}")
                    else:
                        print(f"No ORFs found in sequence '{seq_id_input}' at frame {frame}.")
                else:
                    print(f"Sequence identifier '{seq_id_input}' not found.")

        except ValueError:
            print("Invalid input, please enter 1, 2, 3 or q to quit.")

    # Interactive part: Repeat analysis
    while True:
        repeat_input = input("\nEnter repeat length (n) to find repeats or 'q' to quit: ").strip()
        if repeat_input.lower() == 'q':
            break
        try:
            n = int(repeat_input)
            if n <= 0:
                print("Repeat length must be a positive integer.")
                continue

            repeat_counts, positions, most_frequent_repeats, max_count = find_repeats(records, n)
            print(f"\nMost frequent repeat(s) of length {n} occurring {max_count} times:")
            for repeat_seq in most_frequent_repeats:
                print(f"  Repeat: {repeat_seq}")
            # Optionally show counts for all repeats:
            # for seq, count in repeat_counts.items():
            #     if count > 1:
            #         print(f"Repeat {seq}: {count} times")

        except ValueError:
            print("Please enter a valid integer for repeat length.")

if __name__ == "__main__":
    main()
