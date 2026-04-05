#To calculate the Hamming distance between two strings

from Bio import SeqIO

def hamming_distance(seq1 , seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences should be of same length")
    return sum(nt1 != nt2 for nt1, nt2 in zip(seq1, seq2))

queries = list(SeqIO.parse("rosalind_hamm.txt", "fasta")) #Seperate query seq into 2 fasta seqs
query1 = queries[0]
query2 = queries[1]

print(f"Hamming distance is: {hamming_distance(query1.seq, query2.seq)}")