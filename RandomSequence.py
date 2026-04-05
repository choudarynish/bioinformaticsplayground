#To genearate a random sequence of DNA/RNA/Polypeptide with given prompt of sequence type and length

from Bio.Seq import Seq
import random

def random_sequence_generate(seq_type , seq_len):
    if seq_type == 'DNA':
        bases = "ATGC"
    elif seq_type == 'RNA':
        bases = "AUGC"
    elif seq_type == 'AA':
        bases = "ACDEFGHIKLMNPQRSTVWY"
    else:
        raise ValueError("Enter the required properly\n")
    
    final_seq = ''.join(random.choices(bases, k=seq_len))
    return Seq(final_seq)

qtype = input("Sequence?(DNA/RNA/AA): ")
qlen = int(input("SequenceLength?: "))
qfinal_seq = random_sequence_generate(qtype, qlen)

print(f"Random {qtype} sequence ({qlen}): {qfinal_seq}") 