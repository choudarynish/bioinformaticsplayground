from Bio import SeqIO
import matplotlib.pyplot as plt

sequences = SeqIO.parse("FASTASeq.txt", "fasta")

sizes = [len(rec) for rec in sequences]


plt.hist(sizes , bins=20) #bins=20 for 20 divisions
plt.title(f"{len(sizes)} Sequences\n Length Range: {min(sizes)} to {max(sizes)}")
plt.ylabel("Number of sequences")
plt.xlabel("Sequence length in bp")
plt.show() 
