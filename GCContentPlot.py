
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

import matplotlib.pyplot as plt #For importing module of plotting graphs

sequences = SeqIO.parse("FASTASeq.txt", "fasta")
gc_values = []

for seq_record in sequences: 
    gc_values.append(gc_fraction(seq_record.seq) * 100)
gc_values = sorted(gc_values)
#Since this version doesnt support GC module, we use gc_fraction to calculate GC content

plt.plot(gc_values)
plt.title(f"{len(gc_values)} Sequences\nGC% Range: {min(gc_values):.2f} to {max(gc_values):.2f}")
plt.ylabel("GC%")
plt.xlabel("Genes")
plt.show()
