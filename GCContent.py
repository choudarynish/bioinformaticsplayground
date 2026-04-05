#To compare the highest GC contents amongst various FASTA sequences
#Add FASTA Sequences in FASTASeq.txt

from Bio import SeqIO

seq_record_max = None
gc_max = 0
for seq_record in SeqIO.parse("FASTASeq.txt" , "fasta"):
    sequence = str(seq_record)
    sequence_id = seq_record.id
    gc = (sequence.count("G") + sequence.count("C"))/len(sequence)*100
    if gc > gc_max:
        seq_record_max = sequence_id
        gc_max = gc

print("Highest GC content in: " ,seq_record.id ,"\nGC Content in %: ", gc_max)

#IN Built function for GC Content 
#from Bio.SeqUtils import GC
#newseq = Seq("GAATTCGAATTC")
#GC(newseq)
