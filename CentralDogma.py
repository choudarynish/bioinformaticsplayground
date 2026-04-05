#To understand the central dogma of genetic information
#Try inputting a FASTA sequence and giving output

from Bio.Seq import Seq

qseq = Seq(input("Enter Sequence: "))

qseq_dnacoding = qseq
qseq_dnanoncoding = Seq(qseq_dnacoding).complement()
qseq_rna = Seq(qseq_dnacoding).transcribe()
qseq_aa = Seq(qseq_dnacoding).translate()

print(f"\nDNA Strands\n 3'-{qseq_dnacoding}-5' \n 5'-{qseq_dnanoncoding}-3'")
print(f"\nTranscripted RNA Strand\n 3'-{qseq_rna}-5'")
print(f"\nTranslated Polypeptide Sequence\n 3'-{qseq_aa}-5'\n")

# * in aa seq means stop codon
#GAATTACAATAACGCGCGCGGGCGGATTTTATTATATATATATCGCGCGCGGC