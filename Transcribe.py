from Bio.Seq import Seq

query_dna = Seq(input("Enter DNA Sequence: "))

query_rna = Seq(query_dna).transcribe()

print(query_rna)


