from Bio.Seq import Seq

query = Seq(input("Enter DNA Sequence: "))

query_rc = Seq(query).reverse_complement()

print(query_rc)
