from Bio.PDB import PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from Bio import SeqIO
 
handle = open("1l2y.pdb", "rU")
file1 = open("file.txt","a")

for record in SeqIO.parse(handle, "pdb-seqres") :
	x = ">" + record.id + "\n" + record.seq
	seq1 = str(record.seq)
	file1.write(seq1)
	print x
handle.close()


f = open('file.txt', 'r')
my_seq = f.read()


analysed_seq = ProteinAnalysis(my_seq)
print analysed_seq.molecular_weight()
print analysed_seq.gravy()
print analysed_seq.get_amino_acids_percent()
print analysed_seq.secondary_structure_fraction()
print analysed_seq.isoelectric_point()
print analysed_seq.instability_index()
print analysed_seq.aromaticity()
