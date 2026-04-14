from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Align
from Bio.Align import MultipleSeqAlignment

"""
birds.aln-clustalw is a multiple sequence alignment of 5 bird species:
I pulled the data off of NCBI, aligned with muscle in clustalw format,
and then used that to verify the tree construction
"""
aln = Align.read(open('birds.aln-clustalw'), 'clustal')
constructor = DistanceTreeConstructor()
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(msa=aln)

# UPGMA method constructor
upgma_tree = constructor.upgma(dm)
print(upgma_tree)

tree = constructor.nj(dm)
print(tree) 