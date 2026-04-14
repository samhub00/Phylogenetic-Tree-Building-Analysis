from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Align
from Bio.Align import MultipleSeqAlignment

"""
birds.aln-clustalw is a multiple sequence alignment of 5 bird species:
I pulled the data off of NCBI, aligned with muscle in clustalw format,
and then used that to verify the tree construction
"""

class DistanceBasedAlignment:

    def __init__(self, calculator='identity'):
        """
        Initializes the Distance Based Alignment class. Required are the type of constructor and the type of 
        calculator you would like to use. 
        """
        self.constructor = DistanceTreeConstructor()
        self. calculator = DistanceCalculator(calculator)

  
    def upgma(self, file):
        """
        Uses the UPGMA algorithm to build a tree from a multiple sequence alignment file located in the same directory with
        the clustal format. 
        """
        aln = Align.read(open(file), 'clustal')
        constructor = self.constructor
        calculator = self.calculator
        dm = calculator.get_distance(msa=aln)
        upgma_tree = constructor.upgma(dm)
        Phylo.draw_ascii(upgma_tree)
    
 
    def nj(self, file):
        """
        Uses the Neighbor Joining tree building algorithm to build a tree from a multiple sequence alignment file located 
        in the same directory with the clustal format.
        """
        aln = Align.read(open(file), 'clustal')
        constructor = self.constructor
        calculator = self.calculator
        dm = calculator.get_distance(msa=aln)
        nj_tree = constructor.nj(dm)
        Phylo.draw_ascii(nj_tree)

tree_builder = DistanceBasedAlignment()

tree_builder.upgma('birds.aln-clustalw')
tree_builder.nj('birds.aln-clustalw')