from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Align
from Bio.Align import MultipleSeqAlignment
import Bio.Phylo.TreeConstruction
import subprocess

#SH
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

    def bme(self, file, output_file = 'bme_tree.nwk'):
        """
        Builds a tree using FastME for balanced minimum evolution.
        """
        pass

    def parsimonyTreeBuilder(alignment,starter_tree, a_format='phylip',  st_format='newick'):
        """
        Max Parsimony Tree building algorithm. Requires alignment file and starting tree file. 
        Alignment file default format is PHYLIP, and starting tree default file is newick.
        """
        aln = Align.read(open(alignment), a_format) 

        starter = Phylo.read(starter_tree, st_format) #newick is tree format default, can change

        scorer = Bio.Phylo.TreeConstruction.ParsimonyScorer()
        searcher = Bio.Phylo.TreeConstruction.NNITreeSearcher(scorer)
        constructor = Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starter)
        pars_tree = constructor.build_tree(aln)

        return pars_tree


tree_builder = DistanceBasedAlignment()

"""
birds.aln-clustalw is a multiple sequence alignment of 5 bird species:
I pulled the data off of NCBI, aligned with muscle in clustalw format,
and then used that to verify the tree construction
"""

tree_builder.upgma('birds.aln-clustalw')
tree_builder.nj('birds.aln-clustalw')