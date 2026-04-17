from Bio import AlignIO, Phylo
import Bio.Phylo.TreeConstruction

#SH
def parsimonyTreeBuilder(alignment,starter_tree, a_format='phylip',  st_format='newick'):
    """
    Max Parsimony Tree building algorithm. Requires alignment file and starting tree file. 
    Alignment file default format is PHYLIP, and starting tree default file is newick.
    """
    aln = AlignIO.read(open(alignment), a_format) 

    starter = Phylo.read(starter_tree, st_format) #newick is tree format default, can change

    scorer = Bio.Phylo.TreeConstruction.ParsimonyScorer()
    searcher = Bio.Phylo.TreeConstruction.NNITreeSearcher(scorer)
    constructor = Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starter)
    pars_tree = constructor.build_tree(aln)

    return pars_tree

