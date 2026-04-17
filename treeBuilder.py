from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Align
from Bio.Align import MultipleSeqAlignment
import Bio.Phylo.TreeConstruction
import subprocess
from piqtree import build_tree
from cogent3 import load_aligned_seqs


#SH


def upgma(file):
    """
    Uses the UPGMA algorithm to build a tree from a multiple sequence alignment file located in the same directory with
    the phylip format. 
    """

    aln = Align.read(open(file), 'phylip')
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator(calculator)
    dm = calculator.get_distance(msa=aln)
    upgma_tree = constructor.upgma(dm)
    return upgma_tree


def nj(file):
    """
    Uses the Neighbor Joining tree building algorithm to build a tree from a multiple sequence alignment file located 
    in the same directory with the phylip format.
    """

    aln = Align.read(open(file), 'phylip')
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator(calculator)

    dm = calculator.get_distance(msa=aln)
    nj_tree = constructor.nj(dm)
    return nj_tree


def bme(file, output_file = 'bme_tree.nwk'):
    """
    Builds a tree using FastME for balanced minimum evolution. Takes in a multiple sequence alignment file
    in phylip format and an optional output file name (default is bme_tree.nwk).
    Will use subprocess call to fastME in cmd line format.
    
    """
    
    print("Building tree using FastME for balanced minimum evolution...")

    #FastME Commands:
    # -i input file
    # -o output file
    # -m tree building method (default B for BME)
    # -B ouput bootstrap trees file (optional)
    # refer to https://gite.lirmm.fr/atgc/FastME/ README for more details on FastME command line options

    cmd = ["fastme", "-i", file, "-o", output_file, "-m", "B"]

    try:
        subprocess.run(cmd, check=True)
        bme_tree = Phylo.read(output_file, 'newick')
        return bme_tree
    except subprocess.CalledProcessError as e:
        print("Error occurred while running FastME: ", e)
        return None


def mp(file, starter_tree, a_format='phylip',  st_format='newick'):
    """
    Max Parsimony Tree building algorithm. Requires alignment file and starting tree file. 
    Alignment file default format is PHYLIP, and starting tree default file is newick.
    """
    aln = Align.read(open(file), a_format) 

    starter = Phylo.read(starter_tree, st_format) #newick is tree format default, can change

    scorer = Bio.Phylo.TreeConstruction.ParsimonyScorer()
    searcher = Bio.Phylo.TreeConstruction.NNITreeSearcher(scorer)
    constructor = Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starter)
    pars_tree = constructor.build_tree(aln)

    return pars_tree

def ml(file, model = 'JC', use_random_seed = False):
    """
    Maximum Likelihood Tree building algorithm. Requires alignment file in phylip format. 
    will use piqtree library for tree building. 
    Model is option from IQ-Tree, our default is JC, can be changed to any model supported by IQ-Tree.
    Models found here: https://iqtree.github.io/doc/Substitution-Models
    """
    seed = 0

    if use_random_seed:
        seed = 1

    aln = load_aligned_seqs(file)
    ml_tree = build_tree(aln, model, rand_seed=seed)
    return ml_tree