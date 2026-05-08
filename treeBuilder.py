from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import SeqIO
import Bio.Phylo.TreeConstruction
import io
from piqtree import build_tree
from cogent3 import load_aligned_seqs
import matplotlib.pyplot as plt

def file_convert(input_file, i_type, output_file, o_type):
    SeqIO.convert(input_file, i_type, output_file, o_type)



#SH
def get_dm(file, file_format="clustal", model="identity"):
    """
    Takes in exactly one multiple sequence alignment file and returns a distance matrix.

    file_format accepted: (phylip interleaved is default)
        {
            "phylip": Phylip Interlaced/Interleaved files,
            "phylip-sequential": Phylip sequential files,
            "phylip-relaxed": Phylip relaxed interpretation which allows for long names,
            "maf": Multiple Alignment Format(MAF) avaliable from UCSC Genome Browser
            "clustal": Clustal F or W
        }

    models accepted: (Default is identity)

    "models" = ['identity', 'blastn', 'trans', 'benner6', 'benner22', 'benner74', 
        'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 
        'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75', 'blosum80', 'blosum85', 
        'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 
        'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250', 'pam30', 'pam300', 
        'pam60', 'pam90', 'rao', 'risler', 'structure']
    """
    aln = AlignIO.read(open(file), file_format)
    calculator = DistanceCalculator(model)
    dm = calculator.get_distance(msa=aln)
    return dm

def upgma(dm):
    """
    Uses the UPGMA algorithm to build a tree from a distance matrix.
    """
    constructor = DistanceTreeConstructor()
    upgma_tree = constructor.upgma(dm)
    return upgma_tree


def nj(dm):
    """
    Uses the Neighbor Joining tree building algorithm to build a tree from a distance matrix.
    """
    constructor = DistanceTreeConstructor()
    nj_tree = constructor.nj(dm)
    return nj_tree


def mp(file, starter_tree, a_format='clustal'):
    """
    Max Parsimony Tree building algorithm. Requires alignment file and starting tree file. 
    Alignment file default format is PHYLIP, and starting tree default file is newick.
    """

    #starter = Phylo.read(starter_tree, st_format) #newick is tree format default, can change
    aln = AlignIO.read(open(file), a_format)
    scorer = Bio.Phylo.TreeConstruction.ParsimonyScorer()
    searcher = Bio.Phylo.TreeConstruction.NNITreeSearcher(scorer)
    constructor = Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor(searcher, starter_tree)
    pars_tree = constructor.build_tree(aln)

    return pars_tree

def ml(file, model = 'JC', rs=None):
    """
    Maximum Likelihood Tree building algorithm. Requires alignment file in phylip format. 
    will use piqtree library for tree building. 
    Model is option from IQ-Tree, our default is JC, can be changed to any model supported by IQ-Tree.
    Models found here: https://iqtree.github.io/doc/Substitution-Models

    Random Seed: For reproducible results, a random seed may be specified. This value None by Default.

        Caution: 0 is a specific random seed. None is equivalent to no random seed being specified.


    Troubleshooting: ensure that your alignment file has sequences of equal lengths, and that sequences are complete
    """
    
    try:
        aln = load_aligned_seqs(file, "clustal", moltype='dna')
        ml_tree = build_tree(aln, model, rand_seed=rs)
        return ml_tree
    except Exception as e:
        print("Unable to build ML tree, check file format...")


def visualize(trees):
    """
    Visualizes a list of trees using matplotlib. Each tree is displayed in a subplot with branch lengths labeled.
    Opens in separate window
    """

    num_trees = len(trees)
    cols = 2
    rows = (num_trees + 1) // cols

    # Use squeeze=False to ensure 'axes' is always a 2D array, 
    # even if there's only one row.
    fig, axes = plt.subplots(rows, cols, figsize=(12, 5 * rows), squeeze=False)

    for i, tree in enumerate(trees):
        row = i // cols
        col = i % cols
        ax = axes[row, col]
        
        # add titile
        method_name = ["UPGMA", "NJ", "ML", "MP"][i] if i < 4 else f"Tree {i+1}"
        ax.set_title(method_name)

        # do_show=False prevents the loop from breaking
        Bio.Phylo.draw(tree, axes=ax, do_show=False, 
                       branch_labels=lambda c: f"{c.branch_length:.2f}" if c.branch_length else "")

    # Hide unused subplots (if you have an odd number of trees)
    for j in range(i + 1, rows * cols):
        fig.delaxes(axes[j // cols, j % cols])

    plt.tight_layout()
    plt.show() # Call show ONCE at the very end

def build_trees(
        file,
        model="identity",
        ml_model='JC',
        ml_random_seed=None
        ):
    """
    file: must be in clustal format, and should be a multiple sequence alignment file of DNA sequences.

    model: is the model used for distance matrix calculation, default is identity, 
        but can be any model supported by Bio.Phylo.TreeConstruction.DistanceCalculator

    ml_model: is the model used for maximum likelihood tree building, default is JC,
        but can be any model supported by IQ-Tree

    ml_random_seed: is the random seed used for maximum likelihood tree building, default is None,
        but can be any integer value. Caution: 0 is a specific random seed. None is equivalent to
        no random seed being specified.

    Builds trees using UPGMA, NJ, ML, and MP algorithms from a given alignment file. 

    Returns a list of trees in the order of UPGMA, NJ, ML, MP.
    """
    get_dm(file, file_format='clustal', model=model)

    trees = []

    upgma_tree = upgma(get_dm(file, file_format='clustal', model=model))
    trees.append(upgma_tree)
    nj_tree = nj(get_dm(file, file_format='clustal', model=model))
    trees.append(nj_tree)
    ml_tree = Bio.Phylo.read(io.StringIO(str(ml(file, model=ml_model, rs=ml_random_seed))), 'newick')
    trees.append(ml_tree)
    mp_tree = mp(file, nj_tree, a_format='clustal')
    trees.append(mp_tree)

    return trees
"""
trees = []

upgma_tree = upgma(get_dm(r'alignments and sequences\cats.aln-clustalw', file_format='clustal'))

trees.append(upgma_tree)

nj_tree = nj(get_dm(r'alignments and sequences\cats.aln-clustalw', file_format='clustal'))

trees.append(nj_tree)

ml_tree = Bio.Phylo.read(io.StringIO(str(ml(r'alignments and sequences\cats.aln-clustalw'))), 'newick')

trees.append(ml_tree)

mp_tree = mp(r'alignments and sequences\cats.aln-clustalw', nj_tree, a_format='clustal')

trees.append(mp_tree)

visualize(trees)

"""

#print("UPGMA tree results:")
#Bio.Phylo.draw_ascii(upgma_tree)

#print("NJ tree results: ")
#Bio.Phylo.draw(nj_tree, branch_labels=lambda c: f"{c.branch_length:.2f}" if c.branch_length else "")

#print("ML tree results: ")
#Bio.Phylo.draw_ascii(ml_tree)
#Bio.Phylo.draw(ml_tree, branch_labels=lambda c: f"{c.branch_length:.2f}" if c.branch_length else "")


#print("MP tree results: ")
#Bio.Phylo.draw_ascii(mp_tree)


"""
We need to hook up indelible or rose for true trees
"""

visualize(build_trees(r'alignments and sequences\cats.aln-clustalw', model='blosum62', ml_model='JC', ml_random_seed=1))