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

#file_convert('BB20001.msf', 'msf', 'BB20001.clustal', 'clustal')

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
        aln = load_aligned_seqs(file, "clustal", moltype='DNA') #chagned from dna to protein for test
        ml_tree = build_tree(aln, model, rand_seed=rs)
        return ml_tree
    except Exception as e:
        print("Unable to build ML tree, check file format...")

#print(ml(r'ROSE_trees_sequences\HKY_high\clustal\HKY_HIGH_100_500_1.clustal'))


def visualize(trees, tree_names=None):
    """
    Visualizes a list of trees using matplotlib. Each tree is displayed in a subplot
    with branch lengths labeled. Opens in a separate window.

    Parameters:
        trees: list
            List of tree objects to visualize.

        tree_names: list or None
            Optional list of names for the trees. Example: ["UPGMA", "NJ"].
            If not provided, generic names will be used.
    """

    num_trees = len(trees)
    cols = 2
    rows = (num_trees + 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(12, 5 * rows), squeeze=False)

    for i, tree in enumerate(trees):
        row = i // cols
        col = i % cols
        ax = axes[row, col]

        if tree_names and i < len(tree_names):
            method_name = tree_names[i]
        else:
            method_name = f"Tree {i + 1}"

        ax.set_title(method_name)

        Bio.Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            branch_labels=lambda c: f"{c.branch_length:.2f}" if c.branch_length else ""
        )

    for j in range(i + 1, rows * cols):
        fig.delaxes(axes[j // cols, j % cols])

    plt.tight_layout()
    plt.show()

def build_trees(
        file,
        model="identity",
        ml_model='JC',
        ml_random_seed=None,
        desired_trees=["UPGMA", "NJ", "ML", "MP"]
        ):
    """
    Parameters

        file: str
        Path to a Multiple Sequence Alignment (MSA) file. Must be in Clustal format and contain DNA sequences.

        model: str (default: "identity")
        The model used for distance matrix calculation. Supports any model compatible with Bio.Phylo.TreeConstruction.DistanceCalculator.

        ml_model: str (default: "JC")
        The model used for Maximum Likelihood tree building. Supports any model supported by IQ-TREE.

        ml_random_seed: int (default: None)
        The random seed for ML tree building.

            0 is treated as a specific seed.

            None is equivalent to no random seed being specified.

        desired_trees: list (default: ["UPGMA", "NJ", "ML", "MP"])
        A list of tree-building methods to execute. Can be any combination of the four supported algorithms.

    Returns list of trees in the order of UPGMA, NJ, ML, MP (if specified in desired_trees)
    """


    get_dm(file, file_format='clustal', model=model)

    trees = []

    #build trees using each method and add to list
    if "UPGMA" in desired_trees:
        upgma_tree = upgma(get_dm(file, file_format='clustal', model=model))
        trees.append(upgma_tree)
    if "NJ" in desired_trees:
        nj_tree = nj(get_dm(file, file_format='clustal', model=model))
        trees.append(nj_tree)
    if "ML" in desired_trees:
        ml_tree = Bio.Phylo.read(io.StringIO(str(ml(file, model=ml_model, rs=ml_random_seed))), 'newick')
        trees.append(ml_tree)
    if "MP" in desired_trees:
        mp_tree = mp(file, nj_tree, a_format='clustal')
        trees.append(mp_tree)
    if not trees:
        print("build_trees() error. No valid tree building method specified, please choose from UPGMA, NJ, ML, MP")
    #return list of trees in the order of UPGMA, NJ, ML, MP
    return trees

"""
TODO:
    hook up indelible or rose for true trees
    add ncbi number to name for display purposes
    export tree lengths for evaluation
"""

#visualize(build_trees(r'bb3_release\RV50\clustal_files\BBS50006.clustal', model='blosum62', ml_model='Blosum62', ml_random_seed=1))