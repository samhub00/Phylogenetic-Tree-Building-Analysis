import time
import csv
from Bio import Phylo
from treeBuilder import upgma, nj


# Input dataset
alignment_file = "wolf_and_pals.aln-phylip_interleaved"
alignment_format = "phylip"


def time_algorithm(name, func, file, file_format):
    """
    Runs a tree-building algorithm and records how long it takes.
    """
    start = time.time()
    tree = func(file, file_format)
    end = time.time()

    runtime = end - start
    print(f"{name} runtime: {runtime:.6f} seconds")

    return tree, runtime


def print_branch_lengths(tree, label):
    """
    Prints every branch length in a tree.
    """
    print(f"\nBranch lengths for {label}:")
    found = False

    for clade in tree.find_clades():
        if clade.branch_length is not None:
            print(f"{clade.name}: {clade.branch_length}")
            found = True

    if not found:
        print("No branch lengths found.")


def count_taxa(tree):
    """
    Counts how many terminal taxa/species are in the tree.
    """
    return len(tree.get_terminals())


def total_branch_length(tree):
    """
    Adds up all branch lengths in the tree.
    """
    total = 0

    for clade in tree.find_clades():
        if clade.branch_length is not None:
            total += clade.branch_length

    return total


def save_results_csv(results, output_file="evaluation_results.csv"):
    """
    Saves evaluation results to a CSV file.
    """
    fieldnames = [
        "algorithm",
        "runtime_sec",
        "num_taxa",
        "total_branch_length",
        "tree_file"
    ]

    with open(output_file, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved evaluation results to {output_file}")


# Build trees and measure runtime
upgma_tree, upgma_time = time_algorithm(
    "UPGMA",
    upgma,
    alignment_file,
    alignment_format
)

nj_tree, nj_time = time_algorithm(
    "NJ",
    nj,
    alignment_file,
    alignment_format
)


# Display trees in the terminal
print("\nUPGMA Tree:")
Phylo.draw_ascii(upgma_tree)

print("\nNJ Tree:")
Phylo.draw_ascii(nj_tree)


# Print branch lengths
print_branch_lengths(upgma_tree, "UPGMA")
print_branch_lengths(nj_tree, "NJ")


# Save trees as Newick files
Phylo.write(upgma_tree, "wolf_upgma_tree.nwk", "newick")
Phylo.write(nj_tree, "wolf_nj_tree.nwk", "newick")

print("\nSaved trees as wolf_upgma_tree.nwk and wolf_nj_tree.nwk")


# Save evaluation metrics to CSV
results = [
    {
        "algorithm": "UPGMA",
        "runtime_sec": upgma_time,
        "num_taxa": count_taxa(upgma_tree),
        "total_branch_length": total_branch_length(upgma_tree),
        "tree_file": "wolf_upgma_tree.nwk"
    },
    {
        "algorithm": "NJ",
        "runtime_sec": nj_time,
        "num_taxa": count_taxa(nj_tree),
        "total_branch_length": total_branch_length(nj_tree),
        "tree_file": "wolf_nj_tree.nwk"
    }
]

save_results_csv(results)

"""
Evaluation metrics that we need:
    True trees from INDELible or Rose
    Parameters: Evolutionary models like Jukes-Cantor, general Time reversible
    Topological Accuracy: nRF distance which measures the distance between the inferred tree and true trees (where do we get the true tree?)
    Branch Length Accuracy: Correlation coefficient between inferred branch lengths and true simulated branch lengths (which ones are simulate?)
    Support values: Evaluate the bootstrap percentages for the ML/NJ or posterior probabilities for BI to assess Node Reliability
    CPU usage? is there a way to track the computing done 
    """

"""
Comparative Analysis
    R packages like ape, phangorn, and tree cluster? any python package avaliability?
    Calculate nRF, bootstrap support and running times for every result
    Simulate specific datasets with two long non-sister branches to test for systematic errors
    (what is a systematic error?)
    
    """

"""
Determine which method are the best at which tasks

Visualization:
    Scalabiilty graph: time taken / number of taxa
    Accuracy heatmap: Seaborn or something? topological representation of accuracy accross different evo rates and tree sizes
    decision tree or flowchart that helps researchers select the best algorithm 
"""