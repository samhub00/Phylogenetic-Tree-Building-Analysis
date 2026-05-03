import time
from Bio import Phylo
from treeBuilder import upgma, nj

alignment_file = "birds.aln-clustalw"

def time_algorithm(name, func, file):
    start = time.time()
    tree = func(file)
    end = time.time()
    runtime = end - start
    print(f"{name} runtime: {runtime:.6f} seconds")
    return tree, runtime

def print_branch_lengths(tree, label):
    print(f"\nBranch lengths for {label}:")
    found = False
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            print(f"{clade.name}: {clade.branch_length}")
            found = True
    if not found:
        print("No branch lengths found.")

upgma_tree, upgma_time = time_algorithm("UPGMA", upgma, alignment_file)
nj_tree, nj_time = time_algorithm("NJ", nj, alignment_file)

print("\nUPGMA Tree:")
Phylo.draw_ascii(upgma_tree)

print("\nNJ Tree:")
Phylo.draw_ascii(nj_tree)

print_branch_lengths(upgma_tree, "UPGMA")
print_branch_lengths(nj_tree, "NJ")

Phylo.write(upgma_tree, "upgma_tree.nwk", "newick")
Phylo.write(nj_tree, "nj_tree.nwk", "newick")

print("\nSaved trees as upgma_tree.nwk and nj_tree.nwk")