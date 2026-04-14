# Phylogenetic-Tree-Building-Analysis
This project aims to systematically evaluate, compare, and characterize the strengths and limitations of prominent phylogenetic tree-building algorithms across varying data scenarios. 

# 1. Objectives
Evaluate Accuracy: Measure how accurately different algorithms recover the "true" tree topology and branch lengths under various evolutionary scenarios.
Characterize Robustness: Test algorithm performance against data limitations, including high substitution rates, unequal evolutionary rates (long-branch attraction), and small sample sizes.

Assess Scalability: Compare the computational time and memory requirements of different methods as the number of taxa increases.
Provide Best Practices: Offer recommendations on which algorithm to use based on dataset characteristics (e.g., small, large, sparse, or dense). 

# 2. Algorithms to be Tested

Distance-based:  Unweighted Pair Group Method with Arithmetic Mean (UPGMA), Neighbor-Joining (NJ), Balanced Minimum Evolution (BME/FastME).

Character-based (Parsimony): Maximum Parsimony (MP).

Character-based (Likelihood): Maximum Likelihood (ML) - specifically RAxML and IQ-TREE.

# 3. Data Generation and Datasets

To overcome the lack of a "gold standard" in real-world data, the project will rely heavily on simulation, supplemented by carefully selected empirical data.

## A. Simulated Data (The Core Framework)
Using tools like INDELible or Rose, generate sequences along known "true" trees. Parameters to vary:

Taxa Number: Small (10-50), Medium (100-500), Large (1000+).

Sequence Length: Short (500 bp), Long (5000+ bp).

Evolutionary Model: Jukes-Cantor (JC), General Time Reversible (GTR) + Gamma (\(\Gamma \)) rate heterogeneity.

Evolutionary Rate: Low divergence (conserved) vs. high divergence (saturated).

## B. Empirical Data
Use established, curated datasets (e.g., BAliBASE, comparative RNA websites).
Select datasets with known "species trees" to compare against gene trees.
