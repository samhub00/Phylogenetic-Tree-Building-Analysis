### Phylogenetic-Tree-Building-Analysis

Created by: Sam Hubler, Andrew Le, and Ryan Cartwright

This project aims to systematically evaluate, compare, and characterize the strengths and limitations of prominent phylogenetic tree-building algorithms across varying data scenarios. 

## 1. Objectives

Evaluate Accuracy: Measure how accurately different algorithms recover the "true" tree topology and branch lengths under various evolutionary scenarios.
Characterize Robustness: Test algorithm performance against data limitations, including high substitution rates, unequal evolutionary rates (long-branch attraction), and small sample sizes.

Assess Scalability: Compare the computational time and memory requirements of different methods as the number of taxa increases.
Provide Best Practices: Offer recommendations on which algorithm to use based on dataset characteristics (e.g., small, large, sparse, or dense). 

## 2. Algorithms to be Tested

Distance-based:  Unweighted Pair Group Method with Arithmetic Mean (UPGMA), Neighbor-Joining (NJ)

Character-based (Parsimony): Maximum Parsimony (MP).

Character-based (Likelihood): Maximum Likelihood (ML) - specifically IQ-TREE.

## 3. Data Generation and Datasets

To overcome the lack of a "gold standard" in real-world data, the project will rely heavily on simulation, supplemented by carefully selected empirical data.

## A. Simulated Data (The Core Framework)

Using tools like ROSE, generate sequences along known "true" trees. Parameters to vary:

Taxa Number: Small (10-50), Medium (100-500), Large (1000+).

Sequence Length: Short (500 bp), Long (5000+ bp).

Evolutionary Model: Jukes-Cantor (JC), Hasegawa (HKY)

Evolutionary Rate: Low divergence (conserved) vs. high divergence (saturated).

## B. Empirical Data
Use established, curated datasets (e.g., BAliBASE, comparative RNA websites).
Select datasets with known "species trees" to compare against gene trees.

## 4. How to run
Please consult the INSTALL file for more information on installing and running the program.

## Code Details

main.py: Main is our text UI that will assist the user in deciding which phylogenetic tree building
algortithm suits them best.

# tree_builder.py: 

Tree builder is our tree building python module. It uses two main libraries, Biopython and piqtree. Biopython handles the UPGMA, NJ, and ML tree building algorithms. piqtree is used to implement the ML tree building algorithm. 

Inside tree builder also resides some file formatting methods. You can use these to change sequence files into the desired format.

Tree builder also handles our tree visualization which is done with matplotlib and biopythons draw method. Biopython also supports an ASCII art tree building although the graphic data is more comprehensive. 

# evaluation.py: 

Evaluation.py is where the evaluation of the trees runs. In the evaluation module users are able to evaluate their trees compartively as well as use the pc performance wrapper function to measure run times and ram usage. A comprehensive analysis of the tree builidng algorithms will result in info on run time, sequence length, taxonomy count (if unknown), and peak ram usage. 

Another thing that users are able to do using the evaluation module is measure their built trees against known true trees. In here users can get the normalized Robinson-Fould distance, patristic correlation and branch correlation scores using the dendropy library. 

There are wrapper methods for these evaluations that can evaluate folders at a time.

# file_manager.py

File manager was created to be able to grab the taxonomy data from the balibase xml files. It can convert files to clustal format which is what we chose to use for all of our analysis and tree building.

# visualization.py

Visualization is a module used to generate the graphs seen in the figures file. It will create scatter plots and heatmaps based on time and ram usage as well as accuracy metrics. Most of the scatter plots and heatmaps have a wrapper that generates graphs for all four of our algorithms used as this was the best and fastest way to compare all four together. The heatmaps and graphs are made using seaborn, matplotlib, and handle data passed to them through pandas dataframes. 

# FINAL_CSVs

This folder contains the final csvs from our tree building and analysis. They can be used to generate visuals to explore the data that we collected

# bash

This folder contains the shell programs I used with the ROSE software to generate the sequences used for the analysis. There are four config files for high and low divergence of the JC and HKY mutation scales. The run_all script allowed me to collect the data faster and I included it for ease of use and data replication. 

# Questions
If questions arise surrounding features, please email me at samuel.hubler@sjsu.edu and I will do my best to answer them.