import time
import psutil
import os
import csv
import dendropy
from dendropy.calculate import treecompare
from treeBuilder import build_trees
from Bio import AlignIO

msa_file = r'alignments and sequences\cats.aln-clustalw'

#trees = build_trees(msa_file, model='blosum62', ml_model='JC', ml_random_seed=1)

#FILE EVALUATION METRICS

def get_msa_length(file, file_format='clustal'):
    """
    Get the length of the multiple sequence alignment (MSA) from a given file.
    """

    alignment = AlignIO.read(file, file_format)
    return alignment.get_alignment_length()


#TREE EVALUATION METRICS MEASURE AND SAVE TO CSV

def get_RF(tree1, tree2):
    """
    Calculate the Robinson-Foulds distance between two trees.
    """
    tns = dendropy.TaxonNamespace()
    tree1 = dendropy.Tree.get(data=tree1.format('newick'), schema='newick', taxon_namespace=tns)
    tree2 = dendropy.Tree.get(data=tree2.format('newick'), schema='newick', taxon_namespace=tns)

    tree1.encode_bipartitions()
    tree2.encode_bipartitions()

    return treecompare.symmetric_difference(tree1, tree2)

#print(get_RF(trees[0], trees[3]))



#PC MONITORING METRICS MEASURE AND SAVE TO CSV

def monitor_performance(function, *args, **kwargs):
    """
    Monitor the performance of a function by measuring its execution time cpu time and memory usage.
    """
    monitoring_results = {}

    process=psutil.Process(os.getpid())
    start_cpu_times = process.cpu_times().user + process.cpu_times().system
    start_time = time.perf_counter()

    result = function(*args, **kwargs)
    
    end_time = time.perf_counter()
    end_cpu_times = process.cpu_times().user + process.cpu_times().system

    mem_info = process.memory_info()
    peak_ram_mb = mem_info.peak_wset / (1024 * 1024)

    monitoring_results['execution_time_seconds'] = end_time - start_time
    monitoring_results['cpu_time_seconds'] = end_cpu_times - start_cpu_times
    monitoring_results['peak_ram_mb'] = peak_ram_mb

    return result, monitoring_results

def save_performance_metrics(msa_file, metrics, output_path='msa_performance_metrics.csv'):
    """
    Saves the computational performance metrics to a csv file. If the file has already been analyzed it will
    append to the end. Requires the msa_file and takes the monitoring results from the monitor_performance 
    function as a dictionary.
    """
    file_exists = os.path.isfile(output_path)
    with open(output_path, mode='a', newline='') as csvfile:
        fieldnames = ['msa_file','seq_len', 'execution_time_seconds', 'cpu_time_seconds', 'peak_ram_mb']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if not file_exists:
            writer.writeheader()

        row = {'msa_file': msa_file}
        row['seq_len'] = get_msa_length(msa_file)
        row.update(metrics)
        writer.writerow(row)



trees, pc_perf = monitor_performance(
    build_trees, 
    msa_file, 
    model='blosum62', 
    ml_model='JC', 
    ml_random_seed=1
    )

save_performance_metrics(msa_file, pc_perf)

"""
print("Performance Metrics:")
print(f"Execution Time: {pc_perf['execution_time_seconds']:.2f} seconds")
print(f"CPU Time: {pc_perf['cpu_time_seconds']:.2f} seconds")
print(f"Peak RAM Usage: {pc_perf['peak_ram_mb']:.2f} MB")
"""



"""
Evaluation metrics that we need:
    True trees from INDELible or Rose
    Parameters: Evolutionary models like Jukes-Cantor, general Time reversible
        These are just parameters for the tree building algorithm
    Topological Accuracy: nRF distance which measures the distance between the inferred tree and true trees (where do we get the true tree?)
        implemented - still need to figure out how to get true trees
    Branch Length Accuracy: Correlation coefficient between inferred branch lengths and true simulated branch lengths (which ones are simulate?)
    Support values: Evaluate the bootstrap percentages for the ML/NJ or posterior probabilities for BI to assess Node Reliability
    CPU usage? is there a way to track the computing done 
        done - using psutil to track cpu time and memory usage
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