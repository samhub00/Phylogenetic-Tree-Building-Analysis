from pathlib import Path
import time
import psutil
import threading
import os
import csv
import dendropy
import gc
from dendropy.calculate import treecompare
from treeBuilder import build_trees
from Bio import AlignIO
from Bio import Phylo
import pandas as pd


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
    Monitors the performance of a given function by tracking its execution time, CPU time, and peak RAM usage.
    """

    monitoring_results = {}
    process = psutil.Process(os.getpid())
    
    baseline_mem = process.memory_info().rss / (1024 * 1024)  # in MB

    # Shared variable to store the peak found by the thread
    peak_mem_during_op = [baseline_mem]
    stop_event = threading.Event()

    def track_memory():
        while not stop_event.is_set():
            # Current memory in MB
            current_mem = process.memory_info().rss / (1024 * 1024)
            if current_mem > peak_mem_during_op[0]:
                peak_mem_during_op[0] = current_mem
            # Small sleep to prevent high CPU usage from the monitor itself
            time.sleep(0.005) 

    # Start the memory monitoring thread
    mem_thread = threading.Thread(target=track_memory)
    mem_thread.start()

    # Performance monitoring start
    start_cpu_times = process.cpu_times().user + process.cpu_times().system
    start_time = time.perf_counter()

    try:
        result = function(*args, **kwargs)
    finally:
        # Performance monitoring end
        end_time = time.perf_counter()
        end_cpu_times = process.cpu_times().user + process.cpu_times().system
        
        # Stop the thread
        stop_event.set()
        mem_thread.join()

    monitoring_results['execution_time_seconds'] = end_time - start_time
    monitoring_results['cpu_time_seconds'] = end_cpu_times - start_cpu_times
    monitoring_results['peak_ram_mb'] = peak_mem_during_op[0] - baseline_mem

    return result, monitoring_results

def save_performance_metrics(msa_file, metrics, tree, output_path='computational_performance.csv'):
    """
    Saves the computational performance metrics to a csv file. If the file has already been analyzed it will
    append to the end. Requires the msa_file and takes the monitoring results from the monitor_performance 
    function as a dictionary.
    """
    file_exists = os.path.isfile(output_path)
    with open(output_path, mode='a', newline='') as csvfile:
        fieldnames = ['msa_file','tree', 'seq_len', 'execution_time_seconds', 'cpu_time_seconds', 'peak_ram_mb']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if not file_exists:
            writer.writeheader()

        msa_file_trimmed = os.path.basename(msa_file)[:-len('.clustal')]

        row = {'msa_file': msa_file_trimmed}
        row['seq_len'] = get_msa_length(msa_file)
        row['tree'] = tree
        row.update(metrics)
        writer.writerow(row)

#this is the overarching tree analysis fucnction
def analyze_file(msa_file):
    output_dir = Path("trees")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # List of algorithms to run
    tree_types = ["UPGMA", "NJ", "ML"] #["UPGMA", "NJ", "ML", "MP"]
    results = []

    for t_type in tree_types:
        # 1. Run and Monitor
        # Note: adjust desired_trees for MP if it requires 'NJ' as a starting tree
        dt = ['NJ', 'MP'] if t_type == "MP" else [t_type]
        
        tree_obj, pc_perf = monitor_performance(
            build_trees, 
            msa_file, 
            desired_trees=dt,
            model='blosum62', 
            ml_model='Blosum62', 
            ml_random_seed=1
        )
        
        # 2. Save Performance Metrics
        save_performance_metrics(msa_file, pc_perf, tree=t_type)

        # 3. Save the Tree File Immediately
        base_name = Path(msa_file).stem
        filename = output_dir / f"{base_name}_{t_type}.nwk"
        
        # Handle dict nesting if build_trees returns a dict
        actual_tree = tree_obj[t_type] if isinstance(tree_obj, dict) else tree_obj
        Phylo.write(actual_tree, filename, "newick")
        
        # 4. CRITICAL: Clear memory
        del tree_obj
        del actual_tree
        gc.collect() # Force the collector to release the RAM now
        
        print(f"Finished {t_type} for {msa_file}")

    print(f"Analysis complete for {msa_file}.")


#Sample analysis of a single file
#analyze_file(r'bb3_release\RV50\clustal_files\BBS50006.clustal')

def analyze_folder(folder_path):
    """
    Analyze all the msa files in a folder and save the performance metrics to a csv file. 
    """
    path = Path(folder_path)
    
    for msa_file in path.glob('*.clustal'):
        
        print(f"Analyzing {msa_file.name}...")
        analyze_file(msa_file)

#analyze_folder(r'msa_for_analysis')

def taxa_count(msa_file):
    with open(msa_file, 'r') as f:
        # Skip the header line (CLUSTAL X...)
        next(f) 
        
        in_block = False
        count = 0
        
        for line in f:
            # strip() removes whitespace and newlines
            clean_line = line.strip()
            
            if clean_line:
                # We found the start of the first sequence block
                in_block = True
                count += 1
            elif in_block:
                # We hit the first empty line AFTER the first block
                # This is our answer.
                return count
                
        return count
#print(taxa_count(r'processed_msa\BB11003.clustal'))

def folder_taxa_counts(folder_path):
    path = Path(folder_path)
    counts = {}
    
    df = pd.read_csv('computational_performance.csv')
    df['num_taxa'] = df['msa_file'].apply(lambda x: taxa_count(path / f"{x}.clustal"))
    df.to_csv('computational_performance_with_taxa.csv', index=False)

folder_taxa_counts(r'processed_msa')

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