import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

eval_metrics = pd.read_csv('computational_performance_with_taxa.csv')

def time_heatmaps(eval_metrics):
    upgma_metrics = eval_metrics[eval_metrics['tree'] == 'UPGMA']
    nj_metrics = eval_metrics[eval_metrics['tree'] == 'NJ']
    ml_metrics = eval_metrics[eval_metrics['tree'] == 'ML']
    mp_metrics = eval_metrics[eval_metrics['tree'] == 'MP']

    # Binning the data for UPGMA to create a heatmap
    upgma_metrics['taxa_bin'] = pd.cut(upgma_metrics['num_taxa'], bins=5)
    upgma_metrics['seq_len_bin'] = pd.cut(upgma_metrics['seq_len'], bins=5)
    upgma_binned = upgma_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['execution_time_seconds'].mean().reset_index()
    pivot_upgma = upgma_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="execution_time_seconds")

    # Binning the data for NJ to create a heatmap
    nj_metrics['taxa_bin'] = pd.cut(nj_metrics['num_taxa'], bins=5)
    nj_metrics['seq_len_bin'] = pd.cut(nj_metrics['seq_len'], bins=5)
    nj_binned = nj_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['execution_time_seconds'].mean().reset_index()
    pivot_nj = nj_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="execution_time_seconds")

    # Binning the data for ML to create a heatmap
    ml_metrics['taxa_bin'] = pd.cut(ml_metrics['num_taxa'], bins=5)
    ml_metrics['seq_len_bin'] = pd.cut(ml_metrics['seq_len'], bins=5)
    ml_binned = ml_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['execution_time_seconds'].mean().reset_index()
    pivot_ml = ml_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="execution_time_seconds")

    # Binning the data for MP to create a heatmap
    mp_metrics['taxa_bin'] = pd.cut(mp_metrics['num_taxa'], bins=5)
    mp_metrics['seq_len_bin'] = pd.cut(mp_metrics['seq_len'], bins=5)
    mp_binned = mp_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['execution_time_seconds'].mean().reset_index()
    pivot_mp = mp_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="execution_time_seconds")


    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()


    pivots = [pivot_upgma, pivot_nj, pivot_ml, pivot_mp]
    titles = [
        "UPGMA Execution Time (seconds)",
        "NJ Execution Time (seconds)",
        "ML Execution Time (seconds)",
        "MP Execution Time (seconds)"
    ]


    for i in range(4):
        
        sns.heatmap(
            pivots[i], 
            ax=axes[i], 
            annot=True, 
            fmt=".1f", 
            cmap="YlGnBu",
            cbar_kws={'label': 'Seconds'} 
        )
        
        axes[i].set_title(titles[i], fontsize=14, fontweight='bold')
        axes[i].invert_yaxis()

        axes[i].set_xlabel('Sequence Length (Binned)')
        axes[i].set_ylabel('Number of Taxa (Binned)')

        axes[i].set_xticklabels([])
        axes[i].set_yticklabels([])

    plt.tight_layout()
    plt.show()