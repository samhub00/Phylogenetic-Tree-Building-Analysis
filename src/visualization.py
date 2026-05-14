import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.ticker as ticker


def time_heatmaps(eval_metrics):
    upgma_metrics = eval_metrics[eval_metrics['tree'] == 'UPGMA']
    nj_metrics = eval_metrics[eval_metrics['tree'] == 'NJ']
    ml_metrics = eval_metrics[eval_metrics['tree'] == 'ML']
    mp_metrics = eval_metrics[eval_metrics['tree'] == 'MP']

    taxa_bins = [0, 15, 65, 120, 400, 1010]
    sequence_bins = [0, 100, 500, 3000, 10000]

    # Binning the data for UPGMA to create a heatmap
    upgma_metrics['taxa_bin'] = pd.cut(upgma_metrics['taxa_count'], bins=taxa_bins)
    upgma_metrics['seq_len_bin'] = pd.cut(upgma_metrics['seq_len'], bins=sequence_bins)
    upgma_binned = upgma_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['cpu_time_seconds'].mean().reset_index()
    pivot_upgma = upgma_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="cpu_time_seconds")

    # Binning the data for NJ to create a heatmap
    nj_metrics['taxa_bin'] = pd.cut(nj_metrics['taxa_count'], bins=taxa_bins)
    nj_metrics['seq_len_bin'] = pd.cut(nj_metrics['seq_len'], bins=sequence_bins)
    nj_binned = nj_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['cpu_time_seconds'].mean().reset_index()
    pivot_nj = nj_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="cpu_time_seconds")

    # Binning the data for ML to create a heatmap
    ml_metrics['taxa_bin'] = pd.cut(ml_metrics['taxa_count'], bins=taxa_bins)
    ml_metrics['seq_len_bin'] = pd.cut(ml_metrics['seq_len'], bins=sequence_bins)
    ml_binned = ml_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['cpu_time_seconds'].mean().reset_index()
    pivot_ml = ml_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="cpu_time_seconds")

    # Binning the data for MP to create a heatmap
    mp_metrics['taxa_bin'] = pd.cut(mp_metrics['taxa_count'], bins=[0,30,60,90,120])
    mp_metrics['seq_len_bin'] = pd.cut(mp_metrics['seq_len'], bins=[0,60,130,400,1500,5000])
    mp_binned = mp_metrics.groupby(['taxa_bin', 'seq_len_bin'], observed=False)['cpu_time_seconds'].mean().reset_index()
    pivot_mp = mp_binned.pivot(index="taxa_bin", columns="seq_len_bin", values="cpu_time_seconds")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    pivots = [pivot_upgma, pivot_nj, pivot_ml, pivot_mp]
    titles = [
        "UPGMA Execution Time (seconds)",
        "NJ Execution Time (seconds)",
        "ML Execution Time (seconds)",
        "MP Execution Time (seconds)"
    ]

    #Displays four heatmaps, one for each of the tree buidling algorithms
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

        #axes[i].set_xticklabels([])
        #axes[i].set_yticklabels([])

    plt.tight_layout()
    plt.show()

def time_scatter(df):
    """
    This method creates a scatter plot of time to process over taxa count with dot size signifying
    sequence length
    """
    plt.figure(figsize=(12, 8))

    sns.scatterplot(
        data=df,
        x='taxa_count',
        y='cpu_time_seconds',
        hue='tree',
        alpha=0.6,
        size='seq_len',
        sizes=(20,200)
    )

    plt.title("Computational Performance: Taxa Count vs. Execution Time", fontsize=15, pad=15)
    plt.xlabel("Number of Taxa", fontsize=12)
    plt.ylabel("CPU Time (Seconds)", fontsize=12)

    plt.yscale('log', base=4)
    plt.xscale('log', base=4)

    plt.legend(title='Tree Type', loc='upper left')
    plt.grid(True, which='both', linestyle='--', alpha=0)
    plt.tight_layout()
    plt.show()


def time_scatter2x2(df):
    """
    This shows  the analysis of all the different algorithms in one table by displaying all algorithsm
    and adding a line of best fit.

    KEEP IN MIND: the scales are logarithmic so the slop of the line is exponential
    """
    # Setup the figure and a 2x2 grid of axes
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
    
    # Flatten the 2D axes array to 1D for easy looping
    axes = axes.flatten()
    
    # Identify unique algorithms (e.g., UPGMA, NJ, ML, MP)
    algorithms = ['UPGMA', 'NJ', 'ML', 'MP'] 

    colors = ['red', 'green', 'blue', 'orange']

    for i, algo in enumerate(algorithms):
        ax = axes[i]
        color = colors[i]

        # Filter data for just this specific algorithm
        subset = df[df['tree'] == algo]
        subset = subset[(subset['taxa_count'] > 0) & (subset['cpu_time_seconds'] > 0)]
        
        # 2. Draw the scatterplot on the specific axis (ax=ax)
        sns.scatterplot(
            data=subset,
            x='taxa_count',
            y='cpu_time_seconds',
            color=color,
            alpha=0.6,
            size='seq_len',
            sizes=(20, 200),
            ax=ax  # This tells seaborn which subplot to use
        )

        if not subset.empty:
            log_x = np.log10(subset['taxa_count'])
            log_y = np.log10(subset['cpu_time_seconds'])

            m, b = np.polyfit(log_x, log_y, 1)

            x_vals = np.linspace(log_x.min(), log_x.max(), 100)
            y_fit = 10**(m * x_vals + b)

            ax.plot(10**x_vals, y_fit, color=color, linestyle='--',
                    linewidth=1, label=f'Best Fit Slope: {m:.2f}')

        # Individual Subplot Formatting
        ax.set_title(f"Algorithm: {algo}", fontsize=14)
        ax.set_xlabel("Number of Taxa", fontsize=10)
        ax.set_ylabel("CPU Time (Seconds)", fontsize=10)

        # apply log scales for data clarity
        ax.set_yscale('log')
        ax.set_xscale('log')

        ax.set_xticks([10, 50, 100, 250, 500])
        
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

        ax.xaxis.get_major_formatter().set_scientific(False)
        ax.yaxis.get_major_formatter().set_scientific(False)

        ax.set_title(f"Algorithm: {algo}", color=color, fontweight='bold')
        ax.legend(loc='upper left', fontsize='small')
        ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.4)

    # global formatting
    fig.suptitle("Computational Performance with Short Sequences: Taxa Count vs. Execution Time", fontsize=18, y=1.02)
    plt.tight_layout()
    plt.savefig('short_seq_performance_2x2_plot.png', bbox_inches='tight')
    plt.show()

def ram_scatter2x2(df):
    """
    This shows  the analysis of all the different algorithms in one table by displaying all algorithsm
    and adding a line of best fit.

    Analyzes the RAM usage over the taxa count
    """
    # Setup the figure and a 2x2 grid of axes
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
    
    # Flatten the 2D axes array to 1D for easy looping
    axes = axes.flatten()
    
    # Identify unique algorithms (e.g., UPGMA, NJ, ML, MP)
    algorithms = ['UPGMA', 'NJ', 'ML', 'MP'] 

    colors = ['red', 'green', 'blue', 'orange']

    for i, algo in enumerate(algorithms):
        ax = axes[i]
        color = colors[i]

        # Filter data for just this specific algorithm
        subset = df[df['tree'] == algo]
        subset = subset[(subset['taxa_count'] > 0) & (subset['cpu_time_seconds'] > 0)]
        
        # 2. Draw the scatterplot on the specific axis (ax=ax)
        sns.scatterplot(
            data=subset,
            x='taxa_count',
            y='peak_ram_mb',
            color=color,
            alpha=0.6,
            size='seq_len',
            sizes=(20, 200),
            ax=ax  # This tells seaborn which subplot to use
        )

        if not subset.empty:
            log_x = np.log10(subset['taxa_count'])
            log_y = np.log10(subset['peak_ram_mb'])

            m, b = np.polyfit(log_x, log_y, 1)

            x_vals = np.linspace(log_x.min(), log_x.max(), 100)
            y_fit = 10**(m * x_vals + b)

            ax.plot(10**x_vals, y_fit, color=color, linestyle='--',
                    linewidth=1, label=f'Best Fit Slope: {m:.2f}')

        # Individual Subplot Formatting
        ax.set_title(f"Algorithm: {algo}", fontsize=14)
        ax.set_xlabel("Number of Taxa", fontsize=10)
        ax.set_ylabel("RAM usage (mb)", fontsize=10)

        # apply log scales for data clarity
        ax.set_yscale('log')
        ax.set_xscale('log')

        ax.set_xticks([10, 50, 100, 250, 500])
        
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

        ax.xaxis.get_major_formatter().set_scientific(False)
        ax.yaxis.get_major_formatter().set_scientific(False)

        ax.set_title(f"Algorithm: {algo}", color=color, fontweight='bold')
        ax.legend(loc='upper left', fontsize='small')
        ax.grid(True, which='major', axis='y', linestyle='--', alpha=0.4)

    # global formatting
    fig.suptitle("RAM Usage with Short Sequences: Taxa Count vs. Execution Time", fontsize=18, y=1.02)
    plt.tight_layout()
    plt.savefig('short_ram_performance_2x2_plot.png', bbox_inches='tight')
    plt.show()


def time_scatter_seq(df):
    plt.figure(figsize=(12, 8))

    sns.scatterplot(
        data=df,
        x='seq_len',
        y='cpu_time_seconds',
        hue='tree',
        alpha=0.6,
        size='taxa_count',
        sizes=(20,300)
    )

    plt.title("Computational Performance: Sequence Length vs. Execution Time", fontsize=15, pad=15)
    plt.xlabel("Sequence Length", fontsize=12)
    plt.ylabel("CPU Time (Seconds)", fontsize=12)

    plt.yscale('log')
    plt.xscale('log', base=7)

    plt.legend(title='Tree Type', loc='upper left')
    plt.grid(True, which='both', linestyle='--', alpha=0)
    plt.tight_layout()
    plt.show()


def acc_tree_heatmaps(df):
    # 1. Calculate Accuracy
    df['accuracy'] = 1 - df['nRF_distance']

    # 2. Create Bins for Taxa Count and Sequence Length
    # Using 'pd.cut' creates equal-width bins. 
    # You can also pass specific lists like [0, 50, 100, 500] if you prefer.
    df['taxa_bin'] = pd.cut(df['taxa_count'], bins=[0, 20, 40, 80, 150, 600])
    df['seq_bin'] = pd.cut(df['seq_len'], bins=[0, 100, 300, 800, 1500, 6000])

    # 3. Setup the visualization grid
    methods = ['ML', 'NJ', 'UPGMA', 'MP']
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    for i, method in enumerate(methods):
        subset = df[df['method'] == method]
        
        if subset.empty:
            axes[i].set_title(f"Method: {method} (No Data)")
            continue

        # 4. Pivot using the BINS instead of raw values
        pivot_df = subset.pivot_table(
            index='taxa_bin', 
            columns='seq_bin', 
            values='accuracy', 
            aggfunc='mean'
        )

        # Sort indices: High taxa at top, low sequence length on left
        pivot_df = pivot_df.sort_index(ascending=False)

        # 5. Generate Heatmap
        sns.heatmap(
            pivot_df, 
            annot=True, 
            fmt=".2f", 
            cmap='RdYlGn', 
            ax=axes[i],
            cbar_kws={'label': 'Mean Accuracy'}
        )
        
        fig.suptitle('Accuracy Scores Across Algorithms High Divergence Trees', fontsize=22)


        axes[i].set_title(f'Algorithm: {method}', fontweight='bold', fontsize=10)
        axes[i].set_xlabel('Sequence Length Range (bp)')
        axes[i].set_ylabel('Number of Taxa Range')

    plt.tight_layout()
    plt.savefig('high_div_accuracy_across_algorithms.png')
    #plt.show()

#acc_tree_heatmaps(low_div_df)

def corr_tree_heatmaps(df):
    # 1. Calculate Accuracy
    df['accuracy'] = df['p_corr']
    # 2. Create Bins for Taxa Count and Sequence Length
    # Using 'pd.cut' creates equal-width bins. 
    # You can also pass specific lists like [0, 50, 100, 500] if you prefer.
    df['taxa_bin'] = pd.cut(df['taxa_count'], bins=[0, 20, 40, 80, 150, 600])
    df['seq_bin'] = pd.cut(df['seq_len'], bins=[0, 100, 300, 800, 1500, 6000])

    # 3. Setup the visualization grid
    methods = ['ML', 'NJ', 'UPGMA', 'MP']
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    for i, method in enumerate(methods):
        subset = df[df['method'] == method]
        
        if subset.empty:
            axes[i].set_title(f"Method: {method} (No Data)")
            continue

        # 4. Pivot using the BINS instead of raw values
        pivot_df = subset.pivot_table(
            index='taxa_bin', 
            columns='seq_bin', 
            values='accuracy', 
            aggfunc='mean'
        )

        # Sort indices: High taxa at top, low sequence length on left
        pivot_df = pivot_df.sort_index(ascending=False)

        # 5. Generate Heatmap
        sns.heatmap(
            pivot_df, 
            annot=True, 
            fmt=".2f", 
            cmap='RdYlGn', 
            ax=axes[i],
            cbar_kws={'label': 'Mean Accuracy'}
        )
        
        fig.suptitle('Patristic Correlation Scores on High Divergence Trees', fontsize=24, fontweight='bold')

        axes[i].set_title(f'Algorithm: {method}', fontsize=10)
        axes[i].set_xlabel('Sequence Length Range (bp)')
        axes[i].set_ylabel('Number of Taxa Range')

    plt.tight_layout()
    plt.savefig('high_divergence_pcorr_heatmap.png')
    plt.show()

