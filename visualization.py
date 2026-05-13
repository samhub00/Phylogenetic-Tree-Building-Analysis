import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

emp_eval = pd.read_csv(r'CSVs/computational_performance_with_taxa.csv')
#print(emp_eval.head())
emp_eval = emp_eval.rename(columns={'num_taxa' : 'taxa_count'})

sim_eval = pd.read_csv('ROSE_performance.csv')

eval_metrics = pd.concat([emp_eval, sim_eval], ignore_index=True)

#print(eval_metrics.head())

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

#time_heatmaps(eval_metrics)
#time_scatter(eval_metrics)


def acc_tree_heatmaps(df):
    # 1. Parse the 'name' column to extract Taxa Amount and Sequence Length
    # Format: HKY_high_100_100_ML -> ['HKY', 'high', '100', '100', 'ML']
    def parse_metadata(name):
        parts = name.split('_')
        return int(parts[2]), int(parts[3])

    df[['taxa', 'length']] = df.apply(
        lambda r: parse_metadata(r['name']), axis=1, result_type='expand'
    )

    # 2. Calculate Accuracy (1 - normalized RF distance)
    df['accuracy'] = 1 - df['nRF_distance']

    # 3. Setup the visualization grid (2x2 for ML, NJ, UPGMA, MP)
    methods = ['ML', 'NJ', 'UPGMA', 'MP']
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for i, method in enumerate(methods):
        # Filter data for specific algorithm
        subset = df[df['method'] == method]
        
        if subset.empty:
            axes[i].set_title(f"Method: {method} (No Data)")
            continue

        # Pivot data: Index = Taxa, Columns = Length, Values = Mean Accuracy
        pivot_df = subset.pivot_table(
            index='taxa', 
            columns='length', 
            values='accuracy', 
            aggfunc='mean'
        )

        # Sort indices to ensure heatmap axes are logical
        pivot_df = pivot_df.sort_index(ascending=False) # High taxa at top
        pivot_df = pivot_df.reindex(sorted(pivot_df.columns), axis=1)

        # 4. Generate Heatmap
        sns.heatmap(
            pivot_df, 
            annot=True, 
            fmt=".2f", 
            cmap='RdYlGn', 
            ax=axes[i],
            cbar_kws={'label': 'Accuracy (1 - nRF)'}
        )
        axes[i].set_title(f'Algorithm Performance: {method}', fontweight='bold')
        axes[i].set_xlabel('Sequence Length (bp)')
        axes[i].set_ylabel('Number of Taxa')

    plt.tight_layout()
    plt.show()

"""csv_paths = [
    r'built_tree_analysis.csv',
    r'CSVs\built_tree_analysis.csv'
]

combined_df = pd.concat([pd.read_csv(f) for f in csv_paths], ignore_index=True)

combined_df.to_csv('tree_analysis.csv', index = False)
"""
# To use:
#df = pd.read_csv('CSVs/BuiltTreeAnalysis.csv')

#df_unique = df.drop_duplicates()
#acc_tree_heatmaps(df_unique)

#mp_current = df_unique[df_unique['method'] == "MP"]
#print(mp_current)
"""
df = pd.read_csv('CSVs/BuiltTreeAnalysis.csv')

# 2. Extract the Model and Rate prefix (e.g., 'HKY_high') from the 'name' column
# We split by '_' and take the first two elements
df['group'] = df['name'].apply(lambda x: '_'.join(x.split('_')[:2]))

# 3. Get unique groups and methods to iterate through
groups = df['group'].unique()      # ['HKY_high', 'HKY_low', 'JC_high', 'JC_low']
methods = df['method'].unique()    # ['ML', 'NJ', 'UPGMA', 'MP']

# 4. Separate and save each combination to a CSV file
for g in groups:
    for m in methods:
        # Filter the dataframe for the specific group and method
        subset = df[(df['group'] == g) & (df['method'] == m)]
        
        # Only save if the subset is not empty
        if not subset.empty:
            filename = f"{g}_{m}.csv"
            subset.to_csv(filename, index=False)
            print(f"Saved: {filename}")



"""
df = pd.read_csv('CSVs/BuiltTreeAnalysis.csv')
df['group'] = df['name'].apply(lambda x: '_'.join(x.split('_')[:2]))
df['ParamA'] = df['name'].apply(lambda x: x.split('_')[2]).astype(int)
df['ParamB'] = df['name'].apply(lambda x: x.split('_')[3]).astype(int)

groups = sorted(df['group'].unique())   # HKY_high, HKY_low, JC_high, JC_low
methods = sorted(df['method'].unique()) # ML, MP, NJ, UPGMA


fig, axes = plt.subplots(len(groups), len(methods), figsize=(22, 18))
fig.suptitle('nRF Distance Heatmaps: Comparative Analysis', fontsize=24, y=1.02)

for i, g in enumerate(groups):
    for j, m in enumerate(methods):
        ax = axes[i, j]
        subset = df[(df['group'] == g) & (df['method'] == m)]
        
        if not subset.empty:
            # Pivot data for the heatmap
            pivot_df = subset.pivot(index='ParamA', columns='ParamB', values='nRF_distance')
            sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt=".3f", ax=ax, cbar=(j==3))
            ax.set_title(f"{g} | {m}", fontsize=14)
        else:
            # Handle the missing HKY_low_MP data
            ax.text(0.5, 0.5, 'NO DATA', ha='center', va='center', color='red')
            ax.set_title(f"{g} | {m} (Missing)")
        
        # Clean up labels to prevent crowding
        if i < len(groups) - 1: ax.set_xlabel('')
        if j > 0: ax.set_ylabel('')

plt.tight_layout()
plt.savefig('all_heatmaps_grid.png', bbox_inches='tight')
plt.show()