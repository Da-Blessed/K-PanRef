#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%%
def style_boxplot(bp, color):
    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        plt.setp(bp[element], color=color, linewidth=1.5)
    plt.setp(bp['fliers'], marker='o', markerfacecolor=color,
             markeredgecolor=color, markersize=4, alpha=0.5)
    for patch in bp['boxes']:
        patch.set_facecolor(color + '99')

#%%
def plot_assembly_stats(input_path, output_path):
    df = pd.read_csv(input_path, sep="\t")

    metrics = [
        "Genome Size (Gb)", "N50 (Mb)", "Number of Scaffolds",
        "QV", "N's per 100kb", "Completeness"
    ]
    ylabels = ["Gb", "Mb", "Count", "Value", "Count", "%"]
    formats = [".2f", "d", "d", "d", ".1f", "d"]
    titles  = [
        "Genome Size", "N50", "Number of Scaffolds",
        "QV", "N's per 100kb", "Completeness"
    ]

    n_cols = 3
    n_rows = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(wspace=0.4, hspace=0.5)

    for idx, (metric, ylabel, fmt, title) in enumerate(zip(metrics, ylabels, formats, titles)):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        ax.set_facecolor('white')

        color_box = "#CD2E3A" if row == 0 else "#0047A0"

        data = df[metric].dropna()
        data = data[np.isfinite(data)]

        bp = ax.boxplot(data, widths=0.5, patch_artist=True)
        style_boxplot(bp, color_box)

        data_min    = data.min()
        data_max    = data.max()
        data_median = data.median()
        margin = data_max - data_median
        ax.set_ylim(data_min - margin, None)

        ax.set_title(title, fontsize=14, fontweight='bold', color='black', pad=15)
        ax.set_xticks([])
        ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.tick_params(colors='black', labelsize=13)
        ax.yaxis.grid(False)

        if fmt == "d":
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _, fmt=fmt: f'{int(round(x))}'))
        else:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _, fmt=fmt: f'{x:{fmt}}'))

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    input_path  = $PATH
    output_path = $PATH
    plot_assembly_stats(input_path, output_path)
