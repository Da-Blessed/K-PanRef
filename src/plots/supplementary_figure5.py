#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%%
def plot_sample_chrom_length(input_path, output_path):
    df = pd.read_csv(input_path, sep="\t")
    df = df[df["Graph"] != "CPC-HPRC"].reset_index(drop=True)

    graphs = ["KPanRef", "CPC", "HPRC"]
    graphs = [g for g in graphs if g in df["Graph"].unique()]

    display_names = {
        "KPanRef": "K-PanRef",
        "CPC":     "CPC",
        "HPRC":    "HPRC",
    }

    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    chroms = [c for c in chrom_order if c in df["Chrom"].unique()]

    colors = {
        "KPanRef": "#0047A0",
        "CPC":     "#CD2E3A",
        "HPRC":    "#000000",
    }

    n_rows = len(chroms)
    n_cols = len(graphs)

    col_widths = []
    for graph in graphs:
        max_s = max(len(df[(df["Graph"] == graph) & (df["Chrom"] == c)]) for c in chroms)
        col_widths.append(max(max_s, 1))

    chrom_max = {chrom: df[df["Chrom"] == chrom]["Length"].max() for chrom in chroms}

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(15, 20),
        gridspec_kw={'width_ratios': col_widths}
    )
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(hspace=0.4, wspace=0.1, bottom=0.04, left=0.07)

    if n_rows == 1:
        axes = [axes]
    if n_cols == 1:
        axes = [[ax] for ax in axes]

    for col_idx, graph in enumerate(graphs):
        df_g = df[df["Graph"] == graph]
        color = colors.get(graph, "#555555")

        for row_idx, chrom in enumerate(chroms):
            ax = axes[row_idx][col_idx]
            df_c = df_g[df_g["Chrom"] == chrom].sort_values("Sample")

            x = np.arange(len(df_c))
            ax.bar(x, df_c["Length"], color=color + '99', edgecolor=color, linewidth=1.0, width=0.6)

            ax.set_facecolor('white')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_color('black')
            ax.spines['bottom'].set_color('black')
            ax.tick_params(axis='y', colors='black', labelsize=9)
            ax.yaxis.grid(False)

            ax.set_xticks([])
            ax.set_xlim(-0.5, len(df_c) - 0.5)
            ax.set_ylim(0, chrom_max[chrom] * 1.1)
            ax.set_ylabel("")

            max_val = chrom_max[chrom]
            if max_val >= 1e6:
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1e6)}'))
            elif max_val >= 1e3:
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1e3)}'))
            else:
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))

            if col_idx != 0:
                ax.set_yticks([])
                ax.spines['left'].set_visible(False)

            if row_idx == 0:
                ax.set_title(display_names.get(graph, graph), fontsize=13, fontweight='bold', color='black')

    fig.canvas.draw()
    for row_idx, chrom in enumerate(chroms):
        ax = axes[row_idx][0]
        bbox = ax.get_position()
        y_center = (bbox.y0 + bbox.y1) / 2

        max_val = chrom_max[chrom]
        unit = "Mb" if max_val >= 1e6 else ("kb" if max_val >= 1e3 else "bp")

        fig.text(0.035, y_center, f"{chrom}\n({unit})", fontsize=9, fontweight='bold',
                 ha='center', va='center', color='black',
                 transform=fig.transFigure)

    fig.text(0.5, 0.02, "Sample", fontsize=13, fontweight='bold',
             ha='center', va='bottom', color='black')

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    input_path  = $PATH
    output_path = $PATH
    plot_sample_chrom_length(input_path, output_path)
