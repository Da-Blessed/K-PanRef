#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

#%%
def plot_chrom_stats(input_path, output_path):
    df = pd.read_csv(input_path, sep="\t")
    df = df[df["Graph"] != "CPC-HPRC"].reset_index(drop=True)

    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    df["Chrom"] = pd.Categorical(df["Chrom"], categories=chrom_order, ordered=True)
    df = df.sort_values(["Chrom", "Graph"]).reset_index(drop=True)

    graphs = ["KPanRef", "CPC", "HPRC"]
    graphs = [g for g in graphs if g in df["Graph"].unique()]

    display_names = {
        "KPanRef": "K-PanRef",
        "CPC":     "CPC",
        "HPRC":    "HPRC",
    }

    chrom_order_no_m = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chroms = [c for c in chrom_order_no_m if c in df["Chrom"].values]
    chroms_m = ["chrM"]

    colors = {
        "KPanRef": "#0047A0",
        "CPC":     "#CD2E3A",
        "HPRC":    "#000000",
    }

    metrics = ["Node", "Edge", "Length"]
    n_graphs = len(graphs)
    width = 0.8 / n_graphs

    fig = plt.figure(figsize=(20, 15))
    fig.patch.set_facecolor('white')

    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.4, wspace=0.4)

    def plot_bars(ax, df, chrom_list, metric, scale, ylabel, graphs, colors, width, show_legend=False, show_xlabel=False, is_main=False):
        x = np.arange(len(chrom_list))
        for i, graph in enumerate(graphs):
            df_g = df[df["Graph"] == graph]
            vals = []
            for chrom in chrom_list:
                row = df_g[df_g["Chrom"] == chrom]
                vals.append(row[metric].values[0] / scale if len(row) > 0 else 0)
            color = colors.get(graph, "#555555")
            offset = (i - n_graphs / 2 + 0.5) * width
            ax.bar(x + offset, vals, width=width * 0.9,
                   color=color + '99', edgecolor=color, linewidth=1.0, label=display_names.get(graph, graph))

        ax.set_title(metric, fontsize=18, fontweight='bold', color='black')
        ax.set_ylabel(ylabel, fontsize=15, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("chr", "") for c in chrom_list], fontsize=11)
        ax.set_facecolor('white')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.tick_params(axis='x', colors='black', labelsize=11)
        ax.tick_params(axis='y', colors='black', labelsize=12)
        ax.yaxis.grid(False)
        if is_main:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))
        else:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
        if show_legend:
            ax.legend(fontsize=11, frameon=True, loc='upper left', bbox_to_anchor=(1.05, 1.0))
        if show_xlabel:
            ax.set_xlabel("Chromosome", fontsize=13, fontweight='bold')

    for row_idx, metric in enumerate(metrics):
        if metric == "Length":
            scale = 1e6
            ylabel = "Mb"
            scale_m = 1e3
            ylabel_m = "kb"
        else:
            max_val = df[df["Chrom"] != "chrM"][metric].max()
            max_val_m = df[df["Chrom"] == "chrM"][metric].max()

            if max_val >= 1e6:
                scale = 1e6
                ylabel = "Count (M)"
            elif max_val >= 1e3:
                scale = 1e3
                ylabel = "Count (k)"
            else:
                scale = 1
                ylabel = "Count"

            if max_val_m >= 1e6:
                scale_m = 1e6
                ylabel_m = "Count (M)"
            elif max_val_m >= 1e3:
                scale_m = 1e3
                ylabel_m = "Count (k)"
            else:
                scale_m = 1
                ylabel_m = "Count"

        ax_main = fig.add_subplot(gs[row_idx, 0:3])
        plot_bars(ax_main, df, chroms, metric, scale, ylabel, graphs, colors, width,
                  show_legend=False, show_xlabel=(row_idx == 2), is_main=True)

        ax_m = fig.add_subplot(gs[row_idx, 3])
        plot_bars(ax_m, df, chroms_m, metric, scale_m, ylabel_m, graphs, colors, width=0.25,
                  show_legend=(row_idx == 0), show_xlabel=(row_idx == 2), is_main=False)

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    input_path  = $PATH
    output_path = $PATH
    plot_chrom_stats(input_path, output_path)
