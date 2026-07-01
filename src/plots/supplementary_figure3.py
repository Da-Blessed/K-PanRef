#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%%
def plot_sample_length(input_path, output_path):
    df = pd.read_csv(input_path, sep="\t")
    df = df[df["Graph"] != "CPC-HPRC"].reset_index(drop=True)
    df["Length_Gb"] = df["Length"] / 1e9

    order = ["KPanRef", "CPC", "HPRC"]
    df["Graph"] = pd.Categorical(df["Graph"], categories=order, ordered=True)
    graphs = order

    display_names = {
        "KPanRef": "K-PanRef",
        "CPC":     "CPC",
        "HPRC":    "HPRC",
    }

    n_graphs = len(graphs)

    fig = plt.figure(figsize=(15, 15))
    fig.patch.set_facecolor('white')

    gs = fig.add_gridspec(n_graphs, 3, hspace=0.5)
    axes = [fig.add_subplot(gs[i, :]) for i in range(n_graphs)]

    colors = {
        "KPanRef": "#0047A0",
        "CPC":     "#CD2E3A",
        "HPRC":    "#000000",
    }

    for ax, graph in zip(axes, graphs):
        df_g = df[df["Graph"] == graph].copy()
        df_g = df_g.sort_values("Sample")

        x = np.arange(len(df_g))
        color = colors.get(graph, "#555555")

        ax.bar(x, df_g["Length_Gb"], color=color + '99', edgecolor=color, linewidth=1.2, width=0.6)

        ax.set_title(display_names.get(graph, graph), fontsize=14, fontweight='bold', color='black')
        ax.set_xticks(x)
        ax.set_xticklabels(df_g["Sample"], rotation=90, fontsize=11)
        ax.set_ylabel("Length (Gb)", fontsize=13, fontweight='bold')
        ax.set_xlim(-0.5, len(df_g) - 0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.tick_params(axis='x', colors='black', labelsize=11)
        ax.tick_params(axis='y', colors='black', labelsize=13)
        ax.yaxis.grid(False)
        ax.set_facecolor('white')

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    input_path  = $PATH
    output_path = $PATH
    plot_sample_length(input_path, output_path)
