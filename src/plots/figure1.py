#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np

#%%
def parse_panacus(panacus_path):
    df = pd.read_csv(panacus_path, sep="\t", skiprows=6, header=None)
    df.columns = ['Sample', 'cov27', 'cov1', 'cov3']
    df[['cov27', 'cov1', 'cov3']] = df[['cov27', 'cov1', 'cov3']].astype(float)
    df = df[['Sample', 'cov27', 'cov3', 'cov1']].copy()
    return df

#%%
def plot_combined(graph_stats_path, panacus_path, output_path):
    df_graph = pd.read_csv(graph_stats_path, sep="\t")
    df_graph = df_graph[df_graph["Graph"] != "CPC-HPRC"].reset_index(drop=True)

    order = ["KPanRef", "CPC", "HPRC"]
    df_graph["Graph"] = pd.Categorical(df_graph["Graph"], categories=order, ordered=True)
    df_graph = df_graph.sort_values("Graph").reset_index(drop=True)

    colors = {
        "KPanRef":  "#0047A0",
        "CPC":      "#CD2E3A",
        "HPRC":     "#000000",
    }

    display_names = {
        "KPanRef": "K-PanRef",
        "CPC":     "CPC",
        "HPRC":    "HPRC",
    }

    metrics = ["Sample Size", "Node", "Edge", "Length"]
    scales  = [1, 1e6, 1e6, 1e9]
    ylabels = ["Count", "Count (M)", "Count (M)", "Count (G)"]

    fig = plt.figure(figsize=(18, 20))
    fig.patch.set_facecolor('white')

    gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.45, wspace=0.35)

    axes_top = []
    for idx, (metric, scale, ylabel) in enumerate(zip(metrics, scales, ylabels)):
        ax = fig.add_subplot(gs[0, idx])
        ax.set_facecolor('white')
        axes_top.append(ax)

        vals = df_graph[metric].astype(float) / scale
        x = np.arange(len(df_graph))
        bar_colors  = [colors.get(g, "#999999") + '99' for g in df_graph["Graph"]]
        edge_colors = [colors.get(g, "#999999") for g in df_graph["Graph"]]

        ax.bar(x, vals, color=bar_colors, edgecolor=edge_colors, linewidth=1.2, width=0.5)

        ax.set_title(metric, fontsize=14, fontweight='bold', color='black', pad=15)
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([display_names.get(g, g) for g in df_graph["Graph"]], fontsize=10, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.tick_params(colors='black', labelsize=11)
        ax.yaxis.grid(False)

        if metric == "Sample Size":
            for bar, val in zip(ax.patches, vals):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + vals.max() * 0.02,
                        f'{int(val)}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        elif metric == "Length":
            for bar, val in zip(ax.patches, vals):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + (vals.max() - vals.min()) * 0.05,
                        f'{val:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        else:
            for bar, val in zip(ax.patches, vals):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + vals.max() * 0.02,
                        f'{val:.1f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        if metric == "Length":
            margin = (vals.max() - vals.min()) * 0.3
            ax.set_ylim(vals.min() - margin, vals.max() + margin)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.2f}'))

    ax_pan = fig.add_subplot(gs[1:4, :])
    ax_pan.set_facecolor('white')

    df_pan = parse_panacus(panacus_path)
    x = np.arange(len(df_pan))
    w = 0.8

    core      = (df_pan['cov27'] - df_pan['cov27'].iloc[0]) / 1e6
    common    = ((df_pan['cov3'] - df_pan['cov27']) - (df_pan['cov3'].iloc[0] - df_pan['cov27'].iloc[0])) / 1e6
    singleton = ((df_pan['cov1'] - df_pan['cov3']) - (df_pan['cov1'].iloc[0] - df_pan['cov3'].iloc[0])) / 1e6

    color_core      = "#08306B"  
    color_common    = "#0047A0"  
    color_singleton = "#4292C6"  

    ax_pan.bar(x, core,      width=w, color=color_core,
               edgecolor=color_core, linewidth=1.5, zorder=2)
    ax_pan.bar(x, common,    width=w, color=color_common + '99',
               bottom=core, edgecolor=color_common, linewidth=1.5, zorder=2)
    ax_pan.bar(x, singleton, width=w, color=color_singleton + '99',
               bottom=core + common, edgecolor=color_singleton, linewidth=1.5, zorder=2)

    legend_elements = [
        mpatches.Patch(facecolor=color_singleton + '99', edgecolor=color_singleton, label='Singleton'),
        mpatches.Patch(facecolor=color_common + '99', edgecolor=color_common, label='Common (≥10%)'),
        mpatches.Patch(facecolor=color_core + '99', edgecolor=color_core, label='Core (≥95%)'),
    ]

    ax_pan.set_title("Pangenome Growth", fontsize=14, fontweight='bold', color='black', pad=15)
    ax_pan.set_ylabel("Added Sequences (Mb)", fontsize=12, fontweight='bold')
    ax_pan.set_xlabel("Samples", fontsize=12, fontweight='bold')
    ax_pan.set_xticks([])
    ax_pan.margins(x=0.01)
    ax_pan.spines['top'].set_visible(False)
    ax_pan.spines['right'].set_visible(False)
    ax_pan.spines['left'].set_color('black')
    ax_pan.spines['bottom'].set_color('black')
    ax_pan.tick_params(colors='black', labelsize=10)
    ax_pan.yaxis.grid(False)
    ax_pan.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))
    ax_pan.legend(handles=legend_elements, fontsize=11, frameon=True, loc='upper left')

    fig.canvas.draw()
    label_axes = [axes_top[0], ax_pan]
    labels     = ['A', 'B']
    for ax, label in zip(label_axes, labels):
        bbox = ax.get_position()
        fig.text(bbox.x0 - 0.02, bbox.y1 + 0.01, label, fontsize=16, fontweight='bold',
                 color='black', va='bottom', ha='left', transform=fig.transFigure)

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    graph_stats_path = $PATH
    panacus_path = $PATH
    output_path  = $PATH
    plot_combined(graph_stats_path, panacus_path, output_path)
