#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib_venn import venn2

#%%
def style_boxplot(bp, color):
    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        plt.setp(bp[element], color=color, linewidth=1.5)
    plt.setp(bp['fliers'], marker='o', markerfacecolor=color,
             markeredgecolor=color, markersize=4, alpha=0.5)
    for patch in bp['boxes']:
        patch.set_facecolor(color + '99')

#%%
def plot_variant_stats(per_sample_path, annot_paths, venn_path, output_path):
    order = ["KPanRef", "CPC", "HPRC"]

    display_names = {
        "KPanRef": "K-PanRef",
        "CPC":     "CPC",
        "HPRC":    "HPRC",
    }

    df_samp = pd.read_csv(per_sample_path, sep="\t")
    df_samp["SNV"] = df_samp["Transitions"] + df_samp["Transversions"]
    df_samp["InDel"] = df_samp["Indels"]
    df_samp = df_samp[df_samp["Graph Name"] != "CPC-HPRC"].reset_index(drop=True)
    df_samp["Graph Name"] = pd.Categorical(df_samp["Graph Name"], categories=order, ordered=True)
    df_samp = df_samp.sort_values("Graph Name").reset_index(drop=True)

    annot_dfs = {}
    for graph, path in annot_paths.items():
        df_annot = pd.read_csv(path, sep="\t", low_memory=False)
        df_annot = df_annot.drop_duplicates(subset="AnnotSV_ID", keep="first").reset_index(drop=True)
        df_annot["SV_length"] = df_annot["SV_length"].abs()
        annot_dfs[graph] = df_annot

    df_venn = pd.read_csv(venn_path, sep="\t")

    colors = {
        "KPanRef": "#0047A0",
        "CPC":     "#CD2E3A",
        "HPRC":    "#000000",
    }

    box_metrics = ["SNV", "InDel", "SV"]
    box_titles  = ["Per-sample SNVs", "Per-sample InDels", "Per-sample SVs"]
    box_ylabels = ["Count (k)", "Count (k)", "Count (k)"]
    box_scales  = [1e3, 1e3, 1e3]

    def fmt_val(v):
        if v >= 1e6:
            return f'{v/1e6:.1f}M'
        elif v >= 1e3:
            return f'{v/1e3:.1f}k'
        else:
            return f'{v:,}'

    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor('white')

    gs_top = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35,
                               top=0.95, bottom=0.42)
    gs_bot = gridspec.GridSpec(1, 3, figure=fig, wspace=0.3,
                               top=0.36, bottom=0.06)

    axes = np.empty((2, 3), dtype=object)
    for i in range(2):
        for j in range(3):
            axes[i, j] = fig.add_subplot(gs_top[i, j])

    ax_venn_small = fig.add_subplot(gs_bot[0, 0])
    ax_venn_sv    = fig.add_subplot(gs_bot[0, 1])

    for idx, (metric, title, ylabel, scale) in enumerate(zip(box_metrics, box_titles, box_ylabels, box_scales)):
        ax = axes[0, idx]
        ax.set_facecolor('white')

        data_list  = []
        graph_list = []
        for graph in order:
            if graph in df_samp["Graph Name"].values:
                vals = df_samp[df_samp["Graph Name"] == graph][metric] / scale
                data_list.append(vals.dropna().values)
                graph_list.append(graph)

        positions = np.arange(1, len(graph_list) + 1)
        for i, (data, graph) in enumerate(zip(data_list, graph_list)):
            color = colors.get(graph, "#555555")
            bp = ax.boxplot(data, positions=[positions[i]], widths=0.5, patch_artist=True)
            style_boxplot(bp, color)

        all_data = np.concatenate(data_list)
        data_min = np.min(all_data)
        data_max = np.max(all_data)
        data_median = np.median(all_data)
        margin = data_max - data_median
        ax.set_ylim(data_min - margin, None)

        ax.set_title(title, fontsize=14, fontweight='bold', color='black', pad=15)
        ax.set_xticks(positions)
        ax.set_xticklabels([display_names.get(g, g) for g in graph_list], fontsize=11, fontweight='bold')
        if idx == 0:
            ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        else:
            ax.set_ylabel("")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.tick_params(colors='black', labelsize=12)
        ax.yaxis.grid(False)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))

    for idx, graph in enumerate(["KPanRef", "CPC", "HPRC"]):
        ax = axes[1, idx]
        ax.set_facecolor('white')

        color = colors.get(graph, "#555555")
        sv_len_orig = annot_dfs[graph]["SV_length"].dropna()
        sv_len_orig = sv_len_orig[sv_len_orig > 0].copy()
        sv_len = sv_len_orig.clip(upper=10000)

        bins = np.arange(0, 10100, 100)
        ax.hist(sv_len, bins=bins, color=color + '99', edgecolor=color, linewidth=0.8)

        stats_text = (
            f"Max: {int(sv_len_orig.max()):,}\n"
            f"Mean: {int(sv_len_orig.mean()):,}\n"
            f"Median: {int(sv_len_orig.median()):,}\n"
            f"Min: {int(sv_len_orig.min()):,}"
        )
        ax.text(0.98, 0.98, stats_text,
                transform=ax.transAxes,
                fontsize=10, va='top', ha='right',
                color='black', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#CCCCCC', alpha=0.8))

        ax.set_title(f"{display_names.get(graph, graph)} SV Length Distribution", fontsize=14, fontweight='bold', color='black', pad=15)
        ax.set_xlabel("kb", fontsize=12, fontweight='bold')
        if idx == 0:
            ax.set_ylabel("Count (k)", fontsize=12, fontweight='bold')
        else:
            ax.set_ylabel("")
        ax.set_xticks(np.arange(0, 10100, 1000))
        ax.set_xticklabels(
            [str(int(b/1e3)) for b in np.arange(0, 10100, 1000)]
        )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.tick_params(colors='black', labelsize=12)
        ax.yaxis.grid(False)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1e3)}'))

    def draw_venn(ax, vtype, vtitle):
        ax.set_facecolor('white')
        row = df_venn[df_venn["Type"] == vtype].iloc[0]
        kpanref_only = int(row["KPanRef"])
        shared       = int(row["Shared"])
        cpchprc_only = int(row["CPC-HPRC"])

        v = venn2(
            subsets=(1, 1, 1),
            set_labels=("K-PanRef", "CPC-HPRC"),
            ax=ax
        )

        v.get_patch_by_id('10').set_color('#0047A099')
        v.get_patch_by_id('10').set_edgecolor('#0047A0')
        v.get_patch_by_id('01').set_color('#CD2E3A99')
        v.get_patch_by_id('01').set_edgecolor('#CD2E3A')
        v.get_patch_by_id('11').set_color('#00000099')
        v.get_patch_by_id('11').set_edgecolor('#000000')

        for label_id, val in [('10', kpanref_only), ('01', cpchprc_only), ('11', shared)]:
            lbl = v.get_label_by_id(label_id)
            if lbl:
                lbl.set_text(fmt_val(val))
                lbl.set_fontsize(12)
                lbl.set_fontweight('bold')

        for lbl in v.set_labels:
            if lbl:
                lbl.set_fontsize(12)
                lbl.set_fontweight('bold')

        ax.set_title(vtitle, fontsize=14, fontweight='bold', color='black', pad=15)

    draw_venn(ax_venn_small, "Small Variant", "Small Variant")
    draw_venn(ax_venn_sv,    "SV",            "SV")

    fig.canvas.draw()
    label_axes = [axes[0, 0], axes[1, 0], ax_venn_small]
    labels     = ['A', 'B', 'C']
    for ax, label in zip(label_axes, labels):
        bbox = ax.get_position()
        fig.text(bbox.x0 - 0.02, bbox.y1 + 0.01, label, fontsize=16, fontweight='bold',
                 color='black', va='bottom', ha='left', transform=fig.transFigure)

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    per_sample_path = $PATH
    annot_paths = {
        "KPanRef": $PATH,
        "CPC":     $PATH,
        "HPRC":    $PATH,
    }
    venn_path   = $PATH
    output_path = $PATH
    plot_variant_stats(per_sample_path, annot_paths, venn_path, output_path)
