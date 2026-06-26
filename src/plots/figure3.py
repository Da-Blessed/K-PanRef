#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
    order = ["Control", "MI"]

    df_samp = pd.read_csv(per_sample_path, sep="\t")
    df_samp["Graph Name"] = df_samp["Graph Name"].replace("control", "Control")
    df_samp["SNV"] = df_samp["Transitions"] + df_samp["Transversions"]
    df_samp["InDel"] = df_samp["Indels"]
    df_samp["Graph Name"] = pd.Categorical(df_samp["Graph Name"], categories=order, ordered=True)
    df_samp = df_samp.sort_values("Graph Name").reset_index(drop=True)

    annot_dfs = {}
    for graph, path in annot_paths.items():
        df_annot = pd.read_csv(path, sep="\t", low_memory=False)
        df_annot = df_annot.drop_duplicates(subset="AnnotSV_ID", keep="first").reset_index(drop=True)
        df_annot["SV_length"] = df_annot["SV_length"].abs()
        annot_dfs[graph] = df_annot

    df_venn = pd.read_csv(venn_path, sep="\t")
    df_venn.columns = df_venn.columns.str.strip()

    colors = {
        "Control": "#3D7A52",
        "MI":      "#C9A227",
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

    fig = plt.figure(figsize=(22, 12))
    fig.patch.set_facecolor('white')

    gs_top = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35,
                               top=0.95, bottom=0.54)
    gs_bot = gridspec.GridSpec(1, 3, figure=fig, wspace=0.3,
                               top=0.40, bottom=0.05)

    axes = np.empty((1, 3), dtype=object)
    for j in range(3):
        axes[0, j] = fig.add_subplot(gs_top[0, j])

    ax_sv_hist    = fig.add_subplot(gs_bot[0, 0])
    ax_venn_small = fig.add_subplot(gs_bot[0, 1])
    ax_venn_sv    = fig.add_subplot(gs_bot[0, 2])

    pos = ax_sv_hist.get_position()
    shrink = 0.86
    new_height = pos.height * shrink
    new_y0 = pos.y0 + (pos.height - new_height) / 2
    ax_sv_hist.set_position([pos.x0, new_y0, pos.width, new_height])

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

        margin = (data_max - data_min) * 6.0
        ax.set_ylim(data_min - margin * 0.5, data_max + margin * 0.5)

        ax.set_title(title, fontsize=14, fontweight='bold', color='black', pad=15)
        ax.set_xticks(positions)
        ax.set_xticklabels(graph_list, fontsize=11, fontweight='bold')
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

        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))

    ax_sv_hist.set_facecolor('white')
    alphas = {"Control": 0.5, "MI": 0.7}
    sv_lens = {}
    for graph in ["Control", "MI"]:
        color = colors.get(graph, "#555555")
        sv_len_orig = annot_dfs[graph]["SV_length"].dropna()
        sv_len_orig = sv_len_orig[sv_len_orig > 0].copy()
        sv_len = sv_len_orig.clip(upper=10000)
        sv_lens[graph] = sv_len

        bins = np.arange(0, 10100, 100)
        ax_sv_hist.hist(sv_len, bins=bins, color=color + '99', edgecolor=color,
                        linewidth=0.8, alpha=alphas[graph])

    legend_elements = [
        mpatches.Patch(facecolor=colors["Control"] + '99', edgecolor=colors["Control"],
                       alpha=0.5, label='Control'),
        mpatches.Patch(facecolor=colors["MI"] + '99', edgecolor=colors["MI"],
                       alpha=0.7, label='MI'),
    ]
    ax_sv_hist.legend(handles=legend_elements, fontsize=10, frameon=True, loc='upper right')

    ax_sv_hist.set_title("SV Length Distribution", fontsize=14, fontweight='bold', color='black', pad=15)
    ax_sv_hist.set_xlabel("kb", fontsize=12, fontweight='bold')
    ax_sv_hist.set_ylabel("Count (k)", fontsize=12, fontweight='bold')
    ax_sv_hist.set_xticks(np.arange(0, 10100, 1000))
    ax_sv_hist.set_xticklabels([str(int(b/1e3)) for b in np.arange(0, 10100, 1000)])
    ax_sv_hist.spines['top'].set_visible(False)
    ax_sv_hist.spines['right'].set_visible(False)
    ax_sv_hist.spines['left'].set_color('black')
    ax_sv_hist.spines['bottom'].set_color('black')
    ax_sv_hist.tick_params(colors='black', labelsize=12)
    ax_sv_hist.yaxis.grid(False)
    ax_sv_hist.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1e3)}'))

    ax_inset = inset_axes(ax_sv_hist, width="48%", height="48%",
                          loc='center',
                          bbox_to_anchor=(-0.02, -0.05, 1, 1),
                          bbox_transform=ax_sv_hist.transAxes)
    ax_inset.set_facecolor('white')

    bins_inset = np.arange(0, 1050, 50)
    for graph in ["Control", "MI"]:
        color = colors.get(graph, "#555555")
        sv_zoom = sv_lens[graph][sv_lens[graph] <= 1000]
        ax_inset.hist(sv_zoom, bins=bins_inset, color=color + '99', edgecolor=color,
                      linewidth=0.6, alpha=alphas[graph])

    ax_inset.set_title("0–1 kb", fontsize=8, fontweight='bold', color='black', pad=4)
    ax_inset.set_xlabel("bp", fontsize=7, fontweight='bold')
    ax_inset.set_xticks(np.arange(0, 1100, 200))
    ax_inset.set_xticklabels([str(int(b)) for b in np.arange(0, 1100, 200)], fontsize=6)
    ax_inset.tick_params(colors='black', labelsize=6)
    ax_inset.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1e3)}'))
    ax_inset.set_ylabel("Count (k)", fontsize=7, fontweight='bold')
    ax_inset.spines['top'].set_visible(False)
    ax_inset.spines['right'].set_visible(False)
    ax_inset.spines['left'].set_color('black')
    ax_inset.spines['bottom'].set_color('black')
    ax_inset.yaxis.grid(False)

    def draw_venn(ax, vtype, vtitle):
        ax.set_facecolor('white')
        row = df_venn[df_venn["Type"] == vtype].iloc[0]
        control_only = int(row["Control"])
        shared       = int(row["Shared"])
        mi_only      = int(row["MI"])

        v = venn2(
            subsets=(1, 1, 1),
            set_labels=("Control", "MI"),
            ax=ax
        )

        v.get_patch_by_id('10').set_color(colors["Control"] + '99')
        v.get_patch_by_id('10').set_edgecolor(colors["Control"])
        v.get_patch_by_id('01').set_color(colors["MI"] + '99')
        v.get_patch_by_id('01').set_edgecolor(colors["MI"])
        v.get_patch_by_id('11').set_color('#9ACD3299')
        v.get_patch_by_id('11').set_edgecolor('#9ACD32')

        for label_id, val in [('10', control_only), ('01', mi_only), ('11', shared)]:
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
    label_axes = [axes[0, 0], ax_sv_hist, ax_venn_small]
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
        "Control": $PATH,
        "MI":      $PATH,
    }
    venn_path   = $PATH
    output_path = $PATH
    plot_variant_stats(per_sample_path, annot_paths, venn_path, output_path)
