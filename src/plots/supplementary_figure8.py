#%%
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

#%%
def plot_sv_length(annot_paths, output_path):

    colors = {
        "Control": "#3D7A52",
        "MI":      "#C9A227",
    }

    annot_dfs = {}
    for graph, path in annot_paths.items():
        df_annot = pd.read_csv(path, sep="\t", low_memory=False)
        df_annot = df_annot.drop_duplicates(subset="AnnotSV_ID", keep="first").reset_index(drop=True)
        annot_dfs[graph] = df_annot

    fig = plt.figure(figsize=(15, 10))
    fig.patch.set_facecolor('white')

    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.4)

    for idx, graph in enumerate(["Control", "MI"]):
        ax = fig.add_subplot(gs[idx, 0:3])
        ax.set_facecolor('white')

        color = colors.get(graph, "#555555")
        sv_len_orig = annot_dfs[graph]["SV_length"].dropna()
        sv_len_orig = sv_len_orig[sv_len_orig != 0].copy()
        sv_len = sv_len_orig.clip(lower=-10000, upper=10000)

        bins = np.arange(-10000, 10100, 100)
        ax.hist(sv_len, bins=bins, color=color + '99', edgecolor=color, linewidth=0.8)

        ax.axvline(x=0, color='black', linewidth=1.0, linestyle='--')

        ax.set_title(f"{graph}", fontsize=14, fontweight='bold', color='black', pad=15)
        ax.set_xlabel("kb", fontsize=13, fontweight='bold')
        ax.set_ylabel("Count (k)", fontsize=13, fontweight='bold')
        ax.set_xticks(np.arange(-10000, 10100, 2000))
        ax.set_xticklabels(
            [str(int(b/1e3)) for b in np.arange(-10000, 10100, 2000)]
        )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_color('black')
        ax.spines['bottom'].set_color('black')
        ax.tick_params(colors='black', labelsize=12)
        ax.yaxis.grid(False)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x/1e3)}'))

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    annot_paths = {
        "Control": $PATH,
        "MI":      $PATH,
    }
    output_path = $PATH
    plot_sv_length(annot_paths, output_path)
