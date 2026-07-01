#%%
import pandas as pd
import matplotlib.pyplot as plt
import re

#%%
def style_boxplot(bp, color):
    for element in ['boxes', 'whiskers', 'caps', 'medians']:
        plt.setp(bp[element], color=color, linewidth=1.5)
    plt.setp(bp['fliers'], marker='o', markerfacecolor=color,
             markeredgecolor=color, markersize=4, alpha=0.5)
    for patch in bp['boxes']:
        patch.set_facecolor(color + '99')

#%%
def plot_ax(ax, data, color, label, col, row_title=None):
    clean = data.dropna().copy()

    if 'rate' in col or 'gc_content' in col:
        clean = clean * 100

    if 'rate' in col or 'gc_content' in col or '(%)' in label:
        unit    = "%"
        display = clean
        ylabel  = "%"
    elif 'QV' in col:
        unit    = ""
        display = clean
        ylabel  = "Value"
    elif 'Reads' in col or 'reads' in col:
        if clean.max() >= 1e9:
            display = clean / 1e9
            unit = "G"
        elif clean.max() >= 1e6:
            display = clean / 1e6
            unit = "M"
        elif clean.max() >= 1e3:
            display = clean / 1e3
            unit = "k"
        else:
            display = clean
            unit = ""
        ylabel = f"Count ({unit})" if unit else "Count"
    else:
        if clean.max() >= 1e9:
            display = clean / 1e9
            unit = "G"
        elif clean.max() >= 1e6:
            display = clean / 1e6
            unit = "M"
        elif clean.max() >= 1e3:
            display = clean / 1e3
            unit = "k"
        else:
            display = clean
            unit = ""
        m = re.search(r'\((.+?)\)', label)
        if m:
            ylabel = m.group(1)
        else:
            ylabel = unit if unit else ""

    title_label = re.sub(r'\s*\(.*?\)', '', label).strip()

    bp = ax.boxplot(display, widths=0.5, patch_artist=True)
    style_boxplot(bp, color)
    ax.set_title(title_label, fontsize=14, fontweight='bold', color='black', pad=15)
    ax.set_ylabel(ylabel, fontsize=13, fontweight='bold')
    ax.set_xticks([])
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(colors='black', labelsize=12)
    ax.yaxis.grid(False)

    if 'QV' in col:
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))
    elif 'rate' in col or 'gc_content' in col:
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))

    if row_title:
        ax.text(0.5, 1.32, row_title, transform=ax.transAxes,
                ha='center', fontsize=16, fontweight='bold', color='black')

#%%
def plot_combined(pangenie_path, hifiasm_path, output_path):
    df_p = pd.read_csv(pangenie_path, sep="\t")
    df_h = pd.read_csv(hifiasm_path, sep="\t")
    df_p = df_p[df_p["Sample Name"] != "KU10K-02029"].reset_index(drop=True)
    df_h = df_h[df_h["Sample Name"] != "KU10K-02029"].reset_index(drop=True)

    df_p["WGS_total_bases_Gb"] = df_p["WGS_total_bases"] / 1e9
    df_h["ONT_Total_Bases_Gb"]  = df_h["ONT_Total_Bases"] / 1e9
    df_h["ONT_Read_N50_kb"]     = df_h["ONT_Read_N50"] / 1e3
    df_h["HiFi_Total_Bases_Gb"] = df_h["HiFi_Total_Bases"] / 1e9
    df_h["HiFi_Read_N50_kb"]    = df_h["HiFi_Read_N50"] / 1e3
    df_h["HiC_total_bases_Gb"]  = df_h["Hi-C_total_bases"] / 1e9

    rows = [
        (df_h, ["HiFi_Total_Bases_Gb", "HiFi_Read_N50_kb", "HiFi_QV"],
               ["Total Bases (Gb)", "Read N50 (kb)", "QV"],
               "#000000", "PacBio HiFi"),
        (df_h, ["ONT_Total_Bases_Gb", "ONT_Read_N50_kb", "ONT_QV"],
               ["Total Bases (Gb)", "Read N50 (kb)", "QV"],
               "#CD2E3A", "ONT-UL"),
        (df_h, ["HiC_total_bases_Gb", "Hi-C_q30_rate", "Hi-C_gc_content"],
               ["Total Bases (Gb)", "Q30 Rate (%)", "GC Content (%)"],
               "#0047A0", "Hi-C"),
        (df_p, ["WGS_total_bases_Gb", "WGS_q30_rate", "WGS_gc_content"],
               ["Total Bases (Gb)", "Q30 Rate (%)", "GC Content (%)"],
               "#0047A0", "WGS"),
    ]

    fig, axes = plt.subplots(4, 3, figsize=(15, 20))
    fig.patch.set_facecolor('white')
    plt.subplots_adjust(wspace=0.4, hspace=0.7)

    for row_idx, (df, cols, labels, color, title) in enumerate(rows):
        for i, (col, label) in enumerate(zip(cols, labels)):
            row_title = title if i == 1 else None
            plot_ax(axes[row_idx, i], df[col], color, label, col, row_title=row_title)

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    plot_combined(
        pangenie_path = $PATH,
        hifiasm_path  = $PATH,
        output_path   = $PATH
    )
