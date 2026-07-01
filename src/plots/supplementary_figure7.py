#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#%%
def plot_go_enrichment(go_path, output_path):
    df = pd.read_csv(go_path)
    df = df.head(20).copy()
    df["-log10(FDR)"] = -np.log10(df["Enrichment FDR"])
    df["Pathway"] = df["Pathway"].str.replace(r"GO:\d+ ", "", regex=True).str.strip()
    df = df.sort_values("-log10(FDR)", ascending=True)

    fig, ax = plt.subplots(figsize=(15, 15))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    y_pos = np.arange(len(df))
    ax.barh(y_pos, df["-log10(FDR)"], color="#0047A099", edgecolor="#0047A0", linewidth=1.0)
    ax.axvline(x=-np.log10(0.05), color='#CD2E3A', linewidth=1.5, linestyle='--')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["Pathway"], fontsize=13)
    ax.set_xlabel("-log10(FDR)", fontsize=15, fontweight='bold')
    ax.set_title("GO Enrichment (Top 20)", fontsize=18, fontweight='bold', color='black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.tick_params(colors='black', labelsize=12)

    plt.savefig(output_path, bbox_inches='tight', format='svg')
    plt.show()

#%%
if __name__ == "__main__":
    go_path     = $PATH
    output_path = $PATH
    plot_go_enrichment(go_path, output_path)
