#%%
import glob
import pickle
from itertools import product
import pandas as pd
import math
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import argparse

#%%
def pickle_to_dict(pickle_path):
    with open(pickle_path, "rb") as frb:
        dict_stat = pickle.load(frb)
    return dict_stat

#%%
def make_kmer_count_dataframe(kmer_count, n_kmer, contexts = "ATCG"):
    n_mer_index = math.ceil(n_kmer/2)
    n_mer_column = n_kmer - n_mer_index
    
    list_context_index = list(map(lambda vals: ''.join(vals), product(contexts, repeat = n_mer_index)))
    list_context_column = list(map(lambda vals: ''.join(vals), product(contexts, repeat = n_mer_column)))
    
    table_kmer_count = pd.DataFrame(index = list_context_index, columns = list_context_column).fillna(0).astype(int)
    
    for kmer_seq, cnt in kmer_count.items():
        context_ind = kmer_seq[:n_mer_index]
        context_col = kmer_seq[-n_mer_column:]
        
        table_kmer_count.loc[context_ind, context_col] = cnt
    
    return table_kmer_count

#%%
def plot_kmer_heatmaps(dict_stat, pickle_path, n_kmer=5, cmap="flare_r"):
    plt.rcParams["font.size"] = 10
    kmer_types = ["front", "end", "intermediate"]

    for kmer_type in kmer_types:
        table = make_kmer_count_dataframe(dict_stat[f"kmer_{kmer_type}"], n_kmer=n_kmer)

        plt.figure(figsize=(15, 13))
        sns.heatmap(data=table, cmap=cmap)
        plt.title(f"Read {kmer_type} k-mer composition", fontsize=12)
        plt.xlabel("Suffix Sequence")
        plt.ylabel("Prefix Sequence")
        plt.tight_layout()
        plt.savefig(f"{('/').join(pickle_path.split('/')[0:-1])}/{pickle_path.split('/')[-1].split('.')[0]}_read_{kmer_type}_k-mer_composition.png", format="png", dpi=300)
        plt.close()

#%%
def get_nvalue(list_length, n = 50):
    total_length = sum(list_length)
    
    val_thresh = total_length * n / 100
    
    sum_val = 0
    
    length_on_thresh = None
    n_sequence_for_thresh = None
    
    for ind, val in enumerate(sorted(list_length, reverse = True), start = 1):
        sum_val += val
        if sum_val >= val_thresh:
            length_on_thresh = val
            n_sequence_for_thresh = ind
            break
    return length_on_thresh, n_sequence_for_thresh
            
#%%    
def write_stat_tsv(dict_stat, pickle_path):
    n50, l50 = get_nvalue(dict_stat["length"])
    n90, l90 = get_nvalue(dict_stat["length"], 90)

    stats = {
        "Total_Reads": len(dict_stat["readID"]),
        "Total_Bases": sum(dict_stat["length"]),
        "Total_N_Bases": sum(dict_stat["num_N_base"]),
        
        "Max_Read_Length": max(dict_stat["length"]),
        "Mean_Read_Length": np.mean(dict_stat["length"]),
        "Min_Read_Length": min(dict_stat["length"]),
        "Read_Length_1_Percentile": np.quantile(dict_stat["length"], 0.01),
        "Read_Length_99_Percentile": np.quantile(dict_stat["length"], 0.99),

        
        "Read_N50": n50,
        "Read_L50": l50,
        "Read_N90": n90,
        "Read_L90": l90,
        
        "Q30_Base_Ratio(%)": sum(dict_stat["num_Q30_base"]) / sum(dict_stat["length"]) * 100,
        "Mean_Read_QV": np.mean(dict_stat["QV"]),
        "Mean_Base_Quality": sum(np.array(dict_stat["mean_base_quality"]) * np.array(dict_stat["length"])) / sum(dict_stat["length"]),
        "QV": -10 * np.log10(
            sum((10 ** (-np.array(dict_stat["QV"]) / 10)) * np.array(dict_stat["length"])) / sum(dict_stat["length"])
        ),
    }
    
    save_path = f"{('/').join(pickle_path.split('/')[0:-1])}/{pickle_path.split('/')[-1].split('.')[0]}_read_stats.tsv"
    pd.DataFrame(list(stats.items())).to_csv(save_path, sep="\t", index=False, header=False)

#%%
def plot_quality_scatterplot(dict_stat, pickle_path):
    table_stat = pd.DataFrame(
        {
            key : dict_stat[key] for key in ['mean_base_quality', 'std_base_quality', 'QV', 'length', 'num_N_base', 'num_Q30_base', 'num_LowQ_base', 'readID']
        }    
    )
    table_stat["perc_Q30_base"] = table_stat["num_Q30_base"] / table_stat["length"] * 100
    table_stat["perc_lowQ_base"] = table_stat["num_LowQ_base"] / table_stat["length"] * 100
    
    dict_key_to_show = {
        "mean_base_quality" : "Mean Base Phred Quality Score",
        "QV" : "Quality Value per Read",
        # "std_base_quality" : "STD Base Quality",
        "perc_Q30_base" : "Q30 Base (%)",
        # "perc_lowQ_base" : "<Q10 Base (%)",
        "length" : "Read Length (bp)"
    }
    
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(
        len(dict_key_to_show), len(dict_key_to_show), figsize = (4 * len(dict_key_to_show), 4 * len(dict_key_to_show))
    )
    plt.subplots_adjust(
        wspace = 0.25, hspace = 0.25
    )
    
    dict_item_to_lim = dict()
    for item in dict_key_to_show.keys():
        valmin = min(table_stat[item])
        valmax = max(table_stat[item])
        valwidth = valmax - valmin
        valadjust = valwidth / 20
        
        dict_item_to_lim[item] = (valmin - valadjust, valmax + valadjust)
    
    for ind_item1, item1 in enumerate(dict_key_to_show.keys()):
        for ind_item2, item2 in enumerate(dict_key_to_show.keys()):
            ax_here = ax[ind_item1, ind_item2]
            if ind_item2 == 0:
                ax_here.set_ylabel(dict_key_to_show[item1])
            if ind_item1 == len(dict_key_to_show)-1:
                ax_here.set_xlabel(dict_key_to_show[item2])
            
            if ind_item1 < ind_item2:
                ax_here.axis("off")
                continue
            if ind_item1 == ind_item2:
                lowerbound, upperbound = dict_item_to_lim[item1]
                ax_here.hist(
                    table_stat[item1],
                    color = "gray",
                    bins = np.linspace(lowerbound, upperbound, 67)
                )    
                ax_here.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                ax_here.set_xlim(lowerbound, upperbound)
                # ax_here.set_ylabel("# Reads")
            else:
                xlim_lowerbound, xlim_upperbound = dict_item_to_lim[item2]
                ylim_lowerbound, ylim_upperbound = dict_item_to_lim[item1]
                
                ax_here.scatter(
                    table_stat[item2],
                    table_stat[item1],
                    s = 5,
                    color = "gray"
                )
                ax_here.set_xlim(xlim_lowerbound, xlim_upperbound)
                ax_here.set_ylim(ylim_lowerbound, ylim_upperbound)
                if item1 == "length":
                    ax_here.set_yticklabels([f"{ytickval/1000:.1f}k" for ytickval in ax_here.get_yticks()])
            # if item2 == "qv_by_read":
            #     ax_here.axvline(20, linewidth = 2, color = "firebrick", linestyle = '--')
            # elif item1 == "qv_by_read":
            #     ax_here.axhline(20, linewidth = 2, color = "firebrick", linestyle = '--')
    plt.savefig(f"{('/').join(pickle_path.split('/')[0:-1])}/{pickle_path.split('/')[-1].split('.')[0]}_Read_stat_plot.png", format="png", dpi=300)
    plt.close()
    
# %%
def get_lowquality_base_count_from_end(dict_stat, check_length):
    cnt_front = np.zeros(check_length)
    cnt_end = np.zeros(check_length)
    
    checkarange = np.arange(check_length)
    
    for length, list_lowqual_ind in zip(dict_stat["length"], dict_stat["idx_lowQ_base"]):
        arr_lowqual_ind_front = np.array(list_lowqual_ind)
        arr_lowqual_ind_end = (length - 1) - arr_lowqual_ind_front
        
        cnt_front += np.isin(checkarange, arr_lowqual_ind_front[:check_length])
        cnt_end += np.isin(checkarange, arr_lowqual_ind_end[-check_length:])
        
    return cnt_front, cnt_end
    

# %%
def draw_lowquality_perc_on_end(dict_stat, pickle_path, check_length):
    cnt_front, cnt_end = get_lowquality_base_count_from_end(dict_stat=dict_stat, check_length=check_length)
    n_read = len(dict_stat["readID"])
    
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))

    perc_front = cnt_front / n_read * 100
    perc_end = cnt_end / n_read * 100

    ymin = min(perc_front.min(), perc_end.min()) * 0.9
    ymax = max(perc_front.max(), perc_end.max()) * 1.1  

    ax[0].plot(range(check_length), perc_front, color="gray")
    ax[0].set_xlabel("Coordinate (bp)")
    ax[0].set_ylabel("Percentage of Low Quality Base (%)")
    ax[0].set_title("Read Front")
    ax[0].set_ylim(ymin, ymax)

    ax[1].plot(range(check_length), perc_end, color="gray")
    ax[1].set_xlabel("Coordinate (bp)")
    ax[1].set_ylabel("Percentage of Low Quality Base (%)")
    ax[1].set_title("Read End")
    ax[1].set_ylim(ymin, ymax)

    plt.tight_layout()
    plt.savefig(
        f"{('/').join(pickle_path.split('/')[0:-1])}/{pickle_path.split('/')[-1].split('.')[0]}_Read_front_end_{check_length}bp_lowQ_base_percent.png",
        format="png", dpi=300
    )
    plt.close()


#%%
def run_plot_read_stats(pickle_path, check_length=[500], n_kmer=None):
    dict_stat = pickle_to_dict(pickle_path)
    if n_kmer != None:
        plot_kmer_heatmaps(dict_stat=dict_stat, pickle_path=pickle_path, n_kmer=n_kmer, cmap="flare_r")
    write_stat_tsv(dict_stat=dict_stat, pickle_path=pickle_path)
    plot_quality_scatterplot(dict_stat=dict_stat, pickle_path=pickle_path)
    for c_length in check_length:
        draw_lowquality_perc_on_end(dict_stat=dict_stat, pickle_path=pickle_path, check_length=c_length)

# %%
if __name__=="__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--input")
    argument_parser.add_argument("--check_length", type=str)
    argument_parser.add_argument("--n_kmer", type=int)
    args = argument_parser.parse_args()
    
    check_length = [int(x) for x in args.check_length.split(",")]
    run_plot_read_stats(pickle_path=args.input, check_length=check_length, n_kmer=args.n_kmer)
