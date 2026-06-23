#%%
from collections import Counter
import gzip
import numpy as np
import subprocess
import glob
from io import SEEK_END
import io
import os
import struct
import zlib
import pickle
from joblib import Parallel, delayed
import argparse

#%%
def estimate_uncompressed_gz_size(filename, read_bytes = 10000000):
    with open(filename, "rb") as gz_in:
        gz_in.seek(-4, SEEK_END)
        lsb = struct.unpack('<I', gz_in.read(4))[0]
        file_size = os.fstat(gz_in.fileno()).st_size
    
    with gzip.open(filename, "rb") as frb:
        sample = frb.read(read_bytes)
    
    compressed_sample = io.BytesIO()
    with gzip.GzipFile(fileobj = compressed_sample, mode = "wb") as fwb:
        fwb.write(sample)
    
    compressed_len = len(compressed_sample.getvalue())
    decompressed_len = len(sample)

    estimate = int(file_size * decompressed_len / compressed_len)

    mask = ~0xFFFFFFFF

    adjusted_estimate = (estimate & mask) | lsb
    return adjusted_estimate

#%%
def get_kmer_count(sequence, nmer):
    lenseq = len(sequence)
    
    list_kmer = list(map(lambda ind: sequence[ind:ind+nmer], range(lenseq-nmer+1)))
    return Counter(list_kmer)

#%%
def get_kmer_count_of_sequence_with_chop(sequence, chop_length, nmer = 5):
    sequence_front = sequence[:chop_length]
    sequence_end = sequence[-chop_length:]
    sequence_intermediate = sequence[chop_length:-chop_length]
    
    count_kmer_front        = get_kmer_count(sequence_front, nmer)
    count_kmer_end          = get_kmer_count(sequence_end, nmer)
    count_kmer_intermediate = get_kmer_count(sequence_intermediate, nmer)
    
    return (count_kmer_front, count_kmer_end, count_kmer_intermediate)

#%%
class FILE_READER:
    def __init__(self, file_path, is_gzip):
        self.is_gzip = is_gzip
        if is_gzip:
            self.file_reader = gzip.open(file_path, "rb")
        else:
            self.file_reader = open(file_path, 'r')
    
    def readline(self):
        if self.is_gzip:
            return self.file_reader.readline().decode()
        else:
            return self.file_reader.readline()
    
    def seek(self, val):
        self.file_reader.seek(val)
    
    def tell(self):
        return self.file_reader.tell()

    def close(self):
        self.file_reader.close()
    
    def __del__(self):
        self.close()

#%%
def calculate_mean_error_rate_as_phred_scale(arr_phredscaled_quality):
    arr_error_rate = 10 ** (-arr_phredscaled_quality/10)
    mean_error_rate = np.mean(arr_error_rate)
    mean_error_rate_phredscaled = -10 * np.log10(mean_error_rate)
    
    return mean_error_rate_phredscaled

#%%
def get_stat_from_file_range(file_path, start, end, is_gzip, kmer_chop_length, n_kmer, index_low_quality, count_low_quality, readid_prefix = None, readid_suffix = None):
    dict_stat = {
        "mean_base_quality":list(),
        "QV":list(),
        "std_base_quality":list(),
        "length":list(),
        "num_N_base":list(),
        "num_Q30_base":list(),
        "num_LowQ_base":list(),
        "idx_lowQ_base":list(),
        "kmer_front":Counter(),
        "kmer_end":Counter(),
        "kmer_intermediate":Counter(),
        "readID":list()
    }
    
    filereader = FILE_READER(file_path, is_gzip)
    filereader.seek(start)
    set_atgc = set("ATGCN")
    while 1:
        name = filereader.readline()
        if name == '':
            break
        if name.startswith('@'):
            name = name.strip("@\n")
            if readid_prefix != None:
                if not name.startswith(readid_prefix):
                    continue
            if readid_suffix != None:
                if not name.endswith(readid_suffix):
                    continue
            # if len(name) != 36:
            #     continue
        else:
            continue
        sequence = filereader.readline().strip('\n').upper()
        _ = filereader.readline()
        quality = filereader.readline().strip('\n')
        
        assert set(sequence) <= set_atgc, name
        
        list_quality_int = list(map(lambda val: ord(val)-33, quality))
        n_q30_base = sum(list(map(lambda val: val >= 30, list_quality_int)))
        
        dict_stat["mean_base_quality"].append(np.mean(list_quality_int))
        dict_stat["QV"].append(calculate_mean_error_rate_as_phred_scale(np.array(list_quality_int)))
        dict_stat["std_base_quality"].append(np.std(list_quality_int))
        dict_stat["length"].append(len(quality))
        dict_stat["num_N_base"].append(Counter(sequence)['N'])
        dict_stat["num_Q30_base"].append(n_q30_base)
        dict_stat["readID"].append(name)
        
        if count_low_quality != None:
            dict_stat["num_LowQ_base"].append(sum(list(map(lambda val: val < count_low_quality, list_quality_int))))
        
        if index_low_quality != None:
            dict_stat["idx_lowQ_base"].append(
                tuple(filter(lambda ind: list_quality_int[ind] < index_low_quality, range(len(list_quality_int))))
            )
        
        if n_kmer != None:
            count_kmer_front, count_kmer_end, count_kmer_intermediate = get_kmer_count_of_sequence_with_chop(sequence, kmer_chop_length, n_kmer)
            for kmer_seq, cnt in count_kmer_front.items():
                dict_stat["kmer_front"][kmer_seq] += cnt
            for kmer_seq, cnt in count_kmer_end.items():
                dict_stat["kmer_end"][kmer_seq] += cnt
            for kmer_seq, cnt in count_kmer_intermediate.items():
                dict_stat["kmer_intermediate"][kmer_seq] += cnt
        
        current_cursor = filereader.tell()
        if end > 0 and current_cursor > end:
            break
    filereader.close()
    return dict_stat

#%%
def get_stat_from_file_with_ncores(file_path, n_cpu, kmer_chop_length, n_kmer, index_low_quality = None, count_low_quality = None, filesize = None, is_gzip = True, readid_prefix = None, readid_suffix = None): 
    if filesize == None:
        if is_gzip:
            filesize = estimate_uncompressed_gz_size(file_path)
        else:
            filesize = os.stat(file_path).st_size
    
    list_ranges = list(map(int, np.linspace(0, filesize, n_cpu+1)))
    list_ranges[-1] = -1
    
    list_starts = list_ranges[:-1]
    list_ends = list_ranges[1:]
    
    with Parallel(n_jobs=n_cpu) as parallel:
        list_stats = parallel(delayed(get_stat_from_file_range)(
            file_path, start, end, is_gzip, kmer_chop_length, n_kmer, index_low_quality, count_low_quality, readid_prefix, readid_suffix
        )for start, end in zip(list_starts, list_ends))
    
    dict_stat = {
        "mean_base_quality":list(),
        "QV":list(),
        "std_base_quality":list(),
        "length":list(),
        "num_N_base":list(),
        "num_Q30_base":list(),
        "num_LowQ_base":list(),
        "idx_lowQ_base":list(),
        "kmer_front":Counter(),
        "kmer_end":Counter(),
        "kmer_intermediate":Counter(),
        "readID":list()
    }
    
    for stat_tmp in list_stats:
        for key in dict_stat.keys():
            if isinstance(dict_stat[key], list):
                dict_stat[key].extend(stat_tmp[key])
            if isinstance(dict_stat[key], Counter):
                for kmer_seq, cnt in stat_tmp[key].items():
                    dict_stat[key][kmer_seq] += cnt
    return dict_stat

#%%
if __name__=="__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--input")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--n_cpu", type=int)
    argument_parser.add_argument("--kmer_chop_length", type=int)
    argument_parser.add_argument("--n_kmer", type=int)
    argument_parser.add_argument("--index_low_quality", type=float)
    argument_parser.add_argument("--count_low_quality", type=float)
    argument_parser.add_argument("--filesize", type=float)
    argument_parser.add_argument("--is_gzip", action="store_true")
    argument_parser.add_argument("--readid_prefix")
    argument_parser.add_argument("--readid_suffix")

    args = argument_parser.parse_args()
    
    dict_stat = get_stat_from_file_with_ncores(file_path=args.input, n_cpu=args.n_cpu, kmer_chop_length=args.kmer_chop_length, n_kmer=args.n_kmer, index_low_quality=args.index_low_quality, count_low_quality=args.count_low_quality, filesize=args.filesize, is_gzip=args.is_gzip, readid_prefix=args.readid_prefix, readid_suffix=args.readid_suffix)
    with open(args.out, "wb") as fwb:
        pickle.dump(dict_stat, fwb)    
