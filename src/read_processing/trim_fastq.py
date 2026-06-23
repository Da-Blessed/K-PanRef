#%%
from collections import Counter
import gzip, subprocess
import numpy as np
from zlib_ng import gzip_ng
import argparse
import os
  
#%%
class FASTQ_READER:
    def __init__(self, path_fasta, is_gzip):
        self.is_gzip = is_gzip
        if is_gzip:
            self.file_reader = gzip_ng.open(path_fasta, "rb")
        else:
            self.file_reader = open(path_fasta, 'r')
    
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
class FASTQ_WRITER:
    def __init__(self, path_fasta, is_gzip):
        self.is_gzip = is_gzip
        if is_gzip:
            self.file_reader = gzip_ng.open(path_fasta, "wb")
        else:
            self.file_reader = open(path_fasta, 'w')
    
    def write(self, value:str):
        if self.is_gzip:
            return self.file_reader.write(value.encode())
        else:
            return self.file_reader.write(value)
        
    def writeline(self, value:str):
        value_end = f"{value}\n"
        self.write(value_end)

    def close(self):
        self.file_reader.close()
    
    def __del__(self):
        self.close()

#%%
def trim_fastq_front_and_tail(path_fastq, path_save, is_gzip, trim_front, trim_tail, min_length):
    import numpy as np
    
    filereader = FASTQ_READER(path_fastq, is_gzip)
    filewriter = FASTQ_WRITER(path_save, is_gzip)
    
    set_atgc = set("ATGCN")
    
    n_passed = 0
    n_filtered_by_length = 0
    
    while 1:
        namevalue = filereader.readline()
        if namevalue == '':
            break
        if namevalue.startswith('@'):
            name = namevalue.strip("@\n")
        else:
            continue
        sequencevalue = filereader.readline()
        tmp = filereader.readline()
        qualityvalue = filereader.readline()
        
        sequence = sequencevalue.strip('\n').upper()
        quality = qualityvalue.strip('\n')
        
        assert set(sequence) <= set_atgc, name
        
        len_read = len(sequence)
        
        write_line = True
        if len_read < (trim_front + trim_tail + min_length):
            write_line = False
            n_filtered_by_length += 1
        if write_line:
            n_passed += 1
            filewriter.write(namevalue)
            filewriter.write(sequence[trim_front:-trim_tail]+'\n')
            filewriter.write(tmp)
            filewriter.write(quality[trim_front:-trim_tail]+'\n')
        
        n_processed = n_passed+n_filtered_by_length
        if (n_processed) % 10000 == 0:
            print(n_processed)
        
    filereader.close()
    filewriter.close()
    return (n_passed, n_filtered_by_length)

#%%
def run(path_fastq, path_save, is_gzip, trim_front, trim_tail, min_length):    
    n_passed, n_filtered_by_length = trim_fastq_front_and_tail(
        path_fastq, path_save, is_gzip, trim_front, trim_tail, min_length
    )
    
    with open(f"{path_save}.stat", 'w') as fw:
        fw.write(f"#Total : {n_passed+n_filtered_by_length}\n")
        fw.write(f"#Passed : {n_passed}\n")
        fw.write(f"#Filtered_by_Length : {n_filtered_by_length}\n")

    print("Done:", path_fastq, sep = ' ')
    
#%%
if __name__ == "__main__":    
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--input", type=str)
    argument_parser.add_argument("--out", type=str)
    argument_parser.add_argument("--is_gzip", action="store_true")
    argument_parser.add_argument("--trim_front", type=int)
    argument_parser.add_argument("--trim_end", type=int)
    argument_parser.add_argument("--min_length", type=int)
    
    args = argument_parser.parse_args()
    
    run(args.input, args.out, args.is_gzip, args.trim_front, args.trim_end, args.min_length)
