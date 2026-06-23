#%%
import gzip,re,json,os,sys
from collections import Counter
from zlib_ng import gzip_ng
import argparse

#%%
def get_fastq_reader(path_file):
    if not os.path.exists(path_file):
        return iter([])
    if path_file.endswith("gz"):
        reader = gzip_ng.open(path_file, 'rb')
    else:
        reader = open(path_file, 'r')
    return reader

#%%
def get_fastq_writer(path_file):
    if path_file.endswith("gz"):
        reader = gzip_ng.open(path_file, 'wb')
    else:
        reader = open(path_file, 'w')
    return reader

#%%
def close_buffer(file_buffer):
    if hasattr(file_buffer, "closed") and not file_buffer.closed:
        file_buffer.close()

#%%
def read_line_from_reader(reader, replace):
    line = next(reader, replace)
    if line != None and hasattr(line, "decode"):
        line = line.decode()
    return line   

#%%
def get_single_read_set(reader):
    seqid = read_line_from_reader(reader, None)
    seq = read_line_from_reader(reader, None)
    optional = read_line_from_reader(reader, None)
    quality = read_line_from_reader(reader, None)
    return (seqid, seq, optional, quality)

#%%
def check_read_set_valid(read_set):
    list_is_valid = list(map(lambda x : x != None, read_set))
    is_set_valid = (sum(list_is_valid) == 4)
    is_set_end = (sum(list_is_valid) == 0)
    return is_set_valid, is_set_end

#%%
def get_sequence_identifier_from_read_set(read_set):
    seqid = read_set[0].strip('\n')
    assert seqid.startswith('@'), f"Sequence ID must start with '@' : {seqid}"
    seqid_info = seqid[1:]
    return seqid_info

#%%
def collect_header_info(reader):
    curr_num = 0
    read_identifier = list()
    while True:
        read_set = get_single_read_set(reader)
        is_set_valid, is_set_end = check_read_set_valid(read_set)
        if is_set_end: 
            print(f"Final Read Count : {curr_num}", flush = True, file = sys.stderr)
            break
        if not is_set_valid:
            print(f"{reader.name} Error: File has unrecognizable read. Pass this read.\n\t{read_set[0]}\n\t{read_set[1]}\n\t{read_set[2]}\n\t{read_set[3]}\n", flush = True, file = sys.stderr)
            continue
        seqid_info = get_sequence_identifier_from_read_set(read_set)
        read_identifier.append(seqid_info)
        curr_num += 1
        if curr_num % 1000000 == 0:
            print(f"Current Read Count : {curr_num}", flush = True, file = sys.stderr)
    return read_identifier
        
#%%
def check_header_unique(list_read_identifier):
    counter_read_id = Counter(list_read_identifier)
    duplicated_read_names = list(filter(lambda x : counter_read_id[x]>1, counter_read_id.keys()))
    num_duplicated_reads = list(map(lambda x : counter_read_id[x], duplicated_read_names))
    duplicated_reads = dict(zip(duplicated_read_names, num_duplicated_reads))
    print(f"Duplicated Reads Detected : {sum(num_duplicated_reads)} / {len(list_read_identifier)}", flush = True)
    return duplicated_reads

#%%
def write_line_by_writer(writer, line):
    line = line.strip('\n') + '\n'
    if hasattr(writer, "compress"):
        line = line.encode()
    writer.write(line)

##%
def write_read_set(writer, read_set):
    for line in read_set:
        write_line_by_writer(writer, line)
    
#%%
def write_deduplicated_fastq_file(fastq_reader, fastq_writer, save_unsafe, list_read_identifier, duplicated_reads):
    curr_num = 0
    dict_duplicated_read_written = dict()
    
    def is_read_okay_to_write(seqid_info, seq, qual):
        if duplicated_reads.get(seqid_info) == None:
            return True
        elif dict_duplicated_read_written.get(seqid_info) == None:
            dict_duplicated_read_written[seqid_info] = [(seq,qual)]
            return True
        else:
            return False    
        
    def is_safe_duplication(seqid_info, seq, qual):
        return (seq, qual) in dict_duplicated_read_written[seqid_info]
    
    unsafe_writer = None if save_unsafe == None else get_fastq_writer(save_unsafe)
    while True:
        read_set = get_single_read_set(fastq_reader)
        is_set_valid, is_set_end = check_read_set_valid(read_set)
        if is_set_end: 
            print(f"Final Writing Count : {curr_num} / {len(list_read_identifier)}", flush = True, file = sys.stderr)
            break
        if not is_set_valid:
            continue
        seqid_info = get_sequence_identifier_from_read_set(read_set)
        seq = read_set[1]
        qual = read_set[3]
        
        write_this_set = is_read_okay_to_write(seqid_info, seq, qual)
        if not write_this_set:
            safe_duplication = is_safe_duplication(seqid_info, seq, qual)
            print(f"Duplicated Read Detected and Pass this Read : {seqid_info}", flush = True, file = sys.stderr)
            if not safe_duplication:
                print(f"\nUnsafe Duplicated Read Detected.\nPlease check the raw fastq file!!! : {seqid_info}\n", flush = True)
                if unsafe_writer != None:
                    write_read_set(unsafe_writer, read_set)
        else:
            write_read_set(fastq_writer, read_set)
        
        curr_num += 1
        if curr_num % 1000000 == 0:
            print(f"Current Writing Count : {curr_num} / {len(list_read_identifier)}", flush = True, file = sys.stderr)
    close_buffer(unsafe_writer)
            
#%%
def run(path_fastq, path_save, save_unsafe):
    fastq_reader = get_fastq_reader(path_fastq)
    print(f"Start reading {path_fastq}...", flush = True)
    list_read_identifier = collect_header_info(fastq_reader)
    print(f"Reading {path_fastq} Finished", flush = True)
    close_buffer(fastq_reader)
        
    duplicated_reads = check_header_unique(list_read_identifier)
    
    if len(duplicated_reads.keys()) > 0:
        fastq_reader_to_write = get_fastq_reader(path_fastq)
        fastq_writer = get_fastq_writer(path_save)
            
        write_deduplicated_fastq_file(fastq_reader_to_write, fastq_writer, save_unsafe, list_read_identifier, duplicated_reads)
        print("Deduplicated Fastq File saved.", flush = True)
        close_buffer(fastq_reader_to_write)
        close_buffer(fastq_writer)
        
    else:
        print("No Duplicated Reads.", flush = True)
        try:
            os.symlink(path_fastq, path_save)
        except Exception as e:
            print("Symbolic Link Failed : {e}", flush = True)

#%%
if __name__ == "__main__":   
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--input")
    argument_parser.add_argument("--out")
    argument_parser.add_argument("--unsafe")
    
    args = argument_parser.parse_args()
    
    run(args.input, args.out, args.unsafe)
  
