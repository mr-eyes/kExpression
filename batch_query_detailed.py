import kProcessor as kp
import sys
import subprocess
from tqdm import tqdm


def CountFrequency(my_list): 
      
    # Creating an empty dictionary  
    freq = {}
    length = len(my_list)
    for items in my_list: 
        freq[items] = my_list.count(items)
    
    for k, v in freq.items():
        freq[k] = round((v/length) * 100)

    return(freq)

def freq_to_line(freq_dict):
    new_freq = {k: v for k, v in sorted(freq_dict.items(), key=lambda item: item[1])}
    line = list()
    for tr, per in new_freq.items():
        line.append(f"{per}:{tr}")
    
    return ",".join(line)


if len(sys.argv) != 4:
    sys.exit("run: python batch_query_detailed.py <index_prefix> <unitigs_file> <chunk_size>")

index_prefix = sys.argv[1].replace(".mqf", '')
unitigs_file = sys.argv[2]
chunk_size = int(sys.argv[3])

process = subprocess.Popen(['wc', '-l', unitigs_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = process.communicate()
seqs_no = int(str(out, encoding='utf8').split()[0]) // 2


print("Loading index...")
idx = kp.colored_kDataFrame.load(index_prefix)
names_map = idx.names_map()

chunks_no = seqs_no // chunk_size

print("Processing...")

batchQuery = kp.ckf_batchQuery(idx, unitigs_file, {"mode" : 1, "k_size" : 25} , chunk_size)

all_results = dict()

with open(unitigs_file + ".tsv", 'w') as result_writer:
    header = "unitig\tscore[%:tr_id]"
    result_writer.write(header + '\n')


    for _ in tqdm(range(chunks_no)):
        batchQuery.next()
        result = batchQuery.get_transcripts()
        for read_id, kmers in result.items():
            score = dict()
            flattened = []
            
            for kmer in kmers:
                flattened += [x for x in kmer]
            
            freq_dict = CountFrequency(flattened)
            result_writer.write(read_id + '\t' + freq_to_line(freq_dict) + '\n')