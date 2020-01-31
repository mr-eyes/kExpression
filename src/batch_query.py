import subprocess

import click
import kProcessor as kp
from tqdm import tqdm

from src.click_context import cli


class Batch_Query:
    index_prefix = str()
    fasta_file = str()
    chunk_size = int()
    seqs_no = int()
    kSize = int()

    def __init__(self, logger_obj, index_prefix, fasta_file, chunk_size, kSize):
        self.Logger = logger_obj
        self.index_prefix = index_prefix.replace(".mqf", '')
        self.chunk_size = int(chunk_size)
        self.fasta_file = fasta_file
        self.kSize = kSize

        process = subprocess.Popen(['wc', '-l', fasta_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        self.seqs_no = int(str(out, encoding='utf8').split()[0]) // 2

    @staticmethod
    def CountFrequency(my_list):

        # Creating an empty dictionary
        freq = {}
        length = len(my_list)
        for items in my_list:
            freq[items] = my_list.count(items)

        for k, v in freq.items():
            freq[k] = round((v / length) * 100)

        return (freq)

    @staticmethod
    def freq_to_line(freq_dict):
        new_freq = {k: v for k, v in sorted(freq_dict.items(), key=lambda item: item[1])}
        line = list()
        for tr, per in new_freq.items():
            line.append(f"{per}:{tr}")

        return ",".join(line)

    def process(self):
        self.Logger.INFO("Loading index...")
        idx = kp.colored_kDataFrame.load(self.index_prefix)
        names_map = idx.names_map()
        chunks_no = self.seqs_no // self.chunk_size
        self.Logger.INFO("Processing...")

        batchQuery = kp.ckf_batchQuery(idx, self.fasta_file, {"mode": 1, "k_size": self.kSize}, self.chunk_size)

        with open(self.fasta_file + ".tsv", 'w') as result_writer:
            header = "unitig\tscore[%:tr_id]"
            result_writer.write(header + '\n')

            for _ in tqdm(range(chunks_no)):
                batchQuery.next()
                result = batchQuery.get_transcripts()
                for read_id, kmers in result.items():
                    flattened = []

                    for kmer in kmers:
                        flattened += [x for x in kmer]

                    freq_dict = self.CountFrequency(flattened)
                    result_writer.write(read_id + '\t' + self.freq_to_line(freq_dict) + '\n')

        self.Logger.SUCCESS(f"Result saved in {self.fasta_file}.tsv")


def validate_kSize(ctx, param, value):
    if not value % 2:
        raise click.BadParameter(f"kmer size: {value} is even, please enter an odd value.")
    return value


@cli.command(name="batchQuery", help_priority=1)
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING,
              help="kProcessor index file prefix")
@click.option('-s', '--chunk-size', 'chunk_size', required=False, type=click.INT, default=100, show_default=True,
              help="chunk size, how many reads to process each iteration?")
@click.option('-k', '--kmer-size', "kSize", callback=validate_kSize, required=True,
              type=click.IntRange(7, 31, clamp=False), help="kmer size")
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.pass_context
def batchQuery(ctx, index_prefix, chunk_size, fasta_file, kSize):
    """Query multiFasta file."""

    BQ = Batch_Query(logger_obj=ctx.obj, index_prefix=index_prefix, fasta_file=fasta_file, chunk_size=chunk_size,
                     kSize=kSize)
    BQ.process()
