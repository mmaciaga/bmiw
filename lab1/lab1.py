from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import gzip
import statistics
import matplotlib.pyplot as plt


# zalozenie takie ze pliki sa spakowane
def read_file(filename):
    with gzip.open(filename, 'rt') as file:
        sequences = [seq_record.seq for seq_record in SeqIO.parse(file, 'fastq')]

    return sequences


def gc_content(sequences):
    gc_values = sorted(100 * gc_fraction(seq) for seq in sequences)
    mean_gc = sum(gc_values) / len(gc_values)
    return mean_gc

def count_reads(sequences):
    return len(sequences)


def count_nt(sequences):
    return sum(len(seq) for seq in sequences)


def mean_len(sequences):
    return statistics.mean(len(seq) for seq in sequences)


def median_len(sequences):
    return statistics.median(len(seq) for seq in sequences)


def histogram(sequences):
    sizes = [len(seq) for seq in sequences]
    plt.hist(sizes, bins=50, color='#720026')
    plt.xlim([0, max(sizes)])
    plt.title('read length distribution')
    plt.xlabel("sequence length (bp)")
    plt.ylabel("count")
    plt.show()

# no paskudny ten main ale co zrobię
def main():
    illumina_1 = read_file('IL_1.fastq.gz')
    illumina_2 = read_file('IL_2.fastq.gz')
    nanopore = read_file('NP.fastq.gz')
    mean_gc_IL_1 = gc_content(illumina_1)
    mean_gc_IL_2 = gc_content(illumina_2)
    mean_gc_IL = (mean_gc_IL_1 + mean_gc_IL_2) / 2
    mean_gc_NP = gc_content(nanopore)
    illumina_reads = count_reads(illumina_1) + count_reads(illumina_2)
    print(f'number of illumina reads: {illumina_reads}')
    print(f'number of nt in nanopore reads: {count_nt(nanopore)}')
    print(f'mean length of nanopore reads: {mean_len(nanopore):.2f}')
    print(f'median length of nanopore reads: {median_len(nanopore):.2f}')
    print(f'mean gc in illumina: {mean_gc_IL:.2f}%')
    print(f'mean gc in nanopore: {mean_gc_NP:.2f}%')
    histogram(nanopore)


if __name__ == '__main__':
    main()
