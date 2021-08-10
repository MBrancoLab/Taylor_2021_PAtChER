import sys
from Bio import SeqIO
import re
import argparse


def get_args(argv):
    """
    Parse arguments from command line
    """
    args = argparse.ArgumentParser()
    args.add_argument('-r1', required=True,
                      metavar='<fastq>',
                      help='R1 fastq file')
    args.add_argument('-r2', required=True,
                      metavar='<fastq>',
                      help='R2 fastq file')
    return args.parse_args()


def process(r1, r2):
    """
    Convert read IDs into bed entries
    """
    ##Process R1
    out1 = open(r1.replace('fastq','bed'), mode='w')
    for read in SeqIO.parse(r1,'fastq'):
        rid = read.description.split()
        if 'HIC' in rid:
            coord = rid[2].split(':')
            chrom = coord[0]
            start = coord[1]
            end = int(start) + len(read.seq)
            strand = '.'
        else:
            coord = re.split(':|\.\.', rid[2])
            chrom = coord[0]
            start = coord[1]
            end = int(start) + len(read.seq)
            if coord[3]=='F':
                strand='+'
            else:
                strand='-'
        out1.write(f'{chrom}\t{start}\t{end}\t{rid[0]}\t.\t{strand}\n')
    out1.close()

    ##Process R2
    out2 = open(r2.replace('fastq','bed'), mode='w')
    for read in SeqIO.parse(r2,'fastq'):
        rid = read.description.split()
        if 'HIC' in rid:
            coord = rid[3].split(':')
            chrom = coord[0]
            end = coord[1]
            start = int(end) - len(read.seq)
            strand = '.'
        else:
            coord = re.split(':|\.\.', rid[2])
            chrom = coord[0]
            end = coord[2]
            start = int(end) - len(read.seq)
            if coord[3]=='F':
                strand='-'
            else:
                strand='+'
        out2.write(f'{chrom}\t{start}\t{end}\t{rid[0]}\t.\t{strand}\n')
    out2.close()


if __name__ == '__main__':
    args = get_args(sys.argv[1:])
    process(args.r1, args.r2)

