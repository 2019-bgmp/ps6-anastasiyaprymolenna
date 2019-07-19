#!/usr/bin/env python3
import argparse
import re

def get_args(): # Usage
        parser = argparse.ArgumentParser(description='Take no more than 3 FASTQ files of sequences to calculate expected coverage')
        parser.add_argument("-f1", "--file1", help = "This needs to be a FASTQ file of sequences. Sequences must be contained on one line.", required = True)
        parser.add_argument("-f2", "--file2", help = "Extra FASTQ file")
        parser.add_argument("-f3", "--file3", help = "Extra FASTQ file")
        parser.add_argument("-g", "--genome", help = "Input the size of the genome length")
        parser.add_argument("-k", "--kmer", help = "Input the kmer size to calculate expected kmer cover based off of it")
        return parser.parse_args()#!/usr/bin/env python3
args = get_args()

file1 = args.file1
file2 = args.file2
file3 = args.file3
genome_len = int(args.genome)

def seq_len(read_file):
    total_sequence_length = 0
    with open(read_file,"r") as file: #opens input file for reading
        ln = 0 #start count to pull out only sequence lines from FASTQ file
        for line in file:
            line = line.strip()
            if ln%4==1: #only work with sequence lines
                length = len(line)
                total_sequence_length += length
            ln+=1 #increment line counter
        NR = ln/4
    return (NR, total_sequence_length)

reads1, f1_len = seq_len(file1)
reads2, f2_len = seq_len(file2)
reads3, f3_len = seq_len(file2)

length_of_reads = f1_len+f2_len+f3_len
avg_read_len = length_of_reads/(reads1+reads2+reads3)
coverage = length_of_reads / genome_len

print("The expected coverage is ", coverage)
k_size = int(args.kmer)
kmer_cov = coverage * (avg_read_len-k_size) / avg_read_len

print("The k-mer coverage is ", kmer_cov)
