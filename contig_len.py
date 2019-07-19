#!/usr/bin/env python3
import argparse
import re

def get_args(): # Usage
	parser = argparse.ArgumentParser(description='Parse Fasta file for the longest protein sequence')
	parser.add_argument("-f", "--file", help = "This needs to be a FASTA file of contigs. Sequences must be contained on one line.", required = True)
	parser.add_argument("-k", "--kmer", help = "This is the size of the kmer used in the contig assembly", required = True)
	parser.add_argument("-o", "--out", help = "Specify the name of the file to write out statistical information to", required = True)
	return parser.parse_args()
args = get_args()

con_file = args.file
kmer_size = int(args.kmer)
out_file = args.out
wFile = open(out_file, "w") # open a file to write output to

contig_count = 0 #initialize contig count to total contigs
longest_contig = ""
contig_len_list = []
coverage_count = 0

rfile = open(con_file, "r")
for line in rfile:
	if ">" in line:
		line = line.strip()
		kmer_len = re.findall("length_[0-9]*", line)
		kmer_len = kmer_len[0].strip("length_") #get only number for kmer length
		kmer_cov = re.findall("[0-9]*.[0-9]*$", line) #get kmer cover
		kmer_cov = float(kmer_cov[0])
		phys_len = int(kmer_len)-kmer_size+1 # caculate string length from kmer length
		contig_count+=1 # sum number of contig records to get the toal number of contigs in the file
		contig = rfile.readline() # read line to get the contig sequence associated with the defline
		contig = contig.strip() #strip new line off of the sequence
		contig_len_list.append(len(contig))
		if len(contig) > len(longest_contig):
			longest_contig = contig
		coverage_count += kmer_cov
longest_contig_len = len(longest_contig) # calculate maximum contig length
genome_ass_len = sum(contig_len_list) # calculate the total length of the genome assembly across the contigs
mean_cont_len = (sum(contig_len_list))/len(contig_len_list) #calculate mean contig length
mean_cov_depth = coverage_count/len(contig_len_list)

#Write Stats to out file
wFile.write("Input file for the Velvet run and statistics:\t{}\n".format(con_file))
wFile.write("K-mer size used in contig assembly:\t{}\n".format(kmer_size))
wFile.write("Number of contigs:\t{}\n".format(len(contig_len_list)))
wFile.write("Maximum contig length:\t{}\n".format(longest_contig_len))
wFile.write("Mean contig length:\t{}\n".format(mean_cont_len))
wFile.write("Total length of the genome assembly across the contigs:\t{}\n".format(genome_ass_len))
wFile.write("Mean depth of coverage for the contigs:\t{}\n".format(mean_cov_depth))

# calculate N50 of the contig assembly
def calculate_N50(list_of_contig_len):
	total_genome_len = sum(list_of_contig_len)
	sorted_len = sorted(list_of_contig_len, reverse = True)
	halfway = total_genome_len/2
	sum_of_len = 0
	break_index = 0
	N50 = 0
	while sum_of_len <= halfway:
		sum_of_len += sorted_len[break_index]
		break_index += 1
	if sum_of_len == halfway:
		N50 = sorted_len[break_index]
	elif sum_of_len > halfway:
		N50 = sorted_len[break_index-1]
	return(N50)
N50 = calculate_N50(contig_len_list)
wFile.write("N50 of the contif assembly:\t{}\n\n".format(N50))
#Binning
sorted_list = sorted(contig_len_list, reverse = True) #work with all of the contig lengths in a sorted fashion so we can use the 0th index as the largest bin
bin_counts = {} # define a dictionary to store the tally of contig lengths that fit in each bin
bin_range_list = [] # make a list of possible limits with the right cut off numbers
limits = [] #make a list of tuples containing bin limits
for limit in range(-1, sorted_list[0], 100):
	bin_range_list.append(limit)
for index in range(len(bin_range_list)):
	limits.append((bin_range_list[index]+1, bin_range_list[index]+100))
for length in sorted_list: # work with the lengths of 
	for limit in limits :
		if (limit[0] <= length <= limit[1]):
			if limit not in bin_counts:
				bin_counts[limit] = 1
			if limit in bin_counts:
				bin_counts[limit] += 1
wFile.write("# Contig length\tNumber of contigs in this category\n")
for limits, counts in sorted(bin_counts.items()):
	wFile.write("{0}\t{1}\n".format(limits[0],counts))
