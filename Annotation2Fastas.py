#!/usr/bin/python
import pandas as pd
import glob
import os
from collections import OrderedDict
import argparse
import fnmatch
import re
import random

# produce filtered files containing only trancript id, blastx hit and GO of hit

def filter_spermatophytes(infile, outfile, countfile):
#filter files further to contain only lines with "Spermatophyta" in BLAST Hit
#for all others write a separate table with counts for category

	aliens = OrderedDict()
	aliens["Embryophyta"] = 0
	aliens["Fungi"] = 0
	aliens["Rodentia"] = 0
	aliens["Bacteria"] = 0
	aliens["Hominidae"] = 0
	aliens["Amphibia"] = 0
	aliens["Mollusca"] = 0
	aliens["Insecta"] = 0

	taxon = os.path.basename(infile).split('.')[0]
	with open(infile, 'r') as f:
		with open(outfile, 'a') as filtered:
			f.readline()
			spermcount = 0
			count = 0
			for line in f:
				count += 1
				if "Spermatophyta" in line:
					spermcount += 1
					filtered.write('{}\n'.format(line.strip()))
				else:
					for key in aliens:
						if key in line:
							aliens[key] += 1

			rest = count - (spermcount + sum(aliens.values()))


	with open(countfile, 'a') as out:
		out.write('{T};{S};{Embryophyta};{Fungi};{Rodentia};{Bacteria};{Hominidae};\
			{Amphibia};{Mollusca};{Insecta};{R}\n'.format(T=taxon, S=spermcount, **aliens, R=rest))

def revc(seq):
## helper function: produce reverse compliment
## of a DNA sequence

	baseComplement = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A' }
	return "".join([baseComplement[base] for base in seq[::-1]])

def get_query_coords(blast_hit_entire):
## helper function: get the position of the query hit
##correcting for direction from the blastp output

	coords = blast_hit_entire.split('^')[2].split(',')[0][2:]
	if int(coords.split('-')[0]) > int(coords.split('-')[1]):
		return coords + '[-]'
	else:
		return coords + '[+]'

def cut2ORF(prot_coords, sequence):
## cut the transcript sequence to the query hit position only
## this is to make sure that transcripts can be aligned, which
## can be difficult otherwise (e.g. when there are introns)

	numbers = prot_coords[:-3]
	orientation = prot_coords[-2]

	if orientation == '+':
	   start = int(numbers.split('-')[0])
	   end = int(numbers.split('-')[1])
	   short_seq = sequence[start -1:end ]

	   return short_seq
	else:
	   start = int(numbers.split('-')[1])
	   end = int(numbers.split('-')[0])

	   short_seq = sequence[start -1:end ]
	   rev_comb = revc(short_seq)
	   return rev_comb


def count_abundance(directory):
## Count abundance of gene symbols in all samples
## Idea: create dictionary with gene symbols as key and nested dictionaries as
## values, containing sample name as key and number of occurrences as value
## adding to the dict if I find the gene symbol in the second column of the
## *.spermatophytes* file


	abundance = {}
	for item in glob.glob(os.path.join(directory, "*.spermatophytes.csv")):
		species = os.path.basename(item).split('.')[0]
		with open(item, 'r') as f:
			for line in f:
				symbol = line.split('|')[1].split('_')[0]
				if symbol[0] == '"':
					symbol = symbol[1:]
				#print symbol
				if symbol in abundance.keys():
					if species in abundance[symbol].keys():
						abundance[symbol][species] += 1
					else:
						abundance[symbol][species] = 1
				else:
					abundance[symbol] = {species : 1}
	return abundance


def minimum_depth(abundance_dict, cutoff):
## reports only genes found in at least >cutoff< number of samples
	symbols = abundance_dict.keys()
	reduced_ab = {}
	for symbol in symbols:
		if len(abundance_dict[symbol]) >= cutoff:
			reduced_ab[symbol] = abundance_dict[symbol]
	return reduced_ab


def exclude_putative_repetitives(reduced_dict, full_dict, cutoff):
## some gene symbols show very high species counts - filter these out as they likely represent repetitive stuff
## cutoff = user-specified
## first step: get list of gene symbols this pertains to
## steps:
## 1 lopp through redab dictionary
## 2 for each symbol (key) in dict - loop through sub-dictionary contained in value
## 3 in each sub-dict: for each species (key), check the value - if it is > 15, break loop, store symbol in list, move to next symbol

# putative repetitives
	put_repet = []

	genes = sorted(reduced_dict.keys())
	for symbol in genes:
		for species in full_dict[symbol]:
			if full_dict[symbol][species] > cutoff:
				put_repet.append(symbol)
				break
	#print(put_repet)
	#print(len(reduced_dict))

	for symbol in put_repet:
		del reduced_dict[symbol]
	print('removed %s genes as putatively repetitive' % len(put_repet))
	return reduced_dict



def main():

###########################
#parse arguments

	parser = argparse.ArgumentParser(description=
		" Filters Trinotate annotation report according to "
		" BlastP results in two steps: 1) keep only transcripts with blast hits, "
		" then 2) remove transcripts mapping to non-spermatophyte references."
		" It will create a file 'aliens' listing the transcripts filtered in step 2."
		" Please make sure to include the transcripts in the Trinotate report.")
	parser.add_argument('PATH', type=str, help='path to directory containing Trinotate results')
	parser.add_argument('-n', '--Num_Samples_cutoff', type=int, required=True,
		help='<number> genes are reported when found in more that this number of samples')
	parser.add_argument('-r', '--repetitive_cutoff', type=int,
		help='<number> if a gene is found in any one species more than this number of times '
		'it is regarded as a putatively repetitive region and excluded. Default: 0 '
		'(no exclusion)', default=0 )
	parser.add_argument('-b', '--blast_algorithm', type=str, help='<P or X> choose which BLAST results '
		'will be the basis for filtering. Default (blast)P',
		choices=('P', 'X'), default='P', const='P', nargs='?')
	parser.add_argument('-d', '--New_Directory', type=str,
		help='<string> name for the new directory containing the output *fasta files. '
		'Default: NewFastas' , default='NewFastas' )
	args = parser.parse_args()


###########################
	# create countfile
	countfile = os.path.join(args.PATH, "aliens.csv")
	with open(countfile, 'a') as out:
		out.write(';{};{};{};{};{};{};{};{};{};{}\n'.format('Spermatophyta',
				'Embryophyta', 'Fungi', 'Rodentia', 'Bacteria', 'Hominidae',
				'Amphibia', 'Mollusca', 'Insecta', 'Other'))


	# remove any line without blastp hit
	for f in glob.glob(os.path.join(args.PATH, '*.xls')):
		taxon = os.path.basename(f).split('_')[0]
		sample = pd.read_table(f, sep='\t')
		if args.blast_algorithm == 'P':
			sample_red = sample[["transcript_id", "sprot_Top_BLASTP_hit", "prot_coords", "transcript"]]
			export = sample_red[sample_red['sprot_Top_BLASTP_hit'] != "."]
		else:
			sample_red = sample[["transcript_id", "sprot_Top_BLASTX_hit", "prot_coords", "transcript"]]
			export = sample_red[sample_red['sprot_Top_BLASTX_hit'] != "."]
		filtered_csv = (os.path.join(os.path.dirname(f), taxon + ".filtered.csv"))
		# helper: name of first output file containing only entries with blastp hits

		export.to_csv(filtered_csv, sep='|', index=False)

		# filter out any non-spermatophyty blast hit
		spermatophyte = os.path.join(os.path.dirname(f), taxon + '.spermatophytes.csv')
		filter_spermatophytes(filtered_csv, spermatophyte, countfile)


	full_dict = count_abundance(args.PATH)

	# now sort the dictionary by the length of the values (=nested dictionaries),
	#brauch ich das? wirklich?
	#sorted_list = sorted(full_dict, key=lambda k: len(full_dict[k]), reverse=True)

	# as preparation for possible sorting (i.e. excluding) of loci due to high
	# hit count of genes per sample
	reduced = minimum_depth(full_dict, args.Num_Samples_cutoff)

	final = exclude_putative_repetitives(reduced, full_dict, args.repetitive_cutoff)

	# create new directory to put the resulting Fasta files in
	fastadir = os.path.join(args.PATH, args.New_Directory)
	if not os.path.isdir(fastadir):
		os.mkdir(fastadir)


	for csv in (glob.glob(os.path.join(args.PATH, "*.spermatophytes.csv"))):
		sample = os.path.basename(csv).split('.')[0]
		for symbol in final:
			with open(csv, 'r') as fas:
				for line in fas:
					if symbol in line:
						contigname = line.split('|')[0]
						header = ' | '.join([sample, contigname])
						sequence = line.strip().split('|')[-1]
						query_coords = line.split('|')[-2]
						if len(query_coords) < 5:
							query_coords = get_query_coords(line.split('|')[1])
						orf = cut2ORF(query_coords, sequence)
						fastafile = os.path.join(fastadir, (symbol + '.fasta'))
						with open(fastafile, 'a') as f:
							f.write('>{}\n{}\n'.format(header, orf.upper()))

	print('done')


if __name__ == '__main__':
	main()
