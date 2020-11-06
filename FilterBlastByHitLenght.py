#!/bin/python

import os
import glob
import argparse


def filter_blast(blastfile, minlen, column):
    # input blastoutput-file, minimum length required of hit alignment - AN: uses standard blast outfmt6 with alnlen=column4
    spec = os.path.basename(blastfile).split('.')[0]
    parent = os.path.dirname(blastfile)
    filtered_blast = os.path.join(parent, spec + '_filtered.outfmt6')
    with open(blastfile) as blast:
        with open(filtered_blast, 'a') as out:
            for line in blast:
                alnlen = line.split()[column - 1]
                if int(alnlen) >= minlen :
                    out.write('{}\n'.format(line.strip()))


def main():

###########################
#parse arguments

	parser = argparse.ArgumentParser(description='filters blast outfmt6 files for entries with a user specified minimum hit alignment length and writes resutls to new file')
	parser.add_argument('Path', type=str, help='directory path to blast result files - there should be no other files in this')
	parser.add_argument('-m', '--minLength', type=float, help='minimum length of the alignment of blast hit. Default: 120', default=120)
	parser.add_argument('-c', '--column', type=int, help='Number of the column containing the blast alnlen metric. Default: 4 (=standard outfmt6 format)', default=4)
	
	args = parser.parse_args()
	
	
###########################					
	for f in glob.glob(os.path.join(args.Path, '*')):
		filter_blast(f, args.minLength, args.column)
		
if __name__ == '__main__':
	main()