#! /usr/bin/python3

from Bio import AlignIO
import glob
import os
from amas import AMAS
import argparse


def columns2keep(aln):
	# helper function, finding the column numbers of non-empty columns
	# input: a parsed AlignIO object
    # helper variable: number of species
    num_spec = len(aln)
    
    #counter variables
    start_gap = 0
    end_gap = 0
    #c = 0
    
    # helper lists for start end end of gaps positions
    skeep = []
    ekeep = []


    for col in range(aln.get_alignment_length()-1):
        column = aln[:, col]
        next_col = aln[:, col + 1]
        if col == 0:
            if column.count('-') != num_spec:
                
                ekeep.append(col)
                
        elif col == aln.get_alignment_length() -1:
            if column.count('-') != num_spec:
                
                skeep.append(col)
        else : 
            prev_col = aln[:, col - 1]
            if column.count("-") == num_spec:
               
                if prev_col.count('-') < num_spec:
                    skeep.append(col)
                if next_col.count('-') < num_spec:
                    ekeep.append(col)

    skeep.append(aln.get_alignment_length())

    keeplist = list(zip(ekeep, skeep))
    
    return keeplist


def write_ungapped_fasta(aln, outname):
	# write fasta less the ghost gaps column 
	# will create temporary alignment files by splitting the original at the gap-only-columns...
	# ... then merge the temps into one alignment and delete the temp files
    keep = columns2keep(aln)

    tmplist = []
    if len(keep) > 0:
        for i in keep:
            tempfile = 'temp_' + str(keep.index(i)) +'.fasta'
            tmplist.append(tempfile)
            tempaln = aln[:, i[0] +1 :i[1]]
            AlignIO.write(tempaln, tempfile, 'fasta')
    else:

        tempfile = 'temp_0.fasta'
        tmplist.append(tempfile)
        AlignIO.write(aln[:,:], tempfile, 'fasta')

    #now read in tmp files and concatenate to output files using AMAS.py
    meta_aln = AMAS.MetaAlignment(in_files=tmplist, data_type="dna",in_format="fasta", cores=1)
    parsed_alns = meta_aln.get_parsed_alignments()
    concat_tuple = meta_aln.get_concatenated(parsed_alns)
    concatenated_alignments = concat_tuple[0]
    new_aln = meta_aln.print_fasta(concatenated_alignments)
    with open(outname, 'w') as out:
        out.write(new_aln)
    for file in glob.glob('temp*.fasta'):
        os.remove(file)

		
		
#######################################################

def main():
#parse arguments
##############################
	parser = argparse.ArgumentParser(description='trims all fasta files in directory from gap only colums')
	parser.add_argument('Path', type=str, help='path to directory containing fasta files to degap') 
		
	args = parser.parse_args()
#############################	
	for file in glob.glob(os.path.join(args.Path, '*.fasta')):
    
		gene =  os.path.basename(file).split('.')[0]
		parent = os.path.dirname(file)
					
		write_ungapped_fasta(AlignIO.read(file, 'fasta'), os.path.join(parent, gene + '_ungapped.fasta'))

if __name__ == '__main__':
	main()
    

