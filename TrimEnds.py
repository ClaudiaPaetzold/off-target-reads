#!/usr/bin/python
import glob
import os
import argparse
from Bio import AlignIO




def ends4trimming(aln, fraction, gapchar):
	### helper function; finds the position(s) in the alignment where the alignment will be trimmed
    # aln = handle of parsed AlignIO object
    # fraction = fraction of species covered at alignment ends - float
    num_spec = len(aln)
    frac = int(num_spec - (num_spec * fraction))
    start_gap = []
    end_gap = []

	# loop though alignment and compare coverage of current position with the next
    for col in range(aln.get_alignment_length() -1):
        column = aln[:, col]
        next_col = aln[:, col + 1]
        if col == 0 and column.count(gapchar) < frac:
            end_gap.append(col)
        elif col == aln.get_alignment_length() -2:
            if column.count(gapchar) < frac:
                start_gap.append(col)
        else:
            if column.count(gapchar) > frac and next_col.count(gapchar) <= frac:
                end_gap.append(col)
            elif column.count(gapchar) < frac and next_col.count(gapchar) >= frac:
                start_gap.append(col)
    
    # check if last position in alignment exceeds coverage requirements -> then it needs to be the last element of startgap
    last = aln.get_alignment_length() 
    if aln[:,last-1].count(gapchar) < fraction:
        start_gap.append(last-1)
    
    return [end_gap[0], start_gap[-1]]

def trim_alignment(aln, outname, fraction, gapchar):
    # aln: handle of AlignIO-pased file
    # outname: file to write to
    # fraction = fraction of species covered at alignment ends - float

    pos = ends4trimming(aln, fraction, gapchar)
    trimmed_aln = aln[:, pos[0] + 1 : pos[1] -1]
    AlignIO.write(trimmed_aln, outname, 'fasta')


	
def main():

###########################
#parse arguments

	parser = argparse.ArgumentParser(description='trims ends of all fasta files in directory by user specified sequence coverage and writes the output to a new file')
	parser.add_argument('Path', type=str, help='path to directory containing fasta files to trim')
	parser.add_argument('-c', '--coverage', type=float, help='required fraction of sequences covering ends; if less, positions will be trimmed. Default: 0.5', default=0.5)
	parser.add_argument('-g', '--gap_character', type=str, help='character denoting a gap. Default: -', default='-')
	
	args = parser.parse_args()
	
	
###########################

	for file in glob.glob(os.path.join(args.Path, '*.fasta')):
		
	
		gene =  os.path.basename(file).split('.')[0]

		parent = os.path.dirname(file)

		trim_alignment(AlignIO.read(file, 'fasta'), os.path.join(parent, gene + '_trimmed.fasta'), args.coverage, args.gap_character)

if __name__ == '__main__':
	main()