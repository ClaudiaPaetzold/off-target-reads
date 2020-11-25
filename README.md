# off-target-reads

A workflow aiming to assemble off-target reads from target enrichment alignments, annotate the contigs, exclude putatively contaminant sequences, find putative single or low copy genes recovered in a user-specified number of samples and create alignments from them. Starting point is a \*.sam or \*.bam file of reads from a target enrichment approach mapped to the respective target regions or baits. The workflow comprises steps running in several established OpenSource programs (e.g. [Trinotate](#https://github.com/Trinotate/Trinotate.github.io/wiki) or [samtools](#https://github.com/samtools/samtools)) and some python scripts to filter their output. This repository contains the python script and the individual steps through the workflow.
<br>
If you use the workflow or any of its parts, please cite the repository as well as Reichelt et al., 2020 (Link will be provided upon publication of the paper).

## Summary

  - [Getting Started](#getting-started)
  - [Workflow](#workflow)
  - [Built With](#built-with)
  - [Contributing](#contributing)
  - [Authors](#authors)
  - [License](#license)
  - [Acknowledgments](#acknowledgments)

## Getting Started

These instructions will get you a copy of the project up and running on
your local machine. However depending on the number and size of your input files,
running some steps on a HPC is strongly advised. See [Workflow](#workflow) for details.

### required Additional Software

  - [samtools](#https://www.htslib.org)
  - [Picard](#https://broadinstitute.github.io/picard/)
  - A de-novo assembly program suited to your reads (e.g. [Bowtie2](#https://github.com/BenLangmead/bowtie2), [segemehl](#http://www.bioinf.uni-leipzig.de/Software/segemehl/), [SPAdes](#https://github.com/ablab/spades), etc.)
  - [BLAST+](#https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  - [Trinotate](#https://github.com/Trinotate/Trinotate.github.io/wiki)
  - A multiple sequence alignment program (e.g. [MAFFT](#https://mafft.cbrc.jp/alignment/software/) , [MUSCLE](#https://www.drive5.com/muscle/), etc.)


### Prerequisites

   - Python 3
   - Pandas
   - Biopython
   - AMAS


### Installing

The most convenient way to install the entirety of the required packages and software, is using the Conda (Python 3 distribution); either [Anaconda](#https://www.anaconda.com/) or [Miniconda](#https://docs.conda.io/en/latest/miniconda.html).

However, depending on the size of your dataset, some steps (e.g. the assembly) might perform better in a high performance environment.


## Workflow

The Workflow consists of 10 individual steps, some of which are optional depending on the intended use of the data.
<br><br>

#### Overview

1) Samtools to extract unmapped sequences<br>
2) Convert sam to fastq files with Picard<br>
3) De-novo assembly of off-target reads<br>
4) Blast <br>
5) Filter BLAST results using **FilterBlastByHitLenght.py**<br>
6) Annotation using Trinotate<br>
7) **Annotation2Fastas.py** to filter for putative contamination contigs and write genes meeting user-specified criteria in terms of sample coverage and hit number into individual fasta files<br>
8) Multiple Sequence alignment of each fasta file<br>
9) Hands-on, visual check of alignments - delete sequences not overlapping or poor alignments<br>
10) Trim alignments - **TrimGhostGaps.py** and/or **TrimEnds.py**<br>

#### Details

**1) Samtools to extract unmapped sequences**
<br><br>
Here is an example for paired end reads:
<br><br>
        ```samtools view -f 4 {filename}.bam > {filename}.unmapped.sam```
<br><br>
**2) Convert \*.sam file to \*.fastq file**
<br><br>
Here is the code for extracting the two paired end reads from the unmapped.sam file into individual \*.fastq files using the Picard suite.
<br><br>
```
        java -jar picard.jar SamToFastq I={filename}.unmapped.sam FASTQ={filename}_unmapped.1.fastq SECOND_END_FASTQ={filename}_unmapped.2.fastq
```
<br><br>
**3) De-novo assembly of off-target reads**
<br><br>
Use an assembler fitting your data. Example here SPAdes for PE reads:
<br><br>
```
python /usr/product/bioinfo/SPADES/3.13.2/spades.py -1 [unmapped1].fastq.gz -2 [unmapped2].fastq -k 21,33,55,77 -t 6 --careful --only-assembler -o [outfilename]_spades
```
<br><br>
**4) Blast**
<br><br>
Blast search to identifiy the contigs. Use either BlastX or BlastP and report to BLAST-outfmt6. If you use BlastP, first translate your input sequences using [TransDecoder](#https://github.com/TransDecoder/TransDecoder/wiki). Download the current verison of [uniprot](#https://uniprot.org) as a reference. <br>
Example:

Identify longest open reading frames:
```
TransDecoder.LongOrfs -t [samplename_contigs].fasta

```
Prepare the BLAST database from the uniprot archive:
```
makeblastdb -in uniprot_sprot.pep -dbtype prot

```
Perform a BLASTP search using the identified ORFs as queries:
```
blastp -query [samplename]_LongestOrfs.pep -db uniprot_sprot.pep -num_threads [num] -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > [samplename]_blastp.outfmt6

```
<br><br>
**5) Filter BLAST results**
<br><br>
This step is optional. The script **FilterBlastByHitLenght.py** will filter BLAST results and retain only hits with a user specified minimum required length for the hit alignment. This step might remove more spurious results, but its main purpose is to facilitate a possible bait or primer design from final alignments.
<br><br>
```
usage: FilterBlastByHitLenght.py [-h] [-m MINLENGTH] [-c COLUMN] Path

filters blast outfmt6 files for entries with a user specified minimum hit
alignment length and writes resutls to new file

positional arguments:
  Path                  directory path to blast result files - there should be
                        no other files in this

optional arguments:
  -h, --help            show this help message and exit
  -m MINLENGTH, --minLength MINLENGTH
                        minimum length of the alignment of blast hit. Default:
                        120
  -c COLUMN, --column COLUMN
                        Number of the column containing the blast alnlen
                        metric. Default: 4 (=standard outfmt6 format)
```
<br><br>
**6) Annotation using Trinotate**
<br><br>
Build a Trinotate database for each sample and export it (see [Manual](#https://github.com/Trinotate/Trinotate.github.io/wiki/Loading-generated-results-into-a-Trinotate-SQLite-Database-and-Looking-the-Output-Annotation-Report) if necessary). Use the (filtered) BLAST results.<br>
**Important**: Make sure to include the transcripts in the annotation reports by setting the ``` --incl_trans``` flag.
<br><br>
**7) Annotation results to Fasta files**
<br><br>
This script will parse all annotation results and produce per gene fasta files containing the transcript-sequences for identified genes found in a user-specified minimum number of samples and with a maximum number of occurrences per sample.
<br><br>
```
usage: Annotation2Fastas.py [-h] -n NUM_SAMPLES_CUTOFF [-r REPETITIVE_CUTOFF]
                            [-b [{P,X}]]
                            PATH

Filters Trinotate annotation report according to BlastP results in two steps:
1) keep only transcripts with blast hits, then 2) remove transcripts mapping
to non-spermatophyte references. It will create a file 'aliens' listing the
transcripts filtered in step 2. Please make sure to include the transcripts in
the Trinotate report.

positional arguments:
  PATH                  path to directory containing Trinotate results

optional arguments:
  -h, --help            show this help message and exit
  -n NUM_SAMPLES_CUTOFF, --Num_Samples_cutoff NUM_SAMPLES_CUTOFF
                        <number> genes are reported when found in more that
                        this number of samples
  -r REPETITIVE_CUTOFF, --repetitive_cutoff REPETITIVE_CUTOFF
                        <number> if a gene is found in any one species more
                        than this number of times it is regarded as a
                        putatively repetitive region and excluded. Default: 0
                        (no exclusion)
  -b [{P,X}], --blast_algorithm [{P,X}]
                        <P or X> choose which BLAST results will be the basis
                        for filtering. Default (blast)P

```
<br><br>
**8) Multiple Seuqence Alignment**
<br><br>
Build multiple sequence alignments (MSAs) using the previously generated \*.fasta files as input. If possible choose attempt to align either direction of a given sequence (``` --adjustdirection``` in MAFFT). See this example:
```
mafft --maxiterate 100 --adjustdirection --localpair --thread 5 [genesymbol].fasta > [genesymbol]_mafft.fasta
```
<br><br>
**9) Check alignments**
<br><br>
Visually check your alignments. Currently, **Annotation2Fasta.py** does not ensure that all transcripts overlap. Hence, it is possible that you have sequences that do not match the rest of the alignments, especially if you allow several transcripts per sample. This does not need to be the case, but should be inspected and edited if neccessary.
<br><br>
**10) Trim alignments**
<br><br>
Regardless of any editing, most alignments will have regions of low sample coverage in the beginning and end. **TrimEnds.py** allows you to cut the columns of an alignment start and end until a minimum sample coverage is reached.
```
usage: TrimEnds.py [-h] [-c COVERAGE] [-g GAP_CHARACTER] Path

trims ends of all fasta files in directory by user specified sequence coverage
and writes the output to a new file

positional arguments:
  Path                  path to directory containing fasta files to trim

optional arguments:
  -h, --help            show this help message and exit
  -c COVERAGE, --coverage COVERAGE
                        required fraction of sequences covering ends; if less,
                        positions will be trimmed. Default: 0.5
  -g GAP_CHARACTER, --gap_character GAP_CHARACTER
                        character denoting a gap. Default: -
```
<br>
If you have edited one or several of your alignments, you might have now columns containing gaps only. Some tools in downstream analysis might report errors when encountering these gappy columns. **TrimGaps.py** allows you to cut columns containing only the gap character.

```
usage: TrimGhostGaps.py [-h] Path

trims all fasta files in directory from gap only colums

positional arguments:
  Path        path to directory containing fasta files to degap

optional arguments:
  -h, --help  show this help message and exit
```
<br><br>



## Disclaimer

The workflow and scripts are the result of an investigation of off-target reads. They were conceived during the writing of a paper on Target Enrichment in
<br><br> *Reichelt et al (2020)* **Target enrichment improves phylogenetic resolution in the genus *Zanthoxylum* (Rutaceae) and indicates both incomplete lineage sorting and hybridization events**. <br><br>
They do not represent a fully matured, wrapped, parallelized pipeline, but a basic working concept, that still requires some hands-on work.
That said, I do encourage any kind of further development. Please contact me if you are interested in advancing the workflow.

## Built With

  - [Creative Commons](#https://creativecommons.org/) - Used to choose
    the license


## Authors

  - **Claudia Paetzold** - *conceived the workflow and wrote the scripts*
    [ClaudiaPaetzold](#https://github.com/ClaudiaPaetzold)
  - **Billie Thompson** - *Provided README Template* -
    [PurpleBooth](#https://github.com/PurpleBooth)
  - **Marc Appelhans** - *provided help and feedback conceptualising the workflow*
  - **Niklas Reichelt**(#https://github.com/NiklasReichelt) - *provided feedback on the repo*


## License

This project is licensed under the [GNU General Public License v3.0](LICENSE.md)
Creative Commons License - see the [LICENSE.md](LICENSE.md) file for
details
