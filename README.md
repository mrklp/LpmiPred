# LpmiPred



## About

LpmiPred is a pipeline that used to identify large sets of candidate LP-miRNAs from small RNA (sRNA) sequencing data.

## Requirements

[dnapi.py](https://github.com/jnktsj/dnapi): De novo adapter prediction algorithm.

[fastp](https://github.com/OpenGene/fastp): A tool designed to provide fast all-in-one preprocessing for FastQ files.

[bowtie](https://bowtie-bio.sourceforge.net/index.shtml): An ultrafast, memory-efficient short read aligner. 

[RNAfold](http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi): Predict secondary structures of single stranded RNA or DNA sequences.

## Pipeline 

![](https://github.com/mrklp/LpmiPred/blob/main/LpmiPred.png)

### The first step:

Sequencing data of 37,820 small RNA libraries from six representative species (silkworm, fruit fly, nematode, zebrafish, mouse, and humans)  were downloaded and preprocessed, including quality control, removal of barcodes and adapters, and collapsing consistent reads. Finally, "collapse_reads_all.fasta" files were obtained for subsequent steps.

```shell
# removal of barcodes and adapters, quality control
python de_barcode_adapter.py  -q file.fastq -m  miRNAs.txt -s species  -o clean.fastq
# collapsing consistent reads for each small RNA libraries
python collapse_reads_single.py -a clean.fasta -o collapse_reads.fasta
# merge all librarys
python collapse_reads_merge.py -f ./collapse_reads_single/ -o collapse_reads_all.fasta
```

### The second step:

Excising potential pre-miRNA from genome.

```shell
# build index
bowtie-build --threads 10 genome.fa genome
# Reads mapped against genome
bowtie -p 10 -f -n 0 -a -m 5 --best --strata -x genome collapse_reads_all.fasta reads_mappings.bwt

# Excising potential pre-miRNA 
python excise_precursors.py -a 20  -b reads_mappings.bwt -g genome.fa -c chr1 -o precursors_chr1
```

### The third step:

Discovering high confidence miRNA using the signature and structures as guidelines.

```shell
# prepare signature file
bowtie-build --threads 10 precursors.fa precursors
bowtie -p 10 -f -v 1 -a --best --strata --norc -x precursors  collapse_reads_all.fasta signature.bwt
# prepare structure file
RNAfold < precursors.fa --noPS > precursors.str
# Convert the signal file from '.bwt' to '.h5'
python bwt_to_h5.py -s signature.bwt -o signature.h5
# Discovering high confidence miRNA
python filter_precursors.py  -n 20 -f precursors.fa -c precursors.coords  -s signature.h5 -p precursors.str  -o result.txt
```

