# LpmiPred



## About

LpmiPred is a pipeline that used to identify large sets of candidate LP-miRNAs in from small RNA (sRNA) sequencing data.

## Requirements

[dnapi.py](https://github.com/jnktsj/dnapi): De novo adapter prediction algorithm.

[fastp](https://github.com/OpenGene/fastp): A tool designed to provide fast all-in-one preprocessing for FastQ files.

## Pipeline 

![](https://github.com/mrklp/LpmiPred/blob/main/LpmiPred.png)

**The first step:**

Sequencing data of 37,820 small RNA libraries from six representative species (silkworm, fruit fly, nematode, zebrafish, mouse, and humans)  were downloaded and preprocessed, including quality control, removal of barcodes and adapters, and collapsing consistent reads. Finally, "reads_collapsed.fa" files were obtained for subsequent steps.

```shell
# removal of barcodes and adapters, quality control
python de_barcode_adapter.py  -q file.fastq -m  miRNAs.txt -s species  -o clean.fastq
# 
python collapse_reads_single.py -a clean.fasta -o collapse_reads.fasta
#
python collapse_reads_merger.py -f ./collapse_reads_single/ -o collapse_reads_all.fasta
```



```
# 在reads的引导下截取前体序列
python excise_precursors.py -a 20  -b reads_mappings.bwt -g genome.fa -c chr1 

# -a: 作为引导的reads的最低数量；
# -b: reads比对到基因组文件；
# -g: 基因组文件；
# -c: 待截取的染色体；
```