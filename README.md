# LpmiPred

```shell
# 在reads的引导下截取前体序列
python excise_precursors.py -a 20  -b reads_mappings.bwt -g genome.fa -c chr1 

# -a: 作为引导的reads的最低数量；
# -b: reads比对到基因组文件；
# -g: 基因组文件；
# -c: 待截取的染色体；
```
