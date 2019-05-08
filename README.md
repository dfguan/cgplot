# cgplot
a script for viewing sequence aligments and read-depth

## Dependency
1. python3 

## Quick Start 

```
python cgplot.py -c 3L:23000000-28000000 -q000050F_arrow_arrow,000057F_arrow_arrow,000021F_arrow_arrow tests/Anoph_coluzzii2chrom.paf tests/PB.cov.wig
```
This will generate a image showing how the query contigs 50F, 57F and 21F are mapped to 23M-28M region of chromsome 3L, and also read depth plot for the contigs are shown on the bottom. 
![plot.png](https://github.com/dfguan/cgplot/blob/master/tests/plot.png) 

## Synopsis

```
usage: cgplot.py [options] paf_file wig_file

Genome Comparison plot

positional arguments:
  paf_file              a paf file
  wig_file              a wig file

optional arguments:
  -h, --help            show this help message and exit
  -c                    chromsome region in chr:start-end or chr format
  -q                    query name(s) mapped into the chromsome region,add comma to join multiple query names, support a maximum of 5 contigs
  -l                    minimum mapped length
  -o                    output file
  -t                    figure title
  --version             show program's version number and exit

```





 



