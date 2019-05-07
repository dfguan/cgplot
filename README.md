# cgplot
a script for viewing sequence aligments and read-depth

## Dependency
1. python3 

## Quick Start 

```
python cgplot.py -c 3L:23000000-28000000 
-q000050F_arrow_arrow,000021F_arrow_arrow 
tests/Anoph_coluzzii2chrom.paf tests/PB.cov.wig
```
This will generate a [plot.png](https://github.com/dfguan/cgplot/tree/master/tests/plot.png) file showing how the query contigs 50F and 21F are mapped to 23M-28M region of chromsome 3L, and also read depth plot are shown in the bottom. 

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
  -q                    query name(s) mapped into the chromsome region,add comma to join multiple query names
  -l                    minimum mapped length
  -o                    output file
  -t                    figure title
  --version             show program's version number and exit

```





 



