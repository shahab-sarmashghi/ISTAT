# ISTAT
ISTAT is a command-line tool for computing the significance of overlap between two sets of intervals. The paper where we have described the methods and tested ISTAT on multiple examples of simualted and biological datasets is available online:
  - [Sarmashghi, S., & Bafna, V. (2019). Computing the Statistical Significance of Overlap between Genome Annotations with iStat. Cell systems. https://doi.org/10.1016/j.cels.2019.05.006][1]

ISTAT is written in C++. You need a modern C/C++ compiler and CMake 2.8+ installed on your machine before proceeding with the installation and building the project.

Installation
------------
1. Clone the github repository by running (or you can download and unzip the repo)
```
    git clone https://github.com/shahab-sarmashghi/ISTAT.git
```
2. Change to the ISTAT directory and make an empty directory `build`
3. Change to `build` directory and build the project by running
```
    cmake ..
    make && make install
```
4. The executable is created and placed in `ISTAT/bin`. Change to `ISTAT/bin` and run `./istat` to get discription and usage information.
5. (Optional) add `ISTAT/bin` to the system path. 

Using ISTAT
------------
ISTAT gets 6 input arguments:
```
istat   QUERY_INTERVALS_FILE   REFERENCE_INTERVALS_FILE   CHROMOSOME_SIZES_FILE   SCALING_FACTOR   OUTPUT_FILE   THREADING_MODE
```
The intervals in `QUERY_INTERVALS_FILE` and `REFERENCE_INTERVALS_FILE` should be in bed format, and the intervals should also be sorted based on their start coordinate. An example of two sets of intervals and chromosome sizes are provided in `ISTAT/data` directory. 

A larger than 1 `SCALING_FACTOR` can be given as the input to scale down the running time and memory usage for the case that a large number of intervals are considered (set `SCALING_FACTOR` to 1 for no scaling). Both running time and memory usage scale linearly with the scaling factors. For example, if it takes 10 minutes to run using scaling factor of 1000, it will take about 100 minutes if the scaling factor is reduced to 100. We recommend starting with larger scaling factors (10000 or 1000) to get an estimate of how much time and memory are spent if a smaller scaling factor is used. Based on our experiments, the results are accurate when the scaling factor is the same order of smallest intervals you have in your intervals.

`THREADING_MODE` argument can be either `s` for no parallelization, or `p` to use all available cores. P-value is computed using two methods, DP (dynamic programming) and PB (Poisson binomial), and the summary of results is printed to the standard output. A complete list of p-values (for all possible overlaps) is written to `OUTPUT_FILE`.

For example, to run ISTAT on the intervals provided in `ISTAT/data` with scaling factor of 1000 and using parallelization:
```
istat hirt_chr.txt tcga_chr.txt  hg19_chrom_sizes.txt 1000 hirt-1000-pval.txt p
```
and the following is printed to your standard output:
```
scale is 1000
Number of query intervals (n) is 116
Number of reference intervals (m) is 101
Overlap is 54
eta is 0.00102011
DP p-value is 8.66247e-06
PB p-value is 2.63722e-10
```
Here, `54` out of `101` reference intervals are overlapped by query intervals, and the computed p-values using DP and PB methods for the observed overlap are reported. Further, all possible p-values for overlap ranging from `0` to `m` (number of reference intervals) is written to `hirt-1000-pval.txt`.

[1]: https://www.cell.com/cell-systems/pdf/S2405-4712(19)30187-5.pdf