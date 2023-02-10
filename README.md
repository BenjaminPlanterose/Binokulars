# Binokulars
## Binomial likelihood function-based bootstrap hypothesis test for co-methylation within reads


#### Benjamin Planterose Jiménez, Brontë Kolar, Manfred Kayser, Athina Vidaki
Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands.

## Requirements

    Operating system: tested on Ubuntu 18.04.6 LTS
    R: tested on R version 4.1.2 (2021-11-01) -- "Bird Hippie" & R version 4.2.2 Patched (2022-11-10 r83330) -- "Innocent and Trusting"

For issues/questions on Binokulars, either [report an issue](https://github.com/BenjaminPlanterose/Binokulars/issues) or simply contact me via b.planterosejimenez at erasmusmc.nl

## Dependencies

#### [R](https://cran.r-project.org/)

To install R, click on the hyperlink above and follow instructions. To download dependency R-packages, do as follows. Open R by running:
```bash
R
```
And run the following R commands:

```r
# To install
install.packages('data.table')
install.packages('parallel')
install.packages('MASS')

# To verify that installation was succesful
library(data.table)
library(parallel)
library(MASS)
```

#### Binokulars

Obtain manual page by running:
```bash
cd <path_to_binokulars>/src/
bash binokulars --h
# Usage: binokulars -t FILE -i FILE [-l INT -N INT -R INT -o CHAR -c INT -f INT]
#
#   -t           A file containing target chromosomic locations (one per row) in the following format chr1:1234-3456.
#   -i           An input directory that includes sorted methylated/unmethylated cytosine counts per chromosome (see Github tutorial for details)."
#   -l           Bin length. Default value is 200 bp.
#   -N           Number of iterations for the permutation test; default value is 1000.
#   -R           Pseudo-random number generator seed; default value is 1.
#   -o           output directory; default value is results.
#   -c           Number of cores to employ; default value is 1.
#   -f           Flank length added left and right to each region; default value is 500.
```

## Test run

Under ```/test_run/``` you may find example data: 

* ```im_regions.txt``` - List of intermediately-methylated (im) regions
* ```/CHR/``` directory - includes a file per chromosome (for the sake of debugging, solely chr12 is provided).
	* ```chr12``` - Each line corresponds to a read pair from a paired-end WGBS experiment and contains the following fields:
		* (chr): on which chromosome the given read was aligned to.
		* (start): Start position of the read.
		* (C_Fw): Number of methylated cytosines in the forward read (i.e. paired-end WGBS sequencing).
		* (T_Fw): Number of unmethylated cytosines in the forward read (i.e. paired-end WGBS sequencing).
		* (C_Rv): Number of methylated cytosines in the reverse read (i.e. paired-end WGBS sequencing).
		* (T_Rv): Number of unmethylated cytosines in the reverse read (i.e. paired-end WGBS sequencing).

You may inspect the content of ```chr12``` file by running:
```bash
head <path_to_binokulars>/test_run/CHR/chr12
# (chr)	(start) (C_Fw)	(T_Fw)	(C_Rv)	(T_Rv)
# chr12	500004	1	0	0	0
# chr12	500004	1	0	2	0
# chr12	500004	1	1	0	0
# chr12	500004	1	1	0	0
# chr12	500004	1	2	0	0
# chr12	500004	2	1	0	0

```

To run the example:

```bash
bash <path_to_binokulars>/src/binokulars --h -t <path_to_binokulars>/test_run/im_regions.txt \
-i <path_to_binokulars>/test_run/CHR -l 200 -N 1000 -f 500 -R 4 -o test_results -c 1
```

## Run on your own data

In order to obtain the ```/CHR/``` directory, a substantial amount of parsing is required. You may find at ```src/opt/``` an example script on how to parse into the right format with [BISCUIT](https://huishenlab.github.io/biscuit/) and several bash commands, starting from fastq files.

In any case, if you plan to run Binokulars genome-wide, please visit our [JRC_seeker repository](https://github.com/BenjaminPlanterose/JRC_seeker).



