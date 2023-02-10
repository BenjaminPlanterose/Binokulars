# Binokulars
## Binomial likelihood function-based bootstrap hypothesis test for co-methylation within reads


#### Benjamin Planterose Jiménez, Brontë Kolar, Manfred Kayser, Athina Vidaki
Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands.


## Requirements

    Operating system: tested on Ubuntu 18.04.6 LTS
    R: tested on R version 4.1.2 (2021-11-01) -- "Bird Hippie" & R version 4.2.2 Patched (2022-11-10 r83330) -- "Innocent and Trusting"


## Dependencies

#### [R](https://cran.r-project.org/)

To install R, click on the hyperlink above and follow instructions. To download dependency R-packages, run the following:

Open R by running:
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
```

## Test run

Under ```/test_run/``` you may find example data. This consists of a folder /CHR/ that includes a file per each chromosome, containing the following fields (without header):


```
head <path_to_binokulars>/test_run/CHR/chr12
# (chr)	(start) (C_Fw)	(T_Fw)	(C_Rv)	(T_Rv)
# chr12	500004	1	0	0	0
# chr12	500004	1	0	2	0
# chr12	500004	1	1	0	0
# chr12	500004	1	1	0	0
# chr12	500004	1	2	0	0
# chr12	500004	2	1	0	0

```




To run the example, a test file regions.txt is given under /example/. To run:

```bash
cd example
binokulars -t absolute_path_to_regions.txt -i absolute_path_to_folder_CHR -l 200 -N 1000 -f 500 -R 4 -o test_results -c 4
```


Please contact me at b.planterosejimenez at erasmusmc.nl for any questions or issues concerning the scripts.

