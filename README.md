# Binokulars
## Binomial likelihood function-based bootstrap hypothesis test for co-methylation within reads


#### Benjamin Planterose Jiménez<sup>1</sup>, Brontë Kolar<sup>1</sup>, Manfred Kayser<sup>1</sup>, Athina Vidaki<sup>1</sup>

<sup>1</sup> Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands.


## Requirements

    Operating system: tested on Ubuntu 18.04.6 LTS
    R: tested on R version 4.1.2 (2021-11-01) -- "Bird Hippie"


## Installation

To install R dependencies, run the following R commands:

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

Also, make sure the right permissions are given to binokulars

```bash
chmod 777 binokulars
```

Finally, export path in .bashrc. Verify this step by running

```bash
binokulars --h
```

## Example

An example regions.txt file is given under /example/. To run 

```bash
binokulars -t absolute_path_to_regions.txt -i absolute_path_to_folder_CHR -l 200 -N 1000 -f 500 -R 4 -o test_results -c 4
```


Please contact me at b.planterosejimenez at erasmusmc.nl for any questions or issues concerning the scripts.

### References and Supporting Information
B. Planterose *et al* (**2022**). Identifying jointly regulated CpGs in the human methylome. *TBD*





