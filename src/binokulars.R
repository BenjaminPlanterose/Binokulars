# Call
# Rscript --vanilla BinomL.R --regions regions.txt --window.length 200 --offset 100 --epiread.file sd_epiread --plot F

# Read arguments
options = commandArgs(trailingOnly = TRUE)
index_names = startsWith(options, '--')
if(sum(index_names) != length(index_names)/2) stop('All arguments need to be specified like Rscript --vanilla BinomL.R --regions regions.txt --window.length 200 --offset 100 --epiread.file sd_epiread --plot F')
arguments = options[!index_names]
names(arguments) = options[index_names]
print(arguments)

# Store arguments as variables
regions = fread(arguments['--regions'])
filenames = arguments['--epiread.file']
window.length = as.numeric(arguments['--window.length'])
offset = as.numeric(arguments['--offset'])
plot_indicator = as.logical(arguments['--plot'])

# Load dependencies
library(data.table)
library(parallel)


# Define functions
# Binomial likelihood approach for ASM calling

s2c_bin <- function(seq)
{
  seq = strsplit(seq, split = '')[[1]]
  seq[seq == 'C'] = 1
  seq[seq == 'T'] = 0
  return(as.integer(seq))
}


# Better from old epiread format rather than .bam
# Single Fragment Epireads 
# chr19  -  3083513  CCCCCCC        3083495  ATT  3083513  CCCCCCC      3083495  ATT
# chr19  -  3083545  CCTCCCCCCT     .        .    3083527  CCCCTCCCCCC  3083523  AG

# Reg: chr9:124987720-124990776
read.location <- function(region)
{
  region = 'chr9:124987720-124990776'
  Chr = strsplit(region, split = ':')[[1]][1]
  region = strsplit(region, split = ':')[[1]][2]
  Start = as.integer(strsplit(region, split = '-')[[1]][1])
  End = as.integer(strsplit(region, split = '-')[[1]][2])
  command <- paste("grep ", "^", Chr, " ", file, " | ", "awk ", "'", "$3>", Start, " && ",
                   "$3<", End, "'", sep = "")
  fread(cmd = command, sep = "\t")
}

region_list = mclapply(X = regions, FUN = read.location)

# Window.length + offset
region = data.frame(V1 = c('chr19', 'chr19'), V2 = c('-', '-'), V3 = c('3083513', '3083545'), V4 = c('CCCCCCC', 'CCTCCCCCCT'), 
                    V5 = c('3083495', '.'), V6 = c('ATT', '.'), V7 = c('3083513', '3083527'), V8 = c('CCCCCCC', 'CCCCTCCCCCC'),
                    V9 = c('3083495', '3083523'), V10 = c('ATT', 'AG'))
setA = sapply(region$V4, s2c_bin)
setB = sapply(region$V8, s2c_bin)



# Model 1

deltaL = function(R_list)
{
  L_model2(R_list) - L_model1(R_list)
}

L_model1 = function(R_list)
{
  likelihood <- function(seq, p)
  {
    dbinom(x = sum(seq), size = length(seq), prob = p, log = T)
  }
  
  p = mean(unlist(R_list))
  sum(sapply(R_list, likelihood, p))
}

L_model2 = function(R_list)
{
  p = mean(sapply(R_list, function(x) mean(x)> 0.5))
  n = length(R_list)
  dbinom(x = n*p, size = n, prob = p, log = T)
}


set0 = list(c(1,1,1,1), c(1,1), c(1,1,1,1,1,1,1,1,1,1), c(1,1,1))
set1 = list(c(1,1,1,1), c(1,1), c(1,1,1,1,1,1,1,1,1,0), c(1,1,1))
set2 = list(c(1,0,1,1,0,0,1,1), c(0,1,1,1,1,0,0,0,1,1), c(1,1,0,0), c(1,0,0,0), c(1,0,0,0,0), c(1,0,0,1,1,0,1,1))
set3 = list(c(1,1,1,1), c(0,0,0,0,0,0,0), c(1,1,1,1,1,1,1,1,1,0), c(1,1,1), c(0,0,0,0,0), c(1,1,1,1,1,1), c(0,0,0,0))
set4 = list(c(1,1,1,1,1,1,1), c(0,0,0,0,0,0,0,1), c(1,1,1,1,1,1,1,1,1,0), c(1,1,1,0,1), c(0,0,0,0,0), c(1,1,1,1,1,1), c(0,0,0,0,1,1))
set5 = list(c(1,0,1,1,0,0,1,1), c(0,1,1,1,1,0,0,0,1,1), c(1,1,0,0), c(1,0,0,0,1,0), c(1,0,0,0,0,1,1), c(1,0,0,1,1,0,1,1))
set6 = list(c(0,0,0,0), c(0,1), c(0,0,0,0,0,0,0,0,0,0), c(0,1,0))
set7 = list(c(0,0,0,0), c(0,1), c(0,0,0,0,0,0,0,0,0,0), c(0,1,0), c(0,1,0,0,0,0,0), c(0,0,0,0,0,0,0,1), c(0,0,0,0,0,0,0,1))
set8 = list(c(0,0,0,0), c(0,0), c(0,0,0,0,0,0,0,0,0,0), c(0,0,0))
sets = list(set0, set1, set2, set3, set4, set5, set6, set7, set8)

plot(sapply(sets, L_model1), type = 'o')
plot(sapply(sets, L_model2), type = 'o')

plot(0:8, -sapply(sets, L_model1), type = 'o')
lines(0:8, sapply(sets, deltaL), type = 'o', lty = 2)






