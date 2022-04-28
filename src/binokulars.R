# Rscript --vanilla $MY_PATH/binokulars.R --target $target_file --epiread_chr_dir $input_dir 
# --bin_length $bin_length --N_iter $N_iter --seed $SEED --output_dir $output_dir 
# --Lflank $Lflank --nCores $nCores


# Read arguments
options = commandArgs(trailingOnly = TRUE)
index_names = startsWith(options, '--')
arguments = options[!index_names]
names(arguments) = options[index_names]
target_file = arguments['--target']
epiread_chr_dir = arguments['--epiread_chr_dir']
bin_length = as.numeric(arguments['--bin_length'])
N_iter = as.numeric(arguments['--N_iter'])
seed = as.integer(arguments['--seed'])
nCores = as.numeric(arguments['--nCores'])
Lflank = as.numeric(arguments['--Lflank'])
res_folder = arguments['--output_dir']
res_dir = paste(getwd(), res_folder, sep = '/')


setwd(res_dir)
capture.output(
print(paste('target_file = ', target_file, '; epiread_chr_dir = ', epiread_chr_dir, '; bin_length = ', bin_length, 
            '; N_iter = ', N_iter, '; seed = ', seed, 
            '; res_dir = ', res_dir, '; Lflank = ', Lflank, '; nCores = ', nCores, sep = '')), file = 'params.txt')

# Load dependencies
library(data.table)
library(parallel)
library(MASS)

permute <- function(P, L)
{
  M = vapply(1:length(P), function(i) rbinom(n = 1, size = L[i], prob = P[i]), as.integer(1L))
  -sum(stats::dbinom(x = M, size = L, prob = P, log = T))
}

permutation_test <- function(X, breaks, N_iter)
{
  zero_count = sapply(1:nrow(X), function(x) sum(as.numeric(X[x, c('V3', 'V4', 'V5', 'V6'), ]) == 0))
  X = X[zero_count != 4,]
  X$bins = cut(X$V2, breaks = breaks, labels = F)
  M_fw = aggregate(V3 ~ bins, data = X, FUN = sum)$V3
  n_fw = aggregate(V3+V4 ~ bins, data = X, FUN = sum)$`V3 + V4`
  M_rv = aggregate(V5 ~ bins, data = X, FUN = sum)$V5
  n_rv = aggregate(V5+V6 ~ bins, data = X, FUN = sum)$`V5 + V6`
  p = (M_fw + M_rv)/(n_fw + n_rv)
  l_bins = table(X$bins)
  p_vec = unlist(sapply(1:length(l_bins), function(i) rep(p[i], l_bins[i])))
  l_fw = X$V3+X$V4
  l_rv = X$V5+X$V6
  
  # Combine Fw and Rv in a single vector
  L = c(l_fw, l_rv)
  P = c(p_vec, p_vec)
  
  # Compute logL under H0
  logL_H0 = vapply(1:N_iter, function(i) permute(P, L), numeric(1))
  
  # Compute logL on data
  #logL = permute(P, L) # NULL MODE
  logL = -sum(stats::dbinom(x = c(X$V3, X$V5), size = L, prob = P, log = T)) # test
  if(logL == Inf) # Fwd methylation of 1, Rv 11110. Probability 0, likelihood infinity.
  {
    return(NA)
  }
  
  # Estimate p-values
  pval = 1 - mean(logL > logL_H0)
  if(pval == 0)
  {
    params = fitdistr(logL_H0, 'Gamma')
    pval = pgamma(q = logL, shape = params$estimate['shape'], rate = params$estimate['rate'], lower.tail = F)
  }
  return(pval)
}

process_region <- function(input_dir, regions_i, Chr_i, Start_i, End_i, bin_length, N_iter)
{
  #print(regions_i)
  setwd(input_dir)
  command <- paste("awk ", "'", "$2>=", Start_i, " && ", "$2<", End_i, "' ", Chr_i, sep = "")
  sub_data = fread(cmd = command, sep = "\t")
  breaks = seq(Start_i-1, End_i, bin_length)
  permutation_test(sub_data, breaks, N_iter)
}



# Read target_file
regions = fread(input = target_file, header = F)$V1
Chr = sapply(strsplit(regions, split = ':'), function(x) x[1])
pos = sapply(strsplit(regions, split = ':'), function(x) x[2])
Start = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[1])) - Lflank # Add flanks
End = as.integer(sapply(strsplit(pos, split = '-'), function(x) x[2])) + Lflank # Add flanks

# Run statistics
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)
pval_vec = mclapply(1:length(regions), function(i) tryCatch(process_region(epiread_chr_dir, regions[i], Chr[i], 
                                                                  Start[i], End[i], bin_length, N_iter), error = function(x) NA), 
                    mc.cores = nCores)
pval_vec = unlist(pval_vec)
names(pval_vec) = regions
write.table(x = pval_vec, file = 'p_values.txt', quote = F, col.names = F, sep = '\t')
