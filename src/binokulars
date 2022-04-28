#!/bin/bash

# Manual page
if [ "$1" == "-h" ] || [ "$1" == "--h" ] || [ "$1" == "--help" ] || [ "$1" == "-help" ]; then
  echo "Usage: $0 -t FILE -i FILE [-l INT -N INT -R INT -o CHAR -c INT -f INT]" >&2
  echo
  echo "   -t           A file containing target chromosomic locations (one per row) in the following format chr1:1234-3456."
  echo "   -i           An input directory that includes sorted single_fragment.epiread files per chromosome (see Biscuit manual page for more details: https://huishenlab.github.io/biscuit/epiread_format/#single-fragment-epireads)."
  echo "   -l           Bin length. Default value is 200 bp."
  echo "   -N           Number of iterations for the permutation test; default value is 1000."
  echo "   -R           Pseudo-random number generator seed; default value is 1."
  echo "   -o           output directory; default value is results."
  echo "   -c           Number of cores to employ; default value is 1."
  echo "   -f           Flank length added left and right to each region; default value is 500."
  echo  
  exit 0
fi

# Parse arguments
while getopts ":t:i:l:N:R:o:c:f:" opt; do
  case $opt in
    t) target_file="$OPTARG"
    ;;
    i) input_dir="$OPTARG"
    ;;
    l) bin_length="$OPTARG"
    ;;
    N) N_iter="$OPTARG"
    ;;
    R) SEED="$OPTARG"
    ;;
    o) output_dir="$OPTARG"
    ;;
    c) nCores="$OPTARG"
    ;;
    f) Lflank="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) echo "Option $opt needs a valid argument"
    exit 1
    ;;
  esac
done

# Make sure all mandatory arguments have been input
if [ -z "$N_iter" ]; then
  N_iter=1000
fi

if [ -z "$bin_length" ]; then
  bin_length=200
fi

if [ -z "$SEED" ]; then
  SEED=1
fi

if [ -z "$output_dir" ]; then
  output_dir="results"
fi

if [ -z "$nCores" ]; then
  nCores=1
fi

if [ -z "$Lflank" ]; then
  Lflank=500
fi



if [ -z "$target_file" ]; then
  printf "***********************************\n"
      printf "* Error: target_file not supplied.*\n"
      printf "***********************************\n"
      exit 1
fi

if [ -z "$input_dir" ]; then
  printf "**********************************\n"
      printf "* Error: input_dir not supplied.*\n"
      printf "**********************************\n"
      exit 1
fi

# Find path to script (also, it is where the Rscript is located)
MY_PATH=$(dirname "$0")
echo "$MY_PATH"

# Print selected arguments
printf "Argument input_dir is %s\n" "$input_dir"
printf "Argument target_file is %s\n" "$target_file"
printf "Argument bin_length is %s\n" "$bin_length"
printf "Argument N_iter is %s\n" "$N_iter"
printf "Argument SEED is %s\n" "$SEED"
printf "Argument output_dir is %s\n" "$output_dir"
printf "Argument nCores is %s\n" "$nCores"
printf "Argument Lflank is %s\n" "$Lflank"
printf "Argument MY_PATH is %s\n" "$MY_PATH"


# Run R script
mkdir "$output_dir"
cp $target_file $output_dir
Rscript --vanilla $MY_PATH/binokulars.R --target $target_file --epiread_chr_dir $input_dir --bin_length $bin_length --N_iter $N_iter --seed $SEED --output_dir $output_dir --Lflank $Lflank --nCores $nCores

