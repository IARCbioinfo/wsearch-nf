#!/usr/bin/env bash

chimera="http://drive5.com/uchime/rdp_gold.fa"
silva_seq="https://data.qiime2.org/2021.2/common/silva-138-99-seqs.qza"
silva_tax="https://data.qiime2.org/2021.2/common/silva-138-99-tax.qza"

# Here the script to create the database from the taxonomy and sequence file
# https://github.com/torognes/vsearch/issues/438

wget $silva_seq
wget $silva_tax

7za e silva-*-seqs.qza */data/*.fasta
7za e silva-*-tax.qza */data/*.tsv

awk 'NR>1' taxonomy.tsv \
  | sed 's/; /,/g' | sed 's/__/:/g' | sed 's/ /_/g' | sed 's/\t/;tax=/' | sed 's/^/>/' \
  > txx

awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' dna-sequences.fasta \
  | awk 'NR%2==0' \
  | paste -d'\n' txx - \
  | gzip > SILVA138_RESCRIPt.fasta.gz

rm txx


if [[ $chimera =~ ^http* ]]
then
        wget $chimera
        chimera=$(basename $chimera)
elif [[ ! -f "$chimera" ]]
then
     printf '%s\n' "$chimera not found !" >&2
     exit 1
else
     printf '%s\n' "found file $chimera" >&1
fi

