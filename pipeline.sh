# Pipeline to add BOLD and GenBank sequences to lab mitogenomes

# -------------------
# Programs to install
# -------------------

# python, biopython, requests
# BLAST
# r-base, r-getopt, r-dplyr
# MAFFT
# mPTP
# RAxML-NG
# treeshrink

# -------------------
# Github scripts uesd
# -------------------

# https://github.com/Aileen-S/voglerlab
# bold_cli.R
# get_lab_gbs.py 
# extract_lab_genes.py
# filter_fasta.py
# filter_tips.R
# partitions.py
# process_blast.py
# ptp_filter_output.py
# remove_dups.py
# rename_fasta.py
# ry_code.py
# supermatrix_count.py


# https://github.com/tjcreedy/biotools
# findframe.py
# translate.py
# backtranslate.py

# https://github.com/nylander/catfasta2phyml
# catfasta2phyml.pl

# https://github.com/nylander/fastagap.git
# degap_fasta_alignment.pl

# -------------------
# Files to have ready
# -------------------

# ID list of lab sequences for chosen taxa
# ID list of lab sequences for outgroup
# Comma separted list of outgroup taxa as they are named in the supermatrix (for mPTP)

# Fastas for:
  # Profile: unaligned nucleotides for COX1a and COX1b for BLAST
  # Profile: nucleotide alignments
  # Profile: protein alignments

  # Fasta file name format:
  #   12S.fasta 16S.fasta ATP6.fasta ATP8.fasta COX1a.fasta COX1b.fasta COX2.fasta COX3.fasta CYTB.fasta
  #   ND1.fasta ND2.fasta ND3.fasta ND4.fasta ND4L.fasta ND5.fasta ND6.fasta
  #   18S.fasta 28S.fasta AK.fasta CAD.fasta EF1A.fasta H3.fasta RNApol.fasta Wg.fasta




# Remember to modify each command to add the path to each file and each github script

cat << p
------------
Define taxon
------------
p

# Eg:
taxon=Dytiscidae

mkdir $taxon
cd $taxon


cat << p
---------------
Lab mitogenomes
---------------
p

# Get lab mitogenome sequences from Valentine database
# First write input file with list of required mitogenome IDs - below called 'lab_ids.txt'
# Remember to include outgroup
# Add latest MMG database version to command, eg 'gbmaster_2024-04-20'

mkdir mmg
cd mmg

python3 github/voglerlab/get_lab_gbs.py -i lab_ids.txt -o $taxon.gb -v version 
python3 github/voglerlab/extract_lab_genes.py -g $taxon.gb


cat << p

---------------
Get NCBI TaxIDs
---------------
p

mkdir blast
cd blast

# Define taxon, eg:
taxon=Dytiscidae

# Give path to BLAST database, eg. on Crop Diversity:
export BLASTDB=/mnt/shared/datasets/databases/ncbi/
# or on Valentine:
export BLASTDB=/mbl/share/workspaces/groups/database/nt-2024-05-29_blastdb

# Add your personal API key to speed up GenBank search
export NCBI_API_KEY=

# Get NCBI taxonomy ID for chosen taxon with script from BLAST package
# If multiple taxa are required, run this process for each then concatenate the output ID lists
output=$(get_species_taxids.sh -n $taxon)
txid=$(echo "$output" | grep 'Taxid' | awk '{print $3}')
echo "Taxon ID: $txid"

# Get list of NCBI taxon IDs for chosen taxon
echo "Getting NCBI txid list"
get_species_taxids.sh -t $txid > $taxon.txids
# Get outgroup NCBI taxids from metadata and append to $taxon.txids


cat << p

------------
BLAST Search
------------
p

gene=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' ../config.genes)

# BLAST search for each gene, using profile fasta as query
# Output gives NCBI taxon ID, GenBank accession and match start and stop positions
# Run processes in parallel if desired using a slurm array job.
blastn -db nt -query $gene.fasta -taxidlist $taxon.txids -out $gene.blast -max_target_seqs 100000000 -outfmt '6 staxids sacc sstart send' -evalue 1e-5 -num_threads 4 

# For each gene, get coordinates of longest sequence for each TXID
python3 github/voglerlab/process_blast.py -i $gene.blast -o $gene.out -e email

# For each gene, extract sequences from BLAST database
> "$gene.fasta"
while IFS= read -r line
do
# Read sequence coordinates from $gene.out file
  accession=$(echo "$line" | awk '{print $1}')
  start=$(echo "$line" | awk '{print $2}')
  end=$(echo "$line" | awk '{print $3}')

# Run blastdbcmd and write fasta file
  blastdbcmd -db nt -entry "$accession" -range "$start-$end" -target_only -outfmt %f >> "$gene.fasta"
done < "$gene.out"

# Strip fasta IDs to leave only accession
for file in *fasta
do
  sed -i 's/^\(>[^.]*\).*/\1/' $file
done

cd ..


cat << p

--------------
GenBank Search
--------------
p


file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' ../config.blast)
echo $file

# Direct to BLAST database
export NCBI_API_KEY=

python3 ~/apps/github/voglerlab/get_genbanks.py -f $file.txids -r txid -e mixedupvoyage@gmail.com -t Coleoptera


cat << p

-----------
BOLD search
-----------
p


mkdir bold
cd bold

# Retrive BOLD sequences and metadata for taxon
python3 github/voglerlab/search_bold.py -t $taxon -o bold_metadata.tsv
# Multiple taxa can be seperated by a comma
python3 github/voglerlab/search_bold.py -t taxon1,taxon2 -o bold_metadata.tsv


# Filter output to get one barcode per BIN and write to fasta
Rscript github/voglerlab/bold_cli.R -c bold_metadata.tsv -b -g

cd ..


cat << p

-----------------
Merge gene fastas
-----------------
p

mkdir frame
cd frame

# Define the list of possible file names
file_names=("12S.fasta" "16S.fasta" "18S.fasta" "28S.fasta" "AK.fasta" "ATP6.fasta" "ATP8.fasta" "CAD.fasta" "COX1.fasta" "COX1a.fasta" "COX1b.fasta" "COX2.fasta" "COX3.fasta" "CYTB.fasta" "EF1A.fasta" "H3.fasta" "ND1.fasta" "ND2.fasta" "ND3.fasta" "ND4.fasta" "ND4L.fasta" "ND5.fasta" "ND6.fasta" "RNApol.fasta" "Wg.fasta")

# List of directories to check
directories=(
  "../mmg"
  "../blast"
  "../bold"
)

# Merge fastas for each gene
for f in "${file_names[@]}"; do
  output="$f"
  : > "$output"  # Create or truncate the output file

  for dir in "${directories[@]}"; do
    if [ -e "$dir/$f" ]; then
      cat "$dir/$f" >> "$output"
    fi
  done
done

# Separate into mitochondrial protein coding, nuclear protein coding and RNA
mkdir m1_raw n1_raw r1_raw
mv 1*.fasta 2*.fasta r1_raw
mv COX*fasta CYTB.fasta ATP*fasta ND*fasta m1_raw
mv RNApol.fasta Wg.fasta EF1A.fasta H3.fasta CAD.fasta AK.fasta n1_raw


cat << p

--------------
COX1 sequences
--------------
p
# Merge COX1a and COX1b with complete COX# Mitogenome frame tags
sed -i -E "s/;frame=[0-9]*(;$)?//" m1_raw/$fgene.fasta
mafft --add m1_raw/$fgene.fasta --maxiterate 2 --adjustdirection --thread 8 ~/ascott/profiles/0_NT_profiles/$base >  m2_align/$gene.fasta
sed "/>/! s/-//g" m2_align/$gene.fasta > m3_nogaps/$gene.fasta
~/apps/github/biotools/findframe.py -r ~/ascott/profiles/0_AA_profiles/$gene.fasta -t 5 -s -u < m3_nogaps/$gene.fasta > m4_frame/$gene.fasta

1 sequences. Seperate manually after alignment
mkdir COI
mv m1_raw/COX1* COI
cat COI/COX1* > COI/COX1all.fasta
cp COI/COX1all.fasta m1_raw/COX1.fasta

# Alternatively split complete COX1 sequences into COX1a and COX1b regions at this stage


cat << p

-----------------
Get reading frame
-----------------
p

mkdir m2_align n2_align m3_nogaps n3_nogaps m4_frame n4_frame

# Mitogenome frame tags
for file in m1_raw/*
do
  base=$(basename $file)
  echo $base
  sed -i -E "s/;frame=[0-9]*(;$)?//" $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 8 profiles/0_NT_profiles/$base >  m2_align/$base
  sed "/>/! s/-//g" m2_align/$base > m3_nogaps/$base
  github/biotools/findframe.py -r profiles/0_AA_profiles/$base -t 5 -s -u < m3_nogaps/${base#*/} > m4_frame/${base#*/}
done

# Nuclear frame tags
for file in n1_raw/*
do
  base=$(basename $file)
  echo $base
  sed -i -E "s/;frame=[0-9]*(;$)?//" $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 8 profiles/0_NT_profiles/$base >  n2_align/$base
  sed "/>/! s/-//g" n2_align/$base > n3_nogaps/$base
  github/biotools/findframe.py -r profiles/0_AA_profiles/$base -t 5 -s -u < n3_nogaps/${base#*/} > n4_frame/${base#*/}
done

cd ..

cat << p

---------
Alignment
---------
p

mkdir alignment

# Copy fastas with reading frames to alignment directory
cp -r frame/m4_frame alignment/m1_raw
cp -r frame/n4_frame alignment/n1_raw
cp -r frame/r1_raw/ alignment

cd alignment

# Add directories for further processing
mkdir n2_aa n3_aaal n4_ntal m2_aa m3_aaal m4_ntal r4_ntal


cat << p

--------------------
Translate to Protein
--------------------
p

# Mitochondrial genes
for file in m1_raw/*
do
  echo 'translate' $file
   github/biotools/translate.py 5 < $file > m2_aa/${file#*/}
done

# Nuclear genes
for file in n1_raw/*
do
  echo 'translate' $file
   github/biotools/translate.py 1 < $file > n2_aa/${file#*/}
done


cat << p

----------------
Align to Profile
----------------
p

# Mitochondrial
for file in m2_aa/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 8 profiles/0_AA_profiles/${file#*/} >  m3_aaal/${file#*/}
done

# Nuclear
for file in n2_aa/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 8 profiles/0_AA_profiles/${file#*/} >  n3_aaal/${file#*/}
done

# RNA
for file in r1_raw/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 8 profiles/0_AA_profiles/${file#*/} >  r4_ntal/${file#*/}
done


# Remove profile

cd m3_aaal
for file in *
do
  echo $file
  mv $file ${file}_bak
  perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${file}_bak | grep -A 1 -P '^>(?!PROFILE::).*' > $file
  rm ${file}_bak
done
cd ..

cd n3_aaal
for file in *
do
  echo $file
  mv $file ${file}_bak
  perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${file}_bak | grep -A 1 -P '^>(?!PROFILE::).*' > $file
  rm ${file}_bak
done
cd ..

cd r4_ntal
for file in *
do
  echo $file
  mv $file ${file}_bak
  perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${file}_bak | grep -A 1 -P '^>(?!PROFILE::).*' > $file
  rm ${file}_bak
done
cd ..


cat << p

---------------------------------
Translate to Nucleotide Alignment
---------------------------------
p

# Mitochondrial
for file in m3_aaal/*
do
  echo $fileStothard Road
  github/biotools/backtranslate.py -i $file m1_raw/${file#*/} 5 > m4_ntal/${file#*/}
done

# Nuclear
for file in n3_aaal/*
do
  echo $file
  github/biotools/backtranslate.py -i $file n1_raw/${file#*/} 5 > n4_ntal/${file#*/}
done


###########################################################
# Manual check and cleaning of alignments                 #
# Split COX1 fasta into COX1a and COX1b before proceeding #
###########################################################


cat << p

---------------------------------
Prepare sequences for supermatrix
---------------------------------
p

mkdir ../supermatrix
mkdir ../supermatrix/aa
mkdir ../supermatrix/nt

cp m3_aaal/*fasta n3_aaal/*fasta r4_ntal/*fasta ../supermatrix/aa
cp m4_ntal/*fasta n4_ntal/*fasta r4_ntal/*fasta ../supermatrix/nt

cp m5_aa_cia/*fasta n5_aa_cia/*fasta r6_nt_cia/*fasta ../supermatrix/aa
cp m6_nt_cia/*fasta n6_nt_cia/*fasta r6_nt_cia/*fasta ../supermatrix/nt

cd ../supermatrix


# Remove frame tags
for file in aa/*
do
  sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in aa"

for file in nt/*
do
  sed -i -E "s/;frame=[0-9]*(;$)?//" $file
done
echo "Frame tags removed from files in nt"


# Remove _R_ from reversed sequences
for file in aa/*
do
   sed -i -E "s/^>_R_/>/" $file
done

for file in nt/*
do
   sed -i -E "s/^>_R_/>/" $file
done


cat << p

-----------------------
Rename with NCBI taxids
-----------------------
p

# Rename sequences with ncbi taxids

mkdir ids_aa ids_nt

for file in aa/*
do
  base=$(basename $file)
  python3 github/voglerlab/rename_fasta.py -i $file -c github/voglerlab/metadata/rename_list_ref_taxid.csv -l -r ids_aa/${base%fasta}csv -o aa/$base.r
  mv aa/${base#*/}.r aa/${base#*/}
done

for file in nt/*
do
  base=$(basename $file)
  python3 github/voglerlab/rename_fasta.py -i $file -c github/voglerlab/metadata/rename_list_ref_taxid.csv -l -r ids_nt/${base%fasta}csv -o nt/$base.r
  mv nt/${base#*/}.r nt/${base#*/}
done

# Write gene source CSVs
Rscript github/voglerlab/merge_ptp_csv.R -d ids_aa/ -o ids_aa.csv
Rscript github/voglerlab/merge_ptp_csv.R -d ids_nt/ -o ids_nt.csv

# Write nucleotide supermatrix
github/catfasta2phyml/catfasta2phyml.pl -c -fasta nt/* > 1_nt_supermatrix.fasta 2> 1_cfnt_partitions.txt

# Convert partition files to RAxML format (GTR model)
python3 github/voglerlab/partitions.py -i 1_cfnt_partitions.txt -t nt -o raxml -p 1

# Write gene representation CSV
python3 github/voglerlab/supermatrix_count.py -f 1_nt_supermatrix.fasta -p 1_cfnt_partitions.txt -o nt.csv

cd ..

cat << p

------------------------------------------
Initial phylogeny for species delimitation
------------------------------------------
p

mkdir mptp

# Check required memory/threads
raxml-ng --parse --msa supermatrix/1_nt_supermatrix.fasta --prefix mptp/search1 --model supermatrix/1_partitions_gene.txt

# Ran phylogeny with single starting tree
raxml-ng --search1 --msa supermatrix/1_nt_supermatrix.fasta --prefix mptp/search1 --model supermatrix/1_partitions_gene.txt --threads number_threads

# Run modelfinder at this stage as well?

cat << p

-------------------------
mPTP Species Delimitation
-------------------------
p

cd mptp

# Remove excessively long branches
run_treeshrink.py -t search1.raxml.bestTree

# Add taxonomy to tree/supermatrix at this stage to ease analysis of result
python3 github/voglerlab/rename_fasta.py -i ../supermatrix/1_nt_supermatrix.fasta -o ../supermatrix/1_tax_nt_supermatrix.fasta -c metadata/rename_taxid_taxonomy.csv
Rscript github/voglerlab/rename_tree.R -i search1.raxml_treeshrink/output.bestTree -o tax.treeshrink.tree -c metadata/rename_taxid_taxonomy.csv

# Filter supermatrix to keep taxa remaining after treeshirink
python3 github/voglerlab/filter_fasta.py -i ../../supermatrix/1_tax_nt_supermatrix.fasta -f supermatrix.fasta -s tax.treeshrink.tree -t tree

# Requires file with comma separated list of outgroup taxa
# Get outgroup
outgroup=$(awk 'NR==1' outgroup.txt)

# Calculate minbr value
echo "Getting minbr value"
mptp --minbr_auto supermatrix.fasta --tree_file tax.treeshrink.tree --output_file minbr  --outgroup $outgroup

# Copy minbr value from slurm output
minbr=$(awk 'END {print $NF}' minbr.txt)
echo "minbar value = " $minbr

# Run PTP
echo "Running mPTP"
mptp --ml --single --minbr $minbr --tree_file tax.treeshrink.tree --output_file mptp  --outgroup $outgroup

# Run MCMC sampling to determine statistical significance
echo "
Running MCMC"
mptp --mcmc 1000000 --mcmc_sample 1000 --mcmc_log 1000 --mcmc_runs 2 --multi --minbr $minbr --tree_file tax.treeshrink.tree --output_file mcmc
cd ..

# Identify taxon with most nucleotides for each PTP species
# Use -n and/or -k flags if necessary
# -n option retains all with species level identification by checking if the last string is all lowercase, assuming the delimiter is an underscore
# -k option accepts a text file with a list of taxa to keep (eg mitogenomes, named species)
python3 github/voglerlab/ptp_filter_output.py -i mptp.txt -s supermatrix.fasta -o ptp_chosen.txt -n


# Filter sequence files using PTP output
cd ../supermatrix

mkdir 2_aa_ptp 2_nt_ptp

for file in aa/*
do
  echo $file
  python3 github/voglerlab/filter_fasta.py -i $file -s ../mptp/ptp_chosen.txt -t list -f 2_aa_ptp/${file#*/}
  #Remove gap only columns
  echo 'Removing gap-only columns'
  github/fastagap/degap_fasta_alignment.pl 2_aa_ptp/${file#*/} > 2_aa_ptp/${file#*/}.g
  mv 2_aa_ptp/${file#*/}.g 2_aa_ptp/${file#*/}
done

for file in nt/*
do
  echo $file
  python3 github/voglerlab/filter_fasta.py -i $file -s ../mptp/ptp_chosen.txt -t list -f 2_nt_ptp/${file#*/}
  # Remove gap only columns
  echo 'Removing gap-only columns'
  github/fastagap/degap_fasta_alignment.pl 2_nt_ptp/${file#*/} > 2_nt_ptp/${file#*/}.g
  mv 2_nt_ptp/${file#*/}.g 2_nt_ptp/${file#*/}
done

# Write supermatrices
github/catfasta2phyml/catfasta2phyml.pl -c -fasta 2_aa_ptp/* > 2_aa_ptp_supermatrix.fasta 2> 2_cfaa_partitions.txt
github/catfasta2phyml/catfasta2phyml.pl -c -fasta 2_nt_ptp/* > 2_nt_ptp_supermatrix.fasta 2> 2_cfnt_partitions.txt
python3 github/voglerlab/ry_code.py -i 2_nt_ptp_supermatrix.fasta -o 2_ry_ptp_supermatrix.fasta

# Convert partition files to RAxML format
# Modify AA models if necessary
python3 github/voglerlab/partitions.py -i 2_cfaa_partitions.txt -t aa -o raxml -p 2_ptp
python3 github/voglerlab/partitions.py -i 2_cfnt_partitions.txt -t nt -o raxml -p 2_ptp

# Write gene representation CSV
python3 github/voglerlab/rename_fasta.py -i 2_nt_ptp_supermatrix.fasta -o 2_tax_nt_ptp_supermatrix.fasta -c github/voglerlab/metadata/rename_taxid_taxonomy.csv
python3 github/voglerlab/supermatrix_count.py -f 2_tax_nt_ptp_supermatrix.fasta -p 2_cfnt_partitions.txt -o 2_nt_ptp.csv


cat << p

-------------------------------
Choose taxa for constraint tree
-------------------------------
p

# Probably can skip this section, and instead use subtrees from the 13K tree as constraint


# # Decide on constraint taxa and write backbone.ids file
# mkdir 3_aa_con 3_nt_con

# for file in 2_aa_ptp/*
# do
#   echo $file
#   python3 github/voglerlab/filter_fasta.py -i $file -s backbone.ids -t list -f 3_aa_con/${file#*/}
#   #Remove gap only columns
#   echo 'Removing gap-only columns'
#   github/fastagap/degap_fasta_alignment.pl 3_aa_con/${file#*/} > 3_aa_con/${file#*/}.g
#   mv 3_aa_con/${file#*/}.g 3_aa_con/${file#*/}
# done

# for file in 2_nt_ptp/*
# do
#   echo $file
#   python3 github/voglerlab/filter_fasta.py -i $file -s backbone.ids -t list -f 3_nt_con/${file#*/}
#   # Remove gap only columns
#   echo 'Removing gap-only columns'
#   github/fastagap/degap_fasta_alignment.pl 3_nt_con/${file#*/} > 3_nt_con/${file#*/}.g
#   mv 3_nt_con/${file#*/}.g 3_nt_con/${file#*/}
# done

# # Write supermatrix
# github/catfasta2phyml/catfasta2phyml.pl -c -fasta 3_aa_con/* > 3_aa_con_supermatrix.fasta 2> 3_cfaa_partitions.txt
# github/catfasta2phyml/catfasta2phyml.pl -c -fasta 3_nt_con/* > 3_nt_con_supermatrix.fasta 2> 3_cfnt_partitions.txt

# python3 github/voglerlab/ry_code.py -i 3_nt_con_supermatrix.fasta -o 3_ry_con_supermatrix.fasta
# python3 github/voglerlab/ry_code.py -i 3_nt_con_supermatrix.fasta -o 3_bin_con_supermatrix.fasta -n

# # Convert partition files to RAxML format
# python3 github/voglerlab/partitions.py -i 3_cfaa_partitions.txt -t aa -o raxml -p 3_con
# python3 github/voglerlab/partitions.py -i 3_cfnt_partitions.txt -t nt -o raxml -p 3_con
# python3 github/voglerlab/partitions.py -i 3_cfnt_partitions.txt -t bin -o raxml -p 3_con


cd ..


cat << p

-------------------
Run constraint tree
-------------------
p


mkdir raxml
cd raxml

# Get random number seed for RAxML and save to slurm output
seed=$RANDOM
echo 'RAxML random number seed = '$seed

# Run RaxML
raxml-ng --msa supermatrix/3_nt_con_supermatrix.fasta --model supermatrix/3_con_partitions_codon.txt --prefix raxml/con_codon --threads 8 --seed $seed
raxml-ng --msa supermatrix/3_ry_con_supermatrix.fasta --model supermatrix/3_con_partitions_codon.txt --prefix raxml/con_codon --threads 8 --seed $seed
raxml-ng --msa supermatrix/3_bin_con_supermatrix.fasta --model supermatrix/3_con_partitions_bin_codon.txt --prefix raxml/con_codon --threads 8 --seed $seed
raxml-ng --msa supermatrix/3_aa_con_supermatrix.fasta --model supermatrix/3_con_partitions_aa.txt --prefix raxml/con_aa --threads 8 --seed $seed


cat << p

--------
RAxML full tree
--------
p

# Or run first with onl taxa with barcodes?

seed=$RANDOM
echo "Random Seed:" $RANDOM

mkdir ml_nt_codon 
raxml-ng --msa ../supermatrix/2_nt_ptp_supermatrix.fasta --model ../supermatrix/ptp_partitions_codon.txt --prefix ml_nt_codon --threads 8 --seed $seed --force perf_threads

mkdir ml_aa
raxml-ng --msa ../supermatrix/1_aa_ptp_supermatrixconstrai.fasta --model ../supermatrix/ptp_partitions_aa.txt --prefix ml_aa --threads 8 --seed $seed --force perf_threads

# Make RY supermatrix
mkdir ml_ry_codon
raxml-ng --msa ../supermatrix/3_ry_ptp_supermatrix.fasta --model ../supermatrix/ptp_partitions_codon.txt --prefix ml_ry_codon --threads 8 --seed $seed --force perf_threads

# Remove long branches
outgroup=$(awk 'NR==1' outgroup.txt)
run_treeshrink.py -t input_tree -c -x outgroup_taxa_list -q 0.5 


# Check taxonomic indices

python3 ~/github/voglerlab/supermatrix_count.py -i ../supermatrix/2_nt_coi_ptp_ts_supermatrix.fasta -p ../supermatrix/2_nt_coi_partitions.txt -o ../supermatrix/2_nt_coi_ptp_ts.csv -t

Rscript	~/github/phylostuff/taxonomicindices.R	-p	2_con_aa.raxml.bestTree	-t	../tri_taxonomy.csv	-a	tri_2_con_aa.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	2_con_nt.raxml.bestTree	-t	../tri_taxonomy.csv	-a	tri_2_con_nt.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	2_con_ry.raxml.bestTree	-t	../tri_taxonomy.csv	-a	tri_2_con_ry.csv		-l	tribe

Rscript ~/github/voglerlab/taxind_cat.R * taxind.csv

Rscript	~/github/phylostuff/taxonomicindices.R	-p	aa_ck_coi.raxml_treeshrink/output.bestTree	-t	../tri_taxonomy.csv	-a	tri_aa_ck_coi.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	aa_con10.raxml_treeshrink/output.bestTree	-t	../tri_taxonomy.csv	-a	tri_aa_con10.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	aa_conmt2.raxml_treeshrink/output.bestTree	-t	../tri_taxonomy.csv	-a	tri_aa_conmt2.csv		-l	tribe

Rscript	~/github/phylostuff/taxonomicindices.R	-p	aa_ck_coi.raxml_treeshrink/output.bestTree	-t	../tri_taxonomy.csv	-a	tri_aa_ck_coi.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	aa_con10.raxml_treeshrink/output.bestTree	-t	../tri_taxonomy.csv	-a	tri_aa_con10.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	aa_conmt2.raxml_treeshrink/output.bestTree	-t	../tri_taxonomy.csv	-a	tri_aa_conmt2.csv		-l	tribe

Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.aa_ck_coi.raxml.bestTree	-t	../../rename_tree.csv	-a	aa_ck_coi.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.aa_con10.raxml.bestTree	-t	../../rename_tree.csv	-a	aa_con10.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.aa_conmt2.raxml.bestTree	-t	../../rename_tree.csv	-a	aa_conmt2.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.codon_ck_coi.raxml.bestTree	-t	../../rename_tree.csv	-a	codon_ck_coi.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.fcon_full.raxml.bestTree	-t	../../rename_tree.csv	-a	fcon_full_aa.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.full.raxml.bestTree	-t	../../rename_tree.csv	-a	full_aa.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.nt_ck_coi.raxml.bestTree	-t	../../rename_tree.csv	-a	nt_ck_coi.csv		-l	tribe
Rscript	~/github/phylostuff/taxonomicindices.R	-p	tax.ry_ck_coi.raxml.bestTree	-t	../../rename_tree.csv	-a	ry_ck_coi.csv		-l	tribe

for file in *tbe*support
do
  Rscript ~/github/phylostuff/taxonomicindices.R  -p $file -t ../../taxonomy250416.csv -a subfamily.${file%tbe.support}csv -l subfamily
done

Rscript ~/github/voglerlab/taxind_cat.R -d . -o subfamily_taxindices.csv



# Bootstrap

#!/bin/bash
#SBATCH --job-name=boot_codon
#SBATCH --cpus-per-task=6
#SBATCH --array=1-500
#SBATCH -p short

seed=$RANDOM
echo "Random Seed: $RANDOM_SEED"
3898
mkdir boot_mito_codon
raxml-ng --bootstrap --msa ../supermatrix/nt_con_mito_supermatrix.fasta --model ../supermatrix/con_mito_partitions_codon.txt --prefix boot_mito_codon/$SLURM_ARRAY_TASK_ID.boot --seed $seed --bs-trees 1 --threads 6 --force perf_threads

# Check for convergence
cat boot_mito_codon/*bootstraps > all_codon.bootstraps
raxml-ng --bsconverge --bs-trees all_codon.bootstraps --prefix boot_mito_codon --seed 2
# Add support values to tree
raxml-ng --support --tree ntcon_mito/ntcon_mito.bestTree --bs-trees all_codon.bootstraps --prefix ntcon_mito/boot
# Strict consensus
raxml-ng --consense STRICT  --tree all_codon.bootstraps --prefix ntcon_mito/strict_boot
# Majority rule consensus
raxml-ng --consense MR  --tree all_codon.bootstraps --prefix ntcon_mito/majority_boot

# Check for convergence
cat boot_mito_aa/*bootstraps > all_aa.bootstraps
raxml-ng --bsconverge --bs-trees all_aa.bootstraps --prefix boot_mito_aa --seed 2
# Add support values to tree
raxml-ng --support --tree aacon_mito/aacon_mito.bestTree --bs-trees all_aa.bootstraps --prefix aacon_mito/boot
#TBE
raxml-ng --support --tree aacon_mito/aacon_mito.bestTree --bs-trees all_aa.bootstraps --prefix aacon_mito/boot --bs-metric tbe

# Strict consensus
raxml-ng --consense STRICT  --tree all_aa.bootstraps --prefix aacon_mito/strict_boot
# Majority rule consensus
raxml-ng --consense MR  --tree all_aa.bootstraps --prefix aacon_mito/majority_boot


cat << p

----
MEGA
----
p

# Make MPC supermatrix
cd supermatrix
cp -r nt_coi nt_coi_mpc
cd nt_coi_mpc
rm 1* 2* AK* CA* H* R* W*
cd ../..
github/catfasta2phyml/catfasta2phyml.pl -c -fasta supermatrix/nt_coi_mpc/* > dating/2_nt_mpc_supermatrix.fasta 2> dating/2_nt_mpc_partitions.txt


# Copy best ML tree to dating directory
# Filter supermatrix
cd dating
python3 github/voglerlab/filter_fasta.py -i 2_nt_mpc_supermatrix.fasta -f mpc_supermatrix.fasta -s ml_nt_codon.bestTree -t tree

# Format outgroup and calibration files
sed 's/,/\n/g' ../mptp/outgroup.txt | awk '{print $0 "=outgroup"}' > outgroup.txt
sed -i -E 's/_{2,5}/_/g' outgroup.txt
sed -i -E 's/_{2,5}/_/g' crown.txt

# Run RelTIme
apptainer exec --bind /mnt/shared/apps/ascott/ --no-mount /etc/localtime megacc.sif megacc -a ../../reltime_ml_nucleotide.mao -c crown.txt -d mpc_supermatrix.fasta -f Fasta -o uniform -t *bestTree -g outgroup.txt



cat << p

----
BAMM
----
p

mkdir bamm
cd bamm
Rscript github/voglerlab/bamm_tree_input.R -i ../dating/uniform_exactTimes.nwk -o bamm.tree


