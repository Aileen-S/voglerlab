# Pipeline to add BOLD and GenBank sequences to lab mitogenomes

# -------------------
# Programs to install
# -------------------

# Python, Biopython
# BLAST
# r-base, r-bold, r-getopt, r-dplyr
# MAFFT
# CIAlign
# mPTP
# RAxML-NG

# -------------------
# Github scripts uesd
# -------------------

# https://github.com/Aileen-S/voglerlab
# bold_cli.R
# extract_lab_genes.py
# get_genbanks.py
# fasta_dups_length.py
# filter_fasta.py
# filter_tips.R
# partitions.py
# ptp_filter_output.py
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
  # Profile: nucleotide alignments for aligning RNA genes
  # Profile: protein alignments for all protein coding genes
  # Profiles available at github: ###############################################################################

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
mkdir lab
cd lab

# Add latest MMG database version, eg 'gbmaster_2024-04-20'
python3 get_lab_gbs.py -i lab_ids.txt -o $taxon.gb -v version 
python3 extract_lab_genes.py -g $taxon.gb

# Organise fastas
mkdir frame
mv *fasta frame
mkdir noframe
mv *rf noframe

# Make sure all genes have reading frame. Find reading frame if missing
if [[ -d noframe && ! -z "$(ls -A noframe)" ]]; then
  for file in noframe/*
  do
    base=$(basename ${file%.rf})
    echo $base
    # Add path to profile sequences
    findframe.py -r profiles/aa/$base.fasta -t 5 < $file > noframe/$base.fasta
    cat  noframe/$base.fasta frame/$base.fasta > frame/$base.fasta.a
    mv frame/$base.fasta.a frame/$base.fasta
  done
fi


cat << p

------------
Add Outgroup
------------
p

# Get outgroup, if not already added with lab sequences
mkdir outgroup
cd outgroup

# Add latest MMG database version, eg 'gbmaster_2024-04-20'
python3 get_lab_gbs.py -i outgroup_ids.txt -o $taxon.gb -v version 
python3 extract_lab_genes.py -g $taxon.gb

# Organise fastas
mkdir frame
mv *fasta frame
mkdir noframe
mv *rf noframe

# Find reading frame
if [[ -d noframe && ! -z "$(ls -A noframe)" ]]; then
  for file in noframe/*
  do
    base=$(basename ${file%.rf})
    echo $base
    findframe.py -r profiles/aa/$base.fasta -t 5 < $file > noframe/$base.fasta
    cat  noframe/$base.fasta frame/$base.fasta > frame/$base.fasta.a
    mv frame/$base.fasta.a frame/$base.fasta
  done
fi


cat << p

---------------------
Search BLAST Database
---------------------
p

mkdir genbank
cd genbank

# Direct to BLAST database
# Below is example of location on franklin/HPC: edit as needed
export BLASTDB=/mbl/share/workspaces/groups/database/nt-2024-05-29_blastdb

# Add API key to speed up GenBank search if necessary
# export NCBI_API_KEY=


# Get NCBI taxonomy ID for chosen taxon with script from BLAST package
output=$(get_species_taxids.sh -n $taxon)
txid=$(echo "$output" | grep 'Taxid' | awk '{print $3}')
echo "Taxon ID: $txid"

# Get list of NCBI species taxononmy IDs for chosen taxon
echo "Getting txid list"
get_species_taxids.sh -t $txid > $taxon.txids

# Search BLAST nucleotide database using profile, limit search using taxon ID list, output list of accessions
# Adjust BLAST parameters as desired
blastn -db nt -query profiles/blast/COX1a.fasta -taxidlist $taxon.accs -out COX1a.blast -max_target_seqs 100000 -outfmt '6 sacc'
blastn -db nt -query profiles/blast/COX1b.fasta -taxidlist $taxon.accs -out COX1b.blast -max_target_seqs 100000 -outfmt '6 sacc'


cat << p

---------------------
Get GenBank Sequences
---------------------
p

# Genes available with this script are ATP6, ATP8, COX1, COX2, COX3, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6,
# AK, CAD, EF1A, H3, RNApol, Wg, 12S, 16S, 18S, 28S

# Get COX1a sequences from BLAST search
python3 ~/scratch/github/voglerlab/get_genbanks.py -f COX1a.blast -r gbid -e mixedupvoyage@gmail.com -c
mv COX1.fasta COX1a.fasta
mv metadata.csv meta_COX1a.csv

# Get COX1b sequences from BLAST result
python3 ~/scratch/github/voglerlab/get_genbanks.py -f COX1b.blast -r gbid -e mixedupvoyage@gmail.com -c
mv COX1.fasta COX1b.fasta
mv metadata.csv meta_COX1b.csv

# Get other sequences with keyword search
python3 ~/scratch/github/voglerlab/get_genbanks.py -t $taxon -e mixedupvoyage@gmail.com
rm COX1.fasta
mv metadata.csv meta_full.csv

# Combine metadata
cat meta* > meta_1.csv
sort meta_1.csv | uniq > meta_2.csv
awk '/^ncbi/{print; next} {a[NR]=$0} END{for(i in a) print a[i]}' meta_2.csv > metadata.csv
rm meta_1.csv meta_2.csv

# Get frame tags if necessary
mkdir frame
mv *fasta frame
mkdir noframe
mv *rf noframe

# Find reading frame
if [[ -d noframe && ! -z "$(ls -A noframe)" ]]; then
  for file in noframe/*
  do
    base=$(basename ${file%.rf})
    echo $base
    findframe.py -r profiles/aa/$base.fasta -t 5 < $file > noframe/$base.fasta
    cat  noframe/$base.fasta frame/$base.fasta > frame/$base.fasta.a
    mv frame/$base.fasta.a frame/$base.fasta
  done
fi

cd ..


cat << p

-----------
BOLD search
-----------
p

mkdir bold
cd bold

# Search BOLD database online (or use the same script to search a previously downloaded BOLD file
Rscript bold_cli.R -t $taxon -m ../genbank/metadata.csv -g


cat << p

-----------------------
Get frame tags for BOLD
-----------------------
p

mkdir noframe
mkdir frame

# Find reading frame
for file in profiles/aa/*
do
  base=$(basename $file)
  echo $base
  findframe.py -r $file -t 5 < raw/$base > frame/$base
done

# Copy RNA fastas to 'frame' directory
cp noframe/1* noframe/2* frame
cd ..

cat << p

-------------------
Merge raw sequences
-------------------
p

mkdir alignment
cd alignment

# Define the list of possible file names
file_names=("12S.fasta" "16S.fasta" "18S.fasta" "28S.fasta" "AK.fasta" "ATP6.fasta" "ATP8.fasta" "CAD.fasta" "COX1.fasta" "COX1a.fasta" "COX1b.fasta" "COX2.fasta" "COX3.fasta" "CYTB.fasta" "EF1A.fasta" "H3.fasta" "ND1.fasta" "ND2.fasta" "ND3.fasta" "ND4.fasta" "ND4L.fasta" "ND5.fasta" "ND6.fasta" "Wg.fasta")

# List of directories to check
directories=(
  "../lab/frame"
  "../genbank/frame"
  "../bold/frame"
  "../outgroup/frame"
)

# Join files
for f in "${file_names[@]}"; do
  output="$f"
  : > "$output"
  for dir in "${directories[@]}"; do
    if [ -e "$dir/$f" ]; then
      cat "$dir/$f" >> "$output"
    fi
  done
done

# Delete files with less than 10 sequences

for file in *fasta
do
  sequence_count=$(grep -c "^>" "$file")
  # Check if the sequence count is less than 10
  if (( sequence_count < 10 )); then
    echo "Deleting $file, as it contains only $sequence_count sequences"
    rm "$file"
  fi2533
done

# Split into RNA, mitogenome protein coding and nuclear protein coding
mkdir n1_raw m1_raw r1_raw 
mkdir n2_aa n3_aaal n4_ntal m2_aa m3_aaal m4_ntal r4_ntal
mv 12S.fasta 16S.fasta 18S.fasta 28S.fasta r1_raw
mv COX*fasta CYTB.fasta ATP*fasta ND*fasta m1_raw
mv RNApol.fasta Wg.fasta EF1A.fasta H3.fasta CAD.fasta AK.fasta n1_raw


cat << p

--------------
COX1 sequences
--------------
p
# Merge COX1a and COX1b, removing duplicates. These will be separated again later.
mkdir COI
mv m1_raw/COX1* COI
cat COI/COX1* > COI/COX1all.fasta
python3 fasta_dups_length.py -i COI/COX1all.fasta -o m1_raw/COX1.fasta -d

cat << p

--------------------
Translate to Protein
--------------------
p

# Nuclear genes
for file in n1_raw/*
do
  echo 'translate' $file
  translate.py 1 < $file > n2_aa/${file#*/}
done

# Mitochondrial genes
for file in m1_raw/*
do
  echo 'translate' $file
  translate.py 5 < $file > m2_aa/${file#*/}
done

cat << p

-----------------
Align to profiles
-----------------
p
# Adjust mafft parameters and thread count as necessary

# Mitochondrial
for file in m2_aa/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 12 profiles/aa/${file#*/} >  m3_aaal/${file#*/}
done

# Nuclear
for file in n2_aa/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 12 profiles/aa/${file#*/} >  n3_aaal/${file#*/}
done

# RNA
for file in r1_raw/*
do
  echo $file
  mafft --add $file --maxiterate 2 --adjustdirection --thread 12 profiles/nt/${file#*/} >  r4_ntal/${file#*/}
done


cat << p

---------------
Remove Profiles
---------------
p

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
  echo $file
  backtranslate.py -i $file m1_raw/${file#*/} 5 > m4_ntal/${file#*/}
done

# Nuclear
for file in n3_aaal/*
do
  echo $file
  backtranslate.py -i $file n1_raw/${file#*/} 5 > n4_ntal/${file#*/}
done

# Tidy up sequence ends and remove duplicates
for file in m4_ntal/*
do
  base=$(basename $file)
  fasta_dups_length.py -i $file -o m4_ntal/$base.f -l -d
  mv m4_ntal/$base.f m4_ntal/$base
done

for file in n4_ntal/*
do
  base=$(basename $file)
  fasta_dups_length.py -i $file -o n4_ntal/$base.f -l -d
  mv n4_ntal/$base.f n4_ntal/$base
done

for file in r4_ntal/*
do
  base=$(basename $file)
  fasta_dups_length.py -i $file -o r4_ntal/$base.f -l -d
  mv r4_ntal/$base.f r4_ntal/$base
done


cat << p

----------------
Clean Alignments
----------------
p

# Adjust CIAlign parameters as necessary

mkdir m5_aa_cia m6_nt_cia n5_aa_cia n6_nt_cia r6_nt_cia

for file in m3_aaal/*; do
  echo $file
  base=$(basename $file)
  if [[ $base == ATP8.fasta ]]; then
    CIAlign --infile $file --outfile_stem m5_aa_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.5
  else
    CIAlign --infile $file --outfile_stem m5_aa_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.5 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
  fi
done

for file in m4_ntal/*; do
  echo $file
  base=$(basename $file)
  if [[ $base == ATP8.fasta ]]; then
    CIAlign --infile $file --outfile_stem m6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.6 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
  else
    CIAlign --infile $file --outfile_stem m6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.65 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
  fi
done
 
for file in n3_aaal/*
do
  echo $file
  base=$(basename $file)
  CIAlign --infile $file --outfile_stem n5_aa_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.5 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
done

for file in n4_ntal/*
do
  echo $file
  base=$(basename $file)
  CIAlign --infile $file --outfile_stem n6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.65 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
done

for file in r4_ntal/*
do
  echo $file
  base=$(basename $file)
  CIAlign --infile $file --outfile_stem r6_nt_cia/${base%.fasta} --remove_divergent --remove_divergent_minperc 0.65 --crop_divergent --crop_divergent_min_prop_nongap 0.1 --remove_insertions
done


# Split COX1 into COX1a and COX1b
# Do this manually instead if you prefer

cd m6_nt_cia
awk '
/^>/ {  # When encountering a header line
    if (NR > 1) {  # For all except the first header line
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only process non-empty sequences that are not all gaps
            print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
        }
    }
    header = $0  # Save the header
    seq = ""  # Reset the sequence
    next
}
{
    seq = seq $0  # Append sequence lines
}
END {
    if (length(seq) > 0 && !is_all_gaps(seq)) {  # Process the last sequence if it is not all gaps
        print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
    }
}

function print_seq(file1, file2, header, seq,    len1, len2) {
    len1 = 720
    len2 = length(seq)
    if (len2 > len1) {
        if (length(substr(seq, 1, len1)) > 0 && !is_all_gaps(substr(seq, 1, len1))) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print substr(seq, 1, len1) > file1
        }
        if (length(substr(seq, len1 + 1)) > 0 && !is_all_gaps(substr(seq, len1 + 1))) {  # Only write non-empty sequences that are not all gaps
            print header > file2
            print substr(seq, len1 + 1) > file2
        }
    } else {
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print seq > file1
        }
    }
}

function is_all_gaps(seq,    i) {
    for (i = 1; i <= length(seq); i++) {
        if (substr(seq, i, 1) != "-") {
            return 0  # Not all gaps
        }
    }
    return 1  # All gaps
}
' COX1_cleaned.fasta

cd ../m5_aa_cia
awk '
/^>/ {  # When encountering a header line
    if (NR > 1) {  # For all except the first header line
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only process non-empty sequences that are not all gaps
            print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
        }
    }
    header = $0  # Save the header
    seq = ""  # Reset the sequence
    next
}
{
    seq = seq $0  # Append sequence lines
}
END {
    if (length(seq) > 0 && !is_all_gaps(seq)) {  # Process the last sequence if it is not all gaps
        print_seq("COX1a.fasta", "COX1b.fasta", header, seq)
    }
}

function print_seq(file1, file2, header, seq,    len1, len2) {
    len1 = 250
    len2 = length(seq)
    if (len2 > len1) {
        if (length(substr(seq, 1, len1)) > 0 && !is_all_gaps(substr(seq, 1, len1))) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print substr(seq, 1, len1) > file1
        }
        if (length(substr(seq, len1 + 1)) > 0 && !is_all_gaps(substr(seq, len1 + 1))) {  # Only write non-empty sequences that are not all gaps
            print header > file2
            print substr(seq, len1 + 1) > file2
        }
    } else {
        if (length(seq) > 0 && !is_all_gaps(seq)) {  # Only write non-empty sequences that are not all gaps
            print header > file1
            print seq > file1
        }
    }
}

function is_all_gaps(seq,    i) {
    for (i = 1; i <= length(seq); i++) {
        if (substr(seq, i, 1) != "-") {
            return 0  # Not all gaps
        }
    }
    return 1  # All gaps
}
' COX1_cleaned.fasta
cd ..


# Check alignment files before proceeding


cat << p

-------------------------
Get ready for supermatrix
-------------------------
p


mkdir ../supermatrix
mkdir ../supermatrix/aa
mkdir ../supermatrix/nt

# Copy aligned files to aa and nt directories for futher processing
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


# Remove _R_ prefix from reversed sequences
for file in aa/*
do
   sed -i -E "s/^>_R_/>/" $file
done

for file in nt/*
do
   sed -i -E "s/^>_R_/>/" $file
done


# Combine metadata

cd ..

Rscript combine_metadata.R -b bold/metadata.csv -g genbank/metadata.csv -l lab/lab.ids -o outgroup/outgroup.ids -m ../SITE-100_Database_mitogenomes_240731.csv -c supermatrix/rename.csv

cd supermatrix


# Rename sequences with NCBI taxon IDs and taxonomy

mkdir ids_aa ids_nt
for file in nt/*
do
  base=$(basename $file)
  python3 rename_fasta.py -i $file -c rename.csv -r ids_nt/${base%fasta}csv -o nt/$base.r
  mv nt/${base#*/}.r nt/${base#*/}
done

for file in aa/*
do
  base=$(basename $file)
  python3 rename_fasta.py -i $file -c rename.csv -r ids_aa/${base%fasta}csv -o aa/$base.r
  mv aa/${base#*/}.r aa/${base#*/}
done


# Remove sequences without COI, if desired
# Remove sequences without barcodes by filter using COX1a.fasta instead of COX1_cleaned.fasta

mkdir aa_coi nt_coi

for file in aa/*
do
  echo $file
  python3 filter_fasta.py -i $file -s aa/COX1_cleaned.fasta -t fasta -f aa_coi/${file#*/}
  # Remove gap only columns
  degap_fasta_alignment.pl aa_coi/${file#*/} > aa_coi/${file#*/}.g
  mv aa_coi/${file#*/}.g aa_coi/${file#*/}
done

for file in nt/*
do
  echo $file
  python3 filter_fasta.py -i $file -s aa/COX1_cleaned.fasta -t fasta -f nt_coi/${file#*/}
  # Remove gap only columns
  degap_fasta_alignment.pl nt_coi/${file#*/} > nt_coi/${file#*/}.g
  mv nt_coi/${file#*/}.g nt_coi/${file#*/}
done

rm */COX1_cleaned.fasta


# Delete files with less than 10 sequences

for file in aa_coi/*
do
  sequence_count=$(grep -c "^>" "$file")
  # Check if the sequence count is less than 10
  if (( sequence_count < 10 )); then
    echo "Deleting $file, as it contains only $sequence_count sequences"
    rm "$file"
  fi2533
done

for file in nt_coi/*
do
  sequence_count=$(grep -c "^>" "$file")
  # Check if the sequence count is less than 10
  if (( sequence_count < 10 )); then
    echo "Deleting $file, as it contains only $sequence_count sequences"
    rm "$file"
  fi
done


cat << p

--------------------------------------
Create supermatrix and partition files
--------------------------------------
p

# Supermatrix
catfasta2phyml.pl -c -fasta aa_coi/* > 1_aa_coi_supermatrix.fasta 2> 1_aa_coi_partitions.txt
catfasta2phyml.pl -c -fasta nt_coi/* > 2_nt_coi_supermatrix.fasta 2> 2_nt_coi_partitions.txt


# Change partition files to RAXML compatible format
# Change models as desired
python3 ~/scratch/github/genbank/partitions.py -i 1_aa_coi_partitions.txt -t aa -o raxml
python3 ~/scratch/github/genbank/partitions.py -i 2_nt_coi_partitions.txt -t nt -o raxml

# Check supermatrices and partition files before proceeding to phylogeny

cd ..


cat << p

-----------------
Initial Phylogeny
-----------------
p

mkdir mptp

# Adjust threads as needed
# To check how many threads you need, run:
#    raxml-ng --parse --msa supermatrix/2_nt_coi_supermatrix.fasta --prefix mptp/parse --model supermatrix/partitions_gene.txt

raxml-ng --search1 --msa supermatrix/2_nt_coi_supermatrix.fasta --prefix mptp/search1 --model supermatrix/partitions_gene.txt --threads 8


cat << p

-------------------------
mPTP Species Delimitation
-------------------------
p

# Requires file with comma separated list of outgroup taxa

cd mptp
# Get outgroup
outgroup=$(awk 'NR==1' outgroup.txt)

# Calculate minbr value
echo "Getting minbr value"
mptp --minbr_auto ../supermatrix/2_nt_coi_supermatrix.fasta --tree_file search1.raxml.bestTree --output_file minbr  --outgroup $outgroup

# Copy minbr value from slurm output
minbr=$(awk 'END {print $NF}' minbr.txt)
echo "minbar value = " $minbr

# Run single rate mPTP
echo "Running mPTP"
mptp --ml --single --minbr $minbr --tree_file search1.raxml.bestTree --output_file mptp  --outgroup $outgroup

# Add multi-rate mPTP and MCMC sampling if needed

# Filter taxa

# Choose taxon with most nucleotides for each PTP species
python3 ptp_filter_output.py -i mptp.txt -s ../supermatrix/2_nt_coi_supermatrix.fasta -o ptp_chosen.txt
# Filter supermatrix
python3 filter_fasta.py -i ../supermatrix/2_nt_coi_supermatrix.fasta -f ../supermatrix/2_nt_coi_ptp_supermatrix.fasta -s ptp_chosen.txt -t list

cd ..


cat << p

-----
RAxML
-----
p

mkdir raxml
cd raxml

# Run RaxML
raxml-ng --search1 --msa ../supermatrix/2_nt_coi_ptp_supermatrix.fasta --model ../supermatrix/partitions_gene.txt --prefix ptp_search1 --threads 8

# Remove long branches
run_treeshrink.py -t ptp_search1.raxml.bestTree -o treeshrink

# Check trees
# Remove more long branches if necessary using filter_tips.R, then filter supermatrix

# Filter supermatrix
python3 filter_fasta.py -i ../supermatrix/2_nt_coi_supermatrix.fasta -f ../supermatrix/2_nt_coi_ptp_ts_supermatrix.fasta -s treeshrink/output.bestTree -t tree
python3 filter_fasta.py -i ../supermatrix/1_aa_coi_supermatrix.fasta -f ../supermatrix/1_aa_coi_ptp_ts_supermatrix.fasta -s treeshrink/output.bestTree -t tree
python3 ry_code.py -i ../supermatrix/2_nt_coi_ptp_ts_supermatrix.fasta -o ../supermatrix/3_ry_coi_ptp_ts_supermatrix.fasta

# Proceed with tree search with chosen supermatrices and partition schemes

# Eg:
# Nucleotide, gene partitions
raxml-ng --msa ../supermatrix/2_nt_coi_ptp_ts_supermatrix.fasta --model ../supermatrix/partitions_gene.txt --prefix ml_nt_gene --threads 8
# Nucleotide, codon partitions
raxml-ng --msa ../supermatrix/2_nt_coi_ptp_ts_supermatrix.fasta --model ../supermatrix/partitions_codon123.txt --prefix ml_nt_codon --threads 8
# Protein, gene partitions
raxml-ng --msa ../supermatrix/1_aa_coi_ptp_ts_supermatrix.fasta --model ../supermatrix/partitions_aa.txt --prefix ml_aa --threads 8
# RY code, gene partitions
raxml-ng --msa ../supermatrix/3_ry_coi_ptp_ts_supermatrix.fasta --model ../supermatrix/partitions_gene.txt --prefix ml_ry_gene --threads 8
# RY code, codon partitions
raxml-ng --msa ../supermatrix/3_ry_coi_ptp_ts_supermatrix.fasta --model ../supermatrix/partitions_codon123.txt --prefix ml_ry_codon --threads 8