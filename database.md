# HiMAP database

These notes describe how to construct a HiMAP database from scratch. First, clone the repository:
```sh
git clone https://github.com/taolonglab/himapdb.git
cd himapdb
```

To run these scripts you need:
* Python 3, with modules Biopython and Pandas installed
* blastn: command line NCBI nucleotide blastn installed


## Metadata

Download the NCBI Genome database metadata for Bacterial and Archaeal genome assemblies and select only the best assembly for each unique strain:
```sh
python3 py/dl_genome_meta.py data/

```

This will first download two files:
* `data/archaea_assembly_summary_2018-05-20.txt`
* `data/bacteria_assembly_summary_2018-05-20.txt`

then generate a filtered assembly summary with unique strain names:
* `data/archaea_assembly_summary_filter_2018-05-20.txt`
* `data/bacteria_assembly_summary_filter_2018-05-20.txt`

and files with direct links to download FASTA sequences and GFF annotations:
* `data/archaea_assembly_features_urls.txt`
* `data/archaea_assembly_features_urls.txt`
* `data/bacteria_assembly_features_urls.txt`
* `data/bacteria_assembly_features_urls.txt`

Note: if you want to download only specific kingdom(s) these can be given as a second argument, i.e.
```sh
python3 py/dl_genome_meta.py data/ archaea
```

or both (same as omitting it):

```sh
python3 py/dl_genome_meta.py data/ archaea,bacteria
```

## Download data

Now, download the sequences. This will be a large download and make take a day or so to complete.
```sh
mkdir data/sequences
mkdir data/features

cd data/features
python3 ../../py/download.py ../archaea_assembly_features_urls.txt 10
python3 ../../py/download.py ../bacteria_assembly_features_urls.txt 10

cd ../sequences
python3 ../../py/download.py ../archaea_assembly_sequences_urls.txt 10
python3 ../../py/download.py ../bacteria_assembly_sequences_urls.txt 10
```

where the second argument specifies the number of concurrent download connections to speed up the download. Increasing number of connections can speed up the download significantly but also lead to some failed downloads (files with size 0). These can be re-downloaded, but deleting them and re-running the script with 1 concurrent download option. 

For example to re-download missing sequences, run this from the `sequences` folder:

```sh
for f in $(find . -size 0); do rm $f; done
python3 ../../py/download.py ../archaea_assembly_sequences_urls.txt 1
```

while some downloads will always result in 0 size files (next script will just skip those).


## Extract 16S sequences from full genomes

Extract 16S ribosomal RNA gene sequences from FASTA sequences, based on their GFF annotations. No need to uncompress anything. Switch back to the `himapdb` folder and run `extract.py`. This script takes 3 inputs: (i) filtered assembly summary file, (ii) folder with downloaded GFF features and (iii) folder with downloaded sequences. The last argument is the output FASTA file:

```sh
cd ../../

python3 py/extract.py \
    data/archaea_assembly_summary_filter_2018-05-20.txt \
    data/features/ data/sequences/ \
    data/16s_from_genomes_2018-05-20.fasta

python3 py/extract.py \
    data/bacteria_assembly_summary_filter_2018-05-20.txt \
    data/features/ data/sequences/ \
    data/16s_from_genomes_2018-05-20.fasta

```

The last line will add Bacteria sequences to the existing output FASTA files.

## Generate files for each primer set

```sh

python3 py/primers_fasta.py \
    primers/nucleotide_codes.txt \
    primers/pcr_primers_table.txt \
    data/pcr_primers    
```

The `data/pcr_primers` folder now contains FASTA files for each primer set:
* V1-V8_8F-1392R.fasta
* V3-V4-2_341F-805R.fasta
* V3-V4-V5-V8_337F-1392R.fasta
* V3-V4-V5_337F-926R.fasta
* V3-V4_337F-805R.fasta
* V4_515F-805R.fasta


## Count hypervariable variants

```sh
python3 py/count.py \
    
```

## Add RefSeq hypervariable variants

We go through each strain name from the NCBI RefSeq search results and add sequences if the strain name is not an exact match to something already in the database.