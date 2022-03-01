# 16S Amplicon sequencing of grapevine samples from Central Valley California vineyards 
## Author : Joel F. Swift
## Date : 08/2021

The purpose of this README is to document the steps used in QIIME2 to process the data prior to import into R.
The sequencing was conducted across three runs (i.e. PLATE1|2|3). For plates 1 and 2, the set of 672 samples were
broken into two sets. For plate 3, 345 samples were resequenced to improve the number sequencing depth (for the 
most part these were berry/leaf samples).


#### 1st PLATE ####
Undetermined_S0_L001_I1_001.fastq.gz
MD5: 8eaa6c4614cfb810741e16793f0f4830

Undetermined_S0_L001_R1_001.fastq.gz
MD5: 7e29c56aa000d6424c1b369fe91a412c

Undetermined_S0_L001_R2_001.fastq.gz
MD5: a671e0aa16a41eb3e48e6171b6175a43


##### Importing data into a Qiime2 artifact (.qza file)

Files were renamed to fit the Qiime2 EMPPairedEndSequences, for example forward.fastq, reverse.fastq, and barcodes.fastq

```bash
cp Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
cp Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
cp Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz


qiime tools import --type EMPPairedEndSequences --input-path ~/projects/CHAP2/RAWDATA/ --output-path 16S_CA_PE_PLATE1_seq.qza
```

##### Processing of data

```bash
# Demultiplexing
qiime demux emp-paired \
  --m-barcodes-file metadata_plate1.tsv \
  --m-barcodes-column barcode-sequence \
  --i-seqs 16S_CA_PE_PLATE1_seq.qza \
  --o-per-sample-sequences demux_plate1_barcode_mapping_revcomp.qza \
  --o-error-correction-details error_correction_details_plate1.qza \
  --p-rev-comp-mapping-barcodes \
  --p-rev-comp-barcodes


qiime demux summarize --i-data demux_plate1_barcode_mapping_revcomp.qza --o-visualization demux_plate1_barcode_mapping_revcomp.qzv

# Denoising and merging pair-end reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_plate1_barcode_mapping_revcomp.qza \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --o-table table_plate1.qza \
  --o-representative-sequences rep-seqs_plate1.qza \
  --o-denoising-stats denoising-stats_plate1.qza \
  --verbose \
  --p-n-threads 2

#Visualize dada2 results
qiime metadata tabulate \
  --m-input-file denoising-stats_plate1.qza \
  --o-visualization denoising-stats-trimmed_plate1.qzv
```

#### 2nd PLATE ####
Undetermined_S0_L001_I1_001_plate2.fastq.gz
MD5: 6d691099050044ccaf201122d0b405c6

Undetermined_S0_L001_R1_001_Plate2.fastq.gz
MD5: 15852e370342c53e95e4c496c8c2458f

Undetermined_S0_L001_R2_001_Plate2.fastq.gz
MD5: daf5bf03d4e2b7d1f7675e4cf9a871b5

##### Importing data into a Qiime2 artifact (.qza file)

Files were renamed to fit the Qiime2 EMPPairedEndSequences, for example forward.fastq, reverse.fastq, and barcodes.fastq

```bash
cp Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
cp Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
cp Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz


qiime tools import --type EMPPairedEndSequences --input-path ./ --output-path 16S_CA_PE_PLATE2_seq.qza
```

##### Processing of data

```bash
# Demultiplexing
qiime demux emp-paired \
  --m-barcodes-file metadata_plate2.tsv \
  --m-barcodes-column barcode-sequence \
  --i-seqs 16S_CA_PE_PLATE2_seq.qza \
  --o-per-sample-sequences demux_plate2_barcode_mapping_revcomp.qza \
  --o-error-correction-details error_correction_details_plate2.qza \
  --p-rev-comp-mapping-barcodes \
  --p-rev-comp-barcodes


qiime demux summarize --i-data demux_plate2_barcode_mapping_revcomp.qza --o-visualization demux_plate2_barcode_mapping_revcomp.qzv

# Denoising and merging pair-end reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_plate2_barcode_mapping_revcomp.qza \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --o-table table_plate2.qza \
  --o-representative-sequences rep-seqs_plate2.qza \
  --o-denoising-stats denoising-stats_plate2.qza \
  --verbose \
  --p-n-threads 2

#Visualize dada2 results
qiime metadata tabulate \
  --m-input-file denoising-stats_plate2.qza \
  --o-visualization denoising-stats-trimmed_plate2.qzv
```

#### 3rd PLATE ####
Undetermined_S0_L001_I1_001.fastq.gz
MD5: 008469604166accae1f77db4257d101c

Undetermined_S0_L001_R1_001.fastq.gz
MD5: 481b4f4e5a7954bb6a114297f21a16f8

Undetermined_S0_L001_R2_001.fastq.gz
MD5: 157f03eb3a698703bb429dbe6a5b5c00

```bash
cp Undetermined_S0_L001_I1_001.fastq.gz barcodes.fastq.gz
cp Undetermined_S0_L001_R1_001.fastq.gz forward.fastq.gz
cp Undetermined_S0_L001_R2_001.fastq.gz reverse.fastq.gz


qiime tools import --type EMPPairedEndSequences --input-path /data/projects/joel.swift/projects/CHAP2/RAWDATA/PLATE3/ --output-path 16S_CA_PE_PLATE3_seq.qza
```

##### Processing of data

```bash
# Demultiplexing
qiime demux emp-paired \
  --m-barcodes-file metadata_plate3.csv \
  --m-barcodes-column barcode-sequence \
  --i-seqs 16S_CA_PE_PLATE3_seq.qza \
  --o-per-sample-sequences demux_plate3_barcode_mapping_revcomp.qza \
  --o-error-correction-details error_correction_details_plate3.qza \
  --p-rev-comp-mapping-barcodes \
  --p-rev-comp-barcodes


qiime demux summarize --i-data demux_plate3_barcode_mapping_revcomp.qza --o-visualization demux_plate3_barcode_mapping_revcomp.qzv

# Denoising and merging pair-end reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_plate3_barcode_mapping_revcomp.qza \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --o-table table_plate3.qza \
  --o-representative-sequences rep-seqs_plate3.qza \
  --o-denoising-stats denoising-stats_plate3.qza \
  --verbose \
  --p-n-threads 0

#Visualize dada2 results
qiime metadata tabulate \
  --m-input-file denoising-stats_plate3.qza \
  --o-visualization denoising-stats-trimmed_plate3.qzv
```

#### Merging Plates together ####

```bash
# merge tables 1,2,3 together and add the reads together.
qiime feature-table merge \
  --i-tables table_plate1.qza \
  --i-tables table_plate2.qza \
  --i-tables table_plate3.qza \
  --p-overlap-method 'sum' \
  --o-merged-table table_merged.qza

# Sanity check to see if merge table looks good and also check sample read depth.
# Also added a column to metadata that denotes sequencing run (#seq_run 1|2|3)
qiime feature-table summarize \
  --i-table table_merged.qza \
  --o-visualization table_merged.qzv \
  --m-sample-metadata-file Metadata_all_plates.tsv

#Also merge rep sequences
qiime feature-table merge-seqs \
  --i-data rep-seqs_plate1.qza \
  --i-data rep-seqs_plate2.qza \
  --i-data rep-seqs_plate3.qza \
  --o-merged-data rep-seqs_merged.qza
```

##### Taxonomy Assignment

```bash
# Getting the pretrained classifier
wget https://data.qiime2.org/2021.4/common/silva-138-99-515-806-nb-classifier.qza

# verifing md5sum e05afad0fe87542704be96ff483824d4
md5sum silva-138-99-515-806-nb-classifier.qza

# Obtaining taxonomy for the 16s sequences
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs_merged.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv


# Determining the number of mitochondrial sequences removed
qiime taxa filter-table \
  --i-table table_merged.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria \
  --o-filtered-table table-no-mitochondria.qza

qiime feature-table summarize \
  --i-table table-no-mitochondria.qza \
  --o-visualization table-no-mitochondria.qzv

#Total sequences: 32,208,314
#After filtering: 31,485,950
#Difference: 722,364
#% Diff: 2.24%

# Determining the number of Chloroplast sequences removed
qiime taxa filter-table \
  --i-table table_merged.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude chloroplast \
  --o-filtered-table table-no-chloroplast.qza

qiime feature-table summarize \
  --i-table table-no-chloroplast.qza \
  --o-visualization table-no-chloroplast.qzv

#Total Sequences: 32,208,314
#After filtering: 31,564,440
#Difference: 643,874
#% Diff: 2.00%


# Filtering out unassigned sequences without a phylum 
qiime taxa filter-table \
  --i-table table_merged.qza \
  --i-taxonomy taxonomy.qza \
  --p-include p__ \
  --o-filtered-table filtered_table-with-phyla.qza

# Filtering out mitochondria and chloroplast sequences
qiime taxa filter-table \
  --i-table filtered_table-with-phyla.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table table-no-mitochondria-no-chloroplast.qza


# Sanity check, make sure with mito and chloro gone the total is equal to
# agrees numbers
qiime feature-table summarize \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --o-visualization table-no-mitochondria-no-chloroplast.qzv

# Creating bar graph with low abudance ASVS
qiime taxa barplot \
  --i-table table-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Metadata_all_plates.tsv \
  --o-visualization taxa-bar-plots-no-mito-no-chloro.qzv
```