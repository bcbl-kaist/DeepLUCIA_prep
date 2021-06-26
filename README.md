# DeepLUCIA_prep
creating features for DeepLUCIA ( https://github.com/bcbl-kaist/DeepLUCIA )

# Input data
- ChromInfo file for chromosome sizes. ( retrieved from UCSC goldenpath )
- BED files for gaps. ( retrieved from UCSC goldenpath )
- Bed file for blacklist regions. ( retrieved from ENCODE Project )
- The genome sequence fasta file is used to generate genomic feature file
- histone ChIP-seq -log(p) signal track/ in bigWig format is used to generate epigenomic feature file.
  - The Roadmap Epigenomics provides epigenome of more than 100 human tissues
  - If you have your own data, you can follow the guide from Roadmap Epigenomics ( https://egg2.wustl.edu/roadmap/web_portal/processed_data.html )
- chromatin loops file in bedpe format.

# Parameters
- The length of genomic region to test ( DeepLUCIA used 5kb for default )
- The length of overlap between regions ( DeepLUCIA used 0bp for default )
- The resolution of epigenomic feature ( DeepLUCIA used 25bp for default )

# Output files
- putativa anchor annotations. 
- genomic loop annotations.
- One-hot encoded genomic feature npy file
- 25bp-resolution epigenomic feature npy file.

# Dependencies
- pybedtools>=0.8.0 ( https://github.com/daler/pybedtools )
- pysam>=0.15.2 ( https://github.com/pysam-developers/pysam )
- pyBigWig>=0.3.13 ( https://github.com/deeptools/pyBigWig )
