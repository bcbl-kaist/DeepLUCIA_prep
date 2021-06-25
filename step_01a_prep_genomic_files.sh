
###############################################

 mkdir -p input/genome_fasta/hg19/

 wget http://mitra.stanford.edu/kundaje/projects/dragonn/hg19.genome.fa.gz     -O input/genome_fasta/hg19/ref_genome.fa.gz    
 wget http://mitra.stanford.edu/kundaje/projects/dragonn/hg19.genome.fa.gz.fai -O input/genome_fasta/hg19/ref_genome.fa.gz.fai
 wget http://mitra.stanford.edu/kundaje/projects/dragonn/hg19.genome.fa.gz.gzi -O input/genome_fasta/hg19/ref_genome.fa.gz.gzi
 wget -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz | gunzip | grep -v "_\|chrY\|chrM" | sort -k1,1V > input/genome_fasta/hg19/chromInfo.txt

###############################################

 mkdir -p input/black_list/hg19/

 wget https://www.encodeproject.org/files/ENCFF200UUD/@@download/ENCFF200UUD.bed.gz -O input/black_list/hg19/blacklist.bed.gz
 wget -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz | gunzip | cut -f 2,3,4 | gzip > input/black_list/hg19/gap.bed.gz

###############################################
 
 mkdir -p input/hiccups_loop/hg19/E017/
 mkdir -p input/hiccups_loop/hg19/E116/
 mkdir -p input/hiccups_loop/hg19/E117/
 mkdir -p input/hiccups_loop/hg19/E119/
 mkdir -p input/hiccups_loop/hg19/E122/
 mkdir -p input/hiccups_loop/hg19/E123/
 mkdir -p input/hiccups_loop/hg19/E127/

 wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_IMR90_HiCCUPS_looplist.txt.gz                     -O input/hiccups_loop/hg19/E017/E017_loop.bedpe.gz
 wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt.gz -O input/hiccups_loop/hg19/E116/E116_loop.bedpe.gz
 wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_HeLa_HiCCUPS_looplist.txt.gz                      -O input/hiccups_loop/hg19/E117/E117_loop.bedpe.gz
 wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_HMEC_HiCCUPS_looplist.txt.gz                      -O input/hiccups_loop/hg19/E119/E119_loop.bedpe.gz
 wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_HUVEC_HiCCUPS_looplist.txt.gz                     -O input/hiccups_loop/hg19/E122/E122_loop.bedpe.gz
 wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_K562_HiCCUPS_looplist.txt.gz                      -O input/hiccups_loop/hg19/E123/E123_loop.bedpe.gz
 wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_NHEK_HiCCUPS_looplist.txt.gz                      -O input/hiccups_loop/hg19/E127/E127_loop.bedpe.gz

###############################################

