# 基因组组装 -- pacbio和nanopore序列组装  
## 一、实验目的  
## 二、知识回顾  
## 三、上机操作  
```
软件安装
git clone https://github.com/marbl/canu.git
canu软件已经下载到服务器，存放路径：/bs1/data/genomeLab/lab1.2/soft/
cd canu/src
make -j8

# By default, canu will correct the reads, then trim the reads, then assemble the reads to unitigs.
一步法
canu \
 -p ecoli -d ecoli-auto \
 genomeSize=4.8m \
 -pacbio-raw p6.25x.fastq

分步法
第一步：纠错
canu -correct \
  -p ecoli -d ecoli \
  genomeSize=4.8m \
  -pacbio-raw  p6.25x.fastq
第二步：trim reads
canu -trim \
  -p ecoli -d ecoli \
  genomeSize=4.8m \
  -pacbio-corrected ecoli/correction/ecoli.correctedReads.fasta.gz

第三步：组装
canu -assemble \
  -p ecoli -d ecoli-erate-0.013 \
  genomeSize=4.8m \
  errorRate=0.013 \
  -pacbio-corrected ecoli/trimming/ecoli.trimmedReads.fasta.gz

canu -assemble \
  -p ecoli -d ecoli-erate-0.025 \
  genomeSize=4.8m \
  errorRate=0.025 \
  -pacbio-corrected ecoli/trimming/ecoli.trimmedReads.fasta.gz
  

```
## 四、作业与思考  
## 五、参考文献  
http://canu.readthedocs.io/en/latest/index.html  
https://github.com/PacificBiosciences/GenomicConsensus  

