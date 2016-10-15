# 基因组组装 -- pacbio和nanopore序列组装  
## 一、实验目的  
## 二、知识回顾  
## 三、上机操作  

```
# 登录服务器
# 设置环境变量
module bioinfo
# 准备工作目录和数据
mkdir lab1_2
cd lab1_2
ln -s /bs1/data/genomeLab/lab1.2/data/selfSampleData/pacbio_filtered.fastq ./

```

### 一步法  
By default, canu will correct the reads, then trim the reads, then assemble the reads to unitigs.

```
canu \
 -p ecoli -d ecoli-auto \
 genomeSize=4.8m \
 -pacbio-raw pacbio_filtered.fastq

# 组装结果
ecoli-auto/ecoli.contigs.fasta

```

### 分步法  
```
第一步：纠错
canu -correct \
  -p ecoli -d ecoli \
  genomeSize=4.8m \
  -pacbio-raw pacbio_filtered.fastq
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
```

While Canu corrects sequences and has 99% identity or greater with PacBio or Nanopore sequences, for the best accuracy we recommend polishing with a sequence-specific tool. We recommend [Quiver](https://github.com/PacificBiosciences/GenomicConsensus) for PacBio and [Nanopolish](http://github.com/jts/nanopolish) for Oxford Nanpore data.

If you have Illumina sequences available, [Pilon](http://www.broadinstitute.org/software/pilon/) can also be used to polish either PacBio or Oxford Nanopore assemblies.


## 四、作业与思考  
1. 用miniasm组装同样的数据  
2. 比较miniasm和canu组装的结果  

### miniasm组装  
```
# Overlap
minimap -Sw5 -L100 -m0 -t8 pacbio_filtered.fastq pacbio_filtered.fastq | gzip -1 > reads.paf.gz
# Layout
miniasm -f pacbio_filtered.fastq reads.paf.gz > reads.gfa
# 将gfa转换为fasta
awk '/^S/{print ">"$2"\n"$3}' reads.gfa | fold > miniasm.fa

```

## 五、参考文献  
1. http://canu.readthedocs.io/en/latest/index.html  
2. https://github.com/PacificBiosciences/GenomicConsensus  
3. Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM. [Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation.](http://biorxiv.org/content/early/2016/08/24/071282) bioRxiv. (2016).
