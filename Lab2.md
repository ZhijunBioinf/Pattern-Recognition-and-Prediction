# 实验二 基因组比对  
## 一、实验目的  
1. 理解比对（mapping, alignment）的含义  
2. 理解全局比对和局部比对的区别和应用  
3. 掌握应用bwa, minimap2, samtools的使用  
4. 理解SAM, BAM文件格式  

## 二、知识回顾  

### 比对的两种策略  
1. Global alignment
2. Local alignment

我们熟悉的blast和blat均属于第二类。   
另外，不同长度的reads比对所用的策略也不一样，对于短reads，基于local alignment的软件如blast, blat不适合。  
将短的reads回帖到长的参考基因组上，这一过程称之为mapping。一般reads数目很大，读长短，参考基因组较长，对于mapping软件有两个要求：

> 1. 速度
> 2. 准确性
 
Mapping软件众多，比较有名的包括bwa, soap, bowtie, novoalign  
另外，由于真核生物mRNA不含有内含子，与一般的DNA mapping软件要求不一样，故转录组mapping使用的软件也不一样，转录组mapping软件比较有名的包括：STAR, hisat  
本实验主要介绍一般意义上的DNA mapping软件的使用。  

### What makes mapping challenging?（挑战）
1. Volume of data
2. Garbage reads
3. Errors in reads, and quality scores
4. Repeat elements and multicopy sequence
5. SNPs/SNVs
6. Indels
7. Splicing (transcriptome)

### 几个影响mapping速度的参数  
1. How many mismatches to allow?
2. Report how many matches?
3. Require best match, or first/any that fit criteria?

## 三、上机操作  
### 设置环境变量和准备工作目录  
```
mkdir lab02
cd lab02
mkdir data
mkdir results
```

### 实验数据  
```
/data/lab/genomic/lab02/data/REL606.fa (参考序列)
/data/lab/genomic/lab02/data/SRR098038.fastq.gz （illumina reads）
/data/lab/genomic/lab02/data/pb_ecoli_0001.fastq （pacbio reads）
```
### (一) Mapping using bwa   

#### 1. 准备数据和index参考基因组  
```
$ cd data
$ ln -s /data/lab/genomic/lab02/data/REL606.fa /data/lab/genomic/lab02/data/reads_* ../data/pb_ecoli_0001.fastq ./
$ samtools faidx REL606.fa
$ mkdir index
$ cd index
$ ln -s ../REL606.fa ./
```
#### 2. 建索引文件  
work_bwaIndex.sh  
```
#PBS -N bwaIdx_W
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
bwa index REL606.fa
```

#### 3. Mapping the reads to the reference genome using bwa  
```
cd ../../result
```
work_bwa.sh
```
#PBS -N bwa
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
bwa mem ../data/index/REL606.fa ../data/reads_1.fq.gz ../data/reads_2.fq.gz > mapping.sam
samtools view -b mapping.sam > mapping.bam
samtools sort -o mapping.sort.bam mapping.bam
samtools index mapping.sort.bam
```
work_bwa2.sh (using pipe)  
```
#PBS -N bwa_pipe
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
bwa mem ../data/index/REL606.fa ../data/reads_1.fq.gz ../data/reads_2.fq.gz | \
 samtools view -b - | \
 samtools sort -o mapping.sort.2.bam -
samtools index mapping.sort.2.bam
```
### (二)Mapping the short reads to the reference genome using minimap2  

work_minimap2.sh  
```
#PBS -N minimap2
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
minimap2 -ax sr ../data/REL606.fa ../data/reads_1.fq.gz ../data/reads_2.fq.gz |\
 samtools view -b - |\
 samtools sort -o mapping.sort.mm.bam -
samtools index mapping.sort.mm.bam
```

### (三) Mapping the long reads to the reference genome using minimap2  
work_minimap_pb.sh  
```
#PBS -N minimap2_pb
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
minimap2 -ax map-pb ../data/REL606.fa ../data/REL606.fa ../data/pb_ecoli_0001.fastq |\
 samtools view -b - |\
 samtools sort -o mapping.sort.pb.bam -
samtools index mapping.sort.pb.bam
```
### (四) 显示和比较比对结果  
使用IGV查看比对结果  
![](./igv_snapshot.png) 


## 四、作业与思考  
1. 先组装，得到contigs，然后将contigs用bwa mem比对到参考基因组上  
2. 用igv显示比对结果   
 
```
# 选K值
kmergenie SRR098038.fastq.gz
# 组装
minia -in SRR098038.fastq.gz -kmer-size 23 -out SRR098038
# mapping
bwa mem ../data/index/REL606.fa SRR098038.contigs.fa | samtools view -b - | samtools sort -o contig_mapping.sort.bam -
# 结果文件为contig_mapping.sort.bam

```
## 五、参考资料  
[bwa](https://github.com/lh3/bwa)  
[minimap2](https://github.com/lh3/minimap2)  
[samtools](http://www.htslib.org/)  
[IGV](http://software.broadinstitute.org/software/igv/)  

