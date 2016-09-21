# 实验二 基因组比对  
## 一、实验目的  
1. 理解比对（mapping, alignment）的含义  
2. 理解全局比对和局部比对的区别和应用  
3. 掌握应用bwa, samtools的使用  
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
### 安装软件  
```
# 准备工作目录
mkdir lab2
cd lab2
mkdir data
mkdir soft
mkdir result

cd soft

# 安装bwa
git clone https://github.com/lh3/bwa.git
cd bwa 
make
./bwa

# 安装bowtie2
git clone https://github.com/BenLangmead/bowtie2.git
cd bowtie2
make
./bowtie2

# 安装samtools
curl -O -L https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar jxvf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
make
./samtools

# 安装igv
curl -O http://data.broadinstitute.org/igv/projects/downloads/IGV_2.3.81.zip
unzip IGV_2.3.81.zip
cd IGV_2.3.81
./igv.sh

```
>> 提示：上述软件已经下载到/bs1/data/genomeLab/lab2/soft/，可以直接拷到你的工作目录下。  

### 实验数据  
```
/bs1/data/genomeLab/lab2/data/REL606.fa
/bs1/data/genomeLab/lab2/data/SRR098038.fastq.gz
```
### Mapping and viewer  
```
# 准备数据和index参考基因组
cd data
ln -s /bs1/data/genomeLab/lab2/data/REL606.fa /bs1/data/genomeLab/lab2/data/SRR098038.fastq.gz ./
samtools faidx REL606.fa
mkdir index
cd index
ln -s ../REL606.fa ./
bwa index REL606.fa

cd ../../result
bwa aln ../data/index/REL606.fa ../data/SRR098038.fastq.gz  > SRR098038.sai
bwa samse ../data/index/REL606.fa SRR098038.sai ../data/SRR098038.fastq.gz > SRR098038.sam
samtools view -b SRR098038.sam > SRR098038.bam
samtools sort -o SRR098038.sort.bam SRR098038.bam
samtools index SRR098038.sort.bam

# 显示比对结果  
samtools tview SRR098038.sorted.bam ../data/REL606.fa

```
## 四、作业与思考  
1. 先组装，得到contigs，然后将contigs用bwa mem比对到参考基因组上  
2. 用igv查看比对情况  

## 五、参考资料  


