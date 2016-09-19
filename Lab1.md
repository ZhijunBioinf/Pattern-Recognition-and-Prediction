# 实验一 基因组组装
## 一、实验目的  
1. 熟悉基因组从头组装原理及步骤  
2. 掌握velvet, minia, SPAdes等短序列拼装软件使用 
3. 熟悉用quast评价组装效果  

## 二、基因组组装组装原理与方法  
### 两种策略  
   1. Overlap/layout/consensus
   2. De Bruijn k-mer graphs  

第一种策略主要应用在长reads组装上，如sanger测序数据和第三代测序数据，组装软件包括phrap, cap3等。  第二种策略主要应用于短reads数据组装上，包括velvet, soapdenovo, ABYSS等。  

**Overlap/layout/consensus基本步骤**  
> 1. Calculate all overlaps. 计算重叠片断  
> 2. Cluster based on overlap. 重叠片断聚类  
> 3. Do a multiple sequence alignment. 多序列比对,取一致序列  

**De Bruijin k-mer graphs基本步骤**  
> 1. Building the k-mer graph，构建k-mer图，有的软件会在之前加一步，对reads进行纠错  
> 2. Construct contigs，搜索路径  
> 3. Scaffolding and fill gaps，构建scaffold并填洞  

本实验主要介绍使用第二种策略的组装软件velvet, minia和SPAdes的使用  

### 如何选择合适k  
1. 多试几个k，看组装效果
2. 利用如[KmerGenie](http://kmergenie.bx.psu.edu/)进行辅助选择  
 
### 如何选择组装软件  
1. 选择你熟悉的软件  
2. 选择大家使用多的软件
3. 选择适合项目性质的软件，如基因组大小，染色体倍性，杂合度，数据类型，宏基因组，转录组等

> [浏览软件库](https://omictools.com/genome-assembly-category)

## 二、上机操作  
**登录服务器**
```
# Linux机器上登录
ssh -l public 172.28.137.55
# Windows机器上登录用putty客户端

```
### 软件下载与编译
组装软件：`velvet`, `minia`, `SPAdes`  
评价软件：`quast`
```
#新建一个目录lab1，本实验所有数据和输出都放入该目录中  
mkdir lab1
cd lab1
mkdir data
mkdir soft
mkdir result
cd soft

# 安装velvet
git clone https://github.com/manogenome/velvet.git
cd velvet
make
./velveth
./velvetg

# 安装minia
wget -c https://github.com/GATB/minia/releases/download/v2.0.7/minia-v2.0.7-Source.tar.gz
tar zxvf minia-v2.0.7-Source.tar.gz
cd minia-v2.0.7-Source
mkdir build
cd build
cmake ..
make -j8
./bin/minia

# 安装SPAdes
 wget -c http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0-Linux.tar.gz
 tar zxvf SPAdes-3.9.0-Linux.tar.gz
 cd SPAdes-3.9.0-Linux
 # 测试SPAdes
 ./bin/spades.py --test

# 安装quast
curl -O -L https://downloads.sourceforge.net/project/quast/quast-2.3.tar.gz
tar xzf quast-2.3.tar.gz

# 安装kmergenie
wget -c http://kmergenie.bx.psu.edu/kmergenie-1.7016.tar.gz
tar zxvf kmergenie-1.7016.tar.gz
cd kmergenie-1.7016
make
./kmergenie
```

### 数据下载  
```
数据存放在服务器位置：
/bs1/data/genomeLab/lab1/data/reads_1.fq.gz
/bs1/data/genomeLab/lab1/data/reads_2.fq.gz
#参考基因组
/bs1/data/genomeLab/lab1/data/ref.fa
```
### 组装  
```
# 准备数据
cd data
ln -s /bs1/data/genomeLab/lab1/data/reads_1.fq.gz /bs1/data/genomeLab/lab1/data/reads_2.fq.gz ./

cd ../result
[path to] velveth ecoli.velvet 21 -shortPaired -fastq.gz -separate ../data/reads_1.fq.gz ../data/reads_2.fq.gz
[path to] velvetg ecoli.velvet -exp_cov auto

[path to] minia -in ../data/reads_1.fq.gz,../data/reads_2.fq.gz -kmer-size 21 -out ecoli.minia

[path to] spades.py -t 4 -1 ../data/reads_1.fq.gz -2 ../data/reads_2.fq.gz -o ecoli.spades

#组装效果评价
[path to] quast.py -R ../data/ref.fa ecoli.velvet/contigs.fa ecoli.minia.contigs.fa ecoli.spades/scaffolds.fasta

#查看评价结果
less quast_results/latest/report.txt 

```
###作业  
1. 不同k-mer值对组装的影响，尝试用kmergenie辅助选择k  
2. 熟悉和理解基因组组装一些术语名词，如N50, NG50, contig, scaffold, gap等
3. 理解k-mer频次分布图，如何估算根据k-mer频次分布图估算基因组大小  


