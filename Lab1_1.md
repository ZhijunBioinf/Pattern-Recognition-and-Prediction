# 实验一 基因组组装
## 一、实验目的  
1. 熟悉使用高性能计算机集群
2. 熟悉基因组从头组装原理及步骤  
3. 掌握velvet, minia, SPAdes等短序列拼装软件使用 
4. 熟悉用quast评价组装效果  

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

## 三、上机操作  

### 创建工作目录  
```
#新建一个目录lab1，本实验所有数据和输出都放入该目录中  
$ mkdir lab1
$ cd lab1
$ mkdir data
$ mkdir result
```

### 数据存放位置  
```
数据存放在服务器位置：
/data/lab/genomic/lab01/data/reads_1.fq.gz
/data/lab/genomic/lab01/data/reads_2.fq.gz
#参考基因组
/data/lab/genomic/lab01/data/ref.fa
```
### 组装  
#### 准备数据  
```
$ cd data
$ ln -s /data/lab/genomic/lab01/data/reads_1.fq.gz /data/lab/genomic/lab01/data/reads_2.fq.gz ./
$ cd ../result
```
#### 估算k值  
```
$ ls ../data/reads_* > reads.file
```
新建一个脚本文件，work_kmer.sh，写入以下内容:  
```
#!/bin/bash
#$ -S /bin/bash
#$ -N kmer
#$ -cwd
#$ -j y
/opt/bio/kmergenie-1.7048/kmergenie reads.file
```

```
# 用qsub提交任务至计算节点
$ qsub work_kmer.sh
```
结束后查看结果，选择最优k值57  

#### 1. 用velvet组装
新建一个脚本文件，work_velvet.sh，写下下列内容:  
```
#!/bin/bash
#$ -S /bin/bash
#$ -N velvet
#$ -cwd
#$ -j y
velveth ecoli.velvet 57 -shortPaired -fastq.gz -separate ../data/reads_1.fq.gz ../data/reads_2.fq.gz
velvetg ecoli.velvet -exp_cov auto
```

```
# 用qsub提交任务至计算节点
$ qsub work_velvet.sh
```

#### 2.用minia组装  
新建一个脚本文件，work_minia.sh，写入下列内容:  
```
#!/bin/bash
#$ -S /bin/bash
#$ -N minia
#$ -cwd
#$ -j y
/opt/bio/bin/minia -in ../data/reads_1.fq.gz,../data/reads_2.fq.gz -kmer-size 57 -out ecoli.minia
```

```
# 用qsub提交任务至计算节点
$ qsub work_minia.sh
```

#### 3. 用SPAdes组装  
新建一个脚本文件，work_spades.sh，写入下列内容:  
```
#!/bin/bash
#$ -S /bin/bash
#$ -N spades
#$ -cwd
#$ -j y
spades.py -t 4 -1 ../data/reads_1.fq.gz -2 ../data/reads_2.fq.gz -o ecoli.spades
```

```
# 用qsub提交任务至计算节点
$ qsub work_spades.sh
```

#### 组装效果评价  
```
# 可直接在命令行执行
$ quast -R /data/lab/genomic/lab01/data/ref.fa \
   ecoli.velvet/contigs.fa \
   ecoli.minia.contigs.fa \
   ecoli.spades/scaffolds.fasta
```
#### 查看评价结果  
```
$ less quast_results/latest/report.txt 
```

## 四、作业  
1. 不同k-mer值对组装的影响，对velvet和minia用31和57进行组装，比较组装效果  
2. 熟悉和理解基因组组装一些术语名词，如N50, NG50, contig, scaffold, gap等
3. 理解k-mer频次分布图，如何根据k-mer频次分布图估算基因组大小及杂合度  
 
## 五、参考文献  
1. [velvet手册](https://github.com/dzerbino/velvet/blob/master/Manual.pdf)
2. [mimia手册](https://github.com/GATB/minia#introduction)
3. [SPAdes手册](https://github.com/ablab/spades)

