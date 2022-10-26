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

> [98 Free Whole Genome Assembly (WGA) Analysis Tools](https://bioinformaticshome.com/tools/wga/wga.html)

## 三、上机操作
### 在conda中重新安装R
```shell
$ source /opt/miniconda3/bin/activate
$ conda create -n r_env r-essentials r-base
```
### 如果不成功，需要自己先安装[miniconda](https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh)
```shell
# 按步骤输入yes
$ sh Miniconda3-py37_4.12.0-Linux-x86_64.sh
# 可能需要重启shell终端，在conda中安装R环境
$ conda create -n r_env r-essentials r-base
```

### 创建工作目录
```shell
# 每个用户home目录下限额使用10G硬盘，主要存放代码。做实验需用到大量数据，因此在专门路径中做实验
# 先按自己的学号建立专属文件夹
$ cd /data/stdata/genomic/grade2020
$ mkdir YourStudentID
$ cd # 回到home路径

# 建立工作路径的软链接
$ ln -s /data/stdata/genomic/grade2020/你的学号
$ cd YourStudentID

#新建一个目录lab1，本实验所有数据和输出都放入该目录中  
$ mkdir lab1
$ cd lab1
$ mkdir data
$ mkdir result
```

### 数据存放位置  
DNA测序数据位于：  
> [genomics_lab1_reads.fastq.gz](./genomics_lab1_reads.fastq.gz)  
> 【更靠谱的数据】/data/stdata/genomic/lab01/data/reads_1.fq.gz  
> 【更靠谱的数据】/data/stdata/genomic/lab01/data/reads_2.fq.gz  

参考基因组位于：  
> [genomics_lab1_ref.fa.gz](./genomics_lab1_ref.fa.gz)  
> 【更靠谱的数据】/data/stdata/genomic/lab01/data/ref.fa  

### 组装  
#### 准备数据  
```shell
$ cd data

# *** 如果使用本地上传的数据，用下面的命令 ***
$ gunzip genomics_lab1_reads.fastq.gz # 解压缩paired-end reads数据
$ for i in `seq 1 8 196904`; do let j=$i+3; sed -n "${i},${j}p" genomics_lab1_reads.fastq; done > reads_1.fq
$ for i in `seq 5 8 196904`; do let j=$i+3; sed -n "${i},${j}p" genomics_lab1_reads.fastq; done > reads_2.fq
$ gzip reads_1.fq reads_2.fq
$ rm -f genomics_lab1_reads.fastq
$ gunzip genomics_lab1_ref.fa.gz # 解压缩参考基因组数据
$ mv genomics_lab1_ref.fa ref.fa

# *** 如果使用集群上的数据，用下面的命令 ***
$ ln -s /data/stdata/genomic/lab01/data/reads_* ./
$ ln -s /data/stdata/genomic/lab01/data/ref.fa ./

$ cd ../result
```

#### 估算k值  
```shell
$ ls ../data/reads_* > reads.file
```

#### 准备kmergenie软件（如果未安装kmergenie）
> 1. 下载：[kmergenie](http://kmergenie.bx.psu.edu/)
> 2. 上传到集群并解压缩：
```shell
$ tar -zxfv kmergenie-1.7051.tar.gz
```
> 3. 进入kmergenie-1.7051文件夹并编译
```shell
$ cd kmergenie-1.7051
$ make
```
> 4. 将kmergenie加入到系统路径，方便调用（先回忆熟悉vi的使用）
```shell
$ vi ~/.bash_profile
PATH=$PATH:$HOME/software/kmergenie-1.7051
$ source ~/.bash_profile
```

新建一个脚本文件，work_kmer.sh，写入以下内容:  
```shell
#!/bin/bash
#$ -S /bin/bash
#$ -N kmer
#$ -cwd
#$ -j y
conda activate r_env
kmergenie reads.file
```

```shell
# 用qsub提交任务至计算节点
$ qsub work_kmer.sh
$ qstat # 查看自己用户名下的kmer任务是否运行起来（r状态）
```

> 结束后查看结果，选择最优k值27 (对应本地数据)  
> 结束后查看结果，选择最优k值111 (对应集群数据)  

#### 1. 用velvet组装
新建一个脚本文件，work_velvet.sh，写下下列内容:  
```shell
#!/bin/bash
#$ -S /bin/bash
#$ -N velvet
#$ -cwd
#$ -j y
velveth ecoli.velvet 27 -shortPaired -fastq.gz -separate ../data/reads_1.fq.gz ../data/reads_2.fq.gz # (若使用集群数据，请设置k为111)
velvetg ecoli.velvet -exp_cov auto
```

```shell
# 用qsub提交任务至计算节点
$ qsub work_velvet.sh
```

#### 2.用minia组装  
新建一个脚本文件，work_minia.sh，写入下列内容:  
```shell
#!/bin/bash
#$ -S /bin/bash
#$ -N minia
#$ -cwd
#$ -j y
minia -in ../data/reads_1.fq.gz,../data/reads_2.fq.gz -kmer-size 27 -out ecoli.minia # (若使用集群数据，请设置k为111)
```

```shell
# 用qsub提交任务至计算节点
$ qsub work_minia.sh
```

#### 3. 用SPAdes组装  
新建一个脚本文件，work_spades.sh，写入下列内容:  
```shell
#!/bin/bash
#$ -S /bin/bash
#$ -N spades
#$ -cwd
#$ -j y
source /opt/miniconda3/bin/activate
conda activate genomelab
spades.py -t 4 -1 ../data/reads_1.fq.gz -2 ../data/reads_2.fq.gz -o ecoli.spades
```

```shell
# 用qsub提交任务至计算节点
$ qsub work_spades.sh
```

#### 组装效果评价
[下载quast](https://sourceforge.net/projects/quast/)    
上传到集群home目录  

```shell
$ cd
$ tar -zxvf quast-5.0.2.tar.gz -C ./ # 解压缩quast-5.0.2.tar.gz
$ chmod a+x quast-5.0.2/quast.py # 为quast.py增加执行权限
$ rm -f quast-5.0.2.tar.gz
$ cd YourStudentID/lab1/result/ # 返回工作路径
$ ln -s ecoli.velvet/contigs.fa velvet.contigs.fa # 在当前路径建立组装结果文件的软链接，方便比较
$ ln -s ecoli.spades/scaffolds.fasta spades.scaffolds.fa

# 用quast评价组装结果
$ ~/quast-5.0.2/quast.py -R ../data/ref.fa velvet.contigs.fa ecoli.minia.contigs.fa spades.scaffolds.fa
```

#### 查看评价结果  
```shell
$ less quast_results/latest/report.txt 
```

## 四、作业  
1. 不同k-mer值对组装的影响，对velvet和minia用41和51进行组装，比较组装效果  
2. 熟悉和理解基因组组装一些术语名词，如N50, NG50, contig, scaffold, gap等
3. 理解k-mer频次分布图，如何根据k-mer频次分布图估算基因组大小及杂合度  
 
## 五、参考文献  
1. [velvet手册](./Velvet-Manual.pdf)
2. [mimia手册](https://github.com/ZhijunBioinf/minia)
3. [SPAdes手册](https://github.com/ZhijunBioinf/spades)

