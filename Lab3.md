# 实验三 基因组注释  
## 一、实验目的  
1. 理解基因组注释，了解原核生物与真核生物基因结构的差异
2. 掌握原核生物基因组注释方法(prokka)
3. 了解真核生物基因组注释方法（maker）
4. 了解常见基因组注释文件格式，如gff, bed


## 二、知识回顾与要求  
生物的遗传信息物质是DNA，DNA是基因的载体，基因是生物行使功能的单位，支持着生命的基本构造和性能，将基因组所有基因找出来是基因组注释的第一步。目前一般基因组项目一般流程是首先组装得到基因组草图(draft)，然后对草图进行基因预测和基因功能预测，即所谓的基因组注释。基因组注释结果的好坏会直接影响后续的分析，所以基因组注释对于基因组项目非常关键。   
原核生物与真核生物由于基因结构不同，基因预测方法也不一样。真核生物基因组注释比较复杂，一般由基因组中心或相关专业人员完成。原核生物基因组注释相对比较简单，已有较成熟的基因组注释软件。    

### 原核生物与真核生物基因结构差异
> 原核生物：不含内含子 -》RNA与DNA序列一致  
> 真核生物：含有内含子

### 实验要求  
1. 掌握用prokka注释原核生物基因组  
2. 掌握用maker注释真核生物基因组  

## 三、上机操作  
### 数据存放位置  
```
/data/lab/genomic/lab03/data/
```

### 数据及工作目录准备  
```
mkdir lab03
cd lab03
ln -s /data/lab/genomic/lab03/data/ ./
mkdir results

```

## (一) 原核生物基因组注释--prokka    
```
cd ../results

```

work_prokka.sh  
```
#PBS -N prokka
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
prokka --outdir anno --prefix PROKKA ../data/REL606.fa

```
注释结果存放在anno目录中，查看结果，了解基因组注释常见的几种格式。  

## （二）真核生物基因组溈--maker  
```
# create control files for maker
$ maker -CTL
```
会产生4个参数设置文件：  
-rw-rw-r-- 1 wangys wangys  1479 Oct 22 08:22 maker_bopts.ctl  
-rw-rw-r-- 1 wangys wangys   893 Oct 22 08:22 maker_evm.ctl  
-rw-rw-r-- 1 wangys wangys  1488 Oct 22 08:22 maker_exe.ctl  
-rw-rw-r-- 1 wangys wangys  4765 Oct 22 08:22 maker_opts.ctl  

编辑maker_opts.ctl文件，改变以下几个参数，几他的用默认参数：  
genome=../data/dpp_contig.fasta
est=../data/dpp_est.fasta
protein=../data/dpp_protein.fasta
est2genome=1

work_maker.sh
```
#PBS -N maker
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
maker
```
真核生物基因组注释比较复杂，这里只是向大家介绍了maker的一般使用，如果要使用maker注释新的基因组，建议参阅：[http://gmod.org/wiki/MAKER_Tutorial](http://gmod.org/wiki/MAKER_Tutorial)  
查看结果文件：  
```
$ cd dpp_contig.maker.output/dpp_contig_datastore/05/1F/contig-dpp-500-500
$ ls -l
```
-rw-rw-r-- 1 wangys wangys 65341 Oct 22 08:30 contig-dpp-500-500.gff
-rw-rw-r-- 1 wangys wangys   717 Oct 22 08:30 contig-dpp-500-500.maker.proteins.fasta
-rw-rw-r-- 1 wangys wangys  4443 Oct 22 08:30 contig-dpp-500-500.maker.transcripts.fasta
-rw-rw-r-- 1 wangys wangys  4120 Oct 22 08:30 run.log
drwxrwxr-x 3 wangys wangys  4096 Oct 22 08:30 theVoid.contig-dpp-500-500

### 用Artemis查看注释结果（选做）  
这一部分是在本地台式机上完成。  
下载地址：http://www.sanger.ac.uk/science/tools/artemis  
将prokka注释得到的prokka.gff文件拷到本地电脑上(用scp)  
打开Artemis，装载注释结果  
>    1. Start Artemis  
>    2. Click OK  
>    3. Go to File -> Open File Manager  
>    4. Navigate to the folder  
>    5. Choose the prokka.gff file you copied from the server

## 四、作业与思考  
1. 尝试用其他原核生物基因预测软件进行基因预测，如[GLIMMER](http://ccb.jhu.edu/software/glimmer/index.shtml)，或者[GeneMark](http://topaz.gatech.edu/GeneMark/)，并比较不同软件注释结果有什么不同
2. 了解真核生物基因组注释软件，如[Augustus](http://bioinf.uni-greifswald.de/augustus/), [GlimmerHMM](http://ccb.jhu.edu/software/glimmerhmm/), [maker](http://www.yandell-lab.org/software/maker.html)

## 五、参考文献  
1. [prokka github](https://github.com/tseemann/prokka)
 
