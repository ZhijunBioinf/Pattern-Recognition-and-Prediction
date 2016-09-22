# 实验三 基因组注释  
## 一、实验目的  
1. 理解基因组注释，了解原核生物与真核生物基因结构的差异
2. 掌握原核生物基因组注释方法


## 二、知识回顾  
生物的遗传信息物质是DNA，DNA是基因的载体，基因是生物行使功能的单位，支持着生命的基本构造和性能，将基因组所有基因找出来是基因组注释的第一步。目前一般基因组项目一般流程是首先组装得到基因组草图(draft)，然后对草图进行基因预测和基因功能预测，即所谓的基因组注释。基因组注释结果的好坏会直接影响后续的分析，所以基因组注释对于基因组项目非常关键。   
原核生物与真核生物由于基因结构不同，基因预测方法也不一样。真核生物基因组注释比较复杂，一般由基因组中心或相关专业人员完成。原核生物基因组注释相对比较简单，已有较成熟的基因组注释软件。本实验主要介绍原核生物基因组注释软件prokka的使用。    

### 原核生物与真核生物基因结构差异
> 原核生物：不含内含子 -》RNA与DNA序列一致  
> 真核生物：含有内含子


## 三、上机操作  
### 数据准备
```
# 基因组数据  
/bs1/data/genomeLab/lab3/data/contigs.fasta

# 新建工作目录
mkdir lab3
cd lab3
mkdir data
mkdir soft
mkdir result

```
### 安装软件  
```
cd soft
git clone https://github.com/tseemann/prokka.git
cd 
prokka --setupdb
prokka --version

```
### Run Prokka on the contigs  
```
cd ../result
ln -s ../data/contigs.fasta ./
prokka --outdir anno --prefix prokka contigs.fasta
cat ./anno/prokka.txt

```
### 用Artemis查看注释结果  
这一部分是在本地台式机上完成。  
下载地址：http://www.sanger.ac.uk/science/tools/artemis  
将prokka注释得到的gff文件拷到本地电脑上，用scp  
打开Artemis，装载注释结果  
>    1. Start Artemis  
>    2. Click OK  
>    3. Go to File -> Open File Manager  
>    4. Navigate to the ~/Downloads folder  
>    5. Choose the prokka.gff file yoiu copied from Amazon

## 四、作业与思考  


## 五、参考文献  
