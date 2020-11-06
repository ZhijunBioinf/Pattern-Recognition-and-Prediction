# 实验五 RNA-Seq分析   
## 一、实验目的  
1. 了解生物的RNA种类，真核生物mRNA的特点及如何分离纯化
2. 了解RNA、DNA比对的区别
3. 掌握DESeq2分析差异表达基因(Differentially Expressed Genes, DEGs)方法

## 二、知识回顾  
* 转录组测序(RNA-Seq)应用非常广泛，目前测序市场有超过一半做的是转录组。  
* 转录组测序是对生物体内所有的RNA进行测序，一般可以按测序目标不同而分为mRNA, microRNA, total RNA测序，传统意义上的RNA-Seq指的是mRNA测序。
* 一般RNA-Seq项目分析包括以下几个步骤：

> 1. 序列质控
> 2. 将得到的clean reads 回帖到参考基因组
> 3. 计算每个基因的reads数
> 4. 比较不同组之间的DE基因
> 5. DE基因的功能挖掘
 
本实验主要介绍DESeq2的使用，DESeq2是用来做差异基因分析的一套R软件包，  


### 材料背景  
有两个水稻品种材料，品种BR-IRGA 409耐热，品种IRGA 428不耐热，拟比较这两个品种经热处理后基因表达差异。  

| RUN	| Class	| Cultivar |	Status |	Group	| Description |
| --- | --- | --- | --- | --- | --- |
| SRR7760302	| IRGA409_T	| BR-IRGA 409	| tolerant to heat stress	| BR-IRGA 409 treated_spikelet	| high temperature treatment (38C) for 7 hrs |
| SRR7760301	| IRGA409_T	| BR-IRGA 409	| tolerant to heat stress	| BR-IRGA 409 treated_spikelet	| high temperature treatment (38C) for 7 hrs |
| SRR7760300	| IRGA409_T	| BR-IRGA 409	| tolerant to heat stress	| BR-IRGA 409 treated_spikelet	| high temperature treatment (38C) for 7 hrs |
| SRR7760299	| IRGA409_C	| BR-IRGA 409	| tolerant to heat stress	| BR-IRGA 409 control_spikelet	| control condition (29C) for 7 hrs |
| SRR7760298	| IRGA409_C	| BR-IRGA 409	| tolerant to heat stress	| BR-IRGA 409 control_spikelet	| control condition (29C) for 7 hrs |
| SRR7760297	| IRGA409_C	| BR-IRGA 409	| tolerant to heat stress	| BR-IRGA 409 control_spikelet	| control condition (29C) for 7 hrs |
| SRR7760296	| IRGA428_T	| IRGA 428	| sensitive to heat stress	| IRGA 428 treated_spikelet	| high temperature treatment (38C) for 7 hrs |
| SRR7760295	| IRGA428_T	| IRGA 428	| sensitive to heat stress	| IRGA 428 treated_spikelet	| high temperature treatment (38C) for 7 hrs |
| SRR7760294	| IRGA428_T	| IRGA 428	| sensitive to heat stress	| IRGA 428 treated_spikelet	| high temperature treatment (38C) for 7 hrs |
| SRR7760293	| IRGA428_C	| IRGA 428	| sensitive to heat stress	| IRGA 428 control_spikelet	| control condition (29C) for 7 hrs |
| SRR7760292	| IRGA428_C	| IRGA 428	| sensitive to heat stress	| IRGA 428 control_spikelet	| control condition (29C) for 7 hrs |
| SRR7760291	| IRGA428_C	| IRGA 428	| sensitive to heat stress	| IRGA 428 control_spikelet	| control condition (29C) for 7 hrs |


## 三、上机操作  
### 进入genomelab环境
```
$ source /opt/miniconda3/bin/activate
$ conda activate genomelab
```
### 数据存放位置及工作目录准备  
```
Data: /data/lab/genomic/lab05/data
Ref: /data/lab/genomic/lab05/ref

$ mkdir lab5
$ cd lab5

# 建立data和ref的软链接
$ ln -s /data/lab/genomic/lab05/data/
$ ln -s /data/lab/genomic/lab05/ref/

$ mkdir results
$ cd results
```

### 1. Mapping(较慢，运行约8小时)
work_mapping.sh  
```
#!/bin/bash
#$ -S /bin/bash
#$ -N hisat2
#$ -j y
#$ -cwd

source /opt/miniconda3/bin/activate
conda activate genomelab
for i in $(seq 291 302)
do 
 hisat2 -p 1 -x ../ref/index/osa -q ../data/SRR7760${i}.fastq | \
 samtools view -b - | \
 samtools sort -o SRR7760${i}.sort.bam - > SRR7760${i}.log
done
```

```
# 提交任务
$ qsub work_mapping.sh
```

### 2. Count（运行约1小时）
work_count.sh  
```
#!/bin/bash
#$ -S /bin/bash
#$ -N count
#$ -j y
#$ -cwd

source /opt/miniconda3/bin/activate
conda activate genomelab
TPMCalculator -g ../ref/Oryza_sativa.IRGSP-1.0.41.gtf -d ./ -a
```

```
# 提交任务
$ qsub work_count.sh
```

### 3. Merge the counting matrix（直接运行paste, sed命令）
```
$ paste <(cut -f 1,6 SRR7760291.sort_genes.out) \
<(cut -f 6 SRR7760292.sort_genes.out) \
<(cut -f 6 SRR7760293.sort_genes.out) \
<(cut -f 6 SRR7760294.sort_genes.out) \
<(cut -f 6 SRR7760295.sort_genes.out) \
<(cut -f 6 SRR7760296.sort_genes.out) \
<(cut -f 6 SRR7760297.sort_genes.out) \
<(cut -f 6 SRR7760298.sort_genes.out) \
<(cut -f 6 SRR7760299.sort_genes.out) \
<(cut -f 6 SRR7760300.sort_genes.out) \
<(cut -f 6 SRR7760301.sort_genes.out) \
<(cut -f 6 SRR7760302.sort_genes.out) > counts.tsv

# Replace the header line
$ sed -i '1c\GeneID\tSRR7760291\tSRR7760292\tSRR7760293\tSRR7760294\tSRR7760295\tSRR7760296\tSRR7760297\tSRR7760298\tSRR7760299\tSRR7760300\tSRR7760301\tSRR7760302' counts.tsv 
```

### 4. DE analysis(使用R包实现)
* 如果在genomelab环境，其中的R没有安装以下需要的包。
* 可以退出到登录集群时的默认环境，其中的R版本已经安装好tidyverse, DESeq2包，不需另外安装
```
$ conda activate
# 进入R
$ R
```
```R
if(require(tidyverse) == FALSE)
  install.packages("tidyverse", repos = "https://mirror.lzu.edu.cn/CRAN/")
if(require(DESeq2) == FALSE)
  install.packages("BiocManager", repos = "https://mirror.lzu.edu.cn/CRAN/")
  BiocManager::install("DESeq2")

require(tidyverse)
require(DESeq2)
cts <- read_tsv("counts.tsv") %>% as.data.frame()
rownames(cts) <- cts$GeneID
cts <- cts[,-1]
colData <- data.frame(SampleID=c("SRR7760291","SRR7760292","SRR7760293",
                                 "SRR7760294","SRR7760295","SRR7760296",
                                 "SRR7760297","SRR7760298","SRR7760299",
                                 "SRR7760300","SRR7760301","SRR7760302"), 
                      Group=factor(c(rep("IRGA428_C",3),
                              rep("IRGA428_T",3),
                              rep("IRGA409_C",3),
                              rep("IRGA409_T",3))))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ Group)
dds <- DESeq(dds)
res.IRGA428 <- results(dds, contrast=c("Group","IRGA428_C","IRGA428_T"))
res.IRGA409 <- results(dds, contrast=c("Group","IRGA409_C","IRGA409_T"))
write.csv(res.IRGA428[order(res.IRGA428$pvalue),],"Results_428.csv")
write.csv(res.IRGA409[order(res.IRGA409$pvalue),],"Results_409.csv")
```

## 四、作业与思考  
比较两种水稻品种经热处理后基因表达变化情况有什么差异。  

## 五、参考文献  
1. [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) 
2. [TPMcalculator](https://github.com/ncbi/TPMCalculator)  
3. [HISAT2](http://daehwankimlab.github.io/hisat2/manual/)
