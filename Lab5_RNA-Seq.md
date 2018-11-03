# 实验五 RNA-Seq分析   
## 一、实验目的  
1. 了解生物的RNA种类，真核生物mRNA的特点及如何分离纯化
2. 了解RNA、DNA比对的区别
3. 掌握DESeq2分析DE基因方法

## 二、知识回顾  
转录组测序（RNA-Seq）应用非常广泛，目前测序市场有超过一半做的是转录组。  
转录组没序是对生物体内所有的RNA进行测序，一般可以按测序鞋标不同而分为  
mRNA，microRNA, total RNA测序，传统意义上的RNA-Seq指的是mRNA测序。一般  
RNA-Seq项目分析包括以下几个步骤：

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
### 1. Mapping  
work_mapping.sh  
```
#PBS -N nucmer
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
for i in $(seq 291 302);
do 
hisat2 -p 1 \
-x /data/lab/genomic/lab05/ref/index/osa \
-q /data/lab/genomic/lab05/data/SRR7760${i}.fastq \
-S SRR7760$i.sam >SRR7760${i}.log;
done
```

## 四、作业与思考  

## 五、参考文献  

