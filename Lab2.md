# 实验二 基因组比对  
## 一、实验目的  
1. 理解比对（mapping, alignment）的含义  
2. 理解全局比对和局部比对的区别和应用  
3. 掌握应用bwa, samtools的使用  
4. 理解SAM, BAM文件格式  

## 二、先导知识  

### 比对的两种策略  
1. Global alignment
2. Local alignment


我们熟悉的blast和blat均属于第二类。   
另外，不同长度的reads比对所用的策略也不一样，对于短reads，基于local alignment的软件如blast, blat不适合。  
将短的reads回帖到长的参考基因组上，这一过程称之为mapping。一般reads数目很大，读长短，参考基因组较长，对于mapping软件有两个要求：
1. 速度
2. 准确性
3. 

Mapping软件众多，比较有名的包括bwa, soap, bowtie, novoalign  
另外，由于真核生物mRNA不含有内含子，与一般的DNA mapping软件要求不一样，故转录组mapping使用的软件也不一样，转录组mapping软件比较有名的包括：STAR, hisat  
本实验主要介绍一般意义上的DNA mapping软件的使用。  

