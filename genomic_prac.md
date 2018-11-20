# 基因组学教学实习  


## 一、本地Blast  

1. 准备数据  
数据已经下载，放在```/data/lab/genomic/prac/blast```目录中。  
```
curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.2.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.1.protein.faa.gz

gunzip *.faa.gz
```
2. 建索引  
```
makeblastdb -in zebrafish.1.protein.faa -dbtype prot
```
3. 运行blastp  
我们先取2条序列试一下  
```
head -n 11 mouse.1.protein.faa > mm-first.faa
blastp -query mm-first.faa -db zebrafish.1.protein.faa
blastp -query mm-first.faa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt
blastp -query mm-first.faa -db zebrafish.1.protein.faa -outfmt 6
less mm-first.x.zebrafish.txt
```
blast1.sh  
```
#!/bin/bash
#PBS -S /bin/bash
#PBS -N Blast1
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
blastp -query mouse.1.protein.faa -db zebrafish.1.protein.faa -out mouse.1.zebrafish.txt -outfmt 6
blastp -query mouse.2.protein.faa -db zebrafish.1.protein.faa -out mouse.2.zebrafish.txt -outfmt 6
```
4. Visualizing BLAST score distributions in RStudio  
```
blast_out1 <- read.table('mouse.1.zebrafish.txt', sep='\t')
blast_out2 <- read.table('mouse.2.zebrafish.txt', sep='\t')
colnames(blast_out1) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(blast_out2) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
hist(blast_out1$evalue)
hist(blast_out2$evalue)
hist(blast_out1$bitscore) 
hist(blast_out2$bitscore) 
plot(blast_out1$pident, blast_out1$bitscore)
plot(blast_out2$pident, blast_out2$bitscore)
plot(blast_out1$pident  * (blast_out1$qend - blast_out1$qstart), blast_out1$bitscore)
```
思考题：如果只统计the best hsp evalue，要如何改？  

## 二、根据blast结果对蛋白序列进行聚类 -- 构建基因家族  
任务：构建Brevibacillus基因家族  
目前在GenBank RefSeq中有66个菌株的基因组序列已经释放，[https://www.ncbi.nlm.nih.gov/genome/?term=txid55080[Organism:exp]](https://www.ncbi.nlm.nih.gov/genome/?term=txid55080[Organism:exp])，我们需要对这66个菌株的蛋白进行聚类分析，构建基因家族  

2.1 数据准备：  
请完成以下表格，收集基因组信息：[https://docs.qq.com/sheet/DUEZiWFBEcktGTWRO](https://docs.qq.com/sheet/DUEZiWFBEcktGTWRO)  

```
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/282/015/GCF_000282015.1_Brevibacillus.strCF112_v1.0/GCF_000282015.1_Brevibacillus.strCF112_v1.0_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/282/075/GCF_000282075.1_Brevibacillus.strBC25_v1.0/GCF_000282075.1_Brevibacillus.strBC25_v1.0_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/645/205/GCF_001645205.1_ASM164520v1/GCF_001645205.1_ASM164520v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/673/705/GCF_001673705.1_ASM167370v1/GCF_001673705.1_ASM167370v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/395/GCF_003013395.1_ASM301339v1/GCF_003013395.1_ASM301339v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/445/GCF_003013445.1_ASM301344v1/GCF_003013445.1_ASM301344v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/475/GCF_003013475.1_ASM301347v1/GCF_003013475.1_ASM301347v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/400/265/GCF_003400265.1_ASM340026v1/GCF_003400265.1_ASM340026v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/423/865/GCF_001423865.2_ASM142386v2/GCF_001423865.2_ASM142386v2_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/114/075/GCF_900114075.1_IMG-taxon_2687453682_annotated_assembly/GCF_900114075.1_IMG-taxon_2687453682_annotated_assembly_protein.faa.gz
.
.
.
gunzip *.gz
```
1. 对蛋白序列进行all-to-all blast  

```

```
2. 提取每个hit的e-value或score值，构建一个表征两条序列的相似性的特征值  
3. 用mcl进行聚类  


