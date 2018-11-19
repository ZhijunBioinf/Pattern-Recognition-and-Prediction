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

二、
