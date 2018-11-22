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

head -n 5000 mouse.1.protein.faa > mouse.1_sub5k.faa
head -n 5000 mouse.2.protein.faa > mouse.2.sub5k.faa
```
blast1.sh  
```
#!/bin/bash
#PBS -S /bin/bash
#PBS -N Blast1
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
blastp -query mouse.1_sub5k.faa -db zebrafish.1.protein.faa -out mouse.1.zebrafish.txt -outfmt 6
blastp -query mouse.2_sub5k.faa -db zebrafish.1.protein.faa -out mouse.2.zebrafish.txt -outfmt 6
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
构建基因家族可以使用[OrthoMCL](http://orthomcl.org/orthomcl/)，本实验使用的方法与OrthoMCL类似，目的是为了让大家更清楚背后的原理。  

目前在GenBank中有66个菌株的基因组序列已经释放，[https://www.ncbi.nlm.nih.gov/genome/?term=txid55080[Organism:exp]](https://www.ncbi.nlm.nih.gov/genome/?term=txid55080[Organism:exp])，但在RefSeq中只有61个菌株有基因组序列，我们需要对这61个菌株的蛋白进行聚类分析，构建基因家族  
![](https://micans.org/mcl/img/fa75.png) . 

2.1 数据准备：  
请完成以下表格，收集基因组信息：[https://docs.qq.com/sheet/DUEZiWFBEcktGTWRO](https://docs.qq.com/sheet/DUEZiWFBEcktGTWRO)  

```
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/612/125/GCF_000612125.1_BreAgr1.0/GCF_000612125.1_BreAgr1.0_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/328/345/GCF_000328345.1_Brevibacillus_agri_BAB-2500_v3.1/GCF_000328345.1_Brevibacillus_agri_BAB-2500_v3.1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/885/GCF_003710885.1_ASM371088v1/GCF_003710885.1_ASM371088v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/612/185/GCF_000612185.1_BreBor1.0/GCF_000612185.1_BreBor1.0_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/353/565/GCF_000353565.1_ASM35356v1/GCF_000353565.1_ASM35356v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/785/GCF_000738785.1_ASM73878v1/GCF_000738785.1_ASM73878v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/865/GCF_003710865.1_ASM371086v1/GCF_003710865.1_ASM371086v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/161/835/GCF_002161835.1_ASM216183v1/GCF_002161835.1_ASM216183v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/385/915/GCF_003385915.1_ASM338591v1/GCF_003385915.1_ASM338591v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/039/275/GCF_001039275.2_ASM103927v2/GCF_001039275.2_ASM103927v2_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/255/GCF_000346255.1_FJAT-GLX-0809/GCF_000346255.1_FJAT-GLX-0809_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/649/505/GCF_001649505.1_ASM164950v1/GCF_001649505.1_ASM164950v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/165/GCF_000010165.1_ASM1016v1/GCF_000010165.1_ASM1016v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/748/185/GCF_001748185.1_ASM174818v1/GCF_001748185.1_ASM174818v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/012/835/GCF_003012835.1_ASM301283v1/GCF_003012835.1_ASM301283v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/296/715/GCF_000296715.2_ASM29671v2/GCF_000296715.2_ASM29671v2_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/815/GCF_003710815.1_ASM371081v1/GCF_003710815.1_ASM371081v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/420/695/GCF_001420695.1_ASM142069v1/GCF_001420695.1_ASM142069v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/825/GCF_003710825.1_ASM371082v1/GCF_003710825.1_ASM371082v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/012/775/GCF_001012775.1_ASM101277v1/GCF_001012775.1_ASM101277v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/215/075/GCF_002215075.1_ASM221507v1/GCF_002215075.1_ASM221507v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/405/GCF_003013405.1_ASM301340v1/GCF_003013405.1_ASM301340v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/935/GCF_003710935.1_ASM371093v1/GCF_003710935.1_ASM371093v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/915/GCF_003710915.1_ASM371091v1/GCF_003710915.1_ASM371091v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/733/515/GCF_000733515.1_ASM73351v1/GCF_000733515.1_ASM73351v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/926/995/GCF_002926995.1_ASM292699v1/GCF_002926995.1_ASM292699v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/927/075/GCF_002927075.1_ASM292707v1/GCF_002927075.1_ASM292707v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/927/085/GCF_002927085.1_ASM292708v1/GCF_002927085.1_ASM292708v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/663/745/GCF_003663745.1_ASM366374v1/GCF_003663745.1_ASM366374v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/594/765/GCF_003594765.1_ASM359476v1/GCF_003594765.1_ASM359476v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/374/385/GCF_000374385.1_ASM37438v1/GCF_000374385.1_ASM37438v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/706/795/GCF_002706795.1_ASM270679v1/GCF_002706795.1_ASM270679v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/590/075/GCF_003590075.1_ASM359007v1/GCF_003590075.1_ASM359007v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/237/005/GCF_000237005.1_ASM23700v2/GCF_000237005.1_ASM23700v2_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/412/145/GCF_002412145.1_A5-miseq/GCF_002412145.1_A5-miseq_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/219/535/GCF_000219535.2_ASM21953v3/GCF_000219535.2_ASM21953v3_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/265/735/GCF_003265735.1_ASM326573v1/GCF_003265735.1_ASM326573v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/259/955/GCF_002259955.1_ASM225995v1/GCF_002259955.1_ASM225995v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/472/325/GCF_000472325.2_PE_F_paired_/GCF_000472325.2_PE_F_paired__protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/696/705/GCF_001696705.1_ASM169670v1/GCF_001696705.1_ASM169670v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/865/525/GCF_002865525.1_ASM286552v1/GCF_002865525.1_ASM286552v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/311/785/GCF_000311785.1_ASM31178v1/GCF_000311785.1_ASM31178v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/965/GCF_003710965.1_ASM371096v1/GCF_003710965.1_ASM371096v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/985/GCF_003710985.1_ASM371098v1/GCF_003710985.1_ASM371098v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/503/775/GCF_000503775.1_version_1_for_the_Brevibacillus_panacihumi_W25_genome/GCF_000503775.1_version_1_for_the_Brevibacillus_panacihumi_W25_genome_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/619/605/GCF_001619605.1_ASM161960v1/GCF_001619605.1_ASM161960v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/710/905/GCF_003710905.1_ASM371090v1/GCF_003710905.1_ASM371090v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/187/725/GCF_001187725.1_ASM118772v1/GCF_001187725.1_ASM118772v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/897/295/GCF_002897295.1_Brl_assembly01/GCF_002897295.1_Brl_assembly01_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/282/075/GCF_000282075.1_Brevibacillus.strBC25_v1.0/GCF_000282075.1_Brevibacillus.strBC25_v1.0_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/282/015/GCF_000282015.1_Brevibacillus.strCF112_v1.0/GCF_000282015.1_Brevibacillus.strCF112_v1.0_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/423/865/GCF_001423865.2_ASM142386v2/GCF_001423865.2_ASM142386v2_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/475/GCF_003013475.1_ASM301347v1/GCF_003013475.1_ASM301347v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/395/GCF_003013395.1_ASM301339v1/GCF_003013395.1_ASM301339v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/013/445/GCF_003013445.1_ASM301344v1/GCF_003013445.1_ASM301344v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/114/075/GCF_900114075.1_IMG-taxon_2687453682_annotated_assembly/GCF_900114075.1_IMG-taxon_2687453682_annotated_assembly_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/645/205/GCF_001645205.1_ASM164520v1/GCF_001645205.1_ASM164520v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/400/265/GCF_003400265.1_ASM340026v1/GCF_003400265.1_ASM340026v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/673/705/GCF_001673705.1_ASM167370v1/GCF_001673705.1_ASM167370v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/454/065/GCF_000454065.1_ASM45406v1/GCF_000454065.1_ASM45406v1_protein.faa.gz
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/744/635/GCF_000744635.1_ASM74463v1/GCF_000744635.1_ASM74463v1_protein.faa.gz
gunzip *.gz
```
1. 对蛋白序列进行all-to-all blast  
做Blast之前请改下序列名，在各自序列名后面加上GCF编号，如将```WP_003333770.1```改成```WP_003333770.1:GCF_000010165```，将所有蛋白序列合并到一个文件```all_pro.faa```，建索引```makeblastdb -in all_pro.faa -dbtype prot```。  
blastAll.sh
```
#!/bin/bash
#PBS -S /bin/bash
#PBS -N blast_all
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR
blastp -query all_pro.faa -db all_pro.faa -out allBlast.tsv -outfmt 6 -evalue 1e-10
```
2. 提取每个hit的score值，构建一个表征两条序列的相似性的特征值  

```
cut -f 1,2,12 allBlast.tsv > allBlast.abc
```
3. 用mcl进行聚类  
work_mcl.sh
```
#!/bin/bash
#PBS -S /bin/bash
#PBS -N MCL
#PBS -l nodes=1:ppn=1
#PBS -j oe
cd $PBS_O_WORKDIR

mcxload -abc allBlast.abc --stream-mirror -write-tab data.tab -o data.mci
mcl data.mci -I 1.4
mcl data.mci -I 2
mcl data.mci -I 4
mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40
clm dist --chain out.data.mci.I{14,20,40}
```
统计有多少基因家族，每个基因家族中每个菌株基因数。  

# 3. Reference  
1. [MCL](https://micans.org/mcl)  
2. [OrthoMCL](http://orthomcl.org/orthomcl)  
