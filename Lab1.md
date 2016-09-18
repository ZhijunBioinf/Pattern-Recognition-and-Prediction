# 基因组组装
## 两种组装策略  
   1. Overlap/layout/consensus
   2. De Bruijn k-mer graphs  

第一种策略主要应用在长reads组装上，如sanger测序数据和第三代测序数据，组装软件包括phrap, cap3等。  第二种策略主要应用于短reads数据组装上，包括velvet, soapdenovo, ABYSS等。      
**Overlap/layout/consensus基本步骤**  
> 1. Calculate all overlaps. 计算重叠片断  
> 2. Cluster based on overlap. 重叠片断聚类  
> 3. Do a multiple sequence alignment. 多序列比对,取一致序列  
** De Bruijin k-mer graphs基本步骤 **  
> 1. Building the k-mer graph  
> 2. Construct contigs  
> 3. Scaffolding and fill gaps  

