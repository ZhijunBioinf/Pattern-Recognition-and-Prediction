# 实验五：序列表征/数值化2(以定量构效关系建模为例)

## 实验目的
* 1）了解定量构效关系建模的研究背景
* 2）编程实现肽序列的AA531(531 properties of Amino Acids)特征表征/数值化

## 1. 定量构效关系
* 分子是物质的基本组成单位。分子结构属性决定其生理活性。
* 通过统计学、信息学方法从分子结构中提取、总结分子结构的信息与规律，有助于从理论上指导实验过程。
* 定量构效关系(Quantitative Structure-Acitivity Relationship, QSAR)是以分子的基本理化性质与相应的生理活性为基础，通过数学或统计学手段定量研究有机小分子（如抑制剂等）与生物大分子（如受体、酶等）间的相互作用。
* QSAR可用于高效生物活性分子化合物筛选、药物的环境毒性评价、新药设计与合成等[1-2]。

![](./SchematicOfGeneralApproachInQSAR.jpg)

## 2. ACE抑制剂研究现状
* 血管紧张素转化酶(Angiotensin-Converting Enzyme，ACE)抑制肽是一类从食源性蛋白质中分离得到的具有降高血压活性的多肽。由于其降血压效果好，而且没有降压药物的毒副作用从而引起了广泛关注。
* 近年来，ACE抑制肽的构效关系研究成为研究重点。结构生物信息学研究表明，ACE抑制肽的ACE抑制能力不仅与其分子质量有关，而且与其氨基酸序列以及其立体空间构象之间存在高度相关性。
* ACE抑制肽的抑制类型与ACE抑制活性、构效关系也存在一定相关性。对ACE抑制肽构效关系进行深入研究将有助于指导开发高活性的功能性食品及降血压药物[3]。

## 3. 数据集
* ACE_tri-peptides_150数据集的肽序列及其活性(lg(IC<sub>50</sub>))来源于文献[4].
* [ACEtriPeptidesSequencesActivities.txt](https://github.com/ZhijunBioinf/Pattern-Recognition-and-Prediction/blob/master/Lab5_PeptideSequencesCoding/ACEtriPeptidesSequencesActivities.txt)

## 4. 基于[AA531](http://www.genome.jp/aaindex)的序列表征
* An amino acid index is a set of 20 numerical values representing any of the different physicochemical and biological properties of amino acids. The AAindex1 section of the Amino Acid Index Database is a collection of published indices together with the result of cluster analysis using the correlation coefficient as the distance between two indices. This section currently contains 544 indices.
* 某些氨基酸缺失了部分理化属性，预处理后，符合条件的生理生化属性共531个。
* 对于每条ACE三肽序列，以每个位置氨基酸对应的531个理化属性依次替换序列，可获得531x3 = 1593个特征。
* [AA531properties.txt](https://github.com/ZhijunBioinf/Pattern-Recognition-and-Prediction/blob/master/Lab5_PeptideSequencesCoding/AA531properties.txt)

## 5. 工作目录准备
```sh
# 建立lab_05文件夹
$ mkdir lab_05
$ cd lab_05
```

## 6. 序列表征
* 参考程序：AA531Coding.py
```python3
import numpy as np
import sys

def file2matrix(filename, seqLength, AA531Dict):
    fr = open(filename) # 打开文件
    arrayOLines = fr.readlines() # 读取所有内容
    fr.close() # 及时关闭文件

    numberOfLines = len(arrayOLines) # 得到文件行数
    returnMat = np.zeros((numberOfLines, 531*seqLength)) # 为返回的结果矩阵开辟内存
    Y = np.zeros((numberOfLines, 1))
    lineNum = 0

    for line in arrayOLines:
        line = line.strip() # 删除空白符，包括行尾回车符
        listFromLine = line.split('\t') # 以'\t'为分隔符进行切片
        AASeq = listFromLine[0] # 取出氨基酸序列
        Y[lineNum] = float(listFromLine[1]) # 取出活性值Y
        
        feaVec = []
        for AA in AASeq: # 扫描序列，将每个氨基酸替换为相应的531个理化属性
            if AA in AA531Dict.keys(): # 如果序列中的氨基酸在AA531Dict中有对应的key，才进行替换
                feaVec.extend(AA531Dict[AA])
            else: # 否则以0替换
                print('Warning: nonregular amino acid found! Coding "%s" in "%s"(seqId: %d) with 531 zeros.' % (AA, AASeq, lineNum))
                feaVec.extend([0.0]*531)
                Y[lineNum] = -1

        returnMat[lineNum,:] = np.array(feaVec)
        lineNum += 1
    return Y, returnMat, lineNum

def makeAA531Dict(filename):
    fr = open(filename) # 打开文件
    arrayOLines = fr.readlines() # 读取所有内容
    del(arrayOLines[0]) # 删除head行
    fr.close() # 及时关闭文件

    AA531Dict = {}
    for line in arrayOLines:
        line = line.strip()
        listFromLine = line.split('\t')
        AA = listFromLine[0]
        properties = [float(i) for i in listFromLine[1:]] # 从文件读取的数值默认是字符串类型，需要转换为浮点型
        AA531Dict[AA] = properties
    return AA531Dict


if __name__ == '__main__':
    AASeqFileName = sys.argv[1]
    AA531FileName = sys.argv[2]
    seqLength = int(sys.argv[3])
    outputFileName = sys.argv[4]

    AA531Dict = makeAA531Dict(AA531FileName)
    Y, AA531Mat, SeqNum = file2matrix(AASeqFileName, seqLength, AA531Dict)
    np.savetxt(outputFileName, np.hstack((Y, AA531Mat)), fmt='%g', delimiter='\t')
    print('The number of sequences is %d. Matrix of features is saved in %s' % (SeqNum, outputFileName))

```
```sh
# 运行程序AA531Coding.py实现以AA531表征序列。在命令行指定ACE抑制剂三肽序列文件名、AA531属性文件名、序列长度、输出文件名。
$ python3 AA531Coding.py ACEtriPeptidesSequencesActivities.txt AA531properties.txt 3 result.txt
```

## 作业
自己独立编写序列表征程序。不怕报错，犯错越多，进步越快！

## 参考文献
[1] 代志军. 特征选择与样本选择用于癌分类与药物构效关系研究(湖南农业大学博士学位论文). 2014. <br>
[2] Nongonierma AB, FitzGerald RJ. Learnings from quantitative structure–activity relationship (QSAR) studies with respect to food protein-derived bioactive peptides: a review. RSC advances. 2016, 6(79): 75400-75413. <br>
[3] 王晓丹, 薛璐, 胡志和, 等. ACE抑制肽构效关系研究进展[J]. 食品科学, 2017, 38(5): 305-310. <br>
[4] 刘静, 彭剑秋, 管骁. 基于多元线性回归的血管紧张素转化酶抑制肽定量构效关系建模研究[J]. 分析科学学报, 2012, 28(1): 16-22. <br>

