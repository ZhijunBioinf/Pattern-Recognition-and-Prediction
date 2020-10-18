# 实验三：分类器之KNN、Decision Tree、Naive Bayes

## 实验目的
* 1）使用K近邻（K-nearest neighbor, KNN）完成剪接位点识别
* 2）使用决策树（Decision Tree, DT）完成剪接位点识别
* 3）使用朴素贝叶斯（Naive Bayes）完成剪接位点识别

## 1. 训练集与测试集构建
* 1）编写通用的k-spaced碱基对组分特征表征程序
```bash
# 安装机器学习包sklearn
pip3 install --user sklearn -i https://pypi.tuna.tsinghua.edu.cn/simple
```

参考程序：kSpaceCoding_general.py, 通用程序避免了每次在程序中修改文件名和其他参数的麻烦。
```bash
# 获得供体真实位点序列表征结果：在命令行指定序列文件名'EI_true.seq'、输出结果文件名'EI_true_kSpace.txt'与KMAX值
python3 kSpaceCoding_general.py EI_true.seq EI_true_kSpace.txt 4
```

```python3
import numpy as np # 导入numpy包，并重命名为np
import sys # 导入sys包，用于从命令行传递参数给python程序

def file2matrix(filename, bpTable, KMAX=2): # 为KMAX提供默认参数（与实验二不同）
	fr = open(filename) # 打开文件
	arrayOLines = fr.readlines() # 读取所有内容
	fr.close() # 及时关闭文件

	numberOfLines = len(arrayOLines) # 得到文件行数
	returnMat = np.zeros((numberOfLines, 16*(KMAX+1))) # 为返回的结果矩阵开辟内存
	lineNum = 0

	for line in arrayOLines:
		line = line.strip() # 删除空白符，包括行尾回车符
		listFromLine = line.split(': ') # 以': '为分隔符进行切片
		nt_seq = list(listFromLine[1]) # 取出核酸序列并转换成list
		del(nt_seq[70:72]) # 删除位于第71，72位的供体位点
		
		kSpaceVec = []
		for k in range(KMAX+1): # 计算不同k条件下的kSpace特征
			bpFreq = bpTable.copy() # bpTable是一个字典型变量，一定要用字典的copy函数，Python函数参数使用的址传递

			for m in range(len(nt_seq)-k-1): # 扫描序列，并计算不同碱基对的频率
				bpFreq[nt_seq[m]+nt_seq[m+1+k]] += 1 # 序列的子串会自动在字典中寻找对应的key，很神奇！否则要自己写if语句匹配
			bpFreqVal = list(bpFreq.values()) # 取出bpFreq中的值并转换成list
			kSpaceVec.extend(np.array(bpFreqVal)/(len(nt_seq)-k-1)) # 每个k下的特征，需除以查找的所有子串数

		returnMat[lineNum,:] = kSpaceVec
		lineNum += 1
	return returnMat, lineNum


if __name__ == '__main__':
	filename = sys.argv[1]
	outputFileName = sys.argv[2]
	KMAX = sys.argv[3]
	bpTable = {}
	for m in ('A','T','C','G'):
		for n in ('A','T','C','G'):
			bpTable[m+n] = 0

	kSpaceMat, SeqNum = file2matrix(filename, bpTable, KMAX)
	np.savetxt(outputFileName, kSpaceMat, fmt='%g', delimiter=',')
	print('The number of sequences is %d. Matrix of features is saved in %s' % (SeqNum, outputFileName))
	
```

## 作业
自己独立编写序列表征程序。不怕报错，犯错越多，进步越快！

## 参考文献
[1] Staden R. Computer methods to locate signals in nucleic acid sequences [J]. Nucleic Acids Research. 1984, 12(2):505. <br>
[2] Zhang M Q, Marr T G. A weight array method for splicing signal analysis [J]. Computer applications in the biosciences: CABIOS, 1993, 9(5):499-509. <br>
[3] Pertea M, Lin X Y, Salzberg S L. GeneSplicer: a new computational method for splice site prediction [J]. Nucleic Acids Research. 2001, 29:1185-1190. <br>
[4] Reese M G, Eeckman F H, Kulp D, et al. Improved splice site detection in Genie [J]. Journal of Computational Biology, 1997, 4(3):311-323. <br>
[5] Rogozin I B, Milanesi L. Analysis of donor splice signals in different organisms [J]. Journal of Molecular Evolution, 1997, 45(1):50-59. <br>
[6] Degroeve S, Saeys Y, Baets B D, et al. SpliceMachine: predicting splice sites from high-dimensional local context representations [J]. Bioinformatics. 2005,21:1332-1338. <br>
[7] Meher P K, Sahu T K, Rao A R. Prediction of donor splice sites using random forest with a new sequence encoding approach [J]. Biodata Mining. 2016,9:4. <br>
[8] Pollastro P, Rampone S. HS3D, a dataset of Homo sapiens splice regions and its extraction procedure from a major public database [J]. International Journal of Modern Physics C, 2002, 13(8):1105–1117. <br>
[9] Chen Y Z, Tang Y R, Sheng Z Y, et al. Prediction of mucin-type O-glycosylation sites in mammalian proteins using the composition of k-spaced amino acid pairs [J]. BMC Bioinformatics, 2008, 9(1):101-112.

## 致谢
剪接位点研究背景，部分摘自湖南农业大学博士学位论文《基于卡方决策表的分子序列信号位点预测》(2019)。<br>
感谢曾莹博士提供其学位论文！
