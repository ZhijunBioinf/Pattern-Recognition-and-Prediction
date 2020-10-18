# 实验三：分类器之KNN、Decision Tree、Naive Bayes

## 实验目的
* 1）提取供体真实位点与虚假位点序列的k-space组分特征；构建训练集与测试集
* 2）使用K近邻（K-Nearest Neighbor, KNN）完成剪接位点识别
* 3）使用决策树（Decision Tree, DT）完成剪接位点识别
* 4）使用朴素贝叶斯（Naive Bayes, NB）完成剪接位点识别

## 1. 训练集与测试集构建
* 1）编写更好用的k-spaced碱基对组分特征表征程序（用于HS3D数据的供体真实/虚假位点序列表征）<br>
参考程序：kSpaceCoding_general.py, 该程序避免了每次在程序中修改文件名和其他参数的麻烦。
```python3
import numpy as np # 导入numpy包，并重命名为np
import sys # 导入sys包，其中的argv函数用于从命令行传递参数给python程序

def file2matrix(filename, bpTable, KMAX=2): # 为KMAX提供默认参数(updated)
	fr = open(filename) # 打开文件
	arrayOLines = fr.readlines() # 读取所有内容
	del(arrayOLines[:4]) # 删除头4行（updated, 避免了运行程序之前，另外使用sed删除头4行）
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
				sub_str = nt_seq[m]+nt_seq[m+1+k] # 提出子串(updated)
				if sub_str in bpFreq.keys(): # 如果子串在bpFreq中有对应的key，才统计频次(updated, NOTE:在供体虚假位点序列中存在非正常碱基)
					bpFreq[sub_str] += 1 # 序列的子串会自动在字典中寻找对应的key，很神奇！否则要自己写if语句匹配
			bpFreqVal = list(bpFreq.values()) # 取出bpFreq中的值并转换成list
			kSpaceVec.extend(np.array(bpFreqVal)/(len(nt_seq)-k-1)) # 每个k下的特征，需除以查找的所有子串数

		returnMat[lineNum,:] = kSpaceVec
		lineNum += 1
	return returnMat, lineNum

if __name__ == '__main__':
	filename = sys.argv[1]
	outputFileName = sys.argv[2]
	KMAX = int(sys.argv[3])
	bpTable = {}
	for m in ('A','T','C','G'):
		for n in ('A','T','C','G'):
			bpTable[m+n] = 0
	
	kSpaceMat, SeqNum = file2matrix(filename, bpTable, KMAX)
	np.savetxt(outputFileName, kSpaceMat, fmt='%g', delimiter=',')
	print('The number of sequences is %d. Matrix of features is saved in %s' % (SeqNum, outputFileName))
	
```

```bash
# 建立供体位点序列文件的软链接
ln -s ../lab_02/EI_true.seq ../lab_02/EI_false.seq ./
# 获得供体真实位点序列表征结果：在命令行指定序列文件名为'EI_true.seq'，输出结果文件名为'EI_true_kSpace.txt'，KMAX值为4
python3 kSpaceCoding_general.py EI_true.seq EI_true_kSpace.txt 4
# 获得供体虚假位点序列表征结果：在命令行指定序列文件名为'EI_false.seq'，输出结果文件名为'EI_false_kSpace.txt'，KMAX值为4
python3 kSpaceCoding_general.py EI_false.seq EI_false_kSpace.txt 4
```

* 2）以序列表征文件构建训练集、测试集 <br>
参考程序：getTrainTest.py
```python3
import numpy as np
import sys
from random import sample # 导入sample函数，用于从虚假位点数据中随机抽取样本
from sklearn.model_selection import train_test_split # 用于产生训练集、测试集

trueSiteFileName = sys.argv[1]
falseSiteFileName = sys.argv[2]
trueSitesData = np.loadtxt(trueSiteFileName, delimiter = ',') # 载入true位点数据
numOfTrue = len(trueSitesData)
falseSitesData = np.loadtxt(falseSiteFileName, delimiter = ',') # 载入false位点数据
randVec = sample(range(numOfTrue), len(trueSitesData)) # 随机产生true位点样本个数的随机向量
falseSitesData = falseSitesData[randVec,] # 以随机向量从false位点数据中抽取样本

Data = np.vstack((trueSitesData, falseSitesData)) # 按行将true位点与false位点数据组合
Y = np.vstack((np.ones((numOfTrue,1)),np.zeros((numOfTrue,1)))) # 产生Y列向量
testSize = 0.3 # 测试集30%，训练集70%
X_train, X_test, y_train, y_test = train_test_split(Data, Y, test_size = testSize, random_state = 0)

trainingSetFileName = sys.argv[3]
testSetFileName = sys.argv[4]
np.savetxt(trainingSetFileName, np.hstack((y_train, X_train)), fmt='%g', delimiter=',') # 将Y与X以列组合后，保存到文件
np.savetxt(testSetFileName, np.hstack((y_test, X_test)), fmt='%g', delimiter=',')
print('Generate training set(%d%%) and test set(%d%%): Done!' % ((1-testSize)*100, testSize*100))

```

```bash
# 构建训练集与测试集：在命令行指定true位点数据、false位点数据、train文件、test文件
python3 getTrainTest.py EI_true_kSpace.txt EI_false_kSpace.txt EI_train.txt EI_test.txt
```

## 2. 以KNN进行剪接位点识别
```bash
# 安装机器学习包sklearn
pip3 install --user sklearn -i https://pypi.tuna.tsinghua.edu.cn/simple
```


## 作业
1. 尽量看懂`参考程序`的每一行代码。
2. 参考程序中，供体位点序列的第71、72位保守二核苷酸GT是在程序中指定的，试着改写程序，实现从命令行传递`位置信息`给程序。
3. 参考程序中，测试集的比例是在程序中指定的，试着改写程序，实现从命令行传递`划分比例`给程序。
4. 熟练使用sklearn包中的不同分类器。
不怕报错，犯错越多，进步越快！~_~

