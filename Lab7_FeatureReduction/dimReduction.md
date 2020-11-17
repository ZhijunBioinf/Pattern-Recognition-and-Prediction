# 实验七：特征降维/选择(PCA, MI-filter, SVM-RFE, RF)

## 实验目的
* 1）数据：[实验六](https://github.com/ZhijunBioinf/Pattern-Recognition-and-Prediction/blob/master/Lab6_Regression_MLR-PLSR-SVR/regress1.md)ACE抑制剂的训练集与测试集；基本模型：SVR
* 2）使用主成分分析(Principal Component Analysis, PCA)进行特征压缩降维，再以SVR建模预测，对比`实验六`的预测结果。
* 3）使用基于互信息的单变量过滤法(Mutual Information-based Filter)进行特征选择，后续同(2)。
* 4）使用基于SVM的迭代特征剔除(SVM Recursive Feature Elimination, SVM-RFE)进行特征选择，后续同(2)。
* 5）使用随机森林(Random Forest, RF)进行特征选择，后续同(2)。

## 准备工作目录与数据
```bash
$ mkdir lab_07
$ cd lab_07
# 对lab_06路径中的ACE抑制剂训练集与测试集建立软链接
$ ln -s ../lab_06/ACE_train.txt
$ ln -s ../lab_06/ACE_test.txt
```

## 1. 使用PCA进行特征压缩降维，再以保留主成分建立SVR模型
* 参考程序：myPCA.py
```python3
import sys
import numpy as np
from sklearn.decomposition import PCA # 导入PCA包

trainFile = sys.argv[1]
testFile = sys.argv[2]
train = np.loadtxt(trainFile, delimiter='\t') # 载入训练集
test = np.loadtxt(testFile, delimiter='\t') # 载入测试集
trX = train[:,1:]
trY = train[:,0]
teX = test[:,1:]
teY = test[:,0]

percentVar = float(sys.argv[3]) # 主成分累计解释的百分比
pca = PCA(n_components = percentVar, svd_solver = 'full') # 创建一个PCA实例
trX = pca.fit_transform(trX)
teX = pca.transform(teX)

newTrainFile = sys.argv[4]
newTestFile = sys.argv[5]
np.savetxt(newTrainFile, np.hstack((trY.reshape(-1,1), trX)), fmt='%g', delimiter='\t') # 将Y与X以列组合后，保存到文件
np.savetxt(newTestFile, np.hstack((teY.reshape(-1,1), teX)), fmt='%g', delimiter='\t')
print('New training set is saved into: %s\nNew test set is saved into: %s' % (newTrainFile, newTestFile))
```

```bash
# 设置主成分累计解释的百分比为95%
$ python3 myPCA.py ACE_train.txt ACE_test.txt 0.95 train_pca.txt test_pca.txt
# 以SVR建模预测：规格化、rbf核、10次交叉寻优
$ python3 ../lab_06/mySVR.py train_pca.txt test_pca.txt 1 rbf 10 ObsdYvsPredY_PCA_SVR.pdf
```

## 2. 使用MI-filter进行特征选择，以保留变量建立SVR模型
* 参考程序：myMIFilter.py
```python3

```

* SVM运行时间较长，将命令写到脚本中再用qsub提交任务
* work_mySVR.sh
```bash
#!/bin/bash
#$ -S /bin/bash
#$ -N mySVR
#$ -j y
#$ -cwd

# 1. 数据规格化，径向基核，10次交叉寻优
echo '------ scale: 1; kernel: rbf; numOfCV: 10 --------'
python3 myPCA_SVR.py ACE_train.txt ACE_test.txt 1 rbf 10 -1 ObsdYvsPredY_allX_SVR.pdf
echo

# 2. 数据规格化，径向基核，10次交叉寻优，主成分累计解释的百分比
echo '------ scale: 1; kernel: rbf; numOfCV: 10; percentVar: 0.95 --------'
python3 myPCA_SVR.py ACE_train.txt ACE_test.txt 1 rbf 10 0.95 ObsdYvsPredY_PCA_SVR.pdf
echo

# 规格化，线性核，不参数寻优
echo '------ scale: 1; kernel: linear; numOfCV: 0 --------'
python3 mySVR.py ACE_train.txt ACE_test.txt 1 linear 0 ObsdYvsPredY_SVR2.pdf
echo

# 规格化，径向基核(rbf)，10次交叉寻优
echo '------ scale: 1; kernel: rbf; numOfCV: 10 --------'
python3 mySVR.py ACE_train.txt ACE_test.txt 1 rbf 10 ObsdYvsPredY_SVR3.pdf
echo

# 规格化，径向基核(rbf)，不参数寻优
echo '------ scale: 1; kernel: rbf; numOfCV: 0 --------'
python3 mySVR.py ACE_train.txt ACE_test.txt 1 rbf 0 ObsdYvsPredY_SVR4.pdf
echo
```
```
# qsub提交任务
$ qsub work_mySVR.sh
```

* 尝试更多的选项搭配，看精度变化，比如数据不规格化时，各种核函数、是否寻优、不同交叉验证次数等情形的预测精度。

## 作业
1. 尽量看懂`参考程序`的每一行代码。 <br>
2. 熟练使用sklearn包中的不同特征降维/选择方法。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* MLR手册：[sklearn.linear_model.LinearRegression](https://scikit-learn.org/stable/modules/linear_model.html#ordinary-least-squares)
* PLSR手册：[sklearn.cross_decomposition.PLSRegression](https://scikit-learn.org/stable/modules/cross_decomposition.html)
* SVR手册：[sklearn.svm.SVR](https://scikit-learn.org/stable/modules/svm.html#regression)
