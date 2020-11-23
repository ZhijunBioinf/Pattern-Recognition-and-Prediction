# 实验八：无监督学习之聚类分析(K-Means, Hierarchical clustering, DBSCAN)

## 实验目的
* 1）使用K-means完成聚类分析。
* 2）使用Hierarchical clustering完成聚类分析。
* 3）使用DBSCAN完成聚类分析。
* 4）比较三种聚类分析方法的效果。

## 准备工作目录
```bash
$ mkdir lab_08
$ cd lab_08

# 若python3不可用，需先激活base环境
$ source /opt/miniconda3/bin/activate
$ conda activate
```

## 背景
* [无监督学习](https://en.wikipedia.org/wiki/Unsupervised_learning): is a type of machine learning that looks for previously undetected patterns in a data set with no pre-existing labels and with a minimum of human supervision. In contrast to supervised learning that usually makes use of human-labeled data, unsupervised learning, also known as self-organization allows for modeling of probability densities over inputs. It forms one of the three main categories of machine learning, along with supervised and reinforcement learning. Semi-supervised learning, a related variant, makes use of supervised and unsupervised techniques.
* [聚类分析](https://en.wikipedia.org/wiki/Cluster_analysis): is the task of grouping a set of objects in such a way that objects in the same group (called a cluster) are more similar (in some sense) to each other than to those in other groups (clusters). It is a main task of exploratory data mining, and a common technique for statistical data analysis, used in many fields, including pattern recognition, image analysis, information retrieval, bioinformatics, data compression, computer graphics and machine learning.
* 经典聚类模型包括但不限于：
> * Centroid models: [k-means](https://en.wikipedia.org/wiki/K-means_clustering), k-means++, [Mean Shift](https://scikit-learn.org/stable/modules/clustering.html#mean-shift), etc.
> * Connectivity models: [hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering), [BIRCH(Balanced Iterative Reducing and Clustering using Hierarchies)](https://en.wikipedia.org/wiki/BIRCH), etc.
> * Density models: [DBSCAN(Density-Based Spatial Clustering of Applications with Noise)](https://en.wikipedia.org/wiki/DBSCAN), [OPTICS(Ordering Points To Identify the Clustering Structure)](https://en.wikipedia.org/wiki/OPTICS), etc.

## 1. 使用K-means完成`手写数字`聚类分析
* 参考程序：myKMeansDigits.py
```python3
import sys
from time import time # 函数的计时包
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics # 聚类效果指标包
from sklearn.cluster import KMeans # KMeans包
from sklearn.datasets import load_digits #
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale


```

```bash
# 设置主成分累计解释的百分比为95%
$ python3 myPCA.py ACE_train.txt ACE_test.txt 0.95 train_pca.txt test_pca.txt
# 以SVR建模预测：规格化、rbf核、10次交叉寻优
$ python3 ../lab_06/mySVR.py train_pca.txt test_pca.txt 1 rbf 10 ObsdYvsPredY_PCA_SVR.pdf
```

## 2. 使用MI-filter进行特征选择，再以保留变量建立SVR模型
* 参考程序：myMIFilter.py
```python3
import sys
import numpy as np
from sklearn.feature_selection import SelectPercentile, mutual_info_regression # 导入MI for regression包

trainFile = sys.argv[1]
testFile = sys.argv[2]
train = np.loadtxt(trainFile, delimiter='\t') # 载入训练集
test = np.loadtxt(testFile, delimiter='\t') # 载入测试集
trX = train[:,1:]
trY = train[:,0]
teX = test[:,1:]
teY = test[:,0]

percentile = int(sys.argv[3]) # Percent of features to keep
selector = SelectPercentile(mutual_info_regression, percentile=percentile) # 创建一个基于MI的SelectPercentile实例
trX = selector.fit_transform(trX, trY)
teX = selector.transform(teX)

print('Number of principal components: %d' % trX.shape[1])
newTrainFile = sys.argv[4]
newTestFile = sys.argv[5]
np.savetxt(newTrainFile, np.hstack((trY.reshape(-1,1), trX)), fmt='%g', delimiter='\t') # 将Y与X以列组合后，保存到文件
np.savetxt(newTestFile, np.hstack((teY.reshape(-1,1), teX)), fmt='%g', delimiter='\t')
print('%d features are selected.' % trX.shape[1])
print('New training set is saved into: %s\nNew test set is saved into: %s' % (newTrainFile, newTestFile))
```

```bash
# 保留5%的特征
$ python3 myMIFilter.py ACE_train.txt ACE_test.txt 5 train_mi.txt test_mi.txt
# 以SVR建模预测：规格化、rbf核、10次交叉寻优
$ python3 ../lab_06/mySVR.py train_mi.txt test_mi.txt 1 rbf 10 ObsdYvsPredY_MI_SVR.pdf
```

## 3. 使用SVM-RFE进行特征选择，再以保留变量建立SVR模型
* 参考程序：mySVMRFE.py
```python3
import sys
import numpy as np
from sklearn.feature_selection import RFE # 导入RFE包
from sklearn.svm import SVR # 导入SVR包
from pytictoc import TicToc
from random import sample

trainFile = sys.argv[1]
testFile = sys.argv[2]
train = np.loadtxt(trainFile, delimiter='\t') # 载入训练集
test = np.loadtxt(testFile, delimiter='\t') # 载入测试集
trX = train[:,1:]
trY = train[:,0]
teX = test[:,1:]
teY = test[:,0]

numX = trX.shape[1]
if numX > 200:
    randVec = np.array(sample(range(numX), 100)) + 1 # 考虑到特征数较多，SVM运行时间较长，随机抽100个特征用于后续建模
    print('Note: 100 features are randomly selected to speed up modeling')
    trX = train[:,randVec]
    teX = test[:,randVec]
    
n_features = int(sys.argv[3])
estimator = SVR(kernel="linear")
t = TicToc()
t.tic()
selector = RFE(estimator, n_features_to_select=n_features, step=1) # 创建一个RFECV实例，每次删1个特征
trX = selector.fit_transform(trX, trY)
teX = selector.transform(teX)
print('Time cost in selecting fetures with SVM-RFE: %gs' % t.tocvalue())

newTrainFile = sys.argv[4]
newTestFile = sys.argv[5]
np.savetxt(newTrainFile, np.hstack((trY.reshape(-1,1), trX)), fmt='%g', delimiter='\t') # 将Y与X以列组合后，保存到文件
np.savetxt(newTestFile, np.hstack((teY.reshape(-1,1), teX)), fmt='%g', delimiter='\t')
print('%d features are selected.' % trX.shape[1])
print('New training set is saved into: %s\nNew test set is saved into: %s' % (newTrainFile, newTestFile))
```

* SVM-RFE运行时间较长，将命令写到脚本中再用qsub提交任务。Note: 命令脚本中也需要先激活base环境，才可用python3
* work_mySVMRFE.sh
```bash
#!/bin/bash
#$ -S /bin/bash
#$ -N mySVMRFE
#$ -j y
#$ -cwd

# 激活base环境，保证计算节点上正常使用python3
source /opt/miniconda3/bin/activate
conda activate

# 设置保留10个特征
python3 mySVMRFE.py ACE_train.txt ACE_test.txt 10 train_svmrfe.txt test_svmrfe.txt
```

```bash
# 以SVR建模预测：规格化、rbf核、10次交叉寻优
$ python3 ../lab_06/mySVR.py train_svmrfe.txt test_svmrfe.txt 1 rbf 10 ObsdYvsPredY_RFE_SVR.pdf
```

## 4. 使用Random Forest进行特征选择，再以保留变量建立SVC模型
```bash
# RF只能用于分类问题，对lab_03路径中的剪接位点训练集与测试集建立软链接
$ ln -s ../lab_03/EI_train.txt
$ ln -s ../lab_03/EI_test.txt
```

* 参考程序：myRandomForest.py
```python3
import sys
import numpy as np
from sklearn.feature_selection import SelectFromModel # 导入SelectFromModel包
from sklearn.ensemble import RandomForestClassifier # 导入RF包
from pytictoc import TicToc

trainFile = sys.argv[1]
testFile = sys.argv[2]
train = np.loadtxt(trainFile, delimiter=',') # 载入训练集
test = np.loadtxt(testFile, delimiter=',') # 载入测试集
trX = train[:,1:]
trY = train[:,0]
teX = test[:,1:]
teY = test[:,0]

n_features = sys.argv[3]
t = TicToc()
t.tic()
clf = RandomForestClassifier(max_depth=2, random_state=0, max_features=n_features) # 创建一个RF实例
clf = clf.fit(trX, trY)
selector = SelectFromModel(clf, prefit=True) # 创建一个SelectFromModel实例
trX = selector.transform(trX)
teX = selector.transform(teX)
print('Time cost in selecting fetures with Random-Forest: %gs' % t.tocvalue())

newTrainFile = sys.argv[4]
newTestFile = sys.argv[5]
np.savetxt(newTrainFile, np.hstack((trY.reshape(-1,1), trX)), fmt='%g', delimiter=',') # 将Y与X以列组合后，保存到文件
np.savetxt(newTestFile, np.hstack((teY.reshape(-1,1), teX)), fmt='%g', delimiter=',')
print('%d features are selected.' % trX.shape[1])
print('New training set is saved into: %s\nNew test set is saved into: %s' % (newTrainFile, newTestFile))
```

```bash
# max_features (looking for the best split): set as 'auto'
$ python3 myRandomForest.py EI_train.txt EI_test.txt auto EI_train_rf.txt EI_test_rf.txt
# 以SVR建模预测：规格化、rbf核、10次交叉寻优
$ python3 ../lab_04/mySVC.py EI_train_rf.txt EI_test_rf.txt 1 rbf 1 10
```


## 作业
1. 尽量看懂`参考程序`的每一行代码。 <br>
2. 熟练使用sklearn包中的不同特征降维/选择方法。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* PCA手册：[sklearn.decomposition.PCA](https://scikit-learn.org/stable/modules/decomposition.html#principal-component-analysis-pca)
* MI手册：[sklearn.feature_selection.mutual_info_regression](https://scikit-learn.org/stable/modules/feature_selection.html#univariate-feature-selection)
* RFE手册：[sklearn.feature_selection.RFE](https://scikit-learn.org/stable/modules/feature_selection.html#recursive-feature-elimination)
* RandomForest: [sklearn.ensemble.RandomForestClassifier](https://scikit-learn.org/stable/modules/feature_selection.html#tree-based-feature-selection)
