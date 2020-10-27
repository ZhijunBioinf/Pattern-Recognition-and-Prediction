# 实验四：分类器之NB, SVC

## 实验目的
* 1）使用朴素贝叶斯(Naive Bayes, NB)完成剪接位点识别
* 2）使用支持向量分类(Support Vector Classification, SVC)完成剪接位点识别

## 1. 以NB进行剪接位点识别
参考程序：myNB.py
```python3
import numpy as np
from sklearn import naive_bayes # 导入NB包
import sys

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集，在命令行指定文件名
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

clf = naive_bayes.GaussianNB() # 创建一个NB的实例
trX = train[:,1:]
trY = train[:,0]
clf.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of NB: %g%%(%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))
```

```bash
# KNN分类器：在命令行指定训练集、测试集、近邻数K
python3 myNB.py EI_train.txt EI_test.txt
```

## 3. 以Logistic回归进行剪接位点识别
参考程序：myLR.py
```python3
import numpy as np
from sklearn import linear_model # 导入线性模型包
import sys

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

maxIterations = int(sys.argv[3]) # 在命令行指定最大迭代次数
clf = linear_model.LogisticRegression(max_iter=maxIterations) # 创建一个LR的实例
trX = train[:,1:]
trY = train[:,0]
clf.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of LR: %g%%(%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))
```

```bash
# LR分类器：在命令行指定训练集、测试集、迭代次数
python3 myLR.py EI_train.txt EI_test.txt 1000
```

## 4. 以Decision Tree进行剪接位点识别
参考程序：myDT.py
```python3
import numpy as np
from sklearn import tree # 导入Decision Trees包
import sys
import graphviz # 导入Graphviz包

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

clf = tree.DecisionTreeClassifier() # 创建一个DT的实例
trX = train[:,1:]
trY = train[:,0]
clf.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of DT: %g%%(%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))

# Export the tree in Graphviz format
graphFileName = sys.argv[3] # 从命令行指定图文件名称
dotData = tree.export_graphviz(clf, out_file=None)
graph = graphviz.Source(dotData)
graph.render(graphFileName)
print('The tree in Graphviz format is saved in "%s.pdf".' % graphFileName)
```

```bash
# 安装Graphviz绘图包
pip3 install --user graphviz -i https://pypi.tuna.tsinghua.edu.cn/simple
# DT分类器：在命令行指定训练集、测试集、DT图文件名
python3 myDT.py EI_train.txt EI_test.txt EISplicing_DecisionTreeGraph
```
[获得的DT树](https://github.com/ZhijunBioinf/Pattern-Recognition-and-Prediction/blob/master/Lab3_Classifiers_KNN-LR-DT/EISplicing_DecisionTreeGraph.pdf)

## 作业
1. 尽量看懂`参考程序`的每一行代码。
2. 参考程序kSpaceCoding_general.py中，供体位点序列的第71、72位保守二核苷酸GT是在程序中指定的，试着改写程序，实现从命令行传递`位置信息`给程序。
3. 参考程序getTrainTest.py中，测试集的比例testSize是在程序中指定的，试着改写程序，实现从命令行传递`划分比例`给程序。
4. 熟练使用sklearn包中的不同分类器。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* KNN手册：[sklearn.neighbors.KNeighborsClassifier](https://scikit-learn.org/stable/modules/neighbors.html#nearest-neighbors-classification)
* LR手册：[sklearn.linear_model.LogisticRegression](https://scikit-learn.org/stable/modules/linear_model.html#logistic-regression)
* DT手册：[sklearn.tree.DecisionTreeClassifier](https://scikit-learn.org/stable/modules/tree.html#classification)
