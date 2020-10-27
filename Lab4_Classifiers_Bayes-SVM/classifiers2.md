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
# NB分类器：在命令行指定训练集、测试集
python3 myNB.py EI_train.txt EI_test.txt
```

## 2. 以SVC进行剪接位点识别
[不同核函数的SVC](https://scikit-learn.org/stable/_images/sphx_glr_plot_iris_svc_0011.png)
参考程序：mySVC.py
```python3
import numpy as np
from sklearn import svm # 导入svm包
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
# SVC分类器：在命令行指定训练集、测试集
python3 mySVC.py EI_train.txt EI_test.txt 1000
```


## 作业
1. 尽量看懂`参考程序`的每一行代码。
2. 熟练使用sklearn包中的不同分类器。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* NB手册：[sklearn.naive_bayes.GaussianNB](https://scikit-learn.org/stable/modules/naive_bayes.html#gaussian-naive-bayes)
* SVM手册：[sklearn.svm.SVC](https://scikit-learn.org/stable/modules/svm.html#classification)
