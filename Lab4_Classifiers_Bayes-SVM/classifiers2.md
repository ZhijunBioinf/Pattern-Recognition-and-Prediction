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
from sklearn import preprocessing # 导入数据预处理包
from sklearn.model_selection import GridSearchCV # 导入参数寻优包
from random import sample

train = np.loadtxt(sys.argv[1], delimiter=',') # 载入训练集
test = np.loadtxt(sys.argv[2], delimiter=',') # 载入测试集

train = train[sample(range(len(train)), 200),1:]
trX = train[:,1:]
trY = train[:,0]
teX = test[:,1:]
teY = test[:,0]

if int(sys.argv[3])==1: # 建模前，是否将每个特征归一化到[-1,1]
    max_abs_scaler = preprocessing.MaxAbsScaler()
    trX = max_abs_scaler.fit_transform(trX)
    teX = max_abs_scaler.transform(teX)

if int(sys.argv[4])==1: # 是否寻找最优参数：c, g
    C_range = np.power(2, np.arange(-5,15,2.0)) # 指定C的范围
    gamma_range = np.power(2, np.arange(3,-15,-2.0)) # 指定g的范围
    parameters = dict(gamma=gamma_range, C=C_range) # 将c, g组成字典，用于参数的grid遍历
    numOfFolds = int(sys.argv[5]) # 命令行指定寻优过程的交叉验证次数
        
    clf = svm.SVC() # 创建一个SVC的实例
    grid = GridSearchCV(clf, param_grid=parameters, cv=numOfFolds) # 创建一个GridSearchCV实例
    grid.fit(trX, trY) # grid寻优c, g
    print("The best parameters are %s with a score of %g" % (grid.best_params_, grid.best_score_))
    clf = svm.SVC(C=grid.best_params_['C'], gamma=grid.best_params_['gamma'])
else:
    clf = svm.SVC()
        
clf.fit(trX, trY) # 训练模型
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of SVC: %g%%(%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))
```

```bash
# SVC分类器：在命令行指定训练集、测试集
python3 mySVC.py EI_train.txt EI_test.txt 1 1 5
```


## 作业
1. 尽量看懂`参考程序`的每一行代码。
2. 熟练使用sklearn包中的NB, SVC分类器。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* NB手册：[sklearn.naive_bayes.GaussianNB](https://scikit-learn.org/stable/modules/naive_bayes.html#gaussian-naive-bayes)
* SVM手册：[sklearn.svm.SVC](https://scikit-learn.org/stable/modules/svm.html#classification)
