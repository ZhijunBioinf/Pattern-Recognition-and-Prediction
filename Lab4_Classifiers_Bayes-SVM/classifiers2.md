# 实验四：分类器之NB, SVC

## 实验目的
* 1）使用朴素贝叶斯(Naive Bayes, NB)完成剪接位点识别
* 2）使用支持向量分类(Support Vector Classification, SVC)完成剪接位点识别

## 准备工作目录
```
$ mkdir lab_04
$ cd lab_04
# 建立lab_03路径中供体位点训练集、测试集文件的软链接
$ ln -s ../lab_03/EI_train.txt ../lab_03/EI_test.txt ./

# 集群上若python3不可用，需先激活base环境
$ source /opt/miniconda3/bin/activate
$ conda activate
```

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
$ python3 myNB.py EI_train.txt EI_test.txt
```

## 2. 以SVC进行剪接位点识别
[不同核函数的SVC](https://scikit-learn.org/stable/_images/sphx_glr_plot_iris_svc_0011.png) <br>
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

train = train[sample(range(len(train)), 200),] # 考虑到SVM运行时间较长，从train中随机抽200样本用于后续建模
trX = train[:,1:]
trY = train[:,0]
teX = test[:,1:]
teY = test[:,0]

isScale = int(sys.argv[3]) # 建模前，是否将每个特征归一化到[-1,1]
kernelFunction = sys.argv[4] # {‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’}, default=’rbf’
isChooseCG = int(sys.argv[5]) # 是否寻找最优参数：c, g

if isScale:
    min_max_scaler = preprocessing.MinMaxScaler(feature_range=(-1,1))
    trX = min_max_scaler.fit_transform(trX)
    teX = min_max_scaler.transform(teX)

if isChooseCG:
    numOfFolds = int(sys.argv[6]) # 命令行指定寻优过程的交叉验证次数
    C_range = np.power(2, np.arange(-5,15,2.0)) # 指定C的范围
    gamma_range = np.power(2, np.arange(3,-15,-2.0)) # 指定g的范围
    parameters = dict(gamma=gamma_range, C=C_range) # 将c, g组成字典，用于参数的grid遍历
    
    clf = svm.SVC(kernel=kernelFunction) # 创建一个SVC的实例
    grid = GridSearchCV(clf, param_grid=parameters, cv=numOfFolds) # 创建一个GridSearchCV实例
    grid.fit(trX, trY) # grid寻优c, g
    print("The best parameters are %s with a score of %g" % (grid.best_params_, grid.best_score_))
    clf = svm.SVC(kernel=kernelFunction, C=grid.best_params_['C'], gamma=grid.best_params_['gamma'])
else:
    clf = svm.SVC(kernel=kernelFunction)
    
clf.fit(trX, trY) # 训练模型
predY = clf.predict(teX) # 预测测试集
Acc = sum(predY==teY)/len(teY) # 计算预测正确的样本数
print('Prediction Accuracy of SVC: %g%%(%d/%d)' % (Acc*100, sum(predY==teY), len(teY)))
```

* SVM运行时间较长，将命令写到脚本中再用qsub提交任务 <br>
work_mySVC.sh
```bash
#!/bin/bash
#$ -S /bin/bash
#$ -N mySVC
#$ -j y
#$ -cwd

# 激活base环境，保证计算节点上正常运行python3
source /opt/miniconda3/bin/activate
conda activate

# SVC分类器：在命令行指定训练集、测试集，规格化，线性核，参数寻优，10次交叉
echo '------ scale: 1; kernel: linear; chooseCG: 1; numOfCV: 10 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 linear 1 10
echo

# 规格化，线性核，不参数寻优
echo '------ scale: 1; kernel: linear; chooseCG: 0 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 linear 0
echo

# 规格化，径向基核(rbf)，参数寻优，10次交叉
echo '------ scale: 1; kernel: rbf; chooseCG: 1; numOfCV: 10 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 rbf 1 10
echo

# 规格化，径向基核(rbf)，不参数寻优
echo '------ scale: 1; kernel: rbf; chooseCG: 0 --------'
python3 mySVC.py EI_train.txt EI_test.txt 1 rbf 0
echo
```
```
# qsub提交任务
$ qsub work_mySVC.sh
```

* 尝试更多的选项搭配，看精度变化，比如数据不规格化时，各种核函数、是否寻优、不同交叉验证次数等情形的预测精度。

## 作业
1. 尽量看懂`参考程序`的每一行代码。
2. 熟练使用sklearn包中的NB, SVC分类器。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* NB手册：[sklearn.naive_bayes.GaussianNB](https://scikit-learn.org/stable/modules/naive_bayes.html#gaussian-naive-bayes)
* SVM手册：[sklearn.svm.SVC](https://scikit-learn.org/stable/modules/svm.html#classification)
