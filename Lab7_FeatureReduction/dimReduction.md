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

## 准备SVM参数寻优函数
* 将以下程序保存为：optimise_svm_cv.py
```python3
from sklearn.model_selection import GridSearchCV # 导入参数寻优包
import numpy as np
from sklearn import svm

def optimise_svm_cv(X, y, kernelFunction, numOfFolds):
    C_range = np.power(2, np.arange(-1, 6, 1.0)) # 指定C的范围
    gamma_range = np.power(2, np.arange(0, -8, -1.0)) # 指定g的范围
    epsilon_range = np.power(2, np.arange(-8, -1, 1.0)) # 指定p的范围
    parameters = dict(gamma=gamma_range, C=C_range, epsilon=epsilon_range) # 将c, g, p组成字典，用于参数的grid遍历
    
    reg = svm.SVR(kernel=kernelFunction) # 创建一个SVR的实例
    grid = GridSearchCV(reg, param_grid=parameters, cv=numOfFolds) # 创建一个GridSearchCV实例
    grid.fit(X, y) # grid寻优c, g, p
    print("The best parameters are %s with a score of %g" % (grid.best_params_, grid.best_score_))
    return grid
```

## 准备预测表现作图函数
* 将以下程序保存为：do_plot.py
```python3
import matplotlib.pyplot as plt
import numpy as np

def do_plot(teY, predY, modelName, plotFileName):
    R2 = 1- sum((teY - predY) ** 2) / sum((teY - teY.mean()) ** 2)
    plt.figure()
    plt.scatter(teY, predY,  color='black') # 做测试集的真实Y值vs预测Y值的散点图
    parameter = np.polyfit(teY, predY, 1) # 插入拟合直线
    f = np.poly1d(parameter)
    plt.plot(teY, f(teY), color='blue', linewidth=3)
    plt.xlabel('Observed Y')
    plt.ylabel('Predicted Y')
    plt.title('Prediction performance using %s' % modelName)
    r2text = 'Predicted R2: %g' % R2
    textPosX = min(teY) + 0.2*(max(teY)-min(teY))
    textPosY = max(predY) - 0.2*(max(predY)-min(predY))
    plt.text(textPosX, textPosY, r2text, bbox=dict(edgecolor='red', fill=False, alpha=0.5))
    plt.savefig(plotFileName)
    print('Plot of prediction performance is save into %s' % plotFileName)
```

## 准备PCA函数
* 将以下程序保存为：do_PCA.py
```python3
from sklearn.decomposition import PCA # 导入PCA包
import numpy as np

def do_PCA(trX, teX, percentVar=0.9, solver='full'):
    pca = PCA(n_components = percentVar, svd_solver = solver) # 创建一个PCA实例
    trX = pca.fit_transform(trX)
    teX = pca.transform(teX)
    print("The cummulative percentage of variance explained: %s" % np.cumsum(pca.explained_variance_ratio_))
    print("The estimated number of components: %d" % pca.n_components_)
    return trX, teX
```

## 1. 使用PCA进行特征压缩降维，以保留主成分建立SVR模型
* 参考程序：myPCA_SVR.py
```python3
import numpy as np
from sklearn import svm # 导入svm包
import sys
from sklearn import preprocessing # 导入数据预处理包
from pytictoc import TicToc # 导入程序计时包
from do_PCA import * # 导入自编函数
from optimise_svm_cv import * # 导入自编函数
from do_plot import * # 导入自编函数

if __name__ == '__main__':
    train = np.loadtxt(sys.argv[1], delimiter='\t') # 载入训练集
    test = np.loadtxt(sys.argv[2], delimiter='\t') # 载入测试集
    modelName = 'SVR'

    trX = train[:,1:]
    trY = train[:,0]
    teX = test[:,1:]
    teY = test[:,0]

    isScale = int(sys.argv[3]) # 建模前，是否将每个特征归一化到[-1,1]
    kernelFunction = sys.argv[4] # {‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’}, default=’rbf’
    numOfFolds = int(sys.argv[5]) # 是否寻找最优参数：c, g, p

    percentVar = float(sys.argv[6]) # 主成分累计解释的百分比
    if percentVar < 1 and percentVar > 0:
        trX, teX = do_PCA(trX, teX, percentVar) # 调用自编函数
    else:
        np.warnings.warn('Incorrect percentage of explained variance specified!! The PCA will be omitted!!')

    if isScale:
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(-1,1))
        trX = min_max_scaler.fit_transform(trX)
        teX = min_max_scaler.transform(teX)
    
    t = TicToc() # 创建一个TicToc实例
    t.tic()
    if numOfFolds > 2: # 如果k-fold > 2, 则进行参数寻优
        grid = optimise_svm_cv(trX, trY, kernelFunction, numOfFolds) # 调用自编函数
        print('Time cost in optimising c-g-p: %gs' % t.tocvalue(restart=True))
        bestC = grid.best_params_['C']
        bestGamma = grid.best_params_['gamma']
        bestEpsilon = grid.best_params_['epsilon']
        reg = svm.SVR(kernel=kernelFunction, C=bestC, gamma=bestGamma, epsilon=bestEpsilon)
    else: # 否则不寻优，使用svm默认参数
        reg = svm.SVR(kernel=kernelFunction)
        
    reg.fit(trX, trY) # 训练模型
    print('Time cost in building model: %gs' % t.tocvalue(restart=True))
    predY = reg.predict(teX) # 预测测试集
    print('Time cost in predicting Y of test set: %gs\n' % t.tocvalue(restart=True))

    R2 = 1- sum((teY - predY) ** 2) / sum((teY - teY.mean()) ** 2)
    RMSE = np.sqrt(sum((teY - predY) ** 2)/len(teY))
    print('Predicted R2(coefficient of determination) of %s: %g' % (modelName, R2))
    print('Predicted RMSE(root mean squared error) of %s: %g' % (modelName, RMSE))

    # Plot outputs
    if len(sys.argv) > 7:
        plotFileName = sys.argv[7]
        do_plot(teY, predY, modelName, plotFileName) # 调用自编函数
    
```

## 2. 使用MI-filter进行特征选择，以保留变量建立SVR模型
* 参考程序：myMIFilter_SVR.py
```python3
import numpy as np
from sklearn import svm # 导入svm包
import sys
from sklearn import preprocessing # 导入数据预处理包
from sklearn.model_selection import GridSearchCV # 导入参数寻优包
import matplotlib.pyplot as plt
from pytictoc import TicToc # 导入程序计时包
from sklearn.decomposition import PCA # 导入PCA包

def do_PCA(trX, teX, percentVar=0.9, solver='full'):
    pca = PCA(n_components = percentVar, svd_solver = solver) # 创建一个PCA实例
    trX = pca.fit_transform(trX)
    teX = pca.transform(teX)
    print("The cummulative percentage of variance explained: %s" % np.cumsum(pca.explained_variance_ratio_))
    print("The estimated number of components: %d" % pca.n_components_)
    return trX, teX

def optimise_svm_cv(X, y, kernelFunction, numOfFolds):
    C_range = np.power(2, np.arange(-1, 6, 1.0)) # 指定C的范围
    gamma_range = np.power(2, np.arange(0, -8, -1.0)) # 指定g的范围
    epsilon_range = np.power(2, np.arange(-8, -1, 1.0)) # 指定p的范围
    parameters = dict(gamma=gamma_range, C=C_range, epsilon=epsilon_range) # 将c, g, p组成字典，用于参数的grid遍历
    
    reg = svm.SVR(kernel=kernelFunction) # 创建一个SVR的实例
    grid = GridSearchCV(reg, param_grid=parameters, cv=numOfFolds) # 创建一个GridSearchCV实例
    grid.fit(X, y) # grid寻优c, g, p
    print("The best parameters are %s with a score of %g" % (grid.best_params_, grid.best_score_))
    return grid

def do_plot(teY, predY, modelName, plotFileName):
    plt.figure()
    plt.scatter(teY, predY,  color='black') # 做测试集的真实Y值vs预测Y值的散点图
    parameter = np.polyfit(teY, predY, 1) # 插入拟合直线
    f = np.poly1d(parameter)
    plt.plot(teY, f(teY), color='blue', linewidth=3)
    plt.xlabel('Observed Y')
    plt.ylabel('Predicted Y')
    plt.title('Prediction performance using %s' % modelName)
    r2text = 'Predicted R2: %g' % R2
    textPosX = min(teY) + 0.2*(max(teY)-min(teY))
    textPosY = max(predY) - 0.2*(max(predY)-min(predY))
    plt.text(textPosX, textPosY, r2text, bbox=dict(edgecolor='red', fill=False, alpha=0.5))
    plt.savefig(plotFileName)
    print('Plot of prediction performance is save into %s' % plotFileName)

if __name__ == '__main__':
    train = np.loadtxt(sys.argv[1], delimiter='\t') # 载入训练集
    test = np.loadtxt(sys.argv[2], delimiter='\t') # 载入测试集
    modelName = 'SVR'

    trX = train[:,1:]
    trY = train[:,0]
    teX = test[:,1:]
    teY = test[:,0]

    isScale = int(sys.argv[3]) # 建模前，是否将每个特征归一化到[-1,1]
    kernelFunction = sys.argv[4] # {‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘precomputed’}, default=’rbf’
    numOfFolds = int(sys.argv[5]) # 是否寻找最优参数：c, g, p

    percentVar = float(sys.argv[6]) # 主成分累计解释的百分比
    if percentVar < 1 and percentVar > 0:
        trX, teX = do_PCA(trX, teX, percentVar)
    else:
        np.warnings.warn('Incorrect percentage of explained variance specified!! The PCA will be omitted!!')

    if isScale:
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(-1,1))
        trX = min_max_scaler.fit_transform(trX)
        teX = min_max_scaler.transform(teX)
    
    t = TicToc() # 创建一个TicToc实例
    t.tic()
    if numOfFolds > 2: # 如果k-fold > 2, 则进行参数寻优
        grid = optimise_svm_cv(trX, trY, kernelFunction, numOfFolds)
        print('Time cost in optimising c-g-p: %gs' % t.tocvalue(restart=True))
        bestC = grid.best_params_['C']
        bestGamma = grid.best_params_['gamma']
        bestEpsilon = grid.best_params_['epsilon']
        reg = svm.SVR(kernel=kernelFunction, C=bestC, gamma=bestGamma, epsilon=bestEpsilon)
    else: # 否则不寻优，使用svm默认参数
        reg = svm.SVR(kernel=kernelFunction)
        
    reg.fit(trX, trY) # 训练模型
    print('Time cost in building model: %gs' % t.tocvalue(restart=True))
    predY = reg.predict(teX) # 预测测试集
    print('Time cost in predicting Y of test set: %gs\n' % t.tocvalue(restart=True))

    R2 = 1- sum((teY - predY) ** 2) / sum((teY - teY.mean()) ** 2)
    RMSE = np.sqrt(sum((teY - predY) ** 2)/len(teY))
    print('Predicted R2(coefficient of determination) of %s: %g' % (modelName, R2))
    print('Predicted RMSE(root mean squared error) of %s: %g' % (modelName, RMSE))

    # Plot outputs
    if len(sys.argv) > 7:
        plotFileName = sys.argv[7]
        do_plot(teY, predY, modelName, plotFileName)
    
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
