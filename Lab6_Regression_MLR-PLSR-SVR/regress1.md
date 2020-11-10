# 实验六：回归模型之MLR, PLSR, SVR

## 实验目的
* 1）根据[实验五](https://github.com/ZhijunBioinf/Pattern-Recognition-and-Prediction/blob/master/Lab5_PeptideSequencesCoding/sequence_coding2.md)对ACE抑制剂多肽序列表征的AA531特征，构建训练集与测试集
* 2）使用多元线性回归(Multiple Linear Regression, MLR)完成ACE抑制剂活性预测
* 3）使用偏最小二乘回归(Partial Least Squares Regression, PLSR)完成ACE抑制剂活性预测
* 4）使用支持向量回归(Support Vector Regression, SVR)完成ACE抑制剂活性预测

## 准备工作目录
```
$ mkdir lab_06
$ cd lab_06
# 对lab_05路径中的ACE抑制剂多肽活性及序列的AA531特征(result.txt)建立软链接，并重命名为'ACEtriPeptides_YandAA531.txt'
$ ln -s ../lab_05/result.txt ACEtriPeptides_YandAA531.txt
```

## 1. 训练集与测试集构建
* 参考程序：getTrainTest_regression.py
```python3
import numpy as np
import sys
from sklearn.model_selection import train_test_split # 用于产生训练集、测试集

DataFileName = sys.argv[1]
Data = np.loadtxt(DataFileName, delimiter = '\t') # 载入数据
X = Data[:,1:]
Y = Data[:,0]

testSize = float(sys.argv[2]) # 命令行指定test set比例
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size = testSize)

trainingSetFileName = sys.argv[3]
testSetFileName = sys.argv[4]
np.savetxt(trainingSetFileName, np.hstack((y_train.reshape(-1,1), X_train)), fmt='%g', delimiter='\t') # 将Y与X以列组合后，保存到文件
np.savetxt(testSetFileName, np.hstack((y_test.reshape(-1,1), X_test)), fmt='%g', delimiter='\t')
print('Generate training set(%d%%) and test set(%d%%): Done!' % ((1-testSize)*100, testSize*100))
```

```bash
# 构建训练集与测试集：在命令行指定序列表征数据、测试集比例、train文件名、test文件名
$ python3 getTrainTest_regression.py ACEtriPeptides_YandAA531.txt 0.2 ACE_train.txt ACE_test.txt
```

## 2. 以MLR完成ACE抑制剂活性预测
* 参考程序：myMLR.py
```python3
import numpy as np
from sklearn import linear_model # 导入MLR包
import sys
import matplotlib.pyplot as plt

train = np.loadtxt(sys.argv[1], delimiter='\t') # 载入训练集
test = np.loadtxt(sys.argv[2], delimiter='\t') # 载入测试集
isNormalizeX = bool(int(sys.argv[3])) # 是否标准化每个x
modelName = 'MLR'

reg = linear_model.LinearRegression(normalize = isNormalizeX) # 创建一个MLR的实例
trX = train[:,1:]
trY = train[:,0]
reg.fit(trX, trY) # 训练模型

teX = test[:,1:]
teY = test[:,0]
predY = reg.predict(teX) # 预测测试集

R2 = 1- sum((teY - predY) ** 2) / sum((teY - teY.mean()) ** 2)
RMSE = np.sqrt(sum((teY - predY) ** 2)/len(teY))
print('Predicted R2(coefficient of determination) of %s: %g' % (modelName, R2))
print('Predicted RMSE(root mean squared error) of %s: %g' % (modelName, RMSE))

# Plot outputs
plotFileName = sys.argv[4]
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
```

```bash
# MLR模型：在命令行指定训练集、测试集、是否对特征标准化、图名
$ python3 myMLR.py ACE_train.txt ACE_test.txt 0 ObsdYvsPredY_MLR.pdf
```

## 3. 以PLSR完成ACE抑制剂活性预测
* 参考程序：myPLSR.py
```python3
import numpy as np
from sklearn.cross_decomposition import PLSRegression # 导入PLSR包
import sys
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_predict # 导入交叉验证包

def optimise_pls_cv(X, y, nCompMax): # 以交叉验证技术获得不同潜变量个数情形下的MSE
    MSEVec = []
    nCompVec = np.arange(1, nCompMax)
    for n_comp in nCompVec:
        pls = PLSRegression(n_components=n_comp) # 创建一个PLSR的实例
        y_cv = cross_val_predict(pls, X, y, cv=10) # 获得10次交叉的预测Y
        mse = sum((y - y_cv.ravel()) ** 2)/len(y) # 计算MSE, NOTE: y_cv维度为2，需转成向量
        MSEVec.append(mse)
    bestNComp = np.argmin(MSEVec) # 获得最小MSE对应的下标，即最优潜变量个数
    
    with plt.style.context('ggplot'): # 以潜变量个数为x轴，对应的MSE为y轴，作图
        plt.plot(nCompVec, np.array(MSEVec), '-v', color='blue', mfc='blue') # 带标记点的折线图
        plt.plot(nCompVec[bestNComp], np.array(MSEVec)[bestNComp], 'P', ms=10, mfc='red') # 在图上标记出最小MSE对应的点
        plt.xlabel('Number of PLS components')
        plt.xticks = nCompVec
        plt.ylabel('MSE')
        plt.title('Optimise the number of PLS components')
        plt.savefig('optimizePLSComponents.pdf')
        
    return bestNComp

if __name__ == '__main__':
    train = np.loadtxt(sys.argv[1], delimiter='\t') # 载入训练集
    test = np.loadtxt(sys.argv[2], delimiter='\t') # 载入测试集
    nCompMax = int(sys.argv[3]) # 潜变量个数上限
    modelName = 'PLSR'

    trX = train[:,1:]
    trY = train[:,0]
    bestNComp = optimise_pls_cv(trX, trY, nCompMax)+1 # 得到最优潜变量个数（特征降维思想）
    print('The best number of PLS components: %d' % bestNComp)
    reg = PLSRegression(n_components = bestNComp) # 创建一个PLSR的实例
    reg.fit(trX, trY) # 训练模型

    teX = test[:,1:]
    teY = test[:,0]
    predY = reg.predict(teX) # 预测测试集
    predY = predY.ravel() # NOTE: predY维度为2，需转成向量

    R2 = 1- sum((teY - predY) ** 2) / sum((teY - teY.mean()) ** 2)
    RMSE = np.sqrt(sum((teY - predY) ** 2)/len(teY))
    print('Predicted R2(coefficient of determination) of %s: %g' % (modelName, R2))
    print('Predicted RMSE(root mean squared error) of %s: %g' % (modelName, RMSE))

    # Plot outputs
    plotFileName = sys.argv[4]
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
```

```bash
# PLSR模型：在命令行指定训练集、测试集、最大潜变量个数、图名
$ python3 myPLSR.py ACE_train.txt ACE_test.txt 20 ObsdYvsPredY_PLSR.pdf
```

## 4. 以Decision Tree进行剪接位点识别
* 参考程序：myDT.py
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
$ pip3 install --user graphviz -i https://pypi.tuna.tsinghua.edu.cn/simple
# DT分类器：在命令行指定训练集、测试集、DT图文件名
$ python3 myDT.py EI_train.txt EI_test.txt EISplicing_DecisionTreeGraph
```
[获得的DT树](https://github.com/ZhijunBioinf/Pattern-Recognition-and-Prediction/blob/master/Lab3_Classifiers_KNN-LR-DT/EISplicing_DecisionTreeGraph.pdf)

## 作业
1. 尽量看懂`参考程序`的每一行代码。
2. 参考程序kSpaceCoding_general.py中，供体位点序列的第71、72位保守二核苷酸GT是在程序中指定的，试着改写程序，实现从命令行传递`位置信息`给程序。
3. 参考程序getTrainTest.py中，测试集的比例testSize是在程序中指定的，试着改写程序，实现从命令行传递`划分比例`给程序。
4. 熟练使用sklearn包中的不同分类器。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* MLR手册：[sklearn.linear_model.LinearRegression](https://scikit-learn.org/stable/modules/linear_model.html#ordinary-least-squares)
* PLSR手册：[sklearn.cross_decomposition.PLSRegression](https://scikit-learn.org/stable/modules/cross_decomposition.html)
* DT手册：[sklearn.tree.DecisionTreeClassifier](https://scikit-learn.org/stable/modules/tree.html#classification)
