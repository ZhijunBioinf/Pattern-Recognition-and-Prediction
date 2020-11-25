# 实验八：无监督学习之聚类分析(K-Means, Hierarchical clustering)

## 实验目的
* 1）使用K-means完成聚类分析。
* 2）使用Hierarchical clustering完成聚类分析。
* 3）理解每种聚类方法在不同参数下的聚类表现。

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
> * Centroid models: [K-means](https://en.wikipedia.org/wiki/K-means_clustering), k-means++, [Mean Shift](https://scikit-learn.org/stable/modules/clustering.html#mean-shift), etc.
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
from sklearn.datasets import load_digits # size:1797x64, 10类，每个样本为1个手写数字的8x8维度的图
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

"""
Cluster quality metrics evaluated (see :ref:`clustering_evaluation` for
definitions and discussions of the metrics):
Shorthand    full name
=========== ========================================================
homo         homogeneity score
compl        completeness score
ARI          adjusted Rand index
=========== ========================================================
"""

def bench_k_means(estimator, name, data, labels):
    t0 = time()
    estimator.fit(data)
    t1 = time() - t0
    inertia_ =  estimator.inertia_ # 每个样本到最近聚类中心（质心）的距离的平方和
    homo = metrics.homogeneity_score(labels, estimator.labels_) # 由真实labels和估计的labels计算出homo, 下同
    compl = metrics.completeness_score(labels, estimator.labels_)
    ARI = metrics.adjusted_rand_score(labels, estimator.labels_)
    print('%-9s\t%.2fs\t%i\t%.3f\t%.3f\t%.3f' % (name, t1, inertia_, homo, compl, ARI))
    
def do_plot(data, n_digits, plotFileName): # 基于PCA的结果作图
    reduced_data = PCA(n_components=2).fit_transform(data) # 保留2个主成分
    kmeans = KMeans(init='k-means++', n_clusters=n_digits, n_init=10) # 创建一个KMeans实例
    kmeans.fit(reduced_data)

    x_min, x_max = reduced_data[:, 0].min() - 1, reduced_data[:, 0].max() + 1
    y_min, y_max = reduced_data[:, 1].min() - 1, reduced_data[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.02), np.arange(y_min, y_max, 0.02)) # 网格图
    Z = kmeans.predict(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape) # 获得网格中每个样本点的label.

    plt.figure(1)
    plt.clf()
    plt.imshow(Z, interpolation='nearest', extent=(xx.min(), xx.max(), yy.min(), yy.max()),
               cmap=plt.cm.Paired, aspect='auto', origin='lower')
    plt.plot(reduced_data[:, 0], reduced_data[:, 1], 'k.', markersize=2)
    centroids = kmeans.cluster_centers_     # 得到质心，作图，白色叉标记
    plt.scatter(centroids[:, 0], centroids[:, 1], marker='x', s=169, linewidths=3, color='w', zorder=10)
    plt.title('K-means clustering on the digits dataset (PCA-reduced data)\n'
              'Centroids are marked with white cross')
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xticks(())
    plt.yticks(())
    plt.savefig(plotFileName)
    print('Plot of K-means clustering performance is save into %s' % plotFileName)

if __name__ == '__main__':
    np.random.seed(42)
    X_digits, y_digits = load_digits(return_X_y=True)
    data = scale(X_digits) # 数据标准化处理
    n_samples, n_features = data.shape
    n_digits = len(np.unique(y_digits))
    labels = y_digits
    
    print("n_digits: %d, \t n_samples: %d, \t n_features: %d" % (n_digits, n_samples, n_features))
    print(82 * '_') # 打印横线
    print('init\t\ttime\tinertia\thomo\tcompl\tARI')
    estimator = KMeans(init='k-means++', n_clusters=n_digits, n_init=10) # 创建一个k-means++的实例
    bench_k_means(estimator, name="k-means++", data=data, labels=labels) # 聚类拟合数据，打印信息
    
    estimator = KMeans(init='random', n_clusters=n_digits, n_init=10) # 创建一个经典k-means的实例
    bench_k_means(estimator, name="random", data=data, labels=labels) # 聚类拟合数据，打印信息
    
    pca = PCA(n_components=n_digits).fit(data) # 主成分分析
    estimator = KMeans(init=pca.components_, n_clusters=n_digits, n_init=1) # 传递主成分，质心确定，运行次数n_init设为1
    bench_k_means(estimator, name="PCA-based", data=data, labels=labels) # 聚类拟合数据，打印信息
    print(82 * '_') # 打印横线
    
    do_plot(data, n_digits, 'K-means_clustering.pdf')
    
```

```bash
$ python3 myKMeansDigits.py
```

## 2. 使用Hierarchical clustering完成`手写数字`聚类分析
* 参考程序：myHClusteringDigits.py
```python3
from time import time
import numpy as np
from scipy import ndimage # 导入图形处理包
from matplotlib import pyplot as plt
from sklearn import manifold, datasets # 导入数据降维包manifold，数据集包datasets
from sklearn.cluster import AgglomerativeClustering # 导入HClustering包
from sklearn import metrics

def plot_clustering(X_red, y, labels, title, plotFileName): # 聚类结果可视化
    x_min, x_max = np.min(X_red, axis=0), np.max(X_red, axis=0)
    X_red = (X_red - x_min) / (x_max - x_min)

    plt.figure(figsize=(6, 4))
    for i in range(X_red.shape[0]):
        plt.text(X_red[i, 0], X_red[i, 1], str(y[i]),
                 color=plt.cm.nipy_spectral(labels[i] / 10.),
                 fontdict={'weight': 'bold', 'size': 9})

    plt.xticks([])
    plt.yticks([])
    plt.title(title, size=17)
    plt.axis('off')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(plotFileName)
    print("Plot of hierarchical clustering performance using '%s' is save into '%s'\n" % (title, plotFileName))

if __name__ == '__main__':
    X, y = datasets.load_digits(return_X_y=True) # 导入手写数字数据集
    n_samples, n_features = X.shape
    np.random.seed(0)

    shift = lambda x: ndimage.shift(x.reshape((8, 8)), 0.3*np.random.normal(size=2)).ravel() # Shift an array
    X = np.concatenate([X, np.apply_along_axis(shift, 1, X)]) # 2倍扩充X，便于更好显示聚类结果
    y = np.concatenate([y, y], axis=0) # 2倍扩充y

    print("Spectral embedding for non-linear dimensionality reduction")
    X_red = manifold.SpectralEmbedding(n_components=2).fit_transform(X) # 非线性降维，使用2个主成分
    print("Done!\n")

    for linkage in ('ward', 'average', 'complete', 'single'): # 使用HC聚类算法的不同连接指标
        estimator = AgglomerativeClustering(linkage=linkage, n_clusters=10) # 创建一个HClustering的实例
        t0 = time()
        estimator.fit(X_red)
        homo = metrics.homogeneity_score(y, estimator.labels_) # 由真实labels和估计的labels计算出homo, 下同
        compl = metrics.completeness_score(y, estimator.labels_)
        ARI = metrics.adjusted_rand_score(y, estimator.labels_)
        print("## Use '%s' linkage criterion, time = %.2fs. homoScore: %g, complScore: %g, ARI: %g" % 
            (linkage, time()-t0, homo, compl, ARI))

        title = "%s linkage" % linkage
        plotFileName = "hierarchicalClustering_%sLinkage.pdf" % linkage
        plot_clustering(X_red, y, estimator.labels_, title, plotFileName)
    
```

```bash
$ python3 myHClusteringDigits.py
```


## 作业
1. 尽量看懂`参考程序`的每一行代码。 <br>
2. 熟练使用K-means, Hierarchical clustering完成聚类分析。 <br>
不怕报错，犯错越多，进步越快！

## 参考
* K-means手册：[sklearn.cluster.KMeans](https://scikit-learn.org/stable/modules/clustering.html#k-means)
* Hierarchical clustering手册：[sklearn.cluster.AgglomerativeClustering](https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering)

