---
title: 径向分布函数（RDF）-OVITO
date: 2025-05-12
series: ["OVITO"]
categories: [MolecularDynamics]
---

# 径向分布函数（RDF）计算 by OVITO

## 前言

[OVITO Python Reference — OVITO Python Reference 3.12.3 documentation](https://docs.ovito.org/python/index.html) 是一个开源且功能强大的分子动力学后处理软件包。

本文将介绍如何利用 OVITO python module 计算单个结构以及一段轨迹（多个结构）内的径向分布函数。

适用于无机非晶体，其他体系慎用。

软件：OVITO、matplotlib、numpy

**注意**：本文仅供参考，欢迎指出错误或分享补充。无能力提供任何指导，**求教者切勿留言**。

## The partial RDFs of a single crystal structure 

### 代码展示

```python
#这段代码用于计算 RDF by OVITO
from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier
import numpy as np

#导入一个氧化锆（ZrO2）的cif文件，所有OVITO支持的输入文件格式都可以（确保这个.py文件的路径下有这样一个cif文件，也可以稍微修改指定结构路径）
pipeline = import_file("ZrO.cif")

#施加一个名叫 CoordinationAnalysisModifier 的修饰器，cutoff用于控制截断半径，number_of_bins用于控制网格细分度（大小100-1000内都可以试试）
pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 5.0, number_of_bins = 500,partial=True))

#进行计算
rdf_table = pipeline.compute().tables['coordination-rdf']

#得到用于画图的横纵坐标，默认第一列是x轴数据，其余列是y轴数据
total_rdf = rdf_table.xy()

#记录total_rdf中y轴数据对应是什么类型的pair-wise
#这个例子中，输出为：
#g(r) for pair-wise type combination O-O:
#g(r) for pair-wise type combination O-Zr:
#g(r) for pair-wise type combination Zr-Zr:
#说明total_rdf是一个四列的数据，第一列是x轴坐标（其实是bin），第二列就是不同pair-wise的RDF数据，依次为 O-O,O-Zr,Zr-Zr
rdf_names = rdf_table.y.component_names
for component, name in enumerate(rdf_names):
    print("g(r) for pair-wise type combination %s:" % name)
    
#将total_rdf保存为txt文件，用于后续画图
np.savetxt("total_rdf.txt", total_rdf)

```

```python
#这段代码用于绘图
import numpy as np
import matplotlib.pyplot as plt

rdf_table = np.loadtxt('total_rdf.txt')

#g(r) for pair-wise type combination O-O:
#g(r) for pair-wise type combination O-Zr:
#g(r) for pair-wise type combination Zr-Zr:

#这里取的是 total_rdf.txt 中的第一列（对应[:,0]）和第三列（对应[:,2]），所以绘制的是 Zr-O pair-wise的partial RDF
plt.plot(rdf_table[:,0], rdf_table[:,2])

#matplotlib的常规设置，问问万能的小迪老师吧
title_font = {'fontsize': 24, 'fontfamily': 'Times New Roman'}
xlabel_font = {'fontsize': 22, 'fontfamily': 'Times New Roman'}
ylabel_font = {'fontsize': 22, 'fontfamily': 'Times New Roman'}

plt.title("RDF Zr-O", fontdict=title_font,pad=8)
plt.xlabel(xlabel='distance r',fontdict=xlabel_font,loc='center',labelpad=8)
plt.ylabel(ylabel='g(r)',fontdict=ylabel_font,loc='center',labelpad=8)
plt.tick_params(axis='both', which='major', labelsize=16, direction='in')

ax = plt.subplot()

#因为只有一个静态结构，pair-wise的某些峰很高，所以这里的y轴坐标上限设置大一些，为200，可灵活改变
#x轴设置为6，稍大于截断半径cutoff即可，因为本身也只在截断半径以内统计
ax.set_ylim(0,200)
plt.xlim(0,6)

fig = plt.gcf()

fig.set_size_inches(1200/100, 800/100)
plt.savefig('output.png', dpi=100)

plt.show()
```

### 结果展示

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/output.png"/>

## The partial RDFs of trajectories

### 代码展示

```python
#这段代码用于计算一定时间内（一段轨迹）的平均 RDFs
from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier,TimeAveragingModifier
import numpy as np

#读入轨迹文件，这里是利用 VASP 进行 AIMD 后得到的 XDATCAR 文件
pipeline = import_file("XDATCAR")

#打印轨迹中的结构数
print("Number of MD frames:", pipeline.num_frames)

#添加修饰器，与单个晶体结构相比，多了 TimeAveragingModifier 修饰器
pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 5.0, number_of_bins = 500,partial=True))
pipeline.modifiers.append(TimeAveragingModifier(operate_on='table:coordination-rdf'))

#计算 RDFs 数据
total_rdf = pipeline.compute().tables['coordination-rdf[average]'].xy()

#记录pair-wise类型
rdf_names = pipeline.compute().tables['coordination-rdf[average]'].y.component_names
for name in rdf_names:
    print("g(r) for pair-wise type combination %s:" % name)

#输出数据，用于后续绘图，不再重复
np.savetxt('rdf.txt', total_rdf, delimiter='\t')
```

## OVITO小知识

受限于Python基础和时间精力的限制，以下内容皆为我个人的有限理解，未能严格考究，仅供参考。

### 在 The partial RDFs of a single crystal structure 的计算代码中

```python
pipeline = import_file("ZrO.cif")
```

这段代码将外来的 cif 文件转化为 ovito 的一个 **Pipeline** 类的实例。

然后我们添加了修饰器，代码省略。

再然后是调用 **Pipeline** 类的 **compute()** 方法，得到一个非常重要的 **DataCollection** 类的实例。

再访问这个 **DataCollection** 实例的 **tables** 属性就得到一个**DataTable** 类的实例（在**DataCollection** 类的定义中，用 **@property** 的 装饰器语法将 **tables()** 方法转化为了 **tables** 属性，所以与 **compute()** 方法相比，不需要 **()** 了）

**DataTable** 类用于储存绘制直方和2d图的数据，当我们添加 **CoordinationAnalysisModifier** 修饰器后，经过 **compute()** 方法后，相应的 partial RDFs 数据就会被储存在这里，可以通过名为 `'coordination-rdf'` 的 键(key) 来检索。

我觉得 **DataTable** 类的实例绝非一个简单的字典，但把它当作字典来理解会比较容易。

```python
rdf_table = pipeline.compute().tables['coordination-rdf']
```

即，**rdf_table** 仍然属于 **DataTable** 类。接着用 `.y.component_names` 来获取 多组y轴数据 的表头。这里太复杂了，~~我也没看太懂~~，不展开了。

```python
rdf_names = rdf_table.y.component_names
```

### 在 The partial RDFs of trajectories 的计算代码中

原本 **compute()** 方法只能用于单个结构，对于 trajectories 需要指定是哪个结构。

但添加 **TimeAveragingModifier** 修饰器后，允许我们计算时间相关量的平均值。详见：https://docs.ovito.org/python/modules/ovito_modifiers.html#ovito.modifiers.TimeAveragingModifier



## References

[ovito.modifiers — OVITO Python Reference 3.12.3 documentation](https://docs.ovito.org/python/modules/ovito_modifiers.html#ovito.modifiers.CoordinationAnalysisModifier)
