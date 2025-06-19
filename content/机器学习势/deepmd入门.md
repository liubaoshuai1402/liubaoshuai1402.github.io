---
title: DeepMD入门
date: 2025-06-18
series: ["deepmd"]
categories: [机器学习势]
---

# DeepMD入门

## 前言

之前一直用的mace，现在来了解一下deepmd吧。

### 1.dp的数据集结构

与mace相对简单的数据集结构相比（只需提供打标签后xyz文件），dp组织它的数据集以一种相对不自由的形式，但也更加清楚明了。

dp的数据集是由dpdata处理DFT结果得到的。一个标准的dp数据集路径下，包含了三种文件，分别是：`type.raw`、`type_map.raw`和`set.000`。

`type.raw`其实就是POSCAR中的原子，只不过把元素符号换成了数字。`type_map.raw`则是保留了数字到元素符号的映射关系。

~~很困惑为什么要这样设计~~

set.000文件夹则是存放DFT标签的，具体如坐标、能量、力等。

### 2.如何生成数据集

任务：想象一下，你有多个OUTCAR，分别是：

1000K下：C、O体系；C、O、H体系

2000K下：C、O体系；C、O、H体系

如何把他们整合成训练集和数据集呢？

每个体系AIMD的OUTCAR包含10000个轨迹。但因为相邻轨迹的结构比较接近，我们决定每100步取一个结构作为数据集。

现在你有以下路径，`1000K/CO/OUTCAR`、`1000K/COH/OUTCAR`、`2000K/CO/OUTCAR`、`2000K/COH/OUTCAR`

此外呢，如果把这四体系合并后再划分训练集、验证集，随机性会更大，比如说某个体系的训练集取得很多，而某个体系则几乎全被当作验证集。

为了避免这样的事情发生，我们可以先把每个体系划分成小的训练集、验证区，然后把这些小的训练集、验证集合并。

这就是总体思路，来看看怎么实现吧。

```python
import dpdata
import numpy as np
#生成一个列表用于dpdata取子集的索引，dpdata的system并不支持直接的切片规则，很无语
#从第1个结构到1w个，每100步取一个，很灵活，可以自由改变起始，结束和步长
IndicesInitial = list(range(0,10000,100))
#从100个结构的索引中随机取20个作为验证集，20%
np.random.seed(42)
ValidIndices = list(np.random.choice(IndicesInitial,20,replace=None))
#从数据集中剔除验证集，得到训练集，这里为了实现列表元素进行集合运算，先把他们变成集合
TrainIndices = list(set(IndicesInitial)-set(ValidIndices))
#只读了一个OUTCAR作为例子，适当修改
dpSystem = dpdata.LabeledSystem('1000K/CO/OUTCAR')

dpSystemTrain = dpSystem.sub_system(TrainIndices)
dpSystemTrain.to("deepmd/npy", "trainset/1000K/CO", set_size=dpSystemTrain.get_nframes())

dpSystemValid = dpSystem.sub_system(ValidIndices)
dpSystemValid.to("deepmd/npy", "validset/1000K/CO", set_size=dpSystemValid.get_nframes())


```

`set_size=dpSystemTrain.get_nframes()`,关于set的尺寸为什么要正好等于所有结构，这个问题追溯起来非常远古。

曾经dp的开发者应该希望一个system下有多个set.00x，最后一个set.00x作为测试集。但后来他们放弃了。

所以现在的用法是，只要一个set.000即可，至于这个system下的数据用于训练、验证，在input.json中说明即可。 

顺带一提，dp训练模型时并不要求测试集，mace的话则可以提供。

下图是input.json的一部分：

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/dp1.png" style="zoom:67%;" />

把相应路径替换成`"trainset/1000K/CO"`等等即可

dp中system的概念，同一个system中，结构涉及的元素和原子个数应该完全相同。

### 3.利用dpdata生成扰动结构

dpdata本身已经提供了很好的[示例](https://docs.deepmodeling.com/projects/dpdata/en/master/systems/system.html), 这里简单搬运一下，并添加了一个更深入的应用场景。

```python
import dpdata
perturbed_system = dpdata.System("CONTCAR").perturb(
    pert_num=3,
    cell_pert_fraction=0.05,
    atom_pert_distance=0.6,
    atom_pert_style="normal",
)
print(perturbed_system.data)
```

利用class `system`的`perturb`方法就可以轻松实现。

扰动晶胞大小可以模拟受力变形的情况。随机扰动是种很好的方法。

但根据计算弹性系数的变形模式施加应变，是不是会让数据集更有意义，同时让训练出的势在计算弹性系数时更具有优势？

那么，如何根据晶型确定需要施加的变形模式呢？这需要一定的物理知识，不过好在已经有人帮我们把做成了软件，名叫[elastool](https://github.com/zhongliliu/elastool)。

这个软件是为Linux设计的，如果装到Windows下要改一下源码才能正确运行。

不过我们可以在它的`strain_matrix.py`文件中，找到我们所需的变形模式。~~（当然，可以问问大G老师）~~

最后，这些微扰后生成的结构，是拿去做单点能计算得到能量、力（==个人认为不需要结构优化==，因为本质是提供DFT数据、而不是DFT理论下合理的结构）。

如果生成100个结构，怎么把他们整合到一起呢，总不能在dp的input.json中提供一百个训练集的路径吧！

其实很简单，用一个循环，把他们的`LabeledSystem`对象加起来就行了，是的，class `system`支持直接进行加和。就像ase中的atoms一样。

写一个简单的测试：

```python
import dpdata

#把拓展名写成OUTCAR是为了让dpdata猜它是OUTCAR，不用指定格式
structure1 = dpdata.LabeledSystem('1.OUTCAR')
print(structure1.data['energies'])
structure2 = dpdata.LabeledSystem('2.OUTCAR')
print(structure2.data['energies']) 
StructureSet = structure1 + structure2
print(StructureSet.data['energies'])
```

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/dp2.png"/>

可以看到，很成功。另外，因为dp对于system的规定是相同原子类型且数量相等，所以如果不满足这个条件的数据，不能加在一起，会报错，也是理所当然的。