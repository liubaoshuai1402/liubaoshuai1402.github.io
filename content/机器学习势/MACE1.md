---
title: 如何构建训练集用于训练机器学习势
date: 2025-05-16
series: ["ASE","DeepMD-Kit"]
categories: [机器学习势]
---

# 如何构建训练集用于训练机器学习势

## 前言

本文将介绍如何构建一个训练集，用于训练MACE势以及DP势。需要提前利用AIMD获取DFT数据集。这里的AIMD软件是VASP。

软件：ASE、DeepMD-Kit

**注意**：本文仅供参考，欢迎指出错误或分享补充。无能力提供任何指导，**求教者切勿留言**。

## 构建MACE势的训练集

MACE接受的训练集非常简单，一个`xyz`文件，包含了各种构型和它们对应的DFT数据标签，以及单原子的DFT数据，需要额外的标签`config_type=IsolatedAtom`。

顺带一提，ASE可以输出一种所谓的[Extended XYZ format](https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#extxyz)，会把各种各样的信息（有点类似OVITO中的全局信息），放到`xyz`文件的第二行。这一行会很长很长。MACE所采用的训练集输入格式就是它。

假设我们要构建水分子的MACE势，当前所处的路径下，有两个文件夹。

一个名为`H2O`的文件夹，里面存放着进行第一性原理分子动力学后得到的`OUTCAR.tar.gz`文件。路径为`H2O/OUTCAR.tar.gz`。

一个名为`IsolatedAtoms`的文件夹，里面存放着涉及元素（这里是H、O）的单原子的单点能计算（ISPIN=2）。路径分别为`IsolatedAtoms/H/OUTCAR.tar.gz`和`IsolatedAtoms/O/OUTCAR.tar.gz`。

### 代码展示

```python
from ase.io import read,write
import random

#定义一个简单的函数用于打标签,这里可以自由更改标签的名字
def addlabel(configs,energy_label='energy_dft',forces_label='forces_dft',stress_label='stress_dft',is_isolated=False):
    if is_isolated == False:
        for at in configs:
            at.info[energy_label] = at.get_potential_energy(force_consistent=True)
            at.arrays[forces_label] = at.get_forces()
            at.info[stress_label] = at.get_stress(voigt=True)
    if is_isolated == True:
        for at in configs:
            at.info['config_type'] = 'IsolatedAtom'
            at.info[energy_label] = at.get_potential_energy(force_consistent=True)
            at.arrays[forces_label] = at.get_forces()
            at.info[stress_label] = at.get_stress(voigt=True)

#read()函数，这里，第一个参数是所读文件路径，第二个参数是切片slice
IsolatedH = read('IsolatedAtoms/H/OUTCAR.tar.gz',':')
IsolatedO = read('IsolatedAtoms/O/OUTCAR.tar.gz',':')
IsolatedAtoms = IsolatedH + IsolatedO
addlabel(configs=IsolatedAtoms,is_isolated=True)

#这里的slice的意思是从第一个结构开始到最后一个结构，每100个结构取一个
db = read('H2O/OUTCAR.tar.gz','::100')
addlabel(configs=db)

#打乱训练集，这对训练非常重要
random.seed(42)
random.shuffle(db)

#将打过标签的数据集合并
db = db + IsolatedAtoms
write('trainset.xyz',db)


```

这里有趣的一点是，为什么对于单个结构的`OUTCAR`，也要进行切片：`IsolatedH = read('IsolatedAtoms/H/OUTCAR.tar.gz',':')`，而不是`IsolatedH = read('IsolatedAtoms/H/OUTCAR.tar.gz')`。

这是因为`read()`函数读取只有一个原子的结构是，会返回`atom`类的实例，而非`atoms`类的实例。而`atom`类不支持`info`属性，会很麻烦。

即便是**多原子的单结构**，如果想把他们合并成一个**轨迹**，也一定要用切片的形式`":"`读取，因为两个`atoms`类的实例加和会得到一个新的`atoms`类的实例。

```python
from ase.io import read,write

#假设H2O.cif是一个含有一个水分子的胞
#这样只会输出一个含有两个水分子的胞
db1 = read('H2O.cif')
db2 = read('H2O.cif')
db = db1 + db2
write('H2O.xyz',db)

#使用切片后，read会返回list[atoms],再进行加和得到的是list[atoms1,atoms2],用write()函数写的时候，就能依次形成轨迹了
db1 = read('H2O.cif',':')
db2 = read('H2O.cif',':')
db = db1 + db2
write('H2O.xyz',db)
```

所以，如果我们想让单一结构也加入到我们的数据集时，也要记得用切片的形式进行读取。不过要铭记于心的是，用**切片的形式读取**的其实是一个`list[atoms]`，要使用其中`atoms`实例的方法时，记得用`for`循环遍历其中的`atoms`实例。



