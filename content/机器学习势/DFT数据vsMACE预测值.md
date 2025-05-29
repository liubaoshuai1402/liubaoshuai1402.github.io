---
title: 绘制DFT数据和MACE预测值的关系图
date: 2025-05-29
series: ["ASE","MACE"]
categories: [机器学习势]
---

# 绘制DFT数据和MACE预测值的关系图

## 前言

本文主要搬运一下[MACE教程其一](https://colab.research.google.com/drive/1ZrTuTvavXiCxTFyjBV4GqlARxgFwYAtX#scrollTo=v75wGSfOw9dA)，记录一下如何绘制DFT数据和MACE预测值的关系图，用于评估模型的准确性，加上一些自己的理解。

软件：ASE、[aseMolec](https://github.com/imagdau/aseMolec)、matplotlib、numpy

**注意**：本文仅供参考，欢迎指出错误或分享补充。无能力提供任何指导，**求教者切勿留言**。

## 准备评估函数

MACE官方写好了评估的命令行脚本，我们这里给它一个虚拟参数，方便以编程的方式调用它。

```python
from mace.cli.eval_configs import main as mace_eval_configs_main
import sys

def eval_mace(configs, model, output):
    sys.argv = ["program", "--configs", configs, "--model", model, "--output", output]
    mace_eval_configs_main()
```

这里的`eval_mace`函数接受三个参数，数据集如训练集、测试集的路径、训练好的模型的路径以及输出的文件命。

## 评估数据集

```python
#evaluate the training set
eval_mace(configs="data/solvent_xtb_train_200.xyz",
          model="MACE_models/mace01_run-123_stagetwo.model",
          output="tests/mace01/solvent_train.xyz")

#evaluate the test set
eval_mace(configs="data/solvent_xtb_test.xyz",
          model="MACE_models/mace01_run-123_stagetwo.model",
          output="tests/mace01/solvent_test.xyz")
```

这样，MACE就会对数据集进行评估，因为其实DFT数据是已有的，主要是在输出文件中，补上MACE的预测值。

## 画图

然后就可以用aseMolec的如下代码进行画图了

```python
from aseMolec import pltProps as pp
from ase.io import read
import matplotlib.pyplot as plt
from aseMolec import extAtoms as ea
import numpy as np

def plot_RMSEs(db, labs):
    ea.rename_prop_tag(db, 'MACE_energy', 'energy_mace') #Backward compatibility
    ea.rename_prop_tag(db, 'MACE_forces', 'forces_mace') #Backward compatibility

    plt.figure(figsize=(9,6), dpi=100)
    plt.subplot(1,3,1)
    pp.plot_prop(ea.get_prop(db, 'bind', '_xtb', True).flatten(), \
                 ea.get_prop(db, 'bind', '_mace', True).flatten(), \
                 title=r'Energy $(\rm eV/atom)$ ', labs=labs, rel=False)
    plt.subplot(1,3,2)
    pp.plot_prop(ea.get_prop(db, 'info', 'energy_xtb', True).flatten(), \
                 ea.get_prop(db, 'info', 'energy_mace', True).flatten(), \
                 title=r'Energy $(\rm eV/atom)$ ', labs=labs, rel=False)
    plt.subplot(1,3,3)
    pp.plot_prop(np.concatenate(ea.get_prop(db, 'arrays', 'forces_xtb')).flatten(), \
                 np.concatenate(ea.get_prop(db, 'arrays', 'forces_mace')).flatten(), \
                 title=r'Forces $\rm (eV/\AA)$ ', labs=labs, rel=False)
    plt.tight_layout()
    return

train_data = read('tests/mace01/solvent_train.xyz', ':')
test_data = train_data[:3]+read('tests/mace01/solvent_test.xyz', ':') #append the E0s for computing atomization energy errors

plot_RMSEs(train_data, labs=['XTB', 'MACE'])
plot_RMSEs(test_data, labs=['XTB', 'MACE'])
```

`plot_RMSEs`函数中首先把标签重命名了一下，因为MACE的源代码`mace_eval_configs_main()`部分默认打的标签是`MACE_`+什么什么的。

## 涉及到的比较高级的Python语法

画DFT forces VS MACE forces的代码，有一段如下

```python
def get_prop(db, type, prop='', peratom=False, E0={}):
    if peratom:
        N = lambda a : a.get_global_number_of_atoms()
    else:
        N = lambda a : 1
    if type == 'info':
        return np.array(list(map(lambda a : a.info[prop]/N(a), db)))
    if type == 'arrays':
        return np.array(list(map(lambda a : a.arrays[prop]/N(a), db)), dtype=object)
    if type == 'cell':
        return np.array(list(map(lambda a : a.cell/N(a), db)))
    if type == 'meth':
        return np.array(list(map(lambda a : getattr(a, prop)()/N(a), db)))
    if type == 'atom':
        if not E0:
            E0 = get_E0(db, prop)
        return np.array(list(map(lambda a : (np.sum([E0[s] for s in a.get_chemical_symbols()]))/N(a), db)))
    if type == 'bind':
        if not E0:
            E0 = get_E0(db, prop)
        return np.array(list(map(lambda a : (a.info['energy'+prop]-np.sum([E0[s] for s in a.get_chemical_symbols()]))/N(a), db)))
```

这里，forces在ASE的`atoms`类中属于`arrays`型信息。

他是把`db`，也就是数据集（`list[atoms]`），中的每一个构型(`atoms`)，传递给了**匿名函数lambda**，这个匿名函数会提取相应标签的`arrays`型信息。

显示指定`dtype=object`，避免了不同构型中原子数不同带来的生成数组时的报错。

之后再沿`axis=0`方向拼接，然后展开，就变成了一维数据，方便计算处理。
