---
title: 机器学习势nep
date: 2025-11-13
series: ["nep"]
categories: [机器学习势]
---

# 机器学习势nep

## 前言

nep的优势在于计算资源上的节省，可以模拟十几万个原子体系的分子动力学。这里放一些处理脚本。

## 处理脚本

### 处理loss.out

这个脚本是根据官方提供的matlib脚本丢给gpt改的。

```python
import numpy as np
import matplotlib.pyplot as plt

# 读取数据
loss = np.loadtxt("loss.out")


plt.figure()
plt.loglog(loss[:, 1:6], linewidth=2)
plt.xlabel('Generation/100', fontsize=15)
plt.ylabel('Loss functions', fontsize=15)
plt.tick_params(labelsize=15, length=6)
plt.legend(['Total','L1-Reg','L2-Reg','Energy-train','Force-train'])
plt.tight_layout()

plt.show()
```

### 主动学习流程

这里，在樊老师建议下，我打算用主动学习的完整流程，做一个纯氢气的机器学习势。算是练手，同时也为了采样。

#### 构建初始训练集

这里我没有用AIMD采样，而是采用mace-off的通用势对氢气进行了采样。

在20 * 20 * 15 埃的胞里使用 packmol 随机添加了50个氢分子。

先是10-300K升温，再300K保温0.5ps。保温区采样了100个结构。

然后80个结构正常保留，后20个结构进行扰动生成了40个结构。

一共120个结构，用vasp进行单点计算以后，作为我的初始训练集。

后续可以发现，我的这120的结构的H-H键局限于0.74左右，没有0.6和0.8或者更偏远的数据。

#### 一个turn里要做些什么

##### 训练nep

文件准备：初始训练集，测试集，nep.in，nep.slurm（超算提交任务的脚本）

训练一个nep势

```
version      4
type         1 H
cutoff       6 5
n_max        6 6
basis_size   8 8 
l_max        4 2 1
neuron       60
lambda_1     0
lambda_2     0.04
lambda_e     1
lambda_f     1
lambda_v     0.1
batch        5000
generation   50000
```



##### gpumd

使用nep势采样，MD的初始结构可以变化一些，这里我加大了氢气的密度，选择添加了70个氢分子。

由于势的不完美，MD轨迹跑到后期可能会出现很严重的非物理构型。我们需要截取轨迹

这是一个很简单的脚本，方便复制粘贴

```python
from ase.io import read,write

db = read('dump.xyz','::')
write('waitselect.xyz',db)
```

由于我的初始训练集纰漏很大，所以我觉得第一个turn的可用构型都用来加入训练集，并对他们中的一部分也进行了微扰。

这里是处理轨迹并微扰的脚本，微扰的程度记得根据实际修改，可以通过ovito的RDF计算简单判断一下H-H的范围是否足够多样。

```python
#workflow2
#将前80个结构保存，后20个结构扰动，一共得到120个结构
from ase.io import read,write
from hiphive.structure_generation import generate_mc_rattled_structures


db = read('Hsample5.xyz',':100:')
H2trainset = []
for n,at in enumerate(db):
    if n < 80:
        H2trainset.append(at)
    else:
        H2trainset.extend(generate_mc_rattled_structures(at,2,rattle_std=0.02,d_min=0.65,n_iter=5))

write('H2trainset.xyz',H2trainset)

```

##### neptrain

正常来说，对于gpumd跑出来的构型，应该是使用某种主动学习方法筛选特定结构加入到训练集。

我所知道的方法有，远点采样和委员会法。

而neptrain的select功能可以很方便得实现远点采样。具体来说，就是将gpumd步骤产生的结构按照远点采样原则筛选出一些来。

再加入训练集。

确保当前路径有名为train.xyz的训练集，nep.txt

```
NepTrain select trajectory.xyz -max 100 -o selected.xyz
```

