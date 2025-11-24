---
title: CalculationOfCij
date: 2025-11-24
series: ["gpumd"]
categories: [MolecularDynamics]
---

# 应力应变法求解弹性系数

现在我要用分子动力学求解有限温度下弹性系数，由于势函数是机器学习势。所以首先，我们应该验证势函数的可靠性。

也就是说，首先验证，在0K下，对象是小胞时，通过DFT与MLIP-MD计算得到的Cij的差距。然后在高温下跑大胞的MD，计算有限温度下的Cij

## 原理及编程

这里以四方晶系-1为例，也就是有六个独立的弹性系数。目前，我所见过的所有，应力应变法求解弹性系数的计算思路都是拟合。

我们先看一下所要求解的弹性系数矩阵的样子。

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/Cij.png" style="zoom: 50%;" />

根据线弹性假设，$$\sigma = C * \varepsilon$$，给出应变，通过线性组合，可以表示出应力，即公式

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/Cij2.png" style="zoom:67%;" />

但是，要想使用numpy.linalg.lstsq进行拟合，还需要重新表达一下这关系。

因为lstsq的工作原理是B = A * X，根据已知应变A、应力B，用最小二乘法拟合一个为列向量的系数，X，（这里即Cij）

现在，我们专注于A矩阵的写法。

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/Cijnote.jpg" style="zoom: 50%;" />

可以看到，==应变A矩阵是由所施加的应变，也即是变形模式，以及弹性系数矩阵的形状共同决定的==。（Cij矩阵决定了A矩阵的形状，而变形模式决定了它具体的元素值）

> [!IMPORTANT]
>
> 有几点需要澄清，待求解的Cij可能不止6个。比如说，7个的时候，需要X列向量是7行，A矩阵是七列吗？显然，不需要也不能。因为一种应变模式只能提供六个线性方程组，6个应力最多只能包含6个弹性系数的有效信息。这也是为什么，低对称晶系需要更多变形模式的原因。因此，我们在选择变形模式的时候，要合理。对于大于6个弹性系数的体系，要让弹性系数间尽量解耦——保证一个变形模式内，涉及的弹性系数最多为6个。通过选取多个变形模式，比如说两种，这样有2个A矩阵，12个方程，把两个A矩阵堆叠起来，就能求解7个Cij了

这里我放上一个具体的脚本，提供一些思路

```python
import numpy as np
#针对四方相计算弹性常数C11,C12,C13,C33,C44,C66
#变形幅度
ups = [-0.03,-0.01,0.01,0.03]
#两种变形模式
modes = [0,1]
Ms = []
for mode in modes:
    for up in ups:
        #第一种变形模式
        if mode == 0:
            strain_voigt = np.array([up,0.,0.,0.,0.,0.])
            strain_matrix = strain_matrix_1 = np.array([[up, 0., 0.],
                                                        [0., 0., 0.],
                                                        [0., 0., 0.]])
            #有了应变模式和应变幅度，即可确定M矩阵
            #M矩阵的列数是求解弹性常数的个数，行数则是6
            M_matrix = np.array([[up,0.,0.,0.,0.,0.],
                                [0.,up,0.,0.,0.,0.],
                                [0.,0.,up,0.,0.,0.],
                                [0.,0.,0.,0.,0.,0.],
                                [0.,0.,0.,0.,0.,0.],
                                [0.,0.,0.,0.,0.,0.]])
            Ms.append(M_matrix)
        #第二种变形模式
        if mode == 1:
            strain_voigt = np.array([0.,0.,up,up,0.,up])
            strain_matrix = strain_matrix_1 = np.array([[0., up/2, 0.],
                                                        [up/2, 0., up/2],
                                                        [0., up/2, up]])
            M_matrix = np.array([[0.,0.,up,0.,0.,0.],
                                [0.,0.,up,0.,0.,0.],
                                [0.,0.,0.,up,0.,0.],
                                [0.,0.,0.,0.,up,0.],
                                [0.,0.,0.,0.,0.,0.],
                                [0.,0.,0.,0.,0.,up]])
            Ms.append(M_matrix)
Ms = np.array(Ms)
Ms_2d = Ms.reshape(-1, 6)  # 自动把 8*6=48 行展开
print(Ms_2d.shape)  # (48, 6)
print(Ms_2d)
```

代码中的M，就是上述的A矩阵。除此之外，代码是先把变形模式1的不同幅值的A矩阵堆叠起来，然后再把变形模式2的不同幅值的A矩阵堆叠起来。

顺序是没有影响的。本质上就是拟合2 \* 4 \* 6 = 48个应力应变关系涉及的6个Cij。就是把矩阵行交换，是不影响有效信息的。

最后，有了堆叠后的应变A矩阵和应力向量，简单的使用lstsq拟合即可。

此外，也可以观察到，这些A矩阵中的每一行中，只有一个不为零的值，其实是因为我们的四方晶系比较简单，合适的应变模式选取让Cij之间完全解耦。

这种情况下，如果不嫌弃麻烦，也可以从这个48个方程中挑选出对应的，某个Cij与有些应力的关系，单独拟合。

## 0K下的弹性系数

0K下DFT计算略过，使用MLIP也不需要跑MD，只需要对变形进行评估得出应力即可。

