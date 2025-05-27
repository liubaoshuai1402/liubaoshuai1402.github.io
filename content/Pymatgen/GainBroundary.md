---
title: 晶界建模
date: 2025-05-24
categories: [pymatgen]
---

# 晶界建模

## 前言

基于[共格点阵模型](https://dtrinkle.matse.illinois.edu/MatSE584/kap_7/backbone/r7_1_2.html)（Coincidence Site Lattice），使用pymatgen进行晶界建模。

## 晶界的类型

### 1. 扭转（twist）晶界

旋转轴垂直于晶界面，也就是两者的密勒指数应该成比列。

### 2. 倾斜（tilt）晶界

旋转轴平行于晶界面，也就是两者的密勒指数的点乘为0。

### 3. 混合（mixed）晶界

旋转轴既不垂直也不平行于晶界面。



基于共格点阵模型的晶界命名法，Σ+number(must be odd)+(hkl)/[uvw]。

举一个栗子， Σ13 (510)/[001] symmetric tilt grain-boundary。

这里Σ的大小，是指旋转后重合点阵的单胞的提及是原始晶体单胞的体积的多少倍。通常，这个数字越大，代表两个晶粒的取向相差越远，晶界能也往往越远。

Σ1则代表趋向一致，那些小角度晶界也被认为Σ的值近似于1。

当指明是twist or tilt晶界时，有时可以省略晶向，也不会造成歧义，比如，==the Σ5(310) tilt GB==，这是一个YSZ中典型的低能量晶界。

但它没有给出旋转轴，因为没有必要，tilt GB 要求晶界面与旋转轴平行，所以只能是[001]。

## 使用pymatgen进行晶界建模

首先，我们假设一个应用场景，就是说，我们建模肯定是根据实验来的，实验上对哪些晶界感兴趣，我们就去建模研究。

所以在这个假设的基础上，我们是知道==Σ的值==以及==旋转轴==、==晶界面==的。

这样，用以下代码我们可以得到旋转角。

```python
from pymatgen.core import Structure
from pymatgen.core.interface import GrainBoundaryGenerator

# 1. 读取结构文件
structure = Structure.from_file("ZrO.cif")
structure = structure.to_conventional()

#创建一个晶界生成器，实例化需要一个晶体结构，最好是conventional cell。
gb_gen = GrainBoundaryGenerator(structure)

# 2. 构建 Σ5 晶界，参数分别对应Σ的值、旋转轴、晶格类型，对非立方体系需要指定轴比
#其实这里感觉很奇怪，轴比和晶格类型，pymatgen不应该自己判断吗，感觉这块代码写的不好
rotation_anglen = gb_gen.get_rotation_angle_from_sigma(5,(0,0,1),lat_type='c')
print(rotation_anglen)
```

这里的输出是，[36.86989764584402, 53.13010235415597, 126.86989764584402, 143.13010235415598]。

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/YSZ%E6%99%B6%E7%95%8C.png"/>

第一个值，36.8°与[参考文献]([Sci-Hub | Structure and Chemistry of Yttria-Stabilized Cubic-Zirconia Symmetric Tilt Grain Boundaries. Journal of the American Ceramic Society, 84(6), 1361–1368 | 10.1111/j.1151-2916.2001.tb00842.x](https://sci-hub.st/https://doi.org/10.1111/j.1151-2916.2001.tb00842.x))一致。

同样的，如果我们用`rotation_anglen = gb_gen.get_rotation_angle_from_sigma(13,(0,0,1),lat_type='c')`

则会得到输出，[22.61986494804043, 67.38013505195957, 112.61986494804043, 157.3801350519596]，第一个值与24°差1.4。

```python
#指定旋转轴，旋转角，晶界面
gb = gb_gen.gb_from_parameters(
    rotation_axis=(0,0,1),
    rotation_angle=rotation_anglen[0],
    expand_times=1,
    vacuum_thickness=0,
    plane=(3,1,0)
)
# 3. 获取晶界结构
gb.to("POSCAR", "poscar")

```

再写上这块代码就可以生成晶界了。