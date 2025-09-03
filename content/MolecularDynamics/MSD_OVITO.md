---
title: 均方位移（MSD）-OVITO
date: 2025-05-13
series: ["OVITO"]
categories: [分子动力学]
---

# 均方位移（MSD）计算 by OVITO

## 前言

[OVITO Python Reference — OVITO Python Reference 3.12.3 documentation](https://docs.ovito.org/python/index.html) 是一个开源且功能强大的分子动力学后处理软件包。

本文将介绍如何利用 OVITO python module 计算某类元素原子在一段轨迹内的均方位移。

适用于无机非晶体，其他体系慎用。

软件：OVITO、matplotlib、numpy

**注意**：本文仅供参考，欢迎指出错误或分享补充。无能力提供任何指导，**求教者切勿留言**。

## OVITO版

### 代码展示

```python
from ovito.io import import_file, export_file
from ovito.modifiers import CalculateDisplacementsModifier
from ovito.modifiers import SelectTypeModifier,InvertSelectionModifier,DeleteSelectedModifier,ExpressionSelectionModifier
import numpy as np
import matplotlib.pyplot as plt

#万物起源 import_file ，导入一段要计算的轨迹
pipeline = import_file("1.dump")

#添加 SelectTypeModifier 修饰器
#设置参数 property = 'Particle Type' 指定选择的类型（这里我们指定的是原子类型）
#设置参数 types = {4} 指定具体的原子类型，这里是这个轨迹中的 4 原子（这个值要根据自己的体系修改），这里用数字代表原子是因为我使用的lammps的输出风格没有记录原子的元素符号，如果你的轨迹里记录的是 元素符号 信息，比如说 VASP 输出的 XDATCAR，则需要用类似于 types = {"H"} 的写法
pipeline.modifiers.append(SelectTypeModifier(property = 'Particle Type', types = {4}))

#添加 InvertSelectionModifier 修饰器，进行原子反选，为剔除不需要计算的原子做准备
pipeline.modifiers.append(InvertSelectionModifier())    

#添加 DeleteSelectedModifier 修饰器，删除上一行代码反选的原子，留下需要计算的原子
pipeline.modifiers.append(DeleteSelectedModifier(operate_on= {'particles'}))                    

#添加 CalculateDisplacementsModifier 修饰器，指定计算 MSD 的参考结构，这里 reference_frame = 0 代表初始结构是参考结构
reference_frame = 0
pipeline.modifiers.append(CalculateDisplacementsModifier(reference_frame=0))    #a subclass of ovito.pipeline.ReferenceConfigurationModifier

#自定义一个修饰器函数，用于将 per-particle displacement 转化为相应元素的均方位移
#本文的 OVITO小知识 将简单介绍自定义修饰器是如何工作的
def calculate_msd(frame, data):
    
	#用一个变量 displacement_magnitudes 记录 data.particles['Displacement Magnitude']，简化代码
	displacement_magnitudes = data.particles['Displacement Magnitude']
	#计算 MSD （将所有原子位移的平方加和然后求平均），OVITO 的数据可以直接和 numpy 交互，nice
	msd = np.sum(displacement_magnitudes ** 2) / len(displacement_magnitudes)           
	#把计算的 MSD 传递给 data (DataCollection类)
	data.attributes["MSD"] = msd 
    
#添加自定义 calculate_msd 修饰器
pipeline.modifiers.append(calculate_msd)


#计算 Pipeline, 得到time vs MSD的数据
table = []     #用于存放数据，time vs MSD

for frame,data in enumerate(pipeline.frames):
	if frame >= reference_frame:
		#这里的 *10 一定要根据自己的计算调整，我的轨迹在lammps计算设置：时间步是0.5fs，每20步输出一帧，所以轨迹中每帧其实经历了10fs，所以乘以10
		#我们 time vs MSD 的x横坐标单位是fs，也可以是别的，自己调整
		time = (frame-reference_frame)*10                 
		table.append([time,data.attributes['MSD']])

#.csv文件还是比较高级的，比纯txt好些，delimiter 指定间隔符为 "," ,这样方便直接excel打开
np.savetxt("msd_data.csv",table,delimiter=",")
```

## OVITO小知识

受限于Python基础和时间精力的限制，以下内容皆为我个人的有限理解，未能严格考究，仅供参考。

### 自定义修饰器

这里介绍的自定义修饰器，按官方说法是 Simple programming interface，它是一个函数，接受两个基本输入：`frame` 和 `data`。

即，`def modify(frame: int, data: DataCollection):`

**注意**：不要返回任何值，数据都应保存在 data 这个 **DataCollection** 类中

```python
def calculate_msd(frame, data):
    
	displacement_magnitudes = data.particles['Displacement Magnitude']
	msd = np.sum(displacement_magnitudes ** 2) / len(displacement_magnitudes)           
	data.attributes["MSD"] = msd 
```

`DataCollection` 类的 `attributes` 属性储存了这个实例数据集的所有的`global attributes`（全局信息）。

`attributes` 属性是由经过`@propert`装饰器装饰的方法，这个方法返回一个辅助类，也可以叫功能类（`_AttributesView`）的实例，`_AttributesView`是抽象基类`MutableMapping`的子类。

`_AttributesView`类用于实现类似字典的功能。

每次访问`attributes` 属性时，都会返回一个`_AttributesView`类的实例，这个实例会先通过`__init__`接收`DataCollection`实例的数据，并支持对其进行字典操作。

所以，当我们`data.attributes["MSD"] = msd`时，是通过`_AttributesView`类提供的字典操作，把MSD数据添加到了`DataCollection`实例中。

由于`_AttributesView`类实现了`__repr__`方法，我们可以直接`print(data.attributes)`来查看包含了哪些全局信息。

### `for frame,data in enumerate(pipeline.frames):`语句

`DataCollection`要`pipeline`经过`compute()`方法得到，但这个例子中并没有见到`compute()`方法。其实隐藏在了这句`for`循环中。

`pipeline.frames`会返回一个迭代器，这个迭代器中`yield`产生每一帧的`DataCollection`。

## References

1.[OVITO官方的MSD脚本](https://www.ovito.org/manual/python/introduction/examples/modifiers/msd_calculation.html)