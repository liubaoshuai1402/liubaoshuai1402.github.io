---
title: Python中的文件操作
date: 2025-05-18
categories: [代码]
---

# Python中的文件操作

### 获取当前路径

`os.getcwd()`可以返回当前工作目录的绝对路径。cwd可以理解成current work directory。

### 创建多个用于VASP计算的文件夹

任务：假设你要对多个相似的结构进行计算，即他们所需的INCAR、KPOINTS、POTCAR、vasp.pbs（提交任务的脚本）是一样的。

INCAR、KPOINTS、POTCAR、vasp.pbs在当前文件夹中

我们要在当前文件夹中创建一个名为`'distorted'`的文件夹，并在其中创建若干个计算文件夹。

代码如下：

```python
import os
import shutil

#先记录下四个计算文件的路径
cwd = os.getcwd()
INCARPath = os.path.join(cwd,'INCAR')
KPOINTSPath = os.path.join(cwd,'KPOINTS')
POTCARPath = os.path.join(cwd,'POTCAR')
PBSPath = os.path.join(cwd,'vasp.pbs')
#创建一个总路径，用于存放计算文件夹
os.makedirs('distorted')
path = os.path.join(cwd,'distorted')
#通过循环，在distorted路径下创建若干个计算文件夹，并把计算文件复制过去
for i in range(10):
    CalPath = os.path.join(path,'{:03d}'.format(i))
    os.makedirs(CalPath)
    shutil.copy(INCARPath,CalPath)
    shutil.copy(KPOINTSPath,CalPath)
    shutil.copy(POTCARPath,CalPath)
    shutil.copy(PBSPath,CalPath)
```

