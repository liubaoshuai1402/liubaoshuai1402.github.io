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

### 提交不同深度的VASP任务

实际的任务路径下，POSCAR等输入文件可能位于不同的深度。

比如，有的在./1/POSCAR，有的则在./1/11/111/POSCAR。

这样简单的循环可能难以完成任务提交，需要一个比较高级的循环方式。

sh脚本

```sh
#!/bin/bash

# find POSCAR
find . -type f -name "POSCAR" | while read -r line; do
    # get dir
    dir=$(dirname "$line")
    
    # cd dir
    cd "$dir" || continue
    
    # just an example
    qsub vasp.pbs
    
    #cd last dir
    cd - >/dev/null || exit
done
```

find 找到当前路径以及所有子路径下名为POSCAR的文件，并依次输出这些文件的相对路径

这其实是一个隐藏的循环，传递给while。

然后用dirname获取文件的父路径，进入，投任务。

这是已经提交了前10个任务，提交第11-20个任务。

```sh
#!/bin/bash
count=0
find . -type f -name "POSCAR" | while read -r line; do
    ((count++))
    if ((count <= 10)); then
        continue
    fi
    if ((count > 20)); then
        break
    fi

    dir=$(dirname "$line")
    cd "$dir" || continue
    qsub vasp.pbs
    sleep 2
    cd - >/dev/null || exit
done
```

#### os.walk('.')

```
import os
from stretch import stretch_strycture

prims = []
for root,dirs,files in os.walk('.'):
    if 'POSCAR' in files:
        prims.append(os.path.join(root,'POSCAR'))
    if 'CONTCAR' in files:
        prims.append(os.path.join(root,'CONTCAR'))
for prim in prims:
    stretch_strycture(prim)
```

os.walk('.')可以很方便地遍历一遍当前路径以及所有子路径的文件。

在这个循环中，

root 为某次循环中进入到的路径，

dirs 为某次循环中进入到的路径下的所有文件夹名称列表，这里是单纯的名称而不是路径，

files 为某次循环中进入到的路径下的所有文件名称列表，这里是单纯的名称而不是路径。

所以判断 `if 'POSCAR' in files:` ，如果为真，这用 `os.path.join()` 将 root 和 POSCAR连接起来获得这个POSCAR的路径。

此外，

os.path.basename(str),只保留最后一级的路径的名称

os.path.dirname(prim),只保留最后一级前的路径
