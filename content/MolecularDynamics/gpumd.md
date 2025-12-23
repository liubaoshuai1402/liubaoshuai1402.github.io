---
title: GPUMD的常见run.in
date: 2025-12-22
series: ["gpumd"]
categories: [MolecularDynamics]
---

# GPUMD的run.in

### 初始化配置

```
potential   ./nep.txt
minimize fire 0.02 20
velocity    1500
time_step   0.5
```

告诉gpumd势文件的路径，结构优化的参数，初始速度，以及模拟的时间步长。

### 设置系综

#### npt

```
ensemble npt_mttk temp 1500 1500 aniso 0 0
```

#### nvt

```
ensemble nvt_bdp 1500 1500 100
```

### 设置计算和dump格式

#### msd

```
compute_msd 10 8000 all_groups 0
```

这个命令计算10*8000个时间步的时间内的msd，0表征第1个分类方法（因为gpumd处理group是从0开始的）

计算该分组方法下所有分组的原子的msd。

比如说YSZH体系，会一下子输出H、O、Y、Zr的。

```
compute_msd 5 200 group 1 1
```

这个命令计算第二个分类方法的group 1 原子的msd

#### dump

```
dump_exyz 8000 1 1 1
```

dump_exyz的格式的三个1分别控制速度、力和原子能量的写入

```
dump_xyz -1 1 8000 mydump.xyz velocity potential force virial
```

dump_xyz的前两个参数控制分组，但是一般我们是全写，所以用-1，第二个参数被忽略。

### run 

```
run 1000
```

## 用于结构优化的例子

```
potential ./nep.txt
minimize fire 0.02 1000
velocity    1

ensemble nve
time_step 0
dump_exyz 1 1 1
run 1
```

