---
title: LAMMPS-thermo
date: 2025-08-23
series: ["LAMMPS"]
categories: [分子动力学]
---

# thermo的一些特点

### 出现的位置

虽然fix、thermo、dump、run之间的顺序可能有很多选择，但我更喜欢这个顺序。

### 继承

每个新的run都会默认继承上一个run的thermo配置，如果需要改变，覆盖一下即可。

