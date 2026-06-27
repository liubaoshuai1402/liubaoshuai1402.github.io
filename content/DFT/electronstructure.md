---
title: 电子相关计算总结
date: 2026-04-23
series: ["后处理"]
categories: [DFT]
---



# 电子结构计算

这里总结一下vasp的电子相关的计算流程，包括DOS、band、bader charge、ELF

## 前置流程

- 结构优化

```
#Start Parameters
PREC = A
ISTART = 0
ICHARG = 2
ALGO = Normal
LPLANE = .TRUE.
GGA = PE
ISPIN = 1
LWAVE = .FALSE.
LCHARG = .FALSE.
IVDW = 12
LASPH = .TRUE.
NCORE = 4

#Electronic Relaxation
ENCUT = 400
EDIFF = 1E-6
NELM = 60
NELMIN = 4
NELMDL = -3
LREAL = AUTO


#Ionic Relaxation
EDIFFG = -1E-2
NSW = 100
IBRION = 2
ISIF = 3

ISMEAR = 0
SIGMA = 0.02
```

## DOS

电子结构的计算一般都包括，静态自洽和非静态自洽。DOS计算需要注意的参数有，

NBANDS，这个值一般是体系电子数/2，然后乘以1.3。

ISMEAR，一般取-5，但是图会比较尖锐，需要取更高密度的K点补偿。如果用0配合SIGMA，图会好看一些，但可能不够精准。

静态自洽的输入文件：

```
#Start Parameters
#use an easy k-point grid
ISTART = 1
ICHARG = 2
PREC = A
ALGO = Normal
LPLANE = .TRUE.
GGA = PE
ISPIN = 2
LWAVE = .TRUE.
LCHARG = .TRUE.
IVDW = 12
LASPH = .TRUE.
NCORE = 4
#bader charge
LAECHG=.TRUE.
#bader charge

#Electronic Relaxation
ENCUT = 500
EDIFF = 1E-6
NELM = 60
NELMIN = 4
NELMDL = -3
LREAL = AUTO


#Density of States
#twice the number of electrons is ok
NBANDS = 72
#twice the number of electrons is ok
NEDOS = 1000
LORBIT = 12

#Ionic Relaxation
EDIFFG = -1E-2
NSW = 0
IBRION = -1
ISIF = 2

ISMEAR = 0
SIGMA = 0.02
```

非静态自洽的输入文件：

```
#Start Parameters
#use more dense k-point grid than electron_1
ISTART = 1
#keep the charge density
ICHARG = 11
#keep the charge density
PREC = A
ALGO = Normal
LPLANE = .TRUE.
GGA = PE
ISPIN = 2
LWAVE = .TRUE.
LCHARG = .TRUE.
IVDW = 12
LASPH = .TRUE.
NCORE = 4

#Electronic Relaxation
ENCUT = 500
EDIFF = 1E-6
NELM = 60
NELMIN = 4
NELMDL = -3
LREAL = AUTO


#Density of States
#twice the number of electrons is ok
NBANDS = 72
#twice the number of electrons is ok
NEDOS = 1000
LORBIT = 12

#Ionic Relaxation
EDIFFG = -1E-2
NSW = 0
IBRION = -1
ISIF = 2

ISMEAR = 0
SIGMA = 0.02
```

但实际体验下来，直接第一步，静态自洽的时候，把K设置密集就可以了。

## band

band计算是对于刚入门没有物理基础的人而言，非常容易出错（说的就是我）。

其中，高对称路径的选取让我非常头疼。特别是对于扩胞之后，到底该如何选取。

> [!important]
>
> 结论很简单，原胞对应原胞高对称路径，晶胞对应晶胞高对称路径。

此外，还有一个非常容易出错的点，同一个晶体结构，它的==基矢==可能有多种表达形式。

- 我们对基矢的选择，会影响高对称点的坐标！
- 我们对基矢的选择，会影响高对称点的坐标！
- 我们对基矢的选择，会影响高对称点的坐标！

这也是为什么vaspkit生成高对称路径的功能，会要求把它提供的原胞复制成POSCAR。

我们在MP下载下来的或者通过软件转化得到的晶体结构，未必是seekpath、vaspkit生成的高对称点所对应的标准结构！

如果我们想用这些软件或网站生成的高对称点，就必须使用它们提供的标准结构！

通常电子结构计算是鼓励使用原胞的，但做掺杂或者超晶格时，不得不用晶胞。

这时候需要结合[文献1](https://doi.org/10.1016/j.commatsci.2010.05.010)和[文献2](https://doi.org/10.1016/j.commatsci.2016.10.015)中的标准结构、标准高对称点、标准路径以及自己体系研究中的常用路径，综合考量确定。

注：早期文献上高对称点与seekpath的高对称点在字母的命名上可能有差异，这个无所谓，只要是布里源区内对应的即可。











## bader charge

计算bader charge的文件一般和电子结构计算的一样，只是需要打开参数：

`LAECHG=.TRUE.`

然后用官方网站的脚本处理一下即可。

需要注意的是，bader charge对fft网格的密度有要求，所以要加查一下收敛情况。

设置参数

NGXF、NGYF、NGZF（用于CHGCAR的fft网格），可以变密集一点，直到收敛。

`ADDGRID=.TRUE.`在计算电荷时可以打开。

（PS：NGX、NGY、NGZ是用于WAVECAR的，不用改）

查看，
