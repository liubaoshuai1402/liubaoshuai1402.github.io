---
title: COHP
date: 2025-06-28
series: ["后处理"]
categories: [DFT]
---

# COHP

## 前言



### VASP文件的准备

LOBSTER要求ISYM=0或-1。其他参数没有特殊要求。

- 结构优化

  略过

- 静态自洽

  以较为简单的K点设置计算波函数、电荷密度。

  用于得到WAVECAR、CHGCAR

  ```
  #关键参数
  ISTART  = 0            
  ICHARG  = 2
  NSW = 0
  NBANDS = 500	#这个值要大，多试试
  NEDOS = 1000
  ```

  后续的非静态自洽的NBANDS和NEDOS尽量和这里保持一致。

- 非静态自洽

  不再更新波函数和电荷密度，进行一次单次的求解。

  这里可以设置更密的K点来求态密度，也可以设置特殊路径的K点来求能带结构。

  ```
  #关键参数
  ISTART = 1     # 使用前一轮的波函数
  ICHARG = 11    # 使用前一轮的电荷密度，不再更新
  NSW = 0
  LORBIT = 12
  NBANDS = 500	#这个值要大，多试试
  NEDOS = 1000
  ```

  特别说一下LORBIT参数，这个参数是后处理参数，和VASP的电子计算无关。

  指示VASP生成：DOSCAR and lm-decomposed PROCAR + phase factors (not recommended)

  所以只需要最后非静态自洽的时候设置即可。

另外，如果对精度不高，或者不考虑成本，我看也有人直接用静态自洽的结果（这里直接设置比较高的K点密度）进行LOBSTER。

如果用的`ISMEAR = 0`，后续要在lobsterin文件中指出。

### LOBSTER

LOBSTER程序完全被一个`lobsterin`文件控制。

```
COHPstartEnergy  -10
COHPendEnergy    5
usebasisset pbeVaspFit2015
gaussianSmearingWidth 0.05 # 这里就是前面提到如果你INCAR里面ISMEAR=0或者1，这里就要设置gaussianSmearingWidth
basisfunctions Ga
basisfunctions As
cohpGenerator from 2.483 to 2.485 type Ga type As
```



### LobsterPy

一个用于后处理LOBSTER输出文件的python软件包。

LobsterPy自动键分析的核心算法依赖于LobserEnv。

不同于传统的Env，（如Voronoi，筛选近邻依靠角和距离的截断），LobserEnv只靠ICOHPs来确定原子的近邻。

LobsterEnv needs only one parameter that influences neighbor selection

