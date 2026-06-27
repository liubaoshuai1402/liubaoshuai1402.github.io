---
title: COHP
date: 2025-06-28
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
COHPendEnergy	  6
usebasisset pbeVaspFit2015
gaussianSmearingWidth 0.05 # 这里就是前面提到如果你INCAR里面ISMEAR=0或者1，这里就要设置gaussianSmearingWidth
basisfunctions Fe 3d 4s 4p
basisfunctions Ni 3d 4s 4p
cohpbetween atom 1 and atom 2
cohpGenerator from 2.483 to 2.485 type Ga type As
```

通常，用户只需要指定能量窗口（费米能级以下的是实际占据的态，所以建议正值小点儿，负值大点儿）、`local basis`以及计算方式。

LOBSTER是一个把计算的物理结果（电子密度）转化为化学结果（轨道电子填充）的软件。

- `local basis`就像是在准备杯子，告诉LOBSTER应该如何把电子密度分割，分割到哪里。

指定方法也很简单，就是用类似的语句，`basisfunctions Fe 3d 4s 4p`。

- LOBSTER可以按照pair distance来进行计算，即考虑一定区间内所有pair的cohp，同时也可以指定原子类型。

也可以指定特定两个原子间。

`cohpbetween atom 1 and atom 2 orbitalWise`
`cohpGenerator from 2.483 to 2.485 type Ga type As`

`orbitalWise`关键字用于打开“沿轨道”。

比如这里的atom1是Bi，atom2是O。

打开沿轨道后，我们就可以计算Bi的s、p轨道电子与O的s、p轨道电子的cohp，进一步分析两者的成键由哪些电子主导。

如图，

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/cohp.png"/>

虽然LOBSTER的orbitalWise会把p分解成px、py、pz，但后处理时可以加和起来。

这张图中，不同的颜色代表不同的轨道对儿的COHP，比如说绿色代表Bi_6p与O_2p的COHP。

- 并行

LOBSTER默认尽可能多的使用核数，且不支持多节点并行。

Window版本，vscode中的gitbash终端下运行，终端不能关闭。

nohup.exe lobster-5.1.1-win64.exe > out.log 2>&1 &

- 输出文件-COHPCAR.lobster

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/cohp1.png"/>

从第三行开始，记录了计算的标记以及具体的值

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/cohp2.png"/>

第一列是能量，第二列、第三列分别是第一个标签对应的pCOHP，ICOHP，第四、第五列分别是第二个标签对应的pCOHP，ICOHP。





### LOPOSTER









### LobsterPy

一个用于后处理LOBSTER输出文件的python软件包。

LobsterPy自动键分析的核心算法依赖于LobserEnv。

不同于传统的Env，（如Voronoi，筛选近邻依靠角和距离的截断），LobserEnv只靠ICOHPs来确定原子的近邻。

备忘：目前我个人lobsterpy安装到了名叫ase的conda虚拟环境中，记得把它的脚本路径加入到环境变量中。

- 简单后处理

`lobsterpy.exe plot 1 2`可以画出前两个标签的energy vs pCOHP。也可以不要2。











