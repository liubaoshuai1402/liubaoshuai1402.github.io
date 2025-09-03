---
title: LAMMPS中的常见命令
date: 2025-06-17
series: ["LAMMPS"]
categories: [分子动力学]
---

# LAMMPS中的常见命令

## 前言

虽然说是常见，但也未必常见吧，可能只是我遇到的不懂的、或者觉得有趣的。

### 1.`labelmap`

`lablemap atom 1 H 2 O`

这个命令用于给atom type指定一个映射关系。在用`write_data`写当前帧的data文件时，文件会含有额外的信息，如下：

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/labelmap.png" style="zoom:60%;" />

但只限于data文件，无法影响dump的轨迹文件，所以想要写带有元素符号的轨迹，还是用dump的custom风格比较好。

### 2. `variable`

想要精通LAMMPS中的`variable`总是要费一番功夫的，今日有幸认真研究一番（本身又涉及很多别的命令，真似高中时看牛津字典，遇到一个个新单词，好爽快）。

不同风格的`variable`的==定义==和==使用==，会有一些差异。

```
variable name style args ...
```

这是官方给出的`variable`语法，简明扼要，`variable`由名字、风格、参数组成。

- `equal`风格

  它后接一个公式，可以包含：数字、常数、数学算符、内置函数、原子值（atom values）、原子矢量（atom vectors）以及compute/fix/variable的引用。

  ```
    atom value = id[i], mass[i], type[i], mol[i], x[i], y[i], z[i], vx[i], vy[i], vz[i], fx[i], fy[i], fz[i], q[i]
    atom vector = id, mass, type, mol, radius, q, x, y, z, vx, vy, vz, fx, fy, f
  ```

  这里，原子值是某个原子的信息，而原子矢量则是包含全体原子的信息。

  值得注意的是，如果涉及到原子值，在`atom_style`之后，还要打开`atom_modify map yes`。

  此外，lammps中`i`是从1开始的，和python的从0开始规则不一样。

  `equal`风格存储了一个公式，当被使用时，将会输出一个单值，它的使用场景：

  1. `print`、`fix print`、`run every`

  

  2. 热动力学输出，由`thermo_style`控制

  ```
  variable alpha equal mass[1]
  thermo          100
  thermo_style    custom step temp pe ke etotal press v_alpha
  ```

  输出如下：

  <img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/v1.png"/>

  3. 平均fix的输入

- `atom`风格

  它的定义是和`equal`风格相同的，不同的是它的使用场景：

  1. dump custom 风格时，作为==每原子（per-atom）==的性质输出。

     ```
     variable        alpha atom x+y+z
     dump	    	1 all custom 100 out1 id element x y z v_alpha
     ```

     把变量定义为`atom`风格，==最大的用处==就是可以直接调用一些per-atom性质，比如原子的x、y、z坐标。

     LAMMPS会知道，你写的x，是指原子的x坐标，并且在用dump custom输出per-atom性质时，还会随着原子的变化而变化。

     利用以上两行，我们可以额外输出每个原子的x、y、z坐标和（虽然可能没有物理意义），如下：

     <img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/v2.png" style="zoom:67%;" />

  2. 平均fix的输入（ [fix ave/chunk](https://docs.lammps.org/fix_ave_chunk.html) and [fix ave/atom](https://docs.lammps.org/fix_ave_atom.html)）

- `vector`风格

  它的定义是和`equal`风格相同的，不同的是它的使用场景：

  1. 它可以作为各种平均fix的输入。

  2. 作为一个矢量，也可以看作python中的列表理解，它的元素可以作为热动力学量输出。

     ```
     variable        alpha vector [1,2,3,4]
     thermo          100
     thermo_style    custom step temp pe ke etotal press v_alpha[4]
     ```

     <img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/v3.png"/>

     这里，对于`vector`风格的`variable`，热动力学输出只接受它的元素，而不能直接输出它本身。

     这里也牵出了对于`vector`风格的数据中的元素的调用方式，就是用`[i]`。

- “immediate” variables

  LAMMPS中的$有大学问。

  1. $可以跟{}，此时{}内的内容需要是一个多字母变量名。

  比如说定义了一个`variable`，名叫name，那么`${name}`就可以调用它。

  这里强调，{}内必须且只能是变量名，不能进行计算，不可以是公式。

  2. $可以跟()，这里才是所谓的定义一个“immediate” variables，并且它是`equal`风格的。

  ()内可以进行计算，也可以引用变量（需要用v_name或内置变量）

  比如，在用nvt系综时，热浴的温度阻尼常用100个时间步。，当我们`timestep     0.001`定义完时间步以后，就会有一个内置变量dt。

  我们可以用如下命令

  `fix 1 all nvt temp 300.0 300.0 $(100.0*dt)`

  但如果用：`fix 1 all nvt temp 300.0 300.0 ${100.0*dt}`，就会报错。

  =={}和()是完全不一样的概念，一定不要搞混。==

  3. $可以直接跟一个字母，这个字母需要是已定义的单字母变量的名。


### 3. fix ave/time

`fix ID group-ID ave/time Nevery Nrepeat Nfreq value1 value2 ... keyword args ...`

用于统计平均。

这里以`fix ave all ave/time 1 1 1000 v_tdiff ave running`为例子，

Nfreq控制多少步统计一次，

Nrepeat控制某次（每经历Nfreq步时）会有多少个数据用于计算平均，

Nevery控制某次（每经历Nfreq步时）用于取平均的数据的步长间隔。

在这个例子中，其实是每1000步取一个数据。

如果改成1 3 1000，就是每1000步，取第998，999，1000步的数据的平均值。

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/lammps3.png"/>

ave running的作用是：

最后输出的`f_ave`，是所有Nfreq间隔的平均的平均。

### 4. fix ave/chunk

`fix ID group-ID ave/chunk Nevery Nrepeat Nfreq chunkID value1 value2 ... keyword args ...`

这个命令本质上还是求时间的平均，只是可以很方便的对chunk中所有的区域同时计算。

以`fix 2 all ave/chunk 10 100 1000 layers v_temp file profile.mp`为例子，

这里是每一千步，取100个数据的平均，数据间隔为10.

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/lammps4.png" style="zoom:67%;" />

可以看到，这是输出间隔为1000的profile.mp

第三列是原子数，第四列是温度。

这些值其实是这1000步内的平均，平均方法是10，100，1000.

