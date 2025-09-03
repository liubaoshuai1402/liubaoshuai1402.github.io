---
title: 使用rNEMD方法计算热导率的lammps输入文件
date: 2025-05-16
series: ["LAMMPS"]
categories: [MolecularDynamics]
---

# 使用rNEMD方法计算热导率的LAMMPS输入文件

## 前言

[rNEMD方法](https://doi.org/10.1063/1.473271)，又叫MP方法，计算材料热导率。

LAMMPS官方提供了计算脚本，但是使用的单位却是`lj`单位制，非常不实用，这里是我自己写的`metal`单位制下的脚本。

**注意**：本文仅供参考，欢迎指出错误或分享补充。无能力提供任何指导，**求教者切勿留言**。

## in file

```
# sample LAMMPS input script for thermal conductivity
# Muller-Plathe method via fix thermal_conductivity

# settings      temperature, kB and timestep

variable        t equal 300
variable        k equal 8.6173e-5 
variable        dt equal 0.0005

# convert from LAMMPS metal units to SI
variable        eV2J equal 1.6022e-19          #energy convert
variable        A2m equal 1.0e-10           #distance convert
variable        ps2s equal 1.0e-12          #time convert  
variable        convert equal ${eV2J}/${ps2s}/${A2m}          

# setup problem

units           metal
atom_style      atomic
atom_modify     map yes
newton          on
read_data       222010

pair_style    	mliap unified /home-ssd/Users/nsgm_zcx/macetrain17/YSZH_MACE_model.model-mliap_lammps.pt 0
pair_coeff    	* * H O Y Zr
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

#fix             R all box/relax aniso 0.0 vmax 0.001
minimize        0 1e-8 1000 100000
timestep        0.0005
velocity        all create 10 12345 dist gaussian mom yes rot yes

# npt increase temp
reset_timestep  0
fix 		remove_com all momentum 500 linear 1 1 1
fix             1 all npt temp 10 $t 0.1 x 0.0 0.0 1.0 y 0.0 0.0 1.0 z 0.0 0.0 1.0 couple xy
thermo_style    custom step temp pe etotal enthalpy lx ly lz vol press
thermo          100
dump            1 all custom 100 nvttrj id element x y z fx fy fz
dump_modify     1 sort id element H O Y Zr
run             20000
undump          1
velocity        all scale $t
unfix           1

# npt keep temp
reset_timestep  0
fix             1 all npt temp $t $t 0.1 x 0.0 0.0 1.0 y 0.0 0.0 1.0 z 0.0 0.0 1.0 couple xy
thermo_style    custom step temp pe etotal enthalpy lx ly lz vol press
thermo          100
dump	    	1 all custom 100 yszhtrj id element x y z fx fy fz 
dump_modify	1 sort id element H O Y Zr
run             10000
undump		1
velocity        all scale $t
unfix           1

#nvt keep temp
reset_timestep  0
fix             1 all nvt temp $t $t 0.1 
thermo_style    custom step temp pe etotal enthalpy lx ly lz vol press
thermo          100
dump	    	1 all custom 100 yszhtrj id element x y z fx fy fz 
dump_modify	1 sort id element H O Y Zr
run             10000
undump		1
velocity        all scale $t
unfix           1



# 2nd equilibration run
compute         ke all ke/atom
variable        temp atom c_ke/1.5/${k}

fix             1 all nve

compute         layers all chunk/atom bin/1d z lower 0.05 units reduced
fix             2 all ave/chunk 10 100 1000 layers v_temp file profile.mp
fix             3 all thermal/conductivity 20 z 20

variable        tdiff equal f_2[11][3]-f_2[1][3]
thermo_style    custom step temp epair etotal f_3 v_tdiff
thermo_modify   colname f_3 E_delta colname v_tdiff dTemp_step

thermo          1000
run             80000

# thermal conductivity calculation
# reset fix thermal/conductivity to zero energy accumulation
fix             3 all thermal/conductivity 20 z 20

variable        start_time equal time
variable        kappa equal (f_3/(time-${start_time})/(lx*ly)/2.0)*(lz/2.0)/f_ave

fix             ave all ave/time 1 1 1000 v_tdiff ave running
thermo_style    custom step temp epair etotal f_3 v_tdiff f_ave
thermo_modify   colname f_3 E_delta colname v_tdiff dTemp_step colname f_ave dTemp

run             20000
print           "Running average thermal conductivity units metal: $(v_kappa:%.2f)"
variable        tc  equal ${kappa}*${convert}
print           "Running average thermal conductivity units SI: $(v_tc:%.2f)"
```

## 细节

### fix             3 all thermal/conductivity 20 z 20

这个fix可以实现动量交换，这里`Nstep=100`，代表每100步进行一次交换，这个值越小，交换频率也越大，施加在物体的温差也越大。

可以修改这个Nstep来达到自己满意的温差。

后续也可以添加关键词`swap`，他代表每次交换动能的原子数，默认是1.这个值越大，每次交换的能量也越大，施加在物体的温差也越大。

所以他俩对温差的影响是相反的。

### 平衡问题

要设置合理的模拟时间，尽可能长一些。

MP方法计算热导率，除去最开始的体系平衡，还剩下两个重要阶段。

1. 建立热流

   也就是`fix nve`后的第一个run。

   thermo会输出温差dTemp_step，计算时观察这个值，让他不在单调升高，出现稳定是最好的。

2. 动量归零，开始计算热导率

   也就是`fix nve`后的第二个run。

我们来看一下lammps官方的例子中，两个阶段的特点。

这是一阶段的温差随时间的变化，

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/MP2.png" style="zoom:10%;" />

可以看到，第一阶段的目的主要是为了==建立起逐渐趋于稳定的温差==。

这是第二阶段的温差随时间的变化，

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/MP3.png" style="zoom:10%;" />

可以看到，在官方的演示中，平衡后的二阶段，也会有或者说==允许一些温差的升高==，但总体变化不大。

==个人观点：==

第一阶段的步数可以取大一点，直到观察到了温差“饱和”现象。

第二阶段不必要和第一阶段等长。本质是稳定后，统计一段时间的平均。

### Nstep与温差大小、良好温度梯度的关系



从计算公式中，可以看出，KAPPA * ΔT  正比于 ΔE（动能交换），

那么我们可以假设体系要计算的热导率其实是固定的，我们动能交换的越频繁，ΔE越大，温差也就越大。

如果想要小的温差，只需要降低频率，增大Nstep就可以了。

但根据实践，Nstep取太大会导致slab间的温度梯度建立的效果很差。

其实很好理解，因为温差是通过`fix 3 all thermal/conductivity 20 z 20`建立的，设想一下，如果20换成2w，

那强制交换动能形成温差后，又会在接下来2w步中逐渐通过热平衡消除温度梯度，

所以Nstep不能太大，起码要大于自然热平衡的效果。

现在是Nstep取得20，第一阶段跑8w步，可以实现一阶段温差呈趋近效果，温度梯度建立得也比较良好。





