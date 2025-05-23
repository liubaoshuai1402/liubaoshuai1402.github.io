---
title: 使用rNEMD方法计算热导率的lammps输入文件
date: 2025-05-16
series: ["LAMMPS"]
categories: [分子动力学]
---

# 使用rNEMD方法计算热导率的LAMMPS输入文件

## 前言

[rNEMD方法](https://doi.org/10.1063/1.473271)，又叫MP方法，计算材料热导率。

LAMMPS官方提供了计算脚本，但是使用的单位却是`lj`单位制，非常不实用，这里是我自己写的`metal`单位制下的脚本。

**注意**：本文仅供参考，欢迎指出错误或分享补充。无能力提供任何指导，**求教者切勿留言**。

## in file

```
# sample LAMMPS input script for thermal conductivity of liquid LJ
# Muller-Plathe method via fix thermal_conductivity

# settings      temperature, kB and timestep

variable        t equal 1500
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
read_data       ./333
pair_style      mace no_domain_decomposition
pair_coeff * *  /home-ssd/Users/nsgm_lbs/train/MACE_model/MACE_model_run-123_stagetwo.model-lammps.pt O Zr Y H
neighbor        1.0 bin
neigh_modify    every 500 delay 0 check no
minimize        1e-5 1e-7 1000 1000
timestep        ${dt}
velocity        all create $t 87287


# 1st equilibration run
reset_timestep  0
fix             1 all npt temp $t $t 0.05 iso 0 0 0.5
thermo_style    custom step temp pe etotal enthalpy lx ly lz vol press
thermo          100
run             5000
unfix           1
velocity        all scale $t
fix             1 all nvt temp $t $t 0.05
run             5000
unfix           1

# 2nd equilibration run
compute         ke all ke/atom
variable        temp atom c_ke/1.5/${k}

fix             1 all nve

compute         layers all chunk/atom bin/1d z lower 0.05 units reduced
fix             2 all ave/chunk 10 100 1000 layers v_temp file profile.mp
fix             3 all thermal/conductivity 100 z 20

variable        tdiff equal f_2[11][3]-f_2[1][3]
thermo_style    custom step temp epair etotal f_3 v_tdiff
thermo_modify   colname f_3 E_delta colname v_tdiff dTemp_step

thermo          1000
run             20000

# thermal conductivity calculation
# reset fix thermal/conductivity to zero energy accumulation
fix             3 all thermal/conductivity 100 z 20
variable        start_time equal time
variable        kappa equal (f_3/(time-${start_time})/(lx*ly)/2.0)*(lz/2.0)/f_ave
fix             ave all ave/time 1 1 1000 v_tdiff ave running
thermo_style    custom step temp epair etotal f_3 v_tdiff f_ave
thermo_modify   colname f_3 E_delta colname v_tdiff dTemp_step colname f_ave dTemp
run             20000
print           "Running average thermal conductivity units metal: $(v_kappa)"
variable        tc  equal ${kappa}*${convert}
print           "Running average thermal conductivity units SI: $(v_tc:%.2f)"
```

