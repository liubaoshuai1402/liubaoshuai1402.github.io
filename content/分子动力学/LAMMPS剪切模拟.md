---
title: LAMMPS剪切模拟
date: 2025-06-14
series: ["LAMMPS"]
categories: [分子动力学]

---

# LAMMPS用于剪切模拟的in文件

```
# 3d metal shear simulation

units		metal
boundary	s s p

atom_style	atomic
lattice		fcc 3.52
region		box block 0 16.0 0 10.0 0 2.828427
create_box	3 box

lattice		fcc 3.52 orient	x 1 0 0 orient y 0 1 1 orient z 0 -1 1 &
		origin 0.5 0 0 
create_atoms	1 box

pair_style	eam
pair_coeff	* * Ni_u3.eam

neighbor	0.3 bin
neigh_modify	delay 5

region		lower block INF INF INF 0.9 INF INF
region		upper block INF INF 6.1 INF INF INF
group		lower region lower
group		upper region upper
group		boundary union lower upper
group		mobile subtract all boundary

set		group lower type 2
set		group upper type 3

# void

#region		void cylinder z 8 5 2.5 INF INF
#delete_atoms	region void

# temp controllers

compute		new3d mobile temp
compute		new2d mobile temp/partial 0 1 1

# equilibrate

velocity	mobile create 300.0 5812775 temp new3d
fix		1 all nve
fix		2 boundary setforce 0.0 0.0 0.0

fix		3 mobile temp/rescale 10 300.0 300.0 10.0 1.0
fix_modify	3 temp new3d

thermo		25
thermo_modify	temp new3d

timestep	0.001
run		100

# shear

velocity	upper set 1.0 0 0
velocity	mobile ramp vx 0.0 1.0 y 1.4 8.6 sum yes

unfix		3
fix		3 mobile temp/rescale 10 300.0 300.0 10.0 1.0
fix_modify	3 temp new2d

#dump		1 all atom 100 dump.shear

#dump		2 all image 100 image.*.jpg type type &
#		axes yes 0.8 0.02 view 0 0 zoom 1.5 up 0 1 0 adiam 2.0
#dump_modify	2 pad 4

#dump		3 all movie 100 movie.mpg type type &
#		axes yes 0.8 0.02 view 0 0 zoom 1.5 up 0 1 0 adiam 2.0
#dump_modify	3 pad 4

thermo		100
thermo_modify	temp new2d

reset_timestep	0
run		3000

```

### `compute  ID  group-ID  temp`

计算某个原子群的温度

- compute		new3d mobile temp

  一个名为`new3d`的compute，计算`mobile`的温度。

- compute		new2d mobile temp/partial 0 1 1

  因为lammps的温度与动能挂钩，可以实现只统计特定方向的动能来计算温度。

  这段代码，指的是只统计y、z方向的动能并计算温度（不统计x方向，所以对应 0）

### `velocity`

velocity	mobile ramp vx 0.0 1.0 y 1.4 8.6 sum yes

velocity      group-ID ramp style 

设置速度均匀变化。vx 0.0 1.0 指x方向的速度从0均匀变化到1.0； y 1.4 8.6 是指当y方向坐标从1.4变化到8.6时，速度发生变化。sum yes 表示这个命令产生的是速度会加到之前已有的速度上，而不是取代。

- temp 关键词

  用法：temp value = temperature compute ID

  例子1：velocity	mobile create 300.0 5812775 temp new3d

  因为速度和温度挂钩，使用`velocity`时，可以自定义温度的计算方式，就是使用temp关键词。

  例子2：

  compute		new2d mobile temp/partial 0 1 1

  fix		3 mobile temp/rescale 10 300.0 300.0 10.0 1.0
  fix_modify	3 temp new2d

  这三段代码可以实现对y、z方向的控温。fix涉及温度时，也可以用fix_modify自定义温度。

  