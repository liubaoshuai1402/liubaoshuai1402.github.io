---
title: 机器学习势MACE的输入文件
date: 2025-05-16
series: ["MACE"]
categories: [机器学习势]
---

# 机器学习势MACE的输入文件

## 前言

不同MACE版本的参数设置会有一定的调查，注意查看自己的MACE版本。这里是0.3.13版本

**注意**：本文仅供参考，欢迎指出错误或分享补充。无能力提供任何指导，**求教者切勿留言**。

## 在超算上用slurm提交python任务

```sh
#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -J lmp_job_gpu
#SBATCH --partition=gpu        
#SBATCH -N 1                    
#SBATCH -G 1
#SBATCH --gres=gpu:1
#SBATCH --exclusive


export PATH=/home-ssd/Users/nsgm_zcx/openmpi-5.0.5/build/bin:$PATH
export LD_LIBRARY_PATH=/home-ssd/Users/nsgm_zcx/openmpi-5.0.5/build/lib:$LD_LIBRARY_PATH
export PATH=/home-ssd/Users/nsgm_zcx/lammps-develop/build-mliap:$PATH

source /home-ssd/Users/nsgm_zcx/miniconda3/etc/profile.d/conda.sh
conda activate cuda

mace_run_train --config parameters.yaml
```

前面几行是为了激活lammps和openmpi，在这里没啥用。

重要的是要在sh脚本里激活conda，然后用`mace_run_train`命令行脚本来提交训练任务。

用`--config`参数和一个yaml文件来提供训练参数。

以下是`parameters.yaml`的内容

```yaml
name: YSZH_MACE_model
seed: 123
device: cuda
train_file: train.xyz
valid_fraction: 0.2
test_file: test.xyz
compute_forces: True
compute_stress: True
energy_key: energy_vasp
forces_key: forces_vasp
stress_key: stress_vasp
E0s: 'isolated'
hidden_irreps: '64x0e + 64x1o'
r_max: 4.0
batch_size: 20
max_num_epochs: 600
swa: True
start_swa: 480
ema: True
ema_decay: 0.99
default_dtype: float32
lr: 0.01
scaling: rms_forces_scaling
multiheads_finetuning: False
enable_cueq: True
```

要使用`mliap`，就必须打开参数`enable_cueq: True`，并确保python环境中安装了`cuEquivariance`和`cupy`。
