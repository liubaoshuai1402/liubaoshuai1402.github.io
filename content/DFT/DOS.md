---
title: DOS
date: 2025-09-03
series: ["后处理"]
categories: [DFT]
---

# 电子态密度（DOS）计算

### INCAR

```
Start Parameters for this run
   ISTART = 0
   ICHARG = 2
   PREC = Accurate
   ALGO = Normal
   ISMEAR = -5
   LCHARG = .TRUE.
   LREAL = .FALSE.
   LELF = .TRUE.
   ISYM = 0
   GGA = PS
   LORBIT = 12
   NSW = 0
   NBANDS = 520
   NEDOS = 1000

Electronic minimisation
  ENCUT = 480  
  NELM  = 100
  NELMIN = 4
  EDIFF = 1E-6
  ISPIN = 1   

Ionic relaxation 
  IBRION = -1
```

### DOS的质量

#### k点

k点太稀疏，会导致尖锐的点出现。
