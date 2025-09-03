---
title: 超晶格材料原子尺度建模
date: 2025-08-15
categories: [DFT]
---

根据原子的z轴坐标，指定哪些是自由的，哪些是受限制的。

```python
from ase.constraints import FixAtoms
from ase.io import write, read

fixed_indices = []  # 要固定的原子索引
atoms = read('POSCAR_BiTe1_Pt')
for i in range(len(atoms)):
    position_z = atoms[i].position[2]
    if 40 > position_z > 15.5:
        indice = atoms[i].index
        fixed_indices.append(indice)
    if 3 < position_z < 7:
        indice = atoms[i].index
        fixed_indices.append(indice)

atoms.set_constraint(FixAtoms(indices=fixed_indices))
write('POSCAR_BiTe1_Pt_sel', atoms, format='vasp')




```



