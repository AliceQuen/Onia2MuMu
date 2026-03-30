# Historical Project Notes

This file preserves project-history material from the previous `README.md`. The dated planning items and branch descriptions below are historical notes, not the current maintenance guide.

## 短期工作

### $J/\psi + J/\psi + \Upsilon$ 相关分析

- 2024.11.5 （主负责：王驰）
  - 检查Run2023C-v2系列数据处理结果
  - 应用`CutOptimization`工具对Run2023系列数据筛选进行优化，尝试提高显著度

### $J/\psi + J/\psi + \phi$ 相关分析

- 2024.11.5 （主负责：程幸）
  - 修改代码，使得可以从数据中提取$J/\psi + J/\psi + \phi$过程候选

### $J/\psi + \Upsilon + \phi$ 相关分析

- 2024.11.5 （主负责：石镇鹏）
  - 修改代码，使得可以从数据中提取$J/\psi + J/\psi + \phi$过程候选

## 已完成工作

### 2024.10.31

- 初步开发得到处理$J/\psi + J/\psi + \Upsilon$过程的代码。已经用于处理Run 3在2023年产生的数据。
- 修改README.md文件，增加了对于代码的使用说明。

## `git`结构说明

### `master`分支

仅形式化地保留，不进行更新

### `Dev-J-J-U`分支（主要作者：王驰）

作为开发分支，用于开发$J/\psi + J/\psi + \Upsilon$过程的代码。

相对“保守”，不会进行过多的实验性尝试。主要通过从其他分支合并代码来更新。

此部分已经相对完备，并且已经开始在Run 3的2023年数据集上运行。

### `Dev-J-J-P`分支（主要作者：程幸）

作为开发分支，用于开发$J/\psi + J/\psi + \phi$过程的代码。

相对“保守”，不会进行过多的实验性尝试。主要通过从其他分支合并代码来更新。

此部分尚待进一步完善。

未来将会在Run 2的数据集上进行初步测试，与日内瓦大学S. Leontsinis等人的结果进行对比，以验证代码的正确性。

### `Dev-J-U-P`分支（主要作者：石镇鹏）

作为开发分支，用于开发$J/\psi + \Upsilon + \phi$过程的代码。

相对“保守”，不会进行过多的实验性尝试。主要通过从其他分支合并代码来更新。

此部分尚待进一步完善。可以结合`Dev-J-U-P`和`Dev-J-J-U`分支代码进行开发。

## 主要贡献者

* 王驰(Eric100911)：清华大学致理书院2022级本科生
* 程幸(Endymion)：清华大学致理书院2022级本科生
* 石镇鹏：清华大学物理系2023级本科生

指导教师：胡震副教授
