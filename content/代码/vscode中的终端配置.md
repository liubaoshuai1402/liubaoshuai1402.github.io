---
title: 快速搞定vscode中的终端配置
date: 2026-01-17
series: ["python"]
categories: [代码]
---

# 快速搞定vscode中的终端配置

目标：在vscode中流畅地使用powershell+conda

### 在vscode的设置文件中定义终端

如图，

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/vscode1.png"/>

定义一些终端的名字和路径，然后把默认的终端设置为其中一个，这里我们选择的是`PowerShell7`.

### 设置powershell的`profile.ps1`文件

这个文件夹的路径一般在

<img src="https://xiaoxiaobuaigugujiao.oss-cn-beijing.aliyuncs.com/img/vscode2.png" style="zoom:50%;" />

```
#region conda initialize
# !! Contents within this block are managed by 'conda init' !!
If (Test-Path "D:\miniconda3\Scripts\conda.exe") {
    (& "D:\miniconda3\Scripts\conda.exe" "shell.powershell" "hook") | Out-String | ?{$_} | Invoke-Expression
}
#endregion
conda activate doped
```

最后一句，`conda activate doped` 请把环境名换成自己常用的环境名。

> [!important]
>
> 我们这里采用的终端是powershell，而非microsoftpowershell！

### 最后一步

在powershell终端中运行，

```
conda config --set auto_activate false
```

