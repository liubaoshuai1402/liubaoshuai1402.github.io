---
title: GitHub小贴士
date: 2025-11-13
series: ["command"]
categories: [command]
---

# 命令备忘录

### 网页更新

```shell
git add .
git commit -m 'add a new'
git push
```

### 联系本地仓库和远程仓库

Windows操作系统下。

首先在GitHub上创建一个新的远程仓库，记得勾选readme，自动创建main分支。

然后在本地仓库的路径下，打开git。

```
#本地仓库初始化
git init
#与远程仓库建立连接
git remote add origin https://github.com/liubaoshuai1402/CCkit
#把本地分支名从master改成main
git branch -m master main
#把远程仓库的readme之类的文件拉去下来
git pull --rebase origin main
#把内容加入到本地仓库
git add .
#提交
git commit -m "Initial commit"
#推送
git push -u origin main
```

