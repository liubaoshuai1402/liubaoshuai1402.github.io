---
title: 花朵分类
date: 2025-05-18
categories: [深度学习]
---

# 花朵分类

## 前言

本文主要借助torchvision软件包，简单梳理一下深度学习代码的基本框架。

## 数据集加载

```python
#假设当前工作路径下，存放着一个名为'flower_data'的文件夹，里面存放着训练集和验证集的图片
#用os.path.join()来一级级得获取路径
data_dir = os.path.join(os.getcwd(),'flower_data')
train_dir = os.path.join(data_dir, 'train')
valid_dir = os.path.join(data_dir, 'valid')
```

```python
# Define batch size
batch_size = 32

# Define transforms for the training and validation sets
normalize = transforms.Normalize(mean=[0.485, 0.456, 0.406],
                                 std=[0.229, 0.224, 0.225])

#定义对数据预处理的组合，比如旋转图片、改变尺寸等等，让模型有更强的稳定性
train_data_transforms = transforms.Compose([
        transforms.RandomResizedCrop(size=256, scale=(0.8, 1.0)),
        transforms.RandomRotation(degrees=15),
        transforms.ColorJitter(),
        transforms.RandomHorizontalFlip(),
        transforms.CenterCrop(size=224),
        transforms.ToTensor(),
        normalize,
    ])

validate_data_transforms = transforms.Compose([
        transforms.Resize(256),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        normalize,
    ])


#这里才真正的把图片加载成了二维数据。
train_dataset = datasets.ImageFolder(
    train_dir,
    train_data_transforms)

validate_dataset = datasets.ImageFolder(
    valid_dir,
    validate_data_transforms)



#做了那么多铺垫，其实就是为了把可用于训练的数据(二维数据)放到 DataLoader 里面
train_loader = torch.utils.data.DataLoader(
    train_dataset, batch_size=batch_size, shuffle=True,
    num_workers=4)

validate_loader = torch.utils.data.DataLoader(
    validate_dataset, batch_size=batch_size, shuffle=True,
    num_workers=4)

data_loader = {}
data_loader['train'] = train_loader
data_loader['valid'] = validate_loader
```

`batch_size`，当训练数据很多时，一次性加载全部数据进行训练会是一种挑战。这时就需要用到批次训练。`batch_size`即是一次训练中使用的数据量。注意这里的一次训练，不是指一epoch。只有遍历所有训练集后，才能叫做完成了一代训练。一代训练包含了诸多这样的一次训练。

`shuffle`参数为`True`时，随机采样。



这是CGCNN（晶体卷积图神经网络）的代码中，用于得到训练集、验证集、测试集的`DateLoader`。

```python
def get_train_val_test_loader(dataset, collate_fn=default_collate,
                              batch_size=64, train_ratio=None,
                              val_ratio=0.1, test_ratio=0.1, return_test=False,
                              num_workers=1, pin_memory=False, **kwargs):

    total_size = len(dataset)
    if kwargs['train_size'] is None:
        if train_ratio is None:
            assert val_ratio + test_ratio < 1
            train_ratio = 1 - val_ratio - test_ratio
            print(f'[Warning] train_ratio is None, using 1 - val_ratio - '
                  f'test_ratio = {train_ratio} as training data.')
        else:
            assert train_ratio + val_ratio + test_ratio <= 1
    indices = list(range(total_size))
    if kwargs['train_size']:
        train_size = kwargs['train_size']
    else:
        train_size = int(train_ratio * total_size)
    if kwargs['test_size']:
        test_size = kwargs['test_size']
    else:
        test_size = int(test_ratio * total_size)
    if kwargs['val_size']:
        valid_size = kwargs['val_size']
    else:
        valid_size = int(val_ratio * total_size)
        
    train_sampler = SubsetRandomSampler(indices[:train_size])
    val_sampler = SubsetRandomSampler(
        indices[-(valid_size + test_size):-test_size])
    if return_test:
        test_sampler = SubsetRandomSampler(indices[-test_size:])
    train_loader = DataLoader(dataset, batch_size=batch_size,
                              sampler=train_sampler,
                              num_workers=num_workers,
                              collate_fn=collate_fn, pin_memory=pin_memory)
    val_loader = DataLoader(dataset, batch_size=batch_size,
                            sampler=val_sampler,
                            num_workers=num_workers,
                            collate_fn=collate_fn, pin_memory=pin_memory)
    if return_test:
        test_loader = DataLoader(dataset, batch_size=batch_size,
                                 sampler=test_sampler,
                                 num_workers=num_workers,
                                 collate_fn=collate_fn, pin_memory=pin_memory)
    if return_test:
        return train_loader, val_loader, test_loader
    else:
        return train_loader, val_loader

```

这段代码整体分为两部分，其一是根据参数确定训练集、验证集、测试集的数量，其二是装配`DataLoader`。和花朵分类的代码稍有不同，这样的实现可以让用户自由决定训练集、验证集和测试集。多使用了`DataLoader`的`sampler`参数。

`SubsetRandomSampler`是`torch.utils.data.sampler`模块中的六个类之一。用于随机采样。

这里是常规用法，假设有100个数据，80个用于训练集，15个用于验证集，5个用于测试集。

那么，

训练集的`SubsetRandomSampler`的初始化参数就是`list(range(100))[0:80]`。

验证集的`SubsetRandomSampler`的初始化参数就是`list(range(100))[-20:-5]`。

测试集的`SubsetRandomSampler`的初始化参数就是`list(range(100))[-5:]`。

### 1. Dataset与DataLoader

- 在CGCNN中，作者自己手写了`CIFData`作为自定义的`Dataset`。

  `CIFData`三个重要的魔法方法的功能：

  1.`__init__()`方法

  ```python
      def __init__(self, root_dir, max_num_nbr=12, radius=8, dmin=0, step=0.2,
                   random_seed=123):
  ```

  这里简单初始化了：数据集的路径、晶体的解析细节（比如最大邻居数和截断半径，高斯距离的参数）

  2.`__len__()`方法

  让CIFData的实例可以通过`len()`函数返回数据集大小。

  3.`__getitem__()`方法

  这里先用了`@functools.lru_cache(maxsize=None)`修饰器来缓存数据，避免每次读入相同结构时，都要解析晶体。

- `DataLoader`中的`collate_fn`打包函数

  CGCNN代码中自定义了打包函数，名叫`collate_pool()`

  `DataLoader`实例在初始化同样只是规定了一些加载参数，并没有实际加载数据集。

  只有当遍历`DataLoader`实例时，才开始加载、打包、返回数据。(`for`循环抽打`DataLoader`，`DataLoader`抽打`Dataset`)

  但是遍历`DataLoader`并不会得到一个个的数据点，而是得到一个个batch的数据包，这是由`collate_fn`打包函数完成的，每当经历了batch_size个数据点，打包函数就会把他们合并。

  这里是打包函数返回的内容

  ```python
      return (torch.cat(batch_atom_fea, dim=0),
              torch.cat(batch_nbr_fea, dim=0),
              torch.cat(batch_nbr_fea_idx, dim=0),
              crystal_atom_idx),\
          torch.stack(batch_target, dim=0),\
          batch_cif_ids
  ```

  以及`for`循环的例子， `for i, (input, target, _) in enumerate(train_loader):`

  得到的`input`其实是第一个返回，即一个元组，包含了经过拼接后的同一批次内的原子特征、邻居特征、以及邻居特征索引和全局原子索引。



## 写神经网络

```python
model_input = 'resnet152'
# Build and train network
if model_input == 'vgg16':
    # Build and train network
    model = models.vgg16(pretrained=True)
elif model_input == 'vgg19':
    # Build and train network
    model = models.vgg19(pretrained=True)
elif model_input == 'resnet152':
    model = models.resnet152(pretrained=True)

# Freeze training for all layers
for param in model.parameters():
    param.requires_grad = False


# # Newly created modules have require_grad=True by default
if 'vgg' in model_input:
    num_features = model.classifier[-1].in_features
    model.classifier[6] = nn.Sequential(
                          nn.Linear(num_features, 512), 
                          nn.ReLU(), 
                          nn.Linear(512, len(cat_to_name)),
                          nn.LogSoftmax(dim=1))
elif 'resnet' in model_input:
    num_features = model.fc.in_features
    model.fc = nn.Sequential(
                            nn.Linear(num_features, 512), 
                              nn.ReLU(), 
                              nn.BatchNorm1d(512),
                              nn.Dropout(0.4),
                              nn.Linear(512, len(cat_to_name)),
                              nn.LogSoftmax(dim=1))

print(model.__class__.__name__)
# check to see that your last layer produces the expected number of outputs
# print(model.classifier[-1].out_features)
```

根据参数选择预训练模型并冻结模型参数，然后修改输出层。

对于`vgg`模型，

```python
        self.classifier = nn.Sequential(
            nn.Linear(512 * 7 * 7, 4096),
            nn.ReLU(True),
            nn.Dropout(p=dropout),
            nn.Linear(4096, 4096),
            nn.ReLU(True),
            nn.Dropout(p=dropout),
            nn.Linear(4096, num_classes),
        )
```

他的分类器本来包含了7层，但是作者觉得最后一层的4096个特征太多了，就改写了一下，减少特征为512个。



## `train()`函数

让我们看看一个合格的训练函数都应该做哪些工作。

```python
def train(train_loader, model, criterion, optimizer, epoch, normalizer):
    #每一代epoch都调用一次train函数
    batch_time = AverageMeter()
    data_time = AverageMeter()
    losses = AverageMeter()
    if args.task == 'regression':
        mae_errors = AverageMeter()
    else:
        accuracies = AverageMeter()
        precisions = AverageMeter()
        recalls = AverageMeter()
        fscores = AverageMeter()
        auc_scores = AverageMeter()

    # switch to train mode
    model.train()

    end = time.time()
    #一个训练函数应该要遍历整个训练集
    for i, (input, target, _) in enumerate(train_loader):
        #遍历train_loader，其实是一个个的批次
        # measure data loading time
        data_time.update(time.time() - end)

        if args.cuda:
            input_var = (Variable(input[0].cuda(non_blocking=True)),
                         Variable(input[1].cuda(non_blocking=True)),
                         input[2].cuda(non_blocking=True),
                         [crys_idx.cuda(non_blocking=True) for crys_idx in input[3]])
        else:
            input_var = (Variable(input[0]),
                         Variable(input[1]),
                         input[2],
                         input[3])
        # normalize target
        if args.task == 'regression':
            target_normed = normalizer.norm(target)
        else:
            target_normed = target.view(-1).long()
        if args.cuda:
            target_var = Variable(target_normed.cuda(non_blocking=True))
        else:
            target_var = Variable(target_normed)

        # compute output
        #把这一批次的数据带入模型
        output = model(*input_var)
        #由模型得到的结果和目标值得到损失函数
        loss = criterion(output, target_var)

        # measure accuracy and record loss
        if args.task == 'regression':
            mae_error = mae(normalizer.denorm(output.data.cpu()), target)
            losses.update(loss.data.cpu(), target.size(0))
            mae_errors.update(mae_error, target.size(0))
        else:
            accuracy, precision, recall, fscore, auc_score = \
                class_eval(output.data.cpu(), target)
            losses.update(loss.data.cpu().item(), target.size(0))
            accuracies.update(accuracy, target.size(0))
            precisions.update(precision, target.size(0))
            recalls.update(recall, target.size(0))
            fscores.update(fscore, target.size(0))
            auc_scores.update(auc_score, target.size(0))

        # compute gradient and do SGD step
        optimizer.zero_grad()
        loss.backward()                           #后向传播、计算梯度
        optimizer.step()                          #传递梯度、调整参数

        # measure elapsed time
        batch_time.update(time.time() - end)
        end = time.time()

        if i % args.print_freq == 0:
            if args.task == 'regression':
                print('Epoch: [{0}][{1}/{2}]\t'
                      'Time {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                      'Data {data_time.val:.3f} ({data_time.avg:.3f})\t'
                      'Loss {loss.val:.4f} ({loss.avg:.4f})\t'
                      'MAE {mae_errors.val:.3f} ({mae_errors.avg:.3f})'.format(
                    epoch, i, len(train_loader), batch_time=batch_time,
                    data_time=data_time, loss=losses, mae_errors=mae_errors)
                )
            else:
                print('Epoch: [{0}][{1}/{2}]\t'
                      'Time {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                      'Data {data_time.val:.3f} ({data_time.avg:.3f})\t'
                      'Loss {loss.val:.4f} ({loss.avg:.4f})\t'
                      'Accu {accu.val:.3f} ({accu.avg:.3f})\t'
                      'Precision {prec.val:.3f} ({prec.avg:.3f})\t'
                      'Recall {recall.val:.3f} ({recall.avg:.3f})\t'
                      'F1 {f1.val:.3f} ({f1.avg:.3f})\t'
                      'AUC {auc.val:.3f} ({auc.avg:.3f})'.format(
                    epoch, i, len(train_loader), batch_time=batch_time,
                    data_time=data_time, loss=losses, accu=accuracies,
                    prec=precisions, recall=recalls, f1=fscores,
                    auc=auc_scores)
                )
```



- `def train(train_loader, model, criterion, optimizer, epoch, normalizer):`

这里接受的参数有：训练集的DataLoader、模型、损失函数、优化器、epoch、正则化。

- `model.train()`

  首先，把模型的训练模式打开，这样网格中一些特殊的函数比如Dropout才会被开启。

- `for`循环遍历训练集的`DataLoader`

  把得到的批数据，经过简单的处理，比如把张量上传到CUDA。

  把批数据流入模型，再后向传播。
