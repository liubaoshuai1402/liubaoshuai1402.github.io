---
title: 稍微深入一些Python中的class
data: 2025-05-13
series: ["python"]
categories: [代码]
---

# 稍微深入一些Python中的类（class）

## 前言

**类（class）**在python代码中几乎无处不在，但在近日的学习中发现，我对它真是了解甚少，甚至基础结构都不能熟稔于心，故开此笔记认真学习。和我一起重新认识一下它吧。

## 一个简单的类

```python
#code 1
class Dog:
    # 类属性
    species = "Dog"
    # 初始化方法
    def __init__(self, name, age):
        self.name = name
        self.age = age
    # 实例方法
    def bark(self):
        return print("旺旺")

print(mydog.species)
print(mydog.name,mydog.age)
mydog.bark()

#输出为：
#Dog
#doudou 2
#旺旺
```

**类**是一种对数据进行计算操作的蓝图，离不开**属性**和**方法**。

code 1中，我们定义了一个 Dog 类，并对它进行了实例化，生成了一个对象。

这个类的结构，很简单。

首先是放在第一部分的**类属性**。**类属性**是**直接在类中定义变量**。所有通过这个类生成的对象都具有这些属性。



然后是放在第二部分的诸多方法，其实就是一个个的函数。

```python
    def __init__(self, name, age):
        self.name = name
        self.age = age
```

- `_init_`方法叫**初始化方法**，是**魔法方法**的一种。让实例初始化时就具有`name`和`age`属性和相应的值。具体表现为`mydog = Dog("doudou",2)`，这个类在初始化时就**需要两个参数**才能转变为实例。

- python中所有的实例方法，包括`__init__`方法，第一个参数都必须是`self`，用于区别**普通函数**和**方法**。

- **默认值**，在定义方法的时候，可以传入默认值，这可以保证在不传入参数时，也能生成一个默认实例。

  ```python
      def __init__(self, name="小卡拉米", age=5):
          self.name = name
          self.age = age
  ```

  这样的`__init__`方法就保证了，即使狗主人忘记了填写信息，`mydog = Dog()`也能正常工作，但默认的狗狗是**小卡拉米**，还是应该记得为自己的狗狗正确填写信息呦。（对于其他实例方法也适用）

在一般的实例方法中，也可以设置参数，定义属性。我们看一段新的代码:

```python
#code 2-1
class phone:
    owner = "xiyangyang"
    def __init__(self,type="huawei"):
        self.type = type
    def call(self,number=10086):
        self.number = number
        print("the number is {}".format(self.number))
print(phone.owner)
myphone = phone("apple")
myphone.call(110)
print("what number is called? ",myphone.number)
```

这里在`phone`类的`call`方法下，定义了一个**实例属性**，并需要一个输入`number`指定值。

但实例方法的参数不一定都是实例属性，且看以下代码：

```python
#code 2-2
class phone:
    def __init__(self,type="huawei"):
        self.type = type
    def call(self,number=10086):
        print("the number is {}".format(number))
myphone = phone("apple")
myphone.call(110)
```

这样定义`call`方法，在使用的时候也需要传入一个`number`，但是这个`number`不是这个实例的属性，无法随时查看。

这就是说，对于类中的方法，目前涉及的参数有三类，一个是默认参数`self`，一个是**实例属性参数**，一个是**普通参数**。

1. 截至目前，我们已经了解了一个简单**类的结构**。即：声明类属性，声明实例方法（其中涉及到添加实例属性）。
2. 由类生成实例时的参数接口，由魔法方法`__init__`定义，而一般实例方法的接口，在调用方法的时候才会出现。这也就意味着，`__init__`方法中定义的实例属性，是伴生的，只要生成实例就存在。而随一般实例方法定义的实例属性，则需要首次调用后才存在。

## 魔法方法

在python中，所有被双下划线包围的方法，统称为魔法方法。太多了，我这里只记录我遇到的，不定时更新。

### 1.`__init__`方法

这个魔法方法用于类的初始化。规定了实例化**类**时接受的参数。

### 2.`__len__`方法

```python
#code 3-1
class menu:
    restaurant="hepingfandian"
    def __init__(self,foods):
        self.foods = foods
    def __len__(self):
        count = 0
        for food in foods:
            count = count + 1
        return count
mymenu = menu(["宫保鸡丁","鱼香肉丝","米饭"])
print(len(mymenu))
#输出：3
```

这个魔法方法让实例可以被用`len()`函数统计长度，**核心是要返回一个整数**，至于这个整数是如何计算得到的，这部分内容就由自己定义了。

### 3.`__getitem__`方法

```python
#code 3-2
class menu:
    restaurant="hepingfandian"
    def __init__(self,foods):
        self.foods = foods
    def __getitem__(self,key):
        return 10
mymenu = menu(["宫保鸡丁","鱼香肉丝","米饭"])
print(mymenu["霸王餐"])
#输出：10
```

这个魔法方法可以让实例像字典或列表一样，实现键值对和索引切片的功能。具体来说，就是**根据条件返回不同的值**。不然就像上述代码一样，客人要吃霸王餐，对应的值却是10。

```python
#code 3-3
class menu:
    restaurant="hepingfandian"
    def __init__(self,foods):
        self.foods = foods
    def __getitem__(self,key):
        if not isinstance(key,str):
            raise TypeError("Attribute key must be a string")
        
        if key == "霸王餐":
            return "你吃牛魔"

        for food in self.foods:
            if food == key:
                return self.foods[key]
        raise KeyError(f"Food '{key}' does not exist in menu.")
        
mymenu = menu({"宫保鸡丁":25,"鱼香肉丝":10,"米饭":2})
print(mymenu["霸王餐"])
print(mymenu["米饭"])
#输出：你吃牛魔
#输出：2
```

`__getitem__`方法可以为实例提供独特的接口，就是`[参数]`。这里的参数可以是字符串、数字或者切片，使用起来像字典或者列表。而一般方法的使用则要`.onemethod(参数1，参数2，...)`。

```python
#code 3-4
class top3:
    ip = "harbin"
    def __init__(self,school):
        self.school = school
    def __getitem__(self,index):
        if isinstance(index,int):
            return self.school[index]
        if isinstance(index,slice):
            return self.school[index]
Top3InMyheart = top3(["qinghua","beida","hagongda"])

print("who is top3?",Top3InMyheart[0:3])
print(Top3InMyheart.ip)
```

这段代码定义了一个名叫`top3`类，实现了索引和切片功能，同时可以这个类的实例拥有列表不同的功能，即这个类的实例有一个类属性：`ip`。让我们可以得知这个排名来自`harbin`。

### 4.`__setitem__`方法

```python
#code 3-5
class games:
    def __init__(self):
        self.games={}
    def __setitem__(self,key,value):
        self.games[key]=value
    def __getitem__(self,key):
        return self.games[key]

MyLoveGames = games()
MyLoveGames[1] = "原神"
print(MyLoveGames[1])
#输出：原神
```

这个魔法方法让类的实例可以像字典或者列表一样，添加新的键值对或者添加新的索引和值。

### 5.`__delitem__`方法

```python
#code 3-6
class games:
    def __init__(self):
        self.games={}
    def __setitem__(self,key,value):
        self.games[key]=value
    def __getitem__(self,key):
        return self.games[key]
    def __delitem__(self,key):
        del self.games[key]

MyLoveGames = games()
MyLoveGames[1] = "原神"
MyLoveGames[1] = "DOTA2"
print(MyLoveGames[1])
del MyLoveGames[1]
#输出：DOTA2
```

这个魔法方法可以让类的实例使用`del`关键字来删除键或索引。

### 6.`__iter__`方法

```python
#code 3-7
class games:
    def __init__(self):
        self.games={}
    def __setitem__(self,key,value):
        self.games[key]=value
    def __getitem__(self,key):
        return self.games[key]
    def __delitem__(self,key):
        del self.games[key]
    def __iter__(self):
        for game in self.games:
            yield game
            
            
MyLoveGame = games()
MyLoveGame[0] = "原神"
MyLoveGame[3] = "DOTA2"
MyLoveGame[2] = "明日方舟"
del MyLoveGame[2]
for game in MyLoveGame:
    print(game)
#输出：0
#输出：3
```

这个魔法方法让实例变为可迭代对象，可以进行`for`循环。也可以用`iter()`函数生成迭代器。这里涉及到一个头疼的`yield`关键字，又涉及到**生成器**和**迭代器**。还是另开一文来记录吧。

简单来讲，当一个可迭代对象遇到`for`循环事件时，会自动转到它的`__iter__`方法，又因为这里的`__iter__`方法含有`yield`关键字，所以不会立即执行，而是得到了一个**生成器**对象。

紧接着，`for`循环又根据这个生成器对象的`__next__`方法（注意，这里是**生成器对象的**`__next__`方法，而非是，我们这个自定义类的`__next__`方法，其实，我们这里也没有写`__next__`方法）把得到的值赋值给`for game in MyLoveGame:`中的`game`，然后打印。

在进行第二次循环的时候，同样地，再次进入实例的`__iter__`方法，得到同一个**生成器**对象。并对这个生成器对象再此施加`__next__`方法，并把得到的值赋值给`for game in MyLoveGame:`中的`game`，然后打印。

- 谁是`for game in MyLoveGame:`中的`game`值？

`__iter__`方法本质上定义了一个**生成器**对象，而非函数。`for`循环会不断调用这个生成器对象`__next__`方法并把得到的值传递给`game`。

### 7.`__repr__`方法

```python
#code 3-8
import collections.abc
class games(collections.abc.MutableMapping):
    def __init__(self):
        self.games={}
    def __len__(self):
        return 0
    def __setitem__(self,key,value):
        self.games[key]=value
    def __getitem__(self,key):
        return self.games[key]
    def __delitem__(self,key):
        del self.games[key]
    def __iter__(self):
        for game in self.games:
            yield game
    def __repr__(self):
        return repr(dict(self))
            
MyLoveGame = games()
MyLoveGame[0] = "原神"
MyLoveGame[1] = "DOTA2"
MyLoveGame[2] = "明日方舟"

print(MyLoveGame)
#输出：{0: '原神', 1: 'DOTA2', 2: '明日方舟'}
```

这里让自定义的games继承了一个抽象基类，这样`return repr(dict(self))`才不报错。

这个魔法方法让实例可以直接使用实例时，返回一个官方字符串。

## 继承与抽象基类

- ### 如何继承父类

- ```python
  #code 3-9
  class student:
      def __init__(self,gender=1,age=18):
          self.gender = gender
          self.age = age
          
  class xiaoxiaobuaigugujiao(student):
      def interests(self,data):
          self.interests = data
          print(self.interests)
  me = xiaoxiaobuaigugujiao()
  print(me.gender)
  print(me.age)
  me.interests("games")
      
      
  ```

  在定义子类的时候，在子类的名字后加上`(父类名)`即可。

- ### super()函数

- ```python
  #code 3-9
  class student:
      def __init__(self,gender=1,age=18):
          self.gender = gender
          self.age = age
          
  class xiaoxiaobuaigugujiao(student):
      #子类中，当覆盖父类中的同名方法后，与父类同名方法的默认参数也要重新写一篇，不会继承
      def __init__(self,gender=1,age=18,height=170):
          #调用父类的构造方法
          super().__init__(gender,age)
          #额外补充
          self.height = height
      def interests(self,data):
          self.interests = data
          print(self.interests)
  me = xiaoxiaobuaigugujiao()
  print(me.gender)
  print(me.age)
  print(me.height)
  me.interests("games")
  ```

  `super()`函数通常出现在子类定义方法的时候，用于调用父类的构造方法，并可以加以补充，这里我们在子类初始化时，就调用了父类`student`的构造方法，并额外补充了`height`。