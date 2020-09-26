# 实验一：Python快速入门
参考：[Python基础教程(crossin全60课)](https://github.com/dai0992/Pattern-Recognition-and-Prediction/blob/master/Python基础教程(crossin全60课).pdf)

## 1. 安装并运行Python
> 1）下载[Python](https://www.python.org/downloads/)，按步骤安装，在"Advanced Options"时，勾选"Add Python to environment variables"。<br>
> 2）验证是否安装成功：打开"命令提示符"（桌面按快捷键"Win+r"，输入cmd，回车），在命令行里输入python，如果看到python的版本信息，说明安装成功。<br>
> 3）第一声啼哭：在命令行输入 <br>
```python
print("Hello World!")
```

## 2. 数据结构
### 2.1 基本数据类型
```python
name = 'Tomas' # 字符串变量（单引号、双引号都可）
myInt = 666 # 整数型变量
myFloat = 1.618 # 浮点型变量
myBool = True # 逻辑型变量
```

### 2.2 变量命名规则
> 1) 第一个字符必须是字母或者下划线
> 2) 剩下的部分可以是字母、下划线或数字
> 3) 变量名称对大小写敏感，比如myname和myName不是同一个变量

**几个有效的变量名**
> a, _abc, abc_12, a1b2_c3

**几个无效的变量名**
> 2dog, my-name

### 2.3 列表(list)
```python
myList1 = [4, 2, 3, 2, 5, 1] # 用一对`中括号`创建一个列表型变量，每个元素都是整数
myList2 = ['meat', 'egg', 'fish', 'milk'] # 每个元素都是字符串
myList3 = [365, 'everyday', 0.618, True] # 每个元素也可以是不同的基本数据类型
myList3[1] # 会输出'everyday', python中的元素计数从0开始（不同于R、MATLAB，但和其他大多数语言相同）
myList3[4] # 会报错，4代表'myList3'中的第5个元素，不存在！
myList3[2] = 'happy' # 将其中的第3个元素0.618修改为字符串'happy'
myList3.append(666.666) # 在myList3的尾巴上增加1个元素，这里使用了list的append方法/函数
myList3[4] # 不会报错，会输出 666.666
del(myList3[1]) # 将第2个元素'everyday'删除
# 此时myList3的内容为[365, 'happy', True, 666.666]
myList3[-1] # 表示索引/选取最后1个元素（特别要注意和正向索引的区别，正向从0、负向从1计数）
myList3[1:3] # 输出['happy', True], 这里使用了切片索引（特别特别注意：冒号前的切片包括、冒号后的切片不包括，这是初学python最坑的地方）
myList3[:3] # 输出[365, 'happy', True]，等价于myList3[0:3]
myList3[1:] # 输出['happy', True, 666.666]，等价于myList3[1:4]（虽然索引4代表不存在的第5个元素，但冒号后切片不包括，所以能取到最后1个元素。很奇葩！）
myList3[1:-1] # 输出['happy', True]，等价于myList3[1:3]
```

### 2.4 字符串
在字符串中表示单引号
> "What's your name?"

在字符串中表示双引号
> 'You are a "BAD" man'

用\\'表示单引号，用\\"表示双引号（ \ 被称为转义字符，\n表示换行，\ 还可用在代码中换行）
> 'I\\'m a \\"good\\" student'

字符串拼接
```python
str1 = 'good'
str2 = 'student'
str3 = str1 + ' ' + str2
str4 = 'I\'m a ' + str1 + ' ' + str2
str_and_num = str1 + 666 # 字符串和数字相加会报错
str_and_num = str1 + str(666) # 用str函数把数字转换为字符串，不会报错
str_and_num = 'good %d' % 666 # 用%对字符串进行格式化，另外有%f, %.2f, %s等
"%s's score is %d" % ('Mike', 90) # 同时用多个%对多个变量格式化
# 假如我们有一个字符串列表
str_list = ['apple', 'pear', 'orange']
'-'.join(str_list) # 输出'apple-pear-orange'，以短横线将各字符串元素连接
''.join(str_list) # 输出'applepearorange', 连接符可以是空串
```

字符串分割
```python
sentence = 'I am a sentence'
sentence.split() # split()函数会默认按空格分割字符串，每个子串组成一个list
section = 'Come on. Let\'s go. Go go go.'
section.split('.') # 指定'.'为分隔符
```

字符串的索引和切片
```python
word = 'helloworld'
word[0] # 输出'h'
word[-2] # 输出'l'
word[0] = 'w' # 会报错，字符串不允许通过索引修改其中字符。如果想修改，只有先转换成列表（用list函数），修改后再用空串将列表中每个字符连接成字符串（用''.join(yourList)）
word[:5] # 输出'hello', 切片规则和list相同
'^_^'.join(word) # 输出''h^_^e^_^l^_^l^_^o^_^w^_^o^_^r^_^l^_^d''，（我只是耍个帅 ^_^）
```

### 2.5 字典（dictionary, 类似Perl中的哈希hash）
字典是键/值对(key:value)的集合，每个键值对之间用逗号分隔，整个字典包括在`花括号`中
> 1) 键必须是唯一的
> 2) 键只能是基本数据类型，如整数、浮点数、字符串、逻辑值，list不能作为键
> 3) 值没有要求
> 4) 键值对没有顺序，因此无法用索引访问字典中的内容，需要用键来访问
```python
score = {'萧峰':95, '段誉':97, '虚竹':90, True:'Good brother', 100:'Perfect score'} # 用大括号建立一个字典
score['虚竹'] # 输出键对应的值：90
score[True] # 输出键对应的值：'Good brother'
score['虚竹'] = 99 # 修改“虚竹”的得分（后期有无崖子和天山童姥的百余年功力）
score['慕容复'] = 88 # 增加“慕容复”的得分
del(score[100]) # 删除键为100的键值对（注意这里100不是索引，而是score的一个键）
```

### 2.6 元组(tuple)
和列表list类似，但是元组中的元素不可更改，元组用`小括号`创建。元组和list同样有索引、切片、遍历等操作
```python
position = (147, 258) # 创建一个只包含数字的元组
weather = ('sunny', 'cloudy', 'rainy') # 创建一个只包含字符串的元组
weather_id = (1, 'sunny', 2, 'cloudy', 3, 'rainy') # 创建一个包含数字和字符串的元组
weather_id[:2] # 输出(1, 'sunny')
```

### 2.7 数据类型转换
```python
int('123') # 输出：123, 字符串转整数
float('6.6') # 输出：6.6, 字符串转小数
str(168) # 输出：'168', 数字转字符串
bool(0) # 输出：False, 数字转逻辑值
int('abc') # 会报错，字符串abc不可能转成数字，不符合常识

bool(-3) # 输出：True
bool(False) # 输出：False, 此时的False是python中的特殊关键字，代表0
bool('False') # 输出：True, 此时'False'只是个不为空的字符串
bool('') # 输出：False, 什么都没有
bool(' ') # 输出：True, 看上去空，实际有一个空格
```


## 3. 程序控制
### 3.1 逻辑判断
```python
x1 = 2
x2 = 8
x1 < 3 # True
x1 == x2 # False
x1 != x2 # True
# -------- and ---------
x1 < 10 and x2 < 10 # True
x1 < 10 and x2 > 10 # False
x1 > 10 and x2 < 10 # False
x1 > 10 and x2 > 10 # False
# -------- or ---------
x1 < 10 or x2 < 10 # True
x1 < 10 or x2 > 10 # True
x1 > 10 or x2 < 10 # True
x1 > 10 or x2 > 10 # False
# -------- not ---------
not(x1<10) # False
not(x1>10) # True
```

### 3.2 判断语句
**格式**
> if 判断条件: <br>
>> 执行的内容1 <br>
>> 执行的内容2 <br>

**特别说明：判断条件后面的`冒号`不能少，if内部的语句需要有`统一的缩进`，一般用4个空格或按一次tab键，并且整个文件要统一，不能空格和tab混用**
```python
if x1 < 10:
  print('x1 is less than 10') # 命令行会自动输出3个点，需要按一次tab键，然后再输入print命令，按2次回车，输出结果

# and判断
if x1 < 10 and x2 < 10:
  print('x1 and x2 are both less than 10')

# if...else...语句
if x1 < x2:
  print('x1 is less than x2')
else:
  print('x1 is more than x2')

# if...elif...else语句
if x1 > x2: # 可以通过设置x1, x2的大小，来得到不同的输出
  print('x1 is more than x2')
elif x1 > 10:
  print('x1 is less than x2, but x1 is more than 10')
else:
  print('x1 is less than x2, and x1 is less than 10')

# if的嵌套
if x1 < 10:
  if x2 > 10:
    print('x1 is less than 10, but x2 is more than 10')
  else:
    print('x1 and x2 are both less than 10')
```

### 3.3 循环语句
**while循环语句格式(同样注意判断条件后面的冒号不能丢)**
> while 判断条件: <br>
>> 执行的内容1 <br>
>> 执行的内容2 <br>
```python
iter_m = 0
x1 = 2
while x1 < 20:
  x1 = x1*2
  iter_m = iter_m+1
  print('Iteration %d: x1 is %d' % (iter_m, x1))
```

**for循环语句格式(同样注意循环范围后面的冒号不能丢)**
> for ... in 循环范围: <br>
>> 执行的内容1 <br>
>> 执行的内容2 <br>
```python
x1 = 2
for i in range(5,10): # 这里使用了range函数产生5到10之间的整数，但不包括10
  x1 = x1*i
  print('i in range(5,10) is %d: x1 is %d' % (i, x1))
```

**循环的嵌套**
```python
for i in range(0,5):
  for j in range(0,i+1): # 好好体会这里j循环的范围，是随着外层i的取值变化的
    print('*', end='')
  print()
```

**break: 满足条件则结束`本层`循环**
```python
x1 = 2
iter_m = 0
while 1: # 最粗暴的判断条件，1代表无限循环，即死循环
  x1 = x1+1
  iter_m = iter_m+1
  print('Iteration %d: x1 is %d' % (iter_m, x1))
  if x1 > 10:
    print('x1 is more than 10 and the while loop should be break!')
    break # 如果没有break语句，循环将无限进行下去
```

**continue: 满足条件则结束`本次`循环**
```python
x1 = 2
iter_m = 0
while x1 < 10:
  x1 = x1+1
  iter_m = iter_m+1
  if x1 % 2 == 0: # 如果x1能被2整除，则打印提示，并结束本次循环
    print('x1 is a even number!')
    continue
  print('Iteration %d: x1 is a even number %d' % (iter_m, x1))
```

## 4. 读写文件
比如有一个文件：data.txt，内容如下
> Hi man! <br>
> I am a file. <br>
> Try read my mind and print it onto screen! <br>

一次性读所有内容，并输出到屏幕
```python
f = open('data.txt') # 使用open函数打开文件，并返回文件句柄给变量f
data = f.readlines() # 使用read函数一次性读取所有内容（当文件很大时，慎用！！）
print(data) # 一次性打印所有内容到屏幕
f.close() # 关闭文件（不论多复杂的程序，一旦打开过文件，记得最后一定要关闭文件）

读一行处理一行，并将处理结果写到一个新文件中
```python
f = open('data.txt')
line = f.readline()
iter_num = 1
while line != '':
```



## 5. 函数


## 6. 正则表达式


## 7. 模块


