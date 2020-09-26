# 实验一：Python快速入门
参考：[Python基础教程(crossin全60课)](https://github.com/dai0992/Pattern-Recognition-and-Prediction/blob/master/Python基础教程(crossin全60课).pdf)

## 1. 安装Python
> 1）下载[Python](https://www.python.org/downloads/)，按步骤安装，在"Advanced Options"时，勾选"Add Python to environment variables"。<br>
> 2）验证是否安装成功：打开"命令提示符"（桌面按快捷键"Win+r"，输入cmd，回车），在命令行里输入python，如果看到python的版本信息，说明安装成功。<br>
> 3）第一声啼哭：在命令行输入<br>
>> ```python
>> print("Hello World!")
>> ```

## 2. 数据结构
### 2.1 基本数据类型
```python
name = 'Tomas' # 字符串变量（单引号、双引号都可）
myInt = 666 # 整数型变量
myFloat = 1.618 # 浮点型变量
myBool = True # 逻辑型变量
```

### 2.2 变量命名规则
> 1) 第一个字符必须是字母或者下划线 <br>
> 2) 剩下的部分可以是字母、下划线或数字 <br>
> 3) 变量名称对大小写敏感，比如myname和myName不是同一个变量 <br>

**几个有效的变量名**
> a, _abc, abc_12, a1b2_c3

**几个无效的变量名**
> 2dog, my-name

### 2.3 list
```python
myList1 = [4, 2, 3, 2, 5, 1] # 用一对中括号创建一个列表型变量，每个元素都是整数
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
```

字符串分割
```python
sentence = 'I am a sentence'
sentence.split() # split()函数会默认按空格分割字符串，每个子串组成一个list
section = 'Come on. Let\'s go. Go go go.'
section.split('.') # 指定'.'为分隔符
```



