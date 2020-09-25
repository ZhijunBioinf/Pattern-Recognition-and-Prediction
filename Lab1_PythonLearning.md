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

### 2.3 字符串
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
```



