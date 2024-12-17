# 构造PDA,用秩的概念来优化代码
from goto import with_goto
import copy
import sys
import numpy as np
import time
# 计算高斯系数,k维向量空间中所有m维子空间的个数


def affineCount(k, m, q):
    end = k-m
    sumk = 1
    for i in range(k, end, -1):
        sumk = sumk * (q**i - 1)

    summ = 1
    for i in range(m, 0, -1):
        summ = summ * (q**i - 1)

    return sumk // summ

# 计算有限域q上n维向量空间 k个向量组成的线性无关组的个数(看似不同的组，可以张成同一个空间)


def vectorBase(n, k, q):
    cnt = 1
    for i in range(0, k):
        cnt = cnt * (q**n - q**i)
    return cnt


@with_goto
def inputData():
    label .begin
    print("分别输入维度k ,素数幂q, m(<=k) , n , t (m+n+t+1<=k), :")
    list1 = list(map(int, input().strip().split()))
    if(len(list1) != 5):
        print("输入数字不是5个")
        goto .begin
    return list1


@with_goto
def inputData2():
    label .begin
    print("输入总维度n，素数幂q,子空间维度k:")
    lis = list(map(int, input().strip().split()))
    if (len(lis) != 3):
        print("输入数字不是3个")
        goto .begin
    else:
        return lis

# 一组基向量张成空间
# n 向量维度，k 基向量个数,q 素数幂,layer 递归层数(k层),base1 基向量,coefficient 基的系数,totalSpace 张成的空间


def bVExtendSpace(n, k, q, layer, base1, coefficient, totalSpace):
    if(layer >= k):
        lis = []  # 可能需要lis = list()
        for y in range(n):
            val = 0
            for x in range(k):
                val = val + coefficient[x] * base1[x][y]
            val = val % q  # 结果求模
            lis.append(val)
        if lis not in totalSpace:
            totalSpace.append(lis)
    else:
        for i in range(q):
            coefficient[layer] = i
            bVExtendSpace(n, k, q, layer+1, base1, coefficient, totalSpace)

# 求n维空间的所有向量,q素数幂


def allVector_n(n, q):
    totalSpace = []
    # 求n维空间的基向量
    base1 = []
    for i in range(n):
        item = [0] * n
        item[i] = 1
        base1.append(item)
    totalSpace = []
    coefficient = [0] * len(item)  # 基的系数个数
    bVExtendSpace(len(base1[0]), len(base1), q, 0,
                  base1, coefficient, totalSpace)
    return totalSpace


# 用全部空间向量totalSpace 求出k维空间所有的基base2(包括可生成同一空间的写法不同的基),不包括零向量
def subSpaceBase(k, layer, pos, totalSpace, tmp, base2):
    if(k <= 0):
        base2.append(copy.deepcopy(tmp))  # 重要坑点，tmp改变，直接存tmp base2的内容也跟着改变
    else:
        for i in range(pos+1, len(totalSpace)):
            tmp.append(totalSpace[i])
            if(np.linalg.matrix_rank(tmp) == layer):          # 计算向量组的秩，必须线性无关
                subSpaceBase(k-1, layer+1, i, totalSpace, tmp, base2)
            tmp.pop()

# 求所有k维子空间,q 素数幂,totalSpace n维空间


def allSubSpace(k, q, totalSpace):
    # 求出t维子空间所有的基
    t = k
    tmp, base2 = [], []
    subSpaceBase(t, 1, 0, totalSpace, tmp, base2)
    print("{0}维空间所有的基个数：{1}".format(k, len(base2)))

    '''
    # 从base2中筛选生成的空间大部分不重复的基
    effectBase = []
    flag = [0] * len(base2)
    for i in range(0,len(base2)):
        if flag[i]==0:
            effectBase.append(base2[i])
            for j in range(i+1,len(base2)):
                AB = copy.deepcopy(base2[i])         # 二层引用赋值，牵一发而动全身
                AB.extend(base2[j])
                if np.linalg.matrix_rank(base2[i]) == np.linalg.matrix_rank(AB):
                    flag[j] = 1

    base2 = effectBase
    print("{0}维有效的基{1}".format(t,len(base2)))
    '''
    # 求所有t维空间
    allSubSpace_k = []
    subSpaceBase_k = []
    for item in base2:
        subSpace = []
        coefficient = [0]*t  # 基的系数
        # n 向量维度，k 基向量个数,q 素数幂,layer 递归层数(从0开始k层),base1 基向量,coefficient 基的系数,totalSpace 张成的空间
        bVExtendSpace(len(item[0]), len(item), q, 0,
                      item, coefficient, subSpace)
        subSpace = sorted(subSpace)             # 对子空间内的向量排序
        if len(subSpace) == q**k and subSpace not in allSubSpace_k:       # k维空间的元素有q**k个
            allSubSpace_k.append(subSpace)
            subSpaceBase_k.append(item)
    return allSubSpace_k, subSpaceBase_k


# 输入数据,k q m n t
list1 = inputData()
# sys.stdout = open('F:\desktop\毕设\Python程序\out3.txt','w')
startTime = time.clock()
# 求n维总空间
totalSpace = allVector_n(list1[0], list1[1])
print("{0}维空间".format(list1[0]))

# 求t维子空间及其生成生成空间对应的基
t = list1[4]
allSubSpace_t, subSpaceBase_t = allSubSpace(t, list1[1], totalSpace)

# n维向量空间 基向量个数为t的基的个数(看似不同的组，可以张成同一个空间)
print("{0}维向量空间 {1}个向量组成的线性无关组的公式计算个数:{2}(包括重复)".format(
    list1[0], t, vectorBase(list1[0], t, list1[1])))

# 计算高斯系数,k维向量空间中所有t维子空间的个数
print("{0}维子空间个数 {1},公式计算为{2}".format(
    t, len(allSubSpace_t), affineCount(list1[0], t, list1[1])))

# 输入固定的t-1维子空间L
L = [[1, 0, 0, 0, 0]]
lowwerL = 1
V = []
VBase = []
# 包含L的所有t维子空间
for item in allSubSpace_t:
    flag = 0
    for vect in L:
        if vect not in item:
            flag = 1
            break
    if flag == 0:
        pos = allSubSpace_t.index(item)
        V.append(item)
        VBase.append(subSpaceBase_t[pos])
# 输出包含L的t维子空间
print("包含L的{0}维子空间个数{1} ---,公式计算为{2}".format(t, len(V),affineCount(list1[0]-lowwerL, t-lowwerL, list1[1])))
for item in V:
    print(item)
# for i in range(len(V)):
#     print("---", VBase[i])

# 包含L的n+t维子空间
n_t = list1[3] + list1[4]
allSubSpace_nt, subSpaceBase_nt = allSubSpace(n_t, list1[1], totalSpace)


# n维向量空间 基向量个数为t的基的个数(看似不同(不重复)的组，可以张成同一个空间)
print("{0}维向量空间 {1}个向量组成的线性无关组的个数公式计算为:{2}(包括重复)".format(
    list1[0], n_t, vectorBase(list1[0], n_t, list1[1])))

# 计算高斯系数,k维向量空间中所有t维子空间的个数
print("{0}维子空间个数 {1} 公式计算为{2}".format(
    n_t, len(allSubSpace_nt), affineCount(list1[0], n_t, list1[1])))


R = []
RBase = []
for item in allSubSpace_nt:
    flag = 0
    for vect in L:
        if vect not in item:
            flag = 1
            break
    if flag == 0:
        pos = allSubSpace_nt.index(item)
        RBase.append(subSpaceBase_nt[pos])
        R.append(item)

# 输出包含L的n+t维子空间
print("包含L的{0}维子空间个数{1}---,公式计算为{2}".format(n_t, len(R),
      affineCount(list1[0]-lowwerL, n_t-lowwerL, list1[1])))
for item in R:
    print(item)
# for i in range(len(R)):
#     print("--", RBase[i])

# 包含L的m+t维子空间
m_t = list1[2] + list1[4]
allSubSpace_mt, subSpaceBase_mt = allSubSpace(m_t, list1[1], totalSpace)

# n维向量空间 基向量个数为t的基的个数(看似不同的组，可以张成同一个空间)
print("{0}维向量空间 {1}个向量组成的线性无关组的个数:{2}(包括重复)".format(
    list1[0], m_t, vectorBase(list1[0], m_t, list1[1])))

# 计算高斯系数,k维向量空间中所有t维子空间的个数
print("{0}维子空间个数 {1} 公式计算为{2}".format(
    m_t, len(allSubSpace_mt), affineCount(list1[0], m_t, list1[1])))

S = []
SBase = []
for item in allSubSpace_mt:
    flag = 0
    for vect in L:
        if vect not in item:
            flag = 1
            break
    if flag == 0:
        pos = allSubSpace_mt.index(item)
        SBase.append(subSpaceBase_mt[pos])
        S.append(item)

# 输出包含L的m+t维子空间
print("包含L的{0}维子空间个数{1}---,公式计算为{2}".format(m_t, len(S),
      affineCount(list1[0]-lowwerL, m_t-lowwerL, list1[1])))
for item in S:
    print(item)
# for i in range(len(S)):
#     print("--", SBase[i])


# 包含L的m+n+t+1维子空间
m_n_t = list1[2] + list1[3] + list1[4] + 1
allSubSpace_mnt, subSpaceBase_mnt = allSubSpace(m_n_t, list1[1], totalSpace)

# n维向量空间 基向量个数为t的基的个数(看似不同的组，可以张成同一个空间)
print("{0}维向量空间 {1}个向量组成的线性无关组的个数:{2}(包括重复)".format(
    list1[0], m_n_t, vectorBase(list1[0], m_n_t, list1[1])))

# 计算高斯系数,k维向量空间中所有t维子空间的个数
print("{0}维子空间个数 {1} 公式计算为{2}".format(m_n_t, len(
    allSubSpace_mnt), affineCount(list1[0], m_n_t, list1[1])))

U = []
UBase = []
for item in allSubSpace_mnt:
    flag = 0
    for vect in L:
        if vect not in item:
            flag = 1
            break
    if flag == 0:
        pos = allSubSpace_mnt.index(item)
        UBase.append(subSpaceBase_mnt[pos])
        U.append(item)

# 输出包含L的m+n+t维子空间
print("包含L的{0}维子空间个数{1}---,公式计算为{2}".format(m_n_t, len(U),
      affineCount(list1[0]-lowwerL, m_n_t-lowwerL, list1[1])))
# for item in U:
#     print(item)
for i in range(len(U)):
    print("--", UBase[i])

endTime = time.clock()
print("运行时间:{0}".format(endTime-startTime))  # 暴力搜索用时40.7746059 s

# 生成X Y Z
# 从挑出的空间进行子空间求和,看求和的结果组成的空间是否是R的元素


def judgeAndAdd(VSpace, n, res, p, Vi):
    if(n < 0):
        if res not in Vi:
        #    Vi.append(res)   #Vi 与 lis 关联在一次了
            Vi.append(copy.deepcopy(res))   #只获取res的值，而非简单指向res
        return 
    else :
        lis = [0] * len(VSpace[0][0])
        for item in VSpace[n]:
            for i in range(len(item)):
                lis[i] = (res[i] + item[i]) % p     #python递归时注意全局变量牵一发而动全身
            judgeAndAdd(VSpace, n-1, lis, p, Vi)

# 从V的所有空间中取出n+1,(m+1),(m+n+2)个空间,n+1,(m+1),(m+n+2)个空间的和是R,(S),(U)的元素
def getSpace(X, V, ans, R, pos, VSpace, p):
    if(ans <= 0):
    #    print()
        Vi = []  #若干个向量的和
        res = [0]*len(VSpace[0][0])
        judgeAndAdd(VSpace, len(VSpace)-1, res, p, Vi)  # 计算这些子空间VSpace的和Vi
        Vi = sorted(Vi)
        global ans1
        ans1 = ans1 + 1
        if Vi in R:
        #    print(">> {0}".format(Vi))
            X.append(copy.deepcopy(sorted(VSpace)))   #对若干个空间进行排序 再复制值
        return 

    for i in range(pos, len(V)):
        VSpace.append(V[i])
    #    print("{0}-->".format(i),end=" ")
        getSpace(X, V, ans-1, R, i+1, VSpace,p)
        VSpace.pop()

# 计算列行标


def rowAndCol(k, q,n, t):
    # 计算阶乘
    denominator = 1
    for i in range(1, n+2):
        denominator *= i

    molecule = q ** (n*(n+1)//2)

    multiplier = 1
    for i in range(n+1):
        multiplier = multiplier * affineCount(k-t+1-i, 1, q)
    return multiplier * molecule // denominator


# 生成X list1 = {k,q,m,n,t}
X = []
VSpace = []
start1 = time.clock()
ans1 = 0
getSpace(X, V, list1[3]+1, R, 0, VSpace, list1[1])
end1 = time.clock()
print("生成X 运行时间:{0}".format(end1-start1))

print("X的元素个数 {0} 计数：{1}".format(len(X),ans1))
print("公式计算X的元素个数 {0}".format(rowAndCol(list1[0],list1[1],list1[3],list1[4])))

# 生成Y 
Y = []
VSpace = []
start2 = time.clock()
ans1 = 0
getSpace(Y,V,list1[2]+1,S,0,VSpace,list1[1])
end2 = time.clock()
print("生成Y 运行时间:{0}".format(end2-start2))
print("Y的元素个数 {0} 计数：{1}".format(len(Y),ans1))
print("公式计算Y的元素个数 {0}".format(rowAndCol(list1[0],list1[1],list1[2],list1[4])))

# 生成Z 
Z = []
VSpace = []
start3 = time.clock()
ans1 = 0
getSpace(Z,V,list1[2]+list1[3]+2,U,0,VSpace,list1[1])
end3 = time.clock()
print("生成Z 运行时间:{0}".format(end3-start3))
print("Z的元素个数 {0} 计数：{1}".format(len(Z),ans1))
print("猜测的公式计算Z的元素个数 {0}".format(rowAndCol(list1[0],list1[1],list1[2]+list1[3]+1,list1[4])))

# X Y 子空间合并后属于Z的元素
def judgeBelongToZ(i,j,X,Y,Z):
    contain = []
    for item in X:
        if item not in contain:
            contain.append(item)
    for item in Y:
        if item not in contain:
            contain.append(item)
    contain = sorted(contain)
    global realPDA
    if contain in Z:
        realPDA[i][j] = contain
        return True
    else :
        realPDA[i][j] = "*"
        return False
# 生成PDA,len(X)列,len(Y)行
pda = [[0 for i in range(len(X))] for i in range(len(Y))]
realPDA = [[0 for i in range(len(X))] for i in range(len(Y))]

for i in range(len(Y)):
    for j in range(len(X)):
        if judgeBelongToZ(i,j,X[j],Y[i],Z):
            pda[i][j] = [i,j]
        else :
            pda[i][j] = "*"

for i in range(len(pda)):
    for j in range(len(pda[0])):
        print(pda[i][j],end="  ")
    print()
print("X的空间：")
for i in range(len(X)):
    print("{0}: {1}".format(i,X[i]))

print("Y的空间")
for i in range(len(Y)):
    print("{0}: {1}".format(i,Y[i]))




