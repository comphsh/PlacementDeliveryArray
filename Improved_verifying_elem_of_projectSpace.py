import copy,sys,time
import numpy as np
# 遍历 k维 模q的空间，顺便获取射影空间
def traverseSpace(k,q,n,ans,baseVector,vector):
    if n>=ans:
        global space1,space2
        space1.append(copy.deepcopy(vector))
        
        flag = 1
        for i in range(1,q):
            lis1 =[0 for j in range(k)]
            for j in range(k):
                lis1[j] = (vector[j] * i) % q
            if lis1 in space2:
                flag = 0
                break
        if flag==1 :
            print(vector)
            space2.append(copy.deepcopy(vector))
        return 

    lis = [0 for i in range(k)]
    for i in range(q):
        for j in range(k):
            lis[j] = (vector[j] + i * baseVector[n][j]) % q
        traverseSpace(k,q,n+1,ans,baseVector,lis)

# 模拟多重循环，只适用于k长向量生成k维空间
def multipleLoop(k,q,n,vector):
    if n>=k:
        global allProjectSpace
        allProjectSpace.append(copy.deepcopy(vector))
        return 
    for i in range(q):
        vector[n] = i
        multipleLoop(k,q,n+1,vector);

# 获取射影空间,即第一个非0位是1的所有向量
def getProjectiveSpace(k,q):
    for i in range(k-1,-1,-1):
        vector = [0 for j in range(k)]
        vector[i] = 1
        multipleLoop(k,q,i+1,vector)

# 模拟多重循环,使用于k长向量生成n维空间
def multipleLoopKN(vsize,k,q,n,vector):
    if n>=k:
        lis = [0 for i in range(vsize)]
        # 从全局变量visited 取得空间L的基向量为1的下标, vector的长度等于visited的长度
        global visited
        for i in range(len(visited)):
            lis[visited[i]] = vector[i]
        global L
        L.append(lis)
        return 
    for i in range(q):
        vector[n] = i
        multipleLoopKN(vsize,k,q,n+1,vector)
        

# 获取t-1维子空间L,vsize是向量的长度 ，k是t-1(空间L的维度)，q是素数幂
def getLSpace(vsize,k,q):
    for i in range(k-1,-1,-1):
        vector = [0 for j in range(k)]
        vector[i] = 1
        multipleLoopKN(vsize,k,q,i+1,vector)

def calSpaceElem(k,q):
    s = 0
    for i in range(0,k):
        s = s + q**i
    return s

# 模拟循环，计算a + k1c1 + k2c2 +... ;q是素数，size是c1 c2的个数，direct 是 k1c1 + k2c2 +...的值 
def calRSULoop(q,size,vector,direct,Vin):
    if size <=0:
        global L,allProjectSpace
        temp = [0 for j in range(len(direct))]
        for a in L:
            for j in range(len(a)):
                temp[j] = (a[j] + direct[j]) % q
            if temp in allProjectSpace and temp not in Vin:
                Vin.append(copy.deepcopy(temp))
        return
    for i in range(q):
        temp = []
        for j in range(len(vector[size-1])):
            temp.append((vector[size-1][j] * i + direct[j])%q)  
        calRSULoop(q,size-1,vector,temp,Vin)

# 计算R，从allProjectSpace筛选 n + 1 个向量 分别乘上 [0 - q-1]倍 加上 L中的元素
def RSUFilter(n,q,pos,RSU,vector,allProjectSpace):
    if pos>=n:
        direct = [0 for i in range(len(vector[0]))]
        Vin = [] #接受子空间的向量
        calRSULoop(q,n,vector,direct,Vin)
        Vin = sorted(Vin)   # 此处用于判别空间是否重复
        if Vin not in RSU:
            RSU.append(Vin)
        return
    global Lbase

    for i in range(pos,len(allProjectSpace)):
        vector.append(allProjectSpace[i])
        matrixTemp = []
        matrixTemp.extend(vector)
        matrixTemp.extend(Lbase)        
        if np.linalg.matrix_rank(matrixTemp) == len(Lbase) + 1 + pos:   # 计算机vector是线性无关的
            RSUFilter(n,q,pos+1,RSU,vector,allProjectSpace)
        vector.pop()

# 用递进的思路找出R,在找出 V(t维)的基础上扩展n维度，成R(n+t)维度
def increaseExpansionSpace(n,RSU,V):
    iteratorV = copy.deepcopy(V)
    global allProjectSpace
    for i in range(n):  # 迭代n次
        iteratorRSU = []
        # 列表中向量到字符串的映射
        mapDigitToStr = dict()
        for space in iteratorV:
            # 清空字典
            mapDigitToStr.clear()  
            for item in space:
                mapDigitToStr[''.join(str(i) for i in item)] = 1
            for item in allProjectSpace:
                if not mapDigitToStr.get(''.join(str(i) for i in item)) :
                    Vin = []
                    for i in range(q):
                        qitem = []
                        for j in range(len(item)):
                            qitem.append((i * item[j]) % q)
                        for a in space:
                            lis = []
                            for j in range(len(item)):
                                lis.append((qitem[j] + a[j]) % q)
                            if lis not in Vin and lis in allProjectSpace:
                                Vin.append(lis)
                                mapDigitToStr[''.join(str(i) for i in lis)] = 1
                    Vin = sorted(Vin)
                    if Vin not in iteratorRSU:
                        iteratorRSU.append(Vin)
        iteratorV = copy.deepcopy(iteratorRSU)        
        
    return iteratorV


print("输入向量长度k,素数幂q,基向量个数ans(k>=ans) m n t")
lis = list(map(int,input().strip().split()))
k,q,ans,m,n,t = lis[0],lis[1],lis[2],lis[3],lis[4],lis[5]
baseVector = [[0 for i in  range(k)] for i in range(ans)]
for i in range(ans):
   baseVector[i][i] = 1

print(baseVector)

vector = [0 for i in range(k)]
space1 = []
space2 = []
begin1 = time.time()
begin11 = time.clock()
traverseSpace(k,q,0,ans,baseVector,vector)
end1 = time.time()
end11 = time.clock()
# for item in space1:
#     print(item)
print(len(space1))
print(space2)
print(len(space2))
print("暴力寻找需要的时间：{0}  {1}".format(end1-begin1,end11-begin11))


allProjectSpace = []
allProjectSpace.append([0 for i in range(k)])
begin2 = time.time()
begin22 = time.clock()
getProjectiveSpace(k,q)
end2 = time.time()
end22 = time.clock()
# print(allProjectSpace)
print(len(allProjectSpace))
print("进制寻找需要的时间：{0}    {1}".format(end2-begin2,end22-begin22))

print("公式计算得全射影空间的元素个数：{0}".format(1 + calSpaceElem(ans,q)))

# 生成一个固定的t-1维子空间的基向量
Lbase = [[0 for i in range(k)] for i in range(t-1)]
# L 空间
L = []
L.append([0 for i in range(k)])
# visited 关系到L子空间的查找，务必一动全动
visited = []  
for i in range(t-1):
    Lbase[i][i] = 1
    visited.append(i)
print("空间L的基 ",Lbase)
getLSpace(k,t-1,q)
print("空间L的元素 ",L)

# 第一种生成V
V = []
# 计算V，从allProjectSpace筛选出不在L中的1个向量然后[0 - q-1]倍 加上 L中的元素
for item in allProjectSpace:
    if item not in L:
        Vin = []
        for i in range(q):
            qitem = []
            for j in range(len(item)):
                qitem.append((i * item[j]) % q)
            for a in L:
                lis = []
                for j in range(len(item)):
                    lis.append((qitem[j] + a[j]) % q)
                if lis not in Vin and lis in allProjectSpace:
                    Vin.append(lis)
        Vin = sorted(Vin)
        if Vin not in V:
            V.append(Vin)
print("第一次生成V:")
for item in V:
    print(item)
print(len(V))

# 第二种生成V
V = []
# 列表中向量到字符串的映射
mapDigitToStr = dict()
for item in L:
    mapDigitToStr[''.join(str(i) for i in item)] = 1
for item in allProjectSpace:
    if not mapDigitToStr.get(''.join(str(i) for i in item)) :
        Vin = []
        for i in range(q):
            qitem = []
            for j in range(len(item)):
                qitem.append((i * item[j]) % q)
            for a in L:
                lis = []
                for j in range(len(item)):
                    lis.append((qitem[j] + a[j]) % q)
                if lis not in Vin and lis in allProjectSpace:
                    Vin.append(lis)
                    mapDigitToStr[''.join(str(i) for i in lis)] = 1
        V.append(Vin)
print("第二次计算V：")
for item in V:
    print(item)
print(len(V))

# 生成 R,在L的基础上循环遍历找出所有的n+1 基向量与L线性组合
R = []
RSUFilter(n+1,q,0,R,[],allProjectSpace)
print("第一次计算R的元素")
for item in R:
    print(item)
print("R的元素个数:",len(R))

# 用递进的思路找出R,在找出 V(t维)的基础上扩展n维度，成R(n+t)维度
R = []
# 在V的基础上扩展n维成R
R = increaseExpansionSpace(n,R,V)
print("第二次计算R的元素")
for item in R:
    print(item)
print("R的元素个数:",len(R))


# 计算 S
S = []
# 模拟遍历
RSUFilter(m+1,q,0,S,[],allProjectSpace)
print("第一次计算S的元素")
for item in S:
    print(item)
print("S的元素个数：",len(S))

# 递增扩展,在V(t维)的基础上扩展m维成S (m+t维)
S = []
S = increaseExpansionSpace(m,S,V)
print("第二次计算S的元素")
for item in S:
    print(item)
print("S的元素个数：",len(S))

# 计算U 
U = []
# 在L(t-1维)基础上循环加上 m+n+2维度成 U(m+n+t+1)
RSUFilter(m+n+2,q,0,U,[],allProjectSpace)
print("第一次计算U的元素")
for item in U:
    print(item)
print("U的元素个数：",len(U))

# 递增扩展，在R(n+t) 或 S(m+t) ，两个之中更大的空间的基础上增加 m+1维 或n+1 成U (m+n+t+1)
U = []
if(n>m):
    U = increaseExpansionSpace(m+1,U,R)
else :
    U = increaseExpansionSpace(n+1,U,S)

print("第二次计算U的元素")
for item in U:
    print(item)
print("U的元素个数：",len(U))