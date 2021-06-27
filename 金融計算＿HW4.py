#!/usr/bin/env python
# coding: utf-8

# In[671]:


#import 所需的套件
import math
import numpy as np
from scipy.stats import norm
import random


# In[672]:


#繼承方法
def inherit(numList1,numList2,St,myList):
    #繼承的list
    newlist = []
    for i in numList1:
        #如果 i = 0就結束
        if i == 0:
            break
        #大於 i 就繼承
        if i > St :
            newlist.append(i)
        #小於 i 就補自己
        else:
            newlist.append(St)
            break        
    for i in numList2:
        #如果 i = 0就結束
        if i == 0:
            break
        #大於 i 且不存仔就繼承
        if (i >= St)and(i not in newlist) :
            newlist.append(i)
        else:
            #小於 i 且不存在就補自己
            if (i not in newlist)and(St not in newlist):
                newlist.append(St)
    #寫入答案
    num = 0
    #print(newlist)
    for i in newlist:
        myList[num] = i
        num += 1    


# In[673]:


#寫入選擇價格
def putprice(Smaxprice,Soption,St):
    #寫入的list
    newList = []
    for i in Smaxprice :
        # i = 0 就停止
        if i == 0:
            break
        #寫入價格或補0
        price = i - St
        if price > 0:
            newList.append(price)
        else:
            newList.append(0)
    #輸出答案
    num = 0
    
    for i in newList:
        Soption[num] = i 
        num += 1


# In[674]:


#往後回推
def backprice(St,u,d,P,Myoption,MyMax,Soption1,SMax1,Soption2,SMax2,r,dT,typeof):
    #需先找到正確的上下價格
    upprice = 0.0
    downprice = 0.0
    num = 0
    for i in MyMax:
        if i == 0.0:
            break
        passnext = 0#是否找得到
        for j in range(len(SMax1)):
            if abs(i - SMax1[j]) < 0.0000001:
                upprice = Soption1[j]
                passnext = 1
        #如果找不到的話
        for j in range(len(SMax1)):
            if (abs(St * u - SMax1[j])<0.0000001)and(passnext == 0):
                upprice = Soption1[j]
        
        for j in range(len(SMax2)):
            if abs(i - SMax2[j]) < 0.0000001:
                downprice = Soption2[j]
                break
        optionePrice = (P * upprice + (1-P) * downprice)*math.exp(-r*dT)
        #美式 or 歐式
        if typeof == "American":
            if (i - St) > optionePrice:    
                optionePrice = i -St
        else:
            pass
        Myoption[num] = optionePrice
        num += 1       


# In[675]:


def bonus(St,u,d,P,n,SMax):
    AllList = []
    for i in range(n):
        iList = []
        for j in range(i+1):
            jList = []
            if j == 0:
                if St*u**(i - j )*d**(j) > SMax:
                    jList.append(St*u**(i - j )*d**(j))
                else:
                    jList.append(SMax)
            elif i == j:
                jList.append(SMax)
            
            elif i - j  > j:
                for k in range(j+1):
                    if SMax > St*u**(i - j - k):
                        jList.append(SMax)
                        break
                    jList.append(St*u**(i - j - k))
            else:
                for k in range(i-j+1):
                    if SMax > St*d**(k - i + j):
                        jList.append(SMax)
                        break
                    jList.append(St*d**(k - i + j))
            iList.append(jList)
        AllList.append(iList)
    return AllList


# In[676]:


#計算CRR Binomial
def binomial(St,r,q,sigma,t,T,Smax,n):
    #計算u,d,P
    dT = (T - t)/n
    u = math.exp(sigma*((dT)**(1/2))) 
    d = 1.0/u 
    P = ((math.exp((r-q)*(dT)))-d)/(u-d) 
    #創造所需的矩陣
    SpiceMatrix = np.zeros((n+1,n+1))
    SmaxMatrix = np.zeros((n+1,n+1,2*(n+2)))
    Soption =  np.zeros((n+1,n+1,2*(n+2)))
    numSmax = np.zeros((n+1,n+1,2*(n+2)))
    #寫入Spricematrix，只需最後兩行，其他都用其複製
    for j in range(n + 1):
        if (n - j) > j:
            SpiceMatrix[n][j] = St * (u)** (n - 2*j)
        else:
            SpiceMatrix[n][j] = St * (d)** (2*j - n)
    for j in range(n):
        if (n - 1 - j) > j:
            SpiceMatrix[n-1][j] = St * (u)** (n - 2*j - 1)
        else:
            SpiceMatrix[n-1][j] = St * (d)** (2*j - n + 1)
    #往前複製
    for i in range(n - 2,-1,-1):
        for j in range(i+1):
            SpiceMatrix[i][j] = SpiceMatrix[i+2][j+1]
    #初始化時間t，如果有Smax需要改
    if St > Smax:
        SmaxMatrix[0][0][0] = St
    else:
        SmaxMatrix[0][0][0] = Smax
    #計算個節點的Smax
    for i in range(1,n + 1):
        for j in range(0,i+1):
            #最下面
            if i == j:
                SmaxMatrix[i][j][0] = SmaxMatrix[0][0][0]
            #最上面，且要注意Smax at t
            elif j == 0:
                if SpiceMatrix[i][j] > SmaxMatrix[0][0][0]:
                    SmaxMatrix[i][j][0] = SpiceMatrix[i][j]
                else:
                    SmaxMatrix[i][j][0] = SmaxMatrix[0][0][0]
            #做其他的SmaxMatrix[i][j]
            else:    
                inherit(SmaxMatrix[i-1][j-1],SmaxMatrix[i-1][j],SpiceMatrix[i][j],SmaxMatrix[i][j])
    SmaxMatrix = bonus(St,u,d,P,n+1,Smax)
    #最後一期的定價
    #print(SmaxMatrix)
    for i in range(n+1):
        putprice(SmaxMatrix[n][i],Soption[n][i],SpiceMatrix[n][i])
    #回推價格
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            backprice(SpiceMatrix[i][j],u,d,P,Soption[i][j],SmaxMatrix[i][j],Soption[i+1][j],SmaxMatrix[i+1][j],Soption[i+1][j+1],SmaxMatrix[i+1][j+1],r,dT,"American")
    print("Basic 1 : binomial tree")
    print("美式的價格 ："+str(Soption[0][0][0]))        
    print("-----------------------------------")
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            backprice(SpiceMatrix[i][j],u,d,P,Soption[i][j],SmaxMatrix[i][j],Soption[i+1][j],SmaxMatrix[i+1][j],Soption[i+1][j+1],SmaxMatrix[i+1][j+1],r,dT,"European")
    print("歐式的價格 ："+str(Soption[0][0][0])) 
#------------------------------------------------------------------------------------------
    SmaxMatrix = bonus(St,u,d,P,n+1,Smax)
    #最後一期的定價
    #print(SmaxMatrix)
    for i in range(n+1):
        putprice(SmaxMatrix[n][i],Soption[n][i],SpiceMatrix[n][i])
    #回推價格
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            backprice(SpiceMatrix[i][j],u,d,P,Soption[i][j],SmaxMatrix[i][j],Soption[i+1][j],SmaxMatrix[i+1][j],Soption[i+1][j+1],SmaxMatrix[i+1][j+1],r,dT,"American")
    print("Bonus 1  : binomial tree")
    print("美式的價格 ："+str(Soption[0][0][0]))        
    print("-----------------------------------")
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            backprice(SpiceMatrix[i][j],u,d,P,Soption[i][j],SmaxMatrix[i][j],Soption[i+1][j],SmaxMatrix[i+1][j],Soption[i+1][j+1],SmaxMatrix[i+1][j+1],r,dT,"European")
    print("歐式的價格 ："+str(Soption[0][0][0])) 


# In[677]:


#計算平均數
def Mean(num,N):
    total = 0.0
    for i in range(N):
        total += num[i]
    mean = total/N
    return mean


# In[678]:


#計算Sigma
def Var(num,N):
    total = 0.0
    mean = Mean(num,N)
    for i in range(N):
        total += (num[i] - mean)**2
    var = (total / N)**(1/2)
    return var


# In[679]:


#實作蒙地卡羅法
def Monte_Carlo(St,r,q,sigma,t,T,Smax,n,simu,rep):
    # rep次試驗的list
    meanList = []
    dT = (T-t)/n
    plus = np.array([(r - q - (sigma**2)/2) * dT]*simu)
    for i in range(rep):
        thisNum = np.zeros((n+1,simu))
        #都先取自然對數處理
        for j in range(simu):
            thisNum[0][j] = math.log(St)

        #模擬股價的自然對數
        for j in range(n):
            x = np.random.normal((thisNum[j] + plus),((sigma**2) * dT)**(1/2))
            
            #匯出至下一期
            thisNum[j+1] = x
        #轉正
        thisNum = thisNum.T
        thisNow = []#現在價格
        thisMax = []#最大價格
        for j in range(simu):
            thisMax.append(math.exp(np.max(thisNum[j])))
            thisNow.append(math.exp(thisNum[j][n]))
        #價格
        thisAns = []
        for j in range(simu):
            if(Smax > thisMax[j]):
                thisAns.append((Smax - thisNow[j])*math.exp(-(r)*(T-t)))
            else:
                thisAns.append((thisMax[j] - thisNow[j])*math.exp(-(r)*(T-t)))
        meanList.append(Mean(thisAns,simu))
    
    allmean = Mean(meanList,rep)
    allVar = Var(meanList,rep)
    print("蒙地卡羅法：")
    print("Mean : "+str(allmean)+" Var : "+str(allVar))
    print("95%信賴區間 ： "+str(allmean - 2 * allVar)+" ~ "+str(allmean + 2 * allVar))
    


# In[709]:


def CheukMethod(St,r,q,sigma,t,T,Smax,n):
    #計算u,d,P
    dT = (T - t)/n
    u = math.exp(sigma*((dT)**(1/2))) 
    d = 1.0/u 
    mu = math.exp((r-q)*(dT))
    P = (mu*u-1)/(mu*(u-d)) 
    #創造所需的矩陣   
    Soption =  np.zeros((n+1,n+1))
    for i in range(n + 1):
        Soption[n][i] = (u**(i)-1)
        #print(Soption[n][i])
    #print()
    #back
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            #print(str(i)+" "+str(j))
            if j == 0:
                Soption[i][j] = ((P) * Soption[i+1][0] + (1-P) * Soption[i+1][1]) * math.exp(-(q)*dT)
            else:
                Soption[i][j] = ((P) * Soption[i+1][j-1] + (1-P) * Soption[i+1][j+1]) * math.exp(-(q)*dT)
            if Soption[i][j] < (u**(j)-1):
                Soption[i][j] = (u**(j)-1)
            #print(Soption[i][j])
    print("Bpnus 2 : 美式")
    print(Soption[0][0]*St)
    print("-------------------")
    for i in range(n-1,-1,-1):
        for j in range(i+1):
            #print(str(i)+" "+str(j))
            if j == 0:
                Soption[i][j] = ((P) * Soption[i+1][0] + (1-P) * Soption[i+1][1]) * math.exp(-(q)*dT)
            else:
                Soption[i][j] = ((P) * Soption[i+1][j-1] + (1-P) * Soption[i+1][j+1]) * math.exp(-(q)*dT)
            #print(Soption[i][j])
    print("Bpnus 2 : 歐式")
    print(Soption[0][0]*St)            


# In[692]:


#輸入參數
St = 50
r = 0.1
q = 0.0
sigma = 0.4 
t = 0.2
T = 0.45
Smax = [50.0,60.0,70.0]
n = 100
simu = 10000
rep = 20


# In[695]:


for i in Smax:
    print("Smax,t = "+str(i))
    binomial(St,r,q,sigma,t,T,i,n)
    Monte_Carlo(St,r,q,sigma,t,T,i,n,simu,rep)


# In[711]:


CheukMethod(St,r,q,sigma,0.2,0.45,St,1000)


# In[ ]:




