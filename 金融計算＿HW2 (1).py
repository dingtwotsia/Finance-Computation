#!/usr/bin/env python
# coding: utf-8

# In[653]:


#import 所需的套件
import math
import numpy as np
from scipy.stats import norm
import random
#答案儲存格


# In[654]:


#計算d1、d2
def Nd1(S,K,r,q,T,sigma):
    a = 0
    a = (math.log(S/K)+(r-q+(sigma**2)/2)*T)/(sigma*(T)**(1/2))
    return a
def Nd2(S,K,r,q,T,sigma):
    a = 0
    a = (math.log(S/K)+(r-q-(sigma**2)/2)*T)/(sigma*(T)**(1/2))
    return a


# In[655]:


#Apply 布雷克休斯模型
def black_scholes_Call(S,K,r,q,sigma,T):
    C =  S * math.exp(-q*T)*(norm.cdf(Nd1(S,K,r,q,T,sigma))) - K * math.exp(-r*T)*norm.cdf(Nd2(S,K,r,q,T,sigma))
    print("European Call : "+str(C))
def black_scholes_Put(S,K,r,q,sigma,T):
    P =  K * math.exp(-r*T)*(norm.cdf(-Nd2(S,K,r,q,T,sigma)))- S * math.exp(-q*T)*norm.cdf(-Nd1(S,K,r,q,T,sigma))
    print("European Put : "+str(P))                                                                                        


# In[656]:


#Monte_Carl Call
def Monte_Carlo_Call(S,K,r,q,sigma,T,N,M):
    meanList = []
    allMean = 0
    for i in range(N):
        num = []
        #模擬股價
        for j in range(M):
            x = random.gauss(mu = math.log(S) + (r - q - (sigma**2)/2) * T, sigma = ((sigma**2) * T)**(1/2))
            num.append(math.exp(x))
        for j in range(M):
            if num[j] > K :
                num[j] = (num[j] - K) * math.exp(-r*T)
            else:
                num[j] = 0
        mean  = 0    
        for j in range(M):
            mean += num[j]
        mean = mean / M
        meanList.append(mean)
        allMean += mean
    allMean = allMean/N
    print("European Call")
    print("Mean : "+str(allMean))
    allVar = 0
    for i in range(N):
        allVar = (meanList[i] - allMean)**2
    allVar = (allVar/N)**(1/2)
    print("Var : "+str(allVar))
    print("95%信賴區間 : "+str(allMean-2*allVar)+"~"+str(allMean+2*allVar))


# In[657]:


#Monte_Carl Put
def Monte_Carlo_Put(S,K,r,q,sigma,T,N,M):
    meanList = []
    allMean = 0
    for i in range(N):
        num = []
        #模擬股價
        for j in range(M):
            x = random.gauss(mu = math.log(S) + (r - q - (sigma**2)/2) * T, sigma = ((sigma**2) * T)**(1/2))
            num.append(math.exp(x))
        for j in range(M):
            if num[j] < K :
                num[j] = (K - num[j]) * math.exp(-r*T)
            else:
                num[j] = 0
        mean  = 0    
        for j in range(M):
            mean += num[j]
        mean = mean / M
        meanList.append(mean)
        allMean += mean
    allMean = allMean/N
    print("European Put")
    print("Mean : "+str(allMean))
    allVar = 0
    for i in range(N):
        allVar = (meanList[i] - allMean)**2
    allVar = (allVar/N)**(1/2)
    print("Var : "+str(allVar))
    print("95%信賴區間 : "+str(allMean-2*allVar)+"~"+str(allMean+2*allVar))


# In[658]:


#算Call的價格
def call(S,K):
    C = 0.0
    if(S > K):
        C = S - K
    else:
        C = 0.0
    return C  


# In[659]:


def Put(S,K):
    P = 0.0
    if(S < K):
        P = K - S
    else:
        P = 0.0
    return P


# In[671]:


def CRR_ECall(S,K,r,q,sigma,T,N,M):
    dT = T/(M) 
    u = math.exp(sigma*((dT)**(1/2))) 
    d = 1.0/u 
    P = ((math.exp((r-q)*(dT)))-d)/(u-d) 
    priceMatrix = [0]*(M+1)
    for i in range(M+1):
        priceMatrix[i] = 0.0
        priceMatrix[i] = S * (u**(M-i))*(d**(i))
        priceMatrix[i] = call(priceMatrix[i],K)
    for i in range(M+1):
        for j in range(M-i):
            priceMatrix [j] = math.exp(-r*dT)*(P*priceMatrix[j]+(1-P)*priceMatrix[j+1])
    print("Eurpean Call : "+str(priceMatrix[0]))


# In[672]:


def CRR_EPut(S,K,r,q,sigma,T,N,M):
    dT = T/(M) 
    u = math.exp(sigma*((dT)**(1/2))) 
    d = 1.0/u 
    P = ((math.exp((r-q)*(dT)))-d)/(u-d) 
    priceMatrix = [0]*(M+1)
    for i in range(M+1):
        priceMatrix[i] = S * (u**(M-i))*(d**(i))
        priceMatrix[i] = Put(priceMatrix[i],K)
    for i in range(M+1):
        for j in range(M-i):
            priceMatrix [j] = math.exp(-r*dT)*(P*priceMatrix[j]+(1-P)*priceMatrix[j+1])
    print("European Put : "+str(priceMatrix[0]))


# In[673]:


def CRR_ACall(S,K,r,q,sigma,T,N,M):
    dT = T/(M) 
    u = math.exp(sigma*((dT)**(1/2))) 
    d = 1.0/u 
    P = ((math.exp((r-q)*(dT)))-d)/(u-d) 
    priceMatrix = [0]*(M+1)

    for i in range(M+1):
        priceMatrix[i] = S * (u**(M-i))*(d**(i))
        priceMatrix[i] = call(priceMatrix[i],K)
    for i in range(M+1):
        for j in range(M-i):  
            priceAfter = math.exp(-r*dT)*(P*priceMatrix[j]+(1-P)*priceMatrix[j+1])
            priceBefore = S * (u**(M-i-1-j))*(d**(j)) - K
            if(priceAfter > priceBefore):
                priceMatrix [j] = priceAfter
            else:
                priceMatrix [j] = priceBefore
    print("Amreican Call : "+str(priceMatrix[0]))


# In[674]:


def CRR_APut(S,K,r,q,sigma,T,N,M):
    dT = T/(M) 
    u = math.exp(sigma*((dT)**(1/2))) 
    d = 1.0/u 
    P = ((math.exp((r-q)*(dT)))-d)/(u-d) 
    priceMatrix = [0]*(M+1)
    for i in range(M+1):
        priceMatrix[i] = S * (u**(M-i))*(d**(i))
        priceMatrix[i] = Put(priceMatrix[i],K)
    for i in range(M+1):
        for j in range(M-i):  
            priceAfter = math.exp(-r*dT)*(P*priceMatrix[j]+(1-P)*priceMatrix[j+1])
            priceBefore = K - S * (u**(M-i-1-j))*(d**(j))
            if(priceAfter > priceBefore):
                priceMatrix [j] = priceAfter
            else:
                priceMatrix [j] = priceBefore
    print("American Put : "+str(priceMatrix[0]))


# In[676]:


def OptionPrice(S,K,r,q,sigma,T,N,M):
    print("By black-Scholes")
    black_scholes_Call(S,K,r,q,sigma,T)
    black_scholes_Put(S,K,r,q,sigma,T)
    print("BY Monte Carlo")
    Monte_Carlo_Call(S,K,r,q,sigma,T,N,M)
    Monte_Carlo_Put(S,K,r,q,sigma,T,N,M)
    print("By CRR Model")
    CRR_ECall(S,K,r,q,sigma,T,N,M)
    CRR_EPut(S,K,r,q,sigma,T,N,M)
    CRR_ACall(S,K,r,q,sigma,T,N,M)
    CRR_APut(S,K,r,q,sigma,T,N,M)


# In[678]:


#參數們
S = 50
K = 50
r = 0.1
q = 0.05
T = 0.5
sigma = 0.4
N = 20
M = 500
OptionPrice(S,K,r,q,sigma,T,N,M)


# In[ ]:




