#import 所需的套件
import math
import numpy as np
from scipy.stats import norm
import scipy.linalg
import random


#參數
K = 100
T = 0.5
N = 5
r = 0.1
times = 10000
rep = 20


#計算平均數
def Mean(num,N):
    total = 0.0
    for i in range(N):
        total += num[i]
    mean = total/N
    return mean


#計算Sigma
def Var(num,N):
    total = 0.0
    mean = Mean(num,N)
    for i in range(N):
        total += (num[i] - mean)**2
    var = (total / N)**(1/2)
    return var


#計算相關係數
def Cov(num1,num2,N):
    mean1 = Mean(num1,N)
    mean2 = Mean(num2,N)
    sigma1 = Var(num1,N)
    sigma2 = Var(num2,N)
    total = 0.0
    for i in range(N):
        total += (num1[i] - mean1)*(num2[i] - mean2)
    cov = (total/N)/(sigma1*sigma2)
    return cov


#個產品的資料
sList = [0.0] * N
qList = [0.0] * N
sigmaList = [0.0] * N
relation = [[0.0] * N for i in range(N)]

sList = [95,95,95,95,95]
qList = [0.05,0.05,0.05,0.05,0.05]
sigmaList = [0.5,0.5,0.5,0.5,0.5]
relation = [[1,0.5,0.5,0.5,0.5],[0.5,1,0.5,0.5,0.5],[0.5,0.5,1,0.05,0.05],[0.5,0.5,0.5,1,0.5],[0.5,0.5,0.5,0.5,1]]


#變成ＣＭatrix
def Cmatrix(relation,sigmaList,cMatrix,N,T):
    for i in range(N):
        for j in range(N):
            cMatrix[i][j] = sigmaList[i] * sigmaList[j] * relation[i][j] *  T 

            
#計算價格
def Price(S,r,q,sigma,R,T,K):
    mu = math.log(S)+(r-q-((sigma)**(2))/2)*T
    price = math.exp(mu+R)
    if price > K:
        return (price - K) 
    else:
        return 0


def RainbowOption_by_Cholesky(K,T,N,r,times,rep,sList,qList,sigmaList,relation):
    #建立ＣMatrix
    cMatrix = np.zeros((N,N))
    Cmatrix(relation,sigmaList,cMatrix,N,T)
    #aMatrix = Cholesky decomposition
    aMatrix =  np.linalg.cholesky(cMatrix).T
    
    #紀錄結果
    meanList = [0.0] * rep
    
    #模擬rep次
    for i in range(rep):
        
        #模擬股價 "r"（times次）
        num = np.zeros((times,N))
        for j in range(times):
            for k in range(N):
                num[j][k] = random.gauss(mu = 0, sigma = 1)
        
        SPrice = [0.0]*times
        rMatrix =  np.matmul(num,aMatrix)
        for j in range(times):
            #計算該次價格
            maxPrice = [0.0] * N
            for k in range(N):  
                maxPrice[k] = Price(sList[k],r,qList[k],sigmaList[k],rMatrix[j][k],T,K)
            SPrice[j] = max(maxPrice)*math.exp(-r*T)
        
        #記錄平均值
        meanList[i] = Mean(SPrice,times)
    
    
    allmean = Mean(meanList,rep)
    allVar = Var(meanList,rep)

    print("Basic Requirement")
    print("Mean : "+str(allmean)+" Var : "+str(allVar))
    print("95%信賴區間 ： "+str(allmean - 2 * allVar)+" ~ "+str(allmean + 2 * allVar))


def RainbowOption_by_Cholesky_momentMatching(K,T,N,r,times,rep,sList,qList,sigmaList,relation):
    #建立ＣMatrix
    cMatrix = np.zeros((N,N))
    Cmatrix(relation,sigmaList,cMatrix,N,T)
    #aMatrix = Cholesky decomposition
    aMatrix =  np.linalg.cholesky(cMatrix).T
    
    #紀錄結果
    meanList = [0.0] * rep
    
    #模擬rep次
    for i in range(rep):
        
        #模擬股價 "r"（times次）
        num = np.zeros((N,times))
        for j in range(N):
            for k in range(int(times/2)):
                num[j][k] = random.gauss(mu = 0, sigma = 1)
                num[j][k + int(times/2)] = -num[j][k]
    
        Zmean = [0.0] * N
        Zvar = [0.0] * N
        for j in range(N):
            Zmean[j] = Mean(num[j],times)
            Zvar [j] = Var(num[j],times)
            #print(str(Zmean[j])+" "+str(Zvar[j]))
    
        Znum = np.zeros((times,N))
        for j in range(times):
            for k in range(N):
                Znum[j][k] = (num[k][j] - Zmean[k])/Zvar[k] 
        
        SPrice = [0.0]*times
        rMatrix =  np.matmul(Znum,aMatrix)
        for j in range(times):
            #計算該次價格
            maxPrice = [0.0] * N
            for k in range(N):  
                maxPrice[k] = Price(sList[k],r,qList[k],sigmaList[k],rMatrix[j][k],T,K)
            SPrice[j] = max(maxPrice)*math.exp(-r*T)
        
        #記錄平均值
        meanList[i] = Mean(SPrice,times)
    
    
    allmean = Mean(meanList,rep)
    allVar = Var(meanList,rep)

    print("Bonus 1 ")
    print("Mean : "+str(allmean)+" Var : "+str(allVar))
    print("95%信賴區間 ： "+str(allmean - 2 * allVar)+" ~ "+str(allmean + 2 * allVar))


def InverseCholesky(Znum,N,times,Zrelation,ZsigmaList):    
    #sigma
    ZsigmaList = [0.0] * N
    for j in range(N):
        ZsigmaList[j] = Var(Znum[j],times)
    #關係
    Zrelation = np.zeros((N,N))
    for j in range(N):
        for k in range(N):
            Zrelation[j][k] = Cov(num[j],num[k],times)
    #轉回去
    Znum = Znum.T
    bMatrix = np.zeros((N,N))
    Cmatrix(Zrelation,ZsigmaList,bMatrix,N,T)
    bMatrix =  np.linalg.cholesky(bMatrix).T
    bMatrix = np.linalg.inv(bMatrix)

    return bMatrix


def RainbowOption_by_Cholesky_InverseCholesky(K,T,N,r,times,rep,sList,qList,sigmaList,relation):
    #建立ＣMatrix
    cMatrix = np.zeros((N,N))
    Cmatrix(relation,sigmaList,cMatrix,N,T)
    #aMatrix = Cholesky decomposition
    aMatrix =  np.linalg.cholesky(cMatrix).T
    
    #紀錄結果
    meanList = [0.0] * rep
    
    for i in range(rep):
        num = np.zeros((N,times))
        #模擬股價(basic)
        
        for j in range(N):
            for k in range(int(times/2)):
                num[j][k] = random.gauss(mu = 0, sigma = 1)
                num[j][k + int(times/2)] = -num[j][k]    
        Zmean = [0.0] * N
        Zvar = [0.0] * N
        for j in range(N):
            Zmean[j] = Mean(num[j],times)
            Zvar [j] = Var(num[j],times)

    
        Znum = np.zeros((N,times))
        for j in range(N):
            for k in range(times):
                Znum[j][k] = (num[j][k] - Zmean[j])/Zvar[j] 
    
        bMatrix = InverseCholesky(Znum,N,times,Zrelation,ZsigmaList)
        rMatrix =  np.matmul(Znum.T,aMatrix)
        SPrice = [0.0]*times
        for j in range(times):
            maxPrice = [0.0] * N
            for k in range(N):  
                maxPrice[k] = Price(sList[k],r,qList[k],sigmaList[k],rMatrix[j][k],T,K)
            SPrice[j] = max(maxPrice)*math.exp(-r*T)
        
        meanList[i] =  Mean(SPrice,times)
    
    allmean = Mean(meanList,rep)
    allVar = Var(meanList,rep)

    print("Bonus 2 ")
    print("Mean : "+str(allmean)+" Var : "+str(allVar))
    print("95%信賴區間 ： "+str(allmean - 2 * allVar)+" ~ "+str(allmean + 2 * allVar))
    
#所有答案
print("Ans Block")
RainbowOption_by_Cholesky(K,T,N,r,times,rep,sList,qList,sigmaList,relation)
RainbowOption_by_Cholesky_momentMatching(K,T,N,r,times,rep,sList,qList,sigmaList,relation)
RainbowOption_by_Cholesky_InverseCholesky(K,T,N,r,times,rep,sList,qList,sigmaList,relation)

