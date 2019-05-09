#!-*-coding:utf-8-*- 
def spectrum(ham, diople):
    import pandas as pd
    import numpy as np
    from math import *
    #diop = pd.read_csv(diople, header=None)
    #print diop
    #harmid = pd.read_csv(ham, header=None)
    #print diop
    hams = ham
    dioples = diople
    #print diople
    print diople.shape, ham.shape
    D, V = np.linalg.eig(hams)
    #排序,调整特征值和特征向量的区别
    sor = np.argsort(D)
    D_i = D[sor]
    D = D_i
    #调整特征向量的序号
    V_i = V[:,sor]
    V = V_i
    V= V.T
    newdipole = np.dot(V, dioples) #distruct a new diploe
    print newdipole
    nstep = 1000
    delt=0.0124*2 #the number of gaussian depression
    n = 27 #amplitude number
    #print D[n-1]
    print D[n-1], D[0]
    step = (D[n-1]-D[0])/nstep
    k = 0
    strengths = [] # y
    energy = [] # x
    lamda = [] # x 
    start = D[0]-300*step
    end = D[n-1]+300*step
    for i in np.arange(start, end , step):
        #print i*2, step
        num_strength = 0.0
        #print i
        for j in range(n):
            #print newdipole[j], np.linalg.norm(newdipole[j])
            osc = np.linalg.norm(newdipole[j])**2
            flag = 0.0245*i*exp(-(D[j]-i)**2/(2*delt**2))*osc/(sqrt(2*3.14)*delt)/n
            #print flag
            num_strength = num_strength + flag
        #print "******"
        #print num_strength
        #print "******"
        strengths.append(num_strength)
        energy.append(i)
        lamda.append(1240/(i/8065.48))
    #print len(strengths)
    print end
    strengths = np.array(strengths).reshape(-1)
    energy = np.array(energy).reshape(-1)
    lamda = np.array(lamda).reshape(-1)
    #np.savetxt('strength.txt', strengths)
    #fp = open("strengths", w)
    #for i in strengths:
          
   # print strengths
    #return strengths
    return D,V,newdipole,strengths, lamda, energy
#spectrum('hamid.csv', 'diople.csv')
