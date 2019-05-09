#!-*-coding:UTF-8-*-
from coupling import E_COUP
import sys
import numpy as np
import pandas as pd
from sympy import *
from GSZK import * 
#******************************读取文件数据*****************************************************
# read filename to list
filename = sys.argv[1]
#print filename
Fname = [i.strip() for i in open(filename)]
# restruct 27*27 hartree matrix
dipole_matrix = np.array([])
hartree_matrix = np.array([])

MG_values = []
dipole_values = []
excited_values = []
for name in Fname:
	func = E_COUP(name)
	dipole, excited = func.dipole_excited(name)	#调用函数中的偶极激发模块，得到激发和偶极
	MG_coord = func.atom_coord('MG', name)   #得到MG的原子坐标 
	N2_coord = func.atom_coord('NB', name)   #得到N2的原子坐标
        N4_coord = func.atom_coord('ND', name)   #得到N4的原子坐标
	new_dipole = func.adjust_dipole(dipole, N2_coord, N4_coord) #通过调整偶极方向与N2指向N4的方向一致
	MG_values.append(MG_coord) #添加坐标
	dipole_values.append(new_dipole) #添加偶极
	excited_values.append(excited) #添加激发能


#********************************计算耦合******************************************
normal = 27
ang = 1.0/0.529177 #将ang单位转换为au单位
au = 27.2116 #将au单位转换为ev
MG_vec, MG_nor = func.MG_vec_nor(normal, MG_values) #得到MG向量和MG的模
ham = []
excited_values = np.array(excited_values)#将列表转化为array，并构造27*27矩阵
MG_vec = np.array(MG_vec).reshape(normal,normal,3)
MG_nor = np.array(MG_nor).reshape(normal,normal)
dipole_values = np.array(dipole_values).reshape(normal,3) #转换偶极矩
for i in range(normal):
    ui = dipole_values[i]
    for j in range(normal):
        print i, j
        uj = dipole_values[j]
        Norm = MG_nor[i][j]*ang #将mg的模单位转换为ang转换为au
        #Rij = [MG_vec[i][j][k]*ang for k in range(3)] #ang转换为au
        Rij = MG_vec[i][j]
        if i == j:
            ham.append(excited_values[i])
        else:
            Cij = ((np.dot(ui,uj)/Norm**3)-3*(np.dot(ui, Rij)*np.dot(uj,Rij))/(Norm**3)/MG_nor[i][j]**2)*au #计算电子耦合
            print ui,'ui', uj,'uj', Norm,'Norm', Rij,'Rij', MG_nor[i][j],'MGnor'
	    ham.append(Cij)
#**************************************************计算光谱***********************************************
ham = pd.DataFrame(np.array(ham).reshape(normal,normal))
dipole_values = pd.DataFrame(dipole_values)
D, V, newdipole, strength,lamba, energy = spectrum(ham, dipole_values) #调用spectrum函数计算光谱强度

#**************************************************转化数据类型***********************************************
strength = pd.DataFrame(strength)
lamba = pd.DataFrame(lamba)
energy = pd.DataFrame(energy)
D = pd.DataFrame(D)
V = pd.DataFrame(V)
newdipole = pd.DataFrame(newdipole)

#***********************************************将数据写入文件***********************************************
with pd.ExcelWriter("TEXT.xlsx") as writer:
    ham.to_excel(writer, sheet_name='Hamilton') #将哈密顿矩阵写入
    dipole_values.to_excel(writer,sheet_name='dipole') #将偶极矩写入
    D.to_excel(writer, sheet_name='D')#将特征矩阵写入
    V.to_excel(writer, sheet_name='V')#将特征向量写入
    newdipole.to_excel(writer,sheet_name='newdipole')#将新偶极写入
    strength.to_excel(writer,sheet_name='strength')#写入强度---纵坐标
    lamba.to_excel(writer,sheet_name='lamba')#横坐标
    energy.to_excel(writer,sheet_name='energy')#横坐标





