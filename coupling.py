#!-*-coding:utf-8-*-
import numpy as np
import os
import re
#通过用一个log文件，统计其中MG，N2，N4坐标和偶极矩阵，以及激发态的值
class E_COUP:
    def __init__(self, INF1):
        self.INF1 = INF1 #初始化文件名
        #self.ev = 27.2116 #1au = 27.2116ev
        #self.cm = 219475  #1au = 219475cm-1
        #self.ang = 1.0/0.529177 #1au = 0.529177ang
    #求文件中的偶极矩和激发能
    def dipole_excited(self,INF ):
        import re
        import os
        INF = INF + ".log"
        #偶极矩
	order = "grep -C 3 'electric dipole moments' %s" %(INF) #将含有electric dipole moments 的行下面三列找出
        return_values = os.popen(order).read() #将找出的值返回
        dipole_line = return_values.split('\n')[5] #选出包含偶极矩的行
        dipole_values = re.findall(r"-?\d+\.\d+", dipole_line)[0:3] #选出偶极矩的值
        
        dip_vec = [float(i) for i in dipole_values] #将偶极矩的值设置为浮点性（字符型转浮点型），并将数据结构设置为numpy

	#激发能
	order = "grep 'Excited State   1' %s" %(INF)
	excited_line = os.popen(order).read()
	excited_values = float(re.findall(r"-?\d+\.\d+", excited_line)[0])	
	return dip_vec, excited_values 
    #求原子坐标, 在pdb文件中找原子坐标， 给出输入文件和原子名称
    def atom_coord(self, ATOM, INF):
        INF = INF + ".pdb"
        ATOM = ATOM
        ATOM_coord_order = "grep '%s' %s" %(ATOM, INF) #在pdb文件中包含原子的行
        ATOM_coord_line = os.popen(ATOM_coord_order).read().split("\n")[0] #将这一行转换为列表形式
        ATOM_coord = re.findall(r"-?\d+\.\d+", ATOM_coord_line)[0:3] #得到原子坐标
        coord = [float(i) for i in ATOM_coord] #将原子坐标转换为浮点型，和array格式
        return coord
    #调整偶极
    def adjust_dipole( self, INF1, INF2, INF3):
        import numpy as np
        #INF1 = self.INF1
        #INF2 = self.INF2
        #文件1中的MG 、N2、N4
        #MG_i = self.coord(INF1, "MG")
        #ND_2_i = self.coord(INF1, "NB")
        #NB_4_i = self.coord(INF1, "ND")
        #文件2中的MG、N2、N4
        #MG_j = self.coord(INF2, "MG")
        #ND_2_j = self.coord(INF2, "NB")
        #NB_4_j = self.coord(INF2, "ND")
        #求MG向量， N2到N4的向量
        #MG_ij = np.array([(MG_j[k]-MG_i[k]) for k in range(3)]) #确定向量,并将aug转化为ev
        #print INF1, INF2, INF3
        INF1 = np.array(INF1)
        N_v = np.array([(INF3[k]-INF2[k]) for k in range(3)]) #确定向量
        #print N_v
	if np.dot(N_v,INF1) < 0 :
		new_dipole = [-k for k in INF1]
	else:
		new_dipole = INF1
	return new_dipole
        #N_j = np.array([(NB_4_j[k]-ND_2_j[k]) for k in range(3)])
        #归一化MG向量
        #nom_MG = np.linalg.norm(MG_ij)
        #print nom_MG
        """
        #镁原子向量是需要除以模
        ZH = []
        for i in range(3):
             k = MG_ij[i]/nom_MG
             ZH.append(k)
        MG_ij_norm =  np.array(ZH)
        """
        #ZH = [(MG_ij[k]/nom_MG for k in range(3))]
        #print ZH 
        #return MG_ij,MG_ij_norm, N_i, N_j
        #return MG_ij, nom_MG, N_i, N_j
    #得到MG片与片之间的mg向量, INF1和INF2分别表示维度和MG总向量
    def MG_vec_nor(self, INF1, INF2):
        import numpy as np
        MG_vec = []
        MG_nor = []
        for i in range(INF1):
	    for j in range(INF1):
                if i == j :
                    MG_vec.append([0.0, 0.0, 0.0])
                    MG_nor.append([0.0])
                else:
                    vec = [INF2[j][k]-INF2[i][k] for k in range(3)]
                    #print vec , "MG"
	            MG_vec.append(vec)
                    nor = [np.linalg.norm(np.array(vec))]
                    MG_nor.append(nor)
        print MG_vec
        return MG_vec, MG_nor
    #得到激发能的值
    def Excited_values(self):
        import re
        INF = self.INF1 + ".log"
        #print INF
        """
        if "702" not in INF:
            order = "grep 'Excited State   1' %s" %(INF)
        else:
            order = "grep 'Excited State   2' %s" %(INF)
        """
        order = "grep 'Excited State   1' %s" %(INF)
        #print order
        line = os.popen(order).read()
        #print line
        #values 表示第一激发能的值，单位为电子福特（ev）
        #values = 8065.49413*float(re.findall(r"-?\d+\.\d+", line)[0])
        values = float(re.findall(r"-?\d+\.\d+", line)[0])
        #print values
        #print "****" 
        #vector 表示电子偶极向量
        vector = self.diople(self.INF1)
        #print vector
        #确定头部N2到N4的指向
        ND_2 = self.coord(self.INF1, "NB") #在文件中找出包含ND的行 
        #print ND_2
        NB_4 = self.coord(self.INF1, "ND") #在文件中找出包含NB的行 
        #print NB_4
        N_vec = np.array([(NB_4[i]-ND_2[i]) for i in range(3)]) #用NB_2的xyz 减去ND_2的xyz。 向量的方向ND指向NB
        #统一偶极的方向
        if np.dot(N_vec, vector) < 0:
               vector_copy = [-k for k in vector]
               vector = vector_copy
        print vector , ND_2, NB_4
        print np.dot(N_vec, vector), INF
        return values,vector
        
    #根据公式来开始计算
    def coup_values(self):
        import numpy as np
        INF1 = self.INF1
        INF2 = self.INF2
        Rij, Norm, Ni, Nj = self.vec(INF1,INF2)
        ui = self.diople(INF1)
        uj = self.diople(INF2)
        if np.dot(Ni, ui) < 0:
             neg_ui = [-k for k in ui]
             ui = np.array(neg_ui)
        if np.dot(Nj, uj) < 0:
             neg_uj = [-k for k in uj]
             uj = np.array(neg_uj) 
        #print np.dot(ui,uj)
        #print np.linalg.norm(Rij)
        #Cij = 219475*((np.dot(ui, uj) - 3*np.dot(ui, Rij_norm)*np.dot(uj, Rij_norm))/np.linalg.norm(Rij)**3)
        #Cij = ((np.dot(ui, uj) - 3*np.dot(ui, Rij_norm)*np.dot(uj, Rij_norm))/np.linalg.norm(Rij)**3) #将所有的单位都转换为了电子福特（ev）进行计算 
        Cij = (np.dot(ui,uj)/((Norm*self.ang)**3) - 3*(np.dot(ui,Rij)*np.dot(uj,Rij))/(Norm*self.ang)**3/Norm**2)*self.ev
        #print Cij
        return Cij
#k = E_COUP('/home4/huangxh/workspace/2fkw/BCL/HIS_BCL/crystal/ALZ_BCL/BCL501_ALZ54/BCL501_OPT', '/home4/huangxh/workspace/2fkw/BCL/HIS_BCL/crystal/ALZ_BCL/BCL502_ALZ59/BCL502_OPT')
#k.coup_values()
