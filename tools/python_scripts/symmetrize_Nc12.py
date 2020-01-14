#! python
import sys
from numpy import *


class symmetrize:

	Kvec = array([[0.00000000e+00, 0.00000000e+00],
				  [2.22044605e-16, 3.14159265e+00],
				  [1.04719755e+00, 2.09439510e+00],
				  [1.04719755e+00, 5.23598776e+00],
				  [2.09439510e+00, 1.04719755e+00],
				  [2.09439510e+00, 4.18879020e+00],
				  [3.14159265e+00, 0.00000000e+00],
				  [3.14159265e+00, 3.14159265e+00],
				  [4.18879020e+00, 2.09439510e+00],
				  [4.18879020e+00, 5.23598776e+00],
				  [5.23598776e+00, 1.04719755e+00],
				  [5.23598776e+00, 4.18879020e+00]])

		   #  MOMENTUM_SPACE k-space symmetries for Nc=12: 

		# 0, 0    |               0, 0    0, 0    0, 0    0, 0
		# 1, 0    |               1, 0    6, 0    6, 0    1, 0
		# 2, 0    |               2, 0    9, 0    4, 0    11, 0
		# 3, 0    |               3, 0    3, 0    10, 0   10, 0
		# 4, 0    |               4, 0    11, 0   2, 0    9, 0
		# 5, 0    |               5, 0    5, 0    8, 0    8, 0
		# 6, 0    |               6, 0    1, 0    1, 0    6, 0
		# 7, 0    |               7, 0    7, 0    7, 0    7, 0
		# 8, 0    |               8, 0    8, 0    5, 0    5, 0
		# 9, 0    |               9, 0    2, 0    11, 0   4, 0
		# 10, 0   |               10, 0   10, 0   3, 0    3, 0
		# 11, 0   |               11, 0   4, 0    9, 0    2, 0


	def __init__(self):
		self.data=[]
		self.Nc=12
		self.iSym=4
		self.setup_symm_tables()


	def setup_symm_tables(self):
		# 0, 0    |               0, 0    0, 0    0 , 0    0 , 0
		# 1, 0    |               1, 0    6, 0    6 , 0    1 , 0
		# 2, 0    |               2, 0    9, 0    4 , 0    11, 0
		# 3, 0    |               3, 0    3, 0    10, 0    10, 0
		# 4, 0    |               4, 0    11, 0   2 , 0    9 , 0
		# 5, 0    |               5, 0    5, 0    8 , 0    8 , 0
		# 6, 0    |               6, 0    1, 0    1 , 0    6 , 0
		# 7, 0    |               7, 0    7, 0    7 , 0    7 , 0
		# 8, 0    |               8, 0    8, 0    5 , 0    5 , 0
		# 9, 0    |               9, 0    2, 0    11, 0    4 , 0
		# 10, 0   |               10, 0   10, 0   3 , 0    3 , 0
		# 11, 0   |               11, 0   4, 0    9 , 0    2 , 0
		# self.symTrK[iK,iSym]
		self.symTrK = zeros((self.Nc, self.iSym),dtype='int')
		self.symTrK[0,0]= 0 ; self.symTrK[0,1]= 0 ; self.symTrK[0,2]=  0 ; self.symTrK[0,3]= 0  
		self.symTrK[1,0]= 1 ; self.symTrK[1,1]= 6 ; self.symTrK[1,2]=  6 ; self.symTrK[1,3]= 1  
		self.symTrK[2,0]= 2 ; self.symTrK[2,1]= 9 ; self.symTrK[2,2]=  4 ; self.symTrK[2,3]= 11 
		self.symTrK[3,0]= 3 ; self.symTrK[3,1]= 3 ; self.symTrK[3,2]=  10; self.symTrK[3,3]= 10 
		self.symTrK[4,0]= 4 ; self.symTrK[4,1]= 11; self.symTrK[4,2]=  2 ; self.symTrK[4,3]= 9  
		self.symTrK[5,0]= 5 ; self.symTrK[5,1]= 5 ; self.symTrK[5,2]=  8 ; self.symTrK[5,3]= 8  
		self.symTrK[6,0]= 6 ; self.symTrK[6,1]= 1 ; self.symTrK[6,2]=  1 ; self.symTrK[6,3]= 6  
		self.symTrK[7,0]= 7 ; self.symTrK[7,1]= 7 ; self.symTrK[7,2]=  7 ; self.symTrK[7,3]= 7  
		self.symTrK[8,0]= 8 ; self.symTrK[8,1]= 8 ; self.symTrK[8,2]=  5 ; self.symTrK[8,3]= 5  
		self.symTrK[9,0]= 9 ; self.symTrK[9,1]= 2 ; self.symTrK[9,2]=  11; self.symTrK[9,3]= 4  
		self.symTrK[10,0]=10; self.symTrK[10,1]=10; self.symTrK[10,2]= 3 ; self.symTrK[10,3]=3  
		self.symTrK[11,0]=11; self.symTrK[11,1]=4 ; self.symTrK[11,2]= 9 ; self.symTrK[11,3]=2  

	def apply_point_group_symmetries_Q0(self,G4):
		# for G4[w1,w2,K1,K2]
		# G4(K,K') = G4(Ra(K),Ra(K')) for all frequencies
		Nc = self.Nc
		nwn = G4.shape[0]
		type=dtype(G4[0,0,0,0])
		for iK1 in range(0,Nc):
			for iK2 in range(0,Nc):
				tmp = zeros((nwn,nwn),dtype=type)
				for iSym in range(0,self.iSym): # Apply every point-group symmetry operation
					iK1Trans = self.symTrK[iK1,iSym]
					iK2Trans = self.symTrK[iK2,iSym]
					tmp[:,:] += G4[:,:,iK1Trans,iK2Trans]
				for iSym in range(0,self.iSym):
					iK1Trans = self.symTrK[iK1,iSym]
					iK2Trans = self.symTrK[iK2,iSym]
					G4[:,:,iK1Trans,iK2Trans] = tmp[:,:]/float(self.iSym)

