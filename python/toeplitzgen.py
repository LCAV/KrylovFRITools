# This class stores a Toeplitz matrix implicitly (generator) and provide fast
# methods to evaluate its column-space and the roots of its generator
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs
import scipy.fftpack as fftp
from scipy.linalg import svd, eig, lstsq
class toeplitzgen:
	def __init__(self, col,row):
		self.col=np.asarray(col).reshape((-1))
		self.row=np.asarray(row).reshape((-1))
		self.M=self.row.size
		assert self.M==self.col.size
		assert self.row[0]==self.col[0]

		self.fc=fftp.fft(np.concatenate((self.col,np.zeros(1),self.row[-1:0:-1])))
		self.fcconj=fftp.fft(np.concatenate((np.conj(self.row),np.zeros(1),np.conj(self.col[-1:0:-1]))))
		self.linop=LinearOperator((self.M,self.M),matvec=self.vmult,rmatvec=self.vmult,matmat=self.mmult)
		self.TOL=1e-10

	def vmult(self,x):
		g=fftp.ifft(self.fc*fftp.fft(np.concatenate((x.reshape((-1)),np.zeros((self.M))))))
		g[self.M:]=0;
		return fftp.ifft(self.fcconj*fftp.fft(g))[0:self.M]

	def mmult(self,A):
		P=A.shape[1]
		G=fftp.ifft(self.fc.reshape((-1,1))*fftp.fft(np.concatenate((A,np.zeros((self.M,P)))),None,0),None,0)
		G[self.M:][:]=0;
		return fftp.ifft(self.fcconj.reshape((-1,1))*fftp.fft(G,None,0),None,0)[0:self.M,:]

	def colSpace(self,K):
		evl, evct = eigs(self.linop,K,tol=self.TOL)
		return evct

	def ESPRIT_TLS(self,K):
		V0=eigs(self.linop,K,tol=self.TOL)[1]
		C=svd(np.concatenate((V0[0:-1,:],V0[1:,:]),1),False)[2].T
		return -

