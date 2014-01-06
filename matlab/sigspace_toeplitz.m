function [V D]=sigspace_toeplitz(Y,K)
% INPUT
%	Y :			A Px(2*M+1) array of measurements (see measurement model above)
%	K :			The dimension of the signal space
%
% OUTPUT
%	V :			An orthonormal basis of the signal space
%	D :			The modulus of the associated singular values
  N=size(Y,2);
  L=floor(N/2);
  T=blockToeplitz(Y,L);
  [~,D,V]=svd(T,0);
  D=diag(D);
  V=V(:,1:K);
  D=D(1:K);
end

function T=blockToeplitz(Y,L)
  [P N]=size(Y);
  M=N-L+1;
  T=zeros(P*M,L);
  for b=1:P
    T(M*(b-1)+(1:M),:)=toeplitz(Y(b,L:end),Y(b,L:-1:1));
  end
end

