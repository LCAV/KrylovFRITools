function [V D termsig]=sigspace_toeplitz_Krylov(Y,K,v0)
% INPUT
%	Y :			A Px(2*M+1) array of measurements (see measurement model above)
%	K :			The dimension of the signal space
%	v0 :		(optional) The initial vector for Lanczos algorithm
%
% OUTPUT
%	V :			An orthonormal basis of the signal space
%	D :			The modulus of the associated singular values
%	termsig :	The termination signal of the Krylov decomposition (see the MATLAB documentation of EIGS)

  [P N]=size(Y);
  M=floor(N/2);
  tc=Y(:,(N-M):end).'; %dim = MxP
  tl=Y(:,(N-M):-1:1).';  %dim = MxP
  sn=mod(N,2);
  
  %% FFT implementation of T'*T*f (T is  stack-toeplitz)
  %tc contains the leading column of each toeplitz block (P columns)
  %tr contains the leading row of each toeplitz block (P rows)
  fc=fft([tc;zeros(1,P);tl(end:-1:2,:)]);
  fcconj=fft([conj(tl);zeros(1,P);conj(tc(end:-1:2,:))]);
  function g=mmult(f)
    g=[f; zeros(M+sn,1)];
    g=ifft(bsxfun(@times,fc,fft(g)));
    g((M+2):end,:,:)=0;
    g=ifft(bsxfun(@times,fcconj,fft(g)));
    g=sum(g(1:(M+1),:),2);
  end
  
  %% Arguments for eigs()
  A=@(f)mmult(f); %pointer to the mat-vec multiply
  if nargin>2
    opts.v0=v0; %initial vector of lanczos, set for tracking
  end
  opts.tol=1e-6;  %play with it
  opts.issym=1;   %We have an hermitian square matrix
  opts.isreal=0;  %may apply if Y is real
  opts.p=max(2*K+1,5);     %# of lanczos vectors
  
  [V D termsig]=eigs(A,M+1,K,'lm',opts);
  D=sqrt(diag(D));
end

