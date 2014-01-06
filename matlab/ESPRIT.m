function z=ESPRIT(Y,K,krylov,TLS)
% INPUT
%	Y :				A Px(2*M+1) array of measurements (see measurement model above)
%	K :				The dimension of the signal space
%	krylov, TLS :	Boolean flags enabling/disabling the use of the Krylov method and TLS
%
% OUTPUT
%	z :				A vector containing the phasors exp(D*omega)

  if krylov
      [V D]=sigspace_toeplitz_Krylov(Y,K);
  else
      [V D]=sigspace_toeplitz(Y,K);
  end

  if TLS
      %ESPRIT-TLS
      [~,~,C]=svd([V(1:(end-1),1:K) V(2:end,1:K)],0);
      Phi=eig(-C(1:K,(K+1):(2*K))/C((K+1):(2*K),(K+1):(2*K)));%literally: eig(Psi)
  else
      %ESPRIT-LS
      Phi=eig(V(1:(end-1),1:K)\V(2:end,1:K));
  end
      z=Phi(:);
end

