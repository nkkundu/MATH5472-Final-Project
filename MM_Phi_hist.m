function [Phi_opt_lmmse,iter,MMSE]= MM_Phi_hist(Phi_0,M,K,C_thth,sigma_sq,epsilon)

W=sigma_sq*eye(M*K);
C_inv=inv(C_thth);
c1=norm(C_thth,1);

idx=zeros(M,K);
for ji=1:M
   idx(ji,:)=[ji:M:ji+M*(K-1)]; 
end

Phi_t=Phi_0;
Phi_til=kron(Phi_t,eye(M));
mse_t=trace( inv(C_inv + Phi_til'*inv(W)*Phi_til) );
mse_p=mse_t;
MMSE=[];
MMSE=[MMSE mse_t];
iter=0;
while(1)
    iter=iter+1;
    ph_t=kron(Phi_t,eye(M));
    A_t=((ph_t*C_thth*ph_t') + W)\(ph_t*C_thth);
    
    lambda=c1*norm(A_t*A_t',1);
    B_t=(lambda*ph_t) -(A_t*A_t')*ph_t*C_thth + A_t*C_thth;
    
    tempB=zeros(K,K);
    for ki=1:M
        tempB = tempB+ B_t(idx(ki,:),idx(ki,:));
    end
    
    Phi_t=exp(1i*angle(tempB));
    mse_t=trace( inv(C_inv + kron(Phi_t'*Phi_t,eye(M))/(sigma_sq)) );
    MMSE=[MMSE mse_t];
    if(mse_p-mse_t<=epsilon)
       break; 
    end
    mse_p=mse_t;
    
     
end

Phi_opt_lmmse=Phi_t;

end