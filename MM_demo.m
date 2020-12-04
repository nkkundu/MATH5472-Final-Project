%This script implements an efficient MM Algorithm to solve
%a non-convex problem approximately in signal processing application
%By Neel Kanth KUNDU, Email: nkkundu@connect.ust.hk

clear; clc;
rng(1);
M=5; K=10; Niter=100;
snr=0;
rho=0.9;

C_thth=zeros(M*K,M*K);

for i=1:M*K
    for j=1:M*K
        C_thth(i,j)= rho^(abs(i-j));
    end
end
MMSE_hist=zeros(2000,Niter);
mmse_final=zeros(Niter,1);
epsilon=1e-6;

sigma_sq=10^(-0.1*snr);

for jk=1:Niter
    init_phi=2*pi*rand(K);
    init_phi=exp(1i*init_phi);
    ph_0=init_phi;
    [phi_opt,iterations,mmse]=MM_Phi_hist(ph_0,M,K,C_thth,sigma_sq,epsilon);
    MMSE_hist(1:iterations+1,jk)=real(mmse)';
    mmse_final(jk)=real(mmse(end));
     fprintf('%dth iteration is done...................\n',jk);
end


MMSE_hist1=MMSE_hist;
MMSE_hist1(MMSE_hist1==0) = NaN;
s1=size(MMSE_hist1);

figure(1)
ite=Niter;
for kk=1:ite
    plot(MMSE_hist1(:,kk),'-');
    hold on
end
xlabel('Iterations','Interpreter','latex')
ylabel('Objective','Interpreter','latex')
title("Objective Value vs iterations  (M= "+M+" ,K="+K+", $\epsilon$ = "+epsilon+" )"+" $\frac{1}{\sigma^2}$="+snr+" dB",'Interpreter','latex')

axes('position',[.65 .175 .25 .25])
box on
for kk=1:ite
    plot(MMSE_hist1(1:4000,kk),'-');
    hold on
end
axis tight


figure(2)
boxplot(mmse_final,'Labels','Objective Value')
h = findobj(gcf,'tag','Outliers');
set(h,'MarkerEdgeColor','k')
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
title("Distribution of the Final Objective Value (M= "+M+" ,K="+K+", $\epsilon$ = "+epsilon+" )"+" $\frac{1}{\sigma^2}$="+snr+" dB",'Interpreter','latex')
