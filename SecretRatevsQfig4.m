%% Secrecy capacity Optimization for Q, Tk, Ts
% Hiba Alwazani
% Sept 2023
clc;
clear all;
close all;
%% Secret Key Capacity with IRS
% Incorporate correlation model for the channels. done
% Obtain Rk average key throughput real and approximate value. done
% Obtain RI information rate real and scaling law.
% Pe obtain probability of attack.
% Initialize all parameters put down objectives.
%% System Parameters
T=100; %number of symbols per coherence interval
Ts=2;
Ptab=[2 5 10];
for tt=1:length(Ptab)
    P=Ptab(tt);
    tt
Tktab=[2:2:T];
Qtab=[4 8 16 32 64 128];
for q=1:length(Qtab)
    Q=Qtab(q);
for t=1:length(Tktab)
    Tk=Tktab(t);
N=10;
    scenarioNum=3;
    [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F]=Scenario(scenarioNum);
    sigma_bar=(2*sigma)/(Ts*P); %estimate noises
F=1;%####
    beta= beta_ar*beta_rb; % alice to rose to bob and vice versa
    rho=rand(N+1,1); % correlations coefficient
    %% Network Simulation
    % Channels
    rho0=0.5; %correlation ris
    rho1=0.42; %correlation bob and eve
    [hab, hae,hbe, har, hrb, hre, R]=channels(N,rho0,rho1, beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re);
    %% Average Key Throughput

    channelgain=beta_ab+beta_ar*beta_rb*trace(R*R');
% Average Key Throughput
P=Ptab(tt);
[p0 p1]=Match_probability(Q,sigma,P,channelgain,Ts); %approximate prob, true probability
 x=p0;
  if p0>=0.5
       x=1-p0;
    else
        x=p0;
    end
     
    Hb=-x*log2(x)-(1-x)*log2(1-x);
   Rk=(1-Hb)*log2(Q)/(Ts/2);
%% Information Rate
L=Tk/Ts;
    for tl=1:L
    theta= 2*pi*rand(N,1);
    Theta=diag(exp(1i*theta)); 
    hphase(tl)=hab+har'*Theta*hrb;
    end
    habs=abs(hphase).^2;
    maxh=max( habs);
RI=log2(channelgain*log(L)*P/sigma); %theory
RItilde=  log2(1+(maxh*P)/sigma); % actual
%% Eve Leak Prob
rhoab=beta_ab+beta_ar*beta_rb*trace(R*R');
rhoae=beta_ae+beta_ar*beta_re*trace(R*R');
Ihabhae=-log2(1-(rho1*sqrt(beta_ab*beta_ae)+trace(R*R')*rho1*beta_ar*sqrt(beta_rb*beta_re))^2/((rhoab+sigma_bar)*(rhoae+sigma_bar)));
Pe=(2^(-2*log2(Q))+sqrt(2*Ihabhae))^(F);

%% Secret Transmission Rate
Rs(q,t)=(1-Pe)*min(Tk/(T)*Rk,(T-Tk)/T*RI);
Rsbar(q)=(1-Pe)*min(Tk/(T)*Rk,(T-Tk)/T*RI);


    end
%% Intersect Tk
Ts=2;
Rk=(1-Hb)*log2(Q);
Tk_intersect=(T*RItilde)/(Rk+RItilde); %intersection point.

%% trial and error method
syms Tk1
RItk=log2(channelgain*log(Tk1/2)*P/sigma);
eqn = Tk1*(Rk+RItk)==(T*RItk);
S = vpasolve(eqn,Tk1);
Tk_trial=double(S);
RI2=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*log((Tk_trial/Ts))*P/sigma);
Rs_Tktrial(tt,q)=(1-Pe)*min(Tk_trial/(T)*Rk,(T-Tk_trial)/T*RI2);
Rs_Tk(tt,q)=(1-Pe)*min(Tk_intersect/(T)*Rk,(T-Tk_intersect)/T*RItilde);

end

end
[M, I]=max(Rs(:));
[I1,I2] = ind2sub(size(Rs),I);
Qstar=Qtab(I1);
Tkstar=Tktab(I2);
[X Y]=meshgrid(Qtab,Tktab);
figure(1)
surf(X,Y,Rs')
zlabel('Secret Transmission Rate (bps)')
xlabel('Q')
ylabel('Tk')
colorbar
grid on
set(gca,'fontsize',16);
figure(2)
hold on
plot(Qtab,Rs_Tktrial(1,:),'k-','LineWidth',2)
plot(Qtab,Rs_Tktrial(2,:),'b-','LineWidth',2)
plot(Qtab,Rs_Tktrial(3,:),'r-','LineWidth',2)
plot(Qtab,Rs_Tk(1,:),'k-.','LineWidth',2)
plot(Qtab,Rs_Tk(2,:),'b-.','LineWidth',2)
plot(Qtab,Rs_Tk(3,:),'r-.','LineWidth',2)
ylabel('Secret Transmission Rate for Tk optimized (bps)')
xlabel('Q')
legend('P=2','P=5', 'P=10')
grid on
set(gca,'fontsize',16);
