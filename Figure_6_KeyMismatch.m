%% Secrecy capacity against N or SNR
% Hiba Alwazani
% Sept 2023
clc;
clear all;
close all;

%% System Parameters
% Channels construction and path loss
T=10; %number of symbols per coherence interval
F=1; %number of frames to consider
Tk=50; % number of symbols reserved for CE
Ts=2; %switching symbol time
Ts2=8; %not optimal
fc=1e9; %carrier frequency
c=3e8;% speed of light
global lambda
 lambda =c/fc; % wavelength
Qtab=[2 4 8 16 32 64 128];
Q= 2;
%Qtab(randi([1 7],1,1))


    scenarioNum=3;
    [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F,dbe]=Scenario(scenarioNum);

T=100
Ts=2;
Tk=T/2;
Tstab=[2:2:Tk];
Ltab=Tk./Tstab;
%Tstab=fliplr(Tstab);
N=100;
M=5000;
for n=1:length(Ltab)
 Ts=Tk/Ltab(n);
 Ts
n
    % monte carlo start 
for m=1:M
      %% Network Simulation
    % Channel
rho0=0.2;
rho1=0.3;
sigma=1e-9;
    [hab, hae,hbe, har, hrb, hre, R]=channels(N,rho0,rho1, beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re);
    theta= 2*pi*rand(N,1);
    Theta=diag(exp(1i*theta)); 
    gab=hab+har'*Theta*hrb;
    z=(2*sigma)/(P*Ts);
zab=sqrt(z)*sqrt(1/2).*(randn+1i.*randn);
zba=sqrt(z)*sqrt(1/2).*(randn+1i.*randn);
ze=sqrt(z)*sqrt(1/2).*(randn+1i.*randn);
gabtilde=gab+zab;
gbatilde=gab+zba;
gbe=hae+har'*Theta*hre+ze;

hF=gabtilde;
hB=gbatilde;
hE=gbe;
theta_hatf=angle(hF);%+delthetaf;
theta_hatb=angle(hB);%+delthetab;
theta_hate=angle(hE);%+delthetab;
%% quantization
for q=1:Q+1
partition(q)=(2*pi*(q-1))/Q-pi;
end
ind(:,m)=quantiz([theta_hatf theta_hatb theta_hate],partition);

end 
     %% Eve Leak Prob
    rhoab=beta_ab+beta_ar*beta_rb*trace(R*R');
    rhoae=beta_ae+beta_ar*beta_re*trace(R*R');
    Ihabhae=-log2(1-(rho1*sqrt(beta_ab*beta_ae)+trace(R*R')*rho1*beta_ar*sqrt(beta_rb*beta_re))^2/((rhoab+z)*(rhoae+z)));
    Pe=(sqrt(2*Ihabhae)+0.25);
    %% OTP

s(n)=sum(ind(1,:)~=ind(2,:))/M;
se1(n)=sum(ind(1,:)~=ind(3,:))/M;
se2(n)=sum(ind(2,:)~=ind(3,:))/M;

%p1=s/M; %actual match
end
               


figure(1)

semilogy(Ltab, s,'b-.','linewidth',2)
hold on
semilogy(Ltab, se1,'r-','linewidth',2)
hold on
semilogy(Ltab, se2,'ko','linewidth',2)
xlabel('L intervals')%epsilon (correlation)')
ylabel('Key Mismatch Probability')
legend('Between Bob and Alice','Between Eve and Alice','Between Eve and Bob')
set(gca,'fontsize',16);
grid on

%% Secret Key Rate based on correlation
%Rskg=log2(1-)
% do for correlation and irs based
% kmr and key length with optimal Q
