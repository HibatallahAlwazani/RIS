%% Validation of SKR lowerbound Journal Paper
% Hiba Alwazani
% Jan 2024
clc;
clear all;
close all;
%% Secret Key Generation LOWER BOUND THEORY with IRS (EFFECT OF N)
% fix the following code by adding iid assumption should be less
% performance 
% add no RIS to check the rate
% create the actual setup from channels
%% System Parameters
% Channels construction and path loss
T=1000; %number of symbols per coherence interval
Tk=T/2; % number of symbols reserved for CE
% IRS params
Ts=2;  %must be even number
N=50;
deltaT=1e-3; %microsecond symbol interval
sigma_b= -96; %dBm

fc=3e9; %carrier frequency
c=3e8;% speed of light
global lambda
 lambda =c/fc; % wavelength
%N=100;
%% Path Loss and Variance
% beta_ab=0.03; % direct channel Alice to Bob and vice versa
% beta_ae=0.01; % direct channel Alice to Eve
% beta_be=0.01; % direct channel Bob to Eve
% % IRS
% beta_ar=0.7; % alice to rose
% beta_br=0.7; % bob to rose
% beta_re=0.7; % rose to eve
%beta= beta_ar*beta_br; % alice to rose to bob and vice versa
sigma=0.01; % noise received assumed common 
P=1; % Transmit power assumed common
sigma_bar=(2*sigma)/(Ts*P); %estimate noises
rho_ris=0.5; %% CORRELATION
rho=0.85;

% 
N_tab=[0:10:100];
 %% Network Simulation
 for s=1:length(N_tab)
     N=N_tab(s);
% 
s
%% Correlation matrix
R=zeros(N,N);
for i=1:N
    for j=1:N
        R(i,j)=rho_ris^abs(i-j);
    end
end


R2=eye(N);
monte=5000;
for m=1:monte

 scenarioNum=1;
    [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F,dbe]=Scenario(3);
   
    T=100; % has to be relative to N for the approx to hold true
%     
%     x_Bob=30;y_Bob=0;
% x_eve=30.5;y_eve=0; % at least 0.1 away
% dBE=sqrt((x_eve-x_Bob)^2+(y_eve-y_Bob)^2);
%    [rho, rho2] = bessel(dBE);%correlation bob and eve
    %rho0=0.25; %correlation ris

h=sqrt(1/2).*(randn(1)+1i.*randn(1));
htilde=sqrt(1/2).*(randn(1)+1i.*randn(1));
hab=sqrt(beta_ab).*h;
hae=sqrt(beta_ae).*(rho*h+sqrt(1-rho^2)*htilde);
%% RIS to Eve and Bob
hris=sqrt(1/2).*(randn(N,1)+1i.*randn(N,1));
h_ristilde=sqrt(1/2).*(randn(N,1)+1i.*randn(N,1));
hrb=sqrt(beta_rb).*hris;
hrb1=sqrtm(R)*hrb;
hre=sqrt(beta_re).*(rho.*hris+sqrt(1-rho^2).*h_ristilde);
hre1=sqrtm(R)*hre;
 har=sqrt(beta_ar)*sqrt(1/2).*(randn(N,1)+1i.*randn(N,1));
 har1=sqrtm(R)*har;
 
Theta=diag(ones(N,1));

 
    sigma_bar=(2*sigma)/(Ts*P); %estimate noises

    beta= beta_ar*beta_rb; % alice to rose to bob and vice versa
    %rho=rand(N+1,1); % correlations coefficient
    %% Network Simulation

   %% A

hcolb1=hab+har1'*Theta*hrb1;
hcole1=hae+har1'*Theta*hre1; 

rhoabe1(m)=hcolb1*hcole1';
rhoab1(m)=hcolb1*hcolb1';
rhoae1(m)=hcole1*hcole1';
%% uncorrelated
hcolb=hab+har'*Theta*hrb;
hcole=hae+har'*Theta*hre; 

rhoabe(m)=hcolb*hcole';
rhoab(m)=hcolb*hcolb';
rhoae(m)=hcole*hcole';

end
%% validation for cross-correlations of Alice, Bob, and Eve channels
rho_abae=real(sum(rhoabe1)/monte);
rho_abae_true=rho*(sqrt(beta_ab*beta_ae)+beta_ar*sqrt(beta_rb*beta_re)*trace(R*R'));
rho_abae_true1=rho*sqrt(beta_ab*beta_ae);


rho_ab=real(sum(rhoab1)/monte);
rho_ab_true=beta_ab+beta_ar*beta_rb*trace(R*R');

rho_ae=real(sum(rhoae1)/monte);
rho_ae_true=beta_ae+beta_ar*beta_re*trace(R*R');

%% uncorrelated
rho_abae1=real(sum(rhoabe)/monte);rho_ab1=real(sum(rhoab)/monte);rho_ae1=real(sum(rhoae)/monte);
%% Theoretical Lowerbound
%rho_ab=beta_ab +N*(beta_ar*beta_br); %correlation between alice and bob 
rho_ab_tilde=rho_ab+sigma_bar;
rho_ae_tilde=rho_ae+sigma_bar;
rho_ab_tildet=rho_ab_true+sigma_bar;
rho_ae_tildet=rho_ae_true+sigma_bar;

rho_ab_tildet1=beta_ab+sigma_bar;
rho_ae_tildet1=beta_ae+sigma_bar;
Rk_lb_noRIS_uncorr(s)=(1/(Ts/2))*log2(1+((beta_ab^2)/(sigma_bar*(2*beta_ab+sigma_bar))));
Rk_lb_noRIS_corr(s)=(1/(Ts/2))*log2(rho_ab_tildet1*(rho_ab_tildet1*rho_ae_tildet1-rho_abae_true1^2)/(rho_ae_tildet1*sigma_bar*(2*beta_ab+sigma_bar)));
Rk_lb_corr(s)=(1/(Ts/2))*log2(rho_ab_tilde*(rho_ab_tilde*rho_ae_tilde-rho_abae^2)/(rho_ae_tilde*sigma_bar*(2*rho_ab+sigma_bar)));
Rk_lb_corr_true(s)=(1/(Ts/2))*log2(rho_ab_tildet*(rho_ab_tildet*rho_ae_tildet-rho_abae_true^2)/(rho_ae_tildet*sigma_bar*(2*rho_ab_true+sigma_bar)));
Rk_lb_uncorr(s)=(1/(Ts/2))*log2(1+((rho_ab^2)/(sigma_bar*(2*rho_ab+sigma_bar))));
Rk_lb_uncorr_true(s)=(1/(Ts/2))*log2(1+((rho_ab_true^2)/(sigma_bar*(2*rho_ab_true+sigma_bar))));
%rho_ab=beta_ab +N*(beta_ar*beta_br); %correlation between alice and bob 
%% uncorrelated
rho_ab_tilde1=rho_ab1+sigma_bar;
rho_ae_tilde1=rho_ae1+sigma_bar;

Rk_lb_corr1(s)=(1/(Ts/2))*log2(rho_ab_tilde1*(rho_ab_tilde1*rho_ae_tilde1-rho_abae1^2)/(rho_ae_tilde1*sigma_bar*(2*rho_ab1+sigma_bar)));
Rk_lb_uncorr1(s)=(1/(Ts/2))*log2(1+((rho_ab1^2)/(sigma_bar*(2*rho_ab1+sigma_bar))));



%% Secrecy Leakage 
Igabgae(s)=-log2(1-rho_abae_true^2/(rho_ab_tildet*rho_ae_tildet));
%Rk_lb_collected(s)=Rk_lb_uncorr_true(s)-Igabgae(s); %validated!!!!!!
%Igabgae(s)=-log2(1-rho_abae^2/(rho_ab_tilde*rho_ae_tilde));
 end
 rholist=0.1:0.1:1
for i=1:length(rholist)
    rho=rholist(i);
    rhoabe=(sqrt(beta_ab*beta_ae)+beta_ar*sqrt(beta_rb*beta_re)*trace(R*R'));
Igabgaerho(i)=-log2(1-(rho*rhoabe)^2/(rho_ab_tilde*rho_ae_tilde));
end
figure(1)
hold on
plot(N_tab,Rk_lb_uncorr_true,'k-','LineWidth',2)
plot(N_tab,Rk_lb_uncorr,'ko','LineWidth',2)
plot(N_tab,Rk_lb_corr_true,'b-','LineWidth',2)
plot(N_tab,Rk_lb_corr,'b^','LineWidth',2)
plot(N_tab,Rk_lb_uncorr1,'r^','LineWidth',2)
plot(N_tab,Rk_lb_corr1,'r-','LineWidth',2)
plot(N_tab,Rk_lb_noRIS_uncorr,'k--','LineWidth',2)
plot(N_tab,Rk_lb_noRIS_corr,'k-.','LineWidth',2)
xlim([0 100])
ylabel('Secret key rate in bits per channel use')
xlabel('N')
legend('Uncorrelated Theory', 'Uncorrelated Monte','Correlated Theory', 'Correlated Monte','Uncorrelated Theory, $\mathbf{R}=\mathbf{I}_N$', 'Correlated Theory, $\mathbf{R}=\mathbf{I}_N$','Uncorrelated Theory, No RIS','Correlated Theory, No RIS','Interpreter','latex')
grid on
set(gca,'fontsize',16);



figure(2)
hold on
plot(N_tab,Igabgae,'k-','LineWidth',2)
ylabel('Secret Leakage in bits per seconds')
xlabel('N')
grid on;
set(gca,'fontsize',16);

figure(3)
hold on
plot(rholist,Igabgaerho,'k--','LineWidth',2)

ylabel('Secret Leakage in bits per seconds')
xlabel('Correlation Coefficient between Bob and Eve')
grid on;
set(gca,'fontsize',16);

%% Functions
d=lambda/10;
function [q, q2] = bessel(d)
global lambda
x=(2*pi*d)/lambda;
q2=besselj(0,x);
fun = @(theta) (2/pi)*cos(x*cos(theta));
q= (integral(fun,0,pi/2));
end

