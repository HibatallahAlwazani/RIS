%% Should we take E[p] or E[SNR] wrt channels and RIS?
%close all;
clc;clear all;
 scenarioNum=3;
%% System Parameters
% Channels construction and path loss
T=60; %number of symbols per coherence interval
Tk=60; % number of symbols reserved for CE
Ts=2; %switching symbol time
fc=3e9; %carrier frequency
c=3e8;% speed of light
global lambda
 lambda =c/fc; % wavelength
 N=100; % RIS elements
 [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F,dbe]=Scenario(scenarioNum);
 [rho1, rho2] = bessel(dbe);%correlation bob and eve
 rho0=0.25; %correlation inter-ris elements
P=0.01;

 sigma_bar=(2*sigma)/(Ts*P); %estimate noises
 beta= beta_ar*beta_rb; % alice to rose to bob and vice versa
 
L=Tk/Ts; % number of keys 
F=300; %number of frames to consider

Qtab=[2,8,16];
snrdb=[85,90,95,100,105,110];
snr=10.^(snrdb./10);
for i=1:length(Qtab)
    Q=Qtab(i)
%Ntab=[1:5:100];
N=100;

%for n=1:length(Ntab)
for n=1:length(snrdb)
    n
    P=snr(n);
  sigma=1;
 sigma_bar=(2*sigma)/(Ts*P); %estimate noises
    m=1;
          Esnr=0; 
 for f=1:F % over all frames %number of key bits
     [hab, hae,hbe, har, hrb, hre, R]=channels(N,rho0,rho1, beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re);
     channelgain=beta_ab+beta_ar*beta_rb*trace(R*R');
         

      for l=1:L %over all switching periods
           theta= 2*pi*rand(N,1); v=exp(1i*theta);  Theta=diag(v);
Esnr= Esnr+abs(hab+ har'*Theta*hrb)^2;
% Noise samples
za=sqrt(sigma_bar).*sqrt(1/2).*(randn(1)+1i.*randn(1));
zb=sqrt(sigma_bar).*sqrt(1/2).*(randn(1)+1i.*randn(1));
% Collective channels
gab= hab+ har'*Theta*hrb; gba=gab;
hB=gab+zb;
hA=gba+za;
%% actual match probability
theta_hatA=angle(hA);%alice
theta_hatB=angle(hB);%bob;
% quantization
for q=1:Q+1
partition(q)=(2*pi*(q-1))/Q-pi;
end
ind(:,m)=quantiz([theta_hatA theta_hatB],partition);
 m=m+1; % index increase
 delta(:,f,l)=quantiz([theta_hatA theta_hatB],partition);
      end

 end
 for l=1:L
     s2=sum(delta(1,:,l)==delta(2,:,l));
     p4(l,n)=s2/F;
 end
        s=sum(ind(1,:)==ind(2,:));
p1(n)=s/m;
 Esnr1= Esnr/m;
 sigmasq=(2*sigma)/( Esnr1*P*Ts);% 
 p3(n)=Match_probability(Q,sigmasq);
sigmasq2=(sigma*2)/(P*channelgain*Ts) ;
 p2(n)=Match_probability(Q,sigmasq2);
 
 %% Entropy
 [Rk1(n) Rkorig(n)]= Rkentropy(p1(n),Q,Ts);
  Rk2(n) = Rkentropy(p2(n),Q,Ts);
   Rk3(n) = Rkentropy(p3(n),Q,Ts);

end
p1q(:,i)=p1;
p2q(:,i)=p2;
end
figure(1)
hold on
plot(snrdb, p1q(:,1),'k-', 'LineWidth',2)
plot(snrdb, p2q(:,1),'k-o','LineWidth',2)
plot(snrdb, p1q(:,2) ,'b--','LineWidth',2)
plot(snrdb, p2q(:,2),'b-o','LineWidth',2)
plot(snrdb, p1q(:,3) ,'r--','LineWidth',2)
plot(snrdb, p2q(:,3),'r-o','LineWidth',2)
%plot(Ntab, p3,'r-.','LineWidth',2)
%plot(Ntab, p4(1,:), 'LineWidth',2)
%plot(Ntab, p4(2,:), 'LineWidth',2)
legend('Actual Q=2', 'Approx. Q=2','Actual Q=8', 'Approx. Q=8','Actual Q=16', 'Approx. Q=16')
 ylabel('Match Probability p_0 ')
xlabel('SNR [dB]')
 grid on;
set(gca,'fontsize',16);


figure(1)
hold on
plot(SINR,p(:,1),'k-','LineWidth',2)
plot(SINR,avg,'k-o','LineWidth',2)

plot(SINR,p(:,2),'b--','LineWidth',2)
plot(SINR,avg4,'b-o','LineWidth',2)

plot(SINR,p(:,3),'r--','LineWidth',2)
plot(SINR,avg8,'r-o','LineWidth',2)

plot(SINR,p(:,4),'m.','LineWidth',2)
plot(SINR,avg16,'m-o','LineWidth',2)
legend('Q=2, th.','Q=2, Sim.','Q=4, th.','Q=4, Sim.','Q=8, th.','Q=8, Sim.','Q=16, th.','Q=16, Sim.','Location','NorthEast')
ylabel('p(SINR,Q)')
xlabel('SINR (dB)')
grid on
set(gca,'fontsize',16);
p=p(2:7,:);
%its not the collective expectation of the SNR
%its the expectation with respect to theta then the expectation with
%respect to the channels.

% figure(2)
% hold on
% plot(Ntab, Rk1,'b--','LineWidth',2)
% %plot(Ntab, Rkorig,'b-','LineWidth',2)
% plot(Ntab, Rk2,'k-','LineWidth',2)
% plot(Ntab, Rk3,'r-.','LineWidth',2)
% legend('Actual Match Probability', 'Approximation Using Expected SNR', 'Approximation Using Theoretical SNR')
% 
% ylabel('Practical Key Rate (bps)')
% xlabel('N')
% grid on;
% set(gca,'fontsize',16);

% 
% Q=128;
% Ts=2;
%   Rkentropy=[];   Rk=[]; Rk_lb=[];
% for t=1:length(Ptab)
%    t
%     channelgain=beta_ab+beta_ar*beta_rb*trace(R*R');
%     P=Ptab(t);
%   [p0 p1]=Match_probability(Q,sigma,P,channelgain,Ts); %approximate prob, true probability 
%     %% consider BSC modelling
%     if p1>0.5
%        x=1-p1;
%     else
%         x=p1;
%     end
% 
%     Hb=-x*log2(x)-(1-x)*log2(1-x);
%     Rkentropy(t)=(1-Hb)*log2(Q)/(Ts/2); %%1-Hb
%     Rk(t)=p1*log2(Q)/(Ts/2);
%     %% Rk Theory LB
%     rho_abae=rho1*(sqrt(beta_ab*beta_ae)+beta_ar*sqrt(beta_rb*beta_re)*trace(R*R'));
%     rho_ab=beta_ab+beta_ar*beta_rb*trace(R*R');
%     rho_ae=beta_ae+beta_ar*beta_re*trace(R*R');
% 
%     %% Theoretical Lowerbound
%     sigma_bar=(2*sigma)/(Ts*P); %estimate noises
%     rho_ab_tilde=rho_ab+sigma_bar;
%     rho_ae_tilde=rho_ae+sigma_bar;
%     Rk_lb(t)=(1/(Ts/2))*log2(rho_ab_tilde*(rho_ab_tilde*rho_ae_tilde-rho_abae^2)/(rho_ae_tilde*sigma_bar*(2*rho_ab+sigma_bar)));
% end
% figure(2)
% hold on
% plot(Ptab,   Rk_lb,'k-','LineWidth',2)
% plot(Ptab, Rk,'r-','LineWidth',2)
% plot(Ptab, Rkentropy,'b--','LineWidth',2)
% 
% legend('Theory LB','Rk with p','Rk entropy')
% ylabel('Secret Transmission Rate (bps)')
% xlabel('Power  (Watt)')
% grid on;
% set(gca,'fontsize',16);
% % 
%  Tk=20;
% Tstab=[2:2:Tk];
% Q=64;
%  for ts=1:length(Tstab)
%      Ts=Tstab(ts)
%  sigma_bar=(2*sigma)/(Ts*P); %estimate noises
%  [p0 p1]=Match_probability(Q,sigma,P,channelgain,Ts); %approximate prob, true probability
% x=p1;
% Rk1(ts)=p1*log2(Q)/(Ts/2);
% Hb=-x*log2(x)-(1-x)*log2(1-x)
% 
%  Rk2(ts)=(1-Hb)*log2(Q)/(Ts/2);
%  end
% 
% figure(3)
% hold on
% plot(Tstab,   Rk2,'k-','LineWidth',2)
% plot(Tstab,   Rk1,'k--','LineWidth',2)
% 
% legend('Rk entropy','Rk')
% ylabel('Key Rate (bps)')
% xlabel('Ts')
% grid on;
% set(gca,'fontsize',16);
% 
function p=Match_probability(Q,sigmasq)
% function to find the match probability using quantization level and noise
% finds actual monte carlo and approximation
%sigmasq=((sigma*2)/(P*channelgain*Ts));%
%sigmasq=0.0001;
dt=2*pi/Q;

  fun= @(theta) tan(theta).^2./(tan(theta).^2+sigmasq);
     int1=integral(fun, 0,dt/2);
     
 fun2= @(theta,x) ((2.*x.*sigmasq)./(x.^2+sigmasq).^2).*(asin(tan(theta)./x).^2);
xmin = @(theta) tan(theta);
 int2=integral2(fun2, 0,dt/2,xmin,tan(dt/2));%;;

     
p=0.5*(tan(dt/2)^2/(tan(dt/2)^2+sigmasq))+(1/dt)*int1+(4/(pi^2*dt))*int2; %approximation for high SNR
end
function [q, q2] = bessel(d)
global lambda
x=(2*pi*d)/lambda;
q2=besselj(0,x);
fun = @(theta) (2/pi)*cos(x*cos(theta));
q= (integral(fun,0,pi/2));
end
function [Rk Rkorig] = Rkentropy(p,Q,Ts)
    Rkorig=p*log2(Q)/(Ts/2);
if p>0.5
       x=1-p;
    else
        x=p;
end
    Hb=-x*log2(x)-(1-x)*log2(1-x);
    Rk=(1-Hb)*log2(Q)/(Ts/2); %%1-Hb

end