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
Q= 4;
%Qtab(randi([1 7],1,1))


    scenarioNum=3;
    [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F,dbe]=Scenario(scenarioNum);
T=100
Tstab=[2:2:T/5]; % consider ts tab power tab, l tab, snr 
N=100;
for n=1:length(Tstab)
     Ts=Tstab(n)
    [rho1, rho2] = bessel(dbe);%correlation bob and eve
    rho0=0.2; %correlation ris
    sigma_bar=(2*sigma)/(Ts*P); %estimate noises
    beta= beta_ar*beta_rb; % alice to rose to bob and vice versa
    rho=0.3; % correlations coefficient
    %% Network Simulation
    % Channel
    [hab, hae,hbe, har, hrb, hre, R]=channels(N,rho0,rho1, beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re);
    %% Average Key Throughput
    channelgain=beta_ab+beta_ar*beta_rb*trace(R*R');
    snr=10*log10((P*channelgain*Ts)/(2*sigma));%dB
    F=1;
    N_quanta=ceil((snr)/6.02); %number of bits optimal -1.76 missing
    [p0 p1]=Match_probability(Q,sigma,P,channelgain,Ts); %approximate prob, true probability
   [p0NR p1NR]=Match_probability(Q,sigma,P,beta_ab,Ts); %approximate prob, true probability
      pNR=p1NR^F;
    p=p1^F;
    p2=p0^F;
   
    pp(n)=p;

    %% consider BSC modelling
      x=p0;
  if p0>=0.5
       x=1-p0;
    else
        x=p0;
    end
     
    Hb=-x*log2(x)-(1-x)*log2(1-x);
   Rk=(1-Hb)*log2(Q)/(Ts/2);   
     
       Hbnr=-p1NR*log2(p1NR)-(1-p1NR)*log2(1-p1NR);
     Rk_NR=(1-Hbnr)*log2(Q)/(Tk/2); % no switchin ts
  
    %% Numerical Optimal Q 
for q=1:length(Qtab)
    Q2=Qtab(q);
     [p0q p1q]=Match_probability(Q2,sigma,P,channelgain,Ts); %approximate prob, true probability
    pq=p1q;
     Hbq=-pq*log2(pq)-(1-pq)*log2(1-pq);
    Rk1(q)=(1-Hbq)*log2(Q2)/(Ts/2); %number of bits optimal   
end
[M, I]=max( Rk1);
Qstar=Qtab(I);
    [p0q p1q]=Match_probability(Qstar,sigma,P,channelgain,Ts); %approximate prob, true probability
    pq=p1q;
     Hbq=-pq*log2(pq)-(1-pq)*log2(1-pq);
    Rk_Qoptim=(1-Hbq)*log2(Qstar)/(Ts/2);
      
    %% Rk Theory LB
    rho_abae=rho1*(sqrt(beta_ab*beta_ae)+beta_ar*sqrt(beta_rb*beta_re)*trace(R*R'));
    rho_ab=beta_ab+beta_ar*beta_rb*trace(R*R');
    rho_ae=beta_ae+beta_ar*beta_re*trace(R*R');
    %% Theoretical Lowerbound
    sigma_bar=(2*sigma)/(Ts*P); %estimate noises
    rho_ab_tilde=rho_ab+sigma_bar;
    rho_ae_tilde=rho_ae+sigma_bar;
    Rk_lb=(1/(Ts/2))*log2(rho_ab_tilde*(rho_ab_tilde*rho_ae_tilde-rho_abae^2)/(rho_ae_tilde*sigma_bar*(2*rho_ab+sigma_bar)));
    %% Finding the right key allocation time
    RItilde=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*P/sigma);
    Tk_trial=tktrial(Rk,  RItilde,T) ; 
    L=Tk/2; 
    Tk_th=tktrial(Rk_lb,  RItilde,T) ;  Lth=  Tk_th/2;
    RItildenr=log2((beta_ab)*P/sigma); %%no RIS
    Tk_nr=tktrial(Rk_NR,   RItildenr,T) ;       Lnr=Tk_nr/2;
    Tk_trialq=tktrial(Rk_Qoptim,   RItilde,T) ;       L2q=Tk_trialq/2;
   % RIotptilde=log2(((abs(hab)+abs(hrb)'*abs(har))^2*P)/sigma); %otp paper ##  competing monte with theory something is wrong here
        RIotp=log2(1+((beta_ab+sqrt(1-rho0)*2*sqrt(pi/4)^(3)*sqrt(beta_ab*beta_ar*beta_rb)*trace(R*R')+(1-rho0)*(pi^2/16)*(beta_ar*beta_rb*trace(R*R')^2))*P)/sigma); % rayleight dist,
       % Tk_trialq2=tktrial(Rk_Qoptim2,   RItilde,T) ;       L2q2=Tk_trialq2/2;
     
        Tk_otp=tktrial(Rk_lb,   RIotp,T) ;   

 %% Information Rate
    L=Tk_trial/2;
    L2=Tk/Ts;
    RI=log2(((beta_ab+beta_ar*beta_rb*trace(R*R'))*P*log(L))/sigma); % Ts=2, tk optimal , q not 
    RIts=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*log(L)*P/sigma);% Ts=random, tk optimal , q not 
    RIq=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*log(L2q)*P/sigma);% Ts,tk,q, optim
   %  RIq2=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*log(L2q2)*P/sigma);% Ts,tk,q, optim
  RIth=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*log(Lth)*P/sigma);% Ts,tk,q, optim
  RInr=log2(1+(beta_ab*P)/sigma); % no RIS
    %% Eve Leak Prob
    rhoab=beta_ab+beta_ar*beta_rb*trace(R*R');
    rhoae=beta_ae+beta_ar*beta_re*trace(R*R');
    Ihabhae=-log2(1-(rho1*sqrt(beta_ab*beta_ae)+trace(R*R')*rho1*beta_ar*sqrt(beta_rb*beta_re))^2/((rhoab+sigma_bar)*(rhoae+sigma_bar)));
    Pe=0
    %% OTP
RIotp=log2(((abs(hab)+abs(hrb)'*abs(har))^2*P)/sigma); % their paper %*log( Lotp)

    %% Secret Transmission Rate
        RsOTP(n)=(1-Pe)*min(Tk_otp/(T)*Rk_lb,(T-Tk_otp)/T*RIotp); %% Otp
        RsNR(n)=(1-Pe)*min(Tk_nr/(T)*Rk_NR,(T-Tk_nr)/T*RInr);
    Rs(n)=(1-Pe)*min(Tk/(T)*Rk,(T-Tk)/T*RI);
        RsTH2(n)=(1-Pe)*min(Tk_th/(T)*Rk_lb,(T-Tk_th)/T*RIth); %Rith tkith
    RsTH(n)=(1-Pe)*min(Tk/(T)*Rk_lb,(T-Tk)/T*RIts); %Rith tkith
       Rs_optimm(n)=(1-Pe)*min(Tk_trial/(T)*Rk,(T-Tk_trial)/T*RI); 
         Rs_optimm2Q(n)=(1-Pe)*min(Tk_trialq/(T)*Rk_Qoptim,(T-Tk_trialq)/T*RIq); 
         
end

figure(1)
hold on
%v'Tk optimal Theory LB',
%plot(Ntab,RsOTP,'k-','LineWidth',2)
plot(Tstab,RsTH2,'k-o','LineWidth',2)
plot(Tstab,RsTH,'k--','LineWidth',2)
plot(Tstab,Rs_optimm2Q,'b-o','LineWidth',2)
%plot(Ntab,Rs_optimmQ1,'b-o','LineWidth',2)
plot(Tstab,Rs_optimm,'b--','LineWidth',2)
plot(Tstab,RsNR,'r-','LineWidth',2)
%all tk optimal
%Q optimal'
legend('SKR Th.,  optimal T_k','SKR Th., Fixed T_k=T/2','Practical SKR, Optimized','Practical SKR, T_k^*, Q=4','No RIS' )
ylabel('Secret Transmission Rate (bps)')
xlabel('N')
xlim([2 20])
grid on;
set(gca,'fontsize',16);

 
 
 

%% Functions
function [q, q2] = bessel(d)
global lambda
x=(2*pi*d)/lambda;
q2=besselj(0,x);
fun = @(theta) (2/pi)*cos(x*cos(theta));
q= (integral(fun,0,pi/2));
end



function Tk1=tktrial(Rk,RI,T) %% may need F here F not 1 is wrong here

    syms Tk

    eqn = Tk*(   Rk+    RI)+(Tk-T)*log2(log((Tk/2)))==(T* RI);
    S = vpasolve(eqn,Tk);
    Tk1=double(S);

end
function Qstar2=optimQ2(sigmasq) %% may need F here F not 1 is wrong here
syms Q

pdiff1=(-(1/sigmasq)*pi*sec(pi/Q)^2*tan(pi/Q))/(Q^2*((1/sigmasq)*tan(pi/Q)^2+1)^2);
dt=2*pi/Q;
fun= @(theta) tan(theta).^2./(tan(theta).^2+sigmasq);
int1=int(fun, 0,pi/Q);
pdiff2= (tan(pi/Q)^2/(tan(pi/Q)^2+sigmasq))*(-1/(2*Q))+(1/(2*pi))*int1;
%% find optimal q by setting equal to 0
pdiff=(pdiff1+pdiff2);


p=(0.5*(tan(dt/2)^2/(tan(dt/2)^2+sigmasq))+(1/dt)*int1);%+(4/(pi^2*dt))*int2;

  
  Hbdiff=-pdiff*log2(p)-pdiff/(log(2)) +(-pdiff)*log2(1-p)+pdiff/log(2);%p*((pdiff)/(p*log(2))
    Hb=-p*log2(p)-(1-p)*log2(1-p);

   
Rprime=1/(Q*log(2))*(1-Hb)+log2(Q)*Hbdiff;
eqn= Rprime==0; % finding optimal Q by setting obj to 0
Qs=vpasolve(eqn);
Qstar2=abs(double((Qs)));

end
 
function Qstar=optimQ(sigmasq) %% may need F here F not 1 is wrong here
syms Q

pdiff1=(-(1/sigmasq)*pi*sec(pi/Q)^2*tan(pi/Q))/(Q^2*((1/sigmasq)*tan(pi/Q)^2+1)^2);
dt=2*pi/Q;
fun= @(theta) tan(theta).^2./(tan(theta).^2+sigmasq);
int1=int(fun, 0,pi/Q);
pdiff2= (tan(pi/Q)^2/(tan(pi/Q)^2+sigmasq))*(-1/(2*Q))+(1/(2*pi))*int1;
%% find optimal q by setting equal to 0
pdiff=(pdiff1+pdiff2);
p=(0.5*(tan(dt/2)^2/(tan(dt/2)^2+sigmasq))+(1/dt)*int1);%+(4/(pi^2*dt))*int2;
Rprime=1/(Q*log(2))*p+log2(Q)*pdiff;
eqn= Rprime==0; % finding optimal Q by setting obj to 0
Qs=solve(eqn);
Qstar=abs(double(vpa(Qs)));

end
