% Ergodic information rate
% check and validate scaling law
clear all;
clc;
N=100;
rho=0.1;

Te=50;
sigma=1
l=50; %Watts
powertab=[1:1:10];
for s=1:length(powertab)  
    s
    scenarioNum=1;
    [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F,dbe]=Scenario(scenarioNum);
    T=100; % has to be relative to N for the approx to hold true

    %[rho1, rho2] = bessel(dbe);%correlation bob and eve
    rho1=0;
    rho0=0.25; %correlation ris
    P=powertab(s);  monte=1000;
 
L=10;
L2=50;
L3=100;
h2=0;
for m=1:monte
    
    [hab, hae,hbe, har, hrb, hre, R]=channels(N,rho0,rho1, beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re);
    channelgain=beta_ab+beta_ar*beta_rb*trace(R*R');

   habexp(m)=(abs(hab).^2);
  maxh=switching(L,hab,har,hrb,N);
    maxh2=switching(L2,hab,har,hrb,N);
    maxh3=switching(L3,hab,har,hrb,N);
   hexp(m)=maxh;

    h2=h2+(abs(hab)+abs(hrb)'*abs(har))^2*P/sigma;
    rate_opt(m)=log2(1+(abs(hab)+abs(hrb)'*abs(har))^2*P/sigma);%maximum choice
ratemax(m)=log2(1+maxh*P/sigma);%maximum choice
ratemax2(m)=log2(1+maxh2*P/sigma);%maximum choice
ratemax3(m)=log2(1+maxh3*P/sigma);%maximum choice
ratenoRIS(m)=log2(1+abs(hab).^2*P/sigma); %simulatedno
end
 hexpmon=sum(hexp)/monte;
 hexpab=sum(habexp)/monte;
 h2=h2/monte;
ratemonUB(s)=log2(1+hexpmon*P/sigma); %expectation of max inside

ratemon(s)=sum(ratemax)/monte; %expectation outisde
ratemon22(s)=sum(ratemax2)/monte; %expectation outisde
ratemon33(s)=sum(ratemax3)/monte; %expectation outisde

rateopt(s)=sum(rate_opt)/monte;
ratemon2(s)=sum(ratenoRIS)/monte; %actual 
ratethlb(s)=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*log(L)*P/sigma); %theoretical scaling
ratetNORISth(s)=log2(1+(beta_ab*P)/(sigma)); % because sqrt(beta_ab).*sqrt(1/2).
ratetNORISmonte(s)=log2(1+( hexpab*P)/(sigma));
RIotp2(s)=log2(1+h2); % best ris

end
figure(1)
hold on
plot(powertab,rateopt,'k-', 'LineWidth',2)
plot(powertab,ratethlb,'k-.', 'LineWidth',2)
%plot(powertab,ratemonUB,'b-', 'LineWidth',2)
plot(powertab,ratemon,'bo', 'LineWidth',2)
plot(powertab,ratemon22,'b--', 'LineWidth',2)
plot(powertab,ratemon33,'b-^', 'LineWidth',2)


%plot(powertab,ratetNORISth,'r-', 'LineWidth',2)
plot(powertab,ratemon2,'r-', 'LineWidth',2)

%plot(powertab,ratetNORISmonte,'ko', 'LineWidth',2)
xlabel('P in (Watts)')
ylabel('Ergodic Information Rate in (bits/ channel use)')
legend('Optimal RIS, Th.','Max Choice, Th. L=10','Max Choice, Sim. L=10','Max Choice, Sim. L=50','Max Choice, Sim. L=100', 'Direct Channel, Sim.')
grid on
xlim([1 10])
set(gca,'fontsize',16);

function maxh=switching(L,hab,har,hrb,N)
    for t=1:L
    theta= 2*pi*rand(N,1);
    Theta=diag(exp(1i*theta)); 
    hphase(t)=hab+har'*Theta*hrb;
    end
    habs=abs(hphase).^2;
    maxh=max( habs);
end

%%%%%%%%%%%%%%%%%%%%%%0
% clear all;
% clc;
% N=100;
% 
% Ntab=[10 50 100 200 300]
% rho=0.1;
% for n=1:length(Ntab)
%     N=Ntab(n)
% %% Correlation matrix
% R=zeros(N,N);
% for i=1:N
%     for j=1:N
%         R(i,j)=rho^abs(i-j);
%     end
% end
% 
% beta_ab=1; % direct channel Alice to Bob and vice versa
% beta_ar=1; % direct channel Alice to Eve
% beta_rb=0.7; % direct channel Bob to Eve
% Te=50;
% sigma=1
% l=50; %Watts
% ratemonUB=[];
% ratemon=[];
% ratethlb=[];
% l=100; 
% Tetab=[10 20 30 40 50 70 90 120 150 200 500]
% sigma=1
% P=1
% for s=1:length(Tetab)
%     s
%    Te=Tetab(s);
% monte=100;
% 
% 
%     scenarioNum=1;
%     [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F,dbe]=Scenario(scenarioNum);
%     T=100; % has to be relative to N for the approx to hold true
% 
%     %[rho1, rho2] = bessel(dbe);%correlation bob and eve
%     rho1=0;
%     rho0=0.25; %correlation ris
% 
% for m=1:monte
%     
%     [hab, hae,hbe, har, hrb, hre, R]=channels(N,rho0,rho1, beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re);
%     channelgain=beta_ab+beta_ar*beta_rb*trace(R*R');
%    for t=1:Te
%    theta= 2*pi*rand(N,1);
%     Theta=diag(exp(1i*theta)); 
% hphase(t)=hab+har'*Theta*hrb;
% 
%    end
%    habexp(m)=(abs(hab).^2);
%    habs=abs(hphase).^2;
%    hexp(m)=max( habs);
%   h= max(habs);
% rate(m)=log2(1+h*P/sigma);
% end
%  hexpmon=sum(hexp)/monte;
%  hexpab=sum(habexp)/monte;
% ratemonUB(s,n)=log2(1+hexpmon*P/sigma);
% ratemon(s,n)=sum(rate)/monte;
% ratethlb(s,n)=log2((beta_ab+beta_ar*beta_rb*trace(R*R'))*log(Te)*P/sigma);
% 
% %ratethub(s)=log2(1+(beta_ab+beta_ar*beta_rb*trace(R*R'))*P*(log(Te)+log(log(Te)))/sigma); % taken from paper
% 
% end
% 
% figure(2)
% hold on
% 
% plot(Tetab,ratethlb(:,n),'r-.', 'LineWidth',2)
% %plot(Tetab,ratethub,'b-', 'LineWidth',2)
% plot(Tetab,ratemon(:,n),'k--', 'LineWidth',2)
% plot(Tetab,ratemonUB(:,n),'b-', 'LineWidth',2)
% 
% xlabel('Switching Rate T_e')
% ylabel('Average Rate')
% legend('Theory Scaling Law','Simulated','Simulated UB')
% 
% grid on
% end
