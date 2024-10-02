%% Simulation scenarios based on practical settings
function [beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re,P, sigma, T,F,dBE]=simScenario(scenarioNum);

%% Scenario 1: RIS's effect is amplified
if scenarioNum==1
%           #Rose
%  *Alice    ^Eve *Bob
T=100; %number of symbols per coherence interval
F=1; %number of frames to consider
Tk=(T-50)/2; % number of symbols reserved for CE
Ts=2;  
N=100;
sigma_b= -96; %dBm
P=25;% Transmit power dBm assumed common
P=10^((P-30)/10); % watt
sigma=10^((sigma_b-30)/10); % % noise received assumed common 


%% Physical Set-up
%% Distance and Simulation Set-up
%  *Alice ^Eve           *Bob

x_IRS= 0; y_IRS=3;
x_alice=0;y_alice=0;
x_Bob=30;y_Bob=0;
x_eve=30;y_eve=1; % at least 0.1 away

%% random eve rotating arond bob
x_eve=30.5+(1 + 2*rand);y_eve=0; % at least 0.1 away
% initial Eve positioning

dAB=abs(x_Bob-x_alice);
d1=y_IRS;
d2=x_Bob-x_IRS;
dRB=sqrt(d1^2+d2^2);
dAR=sqrt(x_IRS^2+y_IRS^2);
dBE=sqrt((x_eve-x_Bob)^2+(y_eve-y_Bob)^2);
dAE=sqrt((x_eve-x_alice)^2+(y_eve-y_alice)^2);
dER=sqrt((x_eve-x_IRS)^2+(y_eve-y_IRS)^2);% Eves are randomly deployed according to a PPP within 1m of Alice


% 
% figure
% hold on
% d =+0.5;
% scatter(x_alice,y_alice,'kx', 'LineWidth',2);
% text(x_alice+d,y_alice+d, 'Alice');
% scatter(x_IRS,y_IRS, 'b^', 'LineWidth',2);
% text(x_IRS+d,y_IRS+d, 'RIS');
% scatter(x_Bob,y_Bob,'bx', 'LineWidth',2);
% text(x_Bob+d,y_Bob+d, 'Bob');
% scatter(x_eve,y_eve,'ro', 'LineWidth',2);
% text(x_eve+d,y_eve+d, 'Eve');
% legend('Alice', 'Rose','Bob','Eve')
% xlabel('distance (m)','Interpreter', 'Latex')
% ylabel('distance (m)', 'Interpreter', 'Latex')
% grid on;
% set(gca,'fontsize',16);

% PATH LOSS EXPONENTS
eta_AB=3.67;
eta_AE=3.67;
eta_EB=3.67;
eta_AR=2.2;
eta_ER=2.2;
eta_RB=2.2;
PL_0=30; %dB
Penetration=0;%dB 20
PLab=PL_0+10*eta_AB*log10(dAB)+Penetration;
PLar=PL_0+10*eta_AR*log10(dAR);
PLrb=PL_0+10*eta_RB*log10(dRB);
PLeb=PL_0+10*eta_EB*log10(dBE);
PLer=PL_0+10*eta_ER*log10(dER);
PLae=PL_0+10*eta_AE*log10(dAE)+Penetration; %%% ::::: note eta_AE unsure could be wrong not mentioned in paper!!!!
beta_ab=10^(-PLab/10); 
beta_ar=10^(-PLar/10);
beta_rb=10^(-PLrb/10);
beta_be=10^(-PLeb/10);
beta_re=10^(-PLer/10);
beta_ae=10^(-PLae/10);

elseif scenarioNum==2
%% Scenario 2: RIS effect minimal no penetration
%           #Rose
%  *Alice    ^Eve *Bob
T=1000; %number of symbols per coherence interval
F=1; %number of frames to consider
Tk=(T-200)/2; % number of symbols reserved for CE
Ts=2;  
N=100;
sigma_b= -96; %dBm
P=20;% Transmit power dBm assumed common
P=10^((P-30)/10); % watt
sigma=10^((sigma_b-30)/10); % % noise received assumed common 


%% Physical Set-up
%% Distance and Simulation Set-up
%  *Alice ^Eve           *Bob

x_IRS= 25; y_IRS=5;
x_alice=0;y_alice=0;
x_Bob=30;y_Bob=0;
x_eve=0.5*rand+30+0.5;y_eve=0.5*rand+0.5; % at least 0.1 away
% initial Eve positioning

dAB=abs(x_Bob-x_alice);
d1=y_IRS;
d2=x_Bob-x_IRS;
dRB=sqrt(d1^2+d2^2);
dAR=sqrt(x_IRS^2+y_IRS^2);
dBE=sqrt((x_eve-x_Bob)^2+(y_eve-y_Bob)^2);
dAE=sqrt((x_eve-x_alice)^2+(y_eve-y_alice)^2);
dER=sqrt((x_eve-x_IRS)^2+(y_eve-y_IRS)^2);% Eves are randomly deployed according to a PPP within 1m of Alice


% 
% figure
% hold on
% d =+0.5;
% scatter(x_alice,y_alice,'kx', 'LineWidth',2);
% text(x_alice+d,y_alice+d, 'Alice');
% scatter(x_IRS,y_IRS, 'b^', 'LineWidth',2);
% text(x_IRS+d,y_IRS+d, 'RIS');
% scatter(x_Bob,y_Bob,'bx', 'LineWidth',2);
% text(x_Bob+d,y_Bob+d, 'Bob');
% scatter(x_eve,y_eve,'ro', 'LineWidth',2);
% text(x_eve+d,y_eve+d, 'Eve');
% legend('Alice', 'Rose','Bob','Eve')
% xlabel('distance (m)','Interpreter', 'Latex')
% ylabel('distance (m)', 'Interpreter', 'Latex')
% grid on;
% set(gca,'fontsize',16);

% PATH LOSS EXPONENTS
eta_AB=4;
eta_EB=3.5;
eta_AR=2;
eta_ER=2.2;
eta_RB=2;
PL_0=30; %dB
Penetration=20;%dB
dAB=10;
PLnew=32.4+20*log10(28)+17.3*log10(dAB);
PLnew2=32.4+20*log10(28)+31.9*log10(dAB);

PLnew=10^(-PLnew/10); 
PLnew2=10^(-PLnew2/10); 
PLab=PL_0+10*eta_AB*log10(dAB);
PLar=PL_0+10*eta_AR*log10(dAR);
PLrb=PL_0+10*eta_RB*log10(dRB);
PLeb=PL_0+10*eta_EB*log10(dBE);
PLer=PL_0+10*eta_ER*log10(dER);
PLae=PL_0+10*eta_ER*log10(dAE); %%% ::::: note eta_AE unsure could be wrong not mentioned in paper!!!!
beta_ab=10^(-PLab/10); 
beta_ar=10^(-PLar/10);
beta_rb=10^(-PLrb/10);
beta_be=10^(-PLeb/10);
beta_re=10^(-PLer/10);
beta_ae=10^(-PLae/10);

else
%  *Alice    ^Eve *Bob
T=1000; %number of symbols per coherence interval
F=1; %number of frames to consider
Tk=(T-50)/2; % number of symbols reserved for CE
Ts=2;  
N=100;
sigma_b= -96; %dBm
P=25;% Transmit power dBm assumed common
P=10^((P-30)/10); % watt
sigma=10^((sigma_b-30)/10); % % noise received assumed common 


%% Physical Set-up
%% Distance and Simulation Set-up
%  *Alice ^Eve           *Bob

x_IRS= 0; y_IRS=3;
x_alice=0;y_alice=0;
x_Bob=30;y_Bob=0;
x_eve=30;y_eve=1; % at least 0.1 away

%% random eve rotating arond bob
x_eve=30.5+(1 + 2*rand);y_eve=0; % at least 0.1 away
% initial Eve positioning
dAB=abs(x_Bob-x_alice);
d1=y_IRS;
d2=x_Bob-x_IRS;
dRB=sqrt(d1^2+d2^2);
dAR=sqrt(x_IRS^2+y_IRS^2);
dBE=sqrt((x_eve-x_Bob)^2+(y_eve-y_Bob)^2);
dAE=sqrt((x_eve-x_alice)^2+(y_eve-y_alice)^2);
dER=sqrt((x_eve-x_IRS)^2+(y_eve-y_IRS)^2);% Eves are randomly deployed according to a PPP within 1m of Alice



fc=3; %Ghz
% PATH LOSS EXPONENTS
eta_AB=3.67; %nlos
eta_AE=3.67; %nlos
eta_EB=3.67; %los
eta_AR=2.2;  %los
eta_ER=2.2;  %los
eta_RB=2.2;  %los
PL_0=30;%32.4+20*log10(3); %dB loss
PLnewlos=17.3*log10(dAB);
PLnewnlos=31.9*log10(dAB);

Penetration=0;%dB 20
PLab=PL_0+10*eta_AB*log10(dAB)+Penetration;
PLar=PL_0+10*eta_AR*log10(dAR);
PLrb=PL_0+10*eta_RB*log10(dRB);
PLeb=PL_0+10*eta_EB*log10(dBE);
PLer=PL_0+10*eta_ER*log10(dER);
PLae=PL_0+10*eta_AE*log10(dAE)+Penetration; %%% ::::: note eta_AE unsure could be wrong not mentioned in paper!!!!
beta_ab=10^(-PLab/10); 
beta_ar=10^(-PLar/10);
beta_rb=10^(-PLrb/10);
beta_be=10^(-PLeb/10);
beta_re=10^(-PLer/10);
beta_ae=10^(-PLae/10);
end
end
