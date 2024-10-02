function [hab, hae,hbe, har, hrb, hre, R]=channels(N, rho0,rho1, beta_ab,beta_ae,beta_be,beta_ar,beta_rb,beta_re);

% Channels simulations
%% RIS Configuration
fc=3e9; %carrier frequency
c=3e8;% speed of light
global lambda
 lambda =c/fc; % wavelength
 
%  Interd=lambda/4;
% 
% 
% 
% % Setup surface
% Nr = 5; %row 
% Nc = 20; %column
% dr = 0.5*lambda;% less than 0.5 
% dc = 0.5*lambda;
% N=Nr*Nc;
% init=[0, 3, 1.5];
% m=1;
% for i=1:Nc;
%     for j=1:Nr
%         u(i,j,1)=1.5+0.25*lambda*j;
%          u(i,j,3)=0+0.25*lambda*i;
%          ss(m,1)=1.5+0.25*lambda*j;
%             ss(m,3)=0+0.25*lambda*i;
%          ss(m,2)=3;
%          m=m+1;
%     end
% end
% u(:,:,2)=3; % remains 3 meters
% B=reshape(u, [100,3]);
% for i=1:N
%     for j=1:N
%         R2(i,j)=sinc((2*pi/lambda)*norm(ss(j)-ss(i)));
%     end
% end
% trace(R2*R2')
%% SINC Model
   for mm1 = 1:N
        m_z = ceil(mm1/N); %% Eq. (3)
        m_x = mod(mm1-1,N)+1; %% Eq. (3)
        for mm2 = 1:N
            n_z = ceil(mm2/N); %% Eq. (3)
            n_x = mod(mm2-1,N)+1; %% Eq. (3)
            d_temp  = sqrt(  (m_x-n_x)^2 +  (m_z-n_z) ^2 )*(lambda/4); %% Eq. (1)
            R2(mm2,mm1) = sinc(2*d_temp/lambda); %% Eq. (14)
        end
    end
%% Correlation matrix
rho0=0.3;
R=zeros(N,N);
for i=1:N
    for j=1:N
        R(i,j)=rho0^abs(i-j);
    end
end

%% Direct Channels
h=sqrt(1/2).*(randn(1)+1i.*randn(1));
htilde=sqrt(1/2).*(randn(1)+1i.*randn(1));
hab=sqrt(beta_ab).*h; %Alice to Bob
hae=sqrt(beta_ae).*(rho1*h+sqrt(1-rho1^2)*htilde); %Alice to Eve
hbe=sqrt(beta_be).*sqrt(1/2).*(randn(1)+1i.*randn(1));

%% RIS to Eve and Bob
hris=sqrt(1/2).*(randn(N,1)+1i.*randn(N,1));
h_ristilde=sqrt(1/2).*(randn(N,1)+1i.*randn(N,1));
hrb=sqrt(beta_rb).*hris;
hrb=sqrtm(R)*hrb;
hre=sqrt(beta_re).*(rho1.*hris+sqrt(1-rho1^2).*h_ristilde);
hre=sqrtm(R)*hre;
har=sqrt(beta_ar)*sqrt(1/2).*(randn(N,1)+1i.*randn(N,1));
har=sqrtm(R)*har;
end
