function [p1]=Match_actual(Q,sigma,P,channelgain,Ts)
% function to find the match probability using quantization level and noise
% finds actual monte carlo
   sigmasq=((sigma*2)/(P*channelgain*Ts));%

M=10000; %number of coefficients
ho=sqrt(1/2).*(randn(M, 1)+1i*randn(M, 1)); % common channel
delhF=sqrt(sigmasq).*sqrt(1/2).*(randn(M, 1)+1i*randn(M, 1)); % phase noise
delhB=sqrt(sigmasq).*sqrt(1/2).*(randn(M, 1)+1i*randn(M, 1)); % phase noise
hF=ho+delhF;
hB=ho+delhB;
theta_hatf=angle(hF);%+delthetaf;
theta_hatb=angle(hB);%+delthetab;

%% quantization
for q=1:Q+1
partition(q)=(2*pi*(q-1))/Q-pi;
end
for m=1:M
ind(:,m)=quantiz([theta_hatf(m) theta_hatb(m)],partition);

end
s=sum(ind(1,:)==ind(2,:));
p1=s/M; %actual match

end


