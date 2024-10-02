function [p p1]=Match_probability(Q,sigma,P,channelgain,Ts)
% function to find the match probability using quantization level and noise
% finds actual monte carlo and approximation

sigmasq=((sigma*2)/(P*channelgain*Ts));%
%sigmasq=0.0001;
dt=2*pi/Q;

  fun= @(theta) tan(theta).^2./(tan(theta).^2+sigmasq);
     int1=integral(fun, 0,dt/2);
     
 fun2= @(theta,x) ((2.*x.*sigmasq)./(x.^2+sigmasq).^2).*(asin(tan(theta)./x).^2);
xmin = @(theta) tan(theta);
 int2=integral2(fun2, 0,dt/2,xmin,tan(dt/2));%;;

     
p=0.5*(tan(dt/2)^2/(tan(dt/2)^2+sigmasq))+(1/dt)*int1+(4/(pi^2*dt))*int2; %approximation for high SNR


%% actual match probability
   
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


