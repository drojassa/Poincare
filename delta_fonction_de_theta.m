clear all;
close(gcf);
close(gcf);
close(gcf);
close(gcf);

n1=1;
n2=1.5;

theta =[45:0.01:60]*pi/180;
kappa=[];
delta=[];
%kappa=[acosh(1/cos(delta/2))-0.5:0.01:acosh(1/cos(delta/2))+0.5];
%vkappa =acosh(1/cos(delta/2))

for k=1:length(theta),
theta2=asin((n1/n2).*sin(theta(k)));

t1s=(2*n1.*cos(theta(k)))/(n1.*cos(theta(k))+n2.*cos(theta2));
t2s=(2*n2.*cos(theta2))/(n1.*cos(theta(k))+n2.*cos(theta2));
t1p=(2*n1.*cos(theta(k)))/(n2.*cos(theta(k))+n1.*cos(theta2));
t2p=(2*n2.*cos(theta2))/(n2.*cos(theta(k))+n1.*cos(theta2));

rs=t1s.*t2s.*t1s.*t2s;
rp=t1p.*t2p.*t1p.*t2p;

kappa(k)=(-1/2).*log(rs/rp);
delta(k)=2*acos(1/cosh(kappa(k)));
end

plot(theta.*180/pi,delta*180/pi)
xlabel('theta')
ylabel('delta')