clear all;
close(gcf);
close(gcf);
close(gcf);
close(gcf);

n1=1;
n2=1.5;

theta =[0:1:70].*pi/180;
kappa=[];
delta = 11.8*pi/180;

for k=1:length(theta),
theta2=asin((n1/n2).*sin(theta(k)));

t1s=(2*n1.*cos(theta(k)))/(n1.*cos(theta(k))+n2.*cos(theta2));
t2s=(2*n2.*cos(theta2))/(n1.*cos(theta(k))+n2.*cos(theta2));
t1p=(2*n1.*cos(theta(k)))/(n2.*cos(theta(k))+n1.*cos(theta2));
t2p=(2*n2.*cos(theta2))/(n2.*cos(theta(k))+n1.*cos(theta2));

rs=-t1s.*t2s.*t1s.*t2s;
rp=t1p.*t2p.*t1p.*t2p;

kappa(k)=(-1/2).*log(rs/rp);


r21=-exp(kappa(k));
r22=exp(-kappa(k));
%matrice du miroir M2
M2=[rp,0;0,rs];

alpha=[0:1:90]*pi/180;
for l=1:length(alpha)

r11=-exp(1i*(delta/2));
r12=exp(-1i*(delta/2));
M1=[r11,0;0,r12]; %de la base x'y' à xy

T=[cos(alpha(l)),sin(alpha(l));-sin(alpha(l)),cos(alpha(l))]; 

J_RT=T*M1*T*M2; %matrice d'un aller-retour

[V,D] = eig(J_RT); %calcul des valeurs propres et vecteurs propres de J_RT
D1(k,l)=D(1,1);D2(k,l)=D(2,2); %valeurs propres des deux modes 

V1=V(:,1);V2=V(:,2);%les deux vecteurs propres


m(k,l)=max(abs(D1(k,l)),abs(D2(k,l)));
end; %fin de la boucle sur alpha
end;


figure(1); %maximum des valeurs propres 
hold on;imagesc(alpha*180/pi,theta*180/pi,m);colormap turbo;colorbar;
xticks(0:10:90)
xlim([0;90])
ylim([0;70])
xlabel('alpha')
ylabel('theta')

figure(2); %rapport entre D1 et D2 echelle logarithmique
rap_vp=log(abs(D1))-log(abs(D2));
hold on;imagesc(alpha*180/pi,theta*180/pi,rap_vp);colormap turbo; colorbar;
xticks(0:10:90)
xlim([0;90])
ylim([0;70])
xlabel('alpha')
ylabel('theta')

figure(3); %differénce de phase entre D1 et D2
dif_phase=abs(angle(D1)-angle(D2));
hold on;imagesc(alpha*180/pi,theta*180/pi,dif_phase);colormap turbo; colorbar;
xticks(0:10:90)
xlim([0;90])
ylim([0;70])
xlabel('alpha')
ylabel('theta')
