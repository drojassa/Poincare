clear all;
close(gcf);
close(gcf);
%IMPORT DATA
mode_theta = xlsread("xyz_fonction_theta.xlsx");

%DEFINE ANGLES IN THE SHEET
theta_exp=(mode_theta(:,1))*pi/180;
theta1_mode1=(mode_theta(:,5)-153)*pi/180;
theta2_mode1=(mode_theta(:,6)-46)*pi/180;
theta1_mode2=(mode_theta(:,7)-153)*pi/180;
theta2_mode2=(mode_theta(:,8)-46)*pi/180;



%PAIR OF ANGLES FOR EACH MODE
xi_mode1=theta1_mode1;chi_mode1=xi_mode1-theta2_mode1-pi/2;
xi_mode2=theta1_mode2;chi_mode2=xi_mode2-theta2_mode2-pi/2;

n1=1;n2=1.5;

theta =0:0.01:73*pi/180;
delta = 11.8*pi/180;

alpha=45*pi/180;
N=13.8e20;
cs=0.15e-20;
d=0.1;

%calcul coordonnés sur sphère de poincarre
X1=zeros(numel(xi_mode1),1);Y1=zeros(numel(xi_mode1),1);Z1=zeros(numel(xi_mode1),1);
error1=zeros(numel(xi_mode1,3));
X2=zeros(numel(xi_mode2),1);Y2=zeros(numel(xi_mode2),1);Z2=zeros(numel(xi_mode2),1);
error2=zeros(numel(xi_mode2,3));

deltaxi=5*pi/180;deltachi=5*pi/180;

for j=1:numel(xi_mode1)   
    X1(j)=cos(2*xi_mode1(j))*cos(2*chi_mode1(j));
    error1(j,1)=((2*sin(2*xi_mode1(j))*cos(2*chi_mode1(j))))*deltaxi+((2*cos(2*xi_mode1(j))*sin(2*chi_mode1(j))))*deltachi;
    Y1(j)=sin(2*xi_mode1(j))*cos(2*chi_mode1(j));
    error1(j,2)=((2*cos(2*xi_mode1(j))*cos(2*chi_mode1(j))))*deltaxi^2+((2*sin(2*xi_mode1(j))*sin(2*chi_mode1(j))))*deltachi;
    Z1(j)=-sin(2*chi_mode1(j));
    error1(j,3)=((2*cos(2*chi_mode1(j))))*deltachi;
end
for j=1:numel(xi_mode2)  
    X2(j)=cos(2*xi_mode2(j))*cos(2*chi_mode2(j));
    error2(j,1)=((2*sin(2*xi_mode2(j))*cos(2*chi_mode2(j))))*deltaxi+((2*cos(2*xi_mode2(j))*sin(2*chi_mode2(j))))*deltachi;
    Y2(j)=sin(2*xi_mode2(j))*cos(2*chi_mode2(j));
    error2(j,2)=((2*cos(2*xi_mode2(j))*cos(2*chi_mode2(j))))*deltaxi+((2*sin(2*xi_mode2(j))*sin(2*chi_mode2(j))))*deltachi;
    Z2(j)=-sin(2*chi_mode2(j));
    error2(j,3)=((2*cos(2*chi_mode2(j))))*deltachi;
end

for k=1:length(theta)
theta2=asin((n1/n2).*sin(theta(k)));

t1s=(2*n1.*cos(theta(k)))/(n1.*cos(theta(k))+n2.*cos(theta2));
t2s=(2*n2.*cos(theta2))/(n1.*cos(theta(k))+n2.*cos(theta2));
t1p=(2*n1.*cos(theta(k)))/(n2.*cos(theta(k))+n1.*cos(theta2));
t2p=(2*n2.*cos(theta2))/(n2.*cos(theta(k))+n1.*cos(theta2));

rs=-t1s.*t2s.*t1s.*t2s;
rp=t1p.*t2p.*t1p.*t2p;

%matrice du miroir M2
M2=[rp,0;0,rs];

ts=t1s*t2s;
tp=t1p*t2p;
Mout=[ts,0;0,tp];
%matrice de passage de la base des axes principaux à celle commune aux 2 miroirs
%matrice du miroir M1
r11=-exp(1i*(delta/2));
r12=exp(-1i*(delta/2));
M1=[sqrt(0.935)*r11,0;0,sqrt(0.930)*r12]; %de la base x'y' à xy

T=[cos(alpha),sin(alpha);-sin(alpha),cos(alpha)]; 

J_RT=T*M1*T*M2; %matrice d'un aller-retour

[V,D] = eig(J_RT); %calcul des valeurs propres et vecteurs propres de J_RT 
D1(k)=D(1,1);D2(k)=D(2,2);%valeurs propres des deux modes 
V1=V(:,1);V2=V(:,2);%les deux vecteurs propres
V1=Mout*V1;V2=Mout*V2;

%paramètres de l'ellipse de polarisation des deux modes 
phi_V1(k)=angle(V1(2)/V1(1));
phi_V2(k)=angle(V2(2)/V2(1));
epsilon_V1(k)=atan(abs(V1(2))/abs(V1(1)));
epsilon_V2(k)=atan(abs(V2(2))/abs(V2(1)));

% calcul des coordonnées des modes sur la sphère de Poincaré
x1(k)=cos(2*epsilon_V1(k));
y1(k)=sin(2*epsilon_V1(k))*cos(phi_V1(k));
z1(k)=sin(2*epsilon_V1(k))*sin(phi_V1(k));
x2(k)=cos(2*epsilon_V2(k));
y2(k)=sin(2*epsilon_V2(k))*cos(phi_V2(k));
z2(k)=sin(2*epsilon_V2(k))*sin(phi_V2(k));


end%fin de la boucle sur theta


figure(1); %coordonnées des vecteurs propres sur la SP (un mode en rouge, l'autre en bleu)
subplot(311); errorbar(theta_exp*180/pi, X1,error1(:,1),'g*');hold on;errorbar(theta_exp*180/pi, X2,error2(:,1),'r*');hold on;plot(theta*180/pi, x1,'r-');hold on; plot(theta*180/pi, x2,'g-');grid on;ylabel('X');axis([theta(1)*180/pi theta(length(theta))*180/pi -1 1]);
subplot(312); errorbar(theta_exp*180/pi, Y1,error1(:,2),'g*');hold on;errorbar(theta_exp*180/pi, Y2,error2(:,2),'r*');hold on;plot(theta*180/pi, y1,'r-');hold on; plot(theta*180/pi, y2,'g-');grid on;ylabel('Y');axis([theta(1)*180/pi theta(length(theta))*180/pi -1 1]);
subplot(313); errorbar(theta_exp*180/pi, Z1,error1(:,3),'g*');hold on;errorbar(theta_exp*180/pi, Z2,error2(:,3),'r*');hold on;plot(theta*180/pi, z1,'r-');hold on; plot(theta*180/pi, z2,'g-');grid on;ylabel('Z');axis([theta(1)*180/pi theta(length(theta))*180/pi -1 1]);
xlabel('theta')

figure(2); %proximité des deux états
errorbar(theta_exp*180/pi,X1.*X2+Y1.*Y2+Z1.*Z2,error1(:,1)+error1(:,2)+error1(:,3)+error2(:,1)+error2(:,2)+error2(:,3),'k*');hold on;plot(theta*180/pi, x1.*x2+y1.*y2+z1.*z2);axis([theta(1)*180/pi 48.7 -1 1]);
xlabel('theta')







