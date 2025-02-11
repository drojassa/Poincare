clear all;
close(gcf);
close(gcf);
close(gcf);
close(gcf);
%résultats donnés dans la base xy 
%valvuls en fonction du déphasage du miroir M2
n1=1;
n2=1.5;

theta =[0:1:90].*pi/180;
kappa=[];
delta = 12.8*pi/180;
%kappa=[acosh(1/cos(delta/2))-0.5:0.01:acosh(1/cos(delta/2))+0.5];

for k=1:length(theta),
theta2=asin((n1/n2).*sin(theta(k)));

t1s=(2*n1.*cos(theta(k)))/(n1.*cos(theta(k))+n2.*cos(theta2));
t2s=(2*n2.*cos(theta2))/(n1.*cos(theta(k))+n2.*cos(theta2));
t1p=(2*n1.*cos(theta(k)))/(n2.*cos(theta(k))+n1.*cos(theta2));
t2p=(2*n2.*cos(theta2))/(n2.*cos(theta(k))+n1.*cos(theta2));

rs=t1s.*t2s.*t1s.*t2s;
rp=t1p.*t2p.*t1p.*t2p;

kappa(k)=(-1/2).*log(rs/rp);


r21=-exp(kappa(k));
r22=exp(-kappa(k));
%matrice du miroir M2
M2=exp(-kappa(k))*[r21,0;0,r22]; %de la base xy à x'y'

alpha=[0:1:90]*pi/180;%valeur de l'angle alpha permettant l'obtention d'un PE

for l=1:length(alpha)

r11=-exp(1i*(delta/2));
r12=exp(-1i*(delta/2));
M1=[r11,0;0,r12]; %de la base x'y' à xy

T=[cos(alpha(l)/2),sin(alpha(l)/2);-sin(alpha(l)/2),cos(alpha(l)/2)]; 

J_RT=T*M1*T*T*M2*T; %matrice d'un aller-retour

[V,D] = eig(J_RT); %calcul des valeurs propres et vecteurs propres de J_RT
D1(k,l)=D(1,1);D2(k,l)=D(2,2); %valeurs propres des deux modes 

V1=V(:,1);V2=V(:,2);%les deux vecteurs propres
V1_cp=T*M2*T*V(:,1); %état propre du mode 1 se propageant en sens contraire
V1_cp=V1_cp/norm(V1_cp); %normalisation de V1_cp
V2_cp=T*M1*T*V(:,2);%état propre du mode 2 se propageant en sens contraire
V2_cp=V2_cp/norm(V2_cp);%normalisation de V2_cp
V1_cp(1)=-V1_cp(1);%réécriture de V1_cp dans la même base que V1 
V2_cp(1)=-V2_cp(1);%réécriture de V2_cp dans la même base que V2 

Ortho1(l)=abs(V1_cp'*V1); %calcul de l'orthogonalité des ondes contra-propagatives du mode 1
Ortho2(l)=abs(V2_cp'*V2); %calcul de l'orthogonalité des ondes contra-propagatives du mode 2



%paramètres de l'ellipse de polarisation des deux modes 
phi_V1(l)=angle(V1(2)/V1(1));
phi_V2(l)=angle(V2(2)/V2(1));
epsilon_V1(l)=atan(abs(V1(2))/abs(V1(1)));
epsilon_V2(l)=atan(abs(V2(2))/abs(V2(1)));

phi_V1_cp(l)=angle(V1_cp(2)/V1_cp(1));
phi_V2_cp(l)=angle(V2_cp(2)/V2_cp(1));
epsilon_V1_cp(l)=atan(abs(V1_cp(2))/abs(V1_cp(1)));
epsilon_V2_cp(l)=atan(abs(V2_cp(2))/abs(V2_cp(1)));


end; %fin de la boucle sur alpha
end;

dif=abs(abs(D1)-abs(D2));
figure(1); %valeurs propres et vecteurs propres
%hold on;surf(theta, alpha,abs(D1));grid on;
hold on;imagesc(alpha*180/pi,theta*180/pi,dif);colormap turbo; colorbar;
xticks(0:10:90)
xlim([0;90])
ylim([0;90])
xlabel('alpha')
ylabel('theta')