clear all;
close(gcf);
close(gcf);
close(gcf);
%IMPORT DATA
mode_alpha = xlsread("xyz_fonction_alpha.xlsx");
%DEFINE ANGLES IN THE SHEET
alpha=mode_alpha(:,2)-(7-abs(mode_alpha(:,2)-90)*7/(90)); %correction d'alpha
alpha=-alpha*pi/180;
theta1_mode1=(mode_alpha(:,8)-153)*pi/180;
theta2_mode1=(mode_alpha(:,9)-46)*pi/180;


%PAIR OF ANGLES FOR EACH MODE
xi_mode1=theta1_mode1;
chi_mode1=xi_mode1-theta2_mode1-pi/2;


%matrice de transmission du miroir de sortie du resonnateur
n1=1;
n2=1.5;

delta = -11.8*pi/180;

theta = 46*pi/180;
theta2=asin((n1/n2).*sin(theta));

t1s=(2*n1.*cos(theta))/(n1.*cos(theta)+n2.*cos(theta2));
t2s=(2*n2.*cos(theta2))/(n1.*cos(theta)+n2.*cos(theta2));
t1p=(2*n1.*cos(theta))/(n2.*cos(theta)+n1.*cos(theta2));
t2p=(2*n2.*cos(theta2))/(n2.*cos(theta)+n1.*cos(theta2));

rs=t1s.*t2s.*t1s.*t2s;
rp=-t1p.*t2p.*t1p.*t2p;

%matrice du miroir M2 (this mirror is rotated the alpha angle)
M2=[rp,0;0,rs];


%matrice du miroir M1 (fixed) 
r11=-exp(1i*(delta/2));
r12=exp(-1i*(delta/2));

M1=[r11,0;0,r12]; %de la base x'y' à xy

ts=t1s*t2s;
tp=t1p*t2p;

Mout=[ts,0;0,tp];

%calcul coordonnés sur sphère de poincarre
X1=zeros(numel(xi_mode1),1);
Y1=zeros(numel(xi_mode1),1);
Z1=zeros(numel(xi_mode1),1);
error1=zeros(numel(xi_mode1,3));


deltaxi=5*pi/180;
deltachi=5*pi/180;

for j=1:numel(xi_mode1)
    
    X1(j)=cos(2*xi_mode1(j))*cos(2*chi_mode1(j));
    error1(j,1)=((2*sin(2*xi_mode1(j))*cos(2*chi_mode1(j))))*deltaxi+((2*cos(2*xi_mode1(j))*sin(2*chi_mode1(j))))*deltachi;
    Y1(j)=sin(2*xi_mode1(j))*cos(2*chi_mode1(j));
    error1(j,2)=((2*cos(2*xi_mode1(j))*cos(2*chi_mode1(j))))*deltaxi^2+((2*sin(2*xi_mode1(j))*sin(2*chi_mode1(j))))*deltachi;
    Z1(j)=-sin(2*chi_mode1(j));
    error1(j,3)=((2*cos(2*chi_mode1(j))))*deltachi;
end


%modele theorique
alpha_teo=-[30:0.001:65]*pi/180;
res=length(alpha_teo);%resolution of the theoric curve
for k=1:res;
    %Alpha Rotation matrix 
    T=[cos(alpha_teo(k)),sin(alpha_teo(k));-sin(alpha_teo(k)),cos(alpha_teo(k))]; 

    J_RT=T*M1*T*M2; %matrice d'un aller-retour

    [V,D] = eig(J_RT); %calcul des valeurs propres et vecteurs propres de J_RT
    D1(k)=D(1,1);D2(k)=D(2,2); %valeurs propres des deux modes 

    V1=V(:,1);V2=V(:,2);%les deux vecteurs propres
    
    V1=Mout*V1;V2=Mout*V2; %on multiplie par la matrice en transmission du miroir de sortie 


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


end; %fin de la boucle sur alpha

figure(1)

subplot(311); errorbar(-alpha*180/pi, X1,error1(:,1),'g*');hold on;plot(-alpha_teo*180/pi,x1,'r-');hold on;plot(-alpha_teo*180/pi,x2,'g-');grid on;ylabel('X');axis([-alpha_teo(1)*180/pi -alpha_teo(length(alpha_teo))*180/pi -1 1]);
subplot(312); errorbar(-alpha*180/pi, Y1,error1(:,2),'g*');hold on;plot(-alpha_teo*180/pi,y1,'r-');hold on;plot(-alpha_teo*180/pi,y2,'g-');grid on;ylabel('Y');axis([-alpha_teo(1)*180/pi -alpha_teo(length(alpha_teo))*180/pi -1 1]);
subplot(313); errorbar(-alpha*180/pi, Z1,error1(:,3),'g*');hold on;plot(-alpha_teo*180/pi,z1,'r-');hold on;plot(-alpha_teo*180/pi,z2,'g-');grid on;ylabel('Z');axis([-alpha_teo(1)*180/pi -alpha_teo(length(alpha_teo))*180/pi -1 1]);
xlabel('alpha')
