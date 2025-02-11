clear all;
close(gcf);
close(gcf);
close(gcf);
close(gcf);

%IMPORT DATA
mode_theta = xlsread("Iseuil_theta.xlsx");
%DEFINE ANGLES IN THE SHEET
theta_exp=(mode_theta(:,1))*pi/180;

IPseuil=mode_theta(:,3);
Puissance=2.8126*IPseuil-2.0613; %Passage de courrent à Puissance

n1=1;n2=1.5;

theta =0:0.01:76*pi/180;
delta = 11.8*pi/180;

N=13.8e20;
cs=0.15e-20;
d=0.1;

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
alpha=45*pi/180;
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

omega_max(k)= max(abs(D1(k)),abs(D2(k)));

end%fin de la boucle sur theta

omega_max_interp = interp1(theta, omega_max, theta_exp, 'linear'); 
% Define the function to fit
func = @(params,theta_exp) params(1)*(N*cs*d - log(1 - params(2)) - 2*log(omega_max_interp));
x0 = [10, 0];  % x0(1) corresponds to C and x0(2) corresponds to A
x = lsqcurvefit(@(params,theta_exp) func(params,theta_exp), x0, theta_exp, Puissance);
% Extract the fitted parameters
C = x(1);A = x(2);

% fprintf('%d\n', C)
% fprintf('%d\n', A)


Pth=C*(N*cs*d-log(1-A)-2*log(omega_max));
figure(1);
errorbar(theta_exp*180/pi, Puissance,0.010*2.8126,0.010*2.8126,1,1,'b*-');hold on; plot(theta*180/pi, Pth);
ylabel('Pseuil (W)');
xlabel('theta (deg)')






