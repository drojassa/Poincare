clear all;
close(gcf);
close(gcf);
%IMPORT DATA
mode_theta = xlsread("Rapport_d_extintion_theta.xlsx");

%DEFINE ANGLES IN THE SHEET
theta_exp=(mode_theta(:,1));
Vmin=(mode_theta(:,7));
Vmax=(mode_theta(:,8));


ext_ratio=Vmax./Vmin;
deg_pol= (ext_ratio-1)./(ext_ratio+1);
plot(theta_exp,deg_pol,'b*-')

ylabel('Degré de polarisation')
xlabel('theta (deg)')
