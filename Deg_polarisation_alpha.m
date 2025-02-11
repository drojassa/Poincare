clear all;
close(gcf);
close(gcf);
%IMPORT DATA
mode_alpha = xlsread("xyz_fonction_alpha.xlsx");

%DEFINE ANGLES IN THE SHEET
alpha_exp=mode_alpha(:,2)-(7-abs(mode_alpha(:,2)-90)*7/(90));
deg_pol=mode_alpha(:,10);

plot(alpha_exp,deg_pol,'b*-')
ylabel('Degré de polarisation')
xlabel('alpha (deg)')
