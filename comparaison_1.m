
nomSerie = 'theta_40_a_48_1_mode';


%position du centre des anneaux
xCenterEst=696;
yCenterEst=508;

for k=1:5
    Visualisation4_JF_version_Theo_1(strcat(nomSerie, '_#', nombreString(k)),xCenterEst,yCenterEst)

end

% % % % boucle manuel pour trouver le centre
% for k =1:3
%     Visualisation4_JF_version_Theo_1((strcat(nomSerie,'_#002')),xCenterEst-2+2*k,yCenterEst)
%     legend
% end


function str = nombreString(nb)
    if(nb < 10)
        str = strcat('00', num2str(nb));
    else
        str = strcat('0', num2str(nb));
    end    
end