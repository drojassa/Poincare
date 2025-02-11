function Visualisation4_JF_version_Theo_1(nomFichier,xCenterEst,yCenterEst)


im=readmatrix(strcat(nomFichier,'.csv'));

%transformation de la matrice (enlève un bruit de font 
S=size(im);
for k=1:S(1)
    for j=1:S(2)
         
        if im(k,j)<50; %la borne d'intensité maximal qu'on peut mettre a zéro pour élliminer le bruit.
            im(k,j)=0;
        end   
    end     

end


radius =800; %rayon à l'intérieur duquel les pixels de l'image sont analysés 


% on prend seulement une section des anneaux, thetamax indique a l'angle
% maximal ou on prend les valeurs d'intensité en référence ou coordonner
% cartésienne usuel
thetamax=56*pi/180;
thetamin=16*pi/180;


imClip= im;

xSize = size(imClip,2); %dimension de l'image en x
ySize = size(imClip,1); %dimension de l'image en y
xAverage = 0;
yAverage = 0; 

sampleSize = 1000; 

 D1=1.00; %indice de distorsion (parfois ça aide à affiner les 'raies') 

D2=D1;

rSpace = linspace(0,radius,sampleSize);  %vecteur allant de 0 à 400 de dimension 1000
deltaRSpace = rSpace(2)-rSpace(1);    %incrément de rSpace
iR = zeros(1,sampleSize);           %vecteur nul de dimension 1000. 
iG = zeros(1,sampleSize); 
iB = zeros(1,sampleSize); 



for xStep = 1:1:xSize
    for yStep = 1:1:ySize      %pour chaque point de l'image...
        
        theta=atan2(yStep-yCenterEst,xStep-xCenterEst);
        if (theta)<thetamax & theta>thetamin 
         
            
        distFromCenter = sqrt(((xStep-xCenterEst)/D2)^2+((yStep-yCenterEst)/D1)^2);   %calcul de la distance du centre par le th. de Pythagore
        
        
        if(distFromCenter>radius)   %si la distance dépasse radius=400, mettre la valeur de l'intensité R, G, et B à 0
            imClip(yStep,xStep) = 0;

 
     
        else              %sinon...   
                        
            [d, ix] = min( abs( rSpace-distFromCenter ) );   %ix est l'indice pour la valeur de rSpace la plus proche de distFromCenter
            iR(ix) = iR(ix) + double(imClip(yStep,xStep))/deltaRSpace;    %%note on enleve le dernier indice quand ce n'est plus des images RGB
%             iG(ix) = iG(ix) + double(imClip(yStep,xStep,2))/deltaRSpace; 
%             iB(ix) = iB(ix) + double(imClip(yStep,xStep,3))/deltaRSpace;
% %             
        end
        
        
    
        
        end
    end
end

S=size(im);


%la focale du la lentille utiliser. variable pour avoir le ISL du FP
f=7.4e-2;%4.795e-2; %distance focale de la lentille (qu'on ajustera manuellement pour obtenir la valeur d'intervalle spectral libre correcte.
n=1.5; %indice dans le Fabry-Perot
c=3e8;
lambda=1030e-9; %longueur d'onde 
nu_0=c/lambda;
pix=(8.77*10^-3/S(2)); %taille d'un pixel
r=rSpace*pix;  %distance physique sur la caméra
theta_out=atan(r/f); %angle d'inclinaison des rayons à la sorie du FP dans l'expérience
theta_in=asin(1/n.*sin(theta_out)); %angle d'inclinaison dans le Fabry-Perot calculé avec la loi de Snell-Descartes
nu_space=(nu_0./cos(theta_in)-nu_0)*1e-9; %relation entre la fréquence, à un facteur additif près, et l'angle dans le FP 


figure(6);
imagesc(imClip); %image brute
hold on
plot(xCenterEst,yCenterEst,'o');%positionnement de notre centre tel que choi aux lignes de code 8 et 9
axis equal
hold off 
colormap hot 


figure(9);
% subplot(3,1,1)
hold on;
plot(nu_space,iR*1e-5,'linewidth',3); %traçage du spectre
ylabel("intensity (arb. u.)",'fontsize',24)
xlabel(" {\it \Delta\nu} (GHz)",'fontsize',24)
set(gca, 'XTick', 0:50:800,'fontsize',24)



% set(gca,'Box', 'off','YTick',[], 'XGrid', 'off','Color','none')
% 
% 
% subplot(3,1,2)
% hold on;
% plot(nu_space,iG*1e-5); %traçage du spectre
% 
% subplot(3,1,3)
% hold on;
% plot(nu_space,iB*1e-5); %traçage du spectre

grid on
% legend

% 
%enregistre les fichiers csv en bmp
% figure('visible' , 'off')
% image=(imagesc(im));
% axis equal
% colormap hot
% saveas(image,strcat(nomFichier,'0.bmp'))

% 
