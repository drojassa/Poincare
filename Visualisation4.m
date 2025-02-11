clear all
cla
close(gcf);
close(gcf);
fichier = 'theta_35_a_48_1_mode_#001.bmp';%nom du fichier de l'image à étudier
im = imread(fichier); %lecture du fichier




% xCenterEst = 375.5; %position du centre x à entrer manuellement
% yCenterEst = 325; %position du centre y à entrer manuellement


xCenterEst=576 ; %position du centre x à entrer manuellement
yCenterEst = 414; %position du centre y à entrer manuellement

radius = 500; %rayon à l'intérieur duquel les pixels de l'image sont analysés 

imClip= im; 

xSize = size(im,2); %dimension de l'image en x
ySize = size(im,1); %dimension de l'image en y
xAverage = 0;
yAverage = 0; 

sampleSize = 1000; 

D1=1.0081; %indice de distorsion (parfois ça aide à affiner les 'raies') 
D2=1/D1;

rSpace = linspace(0,radius,sampleSize);  %vecteur allant de 0 à 400 de dimension 1000
deltaRSpace = rSpace(2)-rSpace(1);    %incrément de rSpace
iR = zeros(1,sampleSize);           %vecteur nul de dimension 1000. 


for xStep = 1:1:xSize
    for yStep = 1:1:ySize      %pour chaque point de l'image...
        distFromCenter = sqrt(((xStep-xCenterEst)/D2)^2+((yStep-yCenterEst)/D1)^2);   %calcul de la distance du centre par le th. de Pythagore
        
        if(distFromCenter>radius)   %si la distance dépasse radius=400, mettre la valeur de l'intensité R, G, et B à 0
            imClip(yStep,xStep,1) = 0;
            imClip(yStep,xStep,2) = 0;
            imClip(yStep,xStep,3) = 0;
 
     
        else              %sinon...   
                        
            [d, ix] = min( abs( rSpace-distFromCenter ) );   %ix est l'indice pour la valeur de rSpace la plus proche de distFromCenter
            iR(ix) = iR(ix) + double(imClip(yStep,xStep,3))/deltaRSpace;
        end
        
        
    
        
    end
end

S=size(im);


f=7.4e-2; %distance focale de la lentille (qu'on ajustera manuellement pour obtenir la valeur d'intervalle spectral libre correcte.
n=1.5; %indice dans le Fabry-Perot
c=3e8;
lambda=1030e-9; %longueur d'onde 
nu_0=c/lambda;
pix=(8.77*10^-3/S(2)); %taille d'un pixel
r=rSpace*pix;  %distance physique sur la caméra
theta_out=atan(r/f); %angle d'inclinaison des rayons à la sorie du FP dans l'expérience
theta_in=asin(1/n.*sin(theta_out)); %angle d'inclinaison dans le Fabry-Perot calculé avec la loi de Snell-Descartes
nu_space=(nu_0./cos(theta_in)-nu_0)*1e-9; %relation entre la fréquence, à un facteur additif près, et l'angle dans le FP 


figure(1);
image(imClip(:,:,:)); %image brute
hold on
plot(xCenterEst,yCenterEst,'o');%positionnement de notre centre tel que choi aux lignes de code 8 et 9
hold off 


figure(2);
hold on;
plot(nu_space,iR*1e-5); %traçage du spectre
%ylim([0,3])
grid on





