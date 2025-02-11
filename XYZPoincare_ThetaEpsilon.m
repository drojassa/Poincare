function [r] = XYZPoincare_ThetaEpsilon(Theta,Epsilon)
%trouve les coordonnées x,y et z sur la sphère de Poincarré
%d'un état de polarisation à partir des paramètre Theta et Epsilon d'une
%ellipse

r(1)=cos(2*Theta)*cos(2*Epsilon);
r(2)=sin(2*Theta)*cos(2*Epsilon);
r(3)=-sin(2*Epsilon);

end

