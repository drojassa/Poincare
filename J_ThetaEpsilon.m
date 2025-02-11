function [J] = J_ThetaEpsilon(Theta,Epsilon)
%trouve le vecteur de Jones à partir des paramètres Thea (angle
%d'inclinaison) et Epsilon (angle du rapport axe_mineur/axe_majeur) 

J(1,1)=cos(Theta)*cos(Epsilon)+i*sin(Theta)*sin(Epsilon);
J(2,1)=sin(Theta)*cos(Epsilon)-i*cos(Theta)*sin(Epsilon);
phi=angle(J(1,1));
J=J*exp(-i*phi);
end

