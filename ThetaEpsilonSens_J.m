function [Theta,Epsilon,Sens] = ThetaEpsilonSens_J(J)
%Trouve l'angle d'inclinaison de l'ellipse (Theta),
%l'angle d'ellipticié (Epsilon), et le sens d'un vecteur de Jones J
%Sens: +1 = droite (horaire)
%      -1 = gauche (anti-horaire)  
%      
% le sens horaire/anti-horaire est pour celui qui observe
% la lumière arriver vers lui, comme sur le plan xy où 
% l'axe z sort de la page.

[Psi, Delta] = PsiDelta_J(J);

Theta=1/2*atan2(sin(2*Psi)*cos(Delta),cos(2*Psi));

Epsilon=-1/2*asin(sin(2*Psi)*sin(Delta));
Sens=sign(Epsilon);


end

