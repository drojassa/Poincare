function [Psi, Delta] = PsiDelta_J(J)
%Trouve les paramètres Psi et Delta de vecteur de Jones J

Eo=dot(J,J);
Eox=norm(J(1));
Eoy=norm(J(2));
Psi=pi/2;
if Eox~=0 
    Psi=atan(Eoy/Eox);
end
Delta=angle(J(2)/J(1));
end

