function c=circleParam(beta,r,cx,cy,cz,phiA,psiA)

%circleParam - Create the equation for a circle
%
%
%c=circleParam(beta,r,cx,cy,cz,phiA,psiA)
%       Inputs:
%               beta - the angle on the circle
%               r - radius of the circle
%               cx,cy,cz - location of the circle's center
%               phiA - rotation of the circle about the global X plane
%               psiA - rotation of the circle about it's Y plane
%
%       Outputs:
%               c - point on the circle (at Beta) in global cs.
%
%        Circle initially is in the X-Y plane and then is rotated about
%        it's X axis (phiA) and then about it's Y axis (psiA)


% Create a vector basis
p=[1 0 0]';

% Make a circle in the X-Plane
B=[cos(beta) sin(beta) 0;
    sin(beta) cos(beta) 0;
    0   0  0];

c_o=r*eye(3)*B*p;


c_o=[c_o;1];  % Append a "1" so that it can be used with Trans Matrix

% Create Transformation Matrix
A1=hmc3_transrotX([cx,cy,cz],0);
A2=hmc3_transrotX([0 0 0],phiA);
A3=hmc3_transrotY([0 0 0],psiA);
Ac=A1*A2*A3;

%Rotate and Translate the circle
c=Ac*c_o;

%The upper 3 then contain the vector
c=c(1:3);