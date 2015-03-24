function c=circleParam(theta,r,cx,cy,cz,phiA,psiA)

%syms theta r cx cy cz phiA psiA
p=[1 0 0]';

B=[cos(theta) sin(theta) 0;
    sin(theta) cos(theta) 0;
    0   0  0];

c_o=r*eye(3)*B*p;

c_o=[c_o;1];

A1=hmc3_transrotX([cx,cy,cz],0);
A2=hmc3_transrotX([0 0 0],phiA);
A3=hmc3_transrotY([0 0 0],psiA);
Ac=A1*A2*A3;

c=Ac*c_o;
%Use only the top 3
c=c(1:3);