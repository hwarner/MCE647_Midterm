
%% This script is used to perfom error estimation of the system

%% Get the Robot transformation matrix
[trans, ~] = rot_jac_mat(5, [1 1 1 0 0]);

%% Get the Equation of the circle
syms thetaC r cx cy cz phiA psiA theta1 theta2 theta3 alphaR xp yp zp
c=circleParam(thetaC,r,cx,cy,cz,phiA,psiA);  %Symbolic calculation of points on circle
c(1)=c(1)-xp;
c(2)=c(2)-yp;
c(3)=c(3)-zp;
c(1)=solve(c(1),cx);
c(2)=solve(c(2),cy);
c(3)=solve(c(3),cz);

%% Tranformation matrix at the laser
A5_0=trans.H50;

%% Equation of the line
l=lineParam(A5_0,alphaR);

%% Substitute line point into circle enter eq
cx=subs(c(1),xp,l(1));
cy=subs(c(2),yp,l(2));
cz=subs(c(3),zp,l(3));

%% Calculate the position of the center
R=(cx^2 + cy^2 + cz^2)^(1/2);
varList=symvar(R);

%Use the Jacobian to get sensitivities
J=jacobian(R,varList);


%Everything above was symbolic, not subs in the values into the jacobian
%to get sensivities
b = .270;
k = .0521;

d1 = .666;
d2 = -.2435;
d3 = .0934;
d4 = k;
d5 = 0;

a1 = 0;
a2 = .4318;
a3 = 0;
a4 = 0;
a5 = b;

alpha1 = -pi/2;
alpha2 = 0;
alpha3 = 0;
alpha4 = 0;
alpha5 = 0;


theta1Array = [-.379; -.227; -.513];
theta2Array = [-.642; -.884; .495];
theta3Array = [.520; .914; -.759];
theta4 = 0;
theta5 = 0;

r=0.330;

theta1 = -.379;
theta2 = -.642;
theta3 = .520;
alphaR=1.4389;
phiA=0.2528;
psiA=-2.3937;
thetaC=0.3923;
eval(J)'
symvar(J).'