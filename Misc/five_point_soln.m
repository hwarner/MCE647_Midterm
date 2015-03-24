% method that solves system of 30 equations and 30 unknowns while
% "constraining the plane"

function main
clear;
clc;
close all

% setup generalized set of rotation matrices with 5 DOF, three revolute and
% two prismatic
[trans, ~] = rot_jac_mat(5, [1 1 1 0 0]);

% define system constants and DH convention parameters
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

% selected five points of our data, spreading them out along the circle's circumference
theta1Array = [-0.74  -0.911 -1.029 -0.919 -0.754];
theta2Array = [-0.04  -0.635 -0.835 -0.796 -1.442];
theta3Array = [-0.349  0.198  0.526  0.704  1.457];
theta4 = 0;
theta5 = 0;

% calculate location of origins four and five for each orientation
for i = 1:length(theta1Array)
    
    theta1 = theta1Array(i);
    theta2 = theta2Array(i);
    theta3 = theta3Array(i);
    
    H40(:,:,i) = double(subs(trans.H40));
    H50(:,:,i) = double(subs(trans.H50));
    
end

% random initial guess
x0 = 10*rand(30,1);

% use unbounded ranges, in general
lb = -inf*ones(30,1);
ub = inf*ones(30,1);

% constrain the plane by constraining the ranges of theta and phi, 
% conventions shown in image: 
% http://en.wikipedia.org/wiki/File:3D_Spherical_2.svg
lb(4) = 4*pi/3;  % 240 degrees
ub(4) = 5*pi/3;  % 300 degrees --> 60 degree range
lb(5) = pi/3;    %  60 degrees
ub(5) = 2*pi/3;  % 120 degrees --> 60 degree range

options = optimset('MaxFunEvals',10000);
[x,fval,exitflag] = fmincon(@cost_func,x0,[],[],[],[],lb,ub,@nonlinconst,options,H40,H50);

% name the results
cx    = x(1)
cy    = x(2)
cz    = x(3)
theta = x(4)
phi   = x(5)
tl    = x(6:10)
tc    = x(11:15)
px    = x(16:20)
py    = x(21:25)
pz    = x(26:30)

% vector normal to circle for plotting
nc = cross([px(1), py(1), pz(1)]-[px(2),py(2),pz(2)],[px(1),py(1),pz(1)]-[px(3),py(3),pz(3)]);

% calculate out the parametric lines for plotting
for tl = -10:1:10
    line1(:,tl+11) = H50(1:3,4,1)+tl*(H50(1:3,4,1)-H40(1:3,4,1));
    line2(:,tl+11) = H50(1:3,4,2)+tl*(H50(1:3,4,2)-H40(1:3,4,2));
    line3(:,tl+11) = H50(1:3,4,3)+tl*(H50(1:3,4,3)-H40(1:3,4,3));
    line4(:,tl+11) = H50(1:3,4,4)+tl*(H50(1:3,4,4)-H40(1:3,4,4));
    line5(:,tl+11) = H50(1:3,4,5)+tl*(H50(1:3,4,5)-H40(1:3,4,5));
    
end

% plot results
plot3(cx, cy, cz, '*r')
hold on;
plot3(px(1),py(1),pz(1),'*b')
plot3(px(2),py(2),pz(2),'*b')
plot3(px(3),py(3),pz(3),'*b')
plot3(px(4),py(4),pz(4),'*b')
plot3(px(5),py(5),pz(5),'*b')
plotCircle3D([cx,cy,cz],nc,.33)
plot3(line1(1,:),line1(2,:),line1(3,:),'--g')
plot3(line2(1,:),line2(2,:),line2(3,:),'--g')
plot3(line3(1,:),line3(2,:),line3(3,:),'--g')
plot3(line4(1,:),line4(2,:),line4(3,:),'--g')
plot3(line5(1,:),line5(2,:),line5(3,:),'--g')
xlabel('x')
ylabel('y')
zlabel('z')
axis([-3 3 -3 3 0 3])
daspect([1 1 1])
grid on
view(-170,20)

% check radius distances
r1 = sqrt((cx-px(1))^2+(cy-py(1))^2+(cz-pz(1))^2);
r2 = sqrt((cx-px(2))^2+(cy-py(2))^2+(cz-pz(2))^2);
r3 = sqrt((cx-px(3))^2+(cy-py(3))^2+(cz-pz(3))^2);
r4 = sqrt((cx-px(4))^2+(cy-py(4))^2+(cz-pz(4))^2);
r5 = sqrt((cx-px(5))^2+(cy-py(5))^2+(cz-pz(5))^2);

% define home position
theta1 = 0;
theta2 = 0;
theta3 = pi/4;

% endpoint coordinates in home position
H50 = double(subs(trans.H50));

% distance from home to circle center
dist_home_to_center =  sqrt(sum((H50(1:3,4) - [cx; cy; cz]).^2))*100/2.54

% final cost
fval

end

function F = cost_func(x,H40,H50)

r = .33;

% parse inputs
cx    = x(1);
cy    = x(2);
cz    = x(3);
theta = x(4);
phi   = x(5);
tl   = x(6:10);
tc    = x(11:15);
px     = x(16:20);
py     = x(21:25);
pz     = x(26:30);

% define center
C = [cx; cy; cz];

% define vector n
n = [cos(phi)*sin(theta); sin(theta)*sin(phi); cos(theta)];

% define vector u
u = [-sin(phi); cos(phi); 0];

n_cross_u = cross(n,u);

% for each point calculate the circle and line equations
for i = 1:5
    
    % circle equation
    F_1(:,i) = r*cos(tc(i))*u + r*sin(tc(i))*n_cross_u + C - [px(i); py(i); pz(i)];
    
    % a point on the line
    r_0 = H50(1:3,4,i);
    
    % a vector parallel to the line
    v_parallel = r_0 - H40(1:3,4,i);
    
    % line equation
    F_2(:,i) = r_0 + tl(i)*v_parallel - [px(i); py(i); pz(i)];
    
end

% sum of the squares cost function
F = sqrt(sum(([reshape(F_1,[],1);reshape(F_2,[],1)]).^2));

end

function [c, ceq] = nonlinconst(x,~,~)

r = .33;

% parse inputs
cx    = x(1);
cy    = x(2);
cz    = x(3);
theta = x(4);
phi   = x(5);
tl   = x(6:10);
tc    = x(11:15);
px     = x(16:20);
py     = x(21:25);
pz     = x(26:30);

% calculate radius constraint |R-R_est|<=1e-10
r1 = abs(r - sqrt((cx-px(1))^2+(cy-py(1))^2+(cz-pz(1))^2))-1e-10;
r2 = abs(r - sqrt((cx-px(2))^2+(cy-py(2))^2+(cz-pz(2))^2))-1e-10;
r3 = abs(r - sqrt((cx-px(3))^2+(cy-py(3))^2+(cz-pz(3))^2))-1e-10;
r4 = abs(r - sqrt((cx-px(4))^2+(cy-py(4))^2+(cz-pz(4))^2))-1e-10;
r5 = abs(r - sqrt((cx-px(5))^2+(cy-py(5))^2+(cz-pz(5))^2))-1e-10;

c = [r1; r2; r3; r4; r5];

% no equality constraints
ceq = 0;

end