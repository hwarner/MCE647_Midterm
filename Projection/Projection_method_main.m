function main

clear;
clc;

runID = strrep(strrep(strrep(datestr(now),'-','_'),':','_'),' ','_');

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

% First three points of Dr. Richter's sample data
theta1Array = [-.379; -.227; -.513];
theta2Array = [-.642; -.884; .495];
theta3Array = [.520; .914; -.759];

% Points 1, 6, and 9 of our collected data
% theta1Array = [-.74 -1.029 -.919];
% theta2Array = [-.04 -.8353 -.796];
% theta3Array = [-.349 .526 .704];

theta4 = 0;
theta5 = 0;

% initial guess for points P
x0 = 2*[rand(3,1); -rand(3,1); rand(3,1)];

% perform optimization to determine points P in world frame relative to end point
options = optimset('MaxFunEvals', 1500);
tic
[x1, fval1, exitflag1] = fmincon(@point_coord_cost_func, x0, [], [], [], [], [], [], @nonlcon, options, trans, d1, d2, d3, d4, d5, a1, a2, a3, a4, a5, alpha1, alpha2, alpha3, alpha4, alpha5, theta1Array, theta2Array, theta3Array, theta4, theta5);
toc
% recalculate endpoints
for i = 1:length(theta1Array)
    
    theta1 = theta1Array(i);
    theta2 = theta2Array(i);
    theta3 = theta3Array(i);
    
    H50 = subs(trans.H50);
    
    endx(i) = double(H50(1,4));
    endy(i) = double(H50(2,4));
    endz(i) = double(H50(3,4));
    
end

% add distances from world frame origin to endpoint to P values
P = [x1(1:3) + endx'; x1(4:6) + endy'; x1(7:9) + endz'];

% define individual points
P1 = [P(1); P(4); P(7)];
P2 = [P(2); P(5); P(8)];
P3 = [P(3); P(6); P(9)];

% calculate radius of final circle
R_final = norm(P1-P2)*norm(P2-P3)*norm(P3-P1)/(2*norm(cross((P1-P2),(P2-P3))))

% calculate center of final circle, http://self.gutenberg.org/articles/Circumscribed\_circle
alpha = norm(P2-P3)^2*dot((P1-P2),(P1-P3))/(2*norm(cross((P1-P2),(P2-P3)))^2);
beta = norm(P1-P3)^2*dot((P2-P1),(P2-P3))/(2*norm(cross((P1-P2),(P2-P3)))^2);
gamma = norm(P1-P2)^2*dot((P3-P1),(P3-P2))/(2*norm(cross((P1-P2),(P2-P3)))^2);
C = alpha*P1 + beta*P2 + gamma*P3;

% calculate vector normal to circle plane
P12 = P1-P2;
P13 = P1-P3;
nV = cross(P12,P13);

% calculate distance from endpoint to each point on the circle
dP1 = sqrt((endx(1)-P1(1))^2+(endy(1)-P1(2))^2+(endz(1)-P1(3))^2);
dP2 = sqrt((endx(2)-P2(1))^2+(endy(2)-P2(2))^2+(endz(2)-P2(3))^2);
dP3 = sqrt((endx(3)-P3(1))^2+(endy(3)-P3(2))^2+(endz(3)-P3(3))^2);

% define home position
theta1 = 0;
theta2 = -pi/4;
theta3 = -pi/2;

% endpoint coordinates in home position
H50 = double(subs(trans.H50));

% distance from home to circle center
dist_home_to_center =  sqrt(sum((H50(1:3,4) - C).^2))*100/2.54

disp('exitflag values of 0, -1, -2, and -3 do not indicate success')

% plot results
h = figure();
plot3(P(1:3), P(4:6), P(7:9), '*')
hold on;
plot3(C(1), C(2), C(3),'g*')
plot3(2.1192, -0.9218, 0.8367,'r*')
plotCircle3D(C', nV', .33) %function credit: Christian Reinbacher, MathWorks file exchange
legend('P_{est}', 'C_{est}', 'C_{given}')
xlabel('x_0, m')
ylabel('y_0, m')
zlabel('z_0, m')
axis([0, 2.5, -2.5, 0, 0, 2.5])
grid on;

save(strcat(runID,'_workspace.mat'),'P1','P2','P3','C','dP1','dP2','dP3','dist_home_to_center','fval1','exitflag1');

saveas(h,strcat(runID,'_plot'),'fig')

end

function cost = point_coord_cost_func(x, trans, d1, d2, d3, d4, d5, a1, a2, a3, a4, a5, alpha1, alpha2, alpha3, alpha4, alpha5, theta1Array, theta2Array, theta3Array, theta4, theta5)

% calculate endpoint and intermediate point (origin 4) in world coordinates
for i = 1:length(theta1Array)
    
    theta1 = theta1Array(i);
    theta2 = theta2Array(i);
    theta3 = theta3Array(i);
    
    H40 = subs(trans.H40);
    H50 = subs(trans.H50);
    
    inter4x(i) = double(H40(1,4));
    inter4y(i) = double(H40(2,4));
    inter4z(i) = double(H40(3,4));
    
    endx(i) = double(H50(1,4));
    endy(i) = double(H50(2,4));
    endz(i) = double(H50(3,4));
       
end

inter4x = inter4x';
inter4y = inter4y';
inter4z = inter4z';
endx = endx';
endy = endy';
endz = endz';

% assign optimization values to P
Px = x(1:3);
Py = x(4:6);
Pz = x(7:9);

% projection onto each plane
cost1 = (endz-inter4z)./(endx-inter4x)-Pz./Px
cost2 = (endz-inter4z)./(endy-inter4y)-Pz./Py
cost3 = (endy-inter4y)./(endx-inter4x)-Py./Px

% calculate cost measure
cost = sqrt(sum(cost1.^2 + cost2.^2 + cost3.^2))

end

function [c, ceq] = nonlcon(x, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~)

R = .33;

% calculate R according to R = abc/(4K) 
% from http://mathworld.wolfram.com/Circumradius.html
P1 = [x(1), x(4), x(7)];
P2 = [x(2), x(5), x(8)];
P3 = [x(3), x(6), x(9)];

d12 = sqrt((P1(1)-P2(1))^2+(P1(2)-P2(2))^2+(P1(3)-P2(3))^2);
d13 = sqrt((P1(1)-P3(1))^2+(P1(2)-P3(2))^2+(P1(3)-P3(3))^2);
d23 = sqrt((P2(1)-P3(1))^2+(P2(2)-P3(2))^2+(P2(3)-P3(3))^2);

s = (d12+d13+d23)/2;

K = sqrt(s*(d12+d13-s)*(d12+d23-s)*(d13+d23-s));

R_calc = d12*d13*d23/(4*K);

% inequality constraint |R-R_calc| <= .0001
c = abs(R - R_calc) - .0001;

% no nonlinear equality constraints
ceq = 0;

end
