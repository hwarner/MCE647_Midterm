
%% PUMA Robot

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


theta1Array = [-.379; -.227; -.513];
theta2Array = [-.642; -.884; .495];
theta3Array = [.520; .914; -.759];
theta4 = 0;
theta5 = 0;

% Subsitute in constants into DH Matrices
namesOfMatrices=fieldnames(trans);
for i=1:length(namesOfMatrices)
   H=trans.(namesOfMatrices{i});
   A.(namesOfMatrices{i})=subs(H,{'b','k','d1','d2','d3','d4','d5','a1','a2','a3','a4','a5',...
       'alpha1','alpha2','alpha3','alpha4','alpha5','theta4','theta5'},...
       [b,k,d1,d2,d3,d4,d5,a1,a2,a3,a4,a5,...
       alpha1,alpha2,alpha3,alpha4,alpha5,theta4,theta5]);  
end

% Make a function out of A5_0
matlabFunction(A.H50,'file','transA5_0','vars',{'theta1','theta2','theta3'});
matlabFunction(A.H40,'file','transA4_0','vars',{'theta1','theta2','theta3'});
matlabFunction(A.H30,'file','transA3_0','vars',{'theta1','theta2','theta3'});
matlabFunction(A.H20,'file','transA2_0','vars',{'theta1','theta2','theta3'});
matlabFunction(A.H10,'file','transA1_0','vars',{'theta1','theta2','theta3'});


