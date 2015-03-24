% Script for the derivation of a general linear and angular velocity 
% Jacobian. Assumes strict Denavit-Hartenberg coordinate systems. Jacobian
% can be reduced for a given robot by substituting any values defined
% within the associated DH table.

function [trans, J] = rot_jac_mat(nDOF, rho)

% define unit vector
k = [0; 0; 1];

for n = 1:nDOF
    
    % define symbolic variables for each link as loop progresses
    a = sym(strcat('a',num2str(n)));
    alpha = sym(strcat('alpha',num2str(n)));
    d = sym(strcat('d',num2str(n)));
    theta = sym(strcat('theta',num2str(n)));
    
    % define A matrices according to DH convention
    A(:,:,n) = HRz(theta)*HTz(d)*HTx(a)*HRx(alpha);
    
    % define naming conventions for homogeneous transformation matrices,
    % both current and fixed
    H_name = strcat('H',num2str(n),num2str(n-1));
    H_name2 = strcat('H',num2str(n),'0');
    
    % define naming convention for extracted rotation matrix, fixed
    R_name = strcat('R',num2str(n),'0');
    
    % define naming convention for extracted origin translation matrix,
    % fixed
    O_name = strcat('O',num2str(n),'0');
    
    % assign homogenous transformation matrix, current
    trans.(H_name) = A(:,:,n);
    
    % combine current homogeneous transformation matrices to form fixed
    % matrices at each link level
    for i = 1:n
        
        if i == 1
            
            trans.(H_name2) = trans.H10;
            
        else
            
            trans.(H_name2) = trans.(H_name2)*trans.(...
                              strcat('H',num2str(i),...
                                     num2str(i-1)));
            
        end
        
    end
    
    % extract and store rotation and translation matrices, fixed
    R.(R_name) = trans.(H_name2)(1:3,1:3);
    O.(O_name) = trans.(H_name2)(1:3,4);
    
    % compute variable z as defined in eq. 4.47 on page 131 of SHV
    if n == 1
        
        z = sym(k);
        
    else
        
        z(:,n) = R.(strcat('R',num2str(n-1),'0'))*k;
        
    end
    
    % define angular velocity Jacobian according to eq. 4.48 on page 131 of
    % SHV
    Jw(:,n) = rho(n)*z(:,n);
    
end

% loop for definition of linear velocity Jacobian was separated because
% access to ON0 required for each value of m
for m = 1:nDOF
    
    if rho(m) == 0
        
        Jv(:,m) = z(:,m);
        
    elseif rho(m) == 1
        
        if m == 1
            
            Jv(:,m) = cross(z(:,m), (O.(strcat('O',num2str(nDOF),'0')) - [0; 0; 0]));
                   
        else
        
            Jv(:,m) = cross(z(:,m), (O.(strcat('O',num2str(nDOF),'0')) - O.(strcat('O',num2str(m-1),'0'))));
        
        end
    end
    
end

% combine linear and angular velocity Jacobians for complete definition of
% J
J = [Jv; Jw];
