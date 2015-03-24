function T = hmc3_transrotX(translation, theta)
% Generates a homogeneous transformation matrix for a 3D translation followed
% by a rotation about the X-axis.
%
% Input:
%	translation		(3 x 1 matrix) translation vector
%	theta			(scalar) amount of rotation in radians (positive is screw towards +X)
%
% Output:
%	T		(4x4 matrix) homogeneous transformation matrix

	T = [	1			0				0			translation(1); ...
			0		cos(theta)		-sin(theta)		translation(2); ...
			0		sin(theta)     	cos(theta)		translation(3); ...
			0			0				0				1		  ];

end