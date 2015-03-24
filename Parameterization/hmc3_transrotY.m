function T = hmc3_transrotY(translation, theta)
% Generates a homogeneous transformation matrix for a 3D translation followed
% by a rotation about the Y-axis.
%
% Input:
%	translation		(3 x 1 matrix) translation vector
%	theta			(scalar) amount of rotation in radians (positive is screw towards +Y)
%
% Output:
%	T		(4x4 matrix) homogeneous transformation matrix

	T = [	cos(theta)		0		sin(theta)		translation(1); ...
				0			1			0			translation(2); ...
			-sin(theta)		0     	cos(theta)		translation(3); ...
				0			0			0				1		  ];

end