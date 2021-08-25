function your_rotated_coord = rotmat( axis_vector,theta,your_coord )
%UNTITLED Summary of this function goes here
%give your coordinates as 3*n matrix for x,y,z
%and the axis vector also 3*1 matrix
%   Detailed explanation goes here
c=cosd(theta);
s=sind(theta);
axis_vector=axis_vector';
Ax=axis_vector(1,1);
Ay=axis_vector(2,1);
Az=axis_vector(3,1);

R=[ c+(1-c)*Ax^2   (1-c)*Ax*Ay-s*Az  (1-c)*Ax*Az+s*Ay     ;...
   (1-c)*Ax*Ay+s*Az  c+(1-c)*Ay^2  (1-c)*Ay*Az-s*Ax ;...
   (1-c)*Ax*Az-s*Ay   (1-c)*Ay*Az+s*Ax  c+(1-c)*Az^2];

your_rotated_coord=R*your_coord';
your_rotated_coord=your_rotated_coord';

end

