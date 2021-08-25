function xyz = xyz_of(coord,A,B,C,spacing,xc,yc,zc )
%XYZ_OF Summary of this function goes here
%   Detailed explanation goes here
xyz(:,1)=xc-(((A/2)+1-coord(:,1))*spacing);
xyz(:,2)=yc-(((B/2)+1-coord(:,2))*spacing);
xyz(:,3)=zc-(((C/2)+1-coord(:,3))*spacing);
xyz(:,4)=coord(:,4);

end

