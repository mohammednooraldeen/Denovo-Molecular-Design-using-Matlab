function R = rotmat3( xyz,abc,uvw,theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
s_xyz=size(xyz,1);
for i=1:s_xyz;
x=xyz(i,1);y=xyz(i,2);z=xyz(i,3);
a=abc(1);b=abc(2);c=abc(3);
u=uvw(1);v=uvw(2);w=uvw(3);
L=u^2+v^2+w^2;


R(i,:)=[ ((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cosd(theta))+L*x*cosd(theta)+sqrt(L)*(-c*v+b*w-w*y+v*z)*sind(theta))/L;
    ((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cosd(theta))+L*y*cosd(theta)+sqrt(L)*(c*u-a*w+w*x-u*z)*sind(theta))/L;
    ((c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cosd(theta))+L*z*cosd(theta)+sqrt(L)*(-b*u+a*v-v*x+u*y)*sind(theta))/L]';

end

end



