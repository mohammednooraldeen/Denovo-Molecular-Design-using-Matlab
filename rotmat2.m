function R = rotmat2( xyz,abc,uvw,theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x=xyz(1);y=xyz(2);z=xyz(3);
a=abc(1);b=abc(2);c=abc(3);
u=uvw(1);v=uvw(2);w=uvw(3);


R=[ (a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cosd(theta))+x*cosd(theta)+(-c*v+b*w-w*y+v*z)*sind(theta);
    (b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cosd(theta))+y*cosd(theta)+(c*u-a*w+w*x-u*z)*sind(theta);
    (c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cosd(theta))+z*cosd(theta)+(-b*u+a*v-v*x+u*y)*sind(theta)];
end

