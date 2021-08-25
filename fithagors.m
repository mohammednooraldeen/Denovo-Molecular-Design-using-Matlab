function  out  = fithagors( a ,b )
%FITHAGORS Summary of this function goes here
%   Detailed explanation goes here

sa=size(a);
sb=size(b);

    
for i=1:sa(1,1)
out(i,1)= roundn(( (a(i,1)-b(1,1))^2+(a(i,2)-b(1,2))^2+ (a(i,3)-b(1,3))^2)^0.5,-4);
end


out;
end

