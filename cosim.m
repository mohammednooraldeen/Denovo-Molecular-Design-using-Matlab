function r = cosim( mat1,mat2 )
%COSIM Summary of this function goes here
%   Detailed explanation goes here
r= dot(mat1,mat2)/(norm(mat1,2)*norm(mat2,2));
 %x=5;for i=1:71; c(i,1)=dot(zinc_60_per(i,:),net_trained.iw{1,1}(x,:))/(norm(zinc_60_per(i,:),2)*norm(net_trained.iw{1,1}(x,:),2));c(i,2)=i;end; sortrows(c,1)
end

