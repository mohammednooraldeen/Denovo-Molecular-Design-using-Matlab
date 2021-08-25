function out = crossover_quatern( A,B,method )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% three methods of crossover
% 1- rand mixture for respective dimensions of A and B 
% 2- rand average for respective dimensions of A and B
% 3- rand mixture+average
% methods


%method=randperm(3,1)
 
out=A;dim=[1 2 3];
rand_num_dimensions_exchange=randperm(3,1)
chose_dimensions_exchange=randperm(3,rand_num_dimensions_exchange)

    if (method==1)
for i=1:rand_num_dimensions_exchange
out(1,chose_dimensions_exchange(i))=B(chose_dimensions_exchange(i));
end

    elseif (method==2)
for i=1:rand_num_dimensions_exchange
out(1,chose_dimensions_exchange(i))= (out(1,chose_dimensions_exchange(i)) + B(chose_dimensions_exchange(i)))/2;
end
        
    elseif (method==3)
for i=1:rand_num_dimensions_exchange
out(1,chose_dimensions_exchange(i))=B(chose_dimensions_exchange(i));
end

remaining_dimensions=dim(~ismember(dim,chose_dimensions_exchange))
s_remaining_dimensions=size(remaining_dimensions); s_remaining_dimensions=s_remaining_dimensions(1,1)

    if (remaining_dimensions)
for i=1:s_remaining_dimensions
out(1,remaining_dimensions(i))= (out(1,remaining_dimensions(i)) + B(remaining_dimensions(i)))/2;
end
    end        

end

