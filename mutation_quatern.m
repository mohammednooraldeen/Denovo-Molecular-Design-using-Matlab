function out = mutation_quatern( A, mutation_rate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
out=A;
rand_num_dimensions=randperm(3,1);
chose_dimensions=randperm(3,rand_num_dimensions);
for i=1:rand_num_dimensions

out(chose_dimensions(i))=out(chose_dimensions(i))+randperm(360,1) * mutation_rate
end


end

