function out = mutation_position( A,lb,ub,mutationRate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% A and lb and ub =[x y z]
out=A;
rand_num_dimensions=randperm(3,1);
chose_dimensions=randperm(3,rand_num_dimensions);
for i=1:rand_num_dimensions;

displacement=randsample(lb(chose_dimensions(i)):ub(chose_dimensions(i)),1)-out(chose_dimensions(i));
out(chose_dimensions(i))=out(chose_dimensions(i))+ mutationRate * displacement;
end


end

