function out = mutation_gene( gene )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

out{1,1}=mutation_mol_2ed(gene{1,1})
out{2,1}=mutation_quatern(gene{2,1})
out{3,1}=mutation_position(gene{3,1},[0 0 0],[100 100 100])


end

