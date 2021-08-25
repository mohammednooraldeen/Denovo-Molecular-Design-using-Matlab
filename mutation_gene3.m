function mutationkids = mutation_gene3(parents ,options,NVARS, ...
    FitnessFcn, state, thisScore,thisPopulation,mutationRate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

num_parents=length(parents);
%mutationkids=cell(num_parents,1)
mutationkids=cell(num_parents,1);

parfor i=1:num_parents
    disp ('############num_parents= ##############'); num_parents
    disp (' ############thisPopulation=#############'); thisPopulation
    
    mutationkids{i,1}{1,1}=mutation_mol_2ed(thisPopulation{parents(i)}{1,1});
    mutationkids{i,1}{2,1}=mutation_quatern(thisPopulation{parents(i)}{2,1});
    mutationkids{i,1}{3,1}=mutation_position(thisPopulation{parents(i)}{3,1},[0 0 0],[100 100 100]);
    
end
mutationkids{:,1}


end

