function mutationkids = mutation_gene5(parents,options,NVARS,FitnessFcn, state, thisScore,thisPopulation,mutationRate_parameters)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

num_parents=length(parents);
%mutationkids=cell(num_parents,1)
mutationkids=cell(num_parents,1);

lb=mutationRate_parameters(1,:);
ub=mutationRate_parameters(2,:);
mutationRate=mutationRate_parameters(3,1);
for i=1:num_parents
    disp ('############num_parents= ##############'); num_parents
    disp (' ############thisPopulation=#############'); thisPopulation
    
    mutationkids{i,1}{1,1}=mutation_mol_4ed(thisPopulation{parents(i)}{1,1});
    if ~isempty(mutationRate)
    mutationkids{i,1}{2,1}=mutation_position(thisPopulation{parents(i)}{2,1},lb,ub,mutationRate);
    mutationkids{i,1}{3,1}=mutation_quatern(thisPopulation{parents(i)}{3,1},mutationRate);
    else
    mutationkids{i,1}{2,1}=mutation_position(thisPopulation{parents(i)}{2,1},lb,ub,0.6);  
    mutationkids{i,1}{3,1}=mutation_quatern(thisPopulation{parents(i)}{3,1},0.6);
    end
end
mutationkids{:,1}


end

