function kids  = crossover_gene7(parents,options,NVARS, ...
    FitnessFcn,thisScore,thisPopulation,energy_weight,C,O,N,H,e,P_C,P_O,P_N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
parents

num_parents=length(parents) %(thisPopulation)
%if ((num_parents/2)==roundn(num_parents/2,0)); num_parents=num_parents/2 ;else num_parents=(num_parents-1)/2; end
num_kids=num_parents/2;
kids=cell(num_kids,1);


count=0;
for i=1:num_kids
  count=count+1;
  %kids{i,1}{1,1}=crossover_mol_2ed(thisPopulation{parents(count)}{1,1},thisPopulation{parents(count+1)}{1,1});
  %kids{i,1}{1,1}=crossover_mol_3ed(thisPopulation{parents(count)}{1,1},thisPopulation{parents(count+1)}{1,1},energy_weight,C,O,N,H,e,P_C,P_O,P_N);
  kids{i,1}{1,1}=crossover_mol_6ed(thisPopulation{parents(count)}{1,1},thisPopulation{parents(count+1)}{1,1},energy_weight,C,O,N,H,e,P_C,P_O,P_N);
  kids{i,1}{2,1}=crossover_quatern(thisPopulation{parents(count)}{2,1},thisPopulation{parents(count+1)}{2,1},2);
  kids{i,1}{3,1}=crossover_quatern(thisPopulation{parents(count)}{3,1},thisPopulation{parents(count+1)}{3,1},2);
  count=count+1
end

    
    
kids{:,:}
end

