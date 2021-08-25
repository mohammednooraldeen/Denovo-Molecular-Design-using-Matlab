function f = fitness_gene_sub_boxes_vectorized_2( chromosome, C,N,O,H,e,P_C,P_O,P_N,keep_original_charge)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

pop=length(chromosome);
%C=importdata(receptor_map,'',6);

for i=1:pop

%xyz=chromosome{i,1}{2,1}
%quat=chromosome{i,1}{3,1}
%mol=chromosome{i,1}{1,1};
mol=translate_mol_gene(chromosome{i,1})

%f(i,1)=xyz(1)+quat(2)*xyz(2)+2*quat(3);
f(i,1)=eval_xyz_sub_surfaces_3(mol,'C',C,N,O,H,e,P_C,P_O,P_N,keep_original_charge);

%cd C:\_malaria_research\matlab_fact
%dos('fred -receptor apo_generated.oeb.gz -dbase lig1.oeb.gz')
%cd d:\MATLAB\ttt
end



end

