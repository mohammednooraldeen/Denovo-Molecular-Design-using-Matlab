function f = simulated_pos_quat(A,gene,receptor_map )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
C=importdata(receptor_map,'',6);
x=A(1); y=A(2); z=A(3);
a=A(4); b=A(5); c=A(6);
[x,y,z]
%[a,b,c]
gene{2,1}=[x,y,z]
%gene{3,1}=[a,b,c]

mol=translate_mol_gene(gene);

f=eval_xyz_2(mol,'C',C)

end

