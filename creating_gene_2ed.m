function out = creating_gene_2ed(num_atoms,position,lb,ub,O_num,N_num )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%lb_position=[-7.9 1.6 -13.3];
%ub_position=[15.2 22.6 7.6];
atom_types=['C' 'O' 'N'];
s_atom_types=size(atom_types,2);
atom_type=atom_types(randsample(s_atom_types,1));
out{1,1}=make_sp2_8ed(num_atoms,O_num,N_num)
if ~isempty(position)
out{2,1}=position;
else
out{2,1}=[randsample(lb(1):ub(1),1),randsample(lb(2):ub(2),1),randsample(lb(3):ub(3),1)]
end
out{3,1}=[randsample(360,1),randsample(360,1),randsample(360,1)]

end

