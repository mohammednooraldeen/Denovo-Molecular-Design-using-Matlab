function child = couple( mol1,mol2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

s_1=size(mol1); s_1=s_1(1,1);
s_2=size(mol2); s_2=s_2(1,1);
child=mol1;
mol1_fusion_part=mol1(ismember(mol1(:,2:4),mol2(:,2:4),'rows'),:)
mol2_fusion_part=mol2(ismember(mol2(:,2:4),mol1(:,2:4),'rows'),:)

mol2_pool=(~ismember(mol2(:,2:4),mol1(:,2:4),'rows'),:)
s_mol2_pool=size(mol2_pool);s_mol2_pool=s_mol2_pool(1,1);

for i=1:s_mol2_pool
choices= mol2_pool(ismember(mol2_pool(:,6),mol2_fusion_part(:,1),'rows'),:)

rand_chosen=choices(find(choices(:,1)==(randsample(choices(:,1),1))),:)
choices= choices(find(choices(:,1)~=rand_chosen(:,1)),:)
child(end+1,:)=rand_chosen;
mol2_fusion_part(end+1,:)=rand_chosen;
end




 



end

