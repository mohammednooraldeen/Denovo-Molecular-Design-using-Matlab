function mol = apply_mol( gene )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



mol=gene{1,1}
xyz=gene{2,1}; x=xyz(1); y=xyz(2); z=xyz(3);
quat=gene{3,1}; angle_about_x=quat(1); angle_about_y=quat(2); angle_about_z=quat(3);

mol_centre= mean(mol(:,2:4));
mol(:,2:4)= bsxfun(@minus,mol(:,2:4),mol_centre)
%view_sp2(mol); hold on;
mol(:,2:4)= rotmat3(mol(:,2:4),[0 0 0],[1 0 0],angle_about_x); % rotation about x axis 
mol(:,2:4)= rotmat3(mol(:,2:4),[0 0 0],[0 1 0],angle_about_y); % rotation about x axis 
mol(:,2:4)= rotmat3(mol(:,2:4),[0 0 0],[0 0 1],angle_about_z); % rotation about x axis
%view_sp2(mol); hold off;

mol(:,2:4)= bsxfun(@plus,mol(:,2:4),xyz)





end

