function mol = view_sp2_full( mol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vdw_r=zeros(10,4);
vdw_r(1,1:2)=['H'  37];        %H
vdw_r(2,1:4)=['C' 77 67 60];  %C
vdw_r(3,1:3)=['O' 66 57];     %O
vdw_r(4,1:3)=['N'  74 65];     %N
vdw_r(5,1:3)=['S'  104 95];    %S
vdw_r(6,1:2)=['P' 110];       %P
vdw_r(7,1:2)=['c'  99];        %Cl
vdw_r(8,1:2)=['F' 64];        %F
vdw_r(9,1:2)=['b' 114];       %Br
vdw_r(10,1:2)=['I' 113] ;     %I
        
s_mol=size(mol,1);
mov=cell(s_mol,s_mol);
for i=1:s_mol
    for j=1:4
        atom_num_2=mol(i,5+j)
         if atom_num_2>0
        bond_type=mol(i,9+j)
        atom_num_1=i
        atom_type_1=mol(i,14)
        atom_vdw_1=vdw_r(ismember(vdw_r(:,1),atom_type_1,'rows'),1+bond_type)
        
        atom_num_2=mol(i,5+j)
        atom_type_2=mol(atom_num_2,14)
        atom_vdw_2=vdw_r(ismember(vdw_r(:,1),atom_type_2,'rows'),1+bond_type)
        
        current_bond_distance=roundn(fithagors(mol(atom_num_1,2:4),mol(atom_num_2,2:4)),-4)
        real_bond_distance=roundn((atom_vdw_1+atom_vdw_2)/100,-4)
        ratio=roundn(real_bond_distance/current_bond_distance,-3)
        
        current_disp=bsxfun(@minus,mol(atom_num_2,2:4),mol(atom_num_1,2:4))
        modified_disp=current_disp*(1-ratio); 
        if (~ismember(atom_num_2,mov{:,1})
            mov{i,1}=atom_num_2;
            mov{i,end}=modified_disp
            modified_coord=bsxfun(@minus,mol(atom_num_2:end,2:4),modified_disp)
          % mol(~ismember( mol(:,2:4),mol(atom_num_1,2:4),'rows'),2:4)=bsxfun(@minus,mol(~ismember( mol(:,2:4),mol(atom_num_1,2:4),'rows'),2:4),modified_disp)
            [a,b,c]=find(ismember(mol(:,6:9),atom_num_2))
            for t=1:size(a,1)
                if (~ismember(a(t),mov{:,1})
                    par_coord=mov(ismember([move{:,1}],t(a),'rows'),2);
                    par_coord=par_coord{1,1};
                    modified_coord=bsxfun(@plus,modified_coord,par_coord);
                end
            end
          
            mol(atom_num_2:end,2:4)=modified_coord
       
        else
           
        end

        
         else j=5;
         end
         
    end
end
        
  mol      

scatter3 (mol(:,2),mol(:,3),mol(:,4),5,'fill','red');hold on;
text(mol(:,2),mol(:,3),char(mol(:,14)))
for i=1:s_mol
    for j=1:4
        if (mol(i,5+j)~=0)
            mol2=mol(i,5+j); 
            switch (mol(i,9+j))  case 1 ;op='-b' ; case 2 ;op='-g'; case 3 ;op=':b' ;end
            plot3([mol(i,2);mol(mol2,2)],...
                            [mol(i,3);mol(mol2,3)],...
                            [mol(i,4);mol(mol2,4)],op,'LineWidth',3)
        end
        end
    end
 
for i=1:s_mol
 s=size(mol);
         for i=1:s(1,1); 
             for j=1:s(1,1); 
                 if ( roundn(((mol(i,2)-mol(j,2))^2+(mol(i,3)-mol(j,3))^2+(mol(i,4)-mol(j,4))^2)^0.5,-2) ==1.33)
                     % ((p(i,1)-p(j,1))^2+(p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2)^0.5 ==bond_length_4  )  
                  
                     plot3([mol(i,2);mol(j,2)],...
                            [mol(i,3);mol(j,3)],...
                            [mol(i,4);mol(j,4)],'-g','LineWidth',3) 
                 else
                            end;
                            end
                            
                            ;end
 xlabel('X');
 ylabel('Y');
 zlabel('Z');
 view([0 0 20])
 axis equal
 hold off;

end

