function mol_updated = mutation_mol_4( mol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

s_mol=size(mol); s_mol=s_mol(1,1);
mol_updated=mol;

iterations=1;
reserved_atoms=mol(1,:);
reserved_rand_chosen_root=mol(1,:);

for i=1:iterations;


    
s_mol_updated=size(mol_updated,1);
choices_updated=mol_updated((find(mol_updated(:,1)>1)),:)
choices_updated=choices_updated(~ismember(choices_updated(:,2:4),reserved_atoms(:,2:4),'rows'),:) %the atoms should not be from reserved_atoms pool which include those previously modified

s_choices_updated=size(choices_updated);s_choices_updated=s_choices_updated(1,1);

if (s_choices_updated>0)






link_hir(1,3)=100;
while (find(link_hir(:,3)>1.33))
        rand_chosen=choices_updated(randperm(s_choices_updated,1),:)
        rand_chosen_root=mol_updated(rand_chosen(1,6),:);
        reserved_rand_chosen_root(end+1,:)=rand_chosen_root
        link_hir=linkage(mol_updated(~ismember(mol_updated(:,2:4),rand_chosen(:,2:4),'rows'),2:4),'single');   
end
        mol_updated(ismember(mol_updated(:,2:4),rand_chosen_root(:,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),rand_chosen_root(:,2:4),'rows'),5)+1
        mol_updated=mol_updated(~ismember(mol_updated(:,2:4),rand_chosen(:,2:4),'rows'),:)
        s_mol_updated=size(mol_updated),s_mol_updated=s_mol_updated(1,1);

        %renumbering of remaining atoms

        %rand_chosen_bounded_atoms=mol(ismember(mol_updated(:,6),rand_chosen(1,1),'rows'),:); %substituted rand_chosen at connection column of other atoms with new nearest atom
        %s_rand_chosen_bound_atoms=size(rand_chosen_bounded_atoms);s_rand_chosen_bound_atoms=s_rand_chosen_bound_atoms(1,1)
        %if (s_rand_chosen_bound_atoms>0)
        %for x=1:s_rand_chosen_bound_atoms
         %   possible_attach_list=sortrows(mol_updated(find((fithagors(mol_updated(:,2:4),rand_chosen_bounded_atoms(x,2:4))==1.33) & mol_updated(:,5)>=1),:),5)
           % if (size(possible_attach_list,1))mol_updated(ismember(mol_updated(:,2:4),rand_chosen_bounded_atoms(x,2:4),'rows'),6)=possible_attach_list(1,1); 
          %  mol_updated(ismember(mol_updated(:,2:4),possible_attach_list(1,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),possible_attach_list(1,2:4),'rows'),5)-1;
         %   end

        %end
        end
        mol_mol(:,1)=mol_updated(:,1);
        mol_mol(:,2)=[1:s_mol_updated];
        mol_updated(:,1)=[1:s_mol_updated];
        for i=1:s_mol_updated
            for j=1:s_mol_updated
            if (mol_updated(i,6)==mol_mol(j,1)) mol_updated(i,6)=mol_mol(j,2);j=s_mol_updated+1;end
            end
        end


%replacing of removed atom
mol_updated
addition_choices=mol_updated(mol_updated(:,5)>=1 & mol_updated(:,6)>=1,:)
addition_choices=addition_choices(~ismember(addition_choices(:,2:4),reserved_rand_chosen_root(:,2:4),'rows'),:) %all except the rand_chosen_root
s_addition_choices=size(addition_choices);s_addition_choices=s_addition_choices(1,1);
%try tomorrow INSHA ALLAH make addition choices as available point(x,y,z)
%to avoid repeated calculation of point and overlaping with existed atoms
if(s_addition_choices>=1)
    i=0;addition=1;
while (addition) 
rand_addition_choice=addition_choices(randperm(s_addition_choices,1),:);
rand_addition_choice_root=rand_addition_choice(1,6);
rand_addition_choice
     theta=120;
     theta2=180;
     if(rand(1)>0.5) theta_sign=-1 ;else theta_sign=1; end
    if(rand(1)>0.5) theta2_sign=1;else theta2_sign=0; end
    
     point=roundn(rotmat3(mol_updated(rand_addition_choice_root,2:4)-rand_addition_choice(1,2:4),[0 0 0],[0 0 1],theta_sign*theta)+rand_addition_choice(1,2:4),-4);
     answer=1;
    answer=ismember(point(1,:),mol_updated(:,2:4));answer=answer(1)*answer(2)*answer(3);
    
    if(~answer) mol_updated(rand_addition_choice_root,5)=mol_updated(rand_addition_choice_root,5)-1; 
        mol_updated(s_mol_updated+1,:)=[s_mol_updated+1,point,2,rand_addition_choice_root]; 
        reserved_atoms(end+1,:)=mol_updated(end,:);addition=0;
        mol_updated(ismember(mol_updated(:,2:4),rand_chosen_root(:,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),rand_chosen_root(:,2:4),'rows'),5)+1;
        mol_updated(ismember(mol_updated(:,2:4),rand_addition_choice(:,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),rand_addition_choice(:,2:4),'rows'),5)-1;
        
    else
        %i=i+1;
        point=roundn(rotmat3(mol_updated(rand_addition_choice_root,2:4)-rand_addition_choice(1,2:4),[0 0 0],[0 0 1],-1*theta_sign*theta)+rand_addition_choice(1,2:4),-4);
       answer=1;
        answer=ismember(point(1,:),mol_updated(:,2:4));answer=answer(1)*answer(2)*answer(3);
       if(~answer) mol_updated(rand_addition_choice_root,5)=mol_updated(rand_addition_choice_root,5)-1;  
        mol_updated(s_mol_updated+1,:)=[s_mol_updated+1,point,2,rand_addition_choice_root]; 
        reserved_atoms(end+1,:)=mol_updated(end,:);addition=0; 
       mol_updated(ismember(mol_updated(:,2:4),rand_chosen_root(:,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),rand_chosen_root(:,2:4),'rows'),5)+1;end
       mol_updated(ismember(mol_updated(:,2:4),rand_addition_choice(:,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),rand_addition_choice(:,2:4),'rows'),5)-1;
       
    end

end



%renumbering



;
else i=iterations+1;
 end   
end



%draw original
p=mol;
scatter3 (p(:,2),p(:,3),p(:,4));hold on;
 
 s=size(p);
         for i=1:s_mol; 
             for j=1:s_mol; 
                 if ( roundn(((p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2+(p(i,4)-p(j,4))^2)^0.5,-2) ==1.33)
                     % ((p(i,1)-p(j,1))^2+(p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2)^0.5 ==bond_length_4  )  
                  
                     plot3([p(i,2);p(j,2)],...
                            [p(i,3);p(j,3)],...
                            [p(i,4);p(j,4)],'-b','LineWidth',10) 
                 else
                            end;
                            end
                            
                            ;end

                   
                        
                        
%draw mutated
s_mol_updated=size(mol_updated);s_mol_updated=s_mol_updated(1,1);
p=mol_updated;
scatter3 (p(:,2),p(:,3),p(:,4));hold on;
 
 s=size(p);
         for i=1:s_mol_updated; 
             for j=1:s_mol_updated; 
                 if ( roundn(((p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2+(p(i,4)-p(j,4))^2)^0.5,-2) ==1.33)
                     % ((p(i,1)-p(j,1))^2+(p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2)^0.5 ==bond_length_4  )  
                  
                     plot3([p(i,2);p(j,2)],...
                            [p(i,3);p(j,3)],...
                            [p(i,4);p(j,4)],'-g','LineWidth',5) 
                 else
                            end;
                            end
                            
                            ;end
 xlabel('X');
 ylabel('Y');
 zlabel('Z');
view ([0,0,100])
                         hold off;
                        
                        

end
