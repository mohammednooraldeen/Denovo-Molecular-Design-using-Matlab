function mol_updated = mutation_mol_3( mol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

s_mol=size(mol); s_mol=s_mol(1,1);
mol_updated=mol;

iterations=1;
reserved_atoms=mol(1,:);
reserved_rand_chosen_addition=mol(1,:);

%while (iterations<s_mol);

addition_choices_updated=mol_updated((find(mol_updated(:,5)>=1)),:)
addition_choices_updated=addition_choices_updated(~ismember(addition_choices_updated(:,2:4),reserved_atoms(:,2:4),'rows'),:) %the atoms should not be from reserved_atoms pool which include those previously modified
s_choices_updated=size(addition_choices_updated);s_choices_updated=s_choices_updated(1,1);

if (s_choices_updated>0)




addition=1;
while (addition) 
     
     rand_chosen_addition=addition_choices_updated(randperm(s_choices_updated,1),:)
     rand_chosen_addition_root=mol_updated(rand_chosen_addition(1,6),:)
     theta=120;
     theta2=180;
     if(rand(1)>0.5) theta_sign=-1 ;else theta_sign=1; end
    if(rand(1)>0.5) theta2_sign=1;else theta2_sign=0; end
    s_mol_updated=size(mol_updated,1);
     point=roundn(rotmat3(rand_chosen_addition_root(1,2:4)-rand_chosen_addition(1,2:4),[0 0 0],[0 0 1],theta_sign*theta)+rand_chosen_addition(1,2:4),-4)
     %pole=p(choice,2:4)-origin;
    
     %point=roundn(rotmat3(point-origin,origin-origin,p(choice,2:4)-origin,theta2_sign*theta2)+origin,-4);
    
     answer=1;
    answer=ismember(point(1,:),mol_updated(:,2:4),'rows');
    
    if(~answer) mol_updated(ismember(mol_updated(:,2:4),rand_chosen_addition(1,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),rand_chosen_addition(1,2:4),'rows'),5)-1; 
        mol_updated(s_mol_updated+1,:)=[s_mol_updated+1,point,2,rand_chosen_addition(1,1)]; 
        addition=0;
        
        
    else
     point=roundn(rotmat3(rand_chosen_addition_root(1,2:4)-rand_chosen_addition(1,2:4),[0 0 0],[0 0 1],-1*theta_sign*theta)+rand_chosen_addition(1,2:4),-4) 
     answer=1;
        answer=ismember(point(1,:),mol_updated(:,2:4),'rows');
       if(~answer) mol_updated(ismember(mol_updated(:,2:4),rand_chosen_addition(1,2:4),'rows'),5)=mol_updated(ismember(mol_updated(:,2:4),rand_chosen_addition(1,2:4),'rows'),5)-1; 
      mol_updated(s_mol_updated+1,:)=[s_mol_updated+1,point,2,rand_chosen_addition(1,1)]
        addition=0;
       end

    end

end

%%% look for atom to delete %%%
if (~addition)
reserved_rand_chosen_addition(end+1,:)=mol_updated(end,:)

link_hir(1,3)=100;count=0;
while (mean(link_hir(:,3))>1.33)
        del_choices=mol_updated(~ismember(mol_updated(:,2:4),reserved_rand_chosen_addition(:,2:4),'rows'),:)
        s_del_choices=size(del_choices,1)
        rand_del_choice=del_choices(randsample(s_del_choices,1),:)
        link_hir=linkage(mol_updated(~ismember(mol_updated(:,2:4),rand_del_choice(1,2:4),'rows'),2:4),'single')
        if (mean(link_hir(:,3))>1.33) reserved_rand_chosen_addition(end+1,:)=rand_del_choice(1,:);end
        reserved_rand_chosen_addition
        %count=count+1
        %if(count>10) reserved_rand_chosen_addition=mol(1,:);end
end

% if (rand_del_choice)
        rand_del_choice_root=rand_del_choice(1,6)
        mol_updated(ismember(mol_updated(:,1),rand_del_choice_root,'rows'),5)=mol_updated(ismember(mol_updated(:,1),rand_del_choice_root,'rows'),5)+1
        mol_updated=mol_updated(~ismember(mol_updated(:,2:4),rand_del_choice(:,2:4),'rows'),:)
        s_mol_updated=size(mol_updated,1);

        molold_molnew(:,1)=mol_updated(:,1);
        molold_molnew(:,2)=[1:s_mol_updated];
        mol_updated(:,1)=[1:s_mol_updated];
        for i=1:s_mol_updated
            for j=1:s_mol_updated
            if (mol_updated(i,6)==molold_molnew(j,1)) mol_updated(i,6)=molold_molnew(j,2);j=s_mol_updated+1;end
            end
        end
        %iterations=s_mol;
 
 %end

end

 end   
%end



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
