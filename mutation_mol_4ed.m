function mol_updated = mutation_mol_4ed( mol )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


disp('MUATION FUNCTION mutation_mol_2ed')
s_mol=size(mol); s_mol=s_mol(1,1);
mol_updated=mol;

iterations=1;
reserved_atoms=mol(1,:);
reserved_rand_chosen_addition=mol(1,:);

%while (iterations<s_mol);

addition_choices_updated=mol_updated((find(mol_updated(:,5)<3 )),:) %exclude atoms of 3 connections and oxygen atoms also
addition_choices_updated=addition_choices_updated(~ismember(addition_choices_updated(:,2:4),reserved_atoms(:,2:4),'rows'),:) %the atoms should not be from reserved_atoms pool which include those previously modified
s_choices_updated=size(addition_choices_updated);s_choices_updated=s_choices_updated(1,1);

if (s_choices_updated>0)




addition=1;
while (addition) 
     
     rand_chosen_addition=addition_choices_updated(randperm(s_choices_updated,1),:)
     rand_chosen_addition_roots=rand_chosen_addition(1,6:8); rand_chosen_addition_roots=rand_chosen_addition_roots(rand_chosen_addition_roots~=0)
     rand_chosen_addition_roots=mol_updated(ismember(mol_updated(:,1),rand_chosen_addition_roots),:)
     s_rand_chosen_addition_roots=size(rand_chosen_addition_roots,1);
     theta=120;
     theta2=180;
     if(rand(1)>0.5) theta_sign=-1 ;else theta_sign=1; end
    if(rand(1)>0.5) theta2_sign=1;else theta2_sign=0; end
    s_mol_updated=size(mol_updated,1); 
    answer=1;count=1
    while (answer==1 & count<=s_rand_chosen_addition_roots)
        
        point=roundn(rotmat3(rand_chosen_addition_roots(count,2:4)-rand_chosen_addition(1,2:4),[0 0 0],[0 0 1],theta_sign*theta)+rand_chosen_addition(1,2:4),-4)
        %pole=p(choice,2:4)-origin;
        %point=roundn(rotmat3(point-origin,origin-origin,p(choice,2:4)-origin,theta2_sign*theta2)+origin,-4);
        answer=ismember(point(1,:),mol_updated(:,2:4),'rows');
        if(~answer) 
             mol_updated(s_mol_updated+1,1:4)=[s_mol_updated+1,point]; %,2,rand_chosen_addition(1,1)]; 
             empty_bonds_loc=find(mol_updated(end,6:9)==0,1);
             addition=0;answer=0;
        else
           point=roundn(rotmat3(rand_chosen_addition_roots(count,2:4)-rand_chosen_addition(1,2:4),[0 0 0],[0 0 1],-1*theta_sign*theta)+rand_chosen_addition(1,2:4),-4)
        
                if(~answer) 
             mol_updated(s_mol_updated+1,1:4)=[s_mol_updated+1,point]; %,2,rand_chosen_addition(1,1)]; 
             empty_bonds_loc=find(mol_updated(end,6:9)==0,1);
             addition=0;answer=0;
        
        
        end
        %if( ~p(i,6)>0) p(i,6)=choice; elseif (~p(i,7)>0) p(i,7)=choice; elseif (~p(i,8)>0) p(i,8)=choice;end; 
        %p(i,5)=numel(find(p(i,6:8)>0));
        
        
     
        end
       count=count+1;
        
    end
    end

end

%%% look for atom to delete %%%
if (~addition)
reserved_rand_chosen_addition(end+1,:)=mol_updated(end,:)

link_hir(1,3)=100;count=0;
while (mean(link_hir(:,3))>1.33)
        %%%exclude heteroatoms delete if the root of the added atom is
        %%%heteroatom
        if (rand_chosen_addition(end,14)>67);
            del_choices=mol_updated(~ismember(mol_updated(:,2:4),reserved_rand_chosen_addition(:,2:4) ,'rows' ) & mol_updated(:,14)==67 ,:)
            
        else
            
            del_choices=mol_updated(~ismember(mol_updated(:,2:4),reserved_rand_chosen_addition(:,2:4),'rows'),:)
        end
        
        s_del_choices=size(del_choices,1)
        
        if (s_del_choices<1) rand_del_choice=[0, point];
        else rand_del_choice=del_choices(randsample(s_del_choices,1),:)
        end
        
        link_hir=linkage(mol_updated(~ismember(mol_updated(:,2:4),rand_del_choice(1,2:4),'rows'),2:4),'single')
        if (mean(link_hir(:,3))>1.33) reserved_rand_chosen_addition(end+1,:)=rand_del_choice(1,:);end
        reserved_rand_chosen_addition
        %count=count+1
        %if(count>10) reserved_rand_chosen_addition=mol(1,:);end
end
        mol_updated=mol_updated(~ismember(mol_updated(:,2:4),rand_del_choice(:,2:4),'rows'),:);
        s_mol_updated=size(mol_updated,1);
        
        %%% make the deleted atom type as type for the previously added atom%%
        mol_updated(end,14)=rand_del_choice(1,14);
        mol_updated(end,15)=0.3;
        mol_updated(end,10:13)=1;
        
     


end

%adjust connection tables 6:8

 mol_updated(:,6:8)=0;
 mol_updated(:,1)=1:s_mol_updated;
 for i=1:s_mol_updated
     for j=1:s_mol_updated
         if (fithagors(mol_updated(i,2:4),mol_updated(j,2:4))==1.33) 
             if (~ismember(mol_updated(j,1),mol_updated(i,6:8)))
                 if( ~mol_updated(i,6)>0) mol_updated(i,6)=mol_updated(j,1); elseif (~mol_updated(i,7)>0) mol_updated(i,7)=mol_updated(j,1); elseif (~mol_updated(i,8)>0) mol_updated(i,8)=mol_updated(j,1);end;
             end
         end
     end
 end
         
 for i=1:s_mol_updated; mol_updated(i,5)=numel(find(mol_updated(i,6:8)>0)); end


mol_updated

%%% calcualte charge for mutated mol
 write_mol_xyz(mol_updated,'moh','temp.xyz')
        %copyfile('d:\MATLAB\ttt\temp.xyz','d:\OpenBabel-2.3.2\temp.xyz')
        %cd d:\OpenBabel-2.3.2\
        dos('obabel -i xyz temp.xyz -o mol2 -O temp.mol2 --partialcharge gasteiger')
        %cd d:\MATLAB\ttt
        temp=importdata('temp.mol2','') %d:\OpenBabel-2.3.2\temp.mol2','')
        for i=1:size(mol_updated,1)
        a= sscanf(temp{i+6},'%i %c %f%f%f %*s %i %*s %f');
        charge=a(end);
        mol_updated(i,15)=charge;
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
axis equal;
                         hold off;
                        
                        

end
