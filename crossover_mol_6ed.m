function child = crossover_mol_6ed( mol1,mol2,energy_weight,C,O,N,H,e,P_C,P_O,P_N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% It is recommended for this function to get the parent molecules with
% atomic energy recorded at column 15, So direct selection and crossover
% can be made. without the need for recalucation of atomic charge and
% energy by eval_xyz_sub_surfaces_3 function. However, this mean that the
% fitness function should be able to modify the evaluated gene by inserting
% the charge or energy value at column-15 which is impossible. To take
% over, the creation function can generate mol and calculate charge for it
% to make it available for fitness or crossover functions.



s_mol1=size(mol1,1);
%         %%%adding charge to mol1
%         write_mol_xyz(mol1,'moh','temp.xyz')
%         %copyfile('d:\MATLAB\ttt\temp.xyz','d:\OpenBabel-2.3.2\temp.xyz')
%         %cd d:\OpenBabel-2.3.2\
%         dos('obabel -i xyz temp.xyz -o mol2 -O temp.mol2 --partialcharge gasteiger')
%         %cd d:\MATLAB\ttt
%         temp=importdata('temp.mol2','') %d:\OpenBabel-2.3.2\temp.mol2','')
%         for i=1:s_mol1
%         a= sscanf(temp{i+7},'%i %c %f%f%f %*s %i %*s %f');
%         charge=a(end);
%         mol1(i,15)=charge;
%         end
%%%add column 16 for parent number
mol1(:,15)=1;
%%%change column 15 from charge to energy impact
%for i=1:s_mol1; mol1(i,15)=eval_xyz_sub_surfaces_4(mol1(i,:),'C',C,O,N,H,e,P_C,P_O,P_N,'true');end


s_mol2=size(mol2,1);
%         %%%adding charge to mol2
%         write_mol_xyz(mol2,'moh','temp.xyz')
%         %copyfile('d:\MATLAB\ttt\temp.xyz','d:\OpenBabel-2.3.2\temp.xyz')
%         %cd d:\OpenBabel-2.3.2\
%         dos('obabel -i xyz temp.xyz -o mol2 -O temp.mol2 --partialcharge gasteiger')
%         %cd d:\MATLAB\ttt
%         temp=importdata('temp.mol2','') %d:\OpenBabel-2.3.2\temp.mol2','')
%         for i=1:s_mol2
%         a= sscanf(temp{i+7},'%i %c %f%f%f %*s %i %*s %f');
%         charge=a(end);
%         mol2(i,15)=charge;
%         end

%%%add column 16 for parent number
mol2(:,15)=2;
%%%change column 15 from charge to energy impact
%for i=1:s_mol2; mol2(i,15)=eval_xyz_sub_surfaces_4(mol2(i,:),'C',C,O,N,H,e,P_C,P_O,P_N,'true');end


mol1
mol2

link_map(1,3)=100;
link_map=linkage([mol1(:,2:4); mol2(:,2:4)], 'single')

dissimilarity=numel(find(ismember(mol1(:,2:4),mol2(:,2:4),'rows')==0))
%if (dissimilarity>0) mol1_identical_to_mol2=0; else mol1_identical_to_mol2=1;end
if (dissimilarity==0) mol_updated=mol1;
else
    main_iteration=0;
    repeat_atom_distribution=1;
    O_assigned=0;
    N_assigned=0;
    while (main_iteration <20 && (O_assigned==0 | N_assigned==0))
        "###########################new___Cross#########################"
        link_hir(1,3)=100;identical_to_mol1=1; identical_to_mol2=1;
        num_iteration=0;
        while (mean(link_hir(:,3))>1.33 | identical_to_mol1==1 | identical_to_mol2==1 | num_iteration < 30 )
            num_iteration=num_iteration+1;
            mol1_fusion_part=mol1(1,:);
            s_mol1_fusion_part=size(mol1_fusion_part,1);
            mol1_pool=mol1(~ismember(mol1(:,2:4),mol1_fusion_part(:,2:4),'rows'),:);
            s_mol1_pool=size(mol1_pool);s_mol1_pool=s_mol1_pool(1,1);
            
            min_num_atoms_mol1=roundn(s_mol1*0.25,0)
            max_num_atoms_mol1=roundn(s_mol1*0.75,0)
            num_atoms_mol1=randsample(min_num_atoms_mol1:max_num_atoms_mol1,1)
            
            w=round(1+((max(mol1_pool(:,15))-mol1_pool(:,15))/(max(mol1_pool(:,15))-min(mol1_pool(:,15))))*100*energy_weight)
            if (isnan(w)) ww=[1:s_mol1_pool];
            else
                ww=0; for i=1:s_mol1_pool; ww(end:end+w(i))=i;end  %% creating vector of atom numbers,each atom number is repeated in the vector
            end                                                %% by times equal to its energy impact (i.e. low energy high frequency of repetition)
            
            for i=1:num_atoms_mol1;
                if (numel(ww)>1) rand_ww=randsample(ww,1); else rand_ww=ww;end
                rand_atom_mol1=mol1_pool(rand_ww,:);     % choose atom from vector
                ww=ww(~ismember(ww,rand_ww))             %remove selected atom from vector
                mol1_fusion_part(i+1,:)= rand_atom_mol1;
            end
            
            s_mol1_fusion_part=size(mol1_fusion_part,1);
            
            mol2_pool=mol2(~ismember(mol2(:,2:4),mol1_fusion_part(:,2:4),'rows'),:)
            num_atoms_mol2=s_mol2-s_mol1_fusion_part;
            s_mol2_pool=size(mol2_pool,1);
            rand_atoms_mol2=zeros(1,15);
            w=round(1+((max(mol2_pool(:,15))-mol2_pool(:,15))/(max(mol2_pool(:,15))-min(mol2_pool(:,15))))*100*energy_weight)
            if (isnan(w)) ww=[1:s_mol2_pool];
            else
                ww=0; for i=1:s_mol2_pool; ww(end:end+w(i))=i;end
            end
            for i=1:num_atoms_mol2;
                if (numel(ww)>1) rand_ww=randsample(ww,1); elseif (numel(ww)==1) rand_ww=ww;end
                rand_atom_mol2=mol2_pool(rand_ww,:);
                ww=ww(~ismember(ww,rand_ww))
                %if (size(rand_atoms_mol2,1)==1) rand_atoms_mol2(1,:)= rand_atom_mol2; else rand_atoms_mol2(end+1,:)= rand_atom_mol2;end
                rand_atoms_mol2(i,:)=rand_atom_mol2;
            end
            
            s_rand_atoms_mol2=size(rand_atoms_mol2,1);
            
            mol_updated=mol1_fusion_part
            mol_updated(end+1:end+s_rand_atoms_mol2,:)=rand_atoms_mol2
            
            
            %check for molecule break
            link_hir=linkage(mol_updated(:,2:4),'single')
            
            %check if exactly similar to one parent
            if (find(ismember(mol1(:,2:4),mol_updated(:,2:4),'rows')==0)) identical_to_mol1=0 ; else identical_to_mol1=1; end
            if (find(ismember(mol2(:,2:4),mol_updated(:,2:4),'rows')==0)) identical_to_mol2=0 ; else identical_to_mol2=1; end
            
            
            
            
            
            %          link_heteros(1,3)=1.33
            %          while (mean(link_heteros<=1.33))
            %
            %          link_heteros=linkage(mol_updated(mol_updated(:,14)==78 | mol_updated(:,14)==79,2:4),'single') %check for adjacent heteros
            %
            %          end
            
        end
        
        
    
    
%     repeat_atom_distribution_iteration=1
%     while (repeat_atom_distribution_iteration<200 )
        s_mol_updated=size(mol_updated,1);
        mol_updated(:,1)=1:s_mol_updated;
        mol_updated(:,15:16)=0;
        for i=1:s_mol_updated
            atom1_type= mol1(ismember(mol1(:,2:4),mol_updated(i,2:4),'row'),14)
            if (atom1_type)
                mol_updated(i,14)=atom1_type;
            else
            end
            
            atom2_type= mol2(ismember(mol2(:,2:4),mol_updated(i,2:4),'row'),14)
            if (atom2_type)
                mol_updated(i,15)=atom2_type;
            else
            end
            
            
        end
        % distribute atom types
        mol1_O_num=numel(find(mol1(:,14)==79));
        mol1_N_num=numel(find(mol1(:,14)==78));
        mol1_C_num=numel(find(mol1(:,14)==67));
        
        mol_updated_O_vector= mol_updated(mol_updated(:,14)==79 | mol_updated(:,15)==79,1)
        mol_updated_N_vector= mol_updated(mol_updated(:,14)==78 | mol_updated(:,15)==78,1)
        mol_updated_heteros_vector= [mol_updated_O_vector;mol_updated_N_vector];
      
        %%%choose distal heteros that are non-adjacent heteros possible positions
        link_heteros= linkage(mol_updated(mol_updated_heteros_vector,2:4),'single');
        mol_updated_heteros_vector_linkage_distal=[link_heteros(link_heteros(:,3)>1.33,1); link_heteros(link_heteros(:,3)>1.33,2)];  % remove close atom ranks
        mol_updated_heteros_vector_linkage_distal=mol_updated_heteros_vector_linkage_distal(mol_updated_heteros_vector_linkage_distal(:,1)<=size(mol_updated_heteros_vector,1),:); % remove new rank numbers for dentries i.e. larger than max atoms intered to linakge function
        mol_updated_heteros_vector_distal=mol_updated_heteros_vector(unique(mol_updated_heteros_vector_linkage_distal(:,1)),1);  % convert vector back from ranks to atom numbers as in mol_updated file
        mol_updated_heteros_vector_distal
        
        %%%choose carbons distal to distal heteros to be spared if heteros
        %%%vector is not sufficient for mol1_O_num+mol1_N_num
        mol_updated_C_vector= mol_updated(~ismember(mol_updated(:,1),mol_updated_heteros_vector_distal),1)
        mol_updated_heteros_carbons_vector=[mol_updated_heteros_vector_distal;mol_updated_C_vector];
        link_heteros_carbons= linkage(mol_updated(mol_updated_heteros_carbons_vector,2:4),'complete');
        mol_updated_heteros_carbons_vector_linkage_distal=[link_heteros_carbons(link_heteros_carbons(:,3)>1.33,1); link_heteros_carbons(link_heteros_carbons(:,3)>1.33,2)];
        mol_updated_heteros_carbons_vector_linkage_distal=mol_updated_heteros_carbons_vector_linkage_distal(mol_updated_heteros_carbons_vector_linkage_distal(:,1)<=size(mol_updated_heteros_carbons_vector,1),:);
        mol_updated_carbons_vector_distal=mol_updated_heteros_carbons_vector(unique(mol_updated_heteros_carbons_vector_linkage_distal(:,1)),1);
        mol_updated_carbons_vector_distal=mol_updated_carbons_vector_distal(~ismember(mol_updated_carbons_vector_distal,mol_updated_heteros_vector_distal),:)
        
        mol_updated_heteros_carbons_vector_distal=[mol_updated_heteros_vector_distal; mol_updated_carbons_vector_distal]
        mol_updated_others_vector=mol_updated(~ismember(mol_updated(:,1),[mol_updated_heteros_vector_distal; mol_updated_carbons_vector_distal]),1)
        
        %mol_updated_hetero_sorted= mol_updated(sortrows(mol_updated,16,'descend')
        %%%%%%%%%%%%%%%%%%%%%%sorting vector and assign oxygen atoms
        mol_updated_heteros_vector_distal_O=mol_updated_heteros_vector_distal(mol_updated(mol_updated_heteros_vector_distal,14)==79 | mol_updated(mol_updated_heteros_vector_distal,15)==67 | mol_updated(mol_updated_heteros_vector_distal,14)==67 | mol_updated(mol_updated_heteros_vector_distal,15)==79,:)
        temp1=sortrows(mol_updated(mol_updated_heteros_vector_distal_O,:),[14:15],'descend')
        mol_updated_heteros_vector_distal_O_sorted=temp1(:,1);
        O_uaza=0;
        if (numel(mol_updated_heteros_vector_distal_O_sorted)< mol1_O_num );
            O_uaza=mol1_O_num- size(mol_updated_heteros_vector_distal_O_sorted,1)
            if (numel(mol_updated_heteros_carbons_vector_distal)>=O_uaza)
                
                mol_updated_heteros_vector_distal_O_sorted=[mol_updated_heteros_vector_distal_O_sorted; mol_updated_heteros_carbons_vector_distal(1:O_uaza)]
            else
            end
            
        else
        end
        
        %check for assignment
        if size(mol_updated_heteros_vector_distal_O_sorted,1)>=mol1_O_num;
            mol_updated(mol_updated_heteros_vector_distal_O_sorted(1:mol1_O_num),16)=79  %assign oxygens
            mol_updated_heteros_vector_distal_O_assigned=mol_updated(mol_updated_heteros_vector_distal_O_sorted(1:mol1_O_num),1);
            O_assigned=1;
        else
            O_assigned=0;
        end
        
        
        if (O_assigned==1)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%sorting vector and assign nitrogen atoms
            %mol_updated_heteros_vector_distal_N=mol_updated_heteros_vector_distal_N(mol_updated(mol_updated_heteros_vector_distal,14)==78 | mol_updated(mol_updated_heteros_vector_distal,15)==67 | mol_updated(mol_updated_heteros_vector_distal,14)==67 | mol_updated(mol_updated_heteros_vector_distal,15)==78,:)
            mol_updated_heteros_vector_distal_N=mol_updated_heteros_carbons_vector_distal(~ismember(mol_updated_heteros_carbons_vector_distal,mol_updated_heteros_vector_distal_O_assigned),:) % all heteros distal except those assigned to oxygens
            
            temp1=sortrows(mol_updated(mol_updated_heteros_vector_distal_N,:),[14:15],'descend')
            mol_updated_heteros_vector_distal_N_sorted=temp1(:,1);
            
            %check for assignment
            if (size(mol_updated_heteros_vector_distal_N_sorted,1)>=mol1_N_num && O_assigned==1 );
                
                mol_updated(mol_updated_heteros_vector_distal_N_sorted(1:mol1_N_num),16)=78  %assign nitrogens
                mol_updated_heteros_vector_distal_N_assigned=mol_updated(mol_updated_heteros_vector_distal_N_sorted(1:mol1_N_num),1);
                N_assigned=1;
            else
                N_assigned=0;
            end
            
            mol_updated(mol_updated(:,16)==0,16)=67   %assign carbons
        end
        
        
    %repeat_atom_distribution==0;
      
    main_iteration=main_iteration+1;

    end 
     %make parent as child if no cross over can be made
     if (O_assigned==0 | N_assigned==0) 
         mol_updated=mol1
         
     else
         mol_updated(:,14)=mol_updated(:,16);
         mol_updated(:,15:16)=[];
     end
    
    %adjust connection tables 6:8
    s_mol_updated=size(mol_updated,1);
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

end
child=mol_updated;

%%%%calcualte atomic charge for the child
 write_mol_xyz(child,'moh','temp.xyz')
        %copyfile('d:\MATLAB\ttt\temp.xyz','d:\OpenBabel-2.3.2\temp.xyz')
        %cd d:\OpenBabel-2.3.2\
        dos('obabel -i xyz temp.xyz -o mol2 -O temp.mol2 --partialcharge gasteiger')
        %cd d:\MATLAB\ttt
        temp=importdata('temp.mol2','') %d:\OpenBabel-2.3.2\temp.mol2','')
        for i=1:size(child,1)
        a= sscanf(temp{i+6},'%i %c %f%f%f %*s %i %*s %f');
        charge=a(end);
        child(i,15)=charge;
        end


%scatter3 (mol1(:,2),mol1(:,3),mol1(:,4),20,[1 0 0]);hold on;
%scatter3 (mol2(:,2),mol2(:,3),mol2(:,4),20,[0 1 0]);hold on;
%scatter3 (child(:,2),child(:,3),child(:,4),20,[0 1 0]);hold on;

%%%%color atoms of parents and child
for i=1:size(child,1)  
scatter3 (mol1(i,2),mol1(i,3),mol1(i,4), 300,atom_type_color(char(mol1(i,14))));hold on;
scatter3 (mol2(i,2),mol2(i,3),mol2(i,4),500,atom_type_color(char(mol2(i,14))));hold on;
scatter3 (child(i,2),child(i,3),child(i,4),700,atom_type_color(char(child(i,14))));hold on;
%'MarkerEdgeColor','k'
end


%scatter3 (mol1(:,2),mol1(:,3),mol1(:,4),20,char(mol1(:,14)));hold on;
%scatter3 (mol2(:,2),mol2(:,3),mol2(:,4),20,char(mol1(:,14)));hold on;
%scatter3 (child(:,2),child(:,3),child(:,4),100,char(mol1(:,14)));hold on;

mol1

s=size(mol1);
         for i=1:s(1,1); 
             for j=1:s(1,1); 
                 if ( roundn(((mol1(i,2)-mol1(j,2))^2+(mol1(i,3)-mol1(j,3))^2+(mol1(i,4)-mol1(j,4))^2)^0.5,-2) ==1.33)
                     % ((mol1(i,1)-mol1(j,1))^2+(mol1(i,2)-mol1(j,2))^2+(mol1(i,3)-mol1(j,3))^2)^0.5 ==bond_length_4  )  
                  
                     plot3([mol1(i,2);mol1(j,2)],...
                            [mol1(i,3);mol1(j,3)],...
                            [mol1(i,4);mol1(j,4)],'-b','LineWidth',7) 
                 else
                            end;
                            end
                            
                            ;end ;hold on;
 s=size(mol2);
         for i=1:s(1,1); 
             for j=1:s(1,1); 
                 if ( roundn(((mol2(i,2)-mol2(j,2))^2+(mol2(i,3)-mol2(j,3))^2+(mol2(i,4)-mol2(j,4))^2)^0.5,-2) ==1.33)
                     % ((mol2(i,1)-mol2(j,1))^2+(mol2(i,2)-mol2(j,2))^2+(mol2(i,3)-mol2(j,3))^2)^0.5 ==bond_length_4  )  
                  
                     plot3([mol2(i,2);mol2(j,2)],...
                            [mol2(i,3);mol2(j,3)],...
                            [mol2(i,4);mol2(j,4)],'-r','LineWidth',4) 
                 else
                            end;
                            end
                            
                            ;end ;hold on;
                        
s=size(child);
         for i=1:s(1,1); 
             for j=1:s(1,1); 
                 if ( roundn(((child(i,2)-child(j,2))^2+(child(i,3)-child(j,3))^2+(child(i,4)-child(j,4))^2)^0.5,-2) ==1.33)
                     % ((child(i,1)-child(j,1))^2+(child(i,2)-child(j,2))^2+(child(i,3)-child(j,3))^2)^0.5 ==bond_length_4  )  
                  
                   plot3([child(i,2);child(j,2)],...
                            [child(i,3);child(j,3)],...
                            [child(i,4);child(j,4)],'.-g','LineWidth',2) 
                 else
                            end;
                            end
                            
                            ;end                         
 xlabel('X');
 ylabel('Y');
 zlabel('Z');
 view([0,0,50]);
 axis equal;
                        hold off;



 



end

