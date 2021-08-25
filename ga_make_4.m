function [mol,score]= ga_make_4( num_atoms,O_num, N_num,num_sub_boxes, C_map,O_map,N_map,H_map,e_map )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

C=importdata(C_map,'',6); 
O=importdata(O_map,'',6); 
N=importdata(N_map,'',6); 
H=importdata(H_map,'',6); 
e=importdata(e_map,'',6); 
centre=sscanf(C.textdata{6,1},'%*[CENTER]%f%f%f');
elements=sscanf(C.textdata{5,1},'%*[NELEMENTS]%i%i%i');
spacing=sscanf(C.textdata{4,1},'%*[SPACING]%f');

x_surfaces=[centre(1)-spacing*elements(1)/2, centre(1)+spacing*elements(1)/2]
y_surfaces=[centre(2)-spacing*elements(2)/2, centre(2)+spacing*elements(2)/2]
z_surfaces=[centre(3)-spacing*elements(3)/2, centre(3)+spacing*elements(3)/2]

x_dist=abs(x_surfaces(2)-x_surfaces(1)); x_dist_frac=x_dist/num_sub_boxes(1);
y_dist=abs(y_surfaces(2)-y_surfaces(1)); y_dist_frac=y_dist/num_sub_boxes(2);
z_dist=abs(z_surfaces(2)-z_surfaces(1)); z_dist_frac=z_dist/num_sub_boxes(3);

for i=0:num_sub_boxes(1)
    x_sub_surfaces(1,i+1)=x_surfaces(1)+i*x_dist_frac;
end
x_sub_surfaces
for i=0:num_sub_boxes(2)
    y_sub_surfaces(1,i+1)=y_surfaces(1)+i*y_dist_frac;
end
y_sub_surfaces
for i=0:num_sub_boxes(3)
    z_sub_surfaces(1,i+1)=z_surfaces(1)+i*z_dist_frac;
end
z_sub_surfaces

sub_elements(1)=roundn(elements(1)/num_sub_boxes(1),0)
sub_elements(2)=roundn(elements(2)/num_sub_boxes(2),0)
sub_elements(3)=roundn(elements(3)/num_sub_boxes(3),0)
count=0;

for i=1:num_sub_boxes(3)
    for j=1:num_sub_boxes(2)
        for k=1:num_sub_boxes(1)
            count=count+1;
            x_sub_surface=x_sub_surfaces(1,k:k+1); sub_centre(count,1)= ((x_sub_surface(2)-x_sub_surface(1))/2)+x_sub_surface(1)
            y_sub_surface=y_sub_surfaces(1,j:j+1); sub_centre(count,2)= ((y_sub_surface(2)-y_sub_surface(1))/2)+y_sub_surface(1) 
            z_sub_surface=z_sub_surfaces(1,i:i+1); sub_centre(count,3)= ((z_sub_surface(2)-z_sub_surface(1))/2)+z_sub_surface(1)
        
        
     %position=sub_centre(i,:);
     position_mut_rate_parameters(1:2,1)=x_sub_surface';
     position_mut_rate_parameters(1:2,2)=y_sub_surface';
     position_mut_rate_parameters(1:2,3)=z_sub_surface'; %mutation limits
     position_mut_rate_parameters(3,1)=0.2 %  <-----mutation rate
     lb=position_mut_rate_parameters(1,:)
     ub=position_mut_rate_parameters(2,:)
     
     [x,feval]=ga(@(gene) fitness_gene_sub_boxes_vectorized_3(gene,C,N,O,H,e,[],[],[],'false'),1,[],[],[],[],[],[],[],[],gaoptimset('PopulationType','custom','PopulationSize',150,'CreationFcn',...
      @(a,b,c) creating_gene4 (a,b,c,num_atoms,sub_centre(count,:),ub,lb,O_num,N_num),'SelectionFcn',{@selectiontournament,5},'CrossoverFcn',@ (a1,a2,a3,a4,a5,a6) crossover_gene7(a1,a2,a3,a4,a5,a6,0,C,O,N,H,e,[],[],[]),'MutationFcn',{@mutation_gene5,position_mut_rate_parameters},'Vectorized','off',...
      'PlotFcn',{@gaplotbestf,@gaplotselection,@gaplotstopping,@gaplotscores},'Elitecount',7,'CrossoverFraction',0.8,'Generations',30,'UseParallel','Always'))
    
  mol{count}=x{1,1};
  score(count)=feval;   
            
            
  mol_t=translate_mol_gene(mol{count})
  filename=sprintf('zainab%d.xyz', count);
  mol_name=sprintf('mohammed_mol%d_%f', count, score(count));
  write_mol_xyz(mol_t,mol_name,filename);          
  
            
        end
    end
end        


 


end

