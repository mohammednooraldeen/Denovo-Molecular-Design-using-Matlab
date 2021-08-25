

num_atoms=30;
num_atoms=num_atoms*2;
lb(1:num_atoms/2)=1;
lb(1+num_atoms/2:num_atoms)=4;
ub(1:num_atoms/2)=108;
ub(1+num_atoms/2:num_atoms)=6;

 %[x,fvals]=ga(  @(r) fitness_function(r,crystal_with_e_2),8,[],[],[],[],[1 1 1 1 4 4 4 4],[108 108 108 108 8 8 8 8 ],@constraint_function,[1:8])%, gaoptimset('PlotFcns',{@gaplotbestf}))
 
 A= [1 -1 0 0 0 0 0 0 0 0]; b=[1] 
 %[x,fvals]=ga(  @(r) fitness_function(r,crystal_with_e_2),num_atoms,[],[],[],[],lb,ub,@(r) constraint_function (r,crystal_with_e_2),[1:num_atoms],...
 %    gaoptimset('Display','iter','PlotFcns',{@gaplotbestf,@gaplotstopping,@gaplotbestindiv,@gaplotscores},...
 %    'PopulationSize',150, 'Generations',150,'StallTimeLimit',150,'StallGenLimit',150,'TimeLimit',150,...
 %    'CreationFcn',@gacreationuniform,'MutationFcn',{@mutationgaussian},'CrossoverFcn',{@crossoversinglepoint},...
 %    'UseParallel','always','Vectorized','off'))
 %[x,fvals]=gamultiobj(@(r) multi_fitness_function(r,crystal_with_e_2),num_atoms,[],[],[],[],lb,ub)
 
 hold off;
 s=size(x); s=s(1,2); 
% scatter3(crystal_with_e_2(x(1:s/2),1), crystal_with_e_2(x(1:s/2),2), crystal_with_e_2(x(1:s/2),3),50,'fill') ...
   %  ,crystal_with_e_2(x(1:s/2),x(1+s/2:s)),'fill')
t=1

for i=1:s/2;
    if (x(i+s/2)==4) marker='C'; color='green';
    elseif (x(i+s/2)==5) marker='N'; color='blue';
         elseif (x(i+s/2)==6) marker='O'; color='red';
              elseif (x(i+s/2)==7) marker='H'; color='cyan';
                   elseif (x(i+s/2)==8) marker='P'; color='magenta';end
        
        
%text(crystal_with_e_2(x(i),1), crystal_with_e_2(x(i),2), crystal_with_e_2(x(i),3),marker,'FontSize',12,'Color',color)
end
  
 hold on;
 
 scatter3(crystal_with_e_2(:,1),crystal_with_e_2(:,2),crystal_with_e_2(:,3),40, crystal_with_e_2(:,4),'fill');
 
 for i=1:s/2; 
             for j=1:s/2; 
                 if ( roundn( ((crystal_with_e_2(x(i),1)-crystal_with_e_2(x(j),1))^2+...
                         (crystal_with_e_2(x(i),2)-crystal_with_e_2(x(j),2))^2+...
                         (crystal_with_e_2(x(i),3)-crystal_with_e_2(x(j),3))^2)^0.5,-2) ==1.33)
                     % ((coord(i,1)-coord(j,1))^2+(coord(i,2)-coord(j,2))^2+(coord(i,3)-coord(j,3))^2)^0.5 ==bond_length_4  )  
                  
                    % plot3([crystal_with_e_2(x(i),1); crystal_with_e_2(x(j),1)],...
                     %    [crystal_with_e_2(x(i),2); crystal_with_e_2(x(j),2)],...
                      %  [ crystal_with_e_2(x(i),3); crystal_with_e_2(x(j),3)],'-b','LineWidth',3) 
                 else
                 end
                            end
 end
 
                        
 %x= 10; scatter3 (all(1:x:end,1),all(1:x:end,2),all(1:x:end,3),20,all(1:x:end,4)); 
 hold on; scatter3(receptor(:,2),receptor(:,3),receptor(:,4),10,'fill'); hold off                   
                        
                        
                 