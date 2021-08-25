function p = make_sp2_8ed( num_atoms,num_O,num_N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%switch (atom_type)) case 67; atom_type='C'; case 79; atom_type='O'; case 78; atom_type='N'; case 72; atom_type=H; end

i=0;
 while i<num_atoms
     i=i+1;
     answer=1;
     if (i==1) p(1,1:15)=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0.3];
     else
   
     n_ready=p(p(p(:,5)<3));
     s=size(p);
     choice=randsample(n_ready(:,1),1);
     
     if(choice ~=1) origins=p(ismember(p(:,1),p(choice,6:8)),:) ;else origins=[0 1.33 0 0];end
     %generate random angle sign
     theta=120;
     theta2=180;
     if(rand(1)>0.5) theta_sign=-1 ;else theta_sign=1; end
     if(rand(1)>0.5) theta2_sign=1;else theta2_sign=0; end
     s_origins=size(origins,1);n=0;
     while (n<s_origins)
         n=n+1;
         origin=origins(n,2:4)
     %point=rotmat3(origin-p(choice,2:4),[0 0 0],[0 0 1],theta_sign*theta)+p(choice,2:4);
     point=roundn(rotmat3(origin-p(choice,2:4),[0 0 0],[0 0 1],theta_sign*theta)+p(choice,2:4),-4);
     pole=p(choice,2:4)-origin;
    
    % point=rotmat3(point,origin,pole,theta2_sign*theta2);
    
     %if(find(p(:,2)==point(1))*find(p(:,3)==point(2))*find(p(:,4)==point(3))) i=i-1; else  p(i,1:5)=[i,point,2]; end
    for j=1:s(1,1); if (p(j,2:4)==point) j=s(1,1)+1 ; answer=answer*0; else answer=answer*1;end; end
    if(answer) p(choice,5)=p(choice,5)-1;p(i,1:4)=[i,point] ; 
        if( ~p(i,6)>0) p(i,6)=choice; elseif (~p(i,7)>0) p(i,7)=choice; elseif (~p(i,8)>0) p(i,8)=choice;end; 
        p(i,5)=numel(find(p(i,6:8)>0));
    %else i=i-1; 
    end
     
     end
     if (answer==0) i=i-1; end
     end
     
%      %%%%optinoal draw of atoms addition progress start
%      
%      scatter3 (p(:,2),p(:,3),p(:,4));hold on;
% % text([p(:,2:4)],'String', p(:,1), 'horizontal','left', 'vertical','bottom');
% % text('Position',p(:,2:4),'String','C')
%  s=size(p);
%          for i=1:s(1,1); 
%              for j=1:s(1,1); 
%                  if ( roundn(((p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2+(p(i,4)-p(j,4))^2)^0.5,-2) ==1.33)
%                      % ((p(i,1)-p(j,1))^2+(p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2)^0.5 ==bond_length_4  )  
%                   
%                     
% 
%                      plot3([p(i,2);p(j,2)],...
%                             [p(i,3);p(j,3)],...
%                             [p(i,4);p(j,4)],'-b','LineWidth',3) 
%                  else
%                             end;
%                             end
%                             
%                             ;end
%  xlabel('X');
%  ylabel('Y');
%  zlabel('Z');
%  axis([-5, 5, -5, 5, -5, 5]); % axis([xmin, xmax, ymin, ymax])
%  axis equal;
%                         
%                         
%      
%  end
%  hold off;
%  %%%%optinoal draw of atoms addition progress end
 end 
 
 %renumbering atoms
 s_p=size(p,1)
 p(:,6:8)=0;
 p(:,10:13)=1;
 p(:,15)=0.3;
 p(:,1)=1:s_p;
 for i=1:num_atoms
     for j=1:num_atoms
         if (fithagors(p(i,2:4),p(j,2:4))==1.33) 
             if (~ismember(p(j,1),p(i,6:8)))
                 if( ~p(i,6)>0) p(i,6)=p(j,1); elseif (~p(i,7)>0) p(i,7)=p(j,1); elseif (~p(i,8)>0) p(i,8)=p(j,1);end;
             end
         end
     end
 end
         
 for i=1:num_atoms; p(i,5)=numel(find(p(i,6:8)>0)); end
 
 %%% adding atom types%%%
 nums=1:s_p;
 %num_O= round( O_per*s_p);
 %num_N= round( N_per*s_p);
 num_C=s_p-num_O-num_N;
 
%%% Os=randsample(nums,num_O)
%%% Os=randsample(nums,num_O)
%%% Ns=randsample(nums(~ismember(nums,Os)),num_N);
%%% Cs=nums(~ismember(nums,[Os,Ns]))

% 2021 edition to avoid O have three bonds
Os=[0]
Ns=[0]
Cs=[0]

for i=1:num_O
 %Os=randsample(p(p(p(:,5)<3),1),1)   %%random select of atoms to be Os if it has less than 3 connections   
 if size(Os,2) == 1 
     
   Os(i)=randsample(p(p(p(:,5)<3 ,1)),1) % exclude previous Os and their connected atoms
 else
   Os(i)=randsample(p(p(p(:,5)<3 & ~ismember (p(:,1),Os(i)) & ~ismember (p(:,6),Os) & ~ismember (p(:,7),Os) & ~ismember (p(:,8),Os) & ~ismember (p(:,9),Os) ,1)),1) % exclude previous Os and their connected atoms

 end;
  end

for i=1:num_N
 %Os=randsample(p(p(p(:,5)<3),1),1)   %%random select of atoms to be Os if it has less than 3 connections   
 if size(Ns,2) == 1 
     
   Ns(i)=randsample(p(p(p(:,5)<3 & ~ismember (p(:,1),Os) & ~ismember (p(:,6),Os) & ~ismember (p(:,7),Os) & ~ismember (p(:,8),Os) & ~ismember (p(:,9),Os) ,1)),1) % exclude previous Os and their connected atoms
 else
  
   Ns(i)=randsample(p(p(p(:,5)<3 & ~ismember (p(:,1),Os) & ~ismember (p(:,6),Os) & ~ismember (p(:,7),Os) & ~ismember (p(:,8),Os) & ~ismember (p(:,9),Os) & ~ismember (p(:,1),Ns) & ~ismember (p(:,6),Ns) & ~ismember (p(:,7),Ns) & ~ismember (p(:,8),Ns) & ~ismember (p(:,9),Ns) ,1)),1) % exclude previous Os and their connected atoms
 
 end;
  end  ; 

%Ns=randsample(p(~ismember(p(:,1),Os),1),num_N);
%Ns=randsample(p(~ismember(p(:,1),Os) & ~ismember(p(:,6),Os) & ~ismember(p(:,7),Os)& ~ismember(p(:,8),Os)& ~ismember(p(:,9),Os) ,1),num_N); % exclude oxygen atoms and atoms coonected to oxygen atoms
Cs=p(~ismember(p(:,1),[Os,Ns]),1)
 %for i=1:s_p; p(i,14)= randsample(['C' 'O' 'N'],1);end
 p(Os,14)='O';
 p(Ns,14)='N';
 p(Cs,14)='C';
 
 
 
 
% sort atoms according to distance from atom 1 (0,0,0)
% for i=1:num_atoms
%     f(i,1)=i;
%     f(i,2)=fithagors(p(1,2:4),p(i,2:4));
% end
% f_sorted=sortrows(f,2);
% for i=1:num_atoms
%     p_sorted(i,:)=p(f_sorted(i,1),:);p_sorted(i,1)=i;
% end
% p=p_sorted;
     

 % renumbering according to fithagors
 %for i=1:num_atoms
 %    n=0; p(i,5)=3;
 %    for j=1:num_atoms
 %        if(fithagors(p(i,2:4),p(j,2:4))==1.33) p(i,5)=p(i,5)-1;if(n==0) n=j;end;end
 %    end
 %    p(i,6)=p(n,1);
 %end

 %%% calcuation of atomic charge
        write_mol_xyz(p,'moh','temp.xyz')
        %copyfile('d:\MATLAB\ttt\temp.xyz','d:\OpenBabel-2.3.2\temp.xyz')
        %cd d:\OpenBabel-2.3.2\
        dos('obabel -i xyz temp.xyz -o mol2 -O temp.mol2 --partialcharge gasteiger')
        %cd d:\MATLAB\ttt
        temp=importdata('temp.mol2','') %d:\OpenBabel-2.3.2\temp.mol2','')
        for i=1:size(p,1)
        a= sscanf(temp{i+6},'%i %c %f%f%f %*s %i %*s %f')
        charge=a(end);
        p(i,15)=charge;
        end
p
 scatter3 (p(:,2),p(:,3),p(:,4));hold on;
% text([p(:,2:4)],'String', p(:,1), 'horizontal','left', 'vertical','bottom');
% text('Position',p(:,2:4),'String','C')
 s=size(p);
         for i=1:s(1,1); 
             for j=1:s(1,1); 
                 if ( roundn(((p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2+(p(i,4)-p(j,4))^2)^0.5,-2) ==1.33)
                     % ((p(i,1)-p(j,1))^2+(p(i,2)-p(j,2))^2+(p(i,3)-p(j,3))^2)^0.5 ==bond_length_4  )  
                  
                     plot3([p(i,2);p(j,2)],...
                            [p(i,3);p(j,3)],...
                            [p(i,4);p(j,4)],'-b','LineWidth',3) 
                 else
                            end;
                            end
                            
                            ;end
 xlabel('X');
 ylabel('Y');
 zlabel('Z');
 axis equal;
                        hold off;

end

