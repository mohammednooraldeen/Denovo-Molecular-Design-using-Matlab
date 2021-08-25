function [ coord ] = make_crystal( num_atoms_x,num_atoms_y, num_atoms_z)
%MAKE_CRYSTAL Summary of this function goes here
%   Detailed explanation goes here
bond_length_1=1.33;
bond_length_2=(bond_length_1*cosd(30));
bond_length_3=bond_length_1-bond_length_2;
bond_length_4=bond_length_1*sind(30);
count=0;
for m=0:num_atoms_z
for i=1:num_atoms_x
       incr=0;k=-2;t=-1;co=0;
    for j=1:num_atoms_y
        k=k+1; count=count+1;t=t+1;
         if(k>0)  incr=1;;else incr=0 ;end; 
         if(roundn(t/2,0)==roundn(t/2,-1))incr_2=1; co=co+1; else incr_2=0; end
         if (roundn(m/2,0)==roundn(m/2,-1)) mult=1; else mult=0; end
       coord(count,:)= [bond_length_2*i*2+incr*bond_length_2,... %mult*bond_length_2,...
           bond_length_1*j-bond_length_3*co  ...
           ,     bond_length_2*m];
        
        if (k==2)k=-2; end;
        end
        
end
end

        
        
        
coord_centre= median(coord);
x_shift=coord_centre(1,1)-3.669;
y_shift=coord_centre(1,2)-12.171;
z_shift=coord_centre(1,3)- -2.829;
coord(:,1)=coord(:,1)-x_shift;
coord(:,2)=coord(:,2)-y_shift;
coord(:,3)=coord(:,3)-z_shift;


scatter3(coord(:,1),coord(:,2),coord(:,3),40,[1 0 0],'fill'); hold on;

s=size(coord);
         for i=1:s(1,1); 
             for j=1:s(1,1); 
                 if ( ((coord(i,1)-coord(j,1))^2+(coord(i,2)-coord(j,2))^2+(coord(i,3)-coord(j,3))^2)^0.5 ==1.33 &...
                      ((coord(i,1)-coord(j,1))^2+(coord(i,2)-coord(j,2))^2+(coord(i,3)-coord(j,3))^2)^0.5 ==1.33 )  
                  
                     plot3([coord(i,1);coord(j,1)],...
                            [coord(i,2);coord(j,2)],...
                            [coord(i,3);coord(j,3)],'-b') 
                 else
                            end;
                            end
                            
                            ;end 

hold off;

end

