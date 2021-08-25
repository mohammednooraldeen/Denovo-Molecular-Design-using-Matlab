function [ coord ] = make_crystal_2( num_atoms_x,num_atoms_y, num_atoms_z,all)
%MAKE_CRYSTAL Summary of this function goes here
%   Detailed explanation goes here
bond_length_1=1.33;
count=1;
 hex=  [0 0;
     1.33*sind(60)    1.33*cosd(60);
     1.33*sind(60)    1.33*cosd(60)+1.33;
     0       1.33*cosd(60)+1.33*cosd(60)+1.33];
for m=0:num_atoms_z
    for i=0:num_atoms_x
      for j=0:num_atoms_y
          if(roundn(m/2,0)==roundn(m/2,-1)) shift=1 ;else shift=0;end
           coord(count:count+3,1:2)=[hex(1:4,1)+i*2*1.33*sind(60)+shift*1.33*sind(60) , hex(1:4,2)+j*2*(1.33*cosd(60)+1.33)];
           coord(count:count+3,3)=m*1.33*sind(60) 
           count=count+4;
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




 

hold off;
s=size(coord)

for i=1:s(1,1)
       a= fithagors (all(:,1:3),coord(i,1:3))
       a(:,2:7)=all(:,4:9);
       a=sortrows(a,-1);
       coord(i,4:9)=a(end,2:7);   
       
end

    
scatter3(coord(:,1),coord(:,2),coord(:,3),100,coord(:,4),'fill'); hold on;

s=size(coord);
         for i=1:s(1,1); 
             for j=1:s(1,1); 
                 if ( roundn(((coord(i,1)-coord(j,1))^2+(coord(i,2)-coord(j,2))^2+(coord(i,3)-coord(j,3))^2)^0.5,-2) ==1.33)
                     % ((coord(i,1)-coord(j,1))^2+(coord(i,2)-coord(j,2))^2+(coord(i,3)-coord(j,3))^2)^0.5 ==bond_length_4  )  
                  
                     plot3([coord(i,1);coord(j,1)],...
                            [coord(i,2);coord(j,2)],...
                            [coord(i,3);coord(j,3)],'-b','LineWidth',3) 
                 else
                            end;
                            end
                            
                            ;end

end

