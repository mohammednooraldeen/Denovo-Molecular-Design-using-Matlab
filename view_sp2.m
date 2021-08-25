function [ output_args ] = view_sp2( p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
scatter3 (p(:,2),p(:,3),p(:,4));hold on;
 
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
 view([0 0 20])
 axis equal
 hold off;

end

