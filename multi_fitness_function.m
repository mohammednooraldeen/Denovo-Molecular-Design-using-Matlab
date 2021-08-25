function e_map = multi_fitness_function(r,crystal_with_e)
%FITNESS_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
total_fitha=0;
total_energy=0;
s=size(r);
e_map=zeros(2)
for i=1:s(1,2)/2
total_fitha= total_fitha+sum(fithagors(crystal_with_e(r(i),1:3),crystal_with_e(r(1:s(1,2)/2),1:3)));
end

connectivity=0;
for i=1:(s(1,2)/2)
    for j=1:(s(1,2)/2)
   if(roundn(fithagors(crystal_with_e(r(i),1:3),crystal_with_e(r(j),1:3)),-2)==1.33) connectivity =connectivity+1; j=(s(1,2)/2)+1;
   
   ; end ; end
 if  (connectivity<1) i=1+s(1,2)/2; end;
    
end;
connectivity=connectivity-s(1,2)/2;

for i=1+s(1,2)/2:s(1,2)
total_energy=total_energy+ crystal_with_e(r(i-s(1,2)/2),r(i));
end;



e_map=[total_fitha; total_energy;connectivity];



end

