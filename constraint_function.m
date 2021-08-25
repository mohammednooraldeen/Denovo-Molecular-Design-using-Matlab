function [c, cequ] = constraint_function(r,crystal_with_e)
%CONSTRAINT_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
%atom types are 4-->8   C-N-O-H-P
s=size(r);count=0;diversity=0;composition=0;total_fitha=0;
for i=1:s(1,2)/2
    for j=1:s(1,2)/2
        if (r(i)==r(j)) count=count+1; end;
    end
end

for i=1+s(1,2)/2:s(1,2)
    for j=1+s(1,2)/2:s(1,2)
        if (r(i)==r(j)) diversity=diversity+1; end;
    end
end



c=0;n=0;o=0;h=0;p=0;
for i=1+s(1,2)/2:s(1,2)
    if (r(i)==4) c=c+1;
    elseif (r(i)==5) n=n+1;
     elseif (r(i)==6) o=o+1;  
         elseif (r(i)==7) h=h+1;
             elseif (r(i)==8) p=p+1;
    end;
end;
total_types=0;
total_types=c+n+o+h+p;
composition=-1*(2*c/total_types +0.01* n/total_types+ 0.01* o/total_types+ 0.01*h/total_types+0.01*p/total_types) 
0.8-(c/total_types)

connectivity=100; coun=0;connec=[0 0];
for i=1:(s(1,2)/2)
   
    for j=1:(s(1,2)/2)
        coun=coun+1;
    if(roundn(fithagors(crystal_with_e(r(i),1:3),crystal_with_e(r(j),1:3)),-2)==1.33); connec(coun,1)=r(i); connec(coun,2)=r(j);end
   
    ; end
    
end;
% now analyse the connect matrix
s_cx=size(connec); s_cx=s_cx(1,1);
s_cy=size(connec); s_cy=s_cy(1,2);

if (s_cx == 2*(-1+s(1,2)/2) || s_cx == 2*(s(1,2)/2)) connectivity=0; else connectivity=10000000; end


%connectivity=connectivity-(s(1,2)/2)+1


for i=1:s(1,2)/2
total_fitha= total_fitha+sum(fithagors(crystal_with_e(r(i),1:3),crystal_with_e(r(1:s(1,2)/2),1:3)));
end

c=[ count;total_fitha;connectivity; composition; abs((0.8- c)/total_types); abs((0.1-n)/total_types);(h/total_types); (p/total_types); abs((0.1-o)/total_types)];

cequ=[];
end

