function energy = eval_xyz_sub_surfaces_4( mol,Atom_Type,C,O,N,H,e,P_C,P_O,P_N,keep_original_charge)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
num_atoms=size(mol,1);energy=0;
%C=importdata('receptor.C.map','',6); %all(:,1)=C.data

%  C=importdata(C_map,'',6);
%  O=importdata(O_map,'',6);
%  N=importdata(N_map,'',6);
%  H=importdata(H_map,'',6);
%  e=importdata(e_map,'',6);

centre=sscanf(C.textdata{6,1},'%*[CENTER]%f%f%f');
elements=sscanf(C.textdata{5,1},'%*[NELEMENTS]%i%i%i');
spacing=sscanf(C.textdata{4,1},'%*[SPACING]%f');

%get_map_point=reshape (C.data,elements(1)+1,elements(2)+1,elements(3)+1)

x_surfaces=[centre(1)-spacing*elements(1)/2, centre(1)+spacing*elements(1)/2]
y_surfaces=[centre(2)-spacing*elements(2)/2, centre(2)+spacing*elements(2)/2]
z_surfaces=[centre(3)-spacing*elements(3)/2, centre(3)+spacing*elements(3)/2]

%%%calculate charge
if (strcmp(keep_original_charge,'false'))
    
    if (isempty(P_C) & isempty(P_O) & isempty(P_N))
%{
        write_mol_xyz(mol,'moh','temp.xyz')
        copyfile('d:\MATLAB\ttt\temp.xyz','d:\OpenBabel-2.3.2\temp.xyz')
        %delete d:\MATLAB\ttt\temp.xyz
        cd d:\OpenBabel-2.3.2\
        dos('obabel -i xyz temp.xyz -o mol2 -O temp.mol2 --partialcharge gasteiger')
        %dos('obabel -i xyz temp.xyz -o mol2 -O temp.mol2 --minimize --ff MMFF94 ')
        %dos('obabel -i mol2 temp.mol2 -o mol2 -O temp2.mol2 -h --minimize --ff MMFF94')
       % dos('obabel -i mol2 temp.mol2 -o mol2 -O temp2.mol2 -h') %-partialcharge gasteiger')
        cd d:\MATLAB\ttt
        temp=importdata('d:\OpenBabel-2.3.2\temp.mol2','')
        for i=1:num_atoms
        a= sscanf(temp{i+7},'%i %c %f%f%f %*s %i %*s %f');
        charge=a(end);
        mol(i,15)=charge;
        end
%}
    else



        for i=1:num_atoms
            switch (mol(i,14)) ;case 'C' ;charge=P_C; case 'O'; charge=P_O; case 'N'; charge=P_N;end
            mol(i,15)=charge;
        end
    end

end

for i=1:num_atoms
 mol(i,:)   
X=mol(i,2);Y=mol(i,3);Z=mol(i,4);
switch (mol(i,14)) ; case 67; map=C; case 79; map=O; case 78; map=N; case 72; map=H; end

if (X<=x_surfaces(1)+spacing | X>=x_surfaces(2)-spacing | Y<=y_surfaces(1)+spacing | Y>=y_surfaces(2)-spacing | Z<=z_surfaces(1)+spacing | Z>=z_surfaces(2)-spacing )
    energy=energy+1000;
else

distances(1,:)=abs([(x_surfaces(1)-X)/spacing, round((x_surfaces(1)-X)/spacing)-1, round((x_surfaces(1)-X)/spacing)+1]);
distances(2,:)=abs([(y_surfaces(1)-Y)/spacing, round((y_surfaces(1)-Y)/spacing)-1, round((y_surfaces(1)-Y)/spacing)+1]);
distances(3,:)=abs([(z_surfaces(1)-Z)/spacing, round((z_surfaces(1)-Z)/spacing)-1, round((z_surfaces(1)-Z)/spacing)+1]);

x=distances(1,1)+1;
y=distances(2,1)+1;
z=distances(3,1)+1;
x1=distances(1,2)+1; x2=distances(1,3)+1;
y1=distances(2,2)+1; y2=distances(2,3)+1;
z1=distances(3,2)+1; z2=distances(3,3)+1;
x1
x2
y1
y2
z1
z2
if (x==0 | x1==0 | x2==0 |y==0 | y1==0 | y2==0 |z==0 | z1==0 | z2==0)% energy=energy+100; 
else
    
    %%%%% VDW energy %%%%%%
xd=(x-x1)/(x2-x1); yd=(y-y1)/(y2-y1); zd=(z-z1)/(z2-z1);

e_x1y1z1_x2y1z1= xd*map.data(get_map_point([x2,y1,z1],elements))  +  map.data(get_map_point([x1,y1,z1],elements))*(1-xd)
e_x1y2z1_x2y2z1= xd*map.data(get_map_point([x2,y2,z1],elements))  +  map.data(get_map_point([x1,y2,z1],elements))*(1-xd)
e_x1y1z2_x2y1z2= xd*map.data(get_map_point([x2,y1,z2],elements))  +  map.data(get_map_point([x1,y1,z2],elements))*(1-xd)
e_x1y2z2_x2y2z2= xd*map.data(get_map_point([x2,y2,z2],elements))  +  map.data(get_map_point([x1,y2,z2],elements))*(1-xd)

e_y1z1_y2z1= yd* e_x1y2z1_x2y2z1 + e_x1y1z1_x2y1z1 *(1-yd);
e_y1z2_y2z2= yd* e_x1y2z2_x2y2z2 + e_x1y1z2_x2y1z2 *(1-yd);


e_z1_z2= zd * e_y1z2_y2z2 + e_y1z1_y2z1 *(1-zd);

energy=energy+e_z1_z2

        %%%%% Electrostatic energy %%%%%%
xd=(x-x1)/(x2-x1); yd=(y-y1)/(y2-y1); zd=(z-z1)/(z2-z1);

e_x1y1z1_x2y1z1= xd*e.data(get_map_point([x2,y1,z1],elements))  +  e.data(get_map_point([x1,y1,z1],elements))*(1-xd)
e_x1y2z1_x2y2z1= xd*e.data(get_map_point([x2,y2,z1],elements))  +  e.data(get_map_point([x1,y2,z1],elements))*(1-xd)
e_x1y1z2_x2y1z2= xd*e.data(get_map_point([x2,y1,z2],elements))  +  e.data(get_map_point([x1,y1,z2],elements))*(1-xd)
e_x1y2z2_x2y2z2= xd*e.data(get_map_point([x2,y2,z2],elements))  +  e.data(get_map_point([x1,y2,z2],elements))*(1-xd)

e_y1z1_y2z1= yd* e_x1y2z1_x2y2z1 + e_x1y1z1_x2y1z1 *(1-yd);
e_y1z2_y2z2= yd* e_x1y2z2_x2y2z2 + e_x1y1z2_x2y1z2 *(1-yd);


e_z1_z2= zd * e_y1z2_y2z2 + e_y1z1_y2z1 *(1-zd);

energy=energy+ (e_z1_z2*mol(i,15)); %multiply energy form e.map by atom charge because autodock4 uses +1 probe in e.map calculation


end

end



%energy=get_map_point(grid_x,grid_y,grid_z)
end

end
