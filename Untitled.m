


%import data from .map files
C=importdata('receptor.C.map','',6); %all(:,1)=C.data
N=importdata('receptor.N.map','',6); %all(:,2)=N.data
O=importdata('receptor.O.map','',6); %all(:,3)=O.data
H=importdata('receptor.H.map','',6); %all(:,4)=H.data
P=importdata('receptor.P.map','',6); %all(:,5)=P.data
E=importdata('receptor.e.map','',6); %all(:,5)=P.data

centre=sscanf(C.textdata{6,1},'%*[CENTER]%f%f%f');
elements=sscanf(C.textdata{5,1},'%*[NELEMENTS]%i%i%i');
spacing=sscanf(C.textdata{4,1},'%*[SPACING]%f');

%dddd=index_less_than_zero(data2_box)

C_box=reshape (C.data,elements(1)+1,elements(2)+1,elements(3)+1);
N_box=reshape (N.data,elements(1)+1,elements(2)+1,elements(3)+1);
O_box=reshape (O.data,elements(1)+1,elements(2)+1,elements(3)+1);
H_box=reshape (H.data,elements(1)+1,elements(2)+1,elements(3)+1);
P_box=reshape (P.data,elements(1)+1,elements(2)+1,elements(3)+1);
E_box=reshape (E.data,elements(1)+1,elements(2)+1,elements(3)+1);

index=index_between(C_box,Inf,-Inf)
all=xyz_of(index,elements(1),elements(2),elements(3),spacing,centre(1),centre(2),centre(3));
all(:,5)=N.data(:,1);
all(:,6)=O.data(:,1);
all(:,7)=H.data(:,1);
all(:,8)=P.data(:,1);
all(:,9)=E.data(:,1);

all

%ddddd=xyz_of(dddd(:,1:4),62,56,56,0.375,3.669,12.171,-2.829)
 x= 10; scatter3 (all(1:x:end,1),all(1:x:end,2),all(1:x:end,3),20,all(1:x:end,4)); hold on; scatter3(receptor(:,2),receptor(:,3),receptor(:,4),10,'fill'); hold off