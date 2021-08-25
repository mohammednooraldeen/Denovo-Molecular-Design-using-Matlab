function [ output_args ] = hexgon( numx,numy,numz )
%HEXGON Summary of this function goes here
%   Detailed explanation goes here

x = [0:1/numx:1];
y = [0:1/numy:1];
z = [0:1/numz:1];

%Where numx, numy and numz are the number of basic hex elements u want in x, y and z %Direction.

%then follow the following step

[X Y Z] = meshgrid(x,y,z);

X = X(:);Y=Y(:);Z=Z(:);
%and then to plot the mesh you need following lines of code

node = [X Y Z];

figure(1);

if sloped==1
patch('Vertices',node,'Faces',faces,...
'FaceVertexCData',hsv(1),'FaceColor','none')
else
patch('Vertices',node,'Faces',faces,...
'FaceVertexCData',hsv(1),'FaceColor','none')
end
view(3); axis square

title(['Cartesian Mesh ', num2str(numx,3),'x',num2str(numy,3),'x',num2str(numz,3)])

% Now here one thing which you need to Compute is the Face connectivites which shoud %be fed into the function Patch which basically patches the different faces of the %hex and there by makes a complete Hexa Hedra. Now to get the Face connectivties you %need to use the following piece of code.

function faces = face_connectivity(num_u,num_v,num_w)


numx = num_u;
numy = num_v;
numz = num_w;

nnodex = numx+1;
nnodey = numy+1;
nnodez = numz+1;

face_pattern = [1 2 nnodex+2 nnodex+1]; % This is your face connectivity Pattern

nnx = numx+1 ;
nny = (nnodex)*(nnodey) ;
inc_u = 1;
inc_v = nnx;
inc_w = nny;
node_pattern=[ 1 2 nnx+2 nnx+1 nny+1 nny+2 nny+nnx+2 nny+nnx+1 ]; % Node connectivity
element = zeros(numx*numy*numz,8);
element = make_elem_hexa(node_pattern,numx,numy,numz,inc_u,inc_v,inc_w,nnx);

% ThisFunction gives the element connectivity.

faces = zeros(1,4);
face = zeros(6,4);
face1 = [1 2 3 4];
face2 = [4 3 7 8];
face3 = [5 6 7 8];
face4 = [2 6 7 3];
face5 = [1 5 8 4];
face6 = [1 2 6 5];

[m,n] = size(element);

for i = 1:size(element,1)

face = [element(i,face1);element(i,face2);element(i,face3);element(i,face4);element(i,face5);element(i,face6)];
faces = cat(1,faces,face);
end

faces(1,:) = [];
faces = faces;
end
end

