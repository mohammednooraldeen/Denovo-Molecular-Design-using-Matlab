function [ output_args ] = write_xyz_file( mol,mol_name, filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
s_mol=size(mol,1);

fid = fopen(filename, 'w');
fprintf(fid,'%i\n',s_mol);
fprintf(fid,'%s\n',mol_name);


for i=1:s_mol
    fprintf(fid, '%c \t %6.4f \t %6.4f \t %6.4f \n', mol(i,14), mol(i,2:4));
end

fclose(fid);

end

