function [atom_color] = atom_type_color(atom_type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
switch atom_type
    case 'O'
        atom_color='red';
    case 'N'
        atom_color='blue';
    case 'C'
        atom_color='black';
    otherwise
        atom_color='black';
end
end

