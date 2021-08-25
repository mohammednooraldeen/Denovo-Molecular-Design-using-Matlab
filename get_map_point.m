function value= get_map_point( xyz,elements )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



x=xyz(1); y=xyz(2); z=xyz(3);

value=( x + (y-1)*(elements(1)+1) + (z-1)*(   (elements(1)+1)  *(elements(2)+1)  ));
end

