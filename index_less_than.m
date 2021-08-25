function [index ] = index_less_than( matrix,value )
%INDEX_LESS_THAN_ZERO Summary of this function goes here
%   Detailed explanation goes here
% colormap('jet') small value= blue, large value =  red
[i]=find(reshape(matrix,1,numel(matrix))<value);
[index(:,1),index(:,2),index(:,3)]=ind2sub(size(matrix),i);
[t]=size(index)
for i=1:t(1,1) 
    index(i,4)=matrix(index(i,1),index(i,2),index(i,3));
end

end

