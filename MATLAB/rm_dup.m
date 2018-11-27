function [vec, ix] = rm_dup(vec, fitness)
%RM_DUP Remove low-fitness duplicate values from vec

ix = 1:numel(vec);

for val=reshape(unique(vec), 1, [])
    maxfit = max(fitness(vec==val));
    ix(vec==val & fitness~=maxfit) = 0;
end
ix = ix~=0;
% ix(ix==0) = [];

vec = vec(ix);

end