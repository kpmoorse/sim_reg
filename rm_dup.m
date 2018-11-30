function [vec, ix] = rm_dup(vec, fitness, SR)
%RM_DUP Remove low-fitness duplicate values from vec

ix = 1:numel(vec);

for val=reshape(unique(vec), 1, [])
    plist = find(vec==val);
    npts = numel(plist);
    if npts > 1
        if exist('SR', 'val')
            for i=1:npts
                subplot(1,npts+1,i)
                imagesc(submask(SR.mdi, SR.mdp(plist(i),:)))
                hold on, plot(21,21,'k.')
                daspect([1 1 1])
                title(sprintf('iModel=%i, sim=%.3f', plist(i), SR.cl(plist(i))))
            end
            subplot(1,npts+1,npts+1)
            imagesc(submask(SR.sci, SR.scp(val,:)))
            hold on, plot(21,21,'k.')
            daspect([1 1 1])
            title(sprintf('jScene=%i', val))
            pause
        end
        
        maxfit = max(fitness(vec==val));
        ix(vec==val & fitness~=maxfit) = 0;
    end
end
ix = ix~=0;
% ix(ix==0) = [];
vec = vec(ix);

end