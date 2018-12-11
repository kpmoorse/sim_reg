function ix = rm_dup(vec, fitness, SR)
%RM_DUP Remove low-fitness duplicate values from vec

ix = 1:numel(vec);

% Loop over unique values in vec
for val=reshape(unique(vec), 1, [])
    if val==0
        continue
    end
    plist = find(vec==val);
    npts = numel(plist);
    % If a value is not unique, remove duplicates
    if npts > 2
        % Plot visualization if data structure is provided
        if exist('SR', 'var')
            for i=1:npts
                subplot(1,npts+1,i)
                imagesc(submask(SR.mdi, SR.mdp(plist(i),:)))
                hold on, plot(21,21,'wo')
                daspect([1 1 1])
                title(sprintf('iModel=%i, sim=%.3f', plist(i), SR.cl(plist(i))))
                set(gca, 'XTick', [])
                set(gca, 'YTick', [])
            end
            subplot(1,npts+1,npts+1)
            imagesc(submask(SR.sci, SR.scp(val,:)))
            hold on, plot(21,21,'wo')
            daspect([1 1 1])
            title(sprintf('jScene=%i', val))
            set(gca, 'XTick', [])
            set(gca, 'YTick', [])
            pause
        end
        
        % Remove elements whose fitness is not maximal
        maxfit = max(fitness(vec==val));
        ix(vec==val & fitness~=maxfit) = 0;
    end
end
ix = ix==0;
% vec = vec(ix);

end