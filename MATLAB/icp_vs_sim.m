function icp_vs_sim(scene, scn_im, model, mdl_im, num_iters, l)

t_icp = icp2(scene, model, num_iters);
t_sim = simreg2(scene, scn_im, model, mdl_im, num_iters);
t_sim2 = simreg2(scene, scn_im, model, mdl_im, num_iters, l);

T_icp = projective2d(t_icp);
T_sim = projective2d(t_sim);
T_sim2 = projective2d(t_sim2);

subplot(2,4,1)
imgformat(scn_im, mdl_im)
title("Original")

subplot(2,4,5)
plot(scene(:,1), scene(:,2), 'g.', model(:,1), model(:,2), 'r.')
pltformat()

subplot(2,4,2)
imgformat(scn_im, imwarp(mdl_im, T_icp))
title("Naive ICP")

subplot(2,4,6)
mdl_icp = append_ones(model)*t_icp;
plot(scene(:,1), scene(:,2), 'g.', mdl_icp(:,1), mdl_icp(:,2), 'r.')
pltformat()

subplot(2,4,3)
imgformat(scn_im, imwarp(mdl_im, T_sim))
title("Hu-Sim ICP (\lambda=0)")

subplot(2,4,7)
mdl_sim = append_ones(model)*t_sim;
plot(scene(:,1), scene(:,2), 'g.', mdl_sim(:,1), mdl_sim(:,2), 'r.')
pltformat()

subplot(2,4,4)
imgformat(scn_im, imwarp(mdl_im, T_sim2))
title(sprintf("Hu-Sim ICP (\\lambda=10^{%i})",round(log(l)/log(10))))

subplot(2,4,8)
mdl_sim2 = append_ones(model)*t_sim2;
plot(scene(:,1), scene(:,2), 'g.', mdl_sim2(:,1), mdl_sim2(:,2), 'r.')
pltformat()

end

% Append a column of ones to the input matrix if it is not already Nx3
function mat = append_ones(mat)

if size(mat, 2) < 3
    mat = [mat, ones(size(mat, 1),1)];
end

end

function imgformat(scn_im, mdl_im)

imshow ...
    ( ...
    imadjust ...
        ( ...
        imfuse ...
            ( ...
            scn_im, ...
            mdl_im, ...
            'falsecolor', ...
            'scaling', 'joint', ...
            'ColorChannels',[1,2,0] ...
            ), ...
        [0 0 0; 0.5 0.5 0.5] ...
        ) ...
    )

end

function pltformat()

axis([1,256,1,256])
set(gca,'Ydir','reverse')
set(gca,'Color','k')
daspect([1 1 1])
l = legend('Scene (Anatomical)', 'Model (Functional)');
l.TextColor = 'w';

end

