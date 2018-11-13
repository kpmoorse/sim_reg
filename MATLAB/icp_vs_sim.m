function icp_vs_sim(scene, scn_im, model, mdl_im, num_iters)

t_icp = icp2(scene, model, num_iters);
t_sim = simreg2(scene, scn_im, model, mdl_im, num_iters);

T_icp = projective2d(t_icp);
T_sim = projective2d(t_sim);

subplot(2,3,1)
imshow(imfuse(scn_im, mdl_im, 'falsecolor','scaling','joint','ColorChannels',[1,2,0]),[1,500])
title("Original")

subplot(2,3,4)
plot(scene(:,1), scene(:,2), '.', model(:,1), model(:,2), '.')
axis([1,256,1,256])
set(gca,'Ydir','reverse')

subplot(2,3,2)
imshow(imfuse(scn_im, imwarp(mdl_im, T_icp), 'falsecolor','scaling','joint','ColorChannels',[1,2,0]),[1,500])
title("Naive ICP")

subplot(2,3,5)
mdl_icp = append_ones(model)*t_icp;
plot(scene(:,1), scene(:,2), '.', mdl_icp(:,1), mdl_icp(:,2), '.')
axis([1,256,1,256])
set(gca,'Ydir','reverse')

subplot(2,3,3)
imshow(imfuse(scn_im, imwarp(mdl_im, T_sim), 'falsecolor','scaling','joint','ColorChannels',[1,2,0]),[1,500])
title("Hu Similarity ICP")

subplot(2,3,6)
mdl_sim = append_ones(model)*t_sim;
plot(scene(:,1), scene(:,2), '.', mdl_sim(:,1), mdl_sim(:,2), '.')
axis([1,256,1,256])
set(gca,'Ydir','reverse')

end

% Append a column of ones to the input matrix if it is not already Nx3
function mat = append_ones(mat)

if size(mat, 2) < 3
    mat = [mat, ones(size(mat, 1),1)];
end

end

