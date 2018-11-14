function t_est = simreg2(scene, scn_im, model, mdl_im, num_iters, l)
%SIMREG2 Register model to scene, weighting correspondence by image similarity
%   Detailed explanation goes here

% Ensure proper point cloud formatting
scene = append_ones(scene);
model = append_ones(model);

% Initialize weight matrix
W = eye(size(model,1));

if ~exist('l', 'var')
    l = 0;
end

% Precalculate Hu moments for each point
scn_hu = zeros(size(scene,1),7);
mdl_hu = zeros(size(model,1),7);
for i=1:size(scene,1)
    scn_hu(i,:) = hu(submask(scn_im, scene(i, 1:2))).';
end
for i=1:size(model,1)
    mdl_hu(i,:) = hu(submask(mdl_im, model(i, 1:2))).';
end

% Normalize Hu moment magnitudes
mean_hu = mean([scn_hu;mdl_hu]);
scn_hu = scn_hu./mean_hu;
mdl_hu = mdl_hu./mean_hu;

% Initialize transform estimate as identity matrix
t_est = eye(3);

for i = 1:num_iters

    % Apply transform to model and find nearest neighbors
    model_current = model*t_est;
    indices = knnsearch(scene, model_current);

    % Generate weights using image similarity via Hu moments
    for j=1:size(model, 1)
        k = indices(j);
        a = scn_hu(k,:);
        b = mdl_hu(j,:);
        W(j,j) = sim(a, b);
    end
    
    % Iterate transform by linear regression
    t_est_new = t_est*linreg(model_current, scene(indices, :), W, l);

    % Break if there is no change in t_est
    if ~all(t_est_new == t_est, 'all')
        t_est = t_est_new;
    else
        model_current = model*t_est;
        break
    end 

end

% Display diagnostic plots for each neuron
T = projective2d(t_est);
for i=1:numel(model)
    
    j = indices(i);
    
    subplot(1,2,1)
    imshow(imadjust(imfuse(imwarp(scale(mdl_im), T), scale(scn_im), 'falsecolor', 'scaling', 'joint', ...
        'ColorChannels',[1,2,0] ), [0 0 0; 0.5 0.5 0.5])) % mdl=red, scn=grn
    hold on
    plot([scene(j,1), model(i,1)], [scene(j,2), model(i,2)], '.-w')
    if l==0
        title("Hu-Sim ICP (\lambda=0)")
    else
        title(sprintf("Hu-Sim ICP (\\lambda=10^{%.1f})",round(log(l)/log(10))))
    end
    
    subplot(2,4,3)
    imshow(submask(scn_im, scene(j, 1:2)),[0,500])
    hold on
    plot(21, 21, 'g.', 'MarkerSize', 15)
    title(sprintf('Scene[%i] (Anatomical)', i))
    
    subplot(2,4,4)
    imshow(submask(mdl_im, model(i, 1:2)),[0,500])
    hold on
    plot(21, 21, 'r.', 'MarkerSize', 15)
    title(sprintf('Model[%i] (Functional)', i))
    
    subplot(2,2,4)
    a = scn_hu(j, :);
    b = mdl_hu(i, :);
    bar([a; b].')
    legend('Scene', 'Model')
    title(sprintf("Hu Moments (sim = %.04f)", sim(a, b)))
    
    pause
end


% % Plot results
% plot(scene(:,1), scene(:,2), 'o')
% hold on
% plot(model(:,1), model(:,2), 'o')
% for i = 1:numel(indices)
%     plot([scene(indices(i),1), model(i,1)], [scene(indices(i),2), model(i,2)], 'k-')
%     plot([scene(indices(i),1), model_current(i,1)], [scene(indices(i),2), model_current(i,2)], 'k-')
% end
% plot(model_current(:,1), model_current(:,2), 'x')
% hold off

end

% Custom linear regression with weight matrix and lambda constraint
function B = linreg(X, Y, W, l)
if ~exist('W', 'var')
    W = eye(size(X,1));
end
if ~exist('l', 'var')
    l = 0;
end
B = inv(X.'*W*X + l*eye(3))*X.'*W*Y;
end

% Append a column of ones to the input matrix if it is not already Nx3
function mat = append_ones(mat)
if size(mat, 2) < 3
    mat = [mat, ones(size(mat, 1),1)];
end
end

% Calculate a basic correlation coefficient between arrays
function r = corr(a, b)
    r = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
%     r = mean(abs(a-b)/(a+b));
end

function s = sim(a, b)
    s = dot(a,b)/(norm(a)*norm(b));  
end
function mat = scale(mat)
    mat = mat/max(mat,[],'all');
end