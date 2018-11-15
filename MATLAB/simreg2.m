function [t_est, indices] = simreg2(scene, scn_im, model, mdl_im, num_iters, varargin)
%SIMREG2 Register model to scene, weighting correspondence by image similarity
%   Detailed explanation goes here

% Set default varargin values
l = 0;
diagplot = true;
t_init = eye(3);

% Loop over varargin elements
while ~isempty(varargin)
    switch lower(varargin{1})

        case 'lambda'
            l = varargin{2};
            assert(and(isnumeric(l), numel(l)==1, mod(l,1)==0), ...
                "Lambda must be a scalar integer")
        case 'diagnostic'
            diagplot = varargin{2};
        case 'initial'
            t_init = varargin{2};
            assert(and(isnumeric(t_init), all(size(t_init)==[3,3])), ...
                "Initial guess must be a 3x3 numeric array")

    end

    varargin(1:2) = [];
end

% Ensure proper point cloud formatting
scene = append_ones(scene);
model = append_ones(model);

% Initialize weight matrix
W = eye(size(model,1));

if ~exist('l', 'var')
    l = 0;
end

% Precalculate Hu moment vectors for each point
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
t_est = t_init;

for i = 1:num_iters

    % Apply transform to model and find nearest neighbors
    model_current = model*t_est;
    indices = knnsearch(scene, model_current);

    % Generate weights using image similarity via Hu moments
    for j=1:size(model, 1)
        k = indices(j);
        a = scn_hu(k,:);
        b = mdl_hu(j,:);
        W(j,j) = corr(a, b);
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

if ~diagplot
    return
end

% Display diagnostic plots for each neuron
T = projective2d(t_est);
rad = ceil(size(submask(scn_im, scene(1, 1:2)), 1)/2);

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
    plot(rad, rad, 'g.', 'MarkerSize', 15)
    title(sprintf('Scene[%i] (Anatomical)', i))
    
    subplot(2,4,4)
    imshow(submask(mdl_im, model(i, 1:2)),[0,500])
    hold on
    plot(rad, rad, 'r.', 'MarkerSize', 15)
    title(sprintf('Model[%i] (Functional)', i))
    
    subplot(2,2,4)
    a = scn_hu(j, :);
    b = mdl_hu(i, :);
    bar([a; b].')
    legend('Scene', 'Model')
    title(sprintf("Hu Moments (corr = %.04f)", corr(a, b)))
    
    flag = true;
    while flag
        uin = input(sprintf('[#%03i] [Y:Accept] / N:Reject / A:Adjust: ', i), 's');
        switch lower(uin)
            case 'y'
                flag = false;
            case 'n'
                indices(i) = 0;
                flag = false;
            case 'a'
                subplot(2,4,3)
                pts = scene(vecnorm(scene - scene(j,:),2, 2)<rad/2, 1:2);
                pts = pts - scene(j,1:2) + [rad rad];
                plot(pts(:,1), pts(:,2), 'go', 'MarkerSize', 10)
                
                gin = myginput(1, 'arrow');
                [~,adj_ix] = min(vecnorm(pts-gin, 2, 2));
                adj_pt = pts(adj_ix, :);
                adj_pt = adj_pt + scene(j,1:2) - [rad rad];
                scene_ix = find(and(scene(:,1)==adj_pt(1), scene(:,2)==adj_pt(2)));
                
                indices(i) = scene_ix;
                j = indices(i);
                
                subplot(2,4,3)
                imshow(submask(scn_im, scene(j, 1:2)),[0,500])
                hold on
                plot(rad, rad, 'g.', 'MarkerSize', 15)
                title(sprintf('Scene[%i] (Anatomical)', i))
            otherwise
                flag = false;
        end
    end
end

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

% function s = sim(a, b)
%     s = dot(a,b)/(norm(a)*norm(b));  
% end

function mat = scale(mat)
    mat = mat/max(mat,[],'all');
end