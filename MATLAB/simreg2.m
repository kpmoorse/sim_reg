function t_est = simreg2(scene, scn_im, model, mdl_im, num_iters)
%SIMREG2 Register model to scene, weighting correspondence by image similarity
%   Detailed explanation goes here

% Ensure proper point cloud formatting
scene = append_ones(scene);
model = append_ones(model);

% Initialize weight matrix
W = eye(size(model,1));

% Precalculate Hu moments for each point
scn_hu = zeros(size(scene,1),7);
mdl_hu = zeros(size(model,1),7);
for i=1:size(scene,1)
    scn_hu(i,:) = hu(submask(scn_im, scene(i, 1:2)));
end
for i=1:size(model,1)
    mdl_hu(i,:) = hu(submask(mdl_im, model(i, 1:2)));
end

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
        W(j,j) = corr(a, b);
    end
    
    % Iterate transform by linear regression
    t_est_new = t_est*linreg(model_current, scene(indices, :), W);

    % Break if there is no change in t_est
    if ~all(t_est_new == t_est, 'all')
        t_est = t_est_new;
    else
        model_current = model*t_est;
        break
    end 

end

% Plot results
plot(scene(:,1), scene(:,2), 'o')
hold on
plot(model(:,1), model(:,2), 'o')
for i = 1:numel(indices)
    plot([scene(indices(i),1), model(i,1)], [scene(indices(i),2), model(i,2)], 'k-')
    plot([scene(indices(i),1), model_current(i,1)], [scene(indices(i),2), model_current(i,2)], 'k-')
end
plot(model_current(:,1), model_current(:,2), 'x')
hold off

end

% Custom linear regression allows non-identity weight matrix
function B = linreg(X, Y, W)

if ~exist('W', 'var')
    W = eye(size(X,1));
end
B = inv(X.'*W*X)*X.'*W*Y;

end

% Append a column of ones to the input matrix if it is not already Nx3
function mat = append_ones(mat)

if size(mat, 2) < 3
    mat = [mat, ones(size(mat, 1),1)];
end

end

% Calculate a basic correlation coefficient between arrays
function r = corr(a, b)

    r = sum(a.*b)/(sum(a)*sum(b));
%     r = mean(abs(a-b)/(a+b));

end