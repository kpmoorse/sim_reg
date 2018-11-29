function t_est = icp2(scene, model, num_iters, W)
%ICP2 Register model to scene with the Iterative Closest Point algorithm

if ~exist('W', 'var')
    W = eye(size(model,1));
end

scene = append_ones(scene);
model = append_ones(model);

% Initialize transform estimate as identity matrix
t_est = eye(3);

for i = 1:num_iters

    % Apply transform to model and find nearest neighbors
    model_current = model*t_est;
    indices = knnsearch(scene, model_current);

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