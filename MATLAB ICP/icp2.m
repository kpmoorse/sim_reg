function t_est = icp2(scene, model, num_iters)

    % Initialize transform estimate as identity matrix
    t_est = eye(3);
    
    for i = 1:num_iters
    
        % Apply transform to model and find nearest neighbors
        model_current = model*t_est;
        indices = knnsearch(scene, model_current);
        
        % Iterate transform by linear regression
        t_est_new = t_est*linreg(model_current, scene(indices, :));
        
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
    end
    plot(model_current(:,1), model_current(:,2), 'x')
    hold off
    
end

function B = linreg(X, Y)

    B = inv(X.'*X)*X.'*Y;
    
end