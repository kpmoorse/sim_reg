function t_est = icp3(scene, model, num_iters)

    % Initialize transform estimate as identity matrix
    t_est = eye(4);
    
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
    plot3(scene(:,1), scene(:,2), scene(:,3), '.')
    hold on
    plot3(model(:,1), model(:,2), model(:,3), '.')
    plot3(model_current(:,1), model_current(:,2), model_current(:,3), 'x')
    for i = 1:numel(indices)
        plot3([scene(indices(i),1), model(i,1)],...
              [scene(indices(i),2), model(i,2)],...
              [scene(indices(i),3), model(i,3)], 'color', [0 0 0]+0.5)
        plot3([scene(indices(i),1), model_current(i,1)],...
              [scene(indices(i),2), model_current(i,2)],...
              [scene(indices(i),3), model_current(i,3)], 'color', [0 0 0]+0.5)
    end
%     legend('Scene', 'Model', 'Transformed Model')
    daspect([1,1,1])
    hold off
    
end

function B = linreg(X, Y)

    B = inv(X.'*X)*X.'*Y;
    
end