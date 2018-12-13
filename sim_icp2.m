function SR = sim_icp2(scene, scn_im, model, mdl_im, varargin)
%SIMREG2 Register model to scene, weighting correspondence by image similarity
%   Detailed explanation goes here

SR = struct('scp', scene, ...
    'sci', scn_im, ...
    'mdp', model, ...
    'mdi', mdl_im);

% Set default varargin values
num_iters = 10;
l = 0;
diagplot = true;
t_init = eye(3);
trim = 0.75;
order = 20;
simfun = @vcorr;

varargin = varargin{1};
% Loop over varargin elements
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'num_iters'
            num_iters = varargin{2};
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
        case 'trim'
            trim = varargin{2};
            assert(and(0<trim, trim<=1), ...
                "Least squares trim value must be on (0,1]")
        case 'order'
            order = varargin{2};
            assert(and(mod(order,1)==0, order>0), ...
                "Zernike order must be a positive integer")
        case 'simfun'
            simfun = str2func(varargin{2});
        otherwise
            error("''%s'' is not a recognized input", varargin{1})
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
disp("Calculating Zernike moments...")
% Precalculate moment vectors for each point
% order = 25; % Maximum order for moment calculation
nm = string(); % List of subscripts
for n=0:order
    for m=0:n
        if mod(n-abs(m), 2)==0
            nm = [nm; [num2str(n),num2str(m)]];
        end
    end
end
nm = nm(2:end);

% Pre-calculate Zernike moments for all scene and model points
scn_mmnt = zeros(size(scene,1),numel(nm));
mdl_mmnt = zeros(size(model,1),numel(nm));
for i=1:size(scene,1)
    scn_mmnt(i,:) = zern2(submask(scn_im, scene(i, 1:2)),order).';
end
for i=1:size(model,1)
    mdl_mmnt(i,:) = zern2(submask(mdl_im, model(i, 1:2)),order).';
end

% Normalize moment magnitudes
% mean_mmnt = mean([scn_mmnt;mdl_mmnt]);
% scn_mmnt = scn_mmnt./mean_mmnt;
% mdl_mmnt = mdl_mmnt./mean_mmnt;

% Initialize transform estimate as identity matrix
t_est = t_init;

for i = 1:num_iters

    % Apply transform to model and find nearest neighbors
    model_current = model*t_est;
    indices = knnsearch(scene, model_current);

    % Generate weights using image similarity via Hu moments
    wvec = zeros(size(model, 1), 1);
    for j=1:size(model, 1)
        k = indices(j);
        a = scn_mmnt(k,:);
        b = mdl_mmnt(j,:);
        wvec(j) = simfun(a, b);
        W(j,j) = simfun(a, b);
    end
    
    % Calculate error distances and trim high-error points
    dvec = vecnorm(model_current - scene(indices, :), 2, 2);
    [~, didx] = sort(dvec);
    didx = didx(1:floor(end*trim));
    
    % Generate trimmed weight matrix
    scn_trim = scene(indices(didx), :);
    mdl_trim = model_current(didx, :);
    W_trim = diag(wvec(didx));
    
    % Iterate transform by linear regression
    t_est_new = t_est*linreg(mdl_trim, scn_trim, W_trim, l);

    % Break if there is no change in t_est
    if ~all(t_est_new(:) == t_est(:))
        t_est = t_est_new;
    else
        model_current = model*t_est;
        break
    end 

end

if ~diagplot
    return
end

% Save initial results to struct
SR.ixi = indices;
corrlist = zeros(numel(indices),1);
for i=1:numel(indices)
    if indices(i) > 0
        corrlist(i) = simfun(mdl_mmnt(i,:), scn_mmnt(indices(i),:));
    else
        corrlist(i) = 0;
    end
end
SR.cli = corrlist;

ixrm = rm_dup(SR.ixi, SR.cli);
SR.ixi(ixrm) = 0;
SR.cli(ixrm) = 0;

% Display diagnostic plots for each neuron
T = projective2d(t_est);
rad = ceil(size(submask(scn_im, scene(1, 1:2)), 1)/2);

disp('Check registration results')
disp('[Y:Accept] / N:Reject / A:Adjust / X:Exit / #:GoTo')

% Loop through registration results
i = 1;
while i<=size(model,1)
    
    j = indices(i);
    
    if SR.ixi(i)==0
        i = i+1;
        continue %skip to next iteration
    end
    if i>1
        clf %clear figure
    end
    
    % Plot overlaid scene and transformed model images
    subplot(1,2,1)
    imshow(imadjust(imfuse(imwarp(scale(mdl_im), T), scale(scn_im)), [0 0 0; 0.5 0.5 0.5])) % mdl=red, scn=grn
    hold on
    plot([scene(j,1), model_current(i,1)], [scene(j,2), model_current(i,2)], 'o-w')
    if l==0
        title("Similarity-Weighted ICP (\lambda=0)")
    else
        title(sprintf("Similarity-Weighted ICP (\\lambda=10^{%.1f})",round(log(l)/log(10))))
    end
    
    % Plot model submask
    subplot(2,4,3)
    imshow(submask(mdl_im, model(i, 1:2)),[0,500])
    hold on
    plot(rad, rad, 'm.', 'MarkerSize', 15)
    title(sprintf('Model[%i] (Functional)', i))
    
    % Plot scene submask
    subplot(2,4,4)
    imshow(submask(scn_im, scene(j, 1:2)),[0,500])
    hold on
    plot(rad, rad, 'g.', 'MarkerSize', 15)
    title(sprintf('Scene[%i] (Anatomical)', j))
    
    % Plot comparative moments
    subplot(2,2,4)
    a = mdl_mmnt(i, :);
    a = reshape([real(a);imag(a)],1,[]);
    b = scn_mmnt(j, :);
    b = reshape([real(b);imag(b)],1,[]);
    
    scatter(a,b,400,(1:numel(a))/2,'.')
    h = colorbar;
    ylabel(h, 'Moment index');
    
    xrng = prctile(a, [5, 95]);
    yrng = prctile(b, [5, 95]);
    xlim(xrng + [-1,1]*0.05*diff(xrng));
    ylim(yrng + [-1,1]*0.05*diff(yrng));
    xlabel('Model (|A_{mn}|)')
    ylabel('Scene (|A_{mn}|)')
    title(sprintf("Normalized Zernike Moments (sim = %.04f)", simfun(a, b)))
    
    flag = true;
    while flag
        % Prompt user evaluation of results
        uin = input(sprintf('(#%03i): ', i), 's');
        switch true
            % On "Accept," move to next point
            case isempty(uin) || strcmpi(uin, 'y')
                flag = false;
                i = i + 1;
            % On "Reject," delete correspondance and move to next point 
            case strcmpi(uin, 'n')
                indices(i) = 0;
                flag = false;
                i = i + 1;
            % On "Adjust," launch adjustment GUI
            case strcmpi(uin, 'a')
                % Display nearby scene points
                subplot(2,4,4)
                pts = scene(vecnorm(scene - scene(j,:),2, 2)<rad/2, 1:2);
                pts = pts - scene(j,1:2) + [rad rad];
                plot(pts(:,1), pts(:,2), 'go', 'MarkerSize', 10)
                
                % Prompt user selection of scene point
                gin_accept = false;
                while ~gin_accept
                    gin = myginput(1, 'arrow');
                    if isempty(gin)
                        gin = [rad, rad];
                    end
                    [~,adj_ix] = min(vecnorm(pts-gin, 2, 2));
                    adj_pt = pts(adj_ix, :);
                    % Only accept nearby selection clicks
                    if norm(gin-adj_pt) <= 2
                        gin_accept = true;
                    end
                end
                adj_pt = adj_pt + scene(j,1:2) - [rad rad];
                scene_ix = find(and(scene(:,1)==adj_pt(1), scene(:,2)==adj_pt(2)));
                
                % Update correspondance
                indices(i) = scene_ix;
                j = indices(i);
                
                flag = false;

            case strcmp(uin, 'x')
                flag = false;
                i = size(model,1) + 1;
            case ~isempty(str2double(uin)) && ~isnan(str2double(uin))
                if any(str2double(uin)==1:size(model,1))
                    flag = false;
                    i = str2double(uin);
                else
                    warning('Numerical input is outside model index range (1:%i)', size(model,1))
                end
            otherwise
                warning('Input not recognized')

        end
    end
end

% Save final results
SR.te = t_est;
SR.ix = indices;
corrlist = zeros(numel(indices),1);
for i=1:numel(indices)
    if indices(i) > 0
        corrlist(i) = simfun(mdl_mmnt(i,:), scn_mmnt(indices(i),:));
    else
        corrlist(i) = 0;
    end
end
SR.cl = corrlist;

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
% function r = corr(a, b)
%     r = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
% end

function mat = scale(mat)
    mat = mat/max(mat(:));
end
function mat = normalize(mat)
    mat = mat/sum(mat(:));
end