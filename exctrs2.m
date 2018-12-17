function [centroids, labels] = exctrs2(img, varargin)
% EXCTRS2 Extract centroids from 2D label matrix
%   For every unique value in the 2D matrix img, find all pixels of that
%   value and calculate their average position (centroid).
%
%   C = EXCTRS2(img) is an N-by-3 matrix where N is the number of unique
%   nonzero values in img and each row contains [value mean(x) mean(y)]
%
%   EXCTRS2(img, 'xdata', x, 'ydata', y) calculates centroids using
%   provided coordinates instead of pixel indices

% Loop over varargin elements
while ~isempty(varargin)
    switch lower(varargin{1})

        case 'xdata'
            x = varargin{2};
        case 'ydata'
            y = varargin{2};

    end

    varargin(1:2) = [];
end

assert(exist('x', 'var') == exist('y', 'var'),...
    'X and Y coordinates must be both defined or both undefined')

if exist('x', 'var')
    assert(and(numel(x)==size(img, 1), numel(y)==size(img, 2)),...
        'XData and YData lengths must match img size')
% Define x and y based on pixel index if not provided as input
else
    shape = size(img);
    x = 1:shape(1);
    y = 1:shape(2);
end

% Convert x and y arrays to ndgrid matrices
[y, x] = ndgrid(x, y);

% Extract label values
labels = unique(img);

% Loop over label values and extract centroids
centroids = zeros(numel(labels), 3);
for i = 1:numel(labels)
    label = labels(i);
    bin = img==label;
    centroid = mean([x(bin), y(bin)], 1);
    centroids(i,:) = [label, centroid];
end

end
