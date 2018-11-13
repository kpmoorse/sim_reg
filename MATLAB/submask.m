function sub = submask(img, ctr, varargin)
%SUBMASK Extract a masked subset of an image
%   Detailed explanation goes here

% Loop over varargin elements
while ~isempty(varargin)
    switch lower(varargin{1})

        case 'mask'
            mask = varargin{2};
        case 'radius'
            radius = varargin{2};

    end

    varargin(1:2) = [];
end

if ~exist('radius', 'var')
    radius = 20;
end
r = radius;
if ~exist('mask', 'var')
    [X,Y] = ndgrid(-r:r, -r:r);
    mask = double(sqrt(X.^2+Y.^2)<=r);
end

sub = double(img(round(ctr(1))+(-r:r), round(ctr(2))+(-r:r)));
sub = sub.*mask;

end

