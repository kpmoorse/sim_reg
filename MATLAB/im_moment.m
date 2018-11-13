function mij = im_moment(img, i, j)
%IMMNT Calculate the moment M_ij of the input image
%   Detailed explanation goes here

% Generate centered coordinates
shape = size(img);
[X, Y] = ndgrid(1:shape(1), 1:shape(2));
x0 = (1+shape(1))/2;
y0 = (1+shape(2))/2;
X = X - x0;
Y = Y - y0;

%Calculate central moment
mij = sum((X.^i).*(Y.^j).*img, 'all');

end

