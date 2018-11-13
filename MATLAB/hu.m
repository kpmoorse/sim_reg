function I = hu(img, varargin)
%EXTHU Calculate Hu moments for an image
%   Detailed explanation goes here

% Calculate central moment
n00 = im_moment(img,0,0);

% Calculate normalized moments
n = @(i,j) im_moment(img,i,j) / n00^(1+(i+j)/2);
n01 = n(0,1);
n11 = n(1,1);
n20 = n(2,0);
n02 = n(0,2);
n12 = n(1,2);
n21 = n(2,1);
n30 = n(3,0);
n03 = n(0,3);

% Calculate Hu moments
I = zeros(4, 1);
I(1) = n20 + n02;
I(2) = (n20 + n02)^2 + 4*n11^2;
I(3) = (n30 - 3*n12)^2 + (3*n21 - n03)^2;
I(4) = (n30 + n12)^2 + (n21 + n03)^2;
I(5) = (n30 - 3*n12)*(n30 + n12)*((n30 + n12)^2 - 3*(n21 + n01)^2) + ...
    (3*n21 - n03)*(n21 + n03)*(3*(n30 + n12)^2 - (n21 + n03)^2);
I(6) = (n20 - n02)*((n30 + n12)^2 - (n21 + n03)^2) + 4*n11*(n30 + n12)*(n21 + n03);
I(7) = (3*n21 - n03)*(n30 + n12)*((n30 + n12)^2 - 3*(n21 + n03)^2) - ...
    (n30 - 3*n12)*(n21 + n03)*(3*(n30 + n12)^2 - (n21 + n03)^2);


end

