function r = vcorr(x, y)
% VCORR Calculate the correlation value between complex vectors x and y
%   NOTE #1: abs() forces real-valued output and r=1 for x=y BUT it also
%       forces r>=0. In other words, vcorr assumes positive correlation.
%   NOTE #2: MATLAB dot() automatically takes the conjugate of one
%       argument, per the definition of a complex inner product

% r = abs(dot(x', y) / sqrt(dot(x, x') * dot(y, y')));
r = abs(dot(x, y) / (norm(x)*norm(y)));

end