function I = minfo(u, v)
%MINFO Calculate the mutual information between two vectors u and v
%
%   Approximate mutual information is calculated according to methods
%   defined in:
%   Kraskov, A., Stogbauer, H., and Grassberger, P., "Estimating Mutual
%       Information" (2008)
%
%   Bins are alloted such that an average of 20 points fall into each px
%   and py bin (less in each pxy bin)


assert(numel(u)==numel(v), ...
    'u and v must be vectors of equal length');

% Calculate either singular or summed MI depending on whether u and v have
% a complex component
if isreal(u) && isreal(v)
    I = mi_comp(u, v);
else
    % Complex summed MI assumes Re(u) and Im(u) are independent
    I = mi_comp(real(u), real(v)) + mi_comp(imag(u), imag(v));
end

end

% Calculate mutual information for real vectors u and v
function Ii = mi_comp(u, v)

nbins = floor(numel(u)/20);
hgram = histcounts2(real(u), real(v), nbins);

pxy = hgram / sum(hgram, 'all');
px = sum(pxy, 2);
py = sum(pxy, 1);
px_py = px * py;
nzs = pxy ~= 0;

Ii = sum(pxy(nzs) .* log(pxy(nzs) ./ px_py(nzs)));

end


