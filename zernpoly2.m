function Vnm = zernpoly2(n,m,dim)
%ZERNIKE2 Calculate Zernike polynomials of a given input size
%
%   Implementation of methods from the following paper:
%   A. Khotanzad and Y. H. Hong, "Invariant Image Recognition by Zernike
%       Moments" (1990)

x = ctrrng(dim(1));
y = ctrrng(dim(2));
[X,Y] = ndgrid(x,y);
r = sqrt(X.^2+Y.^2);
th = mod(atan2(Y,X),2*pi);

Rnm = R(n,m,r);
Vnm =  Rnm .* exp(1j * m * th);

end

% Define radial polynomial Rnm
function Rnm = R(n,m,r)

ss = 0:((n-abs(m))/2);
Rnm = 0*r;
for s = ss
    Rnm = Rnm + (-1).^s .* ...
            ff(n-s) / (ff(s) .* ff((n+abs(m))/2 - s) .* ff((n-abs(m))/2 - s)) .* ...
            r.^(n - 2*s);
end
Rnm(r>1) = 0;

end

% Normalize image coordinates to the unit circle
function y = ctrrng(x)
y = 1:x;
y = y - mean(y);
y = y / max(y);
end

% Shorthand factorial function for readability
function fx = ff(x)
fx = factorial(x);
end