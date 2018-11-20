function [Z, nm] = zern2(img, order, varargin)
%ZERN2 Calculate rotation-invariant Zernike moments for an input image
%
%   Implementation of methods from the following paper:
%   A. Khotanzad and Y. H. Hong, "Invariant Image Recognition by Zernike
%       Moments" (1990)


% Set default varargin values
outform = "complex";

% Loop over varargin elements
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'outform'
            outform = varargin{2};
        otherwise
            warning("%s in not a recognized argument", varargin{1})

    end

    varargin(1:2) = [];
end

% Loop over allowed indices
Z = [];
nm = [];
for n=0:order
    for m=0:n
        if mod(n-abs(m), 2)==0
            % Append Zernike moments and n,m values
            Anm = A(n,m,img);
            Z = [Z; Anm];
            nm = [nm; [n,m]];
        end
    end
end

switch outform
    case 'mag-phase'
        Z = [abs(Z), angle(Z)];
end
end

% Calculate Zernike moment Anm (dependent on Vnm)
function Anm = A(n,m,img)

y = ctrrng(size(img, 1));
x = ctrrng(size(img, 2));
[X,Y] = ndgrid(x,y);
r = sqrt(X.^2+Y.^2);
th = mod(atan2(Y,X),2*pi);

sumterms = double(img).*conj(V(n,m,r,th));
Anm = (n+1)/pi * sum(sumterms(:));

end

% Define Zernike polynomial Vnm (dependent on Rnm)
function Vnm = V(n,m,r,th)

Rnm = R(n,m,r);
Vnm = Rnm .* exp(1j * m * th);

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