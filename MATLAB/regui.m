function indices = regui(num_iters)
%REGUI Select files via UI and call registration function
%   Detailed explanation goes here

if ~exist('num_iters', 'var')
    num_iters = 10;
end

z = 3;

lmat = uigetfile({'*.tif';'*.tiff'}, "Select label matrix (anatomical)");
img = uigetfile({'*.tif';'*.tiff'}, "Select image stack (anatomical)");
cnmf = uigetfile({'*.mat'}, "Select CNMF data (functional)");

lmat = bigread2(lmat);
lmat = lmat(:,:,z);
img = bigread2(img);
img = img(:,:,z);
load(cnmf);

scene = exctrs2(lmat);
scene = scene(:, 2:3);
scn_im = img;

model = fliplr(CNM.cm);
model = model(:, 1:2);
mdl_im = mean(CNM.Y, 3);

[~, indices] = simreg2(scene, scn_im, model, mdl_im, num_iters);

end

