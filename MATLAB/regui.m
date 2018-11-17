function indices = regui(num_iters)
%REGUI Select files via UI and call registration function
%   Detailed explanation goes here

if ~exist('num_iters', 'var')
    num_iters = 10;
end

z = 3;

% Select data files via user input
lmat = uigetfile({'*.tif';'*.tiff'}, "Select label matrix (anatomical)");
img = uigetfile({'*.tif';'*.tiff'}, "Select image stack (anatomical)");
cnmf = uigetfile({'*.mat'}, "Select CNMF data (functional)");

% Read data files
lmat = bigread2(lmat);
img = bigread2(img);
load(cnmf, 'CNM');

% Reformat data for function input
scene = exctrs2(lmat(:,:,z));
scene = scene(:, 2:3);
scn_im = img(:,:,z);

model = fliplr(CNM.cm);
model = model(:, 1:2);
mdl_im = mean(CNM.Y, 3);

% Run similarity registration
[~, indices] = simreg2(scene, scn_im, model, mdl_im, num_iters, 'moment', 'zern2');

end

