function zerntest2(img, order)
%ZERNTEST2 Compare Zernike reconstructions with source image

[Z, nm] = zern2(img, order);
global recon
recon = zeros([size(img), order]);
recon(:, :, 1) = Z(1)*zernpoly2(nm(1,1), nm(1,2), size(img));
for i = 2:order
    % Initialize values from order n-1
    recon(:, :, i) = recon(:, :, i-1);
    
    % Extract subscript and coefficient values for order n
    subnm = nm(nm(:,1)==i, :);
    subZ = Z(nm(:,1)==i);
    
    % Add contributions from all relevant m-projections
    for j = 1:size(subnm, 1)
        recon(:, :, i) = recon(:, :, i) + subZ(j)*zernpoly2(subnm(j,1), subnm(j,2), size(img));
    end
end    


figure
imagesc(img)
colorbar
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

f = figure;
ax = axes('Parent', f, 'position', [0.17 0.11  0.73 0.83]);

h = imagesc(ax, real(recon(:, :, order)));
colorbar
ylabel(sprintf('Order: %i', order))
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

b = uicontrol('Parent', f, 'Style', 'slider','value', order, 'min', 1, 'max', order, ...
    'Position',[81,15,419,23], 'Callback', @slidercall, 'SliderStep', [1/(order-1), 1]);
% b.Callback = @(es, ed) imagesc(ax, real(recon(:, :, round(es.Value))));
% b.Callback = @(es, ed) print(es.Test);
% set(b, 'Test', 5)
% handles = guidata(b)

end

function slidercall(es, ~)
    global recon
    set(es, 'Value', round(es.Value))
    imagesc(real(recon(:, :, es.Value)));
    ylabel(sprintf('Order: %i', es.Value))
    colorbar
    set(gca, 'XTickLabel', [])
    set(gca, 'YTickLabel', [])
end