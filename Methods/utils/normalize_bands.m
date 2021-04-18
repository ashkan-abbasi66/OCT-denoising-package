function normalized_imgs = normalize_bands(imgs)
% imgs: Rows x Columns x bands
[~,~,band] = size(imgs);

normalized_imgs = imgs;

for i =1: band
    y = imgs(:, :, i) ;
    max_y = max(y(:));
    min_y = min(y(:));
    y =  (y - min_y)./ (max_y - min_y);
    normalized_imgs(:, :, i) = y;
end
