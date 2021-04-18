function [padded_imgs,valid_rows,valid_cols] = symextend_3d(imgs,d)
% Extends the image volume height and width dimensions such that its
% spatial resolution is divisible by "d".
% 
% IMGS: the input image volume
% D: is an even integer number
%
% USAGE:
%   [padded_imgs,valid_rows,valid_cols] = symextend_3d(imgs,32);
%   frame1 = padded_imgs(valid_rows,valid_cols,1)
% 

% Make the image size even
imgs = make_size_even(imgs);

[R,C,DD] = size(imgs);

% make the size divisible by "d"
row_padsize = (d - mod(R,d))/2;
col_padsize = (d - mod(C,d))/2;

padded_imgs = zeros(R + 2*row_padsize, C + 2*col_padsize, DD);

for frame = 1:size(imgs,3)
    tmp = wextend('addrow','sym',imgs(:,:,frame),row_padsize);
    padded_imgs(:,:,frame) = wextend('addcol','sym',tmp,col_padsize);
end

valid_rows = row_padsize + 1:row_padsize + R;
valid_cols = col_padsize + 1:col_padsize + C;

