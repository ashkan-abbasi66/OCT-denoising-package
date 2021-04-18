function imc = crop_image(im, d)
% 
% Crop the input image so as its dimensions are divisible by "d"
% 
% E,g., [450,900] --- "d = 4" ---> [448,900]
% 
% 

nd = ndims(im);

if nd>2
    [H,W,~] = size(im);
else
    [H,W] = size(im);
end

new_size = [H - mod(H,d), W - mod(W,d)];

bbox = [(H - new_size(1))/2 + 1,...
        (W - new_size(2))/2 + 1,...
        (H + new_size(1))/2,...
        (W + new_size(2))/2];
bbox = floor(bbox);

if nd>2
    imc = im(bbox(1):bbox(3), ...
         bbox(2):bbox(4),:);
else
    imc = im(bbox(1):bbox(3), ...
         bbox(2):bbox(4));
end

% --------------------------------------------------------------
% The original code was found in Deep Image Prior (DIP) Package
% 
% % % def crop_image(img, d=32):
% % %     '''Make dimensions divisible by `d`'''
% % % 
% % %     new_size = (img.size[0] - img.size[0] % d, 
% % %                 img.size[1] - img.size[1] % d)
% % % 
% % %     bbox = [
% % %             int((img.size[0] - new_size[0])/2), 
% % %             int((img.size[1] - new_size[1])/2),
% % %             int((img.size[0] + new_size[0])/2),
% % %             int((img.size[1] + new_size[1])/2),
% % %     ]
% % % 
% % %     img_cropped = img.crop(bbox)
% % %     return img_cropped