function imgs = make_depth_power_2(imgs)
% 
% Make the 3rd dimension (depth) to be a power of 2
% 
% I have CHANGED THIS FUNCTION to make even and greater than 4.
% 

nf = size(imgs,3);
if mod(nf,2) ~= 0 || nf ~= 2^nextpow2(nf) 
    imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
    imgs = imgs(:,:,1:pow2(ceil(log2(nf))));
elseif nf<4
    imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
end


end