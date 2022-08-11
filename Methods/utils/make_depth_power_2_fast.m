function imgs = make_depth_power_2_fast(imgs)
% 
% Make the 3rd dimension (depth) to be a power of 2
% 
% I have CHANGED THIS FUNCTION to make even and greater than 4.
% 

%  MAKE EVEN AND GREATER THAN 4
nf = size(imgs,3);

if mod(nf,2) ~= 0 
    imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
    imgs = imgs(:,:,1:nf+1);
end

if nf<4
    imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
end

end