function imgs = make_size_even(imgs)

[R,C,~] = size(imgs);

if mod(R,2) ~= 0
    imgs = imgs(1:end - 1,:,:);
end

if mod(C,2) ~= 0
    imgs = imgs(:,1:end - 1,:);
end