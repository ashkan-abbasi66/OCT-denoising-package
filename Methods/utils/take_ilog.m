function img_out = take_ilog(img_log)
% Take inverse logarithm transformation
% 
% See "take_log"
    img_out = exp(img_log) - 1;
end