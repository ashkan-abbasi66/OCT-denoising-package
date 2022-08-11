function img_log = take_log(img)
% Take logarithm transformation of the input image
% 
% USAGE:
%   img_log = take_log(img)
%   Do some operation on "img_log"
%   img_out = take_ilog(img_log)
% 
    img_log = log(img+1);
end