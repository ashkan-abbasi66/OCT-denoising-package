function imshow_roi(img,pos)
    figure
    imshow(img)
    
    N = length(pos);
    
    for i=1:N
        rectangle('Position',pos{i},...
        'EdgeColor','r','LineWidth',1)    
    end
    
end