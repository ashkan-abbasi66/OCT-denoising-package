

%% show magnified ROI for an image

addpath('E:\THESIS\Implements\Pedagogical\BEST ROUTINES');


n=7;% image number  --- 7
scale_factor=1;
%===============================================

% ~~~~~~~~~~~~ reconstructed images ~~~~~~~~~~~~
% experiment_name=sprintf('bm4d_bicubic_%dxInterp',scale_factor)
% % output file name (Here, an image is expected)
% outfolder=['outs/' experiment_name];
outfolder=['../../3D_OCT_denoising_toolbox_Results/benchmark_bm4d_log_dt1'];
outfname=fullfile(outfolder,sprintf('%02d.tif',n));
% outfname=fullfile(outfolder,sprintf('%0.6d.tif',(n-1)*5+1));
% outfname=fullfile(outfolder,sprintf('%dSSR_denoising_result.tif',n));
% % ~~~~~~~~~~~~ Original images ~~~~~~~~~~~~
% outfname=fullfile('../Datasets/dt1_Bioptigen_SDOCT/',...
%     num2str(n),'test.tif');


%===============================================
% Where ROIs were saved? (Here, a .mat file is expected)
posfname=['./rois_for_dt1/im' num2str(n)];
load(posfname,'pos');
%
im_out=imread(outfname);
%=============
if ~isempty(strfind(outfname,'test')) && (scale_factor>1)
    im_out=im_out(:,1:scale_factor:end);
    [R,C]=size(im_out);
    C=C*scale_factor;
    valid_cols=uint16(1:scale_factor:C);
    im_out2=zeros(R,C);
    im_out2(:,valid_cols)=mat2gray(im_out);
    im_out=im_out2;% Now, size(imn)==size(imn2);
    clear im_out2
end
%=============
figure,imshow(im_out,[])
%
% Show ROIs on the image (6 ROIs)
selected_roi=1:6;
roi_for_magnification=1:6
for j=selected_roi
    %imrect(fh.CurrentAxes,pos{j});
    rectangle('Position',pos{j},...
        'EdgeColor','r','LineWidth',1)    
end
% Magnify ROIs
for j=roi_for_magnification
    cropped_im=imcrop(mat2gray(im_out),pos{j});
    figure
    imshow(cropped_im,'InitialMagnification',250)
    title(sprintf('ROI #%d',j));
end


%% show outputs with or without ROIs

% % % iptsetpref('ImshowBorder','loose')
% % % show_ROI=1;
% % % testImg_indices=3;
% % % fh=figure;
% % % for i=1:numel(testImg_indices)
% % %     strnumber=num2str(testImg_indices(i));
% % %     % ** the results are saved into this file
% % %     outfolder=['../../3D_OCT_denoising_toolbox_Results/benchmark_bm4d_log_dt1'];
% % %     outfname=fullfile(outfolder,sprintf('%02d.tif',n));
% % %     im_out=imread(outfname);
% % %     % o_sb3D_msomp_im
% % %     % out_myinp_x2_K52_nl3_tukey_im
% % %     %
% % %     pth=['E:\DTSET\Fang2013\For synthetic experiments\',strnumber];
% % %     cleanfile='average.tif';testfile='test.tif';
% % %     im= single(imread(fullfile(pth,cleanfile)));
% % %     imn = single(imread(fullfile(pth,testfile)));
% % %     %
% % %     imshow(im_out,[])
% % %     if show_ROI==1
% % %         posfname=['./rois_for_dt1/im' num2str(n)];
% % %         load(posfname,'pos');
% % %         for j=1:numel(pos)
% % %             imrect(fh.CurrentAxes,pos{j});
% % %         end
% % %     end
% % %     pause
% % % end


%% crop image #3   ---- local PSNR & SSIM
[PSNR,SSIM]=comp_psnr(imcrop(im,pos),imcrop(im_out,pos),imcrop(imn,pos));
pos=[1 160 890 60];% [xmin ymin width height]  --
figure,imshow(imcrop(im_out,pos),[]);
title(sprintf('PSNR=%.4g, SSIM=%.4g',PSNR,SSIM))
%% show outputs of two method
cd('E:\mfiles_acode_thesis\005_ksvd\Inpainting_Experiment\my_inpainting_HHdict');
testImg_indices=[3];%1:18;
for i=1:numel(testImg_indices)
    strnumber=num2str(testImg_indices(i));
    % ** first method
    outfname1=['outs/out_myinp_x2_K52_nl3_tukey_im' strnumber];% o_minp_ompmod_nlm_im
    % ** second method
    outfname2=sprintf('outs/o_sb%s_msomp_im%s','3D',strnumber);% o_sb3D_msomp_im
    % ** Load images
    load(outfname1,'im_out');
    im_out1=im_out;
    load(outfname2,'im_out');
    im_out2=im_out;
    %
    fh1=figure;
    imshow(im_out1,[]);title([outfname1 ': image #' strnumber],'interpreter', 'none')
    fh2=figure;
    imshow(im_out2,[]);title([outfname2 ': image #' strnumber],'interpreter', 'none')
    pause
    fh1.delete;
    fh2.delete;
end
%% show outputs of two methods in ROI
cd('E:\mfiles_acode_thesis\005_ksvd\Inpainting_Experiment\my_inpainting_HHdict');
testImg_indices=3;%1:18;
for i=1:numel(testImg_indices)
    strnumber=num2str(testImg_indices(i));
    % ** first method
    outfname1=['outs/out_myinp_x2_K52_nl3_tukey_im' strnumber];% o_minp_ompmod_nlm_im
    % ** second method
    outfname2=sprintf('outs/o_sb%s_msomp_im%s','3D',strnumber);
    % ** Load images
    load(outfname1,'im_out');
    im_out1=im_out;
    [R,C]=size(im_out1);
    load(outfname2,'im_out');
    im_out2=im_out;
    %
    fh1=figure;
    imshow(im_out1,[]);title([outfname1 ': image #' strnumber],'interpreter', 'none')
    hrec=imrect;
    wait(hrec);
    pos=hrec.getPosition;
    x=pos(1);
    y=pos(2);
    h=pos(4);
    w=pos(3);
    rows=round(y+1:y+h);
    columns=round(x+1:x+w);
    rows(rows<1 & rows>R)=[];
    columns(columns<1 & columns>C)=[];
    imshow(im_out1(rows,columns),[]);title([outfname1 ': image #' strnumber],'interpreter', 'none')
    fh2=figure;
    imshow(im_out2(rows,columns),[]);title([outfname2 ': image #' strnumber],'interpreter', 'none')
    pause
    fh1.delete;
    fh2.delete;
end
