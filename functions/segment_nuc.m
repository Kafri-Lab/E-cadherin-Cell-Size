function [nuc_labelled, nuc_seeds] = fun(nuc, nuc_mask)
    %% Nuc Segment
    % Correct background brightness
    nuc_closed=imclose(nuc,strel('disk',20)); %imshow(nuc_closed,[]);
    nuc_opened=imopen(nuc_closed,strel('disk',100)); %imshow(nuc_opened,[]);
    nuc_corrected=nuc-nuc_opened; %imshow(nuc_corrected,[]);

    % Smooth
    nuc_smooth = imgaussfilt(nuc_corrected,7);

    % Threshold
    % nuc_mask = imread([thresholded_imgs_path img_names{n} '_thresh.tif']);
    % figure; imshow(nuc_mask,[]);

    % Find seeds (using thresh+bwdist)
    nuc_mask = imclose(nuc_mask,strel('disk',5)); figure; imshow(nuc_mask,[]);
    nuc_bwdist = bwdist(~nuc_mask); %figure; imshow(nuc_bwdist,[0 50])
    nuc_bwdist = imgaussfilt(nuc_bwdist,8);
    nuc_seeds = imregionalmax(nuc_bwdist);
    nuc_seeds = nuc_seeds & nuc_mask;

    nuc_labelled = []

    % Debug nuc seeds
    %[xm,ym]=find(nuc_seeds);
    %figure; imshow(nuc_corrected,[prctile(nuc_corrected(:),0) prctile(nuc_corrected(:),95)]); hold on; plot(ym,xm,'or','markersize',2,'markerfacecolor','r')

    % Not watershedding because it causes cut problems!
    % Watershed Nuc
    % nuc_min = imimposemin(-nuc_smooth2,nuc_seeds);  %imshow(nuc_min,[]);
    % nuc_ws=watershed(nuc_min);
    % nuc_ws = nuc_ws & nuc_mask;
    % nuc_labelled=bwlabel(nuc_ws);  %figure; imshow(nuc_labelled,[]);

    % Join Cut Nuclei
    %nuc_labelled = JoinCutNuclei(nuc_labelled); %figure; imshow(nuc_labelled,[]);

    % % NUC STATS
    % nuc_stats=regionprops(nuc_labelled,'Area');
    % nuc_size = cat(1,nuc_stats.Area);

    % % Remove Nucs that are too big or too small
    % nuc_labelled_orig = nuc_labelled;
    % % nuc_labelled = nuc_labelled_orig; % to aid development
    % min_nuc_size = 0;
    % max_nuc_size = 5000;
    % too_small_or_big = min_nuc_size>nuc_size | nuc_size >max_nuc_size;
    % ids = find(too_small_or_big);
    % nuc_labelled(ismember(nuc_labelled,ids))=0; %imshow(nuc_labelled,[]);

    % Debug Nuc
    %segmentation_color_overlay(nuc_corrected, nuc_seeds, nuc_labelled, nuc_labelled);
    %pause
    %continue
end