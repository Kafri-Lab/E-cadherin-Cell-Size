function [cyto_labelled, cyto_seeds] = fun(cyto)
    %% Cyto Section
    %%
    % figure('name',['cyto' img_names{n}],'NumberTitle', 'off');imshow(cyto,[])
    
    %% SMOOTH
    cyto_smooth = imgaussfilt(cyto,7);
    %figure('name',['cyto_smooth' img_names{n}],'NumberTitle', 'off');imshow(cyto_smooth,[])
    
    %% FIND SEEDS
    cyto_smooth=imhmin(cyto_smooth,2); % suppresing local minima
    % figure('name',['imhmin' img_names{n}],'NumberTitle', 'off');imshow(cyto_smooth,[])
    [cyto_seeds]=imregionalmin(cyto_smooth);
 
    
    % Debug cyto seeds
    [xm,ym]=find(cyto_seeds);
    % figure; imshow(cyto,[]); hold on; plot(ym,xm,'or','markersize',2,'markerfacecolor','r')
    
    %% WATERSHED
    cyto_min = imimposemin(cyto,cyto_seeds);
    cyto_ws=watershed(cyto_min);
    cyto_labelled=bwlabel(cyto_ws);
    
    % CLEAR BOARDER
    boarder_cleared = imclearborder(cyto_labelled);
    cyto_labelled = bwlabel(boarder_cleared);
    % figure('name',['boarder_cleared' img_names{n}],'NumberTitle', 'off');imshow(cyto_labelled,[]); colormap(gca, 'jet');
    
    %Debug cyto
    %segmentation_color_overlay(cyto, nuc_seeds, cyto_labelled, nuc_mask);

end