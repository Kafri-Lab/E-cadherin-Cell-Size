function segmentation_color_overlay(cyto, cyto_seeds, labelled_cyto)
    labelled_cyto_rgb = label2rgb(labelled_cyto,'jet', 'k', 'shuffle');

    % % Display seeds and segmentation over original image
    % cyto_rgb = cat(3, cyto, cyto, cyto);
    % cyto_seeds_rgb = cat(3, cyto_seeds, zeros(size(cyto_seeds)), zeros(size(cyto_seeds)));
    % cyto_overlay = uint8(labelled_cyto_rgb./4) + uint8(cyto_rgb) + uint8(cyto_seeds_rgb)*128;
    % figure('name','seedrgb','NumberTitle', 'off'); imshow(uint8(cyto_overlay),[]);
    
    % % Display segmentation over original image
    figure('name','rgb','NumberTitle', 'off'); imshow(cyto,[]);
    hold on
    labelled_cyto_rgb = label2rgb(uint32(labelled_cyto), 'jet', [1 1 1], 'shuffle');
    himage = imshow(labelled_cyto_rgb,[]); himage.AlphaData = 0.3;
end