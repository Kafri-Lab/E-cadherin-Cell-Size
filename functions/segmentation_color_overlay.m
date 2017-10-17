function segmentation_color_overlay(original, seeds, labels, boundary_labels)
    
    % Option 1.
    % Display seeds and segmentation over original image
%     cyto_rgb = cat(3, original, original, original);
%     cyto_seeds_rgb = cat(3, seeds, zeros(size(seeds)), zeros(size(seeds)));
%     cyto_overlay = uint8(labelled_cyto_rgb./4) + uint8(cyto_rgb) + uint8(cyto_seeds_rgb)*128;
%     figure('name','seedrgb','NumberTitle', 'off'); 
%     imshow(uint8(cyto_overlay),[]);
    
    % Option 2.
%     % % Display segmentation over original image
%     figure('name','rgb','NumberTitle', 'off'); imshow(original,[]);
%     hold on
%     labelled_cyto_rgb = label2rgb(uint32(labels), 'jet', [1 1 1], 'shuffle');
%     himage = imshow(labelled_cyto_rgb,[]);
%     himage.AlphaData = 0.3;
%     [xm,ym]=find(seeds);
%     hold on
%     plot(ym,xm,'or','markersize',2,'markerfacecolor','r')
%      
%     % Option 3.
%     % % Display segmentation over original image
%     figure('name','rgb','NumberTitle', 'off'); imshow(cyto,[]);
%     hold on
%     labelled_cyto_rgb = label2rgb(uint32(labelled_cyto), 'jet', [1 1 1], 'shuffle');
%     himage = imshow(labelled_cyto_rgb,[]);
%     himage.AlphaData = 0.3;
    
        
    % Option 4.
    % % Display segmentation over original image and segmented nuclear
    % boundries
    %figure('name','rgb','NumberTitle', 'off'); 
    imshow(imoverlay(uint8(original),bwperim(boundary_labels),'r'),[]);
    hold on
    labelled_cyto_rgb = label2rgb(uint32(labels), 'jet', [1 1 1], 'shuffle');
    himage = imshow(labelled_cyto_rgb,[]);
    himage.AlphaData = 0.3;
    [xm,ym]=find(seeds);
    hold on
    plot(ym,xm,'or','markersize',2,'markerfacecolor','g','markeredgecolor','g')


end