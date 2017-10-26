function fun(ResultsTable,subsetTable,img_names,imgs_path)
    % Save one png per input image with two subplots: one containing all segmentations and one with error segmentations filtered. Colorize by size
    for n=1:size(img_names,1)
        img_subsetTable = subsetTable(find(strcmp(subsetTable.Image,img_names{n})),:);
        if height(img_subsetTable) == 0
            continue
        end
        img = imread([imgs_path img_names{n}]);
        cyto = double(img(:,:,1));
        labelled_by_size = zeros(size(cyto));
        for i=1:height(img_subsetTable)
            PixelIdxList = cell2mat(img_subsetTable{i,{'PixelIdxList'}});
            labelled_by_size(PixelIdxList)=img_subsetTable{i,'CellSize'};
        end
        labelled_by_size(labelled_by_size>2999)=2999; % make colors more beautiful by putting an upper limit
        labelled_by_size_mod_colors = labelled_by_size;
        labelled_by_size_mod_colors(1)=min(subsetTable{:,'CellSize'});
        labelled_by_size_mod_colors(2)=2999;
        % Display RGB overlay
        figure('name',[img_names{n}],'NumberTitle', 'off'); 
        subplot(1,2,2);
        imshow(cyto,[]);
        hold on
        labelled_by_size_rgb = label2rgb(uint16(labelled_by_size_mod_colors), 'jet', [1 1 1]);
        himage = imshow(labelled_by_size_rgb,[]); himage.AlphaData = 0.3;
        
        subplot(1,2,1);
        labelled_by_size = zeros(size(cyto));
        img_ResultsTable = ResultsTable(find(strcmp(ResultsTable.Image,img_names{n})),:);
        for i=1:height(img_ResultsTable)
            PixelIdxList = cell2mat(img_ResultsTable{i,{'PixelIdxList'}});
            labelled_by_size(PixelIdxList)=img_ResultsTable{i,'CellSize'};
        end
        labelled_by_size(labelled_by_size>2999)=2999; % make colors more beautiful by putting an upper limit
        labelled_by_size_mod_colors = labelled_by_size;
        labelled_by_size_mod_colors(1)=min(ResultsTable{:,'CellSize'});
        labelled_by_size_mod_colors(2)=2999;
        % Display RGB overlay
        imshow(cyto,[]);
        hold on
        labelled_by_size_rgb = label2rgb(uint16(labelled_by_size_mod_colors), 'jet', [1 1 1]);
        himage = imshow(labelled_by_size_rgb,[]); himage.AlphaData = 0.3;
        
        %title([char(unique(img_subsetTable.Animal)) ' - Filename: "' img_names{n} '", Cell size: ' int2str(mean(img_subsetTable.CellSize)) 'px'],'Interpreter','none')
        
        export_fig(['segmentation_errors_filtered/' char(unique(img_subsetTable.Animal)) '_' img_names{n} '_cellsize' int2str(mean(img_subsetTable.CellSize)) '.png'],'-m2')
        close all
    end
end