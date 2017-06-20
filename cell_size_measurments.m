set(0,'DefaultFigureWindowStyle','docked')
addpath '\\carbon.research.sickkids.ca\rkafri\Miriam\Matlab function library'
inpath = '\\carbon.research.sickkids.ca\rkafri\DanielS\yuval_blind_measurements\images\';
ResultsTable = table(); % initialize empty table

img_ids = [1 2 3 4 5 6 7 8];

for img_id=img_ids
    progress = ['loop ' int2str(img_id)]
    
    %% LOAD IMAGES
    img = imread([inpath int2str(img_id) '.tif']);
    cyto = double(img(:,:,1));
    insulin = double(img(:,:,2));
    nuc = double(img(:,:,3));
    
    %%
    %% Insulin Section
    %%
    % figure('name',['insulin' int2str(img_id)],'NumberTitle', 'off');imshow(insulin,[])
    
    %% SMOOTH
    ins_smooth = imgaussfilt(insulin,12);
    % figure('name',['ins_smooth' int2str(img_id)],'NumberTitle', 'off');imshow(ins_smooth,[])
    
    %% THRESHOLD
    ins_thresh = ins_smooth>14;
    % figure('name',['ins_thresh' int2str(img_id)],'NumberTitle', 'off');imshow(ins_thresh,[])
    
    %% FILL HOLES
    ins_fill = imfill(ins_thresh, 'holes');
    % figure('name',['ins_fill' int2str(img_id)],'NumberTitle', 'off');imshow(ins_fill,[])
    
    %% ERODE (compensate for aggresive threshold)
    ins_erode = imerode(ins_fill, strel('disk',10));
    % figure('name',['ins_erode' int2str(img_id)],'NumberTitle', 'off');imshow(ins_erode,[])
    
    %% REMOVE SMALL OBJECTS
    ins_open = bwareaopen(ins_erode, 500);
    % figure('name',['ins_open' int2str(img_id)],'NumberTitle', 'off');imshow(ins_open,[])
    
    insulin_mask = ins_open;
    
    
    %%
    %% Nuclei Section
    %%
    % figure('name',['nuc' int2str(img_id)],'NumberTitle', 'off');imshow(nuc,[])
    
    %% SMOOTH
    nuc_smooth = imgaussfilt(nuc,7);
    % figure('name',['nuc_smooth' int2str(img_id)],'NumberTitle', 'off');imshow(nuc_smooth,[])
    
    %% Flat-field correction
    nuc_closed=imclose(nuc_smooth,strel('disk',15));
    nuc_opened=imopen(nuc_closed,strel('disk',150));
    nuc_corrected=nuc-nuc_opened;
    % figure('name',['corrected' int2str(img_id)],'NumberTitle', 'off');imshow(nuc_corrected,[])
    
    %% THRESHOLD
    nuc_thresh = nuc_corrected>15;
    % figure('name',['nuc_thresh' int2str(img_id)],'NumberTitle', 'off');imshow(nuc_thresh,[])
    
    %% OPEN
    nuc_open = imopen(nuc_thresh,strel('disk',3));
    % figure('name',['nuc_open' int2str(img_id)],'NumberTitle', 'off');imshow(nuc_open,[])
    
    %% REMOVE SMALL OBJECTS
    nuc_open = bwareaopen(nuc_open, 100);
    % figure('name',['open' int2str(img_id)],'NumberTitle', 'off');imshow(nuc_open,[])
    
    %% RENAME
    nuc_clean = nuc_open;
    
    %% FIND SEEDS
    nuc_seeds=imregionalmax(nuc_smooth);
    nuc_seeds=nuc_seeds&nuc_clean;
    nuc_seeds=imdilate(nuc_seeds,strel('sphere',7));
    
    % Debug nuc seeds
    nuc_seed_overlay=(RGBOverlay(mat2gray(nuc+100),nuc_seeds));
    % figure('name',['nuc_seeds' int2str(img_id)],'NumberTitle', 'off');imshow(nuc_seed_overlay,[])
    
    %% WATERSHED
    nuc_ws=watershed(imimposemin(imcomplement(nuc),nuc_seeds))&nuc_clean;
    labelled_nuc = bwlabel(nuc_ws);
    
    % Debug nuc
    labelled_nuc_rgb = label2rgb(labelled_nuc,'jet', 'k', 'shuffle');
    nuc_seeds_rgb = cat(3, nuc_seeds, zeros(size(nuc_seeds)), zeros(size(nuc_seeds)));
    nuc_rgb = cat(3, nuc, nuc, nuc);
    nuc_overlay = uint8(labelled_nuc_rgb./4) + uint8(nuc_rgb+50) + uint8(nuc_seeds_rgb)*128;
    % figure;imshow(uint8(nuc_overlay),[]);
    cyto_rgb = cat(3, cyto, cyto, cyto);
    cyto_nuc_overlay = uint8(labelled_nuc_rgb./4) + uint8(cyto_rgb) + uint8(nuc_seeds_rgb)*128;
    % figure;imshow(uint8(cyto_nuc_overlay),[]);
    
    
    %%
    %% Cyto Section
    %%
    % figure('name',['cyto' int2str(img_id)],'NumberTitle', 'off');imshow(cyto,[])
    
    %% SMOOTH
    cyto_smooth = imgaussfilt(cyto,7);
    % figure('name',['cyto_smooth' int2str(img_id)],'NumberTitle', 'off');imshow(cyto_smooth,[])
    
    %% FIND SEEDS
    cyto_smooth=imhmin(cyto_smooth,2); % suppresing local minima
    [cyto_seeds]=imregionalmin(cyto_smooth);
    
    % Debug cyto seeds
    [xm,ym]=find(cyto_seeds);
    % figure; imshow(cyto,[]); hold on; plot(ym,xm,'or','markersize',2,'markerfacecolor','r')
    
    %% WATERSHED
    cyto_min = imimposemin(cyto,cyto_seeds);
    cyto_ws=watershed(cyto_min);
    labelled_cyto=bwlabel(cyto_ws);
    
    % CLEAR BOARDER
    boarder_cleared = imclearborder(labelled_cyto);
    labelled_cyto = bwlabel(boarder_cleared);
    % figure('name',['boarder_cleared' int2str(img_id)],'NumberTitle', 'off');imshow(labelled_cyto,[]); colormap(gca, 'jet');
    
    % Debug cyto
    labelled_cyto_rgb = label2rgb(labelled_cyto,'jet', 'k', 'shuffle');
    cyto_rgb = cat(3, cyto, cyto, cyto);
    cyto_seeds_rgb = cat(3, cyto_seeds, zeros(size(cyto_seeds)), zeros(size(cyto_seeds)));
    cyto_overlay = uint8(labelled_cyto_rgb./4) + uint8(cyto_rgb) + uint8(cyto_seeds_rgb)*128;
    figure('name',['seedrgb' int2str(img_id)],'NumberTitle', 'off'); imshow(uint8(cyto_overlay),[]);
    
    figure('name',['rgb' int2str(img_id)],'NumberTitle', 'off'); imshow(cyto,[]);
    hold on
    labelled_cyto_rgb = label2rgb(uint32(labelled_cyto), 'jet', [1 1 1], 'shuffle');
    himage = imshow(labelled_cyto_rgb,[]); himage.AlphaData = 0.3;
    
    %%
    %% ResultsTable Section
    %%
    
    % EDGE SCORE
    EdgeScore = [];
    for id=1:max(labelled_cyto(:))
        Per=bwperim(labelled_cyto==id);
        DilPer=imdilate(Per,strel('disk',2));
        Tube=imdilate(Per,strel('disk',4))&~DilPer;
        EdgeScore(id)=mean(cyto(DilPer))/mean(cyto(Tube));
    end
    
    % MISC METRICS
    cyto_stats=regionprops(labelled_cyto,'Area','Solidity','PixelIdxList');
    nuc_stats=regionprops(labelled_cyto,nuc_corrected,'MeanIntensity');
    ins_stats=regionprops(labelled_cyto,insulin,'MeanIntensity');
    stats_matrix=[cat(1,cyto_stats.Area)...
        cat(1,cyto_stats.Solidity)...
        cat(1,nuc_stats.MeanIntensity)...
        cat(1,ins_stats.MeanIntensity)...
        EdgeScore'];
    
    newResults=array2table(stats_matrix,'VariableNames',{'CellSize','Solidity','DAPI','Insulin','EdgeScore'});
    
    % IMAGE ID COLUMN
    labels = cell(1, length(cyto_stats));
    labels(:) = {['Image ' int2str(img_id)]};
    newResults.Image = labels';
    
    % INDICES OF THE PIXELS FOR EACH CYTO
    PixelIdxList = cell(1, length(cyto_stats));
    for id=1:max(labelled_cyto(:))
        PixelIdxList{id} = cyto_stats(id).PixelIdxList;
    end
    newResults.PixelIdxList = PixelIdxList';
    
    % STORE RESULTS
    ResultsTable = [ResultsTable; newResults];
end

save('ResultsTable.mat', 'ResultsTable');


%% MEASUREMENTS SECTION

load('ResultsTable.mat');

% Filter by solidity
solidity_thresholds = [0.75 0.75 0.75 0.75 0.85 0.75 0.85 0.75];
subsetTable = table();
for img_id=img_ids
    newSubset = ResultsTable(find(strcmp(ResultsTable.Image,{['Image ', int2str(img_id)]})),:);
    subset_ids=newSubset.Solidity>solidity_thresholds(img_id);
    newSubset=newSubset(subset_ids,:);
    subsetTable = [subsetTable; newSubset];
end

% Filter very big objects
max_cell_size = 10000;
subsetTable=subsetTable(subsetTable.CellSize<10000,:);

% Filter by insulin
insulin_thresholds = [12 4.6 40 14 33 31 9 9];
insulinTable = table();
for img_id=img_ids
    subset_ids=subsetTable.Insulin>insulin_thresholds(img_id);
    newSubset=subsetTable(subset_ids,:);
    newSubset = newSubset(find(strcmp(newSubset.Image,{['Image ', int2str(img_id)]})),:);
    insulinTable = [insulinTable; newSubset];
end


% Filter by insulin
insulin_thresholds = [12 4.6 40 14 33 31 9 9];
noinsulinTable = table();
for img_id=img_ids
    subset_ids=subsetTable.Insulin<insulin_thresholds(img_id);
    newSubset=subsetTable(subset_ids,:);
    newSubset = newSubset(find(strcmp(newSubset.Image,{['Image ', int2str(img_id)]})),:);
    noinsulinTable = [noinsulinTable; newSubset];
end




%% GRAPHICS SECTION
subsetTable = noinsulinTable;

% RGB Segmentation Overlay
for img_id=img_ids
    img = imread([inpath int2str(img_id) '.tif']);
    cyto = double(img(:,:,1));
    labelled_by_size = zeros(size(cyto));
    img_subsetTable = subsetTable(find(strcmp(subsetTable.Image,{['Image ', int2str(img_id)]})),:);
    for i=1:height(img_subsetTable)
        PixelIdxList = cell2mat(img_subsetTable{i,{'PixelIdxList'}});
        labelled_by_size(PixelIdxList)=img_subsetTable{i,'CellSize'};
    end
    labelled_by_size(labelled_by_size>7777)=7777; % make colors more beautiful by putting an upper limit
    labelled_by_size_mod_colors = labelled_by_size;
    labelled_by_size_mod_colors(1)=min(subsetTable{:,'CellSize'});
    labelled_by_size_mod_colors(2)=7777;
    % Display RGB overlay
    figure('name',['rgb' int2str(img_id)],'NumberTitle', 'off'); imshow(cyto,[]);
    hold on
    labelled_by_size_rgb = label2rgb(uint32(labelled_by_size_mod_colors), 'jet', [1 1 1]);
    himage = imshow(labelled_by_size_rgb,[]); himage.AlphaData = 0.3;
end

% Anova table and boxplot
figure
[p,t,stats] = anova1(subsetTable.CellSize,subsetTable.Image);
set(gca,'FontSize',19)


% Bar chart
Means = grpstats(subsetTable.CellSize,subsetTable.Image,'mean');
Stds = grpstats(subsetTable.CellSize,subsetTable.Image,'std');
Lngth = grpstats(subsetTable.CellSize,subsetTable.Image,'numel');
figure
bar(Means)
hold on
errorbar(Means,Stds./sqrt(Lngth),'.r')
ylabel('Cell Area (pixel count)', 'FontSize', 21);
for i=1:length(Means)
    text(i+0.06,Means(i)+130,int2str(Means(i)),'FontSize',20);
end
set(gca,'FontSize',19,'XTickLabel',{'Image 1', 'Image 2', 'Image 3', 'Image 4', 'Image 5', 'Image 6', 'Image 7', 'Image 8'})


% VIOLIN PLOT
addpath '\\carbon.research.sickkids.ca\rkafri\DanielS\Violinplot-Matlab'
figure('Position', [100, 100, 900, 850]);
vs = violinplot(subsetTable.CellSize, subsetTable.Image);
set(gca,'FontSize',19)
ylabel('Cell Area (pixel count)', 'FontSize', 21);


