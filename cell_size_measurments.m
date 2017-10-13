set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath('functions'));

ResultsTable = table(); % initialize empty table

%% Load image names 
imgs_path = '\\carbon.research.sickkids.ca\rkafri\DanielS\Images\zoo_animal\hepatocyte_images\Second Image Set\';
%imgs_path = 'Z:\DanielS\zoo_animal_images\4th set - Miri Stolovich-Rain - animal data from Dors to Kafris lab070617\tif zoo plot\';
img_names = dir([imgs_path '*.tif']);
img_names = {img_names.name}'; 

%% Map from keywords in filenames to pretty animal names
name_map = containers.Map;
name_map('tiger') = 'Tiger';
name_map('fruit bat') = 'Fruit Bat';
name_map('spalax') = 'Blind Mole Rat';
name_map('wild bore') = 'Wild Bore';
name_map('horse') = 'Horse';
name_map('wolf') = 'Grey Wolf';
name_map('giraf') = 'Giraffe';
name_map('cow') = 'Cow';
name_map('cat') = 'Cat';
name_map('cotton') = 'Cotton Tamarin';
name_map('kangaroo') = 'Kangaroo';
name_map('peccary') = 'Peccary';
name_map('pig') = 'Pig';
name_map('oryx') = 'Arabian Oryx';
name_map('nmr') = 'Naked Mole Rat';
name_map('NMR') = 'Naked Mole Rat';
name_map('prairie dog') = 'Prairie Dog';
% name_map('dog') = 'Dog'; % non-unique name (see prairie dog). So we can't use this as is, need to fix
name_map('zvi') = 'Gazelle';
name_map('RW') = 'Psammomys';
name_map('psamon') = 'Psammomys';
name_map('rat') = 'Black Rat';
name_map('mouse') = 'Mouse';
name_map('mou ') = 'Mouse';
name_map('shrew liv') = 'Shrew';
name_map('shrew 11') = 'Shrew';
name_map('human') = 'Human';
name_map('porcupine') = 'Porcupine';
name_map('dorban') = 'Porcupine';
name_map('mon  ') = 'Macaque';
name_map('grey bat') = 'Grey Bat';

% IMAGE SEGMENTATION SECTION
% for n=1
for n=1:size(img_names,1)
 
 %for n=[24 130]
    progress = {img_names{n} 'loop number' n 'out of' size(img_names,1)}  % progress indicator
    
    %% LOAD IMAGES
    img = imread([imgs_path img_names{n}]);
    cyto = double(img(:,:,1));
    insulin = double(img(:,:,2));
    nuc = double(img(:,:,3));
    % figure; imshow(cyto,[])
    % figure; imshow(insulin,[])
    % figure; imshow(nuc,[])
    
 
   
    %% Nuc Segment
    
    % Correct background brightness
    nuc_closed=imclose(nuc,strel('disk',20)); %imshow(nuc_closed,[]);
    nuc_opened=imopen(nuc_closed,strel('disk',100)); %imshow(nuc_opened,[]);
    nuc_corrected=nuc-nuc_opened; %imshow(nuc_corrected,[]);
    % Smooth
    nuc_smooth = imgaussfilt(nuc_corrected,7);
    % Threshold
    nuc_thresh = nuc_smooth>5; figure; %imshow(nuc_thresh,[]);
    % Remove single isolated pixels
    nuc_open = imopen(nuc_thresh,strel('disk',1)); figure; %imshow(nuc_open,[]);
    % Connect remaining pixels
    nuc_close = imclose(nuc_open,strel('disk',6)); figure; %imshow(nuc_close,[]);
    % Remove smallish groups of pixels
    nuc_open = imopen(nuc_close,strel('disk',5)); figure; %imshow(nuc_open,[]);
    nuc_mask = nuc_open;

    % Find seeds
    nuc_smooth2 = imgaussfilt(nuc_smooth,5);  %imshow(nuc_smooth2,[]);
    nuc_seeds = imregionalmax(nuc_smooth2);
    nuc_seeds = nuc_seeds & nuc_mask;
    
    % Debug nuc seeds
    [xm,ym]=find(nuc_seeds);
    figure; %imshow(nuc_corrected,[prctile(nuc_corrected(:),0) prctile(nuc_corrected(:),95)]); hold on; plot(ym,xm,'or','markersize',2,'markerfacecolor','r')
    
    % Watershed Nuc    
    nuc_min = imimposemin(-nuc_smooth2,nuc_seeds);  %imshow(nuc_min,[]);
    nuc_ws=watershed(nuc_min);
    nuc_ws = nuc_ws & nuc_mask;
    labelled_nuc=bwlabel(nuc_ws);  %imshow(labelled_nuc,[]);

    labelled_nuc = JoinCutNuclei(labelled_nuc); figure; imshow(labelled_nuc,[]);
    
    % NUC STATS
    nuc_stats=regionprops(labelled_nuc,'Area');
    nuc_size = cat(1,nuc_stats.Area);
    
    % Remove Nucs that are too big or too small
    labelled_nuc2 = labelled_nuc;
    min_nuc_size = 0;
    max_nuc_size = 5000;
    too_small_or_big = min_nuc_size>nuc_size | nuc_size >max_nuc_size;
    ids = find(too_small_or_big);
    labelled_nuc2(ismember(labelled_nuc2,ids))=0; imshow(labelled_nuc2,[]);
    labelled_nuc = labelled_nuc2;
    
    % Debug Nuc
    segmentation_color_overlay(nuc_corrected, nuc_seeds, labelled_nuc, labelled_nuc);


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
    labelled_cyto=bwlabel(cyto_ws);
    
    % CLEAR BOARDER
    boarder_cleared = imclearborder(labelled_cyto);
    labelled_cyto = bwlabel(boarder_cleared);
    % figure('name',['boarder_cleared' img_names{n}],'NumberTitle', 'off');imshow(labelled_cyto,[]); colormap(gca, 'jet');
    
    %Debug cyto
    segmentation_color_overlay(cyto, cyto_seeds, labelled_cyto, labelled_nuc);

    %%
    %% ResultsTable Section
    %%
    newResults = table();
    
    %% EDGE SCORE - for filtering segmentation errors (SLOW)
    % BWDIST ON PERIM (to help select edge)
    cyto_perim = bwperim(labelled_cyto);
    cyto_dist = bwdist(cyto_perim);
    % CALC EDGE SCORES
    EDGE_WIDTH_PX = 6;
    EdgeScore = zeros(max(labelled_cyto(:)),1);
    for cell_id=1:max(labelled_cyto(:))
        Cell=cyto(labelled_cyto==cell_id);
        Cell_dist=cyto_dist(labelled_cyto==cell_id);
        Cell_edge=Cell_dist<6; % EDGE WIDTH IN PIXELS
        slope=corrcoef(Cell_dist(Cell_edge),Cell(Cell_edge));
        slope=slope(2,1);
    %     plot(Cell_dist(Cell_edge),Cell(Cell_edge),'.')
    %     title(num2str(slope)); drawnow; pause
        EdgeScore(cell_id)=slope;
    end
    newResults.EdgeScore = EdgeScore;

    % Debug edge score
    % display_edge_scole_color_overlay(EdgeScore, cyto, cyto_seeds, labelled_cyto, -0.5);

    % CYTO STATS
    cyto_stats=regionprops(labelled_cyto,'Area','Solidity','PixelIdxList');
    newResults.CellSize = cat(1,cyto_stats.Area);
    newResults.Solidity = cat(1,cyto_stats.Solidity);

    % Debug Solidity Filter
    % display_solidity_filter_color_overlay(newResults.Solidity, cyto, cyto_seeds, labelled_cyto, 0.8);

    % NUC STATS
    % nuc_stats=regionprops(labelled_cyto,nuc_corrected,'MeanIntensity');
    % newResults.DAPI = cat(1,nuc_stats.MeanIntensity);

    % ANIMAL NAME COLUMN
    labels = cell(1, length(cyto_stats));
    animal_name = img_name_to_animal_name(img_names{n},name_map);
    labels(:) = {animal_name};
    newResults.Animal = labels';

    % IMAGE ID COLUMN
    labels = cell(1, length(cyto_stats));
    labels(:) = {img_names{n}};
    newResults.Image = labels';

    % INDICES OF THE PIXELS FOR EACH CYTO
    PixelIdxList = cell(1, length(cyto_stats));
    for id=1:max(labelled_cyto(:))
        PixelIdxList{id} = cyto_stats(id).PixelIdxList;
    end
    newResults.PixelIdxList = PixelIdxList';
    
    labelled_nuc_orig = labelled_nuc;
    % Remove nuclei in multiple cells
    for nuc_id=1:max(labelled_nuc(:))
        if length(unique(labelled_cyto(labelled_nuc==nuc_id)))>1
           labelled_nuc(labelled_nuc==nuc_id)=0;
        end
    end
    
    % Count Number of Nucs in Cell
    NucCount = zeros(length(cyto_stats),1);
    for cell_id=1:max(labelled_cyto(:))
        NucCount(cell_id) = sum(unique(labelled_nuc(labelled_cyto==cell_id))>0) % node zero here is to ignore background (ie non-nuclei pixels)
        %% OR, the above expanded with very descriptive variable names
        % single_cell = labelled_cyto==cell_id;
        % nuc_pixels_in_single_cell = labelled_nuc(single_cell);
        % unique_pixels_values_in_single_cell_nuc = unique(nuc_pixels_in_single_cell);
        % non_zero_unique_pixels_values_in_single_cell_nuc = unique_pixels_values_in_single_cell_nuc>0;
        % total_non_zero_unique_pixels_values_in_single_cell_nuc = sum(non_zero_unique_pixels_values_in_single_cell_nuc);
    end
    newResults.NucCount = NucCount;
    % STORE RESULTS
    ResultsTable = [ResultsTable; newResults];
end

save('ResultsTable.mat', 'ResultsTable');


%% MEASUREMENTS SECTION

load('ResultsTable.mat');



% Filter by solidity, edgescore, outliers
subsetTable = table();
for n=1:size(img_names,1)
%for n=114 % black rat - rat 4m-304.tif
    imageSubset = ResultsTable(find(strcmp(ResultsTable.Image,img_names{n})),:);
    if height(imageSubset) == 0
        warning('Warning: no cells found in image. Skipping')
        continue
    end

    % Calc Automatic Solidity Threshold (the largest peak in ksdenisty)
    [f,xi] = ksdensity(imageSubset.Solidity);
    [pks, peak_xlocs] = findpeaks(f,'SortStr','descend','NPeaks',1); % find the largest peak
    spacing = xi(2)-xi(1); % Store spacing between sliding window jumps of ksdensity
    peak_xlocs = ((peak_xlocs)*spacing)+min(xi)-spacing; % Set correct range of x values.
    solidity_threshold = peak_xlocs(1); % Use the first (largest) peak location (on x-axis) as the edge score threshold
    solidity_threshold = solidity_threshold - std(imageSubset.Solidity)/3; % adjust threshold by 1/3 a std dev from the mean
    
    % Calc Automatic Edge Score Threshold (the largest peak in ksdenisty)
    [f,xi] = ksdensity(imageSubset.EdgeScore);
    [pks, peak_xlocs] = findpeaks(f,'SortStr','descend','NPeaks',1); % find the largest peak
    spacing = xi(2)-xi(1); % Store spacing between sliding window jumps of ksdensity
    peak_xlocs = ((peak_xlocs)*spacing)+min(xi)-spacing; % Set correct range of x values.
    edge_threshold = peak_xlocs(1); % Use the first (largest) peak location (on x-axis) as the edge score threshold
    

    imageSubset=imageSubset(imageSubset.Solidity>solidity_threshold,:);
    imageSubset=imageSubset(imageSubset.EdgeScore<edge_threshold,:);
    
    % prctile_threshold = prctile(newSubset.CellSize,95);
    %newSubset=newSubset(newSubset.CellSize<prctile_threshold,:);
    % prctile_threshold = prctile(newSubset.CellSize,5);
    %newSubset=newSubset(newSubset.CellSize>prctile_threshold,:);
    
    % debug
    %segmentation_color_overlay(cyto, cyto_seeds, labelled_cyto);

    if height(imageSubset) == 0
        warning('Warning: thresholds removed all cells')
        continue
    end
    subsetTable = [subsetTable; imageSubset];
end

%Filter very big objects
max_cell_size = 6000;
subsetTable=subsetTable(subsetTable.CellSize<6000,:);

% Filter by insulin
% % TODO: Remove loop
% insulinTable = table();
% insulin_threshold = 777;
% for n=1:size(img_names,1)
%     subset_ids=subsetTable.Insulin>insulin_threshold;
%     newSubset=subsetTable(subset_ids,:);
%     newSubset = newSubset(find(strcmp(newSubset.Image,img_names{n})),:);
%     insulinTable = [insulinTable; newSubset];
% end
% % Filter by no insulin
% %TODO: Remove loop
% noinsulinTable = table();
% insulin_threshold = 777;
% for n=1:size(img_names,1)
%     subset_ids=subsetTable.Insulin<insulin_threshold;
%     newSubset=subsetTable(subset_ids,:);
%     newSubset = newSubset(find(strcmp(newSubset.Image,img_names{n})),:);
%     noinsulinTable = [noinsulinTable; newSubset];
% end


% ALL ANIMAL NAMES WITH ORDER NUMBER (for debugging)
for n=1:size(img_names,1)
    fprintf('%s: %s\n',int2str(n),img_names{n});
end


%% Load hand collected cell sizes to compare with
animalsTable = load_animals_table_from_google_spreadsheet();



%% GRAPHICS SECTION

% RGB Segmentation Overlay (subsetTable)
% for n=1
for n=1:size(img_names,1)
%for n=114 % black rat - rat 4m-304.tif
    n
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
    %export_fig('1.png','-m2')
    % print(gcf,['filtered/' char(unique(img_subsetTable.Animal)) '_' img_names{n} '_cellsize' int2str(mean(img_subsetTable.CellSize)) '.png'],'-dpng','-r300');
    close all
end


%% PLOTTING SECTION
% Anova table and boxplot
figure
[p,t,stats] = anova1(subsetTable.CellSize,subsetTable.Animal);
set(gca,'FontSize',6)


% Bar chart
Means = grpstats(subsetTable.CellSize,subsetTable.Animal,'mean');
Stds = grpstats(subsetTable.CellSize,subsetTable.Animal,'std');
Lngth = grpstats(subsetTable.CellSize,subsetTable.Animal,'numel');
figure
bar(Means)
labels = unique(subsetTable.Animal,'stable');
set(gca,'XTickLabel',labels,'XTickLabelRotation',45)
set(gca,'XTick',1:length(labels));
hold on
errorbar(Means,Stds./sqrt(Lngth),'.r')
ylabel('Cell Area (pixel count)', 'FontSize', 21);
for i=1:length(Means)
    text(i+0.06,Means(i)+130,int2str(Means(i)),'FontSize',20);
end

% PRINT ANIMAL CELL SIZES
for n=1:size(Means)
    fprintf('%4.0f: %s\n',Means(n),labels{n});
end

%% Add cell sizes computed by computer vision (CV) in this file to the animals table
animalsTable.HepatocyteCV = NaN(height(animalsTable),1);
cv_animal_names = unique(subsetTable.Animal,'stable'); % animals processed by cv
for n=1:length(cv_animal_names)
    animal_index = find(strcmp(animalsTable.ShortName,cv_animal_names{n}));
    animalsTable.HepatocyteCV(animal_index) = Means(n);
end

%% Plot cell size calculated by hand vs computer vision
animalsSubsetTable = animalsTable(~isnan(animalsTable.Hepatocyte) & ...
                          ~isnan(animalsTable.HepatocyteCV),:);
CellSize = animalsSubsetTable.Hepatocyte;
CellSizeCV = animalsSubsetTable.HepatocyteCV;
figure('Position', [400, 400, 300, 250])
scatter(CellSizeCV,CellSize, 'ok')
ylabel('Hepatocyte Cell Volume Calculated by Hand','FontSize',14)
xlabel('Hepatocyte Cell Volume Calculated by Computer Vision','FontSize',14)
textfit(CellSizeCV,CellSize, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
[f,fresult]=fit(CellSizeCV,CellSize,'poly1');
hold on
plot(CellSizeCV,f(CellSizeCV),'r')
[r p] = corr(CellSizeCV,CellSize);
title(['R = ' num2str(r) ', p = ' num2str(p)])
% xlim([500 4500])
% ylim([500 4500])

%% Plot life span versus cell size (Calculated by Hand)
animalsSubsetTable = animalsTable(~isnan(animalsTable.Hepatocyte) & ...
                            ~isnan(animalsTable.Lifespan),:);
CellSize = animalsSubsetTable.Hepatocyte;
LifeSpan = animalsSubsetTable.Lifespan;

[r p] = corr(CellSize, log(LifeSpan));
figure('Position', [400, 400, 300, 250])
semilogy(CellSize,LifeSpan,'ok')
set(gca,'ytick',[0 3 6 12 25 50 100])
[f,fresult]=fit(CellSize,log(LifeSpan),'poly1');
hold on
plot(CellSize,exp(f(CellSize)),'r')
xlabel(['Hepatocyte Cell Volume (um^3) Calculated by Hand'],'FontSize',14)
ylabel('Life Span (yrs)','FontSize',14)
title(['R = ' num2str(r) ', p = ' num2str(p)])
%textfit(CellSize,LifeSpan, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
text(CellSize,LifeSpan, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
% xlim([400 4500])

%% Plot life span versus cell size (Calculated by CV)
animalsSubsetTable = animalsTable(~isnan(animalsTable.HepatocyteCV) & ...
                            ~isnan(animalsTable.Lifespan),:);
CellSize = animalsSubsetTable.HepatocyteCV;
LifeSpan = animalsSubsetTable.Lifespan;

[r p] = corr(CellSize, log(LifeSpan));
figure('Position', [400, 400, 300, 250])
semilogy(CellSize,LifeSpan,'ok')
set(gca,'ytick',[0 3 6 12 25 50 100])
[f,fresult]=fit(CellSize,log(LifeSpan),'poly1');
hold on
plot(CellSize,exp(f(CellSize)),'r')
xlabel(['Hepatocyte Cell Volume (um^3) Calculated by Computer Vision'],'FontSize',14)
ylabel('Life Span (yrs)','FontSize',14)
title(['R = ' num2str(r) ', p = ' num2str(p)])
text(CellSize,LifeSpan, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
% xlim([400 4500])
