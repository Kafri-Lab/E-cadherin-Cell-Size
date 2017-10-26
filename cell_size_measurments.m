set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath('functions'));

ResultsTable = table(); % initialize empty table

%% Load image names 
imgs_path = '\\carbon.research.sickkids.ca\rkafri\DanielS\Images\zoo_animal\hepatocyte_images\Second Image Set\';
thresholded_imgs_path = '\\carbon.research.sickkids.ca\rkafri\DanielS\Images\zoo_animal\hepatocyte_images\old\Phansalkar_threshold_nuc\';
%imgs_path = 'Z:\DanielS\zoo_animal_images\4th set - Miri Stolovich-Rain - animal data from Dors to Kafris lab070617\tif zoo plot\';
img_names = dir([imgs_path '*.tif']);
img_names = {img_names.name}';

%% Map from keywords in filenames to pretty animal names
name_map = get_name_map()

% IMAGE SEGMENTATION SECTION
% for n=1
for n=1:size(img_names,1)
    progress = {img_names{n} 'loop number' n 'out of' size(img_names,1)}  % progress indicator

    %% LOAD IMAGES
    img = imread([imgs_path img_names{n}]);
    cyto = double(img(:,:,1));
    insulin = double(img(:,:,2));
    nuc = double(img(:,:,3));
    % figure; imshow(cyto,[])
    % figure; imshow(insulin,[])
    % figure; imshow(nuc,[])

    %%
    %% Segmentation Section
    %%

    %% Segment Nuc
    nuc_mask = logical(imread([thresholded_imgs_path img_names{n} '_thresh.tif'])); % load precomputed thresholded image
    [nuc_labelled, nuc_seeds] = segment_nuc(nuc, nuc_mask)

    %% Segment Cyto
    [cyto_labelled, cyto_seeds] = segment_cyto(cyto)

    %%
    %% Error Filtering 1
    %%

    % Remove cells which contain no seeds
    for cell_id=1:max(cyto_labelled(:))
        num_seeds = sum(nuc_seeds(cyto_labelled==cell_id));
        if num_seeds == 0
            cyto_labelled(cyto_labelled==cell_id)=0;
        end
    end

    % Remove seeds without cell
    nuc_seeds(~cyto_labelled)=0;
    
    % Remove nuc without cell
    nuc_mask(~cyto_labelled)=0;

    % % Remove nuclei in multiple cells
    % for nuc_id=1:max(labelled_nuc(:))
    %     if length(unique(cyto_labelled(labelled_nuc==nuc_id)))>1
    %        labelled_nuc(labelled_nuc==nuc_id)=0;
    %     end
    % end
    
    % Debug Segmentation
    figure
    subplot(1,2,1);
    segmentation_color_overlay(cyto, nuc_seeds, cyto_labelled, nuc_mask);
    subplot(1,2,2);
    segmentation_color_overlay(nuc, nuc_seeds, cyto_labelled, nuc_mask);
    export_fig(['debug/' img_names{n} '.png'],'-m2')
    close all;

    %%
    %% ResultsTable Section
    %%

    newResults = table();

    % Edge Score
    EdgeScore = calc_edge_score(cyto,cyto_labelled);
    newResults.EdgeScore = EdgeScore
    % Debug edge score
    % display_edge_scole_color_overlay(EdgeScore, cyto, cyto_seeds, cyto_labelled, -0.5);
    
    % Cyto Stats
    cyto_stats=regionprops(cyto_labelled,'Area','Solidity','PixelIdxList');
    newResults.CellSize = cat(1,cyto_stats.Area);
    newResults.Solidity = cat(1,cyto_stats.Solidity);
    % Debug Solidity Filter
    % display_solidity_filter_color_overlay(newResults.Solidity, cyto, cyto_seeds, cyto_labelled, 0.8);

    % Nuc Stats
    % nuc_stats=regionprops(cyto_labelled,nuc_corrected,'MeanIntensity');
    % newResults.DAPI = cat(1,nuc_stats.MeanIntensity);

    % Animal Name Column
    labels = cell(1, length(cyto_stats));
    animal_name = img_name_to_animal_name(img_names{n},name_map);
    labels(:) = {animal_name};
    newResults.Animal = labels';

    % Image Id Column
    labels = cell(1, length(cyto_stats));
    labels(:) = {img_names{n}};
    newResults.Image = labels';

    % Indices Of The Pixels For Each Cyto
    PixelIdxList = cell(1, length(cyto_stats));
    for id=1:max(cyto_labelled(:))
        PixelIdxList{id} = cyto_stats(id).PixelIdxList;
    end
    newResults.PixelIdxList = PixelIdxList';

    % Count Number of Nucs in Cell (By counting number of seeds)
    NucCount = zeros(max(cyto_labelled(:)),1);
    for cell_id=1:max(cyto_labelled(:))
        NucCount(cell_id) = sum(nuc_seeds(cyto_labelled==cell_id));
    end
    newResults.NucCount = NucCount;

    % Measure Nuclear Area per Cell
    NucArea = zeros(max(cyto_labelled(:)),1);
    for cell_id=1:max(cyto_labelled(:))
        NucArea(cell_id) = sum(logical(nuc_mask(cyto_labelled==cell_id)));
        %figure;%nuc_in_cyto=nuc_mask;nuc_in_cyto(cyto_labelled~=cell_id)=0;imshow(nuc_in_cyto);pause
    end
    newResults.NucArea = NucArea;

    % Normalize cell size to number of nuc seeds in cell
    newResults.NormCellSizeSeeds = newResults.CellSize./newResults.NucCount;

    % Normalize cell size to area of nuc(s) in cell
    newResults.NormCellSizeArea = newResults.CellSize./newResults.NucArea;

    % STORE RESULTS
    ResultsTable = [ResultsTable; newResults];
end

% SAVE RESULTS
save('ResultsTable.mat', 'ResultsTable');

% LOAD RESULTS
load('ResultsTable.mat');
subsetTable = ResultsTable; % initialize subset table (more errors will be filtered)


%%
%% Error Filtering 2
%%

% Filter by solidity, edgescore, prctile and more
subsetTable = filter_results_per_image(subsetTable,img_names)

%Filter very big objects
max_cell_size = 6000;
subsetTable=subsetTable(subsetTable.CellSize<6000,:);


%%
%% Debugging Section
%%

% Print All Animal Names With Order Number (for debugging)
for n=1:size(img_names,1)
    fprintf('%s: %s\n',int2str(n),img_names{n});
end

% Save one png per input image with two subplots: one containing all segmentations and one with error segmentations filtered. Colorize by size
segmentation_color_overlay_by_size(ResultsTable,subsetTable,img_names,imgs_path);


%%
%% Analysis Section
%%

%% Load sizes from Yuval's team to compare with
animalsTable = load_animals_table_from_google_spreadsheet();

% % Get Just The Best Images
% images_of_interest = {'Miri Stolovich-Rain - s10-1517600 human 17y.tif', 'Miri Stolovich-Rain - s10-1517601 human 17y.tif', 'Miri Stolovich-Rain - s10-1517603 human 17y.tif', 'Miri Stolovich-Rain - s10-1517604 human 17y.tif', 'Miri Stolovich-Rain - s10-1517605 human 17y.tif', 'Miri Stolovich-Rain - s10-1517606 human 17y.tif', 'Miri Stolovich-Rain - s10-1517607 human 17y.tif', 'Miri Stolovich-Rain - s10-1517608 human 17y.tif', 'Miri Stolovich-Rain - kangaroo liv03.tif', 'Miri Stolovich-Rain - kangaroo liv04.tif', 'Miri Stolovich-Rain - 103 mou 9m.tif', 'Miri Stolovich-Rain - 106 mou 9m.tif', 'Miri Stolovich-Rain - 206 mou 9m.tif', 'Miri Stolovich-Rain - 207 mou 9m.tif', 'Miri Stolovich-Rain - 410 mou 9m.tif', 'Miri Stolovich-Rain - 415 mou 9m.tif'};
% subsetTable = subsetTable(ismember(subsetTable.Image,images_of_interest),:);

% Bar chart
Median = grpstats(subsetTable.CellSize,subsetTable.Animal,'median');
Stds = grpstats(subsetTable.CellSize,subsetTable.Animal,'std');
Lngth = grpstats(subsetTable.CellSize,subsetTable.Animal,'numel');
MedianNormSeeds = grpstats(subsetTable.NormCellSizeSeeds,subsetTable.Animal,'median');
StdNormSeeds = grpstats(subsetTable.NormCellSizeSeeds,subsetTable.Animal,'std');
LngtNormSeeds = grpstats(subsetTable.NormCellSizeSeeds,subsetTable.Animal,'numel');
MedianNormArea = grpstats(subsetTable.NormCellSizeArea,subsetTable.Animal,'median');
StdsNormArea = grpstats(subsetTable.NormCellSizeArea,subsetTable.Animal,'std');
LngthNormArea = grpstats(subsetTable.NormCellSizeArea,subsetTable.Animal,'numel');


%% Add cell sizes computed by computer vision (CV) in this file to the animals table
animalsTable.HepatocyteCV = NaN(height(animalsTable),1);
animalsTable.HepatocyteNormSeeds = NaN(height(animalsTable),1);
animalsTable.HepatocyteNormArea = NaN(height(animalsTable),1);
cv_animal_names = unique(subsetTable.Animal,'stable'); % animals processed by cv
for n=1:length(cv_animal_names)
    animal_index = find(strcmp(animalsTable.ShortName,cv_animal_names{n}));
    animalsTable.HepatocyteCV(animal_index) = Median(n);
    animalsTable.HepatocyteNormSeeds(animal_index) = MedianNormSeeds(n);
    animalsTable.HepatocyteNormArea(animal_index) = MedianNormArea(n);
end

% PRINT ANIMAL CELL SIZES
for n=1:size(Median)
    fprintf('%4.0f: %s\n',Median(n),labels{n});
end

%%
%% Plotting Section
%%

plot_anova(subsetTable);
plot_boxplot(subsetTable,Median);
plot_CellSize_manual_vs_automated(animalsTable);
plot_CellSize_vs_LifeSpan(animalsTable, 'Hepatocyte');
plot_CellSize_vs_LifeSpan(animalsTable, 'HepatocyteCV');


animalsSubsetTable = animalsTable(~isnan(animalsTable.Hepatocyte) & ...
                                  ~isnan(animalsTable.HepatocyteCV) & ...
                                  ~isnan(animalsTable.HepatocyteNormSeeds) & ...
                                  ~isnan(animalsTable.HepatocyteNormArea)  ...
                                 ,:);


% average human data
animalsSubsetTable = average_human_rows(animalsSubsetTable);

% scale between 0 and 1
animalsSubsetTableTEMP = animalsSubsetTable;
animalsSubsetTableTEMP.Acinar = normalize0to1(animalsSubsetTableTEMP.Acinar);
animalsSubsetTableTEMP.Hepatocyte = normalize0to1(animalsSubsetTableTEMP.Hepatocyte);
animalsSubsetTableTEMP.HepatocyteCV = normalize0to1(animalsSubsetTableTEMP.HepatocyteCV);
animalsSubsetTableTEMP.HepatocyteNormArea = normalize0to1(animalsSubsetTableTEMP.HepatocyteNormArea);
animalsSubsetTableTEMP.HepatocyteNormSeeds = normalize0to1(animalsSubsetTableTEMP.HepatocyteNormSeeds);
animalsSubsetTable = animalsSubsetTableTEMP;

