set(0,'DefaultFigureWindowStyle','docked')
addpath 'functions';
ResultsTable = table(); % initialize empty table

%% Load image names
imgs_path = 'images\4th set - Miri Stolovich-Rain - animal data from Dors to Kafris lab070617\tif zoo plot\';
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
name_map('dog') = 'Dog';
name_map('oryx') = 'Arabian Oryx';
name_map('nmr') = 'Naked Mole Rat';
name_map('NMR') = 'Naked Mole Rat';
name_map('prairie dog') = 'Prairie Dog';
name_map('zvi') = 'Gazelle';
name_map('RW') = 'Psammomys';
name_map('rat') = 'Black Rat';
name_map('mouse') = 'Mouse';
name_map('shrew') = 'Shrew';
name_map('human') = 'Human';
name_map('porcupine') = 'Porcupine';
name_map('mon  pan') = 'Macaque';
name_map('grey bat') = 'Grey Bat';

count = 1;

% for n=1
for n=1:size(img_names,1)
    progress = {img_names{n} 'loop number' count 'out of' size(img_names,1)}  % progress indicator
    count = count+1;

    %% LOAD IMAGES
    img = imread([imgs_path img_names{n}]);
    cyto = double(img(:,:,1));
    insulin = double(img(:,:,2));
    % nuc = double(img(:,:,3));
    % figure; imshow(cyto,[])

    %%
    %% Insulin Section
    %%
    % figure('name',['insulin' img_names{n}],'NumberTitle', 'off');imshow(insulin,[])
    
    %% SMOOTH
    ins_smooth = imgaussfilt(insulin,12);
    % figure('name',['ins_smooth' img_names{n}],'NumberTitle', 'off');imshow(ins_smooth,[])
    
    %% THRESHOLD
    thresh = calc_insulin_threshold(ins_smooth);
    ins_thresh = ins_smooth>thresh;
    % figure('name',['ins_thresh' img_names{n}],'NumberTitle', 'off');imshow(ins_thresh,[])
    
    %% FILL HOLES
    ins_fill = imfill(ins_thresh, 'holes');
    % figure('name',['ins_fill' img_names{n}],'NumberTitle', 'off');imshow(ins_fill,[])
    
    %% ERODE (compensate for aggresive threshold)
    ins_erode = imerode(ins_fill, strel('disk',21));
    % figure('name',['ins_erode' img_names{n}],'NumberTitle', 'off');imshow(ins_erode,[])
    
    %% REMOVE SMALL OBJECTS
    ins_open = bwareaopen(ins_erode, 3000);
    % figure('name',['ins_open' img_names{n}],'NumberTitle', 'off');imshow(ins_open,[])

    insulin_mask = ins_open;

    %%
    %% Cyto Section
    %%
    % figure('name',['cyto' img_names{n}],'NumberTitle', 'off');imshow(cyto,[])
    
    %% SMOOTH
    cyto_smooth = imgaussfilt(cyto,7);
    % figure('name',['cyto_smooth' img_names{n}],'NumberTitle', 'off');imshow(cyto_smooth,[])
    
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
    % figure('name',['boarder_cleared' img_names{n}],'NumberTitle', 'off');imshow(labelled_cyto,[]); colormap(gca, 'jet');
    
    % % Debug cyto
    % labelled_cyto_rgb = label2rgb(labelled_cyto,'jet', 'k', 'shuffle');
    % cyto_rgb = cat(3, cyto, cyto, cyto);
    % cyto_seeds_rgb = cat(3, cyto_seeds, zeros(size(cyto_seeds)), zeros(size(cyto_seeds)));
    % cyto_overlay = uint8(labelled_cyto_rgb./4) + uint8(cyto_rgb) + uint8(cyto_seeds_rgb)*128;
    % figure('name',['seedrgb' img_names{n}],'NumberTitle', 'off'); imshow(uint8(cyto_overlay),[]);
    
    % figure('name',['rgb' img_names{n}],'NumberTitle', 'off'); imshow(cyto,[]);
    % hold on
    % labelled_cyto_rgb = label2rgb(uint32(labelled_cyto), 'jet', [1 1 1], 'shuffle');
    % himage = imshow(labelled_cyto_rgb,[]); himage.AlphaData = 0.3;

    %%
    %% ResultsTable Section
    %%
    newResults = table();
    
    % EDGE SCORE (for filtering segmentation errors)
    EdgeScore = [];
    for id=1:max(labelled_cyto(:))
        Per=bwperim(labelled_cyto==id);
        DilPer=imdilate(Per,strel('disk',2));
        Tube=imdilate(Per,strel('disk',4))&~DilPer;
        EdgeScore(id)=mean(cyto(DilPer))/mean(cyto(Tube));
    end
    newResults.EdgeScore = EdgeScore';
    
    % CYTO STATS
    cyto_stats=regionprops(labelled_cyto,'Area','Solidity','PixelIdxList');
    newResults.CellSize = cat(1,cyto_stats.Area);
    newResults.Solidity = cat(1,cyto_stats.Solidity);

    % INSULIN STATS
    ins_stats=regionprops(labelled_cyto,insulin_mask.*7777,'MeanIntensity'); % the .*7777 is a temporary solution to set cells as their beta or not beta, a continuous solution should be implemented when time permits
    newResults.Insulin = cat(1,ins_stats.MeanIntensity);

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
    
    % STORE RESULTS
    ResultsTable = [ResultsTable; newResults];
end

save('ResultsTable.mat', 'ResultsTable');


%% MEASUREMENTS SECTION

load('ResultsTable.mat');



% Filter by solidity, insulin, prctile (keep up to 90 prctile)
%TODO: Remove loop
subsetTable = table();
solidity_threshold = 0.7;
insulin_threshold = 777;
for n=1:size(img_names,1)
    newSubset = ResultsTable(find(strcmp(ResultsTable.Image,img_names{n})),:);
    newSubset=newSubset(newSubset.Solidity>solidity_threshold,:);
    newSubset=newSubset(newSubset.Insulin<insulin_threshold,:);
    prctile_threshold = prctile(newSubset.CellSize,90);
    newSubset=newSubset(newSubset.CellSize<prctile_threshold,:);
    subsetTable = [subsetTable; newSubset];
end

% Filter by edge score
edge_threshold = 1.3;
subsetTable=subsetTable(subsetTable.EdgeScore>edge_threshold,:);

% Filter very big objects
max_cell_size = 10000;
subsetTable=subsetTable(subsetTable.CellSize<10000,:);

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



%% GRAPHICS SECTION

% RGB Segmentation Overlay (subsetTable)
% for n=1
for n=1:size(img_names,1)
    img = imread([imgs_path img_names{n}]);
    cyto = double(img(:,:,1));
    labelled_by_size = zeros(size(cyto));
    img_subsetTable = subsetTable(find(strcmp(subsetTable.Image,img_names{n})),:);
    for i=1:height(img_subsetTable)
        PixelIdxList = cell2mat(img_subsetTable{i,{'PixelIdxList'}});
        labelled_by_size(PixelIdxList)=img_subsetTable{i,'CellSize'};
    end
    labelled_by_size(labelled_by_size>4444)=4444; % make colors more beautiful by putting an upper limit
    labelled_by_size_mod_colors = labelled_by_size;
    labelled_by_size_mod_colors(1)=min(subsetTable{:,'CellSize'});
    labelled_by_size_mod_colors(2)=4444;
    % Display RGB overlay
    figure('name',[img_names{n}],'NumberTitle', 'off'); imshow(cyto,[]);
    hold on
    labelled_by_size_rgb = label2rgb(uint32(labelled_by_size_mod_colors), 'jet', [1 1 1]);
    himage = imshow(labelled_by_size_rgb,[]); himage.AlphaData = 0.3;
    print(gcf,['filtered/' img_names{n} '.png'],'-dpng','-r300');
    close all
end

% RGB Segmentation Overlay (ResultsTable)
% for n=1
for n=1:size(img_names,1)
    img = imread([imgs_path img_names{n}]);
    cyto = double(img(:,:,1));
    labelled_by_size = zeros(size(cyto));
    img_ResultsTable = ResultsTable(find(strcmp(ResultsTable.Image,img_names{n})),:);
    for i=1:height(img_ResultsTable)
        PixelIdxList = cell2mat(img_ResultsTable{i,{'PixelIdxList'}});
        labelled_by_size(PixelIdxList)=img_ResultsTable{i,'CellSize'};
    end
    labelled_by_size(labelled_by_size>4444)=4444; % make colors more beautiful by putting an upper limit
    labelled_by_size_mod_colors = labelled_by_size;
    labelled_by_size_mod_colors(1)=min(ResultsTable{:,'CellSize'});
    labelled_by_size_mod_colors(2)=4444;
    % Display RGB overlay
    figure('name',[img_names{n}],'NumberTitle', 'off'); imshow(cyto,[]);
    hold on
    labelled_by_size_rgb = label2rgb(uint32(labelled_by_size_mod_colors), 'jet', [1 1 1]);
    himage = imshow(labelled_by_size_rgb,[]); himage.AlphaData = 0.3;
    print(gcf,['all/' img_names{n} '.png'],'-dpng','-r300');
    close all
end

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
hold on
errorbar(Means,Stds./sqrt(Lngth),'.r')
ylabel('Cell Area (pixel count)', 'FontSize', 21);
for i=1:length(Means)
    text(i+0.06,Means(i)+130,int2str(Means(i)),'FontSize',20);
end
%set(gca,'FontSize',19,'XTickLabel',{'Image 1', 'Image 2', 'Image 3', 'Image 4', 'Image 5', 'Image 6', 'Image 7', 'Image 8'})


%% Load hand collected cell sizes to compare with
animalsTable = load_animals_table_from_google_spreadsheet();

%% Add cell sizes computed by computer vision (CV) in this file to the animals table
animalsTable.AcinarCV = NaN(height(animalsTable),1);
cv_animal_names = unique(subsetTable.Animal); % animals processed by cv
for n=1:length(cv_animal_names)
    animal_index = find(strcmp(animalsTable.ShortName,cv_animal_names{n}));
    animalsTable.AcinarCV(animal_index) = Means(n);
end


%% Plot cell size calculated by hand vs computer vision
animalsSubsetTable = animalsTable(~isnan(animalsTable.Acinar) & ...
                          ~isnan(animalsTable.AcinarCV),:);
CellSize = animalsSubsetTable.Acinar;
CellSizeCV = animalsSubsetTable.AcinarCV;

figure('Position', [400, 400, 300, 250])
scatter(CellSizeCV,CellSize, 'ok')
ylabel('Acinar Cell Volume Calculated by Computer Vision')
xlabel('Acinar Cell Volume Calculated by Hand')


%% Plot life span versus cell size Calculated by Hand
animalsSubsetTable = animalsTable(~isnan(animalsTable.Acinar) & ...
                            ~isnan(animalsTable.Lifespan),:);
CellSize = animalsSubsetTable.Acinar;
LifeSpan = animalsSubsetTable.Lifespan;

[r p] = corr(CellSize, log(LifeSpan));
figure('Position', [400, 400, 300, 250])
semilogy(CellSize,LifeSpan,'ok')
set(gca,'ytick',[0 3 6 12 25 50 100])
[f,fresult]=fit(CellSize,log(LifeSpan),'poly1');
hold on
plot(CellSize,exp(f(CellSize)),'r')
xlabel(['Acinar Cell Volume (um^3) Calculated by Hand'])
ylabel('Life Span (yrs)')
title(['R = ' num2str(r) ', p = ' num2str(p)])


%% Plot life span versus cell size Calculated by Hand
animalsSubsetTable = animalsTable(~isnan(animalsTable.AcinarCV) & ...
                            ~isnan(animalsTable.Lifespan),:);
CellSize = animalsSubsetTable.AcinarCV;
LifeSpan = animalsSubsetTable.Lifespan;

[r p] = corr(CellSize, log(LifeSpan));
figure('Position', [400, 400, 300, 250])
semilogy(CellSize,LifeSpan,'ok')
set(gca,'ytick',[0 3 6 12 25 50 100])
[f,fresult]=fit(CellSize,log(LifeSpan),'poly1');
hold on
plot(CellSize,exp(f(CellSize)),'r')
xlabel(['Acinar Cell Volume (um^3) Calculated by Computer Vision'])
ylabel('Life Span (yrs)')
title(['R = ' num2str(r) ', p = ' num2str(p)])

