function fun(subsetTable)
  load('ResultsTable.mat');



  %% Load hand collected cell sizes to compare with
  animalsTable = load_animals_table_from_google_spreadsheet();

  % Get Just The Best Images
  %images_of_interest = {'Miri Stolovich-Rain - s10-1517600 human 17y.tif', 'Miri Stolovich-Rain - s10-1517601 human 17y.tif', 'Miri Stolovich-Rain - s10-1517603 human 17y.tif', 'Miri Stolovich-Rain - s10-1517604 human 17y.tif', 'Miri Stolovich-Rain - s10-1517605 human 17y.tif', 'Miri Stolovich-Rain - s10-1517606 human 17y.tif', 'Miri Stolovich-Rain - s10-1517607 human 17y.tif', 'Miri Stolovich-Rain - s10-1517608 human 17y.tif', 'Miri Stolovich-Rain - kangaroo liv03.tif', 'Miri Stolovich-Rain - kangaroo liv04.tif', 'Miri Stolovich-Rain - 103 mou 9m.tif', 'Miri Stolovich-Rain - 106 mou 9m.tif', 'Miri Stolovich-Rain - 206 mou 9m.tif', 'Miri Stolovich-Rain - 207 mou 9m.tif', 'Miri Stolovich-Rain - 410 mou 9m.tif', 'Miri Stolovich-Rain - 415 mou 9m.tif'};
  %subsetTable = ResultsTable(ismember(ResultsTable.Image,images_of_interest),:);

  subsetTable = ResultsTable; % hashtag no filter

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


  animalsSubsetTable = animalsTable(~isnan(animalsTable.HepatocyteNormSeeds) & ...
                                    ~isnan(animalsTable.HepatocyteNormArea)  ...
                                   ,:);


  % average human data
  animalsSubsetTableTEMP = animalsSubsetTable;
  animalsSubsetTableTEMP{1,9:11} = nanmean(animalsSubsetTable{strcmp(animalsSubsetTable.ShortName,'Human'),9:11});
  animalsSubsetTableTEMP([3 9],:) = [];
  animalsSubsetTable = animalsSubsetTableTEMP;

  % % scale between 0 and 1
  % animalsSubsetTableTEMP = animalsSubsetTable;
  % animalsSubsetTableTEMP.Acinar = normalize0to1(animalsSubsetTableTEMP.Acinar);
  % animalsSubsetTableTEMP.Hepatocyte = normalize0to1(animalsSubsetTableTEMP.Hepatocyte);
  % animalsSubsetTableTEMP.HepatocyteCV = normalize0to1(animalsSubsetTableTEMP.HepatocyteCV);
  % animalsSubsetTableTEMP.HepatocyteNormArea = normalize0to1(animalsSubsetTableTEMP.HepatocyteNormArea);
  % animalsSubsetTableTEMP.HepatocyteNormSeeds = normalize0to1(animalsSubsetTableTEMP.HepatocyteNormSeeds);
  % animalsSubsetTable = animalsSubsetTableTEMP;


  AcinarSize = animalsSubsetTable.Acinar;
  HepSize = animalsSubsetTable.Hepatocyte;
  HepSizeCV = animalsSubsetTable.HepatocyteCV;
  HepSizeSeeds = animalsSubsetTable.HepatocyteNormSeeds;
  HepSizeArea = animalsSubsetTable.HepatocyteNormArea;
  LifeSpan = animalsSubsetTable.Lifespan;
  CellSize = [AcinarSize, HepSize, HepSizeCV, HepSizeSeeds, HepSizeArea];
  CellSizeName = {'AcinarSize','HepSize','HepSizeCV','HepSizeSeeds','HepSizeArea'};
  colors = {[0 0 1], [1 0 0], [1 0.7882 0.3529], [0.5 0 0.9], [0 1 0]};

  % figure('Position', [400, 400, 300, 250])
  % for kk = 1:5 
  %     scatter(CellSize(:,kk),log(LifeSpan),120,colors{kk},'LineWidth',2.5); hold on;
  %     [f,fresult]=fit(CellSize(:,kk),log(LifeSpan),'poly1');
  %     ax = plot(f,CellSize(:,kk),log(LifeSpan));
  %     ax(2).Color = colors{kk};    
  %     %errorbar(Median,Stds./sqrt(Lngth),'.r') 
  % end

  figure;
  for ii = 1:5
      subplot(2,3,ii); scatter(CellSize(:,ii),log(LifeSpan),120,colors{ii},'filled'); hold on;
      [f,fresult]=fit(CellSize(:,ii),log(LifeSpan),'poly1');
      plot(f,CellSize(:,ii),log(LifeSpan));
      [r p] = corr(CellSize(:,ii), log(LifeSpan));   
      title(['R = ' num2str(r) ', p = ' num2str(p)])
      xlabel(CellSizeName{ii}); ylabel('log (Lifespan)')
      set(gca,'yticklabel',[0 5 25 50 100]);
  end

  %text(CellSize(:,kk),LifeSpan, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
