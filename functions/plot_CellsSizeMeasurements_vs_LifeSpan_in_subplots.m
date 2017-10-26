function fun(animalsTable)
  animalsSubsetTable = animalsTable(~isnan(animalsTable.Hepatocyte) & ...
                                    ~isnan(animalsTable.HepatocyteCV) & ...
                                    ~isnan(animalsTable.HepatocyteNormSeeds) & ...
                                    ~isnan(animalsTable.HepatocyteNormArea)  ...
                                   ,:);


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
end
