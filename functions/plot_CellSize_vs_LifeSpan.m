function fun(animalsTable,cell_size_column_name)
  %% Plot life span versus cell size (Calculated by Hand)
  animalsSubsetTable = animalsTable(~isnan(eval(['animalsTable.' cell_size_column_name])) & ~isnan(animalsTable.Lifespan),:);
  CellSize = eval(['animalsSubsetTable.' cell_size_column_name]);
  LifeSpan = animalsSubsetTable.Lifespan;

  [r p] = corr(CellSize, log(LifeSpan));
  figure('Position', [400, 400, 300, 250])
  semilogy(CellSize,LifeSpan,'ok')
  set(gca,'ytick',[0 3 6 12 25 50 100])
  [f,fresult]=fit(CellSize,log(LifeSpan),'poly1');
  hold on
  plot(CellSize,exp(f(CellSize)),'r')
  xlabel([cell_size_column_name ' Cell Volume (um^3) Calculated by Hand'],'FontSize',14)
  ylabel('Life Span (yrs)','FontSize',14)
  title(['R = ' num2str(r) ', p = ' num2str(p)])
  %textfit(CellSize,LifeSpan, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
  text(CellSize,LifeSpan, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
  % xlim([400 4500])
end