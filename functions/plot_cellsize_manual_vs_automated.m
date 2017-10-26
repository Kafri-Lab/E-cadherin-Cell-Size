function fun(animalsTable)
  animalsSubsetTable = animalsTable(~isnan(animalsTable.HepatocyteCV) & ...
                            ~isnan(animalsTable.Hepatocyte),:);
  size_manual = animalsSubsetTable.Hepatocyte;
  size_automated = animalsSubsetTable.HepatocyteCV;
  figure('Position', [400, 400, 300, 250])
  scatter(size_automated,size_manual, 'ok')
  ylabel('Hepatocyte Cell Volume Calculated by Hand','FontSize',14)
  xlabel('Hepatocyte Cell Volume Calculated by Computer Vision','FontSize',14)
  textfit(size_automated,size_manual, animalsSubsetTable.ShortName, 'horizontal','left', 'vertical','bottom','FontSize',11)
  [f,fresult]=fit(size_automated,size_manual,'poly1');
  hold on
  plot(size_automated,f(size_automated),'r')
  [r p] = corr(size_automated,size_manual);
  title(['R = ' num2str(r) ', p = ' num2str(p)])
  % xlim([500 4500])
  % ylim([500 4500])
end