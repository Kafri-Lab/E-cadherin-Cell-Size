function fun(subsetTable, Median, Stds, Lngth)
  figure
  bar(Median)
  labels = unique(subsetTable.Animal,'stable');
  set(gca,'XTickLabel',labels,'XTickLabelRotation',45)
  set(gca,'XTick',1:length(labels));
  hold on
  errorbar(Median,Stds./sqrt(Lngth),'.r')
  ylabel('Cell Area (pixel count)', 'FontSize', 21);
  for i=1:length(Median)
      text(i+0.06,Median(i)+130,int2str(Median(i)),'FontSize',20);
  end
end