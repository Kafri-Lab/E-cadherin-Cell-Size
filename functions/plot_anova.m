function fun(subsetTable)
  figure
  [p,t,stats] = anova1(subsetTable.CellSize,subsetTable.Animal);
  set(gca,'FontSize',6)
end