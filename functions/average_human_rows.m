function animalsTable = fun(animalsTable)
  animalsTableTEMP = animalsTable;
  animalsTableTEMP{1,9:11} = nanmean(animalsTable{strcmp(animalsTable.ShortName,'Human'),9:11});
  animalsTableTEMP([2 5],:) = [];
  animalsTable = animalsTableTEMP;
end