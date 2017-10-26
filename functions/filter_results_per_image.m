function subsetTable = fun(ResultsTable,img_names)
    % Filter by solidity, edgescore, prctile and more
    subsetTable = table();
    for n=1:size(img_names,1)
        imageSubset = ResultsTable(find(strcmp(ResultsTable.Image,img_names{n})),:);
        if height(imageSubset) == 0
            warning('Warning: no cells found in image. Skipping')
            continue
        end

        %% Calc Automatic Solidity Threshold (the largest peak in ksdenisty)
        [f,xi] = ksdensity(imageSubset.Solidity);
        [pks, peak_xlocs] = findpeaks(f,'SortStr','descend','NPeaks',1); % find the largest peak
        spacing = xi(2)-xi(1); % Store spacing between sliding window jumps of ksdensity
        peak_xlocs = ((peak_xlocs)*spacing)+min(xi)-spacing; % Set correct range of x values.
        solidity_threshold = peak_xlocs(1); % Use the first (largest) peak location (on x-axis) as the edge score threshold
        solidity_threshold = solidity_threshold - std(imageSubset.Solidity)/3; % adjust threshold by 1/3 a std dev from the mean
        
        %% Calc Automatic Edge Score Threshold (the largest peak in ksdenisty)
        [f,xi] = ksdensity(imageSubset.EdgeScore);
        [pks, peak_xlocs] = findpeaks(f,'SortStr','descend','NPeaks',1); % find the largest peak
        spacing = xi(2)-xi(1); % Store spacing between sliding window jumps of ksdensity
        peak_xlocs = ((peak_xlocs)*spacing)+min(xi)-spacing; % Set correct range of x values.
        edge_threshold = peak_xlocs(1); % Use the first (largest) peak location (on x-axis) as the edge score threshold

        % Filter by solidity
        imageSubset=imageSubset(imageSubset.Solidity>solidity_threshold,:);
        % Filter by Edge Score
        imageSubset=imageSubset(imageSubset.EdgeScore<edge_threshold,:);
        
        %% Filter by prctile
        % prctile_threshold = prctile(newSubset.CellSize,95);
        %newSubset=newSubset(newSubset.CellSize<prctile_threshold,:);
        % prctile_threshold = prctile(newSubset.CellSize,5);
        %newSubset=newSubset(newSubset.CellSize>prctile_threshold,:);
        
        % Debug
        %segmentation_color_overlay(cyto, cyto_seeds, cyto_labelled);

        if height(imageSubset) == 0
            warning('Warning: thresholds removed all cells')
            continue
        end
        subsetTable = [subsetTable; imageSubset];
    end
end