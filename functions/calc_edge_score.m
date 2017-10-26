function EdgeScore = fun(cyto,labelled_cyto)
    %% EDGE SCORE - for filtering segmentation errors (SLOW)
    % BWDIST ON PERIM (to help select edge)
    cyto_perim = bwperim(labelled_cyto);
    cyto_dist = bwdist(cyto_perim);
    % CALC EDGE SCORES
    EDGE_WIDTH_PX = 6;
    EdgeScore = zeros(max(labelled_cyto(:)),1);
    for cell_id=1:max(labelled_cyto(:))
        Cell=cyto(labelled_cyto==cell_id);
        Cell_dist=cyto_dist(labelled_cyto==cell_id);
        Cell_edge=Cell_dist<6; % EDGE WIDTH IN PIXELS
        slope=corrcoef(Cell_dist(Cell_edge),Cell(Cell_edge));
        slope=slope(2,1);
    %     plot(Cell_dist(Cell_edge),Cell(Cell_edge),'.')
    %     title(num2str(slope)); drawnow; pause
        EdgeScore(cell_id)=slope;
    end
end