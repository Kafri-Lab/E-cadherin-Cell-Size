function display_solidity_filter_color_overlay(EdgeScore, cyto, cyto_seeds, labelled_cyto, thresh)
    labelled_edge_filt = labelled_cyto;
    for id=1:max(labelled_cyto(:))
        if EdgeScore(id) > thresh
            labelled_edge_filt(labelled_edge_filt==id)=0;
        end
    end
    segmentation_color_overlay(cyto, cyto_seeds, labelled_edge_filt);

    labelled_edge_not_filt = labelled_cyto;
    for id=1:max(labelled_cyto(:))
        if EdgeScore(id) < thresh
            labelled_edge_not_filt(labelled_edge_not_filt==id)=0;
        end
    end
    segmentation_color_overlay(cyto, cyto_seeds, labelled_edge_not_filt);
end