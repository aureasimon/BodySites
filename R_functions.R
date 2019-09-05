#07/15/2019: Ceylan Tannes Function sto run after QIIME2 first analysis
### ================
###   User defined functions
### ================
filter_low_coverage <- function(props, frac_cutoff=0.6, min_ab=0){
  frac_nonzero <- function (x) sum(x > min_ab) / length(x)
  apply(props, 1, frac_nonzero) >= frac_cutoff
}
###=====
###  make_pcoa_plot <- function(uu, s, shape_by, color_by, title)
###  uu: distance, s: mapping file, shape_by: variable used for shape, color_by: variable used for color
###=====
make_pcoa_plot <- function(dm, s, shape_by, color_by) {
  dm <- usedist::dist_subset(dm, s$SampleID)
  pc <- pcoa(dm)
  pc_df <- merge(s, pc$vectors[, 1:3], by.x="SampleID", by.y="row.names")
  pc_pct <- round(pc$values$Relative_eig * 100)
  
  pcoa_plot = ggplot(pc_df, aes(x=Axis.1, y=Axis.2)) +
    theme_bw() +
    scale_shape_discrete(name=sub("_", " ", shape_by)) + 
    scale_colour_discrete(name=sub("_", " ", color_by)) +
    labs(
      x=paste0("PCoA axis 1 (", pc_pct[1], "%)"),
      y=paste0("PCoA axis 2 (", pc_pct[2], "%)")
    )
  
  if (is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by))))
  } else if (!is.null(shape_by) & !is.null(color_by)) {
    pcoa_plot <- pcoa_plot + geom_point(aes(colour=factor(get(color_by)), shape=factor(get(shape_by))))
  } else {
    pcoa_plot <- pcoa_plot + geom_point()
  }
  return(pcoa_plot)
}
###=====
###  heatmap_grouped <- function(genus_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1)
###  option=1: rows_to_keep <- filter_low_coverage(heatmap_props, perc_cutoff=thre) ## taxa found in at least 80% of samples
###  option=2: rows_to_keep <- apply(heatmap_props,1,max) >= 0.01 ## taxa with abundance in any sample exceeding 1%
###=====
heatmap_grouped <- function(summed_props, heatmap_s, grps = c("study_group", "study_day"), fname=NULL, thre=0.8, option=1, prop_cut=0.01, satu_limit=0.4){
  
  #color = saturated_rainbow(101)
  color = saturated_rainbow(101, saturation_limit=satu_limit)
  breaks = c(0, 1e-10, seq(0.001, 1, length.out = 100))
  
  heatmap_props <- summed_props[,heatmap_s$SampleID]
  
  if (option == 1) {
    rows_to_keep <- filter_low_coverage(heatmap_props, frac_cutoff=thre) 
  } else if (option == 2) {
    rows_to_keep <- apply(heatmap_props,1,max) >= prop_cut 
  }
  heatmap_props <- heatmap_props[rows_to_keep,]
  
  ## group the SampleIDs
  heatmap_s %<>% arrange_(.dots=grps)
  heatmap_props <- heatmap_props[, heatmap_s$SampleID]
  
  ## update the annotation
  annc <- heatmap_s[,grps] %>% as.data.frame()
  rownames(annc) <- heatmap_s$SampleID
  colnames(annc) <- grps
  
  ## heatmap time
  if (!is.null(fname))
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, filename = fname, 
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
  else
    pheatmap(heatmap_props, annotation = annc, color = color, breaks = breaks, 
             fontsize_col = 8, fontsize_row = 8, cluster_cols = FALSE, cluster_rows = FALSE,cellheight = 8, cellwidth = 8)
}
