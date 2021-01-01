#' Visualize cell evfs, colored using the simulated time
#' @param filename name of the saved figure
#' @param evf_res evf parameters
#' @param width width of the figure
#' @param height height of the figure
#' @param units select from "cm", "mm" and "in", default "in"
#' @param dpi dpi value of the figure
#' @return None
#' @import umap ggplot2
#' @export
plotEVF <- function(evf_res, width = 15, height = 15, units = "in", dpi = 1000){
  library(ggplot2)
  param_names <- c("k_on", "k_off", "s")
  kinetic_params <- evf_res[[1]]

  theme<-theme(panel.background = element_blank(),
               panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))

  lapply(param_names, function(i){
    pca_res <- prcomp(kinetic_params[[i]], scale. = F)
    x_dr <- as.data.frame(pca_res$x)
    x_dr$cell_time <- 1:dim(kinetic_params[[i]])[1]
    plt <- ggplot() + geom_point(data = x_dr, mapping = aes(x = PC1, y = PC2, color = cell_time)) +
      theme + scale_color_gradient(low = "orange", high = "purple")
    ggsave(filename = paste0(i, "_kinetic_pca.pdf"), plot = plt, device = "pdf", path = "./", height = height, width = width, units = units, dpi = dpi)
  })
}

#' Visualize the result, colored using the simulated time
#' @param filename name of the saved figure
#' @param result output from the simulator
#' @param dr dimension reduction method, default "pca
#' @param width width of the figure
#' @param height height of the figure
#' @param units select from "cm", "mm" and "in", default "in"
#' @param dpi dpi value of the figure
#' @return None
#' @import umap ggplot2
#' @export
plotPseudotime <- function(filename, result, dr = "pca", width = 15, height = 15, units = "in", dpi = 1000){
  library(ggplot2)

  counts_u <- result$counts_u
  counts_s <- result$counts_s
  kinetic_params <- result$kinet_params
  state_mat <- result$state_mat
  cell_time <- result$cell_time
  velocity <- result$velocity

  if(dr == "pca") {
    pca_res <- prcomp(log1p(t(counts_s)), scale. = F)
    x_dr <- as.data.frame(pca_res$x)
    x_dr$cell_time <- cell_time
  } else if(dr == "umap") {
    umap_res <- umap(d=t(log1p(counts_s)))
    x_dr <- as.data.frame(umap_res$layout)
    colnames(x_dr) <- c("PC1", "PC2")
    x_dr$cell_time <- cell_time
  } else {stop("The dr can only be umap or pca")}

  theme<-theme(panel.background = element_blank(),
               panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))

  plt <- ggplot(data = x_dr) + geom_point(mapping = aes(x = PC1, y = PC2, color = cell_time)) + theme + scale_color_gradient(low = "orange", high = "purple")
  ggsave(filename = filename, plot = plt, device = "pdf", path = "./", height = height, width = width, units = units, dpi = dpi)
}

#' Visualize the simulated velocity, colored using the simulated time
#' @param filename name of the saved figure
#' @param result output from the simulator
#' @param dr dimension reduction method, default "pca
#' @param width width of the figure
#' @param height height of the figure
#' @param units select from "cm", "mm" and "in", default "in"
#' @param dpi dpi value of the figure
#' @return None
#' @import umap ggplot2
#' @export
plotVelo <- function(filename, result, width = 15, height = 15, units = "in", dpi = 1000, arrow.length = 0.01){
  library(ggplot2)

  counts_u <- result$counts_u
  counts_s <- result$counts_s
  kinetic_params <- result$kinet_params
  state_mat <- result$state_mat
  cell_time <- result$cell_time
  velocity <- result$velocity

  pca <- prcomp(log1p(t(counts_s)), scale. = F)
  new_mat <- counts_s + arrow.length * velocity; new_mat[which(new_mat<0)] <- 0
  pca.pred <- predict(pca, newdata = log1p(t(new_mat)))

  pca_df <- as.data.frame(pca$x)
  pca_df$cell_time <- cell_time
  pca.pred_df <- as.data.frame(pca.pred)
  colnames(pca.pred_df) <- lapply(colnames(pca.pred_df), function(x){paste0(x,".pred")})
  pca_res <- cbind(pca_df, pca.pred_df)

  theme<-theme(panel.background = element_blank(),
               panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))

  plt <- ggplot(data = pca_res) + geom_point(mapping = aes(x = PC1, y = PC2, color = cell_time)) + theme + scale_color_gradient(low = "orange", high = "purple") +
    geom_segment(aes(x = PC1, y = PC2, xend = PC1.pred, yend =  PC2.pred), arrow = arrow(length = unit(0.1, "cm")), color = "black", alpha = 0.5)

  ggsave(filename = filename, plot = plt, device = "pdf", path = "./", height = height, width = width, units = units, dpi = dpi)
}
