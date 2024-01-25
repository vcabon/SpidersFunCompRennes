plot.nmds <- function (nmds.obj, 
                       sp.fit, 
                       pval.sp = 0.05,
                       env.fit,
                       pval.env = 0.05,
                       env.var, 
                       ssp.fit = FALSE) 
{
    
    
    site.scrs <- as.data.frame(scores(nmds.obj, display = "sites"))
    site.scrs <- cbind(site.scrs, Site = rownames(site.scrs))
    site.scrs <- cbind(site.scrs, Var = spiders_env[[paste0(env.var)]])
    
    spp.scrs <- as.data.frame(scores(sp.fit, display = "vectors"))
    spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
    spp.scrs <- cbind(spp.scrs, pval = sp.fit[["vectors"]][["pvals"]])
    sig.spp.scrs <- subset(spp.scrs, pval<=pval.sp)
    
    env.scrs <- as.data.frame(scores(env.fit, display = "vectors"))
    env.scrs <- cbind(env.scrs, env.variables = rownames(env.scrs))
    env.scrs <- cbind(env.scrs, pval = env.fit[["vectors"]][["pvals"]])
    sig.env.scrs <- subset(env.scrs, pval<=pval.env)
    
    
    
    if (ssp.fit == TRUE) {
        
        nmds.plot <- ggplot(site.scrs, 
                            aes(x=NMDS1,
                                y=NMDS2)) +
            
            geom_point(size = 1.5, 
                       aes(NMDS1, 
                           NMDS2, 
                           shape = site.scrs$Var)) +
                         #colour = site.scrs$Var, size = 1.5)) +
            
         #   scale_colour_gradient(low = "lightgrey", high = "black") +
            
            guides(size = "none") +
            
            coord_fixed() +
            
            theme_classic() + 
            
            theme(panel.background = element_rect(fill = NA, 
                                                  colour = "black", 
                                                  size = 1, 
                                                  linetype = "solid")) +
            
            labs(colour = env.var) + 
            
            theme(legend.position = "right", 
                  legend.text = element_text(size = 12), 
                  legend.title = element_text(size = 12), 
                  axis.text = element_text(size = 10)) + 
            
            geom_segment(data = sig.spp.scrs, 
                         aes(x = 0, 
                             xend=NMDS1, 
                             y=0, yend=NMDS2), 
                         arrow = arrow(length = unit(0.25, "cm")), 
                         colour = "grey10", 
                         lwd=0.3) + 
            
            ggrepel::geom_text_repel(data = sig.spp.scrs, 
                                     aes(x=NMDS1, 
                                         y=NMDS2,label = Species), 
                                     cex = 3,
                                     direction = "both",
                                     segment.size = 0.25) +
            
            labs(title="NMDS - species")
    }
    
    else {
        
        nmds.plot <- ggplot(site.scrs, 
                            aes(x=NMDS1,
                                y=NMDS2)) +
            
            geom_point(size = 2, 
                       aes(NMDS1, 
                           NMDS2, 
                           shape = site.scrs$Var,
                           colour = site.scrs$Var)) +
            
            scale_colour_manual(values = c("#137C8B", "#B8CBD0", "#595959")) +
            
           # scale_colour_gradient(low = "lightgrey", high = "black") +
            
            guides(size = "none") +
            
            coord_fixed() +
            
            theme_classic() + 
            
            theme(panel.background = element_rect(fill = NA, 
                                                  colour = "black", 
                                                  size = 1, 
                                                  linetype = "solid")) +
            
            labs(colour = env.var) + 
            
            theme(legend.position = "right", 
                  legend.text = element_text(size = 12), 
                  legend.title = element_text(size = 12), 
                  axis.text = element_text(size = 10)) + 
            
            geom_segment(data = sig.env.scrs, 
                         aes(x = 0, xend=NMDS1, 
                             y=0, yend=NMDS2), 
                         arrow = arrow(length = unit(0.25, "cm")),
                         lwd=0.3) + 
            
            ggrepel::geom_text_repel(data = sig.env.scrs, 
                                     aes(x=NMDS1, y=NMDS2, 
                                         label = env.variables), 
                                     cex = 3, direction = "both", 
                                     segment.size = 0.25) +
            
            labs(title=" ", shape = "Cluster")
        
    }
}


