## Figure state transition

### PACKAGES ####
require(diagram)
require(labdsv)
require(vegan)

### DATA ####
path_data <- "../Doctorat/data/"

# Species
sp.mat.ba <- readRDS(paste0(path_data, "species/sp_mat_ba_jul2018.RDS"))
Trans_df <- readRDS(paste0(path_data, "species/Trans_df_1980-2000.rds"))
gr_vec <- readRDS(paste0(path_data, "species/gr_ba3_vec.rds"))


# Spatial
ecoregion <- st_read(paste0(path_data, "map/ecoregion_simple.gpkg"), quiet = T)
ecoregion <- st_transform(ecoregion, 32198)

ecoregion$SOUS_DOM6 <- factor(ecoregion$SOUS_DOM6, rev(c("Spruce-moss",
                                                         "Balsam fir-white birch",
                                                         "Balsam fir-yellow birch",
                                                         "Sugar maple-yellow birch",
                                                         "Sugar maple-basswood",
                                                         "Sugar maple-bitternut hickory")))



xy <- st_read(paste0(path_data, "species/plot_xy32198_may2018.gpkg"))
st_crs(xy) <- 32198
xy <- subset(xy, plot_id %in% unique(Trans_df$plot_id))


# Plotting parameters

gr.name <- c("Boreal", "Mixed", "Pioneer", "Temperate")

pal.gr <- c("#158282", "#A1BD93","#FEAC19", "#D43650")


### MAP ####


# quartz(width = 6.2, height = 3.5)
png(paste0("images/map_plots.png"), res=300, 
    width = 6.2, height = 3.5, units = 'in', bg = "transparent", type='cairo')
par(mar=c(0.1,0,0,0),las=1)
plot(st_geometry(ecoregion), border = "gray60", 
     col = alpha(brewer.pal(6,"Spectral"),.5)[ecoregion$SOUS_DOM6],
     key.pos=NULL, main=NULL)
plot(st_geometry(xy), cex = .2, pch = 19, col = alpha("black",.3), add = T)
legend(0, 4e+05, legend = rev(levels(ecoregion$SOUS_DOM6)), 
       fill = rev(alpha(brewer.pal(6,"Spectral"),.8)),
       bty="n",border = "grey35", cex= .8)
dev.off()


### SCHEMA TRANSITION VS TIME ####



png(paste0("images/multinom_all_trans.png"), res=300, 
    width = 4, height = 1.8, units = 'in', bg = "transparent", type='cairo')
# quartz(width = 4, height = 1.7)
par(mar=c(0,0,0,0))
openplotmat(main = "")

pos <- coordinates(pos = 4, my = -0.1)
labs <- c(expression(bold('État'[0])), 
          expression(bold('État'[1])), 
          expression(bold('État'["..."])), 
          expression(bold('État'["n"])))

for (i in 1:4) {
  textellipse(mid=pos[i,], radx = 0.11, lab = labs[i], cex = 1.1, 
              col = '#404940', shadow.col='#98B283', shadow.size = 0.005)
  }
for (i in 1:3) { 
  curvedarrow(from = pos[1, ] + c(0, .1), to = pos[i+1, ] + c(0, .12),
              curve = -0.64, arr.type = "triangle", arr.length=.25, arr.width=.2,
              arr.pos = 1, lwd = 2, lcol ="grey40", arr.col ="grey40")}

straightarrow(from = pos[1, ] - c(0, .22), to = pos[4, ] - c(0, .22), 
              arr.type = "triangle", arr.pos = 1, lwd = 10, lcol ="grey60", arr.col ="grey60",
              arr.length=.5, arr.width=.4) 
textplain(pos[4, ] - c(0, .35), lab = "Time")
dev.off()

# quartz(width = 4, height = 2)
png(paste0("images/multinom_trans.png"), res=300, 
    width = 4, height = 2, units = 'in', bg = "transparent", type='cairo')
par(mar=c(0,0,0,0))
openplotmat(main = "")

pos <- coordinates(pos = 4, my = -0.2)
labs <- c(expression(bold('État'[0])), 
          expression(bold('État'[1])), 
          expression(bold('État'["..."])), 
          expression(bold('État'["n"])))

for (i in 1:4) 
  textellipse(mid=pos[i,], radx = 0.1, lab = labs[i], cex = 1,  
              col = '#404940', shadow.col='#98B283', shadow.size = 0.005)

curvedarrow(from = pos[1, ] + c(0,0.1), to = pos[4, ] + c(0,0.12),
            curve = -0.55, arr.type = "triangle", arr.length=.25, arr.width=.2,
            arr.pos = 1, lwd = 2, lcol ="grey40", arr.col ="grey40")

dev.off()

## PCA

# Compute Hellinger transformation
# sp.hel <- decostand(sp.mat[,MySpecies], "hel")
# 
# # PCA on transformed species matrix
# pca <- rda(sp.hel)
# quartz()
# plot(pca, display = "sites", type = 'n', scaling = 1)
# points(pca, scaling = 1, col = alpha(col.gr, 0.5), pch = pch.gr, cex = 0.8)
# legend("bottomright", legend = gr.name, col = levels(col.gr), pch = c(17,15,16,18), bty = "n")
# 


### INDVAL ####

indic.sp <- indval(sp.mat[, MySpecies], gr)

# Table of the significant indicator species
pAdj <- p.adjust(indic.sp$pval, "holm")
maxiv <- indic.sp$maxcls #[pAdj <= 0.05]
iv <- indic.sp$indcls #[pAdj <= 0.05]
pv <- indic.sp$pval #[pAdj <= 0.05]
fr <- apply(sp.mat[,MySpecies] > 0, 2, sum)#[pAdj <= 0.05]
fidg <- data.frame(group=maxiv, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(-fidg$group, -fidg$indval),]
fidg$group <- factor(fidg$group, 1:length(gr.name), gr.name)


# Heatmap

indic <- indic.sp$indval %>% 
  tibble::rownames_to_column(var = "Species")%>% 
  arrange(Species) %>% 
  melt(id.vars = "Species") 

indic$Species <- factor(x = indic$Species, levels = rev(row.names(fidg)))
indic$variable <- factor(x = indic$variable, levels = rev(unique(indic$variable)))

ggplot(data = indic, aes(x = variable, y = Species)) +
  geom_tile(aes(fill = value),colour="white") +
  #low= "#CED5F0" , high = "#072175"; low="#E5E5E5", high = "steelblue" #104E8B
  scale_fill_gradient(low="#E5E5E5", high = "red3", limits=c(0.05,0.9), na.value ="white") +
  theme_classic(base_size=9)+
  theme(axis.text.x= element_text(size = 8.5)) +
  labs(fill='IndVal', x= "Community states") 



### TRANSITION DIAGRAM ####


Trans.map1 <- slice(group_by(Trans_df, plot_id), 1)
Trans.map2 <- slice(group_by(Trans_df, plot_id), n())

Trans.map <- Trans.map1 %>%
  ungroup() %>%
  dplyr::select(plot_id, year1, From) %>%
  full_join(dplyr::select(Trans.map2, plot_id, year2, To), by = "plot_id") %>%
  mutate(lineNo = row_number()) %>%
  mutate(col.tr = factor(To, levels = gr.name, 
                         labels = pal.gr)) %>%
  mutate(From = factor(From, levels = gr.name, ordered = T), 
         To = factor(To, levels = gr.name, ordered = T)) 


M <- matrix(nrow = 4, ncol = 4, byrow = TRUE, data = table(Trans.map$From, Trans.map$To))
M.perc <- t(round(t(M)/colSums(M)*100, 0))


# quartz(width=8, height=6)
png(paste0("images/diag_trans_pourcen.png"), res=300, 
    width = 8, height = 6, units = 'in', bg = "transparent", type='cairo')

par(mar=c(0,0,0,0))
pos.box <- cbind (c(0.5, 0.2, 0.8, 0.5), c(0.8, 0.5, 0.5, 0.2))
plotmat(M.perc, pos = pos.box, curve = 0.07, name = gr.name, lwd = 2, relsize=.98,
        box.cex = 1.2, cex.txt = 1.2, txt.col = "white", dtext = 0.15, txt.font = 2,
        box.lwd = 0.1, box.type = "rect", shadow.size = 0.005,
        box.prop = 0.4, box.size = 0.1, box.col = pal.gr,
        arr.length=.2, arr.width=.15,  arr.type ="curved",
        arr.col ="grey40", arr.lcol ="grey40",
        self.cex = 0.6, self.shifty = c(.07,0,0,-.07), self.shiftx = c(0,-.14,.14,0), self.lwd = 2)

dev.off()



# quartz(width=8, height=6)
png(paste0("images/diag_trans_mixed.png"), res=300, 
    width = 8, height = 6, units = 'in', bg = "transparent", type='cairo')

par(mar=c(0,0,0,0))
pos.box <- cbind (c(0.5, 0.2, 0.8, 0.5), 
                  c(0.8, 0.5, 0.5, 0.2))
lwd.mat <- matrix(c(2,2,2,2,
                   4,4,4,4, 
                   2,2,2,2,
                   2,2,2,2), 4, 4)
col.mat <- matrix(c("grey80","grey80","grey80","grey80",
                    "grey40","grey40","grey40","grey40",
                    "grey80","grey80","grey80","grey80",
                    "grey80","grey80","grey80","grey80"), 4, 4) 

plotmat(M.perc, pos = pos.box, curve = 0.07, name = gr.name, lwd = lwd.mat, relsize=.98,
        box.cex = 1.2, cex.txt = 0, txt.col = "white", dtext = 0.15, txt.font = 2,
        box.lwd = 0.1, box.type = "rect", shadow.size = 0.005,
        box.prop = 0.4, box.size = 0.1, box.col = pal.gr,
        arr.length=.2, arr.width=.15,  arr.type ="curved",
        arr.col =col.mat, arr.lcol = col.mat,
        self.cex = 0.6, self.shifty = c(.07,0,0,-.07), self.shiftx = c(0,-.14,.14,0), 
        self.lwd = c(2,4,2,2))

dev.off()

### MAP OF TRANSITION ####

# figures organized in rectangle
quartz(width = 8.5, height = 7)
nf <- matrix(c(1,2,5,5,3,4), 3,2, byrow = T)
layout(nf,  heights= c(1, 0.7, 1))
par(mar=c(0,0,0,0), oma = c(0,0,0,0))
for(i in 1:4) {
  from_to = Trans.map %>% filter(From == gr.name[i])
  
  pch.gr = rep(19, nrow(from_to))
  pch.gr[from_to$To == gr.name[i]] = 4
  plot(ecoregion$geom, border = "grey60", lwd=0.7)
  plot(xy$geom[from_to$lineNo], pch = pch.gr, 
       col = alpha(from_to$col.tr, 0.4), cex = 0.6, add = T)
  text(-76, 53, labels = paste("From", gr.name[i]), cex = 1.5)
}

par(mar=c(0,13,0,13))
curve.mat = matrix(c(.07,.09,.05,.07, 
                     .09,.07,.07,.05, 
                     .05,.07,.07,.09,
                     .07,.05,.09,.07), 4, 4)

plotmat(M, pos = c(2, 2), curve = curve.mat, name = gr.name, lwd = 1, relsize=0.9,
        box.cex = 1.2, cex.txt = 0.9, 
        box.lwd = 0.1, box.type = "rect", shadow.size = 0,
        box.prop = 0.4, box.size = 0.1, box.col = pal.gr,
        arr.length=.1, arr.width=.15,  arr.type ="triangle", 
        arr.col ="grey40", arr.lcol ="grey40",
        self.cex = 0.5, self.shifty = c(.06,.06,-.06,-.06), self.shiftx = c(-0.1, 0.1, -0.1, 0.1))


### TRANSITION FLOW ####

pos <- coordinates(pos = 4, my = -.4)
pos1 <- c(.5,.9)
curv1 <- c(.4,.48,.52,.6)
curv2 <- c(.2,.4,.6,.8)

for(st in 1:4) {
  png(paste0("images/flow_trans_", gr.name[st], ".png"), res=300,
      width = 5, height = 2.8, units = 'in', bg = "transparent", type='cairo')
  par(mar=c(0,0,0,0))
  openplotmat(main = "")
  for(i in 1:4) {
    xspline(x = c(pos1[1], curv1[i], curv2[i], pos[i,1]), 
            y = c(pos1[2]-.02, 0.5, 0.3, pos[i,2]), 
            s = 1, lwd = M.perc[i,st]/4, border = "grey40")
  }
  
  
  textrect(mid=pos1, radx = 0.11, rady = 0.07, lab = gr.name[st], cex = 1, font = 2, 
           col = "white", box.col = pal.gr[st], lcol = pal.gr[st], shadow.size = 0.005)
  for (i in 1:4) {
    textrect(mid=pos[i,], radx = 0.11, rady = 0.07, lab = gr.name[i], cex = 1, font = 2, 
             col = "white", box.col = pal.gr[i], lcol = pal.gr[i], shadow.size = 0.005)

    if(i %in% 1:2) l = -.05
    if(i %in% 3:4) l = .05
    textplain(mid = pos[i,]+c(l,.1), lab = paste0(M.perc[i,st], "%"), 
              cex=.9, font=2, col = "grey25")
  }
  dev.off()
}





