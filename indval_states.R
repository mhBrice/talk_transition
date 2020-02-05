#### Indicator species of forest states ####

library(labdsv)

### DATA ####

states_ba <- readRDS("data/states_envba.RDS") 

sp_ba <- readRDS("data/sp_mat_ba_nov2019.RDS") %>% 
  filter(ID_PE_MES %in% states_ba$ID_PE_MES)


tree_code <- read.csv2("data/ref_spCode.csv")

# Species selection

sp_ba <- sp_ba %>%
  mutate(SORSP = SORAME + SORDEC) %>%
  select(-MALSP, -SORAME, -SORDEC, -ALNCRI, -CRASP, -ALNRUG, -PRUVIR)

MySpecies <- sort(names(which(colSums(sp_ba[,-c(1:4)])>0)))

sp_ba <- sp_ba %>%
  select(ID_PE, ID_PE_MES, plot_id, year_measured, MySpecies)

### INDVAL ####

indic_demo <- indval(as.data.frame(sp_ba[,-c(1:4)]), states_ba$states_ba)

tab <- t(indic_demo$indval)

# Remove species that barely contribute to any state

tab <- tab[,which(colSums(tab)>.04)]

# Order species by group

tab <- tab[, c(which(colnames(tab) %in% boreal),
              which(colnames(tab) %in% pioneer),
              which(colnames(tab) %in% temperate))]

### HEATMAP ####

# Color palette

pal0 <- c("#fcf1c5", "#9eb625", "#045579", "#0b2d40")
pal <- colorRampPalette(pal0)(30)

quantile.range <- quantile(tab, probs = seq(0, 1, 0.05), na.rm = T)

# Parameters

nsp = ncol(tab)
seqx = seq(0, 1, len = nrow(tab))
seqy = seq(0, 1, len = ncol(tab))

sp_name <- tree_code$complete.name[match(colnames(tab),tree_code$spCode)]

png("images/indval_states.png", width = 6, height = 6, units = "in", res = 300)
#quartz(width = 6, height = 6)
par(mar=c(2.5, 8, .2, .2))
image(tab[,nsp:1], col = (pal), axes = F)

text(seqx, rep(-.05, 4), c("Boréal", "Mixte", "Pionnier", "Tempéré"), 
     xpd = NA, cex = .9, font = 2)
mtext("États", 1, line = 1.3, font = 2, cex = 1)

text(rep(-.2, 28), rev(seqy), sp_name, xpd = NA, 
     cex = .8, adj = 1, font = 3)

arrows(y0 = mean(seqy[18:19]), x0 = -0.6, x1 = 1.2, length = 0, xpd = NA)
arrows(y0 = mean(seqy[14:15]), x0 = -0.6, x1 = 1.2, length = 0, xpd = NA)

dev.off()