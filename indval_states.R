#### Indicator species of forest states ####

library(labdsv)
library(dplyr)
library(diagram)
library(graphicsutils)

### DATA ####

states_ba <- readRDS("data/states_envba.RDS") 

sp_ba <- readRDS("data/sp_mat_ba_nov2019.RDS") %>% 
  filter(ID_PE_MES %in% states_ba$ID_PE_MES)


tree_code <- read.csv2("data/ref_spCode.csv")

### FUNCTIONS ####

addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE, # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
  xpd = NA
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate, xpd = NA)
}

# Species selection

sp_ba <- sp_ba %>%
  mutate(SORSP = SORAME + SORDEC) %>%
  select(-MALSP, -SORAME, -SORDEC, -ALNCRI, -CRASP, -ALNRUG, -PRUVIR)

MySpecies <- sort(names(which(colSums(sp_ba[,-c(1:4)])>0)))

sp_ba <- sp_ba %>%
  select(ID_PE, ID_PE_MES, plot_id, year_measured, MySpecies)

### Species group ####
states <- c("Boréal", "Mixte", "Pionnier", "Tempéré")

pioneer <- c("BETPAP", "BETPOP", 
             "POPBAL", "POPDEL", "POPGRA", "POPTRE", 
             "PRUPEN", "SALSP", "SORSP")

temperate <- c("ACENEG", "ACENIG", "ACEPEN", "ACERIN", "ACERUB", "ACESAC", "ACESPI", 
               "AMESP", "BETALL",
               "CARCAR", "CARCOR", "FAGGRA", 
               "FRAAME", "FRANIG", "FRAPEN", 
               "JUGCIN", "OSTVIR",
               "PICRUB", "PINRES", "PINRIG", "PINSTR",
               "PRUSER", 
               "QUEALB", "QUEBIC", "QUEMAC", "QUERUB", 
               "THUOCC", "TILAME", "TSUCAN", 
               "ULMAME", "ULMRUB", "ULMTHO")

boreal <- c("ABIBAL", "LARLAR", "PICGLA", "PICMAR", "PINBAN")


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

# image

acer <- readPNG("images/acer.png")
betula <- readPNG("images/betula.png")
abies <- readPNG("images/sapin.png")

# Color palette

pal0 <- c("#fcf1c5", "#9eb625", "#045579", "#0b2d40")
pal <- colorRampPalette(pal0)(30)
st_col <- c("#158282", "#A1BD93", "#FEAC19", "#D43650")

quantile.range <- quantile(tab, probs = seq(0, 1, 0.05), na.rm = T)

# Parameters

nsp = ncol(tab)
seqx = seq(0, 1, len = nrow(tab))
seqy = seq(0, 1, len = ncol(tab))

sp_name <- tree_code$complete.name[match(colnames(tab),tree_code$spCode)]



png("images/indval_states.png", width = 7, height = 6, units = "in", res = 300)
#quartz(width = 7, height = 6)
layout(matrix(1:2), heights = c(1,.13))
par(mar=c(0, 12, .2, .2))
image(tab[,nsp:1], col = (pal), axes = F)

# text(seqx, rep(-.05, 4), states, 
#      xpd = NA, cex = 1, font = 2, col = st_col)
# mtext("États", 1, line = 1.3, font = 2, cex = 1.2)

text(rep(-.2, nsp), rev(seqy), sp_name, xpd = NA, 
     cex = .9, adj = 1, font = 3)

arrows(y0 = mean(seqy[18:19]), x0 = -0.6, x1 = 1.2, length = 0, xpd = NA)
arrows(y0 = mean(seqy[14:15]), x0 = -0.6, x1 = 1.2, length = 0, xpd = NA)

addImg(abies, x = -.73, y = .92, width = .16)
addImg(betula, x = -.73, y = .7, width = .15)
addImg(acer, x = -.73, y = .5, width = .25)

par(mar = c(0,12,0,.2))
plot0()

pos <- coordinates(pos = 4, my = .4)
for (i in 1:4) {
  textrect(mid=pos[i,], radx = 0.21, rady = 0.4, lab = states[i],
           cex = 1.1, font = 2,
           col = "white", box.col = st_col[i], lcol = st_col[i],
           shadow.size = 0.005, xpd = NA)
}
mtext("États", 1, line = -1.3, font = 2, cex = 1.2)

dev.off()




