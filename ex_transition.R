
st_tr <- list(c("Boréal", "Boréal", "Pionnier", "Boréal"),
              c("Mixte", "Mixte", "Mixte", "Tempéré", "Tempéré"),
              c("Tempéré", "Tempéré", "Tempéré"))

pos <- list(x1=matrix(c(-0.81, -0.27, 0.27, 0.81, rep(.6, 4)), 4),
            x2=matrix(c(-0.84, -0.42, 0, 0.42, 0.84, rep(0, 5)), 5),
            x3=matrix(c(-0.72, 0, 0.72, rep(-.5, 3)), 3))


trans_nb <- statetable.msm(states_num, plot_id, data = states_ba)
rownames(trans_nb) <- colnames(trans_nb) <- states
trans_perc <- (round(trans_nb/rowSums(trans_nb)*100, 1))

png("images/ex_transition.png", width = 9.5, height = 6, units = "in", res = 300)
#quartz(width = 9.5, height = 6)

layout(matrix(1:2,1), widths = c(1, .74))
par(mar=c(0,0.3,0,.3))
plot0()
for(tr in 1:3) {
  posi <- pos[[tr]]
  st_tri <- st_tr[[tr]]
  ntr <- length(st_tri)
  
  for (i in 1:ntr) {
    col = st_col[which(states == st_tri[i])]
    textrect(mid=posi[i,], radx = 0.17, rady = 0.06, lab = st_tri[i],
             cex = 1, font = 2,
             col = "white", box.col = col, lcol = col,
             shadow.size = 0.005, xpd = NA)
  }
  for (i in 1:(ntr-1)) {
    curvedarrow(from = posi[i, ] + c(.01, .08), to = posi[i+1, ] + c(-.01, .08),
                curve = -.07*ntr, arr.type = "triangle", arr.length=.3, arr.width=.4,
                arr.pos = .5, lwd = 3, lcol ="grey40", arr.col ="grey40")}
  
}

text(0, -.7, ". . .", cex = 2.5, font = 2, col ="grey40")

straightarrow(from = c(-.8,-.83), to = c(.8,-.83),
              arr.type = "triangle", arr.pos = 1, lwd = 15, lcol ="grey60", arr.col ="grey60",
              arr.length=.7, arr.width=.7)
textplain(c(.8,-.93), lab = "Temps", cex = 1.2)

CurlyBraces(x0=.95, x1=.95, y0=-.8, y1=.8, pos = 1, direction = 1, depth=.2)

par(mar = c(7, 4, 9, .5))
plot_trans(trans_perc, states_lab = states, labels = T)

dev.off()