library(reshape2)
library(ggplot2)

## height and width of pdf output
gheight <- 10
gwidth<-20
dotsize <- 3
# directory names
pdfdir <- "pdf/"
graphdir <- "resync"
fulldir <- paste(pdfdir,graphdir, sep="")

if (!dir.exists(pdfdir)) {
  dir.create(pdfdir)
}
if (!dir.exists(fulldir)) {
  dir.create(fulldir)
}
scale <- scale_color_hue()

resync_energy_n15_to_n30 = read.csv("results-mirollo-strogatz_resync_energy.csv", 
                                    col.names=c("N","T","eps","RP", "ML", "U", "CD","avg" 
                                    ), header=FALSE)

resync_energy_n10 = read.csv("results-mirollo-strogatz_resync_energy_n10.csv", 
                             col.names=c("N","T","eps","RP", "ML", "U", "CD","avg" 
                             ), header=FALSE)

resync_energy_n35 = read.csv("results-mirollo-strogatz_resync_energy_n35.csv", 
                             col.names=c("N","T","eps","RP", "ML", "U", "CD","avg" 
                             ), header=FALSE)


resync_energy <- rbind(resync_energy_n15_to_n30,resync_energy_n10)
resync_energy <- rbind(resync_energy,resync_energy_n35)
resync_energy$avg_mwh <- resync_energy$avg *1000 / resync_energy$N

for (e in c(0.1 )) {
  for (ml in c(0.2 )) { 
    ranges_data_set=resync_energy[ resync_energy$eps ==e & resync_energy$ML ==ml & resync_energy$RP < 5,]
    pow_vs_ord <- ggplot(ranges_data_set) + 
      geom_line(aes( x=N,y=avg_mwh, colour=factor(RP), linetype=factor(U)))+
      geom_point(aes(x=N, y=avg_mwh, color=factor(RP), shape=factor(RP)),size = dotsize) +
      scale_x_continuous(expand = c(0.01, 0)) + 
      scale+
      ggtitle(substitute("Power Consumption per Node for Restabilisation  (mWh) (T=10"~epsilon*"="*e*" "~mu*"="*ml*")", list(e=bquote(.(e)), ml=bquote(.(ml))))) +
      labs( x = "N", y = "Expected Power Consumption Cons. (mWh)", colour="R", linetype="U", shape="R") 
    if (e == 0.1 & ml == 0.2) {
      pow_vs_ord <- pow_vs_ord  + theme_bw()+ theme(legend.position = c(.90,.75), legend.box="horizontal") 
      
    } else {
      pow_vs_ord <- pow_vs_ord + theme_bw() 
    }
    print(pow_vs_ord)
    ggsave(sprintf("%s/ms_t10_e%f_ml%f_pow_vs_N.pdf",fulldir,e,ml), plot=pow_vs_ord, device="pdf", width=gwidth, height=gheight, units = "cm")
  }
}