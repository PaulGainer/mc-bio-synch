library(reshape2)
library(ggplot2)

scale <- scale_color_hue()

intervals <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)


show_energy <- TRUE
show_sync_times <- TRUE
show_energy_vs_times <- TRUE

## height and width of pdf output
gheight <- 10
gwidth<-20
dotsize <- 3

# directory names
pdfdir <- "pdf/"
energydir <- "energy"
fullenergydir <- paste(pdfdir,energydir, sep="")
timedir <- "time"
fulltimedir <- paste(pdfdir,timedir, sep="")
vsdir <- "energy_vs_time"
fulle_vs_tdir <- paste(pdfdir,vsdir, sep="")

if (!dir.exists(pdfdir)) {
  dir.create(pdfdir)
}
if (!dir.exists(fullenergydir)) {
  dir.create(fullenergydir)
}
if (!dir.exists(fulltimedir)) {
  dir.create(fulltimedir)
}
if (!dir.exists(fulle_vs_tdir)) {
  dir.create(fulle_vs_tdir)
}


raw_sync_energy = read.csv("results-mirollo-strogatz_micaz_energy_n8.csv", 
                           col.names=c("N","T","eps","RP", "ML", "U", "CD",
                                       "avg1", "avg2", "avg3", "avg4", "avg5", "avg6", "avg7", "avg8", "avg9", "avg10", 
                                       "max1", "max2", "max3", "max4", "max5", "max6", "max7", "max8", "max9", "max10" 
                           ), header=FALSE)

sync_energy <- melt(raw_sync_energy, id = c("N","T","eps","RP", "ML", "U", "CD"))

avg_energy <- sync_energy[grep("avg", sync_energy$variable),]
max_energy <- sync_energy[grep("max", sync_energy$variable),]

avg_energy$interval<- (as.numeric(factor(avg_energy$variable)) ) / 10
max_energy$interval<- (as.numeric(factor(max_energy$variable)) ) / 10


ranges_energy <- avg_energy
ranges_energy$avg <- avg_energy$value *1000/ranges_energy$N
ranges_energy$max <- max_energy$value *1000/ranges_energy$N



ranges_energy$variable <- NULL
ranges_energy <- ranges_energy[c("N", "T", "eps", "RP", "ML", "U", "CD","interval", "avg","max")]

sync_energy_time = read.csv("results-mirollo-strogatz_energy_time_n8.csv", 
                            col.names=c("N","T","eps","RP", "ML", "U", "CD",
                                        "avg_energy", "max_energy", "avg_time", "max_time"
                            ), header=FALSE)
sync_energy_time$avg_energy <- sync_energy_time$avg_energy * 1000/sync_energy_time$N
sync_energy_time$max_energy <- sync_energy_time$max_energy * 1000/sync_energy_time$N

if (show_energy) {
  for (e in c(0.1)) {
    for (ml in c(0.2)) {
      ranges_data_set=ranges_energy[ ranges_energy$eps ==e & ranges_energy$ML ==ml & ranges_energy$RP < 5,]
      
      pow_vs_ord <- ggplot(ranges_data_set) + 
        geom_line(aes( x=interval,y=avg, colour=factor(RP), linetype="avg"))+
        geom_line(aes( x=interval,y=max, colour=factor(RP), linetype="max"))+
        geom_point(aes( x=interval,y=avg, colour=factor(RP), shape=factor(RP)),size = dotsize)+
        geom_point(aes( x=interval,y=max, colour=factor(RP), shape=factor(RP)),size = dotsize)+
        
        expand_limits(y=0,x=0.1)+
        coord_cartesian(ylim= c(0,.6))+
        scale_x_continuous(breaks=intervals, expand = c(0.01, 0)) + 
        scale+
        ggtitle(substitute("Power Consumption per Node (mWh) (N=8 T=10"~epsilon*"="*e*""~mu*"="*ml*")", list(e=bquote(.(e)), ml=bquote(.(ml))))) +
        labs( x = "Phase Coherence", y = "Expected Power Cons. (mWh)", colour="R", shape="R", linetype="") 
      if (e == 0.1 & ml == 0.2) {
        pow_vs_ord <- pow_vs_ord +  theme_bw()+ theme(legend.position = c(.1,.75), legend.box="horizontal") 
      } else {
        pow_vs_ord <- pow_vs_ord +  theme_bw()
      }
      print(pow_vs_ord)
      ggsave(sprintf("%s/ms_n8_t10_e%f_ml%f_pow_vs_ord.pdf",fullenergydir,e,ml), plot=pow_vs_ord, device="pdf", width=gwidth, height=gheight, units = "cm")
    }
  }
  for (e in c(0.1)) {
    full_sync_energy <- sync_energy_time[ sync_energy_time$eps == e & sync_energy_time$RP < 5, ]
    pow_vs_ml <- ggplot(full_sync_energy) + 
      geom_line(aes( x=ML,y=avg_energy, colour=factor(RP), linetype="avg"))+
      geom_line(aes( x=ML,y=max_energy, colour=factor(RP), linetype="max"))+
      geom_point(aes( x=ML,y=avg_energy, colour=factor(RP), shape=factor(RP)),size = dotsize)+
      geom_point(aes( x=ML,y=max_energy, colour=factor(RP), shape=factor(RP)),size = dotsize)+
      
      expand_limits(y=0,x=0.1)+
      coord_cartesian(ylim = c(0,1.5))+
      scale_x_continuous(breaks=intervals, expand = c(0.01, 0))+
      scale+
      ggtitle(substitute("Power Consumption per Node (mWh) (N=8 T=10"~epsilon*"="*e*")", list(e=bquote(.(e))))) +
      labs( x = "Broadcast Failure Probability", y = "Expected Energy Cons. (mWh)", colour="R", shape="R", linetype="")  +
      theme_bw()+
      theme(legend.position = c(.1,.75), legend.box="horizontal") 
    print(pow_vs_ml)
    ggsave(sprintf("%s/ms_n8_t10_e%f_pow_vs_ml.pdf",fullenergydir,e,ml), plot=pow_vs_ml, device="pdf", width=gwidth, height=gheight, units = "cm")      
  }
}

raw_sync_times = read.csv("results-mirollo-strogatz_sync_times_n8.csv", 
                          col.names=c("N","T","eps","RP", "ML", "U", "CD",
                                      "avg1", "avg2", "avg3", "avg4", "avg5", "avg6", "avg7", "avg8", "avg9", "avg10", 
                                      "max1", "max2", "max3", "max4", "max5", "max6", "max7", "max8", "max9", "max10" 
                          ), header=FALSE)

sync_times <- melt(raw_sync_times, id = c("N","T","eps","RP", "ML", "U", "CD"))


avg_times <- sync_times[grep("avg", sync_times$variable),]
max_times <- sync_times[grep("max", sync_times$variable),]

avg_times$interval<- (as.numeric(factor(avg_times$variable)) ) / 10
max_times$interval<- (as.numeric(factor(max_times$variable)) ) / 10

ranges_times <- avg_times
ranges_times$avg <- avg_times$value 
ranges_times$max <- max_times$value 
ranges_times$variable <- NULL
ranges_times <- ranges_times[c("N", "T", "eps", "RP", "ML", "U", "CD","interval", "avg","max")]

if (show_sync_times) {
  for (e in c(0.1)) {
    for (ml in c(0.2)) {
      ranges_data_set=ranges_times[ ranges_times$eps ==e & ranges_times$ML ==ml & ranges_times$RP < 5,]
      pow_vs_ord <- ggplot(ranges_data_set) + 
        geom_line(aes( x=interval,y=avg, colour=factor(RP), linetype="avg"))+
        geom_line(aes( x=interval,y=max, colour=factor(RP), linetype="max"))+
        geom_point(aes( x=interval,y=avg, colour=factor(RP), shape=factor(RP)),size = dotsize)+
        geom_point(aes( x=interval,y=max, colour=factor(RP), shape=factor(RP)),size = dotsize)+
        coord_cartesian(ylim = c(0,6))+
        expand_limits(y=0,x=0.1)+
        scale_x_continuous(breaks=intervals, expand = c(0.01, 0))+
        scale+
        ggtitle(substitute("Time for Synchronisation (cycles) (N=8 T=10"~epsilon*"="*e*""~mu*"="*ml*")", list(e=bquote(.(e)), ml=bquote(.(ml))))) +
        labs( x = "Phase Coherence", y = "Expected Sync Time (cycles)", colour="R", shape="R", linetype="") 
      if (e == 0.1 & ml == 0.2) {
        pow_vs_ord <- pow_vs_ord +   theme_bw()+ theme(legend.position = c(.1,.75), legend.box="horizontal") 
      } else {
        pow_vs_ord <- pow_vs_ord +   theme_bw()
      }
      print(pow_vs_ord)
      ggsave(sprintf("%s/ms_n8_t10_e%f_ml%f_sync_vs_ord.pdf",fulltimedir,e,ml), plot=pow_vs_ord, device="pdf", width=gwidth, height=gheight, units = "cm")
    }
  }
}


if (show_energy_vs_times) {
  vs_plot_avg <- ggplot(sync_energy_time) +
    geom_point(aes(x=avg_time, y=avg_energy, color=factor(RP), shape=factor(RP)),size = dotsize) +
    scale +
    geom_smooth(method ="lm", aes(x=avg_time, y=avg_energy, color=factor(RP)), size=.1)+
    ggtitle(substitute(" Relation between Power Consumption and Time (N=8 T=10) [average]")) +
    labs( x = "Sync. Time (cycles)", y = "Power Consumption (mWh)", colour="R", shape="R",  linetype="") +
    expand_limits(y=0,x=0)+
    scale_x_continuous( expand = c(0.01, 0))+
    theme_bw() + theme(legend.position = c(.05,.75), legend.box="horizontal") 
  print(vs_plot_avg)
  vs_plot_max <- ggplot(sync_energy_time) +
    geom_point(aes(x=max_time, y=max_energy, color=factor(RP), shape=factor(RP)) ,size = dotsize) +
    geom_smooth(method ="lm", aes(x=max_time, y=max_energy, color=factor(RP)), size=.1)+
    scale +
    ggtitle(substitute(" Relation between Power Consumption and Time (N=8 T=10) [maximal]")) +
    labs( x = "Sync. Time (cycles)", y = "Power Consumption (mWh)", colour="R", shape="R",  linetype="") +
    expand_limits(y=0,x=0)+
    coord_cartesian(xlim = c(0,10), ylim = c(0,1.3))+
    scale_x_continuous( expand = c(0.01, 0))+
    theme_bw() + theme(legend.position = c(.05,.75), legend.box="horizontal") ;
  print(vs_plot_max)
  
  ggsave(sprintf("%s/ms_n8_t10_energy_vs_time_avg.pdf",fulle_vs_tdir), plot=vs_plot_avg, device="pdf", width=gwidth, height=gheight, units = "cm")
  ggsave(sprintf("%s/ms_n8_t10_energy_vs_time_max.pdf",fulle_vs_tdir), plot=vs_plot_max, device="pdf", width=gwidth, height=gheight, units = "cm")
}