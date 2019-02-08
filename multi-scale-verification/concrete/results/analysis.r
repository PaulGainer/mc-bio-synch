library(ggplot2)

n4 = read.csv("n4_sync.csv", 
              col.names=c("RP","ML", 
                          "prob_sync"),header=TRUE)

# size of the created pdfs
gheight <- 10
gwidth<-20

# font size in created graphs 
gtext <- 18

# constants to define which graphs and tables should be created
create_graphs <- TRUE
create_tables <- TRUE

# greyscale or colour?
grey_scale <- FALSE

if (grey_scale) {
  scale <- scale_color_grey(end=0.7)
  
} else {
  scale <- scale_color_hue()
}

prob_data = n4[n4$ML <1 , ]

# create graphs for visualising the model checking results
        prob_vs_rp <- ggplot(data=prob_data,aes(y = prob_sync, x = RP, color=factor(ML), shape=factor(ML)) ) +
          geom_line() +
          geom_point(size=3)+
          theme_bw()+
          theme(legend.position=c(.95,.8), legend.justification = c(0.95,0.85), text = element_text(size=gtext))+
          scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14))+
          scale+
          ggtitle(substitute("Concrete Model Synchronisation (N=4 T=10"~epsilon*"=.1)")) +
          labs( x = "R", y = "Sync. Probability", colour=expression(mu), shape=expression(mu)) 

  print(prob_vs_rp)
  ggsave("concrete_n4_t10_e01_prob_vs_rp.pdf", plot=prob_vs_rp, device="pdf", width=gwidth, height=gheight, units = "cm")