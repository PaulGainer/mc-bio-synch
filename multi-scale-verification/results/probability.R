library(ggplot2)
#library(xtable)

#
# Load the model checking results. The first 5 columns are the parameters of 
# the model. Then we have 3 sections of 7 columns each. Within these sections, we have 
# the following structure:
# model checking result, model construction time, states of the model, transistions of the model, size of probability matrix, time for model checking, memory for model checking
#
# The first section contains the result for checking P=?[F synchronised], the second for R{"time_to_synch"}=? [F synchronised]
# The third section contains results we did not use for the analysis.

n3 = read.csv("ms_n3.csv",
              col.names=c("N","T","eps","RP", "ML",
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem",
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem" ,
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem"
              ), header=FALSE)

n4 = read.csv("ms_n4.csv",
              col.names=c("N","T","eps","RP", "ML",
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem",
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem"
              ),header=FALSE)

n5 = read.csv("ms_n5.csv",
              col.names=c("N","T","eps","RP", "ML",
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem",
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem"
              ),header=FALSE)
n6 = read.csv("ms_n6.csv",
              col.names=c("N","T","eps","RP", "ML",
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem",
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem"
              ),header=FALSE)
n7 = read.csv("ms_n7.csv",
              col.names=c("N","T","eps","RP", "ML",
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem",
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem"
              ),header=FALSE)

n8_prob = read.csv("ms-prob_n8.csv", 
              col.names=c("N","T","eps","RP", "ML", "U", "CD",
                          "prob_sync"),header=FALSE)

n8_time = read.csv("ms-time_n8.csv", 
              col.names=c("N","T","eps","RP", "ML", "U", "CD",
                          "time_sync"),header=FALSE)

n8 = merge(n8_prob, n8_time, by=c("N","T","eps","RP", "ML", "U", "CD"), all=TRUE)

keep = c( "N","T","eps","RP", "ML",
          "prob_sync", "time_sync")
# create a single dataframe containing the results
 rawdata = rbind(n3, n4)
 rawdata = rbind(rawdata,n5)
 rawdata = rbind(rawdata,n6)
 rawdata = rbind(rawdata,n7)
 rawdata = rawdata[keep]
 n8 = n8[keep]
 rawdata = rbind(rawdata,n8)

# some cleaning of the data: remove missing values (due to failure of prism) 
# and censordata for T > 10 :) (we did not get results for N=7 and T=11 in reasonable time)
prob_data <- rawdata[!is.na(rawdata$prob_sync) & rawdata$T < 11,]
time_data <- rawdata[!is.na(rawdata$time_sync)& rawdata$T < 11,]

# just out of curiosity, we save the missing values in new dataframes
#prob_missingdata <- rawdata[is.na(rawdata$prob_sync),]
#time_missingdata <- rawdata[is.na(rawdata$time_sync),]

# our own colour scheme, if colour is used in the graphs (most palettes do not contain enough colours)
col <- c("blue","red", "green","violet","black","orange",
         "steelblue4", "brown","yellowgreen", "darkgreen", 
         "palevioletred3", "salmon1", "chocolate4", "firebrick2",
         "grey43")


# size of the created pdfs
gheight <- 10
gwidth<-20
dotsize <- 3

# font size in created graphs 
gtext <- 18

# constants to define which graphs and tables should be created
create_graphs <- TRUE
create_tables <- FALSE

fix_n_t_e_prob_vs_ml <- TRUE
fix_n_t_e_prob_vs_rp <- TRUE
#fix_t_r_e_time_vs_N <- TRUE
fix_n_r_ml_time_vs_T <-FALSE
fix_n_r_e_time_vs_T <- FALSE
fix_n_ml_e_rp_perc <- FALSE
fix_t_r_ml_time_vs_N <-FALSE
fix_t_r_e_time_vs_N <-TRUE
fix_r_ml_e_time_vs_N <- FALSE
fix_r_ml_e_time_vs_T <-FALSE
checking_graphs <- FALSE

# greyscale or colour?
grey_scale <- FALSE

if (grey_scale) {
  scale <- scale_color_grey(end=0.7)
  
} else {
  #scale <-  scale_colour_manual(values= col) 
  scale <- scale_color_hue()
  
  }

# create graphs for visualising the model checking results
if(fix_n_t_e_prob_vs_ml & create_graphs) {
  for (n in c(8)) {
    for (t in c(10)) {
      for (e in c(0.1)) {
#        prob_data_set=prob_data[prob_data$N == n & prob_data$T == t & prob_data$eps ==e,]
        time_data_set=time_data[time_data$N == n & time_data$T == t & time_data$eps ==e ,]
        
        # prob_vs_ml <- ggplot(data=prob_data_set,aes(y = prob_sync, x = ML, color=factor(RP), shape=factor(RP)) ) +
        #   geom_line() +
        #   geom_point(size=3)+
        #   theme_bw()+
        #   scale+
        #   scale_shape_manual(values=1:nlevels(factor(prob_data$RP))) +
        #   ggtitle(substitute("M&S Synchronisation Probability (N="*n*", T="*t*","~epsilon*"="*e*")", list(n = bquote(.(n)), t = bquote(.(t)), e=bquote(.(e))))) +
        #   labs( x = expression(mu), y = "Sync. Probability", colour="RP", shape="RP") 
        
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_ml <- ggplot(time_data_set)+ #,aes(y = time_sync, x = ML, color=factor(RP), shape=factor(RP), size=dotsize) )
        geom_line(aes( x=ML,y=time_sync, color=factor(RP)))+
          geom_point(aes(x=ML, y=time_sync, color=factor(RP), shape=factor(RP)), size=dotsize) +
          theme_bw()+
          theme(legend.position=c(.125,.84), legend.justification = c(0.95,0.85), text = element_text(size=gtext))+
          guides(colour= guide_legend(keywidth = 2))+
#          geom_line() +
#          geom_point(size=3)+
          expand_limits(y=0,x=0)+
          scale+
#          scale_shape_manual(values=1:nlevels(factor(prob_data$RP))) +
          ggtitle(substitute("Time for Synchronisation (N="*n*", T="*t*","~epsilon*"="*e*")", list(n = bquote(.(n)), t = bquote(.(t)), e=bquote(.(e))))) +
          labs( x = expression(mu), y = "Sync. Time (in Cycles)", colour="RP",shape="RP") 
#        print(prob_vs_ml)
        print(time_vs_ml)
#        ggsave(sprintf("pdf/perez_n%d_t%d_e%f_prob_vs_ml.pdf",n,t,e), plot=prob_vs_ml, device="pdf", width=gwidth, height=gheight, units = "cm")
        ggsave(sprintf("pdf/perez_n%d_t%d_e%f_time_vs_ml.pdf",n,t,e), plot=time_vs_ml, device="pdf", width=gwidth, height=gheight, units = "cm")
      }
    }
  }
  
}

if(fix_n_t_e_prob_vs_rp & create_graphs) {
  for (n in  c(8)) {
    for (t in c(10)) {
      for (e in c(0.1)) {
        prob_data_set=prob_data[prob_data$N == n & prob_data$T == t & prob_data$eps ==e & (prob_data$ML==0 | prob_data$ML==0.2 |prob_data$ML==0.4 |prob_data$ML==0.6 |prob_data$ML==0.8 ),]
#        time_data_set=time_data[time_data$N == n & time_data$T == t & time_data$eps ==e,]
        
        prob_vs_rp <- ggplot(data=prob_data_set, aes( y = prob_sync, x = RP, color=factor(ML),shape=factor(ML)) ) +
          theme_bw()+
          theme(legend.position=c(.95,.8), legend.justification = c(0.95,0.85), text = element_text(size=gtext))+
          guides(colour= guide_legend(keywidth = 2))+
          geom_line() +
          geom_point(size=3)+
          scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14))+
          scale_y_continuous()+
#          scale_shape_manual(values=1:nlevels(factor(prob_data$ML))) +
          scale+
          ggtitle(substitute("Synchronisation Probability (N="*n*", T="*t*","~epsilon*"="*e*")", list(n = bquote(.(n)), t = bquote(.(t)), e=bquote(.(e))))) +
          labs( x = "RP", y = "Sync. Probability", color=expression(mu), shape=expression(mu))
        
#        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
#        time_vs_rp <- ggplot(data=time_data_set ) +
#          geom_line(aes(y = time_sync, x = RP, color=factor(ML))) +
          # theme_bw()+
          # expand_limits(y=0)+
          # scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14))+
          # scale+
          # ggtitle(substitute("M&S Synchronisation Time (N="*n*", T="*t*","~epsilon*"="*e*")", list(n = bquote(.(n)), t = bquote(.(t)), e=bquote(.(e))))) +
          # labs( x = "RP", y = "Sync. Time (in Cycles)", color=expression(mu))
        
        print(prob_vs_rp)
#        print(time_vs_rp)
        ggsave(sprintf("pdf/perez_n%d_t%d_e%f_prob_vs_rp.pdf",n,t,e), plot=prob_vs_rp, device="pdf", width=gwidth, height=gheight, units = "cm")
#        ggsave(sprintf("pdf/perez_n%d_t%d_e%f_time_vs_rp.pdf",n,t,e), plot=time_vs_rp, device="pdf", width=gwidth, height=gheight, units = "cm")
        
        
      }
    }
  }
}


if (fix_n_r_ml_time_vs_T & create_graphs) {
  for (n in c(7)) {
    for (r in c(1)) {
      for (ml in c(0.5)) {
        time_data_set=time_data[time_data$N == n &time_data$RP==r & time_data$ML ==ml,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_T <- ggplot(data=time_data_set ) +
          geom_line(aes(y = time_sync, x = T, color=factor(eps))) +
          theme_bw()+
          expand_limits(y=0)+
          scale+
          
          scale_y_continuous(limits = c(0,NA))+
          ggtitle(substitute("M&S Synchronisation Time (N="*np*", RP="*r*","~mu*"="* ml*")", list(np = bquote(.(n)), r = bquote(.(r)), ml=bquote(.(ml))))) +
          labs( x = "T", y = "Sync. Time (in Cycles)", color=expression(epsilon))+
          theme_bw()
        
        print(time_vs_T)
        ggsave(sprintf("pdf/perez_n%d_rp%d_ml%f_time_vs_T.pdf",n,r,ml), plot=time_vs_T, device="pdf", width=gwidth, height=gheight, units = "cm")
        
      }
    }
  }
  
}

if (fix_n_r_e_time_vs_T & create_graphs) {
  for (n in c(7)) {
    for (r in c(1)) {
      for (e in c(0.2)) {
        
        time_data_set=time_data[time_data$N == n &time_data$RP==r & time_data$eps ==e,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_T <- ggplot(data=time_data_set ) +
          geom_line(aes(y = time_sync, x = T, color=factor(ML))) +
          theme_bw()+
          expand_limits(y=0)+
          scale+
          scale_y_continuous(limits = c(0,NA))+
          ggtitle(substitute("M&S Synchronisation Time (N="*np*", RP="*r*","~epsilon*"="* e*")", list(np = bquote(.(n)), r = bquote(.(r)), e=bquote(.(e))))) +
          labs( x = "T", y = "Sync. Time (in Cycles)", color=expression(mu))+
          theme_bw()
        print(time_vs_T)
        ggsave(sprintf("pdf/perez_n%d_rp%d_e%f_time_vs_T.pdf",n,r,e), plot=time_vs_T, device="pdf", width=gwidth, height=gheight, units = "cm")
        
      }
    }
  }
  
}

if (fix_n_ml_e_rp_perc & create_graphs) {
  for (n in c(7)) {
    for (ml in c(0.2)) {
      for (e in c(0.2)) {
        
        prob_data_set=prob_data[ prob_data$N == n & prob_data$eps == e & prob_data$ML ==ml,]
        time_data_set=time_data[time_data$N == n &time_data$eps==e & time_data$ML ==ml,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        prob_vs_rp_perc <- ggplot(data=prob_data_set, aes( y = prob_sync, x = (RP/ T * 100)))  +
          theme_bw()+
          geom_line(aes(color=factor(T))) +
          scale+
          ggtitle(substitute("M&S Synchronisation Probability (N="*np*", "*epsilon*"="*e*","~mu*"="* ml*")", list(np = bquote(.(n)), e = bquote(.(e)), ml=bquote(.(ml))))) +
          labs( x = "RP (percentage of T)", y = "Sync. Probability", color="T")
        
        time_vs_rp_perc <- ggplot(data= time_data_set ) +
          geom_line(aes(y = time_sync, x = ( RP/  T * 100), color=factor( T))) +
          expand_limits(y=0)+
          theme_bw()+
          scale+
          ggtitle(substitute("M&S Synchronisation Time (N="*np*", "*epsilon*"="*e*","~mu*"="* ml*")", list(np = bquote(.(n)), e = bquote(.(e)), ml=bquote(.(ml))))) +
          labs( x = "RP (percentage of T)", y = "Sync. Time (in Cycles)", color="T")
        
        
        print(prob_vs_rp_perc)
        print(time_vs_rp_perc)
        ggsave(sprintf("pdf/perez_n%d_ml%f_e%f_prob_vs_rp_perc.pdf",n,ml,e), plot=prob_vs_rp_perc, device="pdf", width=gwidth, height=gheight, units = "cm")
        ggsave(sprintf("pdf/perez_n%d_ml%f_e%f_time_vs_rp_perc.pdf",n,ml,e), plot=time_vs_rp_perc, device="pdf", width=gwidth, height=gheight, units = "cm")
      }
    }
  }
}

if(fix_t_r_ml_time_vs_N & create_graphs) {
  for (t in c(10)) {
    for (r in c(1,4)) {
      for (ml in c(0.2, 0.3)) {
        time_data_set=time_data[time_data$T == t &time_data$RP==r & time_data$ML ==ml,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_N <- ggplot(data=time_data_set ) +
          theme_bw()+
          expand_limits(y=0)+
          geom_line(aes(y = time_sync, x = N, color=factor(eps))) +
          scale+
          ggtitle(substitute("M&S Synchronisation Time (T="*t*", RP="*r*","~mu*"="* ml*")", list(t = bquote(.(t)),r = bquote(.(r)), ml=bquote(.(ml))))) +
          labs( x = "N", y = "Sync. Time (in Cycles)", color=expression(epsilon))
        print(time_vs_N)
        
        ggsave(sprintf("pdf/perez_t%d_rp%d_ml%f_time_vs_N.pdf",t,r,ml), plot=time_vs_N, device="pdf", width=gwidth, height=gheight, units = "cm")
      }
    }
  }
}

if(fix_t_r_e_time_vs_N & create_graphs) {
  
  for (t in c(10)) {
    for (r in c(1)) {
      for (e in c(0.1)) {
        time_data_set=time_data[time_data$T == t &time_data$RP==r & time_data$eps ==e& (time_data$ML==0.1 | time_data$ML==0.2 |time_data$ML==0.4 |time_data$ML==0.6 |time_data$ML==0.8 ) ,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        # we use a special treatment for refractory period of 1 and epsilon=.1, to create nicer graph in the paper
        if (r ==1 & e==0.1 ) {
          time_vs_N <- ggplot(data=time_data_set,aes(y = time_sync, x = N, color=factor(ML), shape=factor(ML)) ) +
            theme_bw()+
            expand_limits(y=c(0,15))+
            geom_line() +
            geom_point(size=3)+
            scale+
            theme(legend.position=c(.95,.16), legend.justification = c(0.95,0.85),legend.direction = "horizontal", text = element_text(size=gtext))+
            guides(colour= guide_legend(keywidth = 2))+
            scale_y_continuous()+
#            scale_shape_manual(values=1:nlevels(factor(prob_data$ML))) +
            ggtitle(substitute("Time for Synchronisation (T="*t*", RP="*r*","~epsilon*"="* e*")", list(t = bquote(.(t)),r = bquote(.(r)), e=bquote(.(e))))) +
            labs( x = "N", y = "Sync. Time (in Cycles)", color=expression(mu), shape=expression(mu))
        } else {
          time_vs_N <- ggplot(data=time_data_set ) +
            theme_bw()+
            expand_limits(y=0)+
            geom_line(aes(y = time_sync, x = N, color=factor(ML))) +
            scale+
            ggtitle(substitute("M&S Synchronisation Time (T="*t*", RP="*r*","~epsilon*"="* e*")", list(t = bquote(.(t)),r = bquote(.(r)), e=bquote(.(e))))) +
            labs( x = "N", y = "Sync. Time (in Cycles)", color=expression(mu))
        }
        print(time_vs_N)
        
        ggsave(sprintf("pdf/perez_t%d_rp%d_e%f_time_vs_N.pdf",t,r,e), plot=time_vs_N, device="pdf", width=gwidth, height=gheight, units = "cm")
      }
    }
  }
}



if(fix_r_ml_e_time_vs_N & create_graphs) {
  
  for (r in c(1,2)) {
    for (ml in c(0.1,0.2,0.9)) {
      for (e in c(0.2)) {     
        time_data_set=time_data[time_data$RP == r &time_data$ML==ml & time_data$eps ==0.1 & time_data$T < 11,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_N <- ggplot(data=time_data_set ) +
          theme_bw()+
          expand_limits(y=0)+
          geom_line(aes(y = time_sync, x = N, color=factor(T))) +
          scale+
          ggtitle(substitute("M&S Synchronisation Time (RP="*r*", "~mu*"="*ml*")", list(r = bquote(.(r)), ml = bquote(.(ml)))))+
          labs( x = "N", y = "Sync. Time (in Cycles)", color="T")
        print(time_vs_N)
      }
    }
  }
}


if(fix_r_ml_e_time_vs_T & create_graphs) {
  
  for (r in c(1)) {
    for (ml in c(0.1,0.2,0.9)) {
      for (e in c(0.2)) {     
        time_data_set=time_data[time_data$RP == r &time_data$ML==ml & time_data$eps ==e & time_data$T < 11,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_N <- ggplot(data=time_data_set,aes(x = T, y =time_sync, colour=factor(N)) ) +
          theme_bw()+
          expand_limits(y=0)+
          geom_line()+
          scale+
          ggtitle(substitute("M&S Synchronisation Time (RP="*r*", "~epsilon*"="*e*"  "~mu*"="*ml*")", list(r = bquote(.(r)), ml = bquote(.(ml)), e=bquote(.(e)))))+
          labs( x = "T", y = "Sync. Time (Cycles)", color="N")
        print(time_vs_N)
        ggsave(sprintf("pdf/perez_rp%d_ml%f_e%f_time_vs_T.pdf",r,ml,e), plot=time_vs_N, device="pdf", width=gwidth, height=gheight, units = "cm")
      }
    }
  }
}

# graphs for meta-data like model checking time, etc.
if(checking_graphs)  {
  data_set=rawdata[rawdata$RP == 1 & rawdata$eps == .1 & rawdata$ML==0.1 & rawdata$T < 11,]
  
  mconst_vs_n <- ggplot(data=data_set  ) +
    geom_line(aes(x = data_set$N, y = (data_set$prob_mconst), color=factor(data_set$T))) +
    theme_bw()+
    scale+
    ggtitle(sprintf("Model Construction Time ")) +
    labs( x = "N ", y = "Time (s)", color="T") 
  
  mcheck_vs_n <- ggplot(data=data_set  ) +
    geom_line(aes(x = data_set$N, y = (data_set$prob_mcheck), color=factor(data_set$T))) +
    theme_bw()+
    scale+
    ggtitle(sprintf("Modelchecking Time ")) +
    labs( x = "N ", y = "Time (s)", color="T") 
  
  mem_vs_n <- ggplot(data=data_set  ) +
    geom_line(aes(x = N, y = (prob_mem/(1024*1024)), color=factor(T))) +
    theme_bw()+
    scale+
    ggtitle(sprintf("Modelchecking Memory Consumption ")) +
    labs( x = "N ", y = "Size (GB)", color="T") 
  
  if (create_graphs) {
    print(mconst_vs_n)
    print(mcheck_vs_n)
    print(mem_vs_n)
    
    ggsave(sprintf("pdf/perez_mconst_vs_N.pdf"), plot=mconst_vs_n, device="pdf", width=gwidth, height=gheight, units = "cm")
    ggsave(sprintf("pdf/perez_mcheck_vs_N.pdf"), plot=mcheck_vs_n, device="pdf", width=gwidth, height=gheight, units = "cm")
    ggsave(sprintf("pdf/perez_mem_vs_N.pdf"), plot=mem_vs_n, device="pdf", width=gwidth, height=gheight, units = "cm")
  }
  # create tables containing memory consumption during model checking
  if (create_tables) {
    mem_table <- data_set[data_set$T == 10, c("N", "prob_mem")]
    mem_table$mem_MB <- mem_table$prob_mem/(1024)
    mem_table <- mem_table[, c("N", "mem_MB")]
    mconst_table <- data_set[data_set$T == 10, c("N", "prob_mconst")]
    
  }
  
}


