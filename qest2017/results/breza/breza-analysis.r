library(ggplot2)
library(xtable)

#
# Load the model checking results. The first 5 columns are the parameters of 
# the model. Then we have 3 sections of 7 columns each. Within these sections, we have 
# the following structure:
# model checking result, model construction time, states of the model, transistions of the model, size of probability matrix, time for model checking, memory for model checking
#
# The first section contains the result for checking P=?[F synchronised], the second for R{"time_to_synch"}=? [F synchronised]
# The third section contains results we did not use for the analysis.



n3 = read.csv("breza_3.csv", 
              col.names=c("N","T","eps","RP", "ML", 
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem", 
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem" ,
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem" 
              ), header=FALSE)

n4 = read.csv("breza_4.csv", 
              col.names=c("N","T","eps","RP", "ML", 
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem", 
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem" 
              ),header=FALSE)

n5 = read.csv("breza_5.csv", 
              col.names=c("N","T","eps","RP", "ML", 
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem", 
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem" 
              ),header=FALSE)
n6 = read.csv("breza_6.csv", 
              col.names=c("N","T","eps","RP", "ML", 
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem", 
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem" 
              ),header=FALSE)
n7 = read.csv("breza_7.csv", 
              col.names=c("N","T","eps","RP", "ML", 
                          "prob_sync","prob_mconst",  "prob_states", "prob_transitions", "prob_matrix","prob_mcheck", "prob_mem", 
                          "time_sync","time_mconst", "time_states", "time_transitions", "time_matrix", "time_mcheck","time_mem",
                          "timef_sync","timef_mconst", "timef_states", "timef_transitions", "timef_matrix", "timef_mcheck","timef_mem" 
              ),header=FALSE)



# create a single dataframe containing the results
rawdata = rbind(n3, n4)
rawdata = rbind(rawdata,n5)
rawdata = rbind(rawdata,n6)
rawdata = rbind(rawdata,n7)

# some cleaning of the data: remove missing values (due to failure of prism) 
# and censordata for T > 10 :) (we did not get results for N=7 and T=11 in reasonable time)
prob_data <- rawdata[!is.na(rawdata$prob_sync) & rawdata$T < 11,]
time_data <- rawdata[!is.na(rawdata$time_sync) & rawdata$T < 11,]

# just out of curiosity, we save the missing values in new dataframes
prob_missingdata <- rawdata[is.na(rawdata$prob_sync),]
time_missingdata <- rawdata[is.na(rawdata$time_sync),]


# our own colour scheme, if colour is used in the graphs (most palettes do not contain enough colours)
col <- c("blue","red", "green","violet","black","orange",
         "steelblue4", "brown","yellowgreen", "darkgreen", 
         "palevioletred3", "salmon1", "chocolate4", "firebrick2",
         "grey43")

# size of the created pdfs
gwidth <-20
gheight <-15

# font size in created graphs 
gtext <- 18

# constants to define which graphs and tables should be created
create_graphs <- TRUE
create_tables <- TRUE

fix_n_t_vs_ML <- TRUE
fix_n_t_vs_RP <- TRUE
fix_t_r_time_vs_N <- TRUE
fix_n_r_time_vs_T <- TRUE
fix_n_ml_perc_rp <- TRUE
fix_r_ml_time_vs_T <- TRUE
checking_graphs <- TRUE

# greyscale or colour?
grey_scale <- TRUE

if (grey_scale) {
  scale <- scale_color_grey(end=0.7)
  
} else {
  scale <-  scale_colour_manual(values= col) 
}

# create graphs for visualising the model checking results

if(fix_n_t_vs_ML & create_graphs) {
  for (n in c(7)) {
    for (t in c(10)) {
      for (e in c(0.1)) {
        prob_data_set=prob_data[prob_data$N == n & prob_data$T == t & prob_data$eps ==e,]
        time_data_set=time_data[time_data$N == n & time_data$T == t & time_data$eps ==e,]
        
        prob_vs_ml <- ggplot(data=prob_data_set, aes(y = prob_sync, x = ML, color=factor(RP), shape=factor(RP))) +
          geom_line(aes(y = prob_sync, x = ML, color=factor(RP))) +
          theme_bw()+
          scale +
          ggtitle(substitute("MP Synchronisation Probability (N="*n*", T="*t*")", list(n = bquote(.(n)), t = bquote(.(t))))) +
          scale_shape_manual(values=1:nlevels(factor(prob_data$RP))) +
          geom_point(size=3)+
          labs( x = expression(mu), y = "Sync. Probability", colour="RP", shape="RP")  
        
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_ml <- ggplot(data=time_data_set , aes(y = time_sync, x = ML, color=factor(RP), shape=factor(RP))) +
          theme_bw()+
          guides(colour= guide_legend(keywidth = 2))+
          theme(legend.position=c(0.03,.94), legend.justification = c(0.05,0.85),text = element_text(size=gtext))+
          geom_line() +
          expand_limits(y=0)+
          scale+
          ggtitle(substitute("MP Synchronisation Time (N="*n*", T="*t*")", list(n = bquote(.(n)), t = bquote(.(t))))) +
          scale_shape_manual(values=1:nlevels(factor(prob_data$RP))) +
          geom_point(size=3)+
          labs( x = expression(mu), y = "Sync. Time (in Cycles)", colour="RP", shape="RP") 
        print(prob_vs_ml)
        print(time_vs_ml)
        ggsave(sprintf("pdf/breza_n%d_t%d_prob_vs_ml.pdf",n,t), plot=prob_vs_ml, device="pdf", width=gwidth, height=gheight, units = "cm")
        ggsave(sprintf("pdf/breza_n%d_t%d_time_vs_ml.pdf",n,t), plot=time_vs_ml, device="pdf", width=gwidth, height=gheight, units = "cm")
      }
    }
  }
  
}

if(fix_n_t_vs_RP & create_graphs) {
  for (n in  c(7)) {
    for (t in c(10)) {
      for (e in c(0.1)) {
        prob_data_set=prob_data[prob_data$N == n & prob_data$T == t & prob_data$eps ==e,]
        time_data_set=time_data[time_data$N == n & time_data$T == t & time_data$eps ==e,]
        
        prob_vs_rp <- ggplot(data=prob_data_set, aes( y = prob_sync, x = RP, color=factor(ML) , shape=factor(ML)) ) +
          theme_bw()+
          theme(legend.position=c(.15,.6), legend.justification = c(0.95,0.85),text = element_text(size=gtext))+
          guides(colour= guide_legend(keywidth = 2))+
          geom_line() +
          geom_point(size=3)+
          scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14))+
          scale_y_continuous(position="right")+
          scale_shape_manual(values=1:nlevels(factor(prob_data$ML))) +
          scale+
          ggtitle(substitute("MP Synchronisation Probability (N="*n*", T="*t*")", list(n = bquote(.(n)), t = bquote(.(t))))) +
          labs( x = "RP", y = "Sync. Probability", color=expression(mu), shape=expression(mu))
        
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        time_vs_rp <- ggplot(data=time_data_set ) +
          geom_line(aes(y = time_sync, x = RP, color=factor(ML))) +
          theme_bw()+
          expand_limits(y=0)+
          scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14))+
          scale+
          ggtitle(substitute("MP Synchronisation Time (N="*n*", T="*t*")", list(n = bquote(.(n)), t = bquote(.(t))))) +
          labs( x = "RP", y = "Sync. Time (in Cycles)", color=expression(mu))
        
        print(prob_vs_rp)
        print(time_vs_rp)
        ggsave(sprintf("pdf/breza_n%d_t%d_prob_vs_rp.pdf",n,t), plot=prob_vs_rp, device="pdf", width=gwidth, height=gheight, units = "cm")
        ggsave(sprintf("pdf/breza_n%d_t%d_time_vs_rp.pdf",n,t), plot=time_vs_rp, device="pdf", width=gwidth, height=gheight, units = "cm")
        
        
      }
    }
  }
}

if(fix_t_r_time_vs_N & create_graphs) {
  
  for (t in c(10)) {
    for (r in c(1,4)) {
      for (e in c(0.1)) {
        time_data_set=time_data[time_data$T == t &time_data$RP==r & time_data$eps ==e,]
        time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
        # we use a special treatment for refractory period of 1, to create nicer graph in the paper
        if (r == 1) {
          time_vs_N <- ggplot(data=time_data_set,aes(y = time_sync, x = N, color=factor(ML), shape=factor(ML)) ) +
            theme_bw()+
            theme(legend.position=c(.45,.95), legend.justification = c(0.95,0.85),legend.direction = "horizontal", text = element_text(size=gtext))+
            guides(colour= guide_legend(keywidth = 2))+
            scale_shape_manual(values=1:nlevels(factor(prob_data$ML))) +
            expand_limits(y=c(0,30))+
            geom_line() +
            geom_point(size=3)+
            scale+
            ggtitle(substitute("MP Synchronisation Time (T="*t*", RP="*r*")", list(t = bquote(.(t)), r = bquote(.(r))))) +
            labs( x = "N", y = "Sync. Time (in Cycles)", color=expression(mu), shape=expression(mu))
        } else {
          time_vs_N <- ggplot(data=time_data_set ) +
            theme_bw()+
            expand_limits(y=0)+
            geom_line(aes(y = time_sync, x = N, color=factor(ML))) +
            scale+
            ggtitle(substitute("MP Synchronisation Time (T="*t*", RP="*r*")", list(t = bquote(.(t)), r = bquote(.(r))))) +
            labs( x = "N", y = "Sync. Time (in Cycles)", color=expression(mu))
          
        }
        
        print(time_vs_N)
        ggsave(sprintf("pdf/breza_t%d_rp%d_time_vs_N.pdf",t,r), plot=time_vs_N, device="pdf", width=gwidth, height=gheight, units = "cm")
        
      }
    }
  }
}

if (fix_n_r_time_vs_T & create_graphs) {
  for (n in c(7)) {
    for (r in c(1)) {
      time_data_set=time_data[time_data$N == n &time_data$RP==r & time_data$e ==0.1,]
      time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
      time_vs_T <- ggplot(data=time_data_set ) +
        geom_line(aes(y = time_sync, x = T, color=factor(ML))) +
        theme_bw()+
        expand_limits(y=0)+
        scale+
        scale_y_continuous(limits = c(0,NA))+
        ggtitle(substitute("MP Synchronisation Time (N="*np*", RP="*r*")", list(np = bquote(.(n)), r = bquote(.(r))))) +
        labs( x = "T", y = "Sync. Time (in Cycles)", color=expression(mu))+
        theme_bw()
      
      print(time_vs_T)
      ggsave(sprintf("pdf/breza_n%d_rp%d_time_vs_T.pdf",n,r), plot=time_vs_T, device="pdf", width=gwidth, height=gheight, units = "cm")
      
      
      
    }
  }
  
}

if (fix_n_ml_perc_rp & create_graphs) {
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
          ggtitle(substitute("MP Synchronisation Probability (N="*np*", "~mu*"="* ml*")", list(np = bquote(.(n)), ml=bquote(.(ml))))) +
          labs( x = "RP (percentage of T)", y = "Sync. Probability", color="T")
        
        time_vs_rp_perc <- ggplot(data= time_data_set ) +
          geom_line(aes(y = time_sync, x = ( RP/  T * 100), color=factor( T))) +
          expand_limits(y=0)+
          theme_bw()+
          scale+
          ggtitle(substitute("MP Synchronisation Time (N="*np*", "~mu*"="* ml*")", list(np = bquote(.(n)), ml=bquote(.(ml))))) +
          labs( x = "RP (percentage of T)", y = "Sync. Time (in Cycles)", color="T")
        
        
        print(prob_vs_rp_perc)
        print(time_vs_rp_perc)
        ggsave(sprintf("pdf/breza_n%d_ml%f_prob_vs_rp_perc.pdf",n,ml), plot=prob_vs_rp_perc, device="pdf", width=gwidth, height=gheight, units = "cm")
        ggsave(sprintf("pdf/breza_n%d_ml%f_time_vs_rp_perc.pdf",n,ml), plot=time_vs_rp_perc, device="pdf", width=gwidth, height=gheight, units = "cm")
      }
    }
  }
}


if(fix_r_ml_time_vs_T & create_graphs) {
  
  for (r in c(1,2)) {
    for (ml in c(0.1,0.2,0.9)) {
      time_data_set=time_data[time_data$RP == r &time_data$ML==ml & time_data$eps ==0.1 & time_data$T < 11,]
      time_data_set <- time_data_set[time_data_set$time_sync != "Inf",]
      time_vs_T <- ggplot(data=time_data_set ) +
        theme_bw()+
        expand_limits(y=0)+
        geom_line(aes(y = time_sync, x = T, color=factor(N))) +
        scale+
        ggtitle(substitute("MP Synchronisation Time (RP="*r*", "~mu*"="*ml*")", list(r = bquote(.(r)), ml = bquote(.(ml)))))+
        labs( x = "T", y = "Sync. Time (in Cycles)", color="N")
      
      print(time_vs_T)
      ggsave(sprintf("pdf/breza_r%d_ml%f_time_vs_T.pdf",r,ml), plot=time_vs_T, device="pdf", width=gwidth, height=gheight, units = "cm")
      
    }
  }
}

# graphs for meta-data like model checking time, etc.
if(checking_graphs) {
  data_set=rawdata[rawdata$RP == 1 & rawdata$eps == .1 & rawdata$ML==0.1 & rawdata$T < 11,]
  
  mconst_vs_n <- ggplot(data=data_set  ) +
    geom_line(aes(x = data_set$N, y = (data_set$prob_mconst), color=factor(data_set$T))) +
    theme_bw()+
    scale+
    ggtitle(sprintf("MP Model Construction Time ")) +
    labs( x = "N ", y = "Time (s)", color="T") 
  
  mcheck_vs_n <- ggplot(data=data_set  ) +
    geom_line(aes(x = data_set$N, y = (data_set$prob_mcheck), color=factor(data_set$T))) +
    theme_bw()+
    scale+
    ggtitle(sprintf("MP Modelchecking Time ")) +
    labs( x = "N ", y = "Time (s)", color="T") 
  
  mem_vs_n <- ggplot(data=data_set  ) +
    geom_line(aes(x = N, y = (prob_mem/(1024*1024)), color=factor(T))) +
    theme_bw()+
    scale+
    ggtitle(sprintf("MP Modelchecking Memory Consumption ")) +
    labs( x = "N ", y = "Size (GB)", color="T") 
  
  if (create_graphs) {
    print(mconst_vs_n)
    print(mcheck_vs_n)
    print(mem_vs_n)
    
    ggsave(sprintf("pdf/breza_mconst_vs_N.pdf",n,t), plot=mconst_vs_n, device="pdf", width=gwidth, height=gheight, units = "cm")
    ggsave(sprintf("pdf/breza_mcheck_vs_N.pdf",n,t), plot=mcheck_vs_n, device="pdf", width=gwidth, height=gheight, units = "cm")
    ggsave(sprintf("pdf/breza_mem_vs_N.pdf",n,t), plot=mem_vs_n, device="pdf", width=gwidth, height=gheight, units = "cm")
  }
  # create tables containing memory consumption during model checking
  
  if (create_tables) {
    mem_table <- data_set[data_set$T == 10, c("N", "prob_mem")]
    mem_table$mem_MB <- mem_table$prob_mem/(1024)
    mem_table <- mem_table[, c("N", "mem_MB")]
    mconst_table <- data_set[data_set$T == 10, c("N", "prob_mconst")]
    
  }
  
  
}

