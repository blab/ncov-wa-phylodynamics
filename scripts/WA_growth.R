
library(ggplot2)
library("coda")
library("colorblindr") 
library(grid)
library(gridExtra)
library(ggtree)

path = "/Users/nmueller/Documents/github/hCoV-19_WA"

end_time = as.Date('2020-01-31')
# read in the mrsi file
mrsi = read.table(paste(path,'/results/mrsi.tsv',sep="/"), sep="\t", header=T)

mrsi$date = as.Date(mrsi$mrsi)

clusters = read.table(paste(path,'/results/cluster_size.tsv',sep="/"), sep="\t", header=T)

# timeframename = c('samples until 2020-03-10', 'samples until 2020-03-24')
timeframename = c('WA w/o Yakima', '614D', '614G', 'Yakima')

# how many days to ignore before the mrsi
# cutoff.mrsi = c(21,28,21,21)
cutoff.mrsi = c(7,7,7,7)

# dispersion parameter (for offspring distribution)
k = 0.3

# R from exponential growth: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1766383/pdf/rspb20063754.pdf
# R = 1+r/b
# set becomin uninfectious rate
becoming_uninf = 52.1429

# translation to I from: https://academic.oup.com/mbe/article/34/11/2982/3952784
# I = (Ne * R *(1 + 1/k))/gen_time
# set the generation time
generation_time = 1/becoming_uninf

# reporting delay of cases after infection https://science.sciencemag.org/content/early/2020/03/24/science.abb3221
reporting_delay = generation_time/2 * 366

# define the colors for bd and coal methods
both_clades_col = "#268457"
d_clade_col = "#0072B2"
g_clade_col = "#D55E00"



coal_col = "#56B4E9"
bdsky_col = "#009E73"
test_col = "#D55E00"

yakima_col = "#2F4F4F"

# define dates for vlines
vline1 = as.Date('2020-02-29') # report of community transmission in the greater Seattle area
vline2 = as.Date('2020-03-11') # intial official actions to combat coronavirus
vline3 = as.Date('2020-03-23') # stay home stay safe order


# combine runs
# logs = list.files(path=paste(path, "out", sep="/"), pattern="*rep0.log", full.names=T)
# for (i in seq(1,length(logs))){
#     in_command <- " -b 20 -log"
#     for (j in seq(0,2)){
#       in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), logs[i]), sep="")
#     }
#     out_command = gsub("_rep0", "", logs[i])
#     # combine the trees
#     system(paste("/Applications/BEAST\\ 2.6.0/bin/logcombiner", in_command, "-o", out_command, "", sep=" "))
# }

# combine runs
# logs = list.files(path=paste(path, "out", sep="/"), pattern="*rep0.*.trees", full.names=T)
# for (i in seq(1,length(logs))){
#     in_command <- " -b 20 -log"
#     for (j in seq(0,2)){
#       in_command = paste(in_command, " ", gsub("rep0", paste("rep", j,sep=""), logs[i]), sep="")
#     }
#     out_command = gsub("_rep0", "", logs[i])
#     # combine the trees
#     system(paste("/Applications/BEAST\\ 2.6.0/bin/logcombiner", in_command, "-o", out_command, "", sep=" "))
# }



# Reads in the analyses


# onset_cases = read.table(paste(path,'/data/PUBLIC_CDC_Event_Date_SARS_20200731.csv',sep=""), sep=",", header=T)
# onset_cases$Day = as.Date(onset_cases$WeekStartDate)
# testing = data.frame(date=unique(onset_cases$Day))
# testing$count=0
# for (i in seq(1,length(testing$date))){
#   testing[i,"count"] = sum(onset_cases[which(onset_cases$Day==testing[i,"date"]), "NewPos_All"])
# }




onset_cases = read.table(paste(path,'/data/PUBLIC_Tests_by_Specimen_Collection.csv',sep=""), sep=",", header=T)
onset_cases$Day = as.Date(onset_cases$Day)
testing = data.frame(date=unique(onset_cases$Day))
testing$count=0
for (i in seq(1,length(testing$date))){
  testing[i,"count"] = sum(onset_cases[which(onset_cases$Day==testing[i,"date"]), "Positive"])
}

first = T
first.growth = T
first.intro = T
for (i in seq(1,length(mrsi$filename))){
  time_diff = mrsi[i,"date"]-end_time
  time = seq(0,as.numeric(time_diff),3.5)
  
  if (!startsWith( as.character(mrsi[i,"filename"]),"multibd")){
    # read in the log file
    t = read.table(paste(path,'/out/', mrsi[i,'filename'], '.log',sep=""), sep="\t", header=T)
    # read in all the Ne's
    for (j in seq(1,length(time))){
      method = strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[2]]
      timeframe = timeframename[as.numeric(strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[3]])]
      name = paste('Ne',j,sep="")
      hpdInt = HPDinterval(exp(as.mcmc(t[,name])))
      hpdInt.m = HPDinterval(exp(as.mcmc(t[,name])), prob=0.5)
      
      date_val = mrsi[i,"date"]-time[j]

      test.numbers = NA
      if (length(testing[which(testing$date==date_val), "count"])==1){
        test.numbers = testing[which(testing$date==date_val), "count"]
      }
      
      
      new.dat = data.frame(time=date_val, 
                           Ne.mean=median(exp(t[,name])), Ne.lower=hpdInt[1,'lower'], Ne.upper=hpdInt[1,'upper'], 
                           Ne.ll=hpdInt.m[1,'lower'], Ne.uu=hpdInt.m[1,'upper'], 
                           nr_cases =test.numbers,
                           method = method, timeframe=timeframe)
      if (first){
        dat = new.dat
        first = F
      }else{
        dat = rbind(dat, new.dat)
      }
    }
    
    # get all the growth rates
    average_over = 1
    for (j in seq(1,length(time)-average_over,1)){
      method = strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[2]]
      timeframe = timeframename[as.numeric(strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[3]])]
      
      name1 = paste('Ne',j,sep="")
      name2 = paste('Ne',j+average_over,sep="")
      
      # get the growth rates
      values = (t[,name1]-t[,name2])/(time[j+average_over]-time[j])*365
      
      doubling = log(2)/values
      R = 1+values/becoming_uninf
      R[which(R<=0)] = 0
      hpdInt = HPDinterval(as.mcmc(values))
      hpdInt.5 = HPDinterval(as.mcmc(values), prob=0.5)
      hpdInt.doubling = HPDinterval(as.mcmc(doubling))
      hpdInt.R = HPDinterval(as.mcmc(R))
      hpdInt.R.5 = HPDinterval(as.mcmc(R), prob=0.5)
      
      # get the percentage of intros
      name.intro = paste('immigrationRate',floor((j-1)/2)+1,sep="")
      values.intro = exp(t[,name.intro])
      
      # get the approximate transmission rate from the growth rate
      transmission = values + becoming_uninf
      # assume that the minimal R0 is .5
      transmission[which(transmission<becoming_uninf/2)] = becoming_uninf/2
      # get the ratio of intros
      ratio.intro = exp(t[,name.intro])/(exp(t[,name.intro]) + transmission)
      
      hpd.intro = HPDinterval(as.mcmc(ratio.intro))
      hpd.intro.5 = HPDinterval(as.mcmc(ratio.intro), prob=0.5)
      
      hpd.force = HPDinterval(as.mcmc(values.intro/exp(t[,name1])))
      hpd.force.5 = HPDinterval(as.mcmc(values.intro/exp(t[,name1])), prob=0.5)

      

      new.dat = data.frame(time=mrsi[i,"date"]-time[j], 
                           growth.mean=median(values), growth.lower=hpdInt[1,'lower'], growth.upper=hpdInt[1,'upper'], 
                           growth.ll=hpdInt.5[1,'lower'], growth.uu=hpdInt.5[1,'upper'], 
                           doubling.mean=median(doubling), doubling.lower=hpdInt.doubling[1,'lower'], doubling.upper=hpdInt.doubling[1,'upper'], 
                           R.lower=hpdInt.R[1,'lower'], R.upper=hpdInt.R[1,'upper'], 
                           R.ll=hpdInt.R.5[1,'lower'], R.uu=hpdInt.R.5[1,'upper'], 
                           intro.l=hpd.intro[1,'lower'], intro.u=hpd.intro[1,'upper'],
                           intro.ll=hpd.intro.5[1,'lower'], intro.uu=hpd.intro.5[1,'upper'],
                           force.l=hpd.force[1,'lower'], force.u=hpd.force[1,'upper'],
                           force.ll=hpd.force.5[1,'lower'], force.uu=hpd.force.5[1,'upper'],
                           method = method, timeframe=timeframe)
      
      new.dat = rbind(new.dat, data.frame(time=mrsi[i,"date"]-time[j+1] +0.01, 
                           growth.mean=median(values), growth.lower=hpdInt[1,'lower'], growth.upper=hpdInt[1,'upper'], 
                           growth.ll=hpdInt.5[1,'lower'], growth.uu=hpdInt.5[1,'upper'], 
                           doubling.mean=median(doubling), doubling.lower=hpdInt.doubling[1,'lower'], doubling.upper=hpdInt.doubling[1,'upper'], 
                           R.lower=hpdInt.R[1,'lower'], R.upper=hpdInt.R[1,'upper'], 
                           R.ll=hpdInt.R.5[1,'lower'], R.uu=hpdInt.R.5[1,'upper'], 
                           intro.l=hpd.intro[1,'lower'], intro.u=hpd.intro[1,'upper'],
                           intro.ll=hpd.intro.5[1,'lower'], intro.uu=hpd.intro.5[1,'upper'],
                           force.l=hpd.force[1,'lower'], force.u=hpd.force[1,'upper'],
                           force.ll=hpd.force.5[1,'lower'], force.uu=hpd.force.5[1,'upper'],
                           method = method, timeframe=timeframe))
  
      if (first.growth){
        growth = new.dat
        first.growth = F
      }else{
        growth = rbind(growth, new.dat)
      }
    }
  }
}

# compute the growth rates from testing
positive = testing
average_over = 7
for (i in seq(2,length(positive$date)-average_over, 1)){
  growt_val = 0
  for (j in seq(1,average_over)){
    growt_val = growt_val + log(positive[i+j,"count"]/positive[i+j-1,"count"])/(1/366)
  }
  time = positive[i,"date"]+average_over/2
  growt_val = growt_val/average_over
  new.dat = data.frame(time=time, growth=growt_val)
  # new.dat = rbind(new.dat,data.frame(time=positive[i,"date"]+1-0.001, growth=growt_val))
  if (i==2){
    growth.testing = new.dat
  }else{
    growth.testing = rbind(growth.testing, new.dat)
  }
}

# read in the mrsi file
first = T
for (i in seq(1,length(mrsi$filename))){
  time_diff = mrsi[i,"date"]-end_time
  time = seq(0,as.numeric(time_diff),3.5)
  
  if (startsWith( as.character(mrsi[i,"filename"]),"multibd")){
    # read in the log file
    t = read.table(paste(path,'/out/', mrsi[i,'filename'], '.log',sep=""), sep="\t", header=T)
    # read in all the Ne's
    for (j in seq(1,length(time))){
      method = strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[2]]
      timeframe = timeframename[as.numeric(strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[3]])]
      name = paste('logReproductiveNumber',length(time)-j+1,sep="")
      # name2 = paste('samplingProportion',length(time)-j+1,sep="")
      
      hpdInt = HPDinterval(exp(as.mcmc(t[,name])))
      hpdInt.5 = HPDinterval(exp(as.mcmc(t[,name])),prob=0.5)
      # hpdInt.samp = HPDinterval(as.mcmc(t[,name2]))
      hpdInt.growth = HPDinterval(as.mcmc((exp(t[,name])-1)*becoming_uninf))
      hpdInt.growth.5 = HPDinterval(as.mcmc((exp(t[,name])-1)*becoming_uninf),prob=0.5)
      
      new.dat = data.frame(time=mrsi[i,"date"]-time[j], 
                           Ne.mean=median(t[,name]), Ne.lower=hpdInt[1,'lower'], Ne.upper=hpdInt[1,'upper'], 
                           Ne.ll=hpdInt.5[1,'lower'], Ne.uu=hpdInt.5[1,'upper'], 
                           # samp.mean=median(t[,name2]), samp.lower=hpdInt.samp[1,'lower'], samp.upper=hpdInt.samp[1,'upper'], 
                           growth.lower=hpdInt.growth[1,'lower'], growth.upper=hpdInt.growth[1,'upper'], 
                           growth.ll=hpdInt.growth.5[1,'lower'], growth.uu=hpdInt.growth.5[1,'upper'], 
                           
                           method = method, timeframe=timeframe)
      new.dat = rbind(new.dat, data.frame(time=mrsi[i,"date"]-time[j]-2+0.001, 
                       Ne.mean=median(t[,name]), Ne.lower=hpdInt[1,'lower'], Ne.upper=hpdInt[1,'upper'], 
                       Ne.ll=hpdInt.5[1,'lower'], Ne.uu=hpdInt.5[1,'upper'], 
                       # samp.mean=median(t[,name2]), samp.lower=hpdInt.samp[1,'lower'], samp.upper=hpdInt.samp[1,'upper'], 
                       growth.lower=hpdInt.growth[1,'lower'], growth.upper=hpdInt.growth[1,'upper'], 
                       growth.ll=hpdInt.growth.5[1,'lower'], growth.uu=hpdInt.growth.5[1,'upper'], 
                       method = method, timeframe=timeframe))

      if (first){
        dat.bdsky = new.dat
        first = F
      }else{
        dat.bdsky = rbind(dat.bdsky, new.dat)
      }
    }
  }
}





# comparison between number of samples and the effective population sizes.

transform = 125
text_height = 950

dat$timeframe <- factor(dat$timeframe, levels = c("614D", "614G", "WA w/o Yakima", "Yakima"))

dat_tmp = dat
growth_tmp = growth
for (i in seq(1,length(levels(mrsi$clade)))){
  ind = which(mrsi$clade==levels(mrsi$clade)[[i]])
  end = mrsi[ind[[1]], "date"]-cutoff.mrsi[[i]]
  clade=timeframename[[i]]
  dat_tmp = dat_tmp[-which(dat_tmp$timeframe==clade & dat_tmp$time>end),]
  growth_tmp = growth_tmp[-which(growth_tmp$timeframe==clade & growth_tmp$time>end),]
}


onset_cases = read.table(paste(path,'/data/PUBLIC_Tests_by_Specimen_Collection.csv',sep=""), sep=",", header=T)
onset_cases$Day = as.Date(onset_cases$Day)

onset_cases_yakima = onset_cases[which(onset_cases$County=="Yakima County"),]

onset_cases = onset_cases[-which(onset_cases$County=="Yakima County"),]
testing = data.frame(date=unique(onset_cases$Day))
testing$count=0
for (i in seq(1,length(testing$date))){
  testing[i,"count"] = sum(onset_cases[which(onset_cases$Day==testing[i,"date"]), "Positive"])
}
testing_yakima = data.frame(date=unique(onset_cases_yakima$Day))
testing_yakima$count=0
for (i in seq(1,length(testing_yakima$date))){
  testing_yakima[i,"count"] = sum(onset_cases_yakima[which(onset_cases_yakima$Day==testing_yakima[i,"date"]), "Positive"])
}



p = ggplot(dat_tmp[which(dat_tmp$method=="skygrowth" & dat_tmp$timeframe!="Yakima"),]) + 
  geom_histogram(data=testing, aes(x=date,y=count), stat="identity", color="grey", fill="grey")+
  geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)+
  # geom_line(aes(x=time, y=Ne.mean*transform, color=timeframe)) +
  geom_ribbon(aes(x=time, ymin=Ne.lower*transform, ymax=Ne.upper*transform, fill=timeframe), alpha=0.2)+
  geom_ribbon(aes(x=time, ymin=Ne.ll*transform, ymax=Ne.uu*transform, fill=timeframe), alpha=0.8)+
  scale_x_date(limits=c(as.Date('2020-01-31'), max(mrsi$date)))  +
  scale_fill_manual(name="Clade:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col))+
  theme_minimal() +
  scale_y_continuous(name="positive tests ", sec.axis = sec_axis(~ 1/transform*., name = "effective population size")) +
  theme( axis.title.y.right = element_text( angle = 90))+
  coord_cartesian(ylim=c(0.01,1000)) +
    geom_text(aes(x=vline1, label="local spread announced", y=text_height), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
    geom_text(aes(x=vline2, label="initial state wide measures", y=text_height*0.85), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
      geom_text(aes(x=vline3, label="state wide lockdown", y=text_height*0.7), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
  xlab("")
plot(p)
ggsave(plot=p, file=paste(path, 'figures/coal_Ne_vs_test_withlegend.pdf', sep='/'), height=2.5, width=6)
ggsave(plot=p + theme(legend.position="none"), file=paste(path, 'figures/coal_Ne_vs_test.pdf', sep='/'), height=2.5, width=6)

transform = 125
text_height = 250


p = ggplot(dat_tmp[which(dat_tmp$method=="skygrowth" & dat_tmp$timeframe=="Yakima"),]) + 
    geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)+
  # geom_line(aes(x=time, y=Ne.mean*transform, color=timeframe)) +
  geom_histogram(data=testing_yakima, aes(x=date,y=count), stat="identity")+
  geom_ribbon(aes(x=time, ymin=Ne.lower*transform, ymax=Ne.upper*transform, fill=timeframe), alpha=0.2)+
  geom_ribbon(aes(x=time, ymin=Ne.ll*transform, ymax=Ne.uu*transform, fill=timeframe), alpha=0.8)+
  scale_x_date(limits=c(as.Date('2020-01-31'), max(mrsi$date)))  +
  theme_minimal() +
  scale_fill_manual(name="Clade:", values = c("Yakima"=yakima_col))+
  scale_y_continuous(name="positive tests ", sec.axis = sec_axis(~ 1/transform*., name = "effective population size")) +
  theme( axis.title.y.right = element_text( angle = 90))+
  coord_cartesian(ylim=c(0.01,300)) +
    geom_text(aes(x=vline1, label="local spread announced", y=text_height), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
    geom_text(aes(x=vline2, label="initial state wide measures", y=text_height*0.85), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
      geom_text(aes(x=vline3, label="state wide lockdown", y=text_height*0.7), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
  xlab("")
plot(p)
ggsave(plot=p+ theme(legend.position="none"), file=paste(path, 'figures/coal_Ne_vs_test_yakima.png', sep='/'), height=2.5, width=6)


p = ggplot(dat_tmp[which(dat_tmp$timeframe==timeframename[[1]]),]) + 
  geom_point(aes(x=nr_cases, y=Ne.mean, color=time)) +
  geom_errorbar(aes(x=nr_cases, ymin=Ne.lower, ymax=Ne.upper, color=time)) +
  theme_minimal() +
  facet_wrap(method~., ncol=1) +
  geom_smooth(aes(x=nr_cases, y=Ne.mean), method="lm")+
  xlab("number of cases") +
  ylab("effective population size")+
  scale_y_log10()+
  scale_x_log10()
  
plot(p)




# Comparison between the different runs

p = ggplot(dat) + 
  geom_ribbon(aes(x=time, ymin=log(Ne.lower), ymax=log(Ne.upper), fill=timeframe), alpha=0.2)+
  geom_ribbon(aes(x=time, ymin=log(Ne.ll), ymax=log(Ne.uu), fill=timeframe), alpha=0.8)+
  facet_wrap(method~., ncol=1) +
  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"=yakima_col))+
  scale_y_continuous(breaks=seq(-6.9078,6.9078,2.3026), labels=c(0.001,0.010,0.1,1,10,100,1000),sec.axis = sec_axis(~ .,breaks=seq(-6.9078,6.9078,2.3026), labels=c(0.001,0.010,0.1,1,10,100,1000), name = "effective population size")) +
  theme_minimal() +
  # theme(legend.position = "none") + 
  xlab("") +
  ylab("effective population size") +
  coord_cartesian(ylim=c(-6.9078,6.9078))
plot(p)
ggsave(plot=p, file=paste(path, 'figures/coal_smoothingcomp.png', sep='/'), height=7,width=12)

text_height = 3

# dat.bdsky = old.dat.bdsky
# 
# dat.bdsky[which(dat.bdsky$timeframe==timeframename[[3]]),"time"] = dat.bdsky[which(dat.bdsky$timeframe==timeframename[[3]]),"time"]-10

dat.bdsky_tmp = dat.bdsky
for (i in seq(1,length(levels(mrsi$clade)))){
  ind = which(mrsi$clade==levels(mrsi$clade)[[i]])
  end = mrsi[ind[[1]], "date"]-cutoff.mrsi[[i]]
  clade=timeframename[[i]]
  dat.bdsky_tmp = dat.bdsky_tmp[-which(dat.bdsky_tmp$timeframe==clade & dat.bdsky_tmp$time>end),]
}


p = ggplot(data=dat.bdsky_tmp) + 
  geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)+
    geom_text(aes(x=vline1, label="local spread announced", y=text_height), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
    geom_text(aes(x=vline2, label="initial state wide measures", y=text_height*0.85), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
      geom_text(aes(x=vline3, label="state wide lockdown", y=text_height*0.7), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
  geom_ribbon(data=dat.bdsky, aes(x=time, ymin=Ne.lower, ymax=Ne.upper, fill=timeframe), alpha=0.2) +
  geom_ribbon(data=dat.bdsky, aes(x=time, ymin=Ne.ll, ymax=Ne.uu, fill=timeframe), alpha=0.8) +
  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"=yakima_col))+
  theme_minimal() +
  # theme(legend.position = "none") + 
  xlab("") +
  ylab("effective reproduction number")+
  coord_cartesian(ylim=c(0,4))
plot(p)
ggsave(plot=p, file=paste(path, 'figures/bdsky_R0.pdf', sep='/'), height=5,width=12)


text_height = 8
p = ggplot(data=growth) + 
  geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)+
    geom_text(aes(x=vline1, label="local spread announced", y=text_height), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
    geom_text(aes(x=vline2, label="initial state wide measures", y=text_height*0.85), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
      geom_text(aes(x=vline3, label="state wide lockdown", y=text_height*0.7), colour="black", angle=0, text=element_text(size=6), size=4,hjust = 0) +
  
  geom_ribbon(aes(x=time, ymin=R.lower, ymax=R.upper, fill=timeframe), alpha=0.2) +
  geom_ribbon(aes(x=time, ymin=R.ll, ymax=R.uu, fill=timeframe), alpha=0.8) +
  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"= yakima_col))+
  theme_minimal() +
  # theme(legend.position = "none") + 
  xlab("") +
  theme_minimal() +
  ylab("effective reproduction number") +
    theme_minimal() +
  facet_wrap(method~., ncol=1) +
  coord_cartesian(ylim=c(0,5))
plot(p)
ggsave(plot=p, file=paste(path, 'figures/coal_R0.png', sep='/'), height=5,width=12)

# Plot the lineages through time for each introduction

  first = T
first.rep = T

for (i in seq(1,4)){
  filename = paste('multicoal_skygrowth_', i,sep='')
  
  # read in the log file 
  t = read.table(paste(path,'/out/', filename, '.log',sep=""), sep="\t", header=T)
  # take a 10% burnin
  t = t[-seq(1,length(t$posterior)/10),]
  
  # get all the intro times
  indices = which(clusters$filename==filename)
  
  # read in the ltt values
  ltt = read.table(paste(path,'/results/ltt_mean.tsv',sep=""), sep="\t", header=T)
  ltt = ltt[which(ltt$dataset==i),]
  
  ltt.clusters = clusters[which(clusters$filename==filename),]
  c = 1;
  current_max = 0;
  plot_offset = 0;
  
  indices = order(ltt[,2])
  
  
  for (l in seq(1,length(ltt$dataset))){
    j=indices[[l]]
    # get the offset in days 
    name = paste('Tree.', ltt[j,2], '.trueHeight',sep='')
    name2 = gsub("trueHeight","height", name) 
    offset  = round((t[1,name]-t[1,name2])*365)
    
    max_val = max(ltt[j,3:length(ltt)])
    # get the clade
    clade_label = clusters[which(clusters$filename==filename & clusters$number==ltt[j,2]),"clade"][[1]]
    
    for (k in seq(3,length(ltt))){
      new.dat = data.frame(time = mrsi[which(mrsi$filename==filename), "date"] - offset - k-1,
                           val = ltt[j,k], val_offset = current_max + max_val/2, cluster = ltt[j,2], clade=paste(clade_label, "clade"),timeframe=timeframename[[i]])
      
      if (first){
        ltt.dat = new.dat
        first = F
      }else{
        ltt.dat = rbind(ltt.dat, new.dat)
      }
    }
    current_max = current_max + max_val + plot_offset
  }

  # uni_times = unique(ltt.dat$time)
  # ltt.dat$rel = 0
  # 
  # ltt.dat$min_val = ltt.dat$val
  # ltt.dat[which(ltt.dat$min_val<1),"min_val"] = 0
  # for (i in seq(1,length(uni_times))){
  #   ind = which(ltt.dat$time==uni_times[i])
  #   ltt.dat[ind,"rel"] = ltt.dat[ind, "min_val"]/sum(ltt.dat[ind, "min_val"])
  # }
  # 
  # ltt.dat[which(ltt.dat$rel==0), "rel"] =0.001
  # p = ggplot(ltt.dat[which(!is.nan(ltt.dat$rel)),]) + 
  #   geom_area(aes(x=time,y=rel,fill=clade, group=cluster),position="stack", stat="identity") +
  #   scale_fill_manual(name="Clade:", values = c("both clades"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col))+
  #   ylab("relative number of lineages") +
  #   xlab("")+
  #   theme_minimal()
  # 
  # plot(p)
  # ggsave(plot=p, file=paste(path, 'figures/coal_ltt_rel.png', sep='/'), height=5,width=9)
}

ltt.dat$newcladeName = "614G"
ltt.dat[which(ltt.dat$clade=="D clade"), "newcladeName"] = "614D"
# ltt.dat[which(ltt.dat$clade=="G clade"), "newcladeName"] = "614G"


segment.dat = data.frame(x=as.Date("2020-06-7"), xend=as.Date("2020-06-7"), y=50, yend=450)
segment.dat = rbind(segment.dat, data.frame(x=as.Date("2020-06-7"), xend=as.Date("2020-06-6"), y=50, yend=50))
segment.dat = rbind(segment.dat, data.frame(x=as.Date("2020-06-7"), xend=as.Date("2020-06-6"), y=250, yend=250))
segment.dat = rbind(segment.dat, data.frame(x=as.Date("2020-06-7"), xend=as.Date("2020-06-6"), y=450, yend=450))
segment.dat$timeframe = levels(ltt.dat$timeframe)[[2]]

segment.text = data.frame(x=as.Date("2020-06-4"), label="0", y=50,angle=0)
segment.text = rbind(segment.text,data.frame(x=as.Date("2020-06-4"), label="200", y=250,angle=0))
segment.text = rbind(segment.text,data.frame(x=as.Date("2020-06-4"), label="400", y=450,angle=0))
segment.text = rbind(segment.text,data.frame(x=as.Date("2020-06-9"), label="nr lineages in cluster", y=250, angle=90))
segment.text$timeframe = levels(ltt.dat$timeframe)[[2]]


p = ggplot(ltt.dat) + geom_ribbon(aes(x=time, ymin=val_offset-val/2, ymax=val_offset+val/2, group=cluster, fill=newcladeName)) +
  scale_x_date(limits=c(as.Date('2020-02-01'), max(mrsi$date)))  +
  theme_minimal() +
  scale_fill_manual(name="Spike protein variant:", values = c("614D"=d_clade_col, "614G"=g_clade_col))+
  xlab("") +
  scale_y_continuous(breaks=NULL) +
  facet_wrap(.~timeframe, ncol=1)+
  geom_segment(data=segment.dat, aes(x=x, xend=xend, y=y, yend=yend)) +
  geom_text(data=segment.text,aes(x=x, label=label, y=y,angle=angle), colour="black", text=element_text(size=3))+
  ylab("")

plot(p)
ggsave(plot=p, file=paste(path, 'figures/coal_ltt.png', sep='/'), height=8,width=12)

# p_coal_growth = ggplot(growth[which(growth$method=="skygrowth"),]) + 
#   geom_ribbon(aes(x=time, ymin=growth.lower, ymax=growth.upper, fill="coalescent skyline"), alpha=0.2) +
#   geom_ribbon(aes(x=time, ymin=growth.ll, ymax=growth.uu, fill="coalescent skyline"), alpha=0.8) +
#   
#   geom_ribbon(data=dat.bdsky[which(dat.bdsky$timeframe=="both clades"),], aes(x=time, ymin=growth.lower, ymax=growth.upper, fill="birth-death skyline"), alpha=0.2) +
#   geom_ribbon(data=dat.bdsky[which(dat.bdsky$timeframe=="both clades"),], aes(x=time, ymin=growth.ll, ymax=growth.uu, fill="birth-death skyline"), alpha=0.8) +
#   geom_path(data=growth.testing, aes(x=time-reporting_delay, y=growth, color="confirmed positive tests\nby time of symptom onset\nminus 5 days")) +
#   scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
#   coord_cartesian(ylim=c(-100,300))+
# 
#   scale_color_manual(name="computed using:", values = c("confirmed positive tests\nby time of symptom onset\nminus 5 days"=test_col, "birth-death skyline"=bdsky_col, "coalescent skyline"=coal_col))+
#   scale_fill_manual(name="inferred using:", values = c("confirmed positive tests\nby time of symptom onset\nminus 5 days"=test_col, "birth-death skyline"=bdsky_col, "coalescent skyline"=coal_col))+
#   xlab("")+
#   theme_minimal() 
# doubling_labels = c(-1,-2,-6,6,2,1)
# p_coal_growth <- p_coal_growth + scale_y_continuous(sec.axis = sec_axis(~ .,breaks=round(log(2)/doubling_labels*365), labels=doubling_labels, name = "doubling times in days")) + ylab("growth rate per year") +
#     theme( axis.title.y.right = element_text( angle = 90))
# 
# plot(p_coal_growth)
# ggsave(plot=p_coal_growth, file=paste(path, 'figures/growth_rates.pdf', sep='/'), height=3,width= 7.8000)


p_R0 = ggplot(growth[which(growth$method=="skygrowth"),]) + 
  geom_ribbon(aes(x=time, ymin=R.lower, ymax=R.upper, fill="coalescent skygrowth"), alpha=0.2) +
  geom_ribbon(aes(x=time, ymin=R.ll, ymax=R.uu, fill="coalescent skygrowth"), alpha=0.8) +
  geom_ribbon(data=dat.bdsky, aes(x=time, ymin=Ne.lower, ymax=Ne.upper, fill="birth-death skyline"), alpha=0.2)+
  geom_ribbon(data=dat.bdsky, aes(x=time, ymin=Ne.ll, ymax=Ne.uu, fill="birth-death skyline"), alpha=0.8)+
  
  ylab("Effective Reproduction Number")+
  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  coord_cartesian(ylim=c(0,5))+
  geom_vline(xintercept=vline1)+
  geom_vline(xintercept=vline2)+
  geom_vline(xintercept=vline3)+
  scale_fill_manual(name="inference using:", values = c("confirmed positive tests"=test_col, "birth-death skyline"=bdsky_col, "coalescent skygrowth"=coal_col))+
  theme_minimal()+
  facet_wrap(.~timeframe, ncol=1) +
  xlab("")

plot(p_R0)

ggsave(plot=p_R0, file=paste(path, 'figures/R0_comparison.png', sep='/'), height=8,width=12)



# ```{r prevalence_coal}
# first = T
# 
# for (i in seq(1,length(mrsi$filename))){
#   time_diff = mrsi[i,"date"]-end_time
#   time = seq(0,as.numeric(time_diff),2)
#   
#   if (!startsWith( as.character(mrsi[i,"filename"]),"multibd")){
#     # read in the log file
#     t = read.table(paste(path,'/out/', mrsi[i,'filename'], '.log',sep=""), sep="\t", header=T)
#     # take a 10% burnin
#     t = t[-seq(1,length(t$posterior)/10),]
#     
#     # get all the prevalence estimates
#     average_over=1
#     for (j in seq(1,length(time)-average_over,1)){
#       method = strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[2]]
#       timeframe = timeframename[as.numeric(strsplit(as.character(mrsi[i,'filename']), split="_")[[1]][[3]])]
#       
#       name1 = paste('Ne',j,sep="")
#       name2 = paste('Ne',j+average_over,sep="")
#       
#       growth_rate_val = (t[,name1]-t[,name2])/(time[j+average_over]-time[j])*365
#       R = 1+growth_rate_val/(becoming_uninf*2)
#       # set all R smaller than 0 to 0
#       R[which(R<0.25)]=0.25
#       # get the Ne at the mid point for this time interval
#       Ne = exp(t[,name1])
#       
#       # compute the prevalence
#       I = (Ne*R*(1+1/k))/generation_time
#       
#       hpdInt = HPDinterval(as.mcmc(I))
#       hpdInt.narrow = HPDinterval(as.mcmc(I), prob=0.5)
# 
#       new.dat = data.frame(time=mrsi[i,"date"] - time[j], 
#                            I.mean=median(I), I.lower=hpdInt[1,'lower'], I.upper=hpdInt[1,'upper'],
#                            I.ll=hpdInt.narrow[1,'lower'], I.uu=hpdInt.narrow[1,'upper'],
#                            method = method, timeframe=timeframe)
#        if (first){
#         prev = new.dat
#         first = F
#       }else{
#         prev = rbind(prev, new.dat)
#       }
#  
#     }
#   }
# }
# 
# 
# p = ggplot()+
#   geom_ribbon(data=prev[which(prev$method=="skygrowth" & prev$timeframe==timeframename[[1]]),],aes(x=time, ymin=I.lower, ymax=I.upper, fill=timeframe), alpha=0.2)+
#   geom_ribbon(data=prev[which(prev$method=="skygrowth" & prev$timeframe==timeframename[[1]]),],aes(x=time, ymin=I.ll, ymax=I.uu, fill=timeframe), alpha=0.8)+
#   scale_x_date(limits=c(as.Date('2020-02-01'), max(mrsi$date)))  +
#   scale_y_continuous(sec.axis = sec_axis(~ .), name = "Prevalence") +
#   theme_minimal() +
#   theme(legend.position = "none")  +
#   scale_color_manual(name="data source", values=c("#56B4E9"))+
#   scale_fill_manual(name="data source", values=c("#56B4E9"))+
#   coord_cartesian(ylim=c(1,2000))+
#   scale_y_log10()
# plot(p)
# ggsave(plot=p, file=paste(path, 'figures/coal_prevalence.png', sep='/'), height=3,width=9)
# 
# 
# ```
# Comparison between the number of intros and samples
# g.perc_intros = ggplot(growth_tmp[which(growth_tmp$method=="skygrowth" & growth_tmp$timeframe==timeframename[[1]]),]) +
#   geom_ribbon(aes(x=time, ymin=intro.l, ymax=intro.u), alpha=0.2, fill=both_clades_col)+
#   geom_ribbon(aes(x=time, ymin=intro.ll, ymax=intro.uu), alpha=0.8, fill=both_clades_col)+
#   theme_minimal()+
#   facet_grid(.~timeframe)+
#   scale_y_log10()+
#   scale_x_date()+
#   ylab("Percentage of cases due to introductions") +
#   xlab("") 
# plot(g.perc_intros)
# ggsave(plot=g.perc_intros, file=paste(path, 'figures/intro_percentage.pdf', sep='/'), height=4,width=7)
# ggsave(plot=g.perc_intros, file=paste(path, 'figures/intro_percentage.png', sep='/'), height=4,width=7)


g.perc_intros = ggplot(growth_tmp[which(growth_tmp$method=="skygrowth" & growth_tmp$timeframe!="Yakima"),]) +
  geom_ribbon(aes(x=time, ymin=intro.l*100, ymax=intro.u*100, fill=timeframe), alpha=0.2, )+
  geom_ribbon(aes(x=time, ymin=intro.ll*100, ymax=intro.uu*100,fill=timeframe), alpha=0.8)+
  theme_minimal()+
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"= yakima_col))+
  scale_y_log10()+
  scale_x_date()+
  coord_cartesian(ylim=c(0.1,100))+
  ylab("Percentage of cases due to introductions") +
  xlab("") 
plot(g.perc_intros)
ggsave(plot=g.perc_intros+theme(legend.position = "none"), file=paste(path, 'figures/intro_percentage_clade.pdf', sep='/'), height=4,width=7)

g.perc_intros = ggplot(growth[which(growth$method=="skygrowth" & growth$timeframe!="WA w/o Yakima"& growth$timeframe!="Yakima"),]) +
  geom_ribbon(aes(x=time, ymin=force.l, ymax=force.u, fill=timeframe), alpha=0.2 )+
  geom_ribbon(aes(x=time, ymin=force.ll, ymax=force.uu,fill=timeframe), alpha=0.8)+
  theme_minimal()+
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"= yakima_col))+
  scale_x_date()+
  coord_cartesian(ylim=c(0,100))+
  ylab("proportional number introductions") +
  xlab("") 
plot(g.perc_intros)
ggsave(plot=g.perc_intros+theme(legend.position = "none"), file=paste(path, 'figures/intros_clade.png', sep='/'), height=4,width=7)


nrs = clusters[which(clusters$filename=="multicoal_skygrid_1"),]

# keeps track of the membership of samples
cluster_member = c()
for (i in seq(1,length(nrs$number))){
  cluster_member = append(cluster_member, rep(nrs[i,"number"], nrs[i,"size"]))
}
first = T
for (i in seq(0,length(cluster_member)/3 - 1,10)){
  nr_intros_samples = c()
  for (r in seq(1,1000)){
    y = sample(cluster_member, size=i+1, replace =F)
    with_last = length(unique(y))
    without_last = length(unique(head(y, -1)))
    nr_intros_samples = append(nr_intros_samples, with_last-without_last)
  }
  # hpd = HPDinterval(as.mcmc(nr_intros_samples))
  # hpd.5 = HPDinterval(as.mcmc(nr_intros_samples), prob=0.5)
  new.dat = data.frame(nr_samples = i, nr_intros = mean(nr_intros_samples))

  if (first){
    nr_intros = new.dat
    first = F;
  }else{
    nr_intros = rbind(nr_intros, new.dat)
  }
}

# get the simulated values
sims_intros = read.table(paste(path,'/results/sim_cluster_size.tsv',sep="/"), sep="\t", header=T)
sims_intros$intro_percentage = as.character(sims_intros$intro_percentage*100)

sims_intros$intro_percentage = factor(sims_intros$intro_percentage, levels = c("1", "3", "10"))

g.nrintros = ggplot() +
  geom_line(data = sims_intros, aes(x=nr_samples, y = intro_prob, group=interaction(run, intro_percentage), color = intro_percentage))+
  geom_line(data=nr_intros, aes(x=nr_samples-1, y = nr_intros, color="observed"))+
  theme_minimal()+
  scale_y_log10()+
  ylab("Probability that new sample reveals new introduction") +
  xlab("Number of samples") +
  scale_color_manual(name="percentage of cases\ndue to introductions", breaks=c("1","3","10","observed"),
                     values=c("1"="#F9C988", "3"="#71BD89", "10"="#24A897", "observed"=both_clades_col))
plot(g.nrintros)
ggsave(plot=g.nrintros, file=paste(path, 'figures/intro_clustersampling.png', sep='/'), height=4,width=7)


# regions given by https://www.doh.wa.gov/Portals/1/Documents/1200/phsd-PHEPR.pdf
region_map = data.frame(county="King", region="King")
region_map = rbind(region_map, data.frame(county="Snohomish", region="North"))
region_map = rbind(region_map, data.frame(county="Island", region="North"))
region_map = rbind(region_map, data.frame(county="Whatcom", region="North"))
region_map = rbind(region_map, data.frame(county="Skagit", region="North"))
region_map = rbind(region_map, data.frame(county="San Juan", region="North"))


region_map = rbind(region_map, data.frame(county="Pierce", region="Pierce"))

region_map = rbind(region_map, data.frame(county="Mason", region="North West"))
region_map = rbind(region_map, data.frame(county="Kitsap", region="North West"))
region_map = rbind(region_map, data.frame(county="Jefferson", region="North West"))
region_map = rbind(region_map, data.frame(county="Clallam", region="North West"))


region_map = rbind(region_map, data.frame(county="Thurston", region="West"))
region_map = rbind(region_map, data.frame(county="Lewis", region="West"))
region_map = rbind(region_map, data.frame(county="Pacific", region="West"))
region_map = rbind(region_map, data.frame(county="Grays Harbor", region="West"))

region_map = rbind(region_map, data.frame(county="Wahkiakum", region="Southwest"))
region_map = rbind(region_map, data.frame(county="Cowlitz", region="Southwest"))
region_map = rbind(region_map, data.frame(county="Skamania", region="Southwest"))
region_map = rbind(region_map, data.frame(county="Clark", region="Southwest"))


region_map = rbind(region_map, data.frame(county="Yakima", region="South Central"))
region_map = rbind(region_map, data.frame(county="Kittitas", region="South Central"))
region_map = rbind(region_map, data.frame(county="Benton", region="South Central"))
region_map = rbind(region_map, data.frame(county="Klickitat", region="South Central"))
region_map = rbind(region_map, data.frame(county="Franklin", region="South Central"))
region_map = rbind(region_map, data.frame(county="Walla Walla", region="South Central"))

region_map = rbind(region_map, data.frame(county="Adams", region="East"))
region_map = rbind(region_map, data.frame(county="Whitman", region="East"))
region_map = rbind(region_map, data.frame(county="Garfield", region="East"))
region_map = rbind(region_map, data.frame(county="Columbia", region="East"))
region_map = rbind(region_map, data.frame(county="Asotin", region="East"))
region_map = rbind(region_map, data.frame(county="Lincoln", region="East"))
region_map = rbind(region_map, data.frame(county="Spokane", region="East"))
region_map = rbind(region_map, data.frame(county="Ferry", region="East"))
region_map = rbind(region_map, data.frame(county="Stevens", region="East"))
region_map = rbind(region_map, data.frame(county="Pend Oreille", region="East"))

region_map = rbind(region_map, data.frame(county="Grant", region="North Central"))
region_map = rbind(region_map, data.frame(county="Okanogan", region="North Central"))
region_map = rbind(region_map, data.frame(county="Chelan", region="North Central"))
region_map = rbind(region_map, data.frame(county="Douglas", region="North Central"))

region_map$county = paste(region_map$county, "County")


mobility = read.delim(paste(path,'/data/Global_Mobility_Report.csv',sep=""), sep=",", header=T,na.strings=c("",".","NA"))

mobility = mobility[which(!is.nan(mobility$sub_region_2) & mobility$sub_region_1=="Washington"),]
mobility$date = as.Date(mobility$date)

mobility$sub_region_2 = factor(as.character(mobility$sub_region_2))

header.labs = labels(mobility)[[2]]




# only use the 4 main regions
mobility$region = ""
for (i in seq(1,length(mobility$region))){
  if (!is.na(mobility[i,"sub_region_2"])){
    mobility[i,"region"] = as.character(region_map[which(region_map$county==mobility[i,"sub_region_2"]),"region"])
  }
}
subset = mobility[which(mobility$sub_region_2=="Skagit County" | mobility$region=="King" | mobility$region=="Pierce" | mobility$sub_region_2=="Snohomish County" | mobility$sub_region_2=="Yakima County"),]

# define time period for which to look into correlation
start = as.Date('2020-02-29')
end = as.Date('2020-03-24')
# remove weekend data
sat = as.Date('2020-01-04')
sun = as.Date('2020-01-05')

mob.val = labels(subset)[[2]]

first = T

for (l in seq(10,length(mob.val)-2)){
  subset$val = NA
  
  for (j in seq(4,length(subset$country_region_code)-3)){
    subset[j,"val"]=mean(subset[seq(j-3,j+3),mob.val[[l]]])
  }

  king_county = subset[which(subset$sub_region_2=="King County"),]
  king_county$sub_region_2 = c()
  
  # compute lag between measures
  counties = unique(subset$sub_region_2)
  a=2
  require("Hmisc")
  for (b in seq(1,length(counties))){
    if (a!=b){
      for (i in seq(0,14)){
        c.vals.a = subset[which(subset$sub_region_2==counties[[a]]),]
        c.vals.b = subset[which(subset$sub_region_2==counties[[b]]),]
        indices = which(c.vals.a$date>=start & c.vals.a$date<=end)
        c.vals.a = c.vals.a[indices,]
        c.vals.b = c.vals.b[indices+i,]
        p=rcorr(c.vals.a$val,c.vals.b$val)
        tmp = abs(c.vals.a$val-c.vals.b$val)
        mean.diff = mean(tmp[which(!is.na(tmp))])
        new.dat = data.frame(lag=i, r=p$r[1,2], county=counties[[b]], diff=mean.diff,mobility = mob.val[[l]])
        if (first){
          correlation = new.dat
          first = F
        }else{
          correlation = rbind(correlation,new.dat)
        }

      }
    }
  }
}

subset$bestlag = NA
for (b in seq(1,length(counties))){
  if (a!=b){
    best.lag = correlation[which.max(correlation[which(correlation$county==counties[[b]]),"r"]),"lag"]
    subset[which(subset$sub_region_2==counties[[b]]), "bestlag"] = best.lag
  }
}

subset[which(subset$region=="North"),"region"] = "north western counties"
p = ggplot(data=dat.bdsky_tmp[which(dat.bdsky_tmp$timeframe!="Yakima"),]) + 
  geom_ribbon(aes(x=time, ymin=Ne.lower, ymax=Ne.upper, fill=timeframe), alpha=0.2) +
  geom_ribbon(aes(x=time, ymin=Ne.ll, ymax=Ne.uu, fill=timeframe), alpha=0.8) +
  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"=yakima_col))+
  scale_y_continuous(sec.axis = sec_axis(~ .*30 - 75, name="workplace mobility")) +
  geom_line(data=subset[which(subset$sub_region_2!="Yakima County"),], aes(x=date, y=val/30+2.5, linetype=sub_region_2,group=sub_region_2))+
  scale_linetype(name="Mobility in County:")+ 
  theme_minimal() +
  theme(legend.position="bottom", legend.box = "vertical") +
  xlab("") +
  ylab("effective reproduction number")+
  coord_cartesian(ylim=c(0,5)) +
    geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)
plot(p)


ggsave(plot=p + ylab("")+scale_y_continuous(sec.axis = sec_axis(~ ., name="")), file=paste(path, 'figures/bdsky_R0_mobility_withlegend.pdf', sep='/'), height=1.5,width=8)
ggsave(plot=p + theme(legend.position="right", legend.box = "horizontal"), file=paste(path, 'figures/bdsky_R0_mobility_withlegend2.pdf', sep='/'), height=4,width=6)
ggsave(plot=p + theme(legend.position="none"), file=paste(path, 'figures/bdsky_R0_mobility.pdf', sep='/'), height=2.5,width=6)

p = ggplot(data=growth_tmp[which(growth_tmp$method=="skygrowth" & growth_tmp$timeframe!="Yakima"), ]) + 

  geom_ribbon(aes(x=time, ymin=R.lower, ymax=R.upper, fill=timeframe), alpha=0.2) +
  geom_ribbon(aes(x=time, ymin=R.ll, ymax=R.uu, fill=timeframe), alpha=0.8) +

  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"=yakima_col))+
  scale_y_continuous(sec.axis = sec_axis(~ .*30 - 75, name="workplace mobility")) +
  geom_line(data=subset[which(subset$sub_region_2!="Yakima County"),], aes(x=date, y=val/30+2.5, linetype=sub_region_2,group=sub_region_2))+
  scale_linetype(name="Mobility in County:")+ 
  theme_minimal() +
  # theme(legend.position = "none") + 
  xlab("") +
  ylab("effective reproduction number")+
  coord_cartesian(ylim=c(0,5)) +
    geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)
plot(p)


# ggsave(plot=p, file=paste(path, 'figures/coal_R0_mobility_withlegend.pdf', sep='/'), height=2.5,width=6)
ggsave(plot=p + theme(legend.position = "none"), file=paste(path, 'figures/coal_R0_mobility.pdf', sep='/'), height=2.5,width=6)


p = ggplot(data=dat.bdsky[which(dat.bdsky$timeframe=="Yakima"),]) + 
  geom_ribbon(aes(x=time, ymin=Ne.lower, ymax=Ne.upper, fill=timeframe), alpha=0.2) +
  geom_ribbon(aes(x=time, ymin=Ne.ll, ymax=Ne.uu, fill=timeframe), alpha=0.8) +
  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"=yakima_col))+
  geom_line(data=subset[which(subset$sub_region_2!="Snohomish County" & subset$sub_region_2!="Skagit County"),], aes(x=date, y=val/30+2.5, linetype=sub_region_2,group=sub_region_2))+
  scale_linetype(name="Mobility in County:")+ 
  scale_y_continuous(sec.axis = sec_axis(~ .*30 - 75, name="workplace mobility")) +
  theme_minimal() +
  # theme(legend.position = "none") + 
  xlab("") +
  ylab("effective reproduction number")+
  coord_cartesian(ylim=c(0,3)) +
    geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)
plot(p)

ggsave(plot=p + theme(legend.position = "none"), file=paste(path, 'figures/bdsky_R0_mobility_yakima.pdf', sep='/'), height=2.5,width=6)
ggsave(plot=p, file=paste(path, 'figures/bdsky_R0_mobility_yakima_legend.pdf', sep='/'), height=2.5,width=4)


p = ggplot(data=growth[which(growth$method=="skygrowth" & growth$timeframe=="Yakima"), ]) + 

  geom_ribbon(aes(x=time, ymin=R.lower, ymax=R.upper, fill=timeframe), alpha=0.2) +
  geom_ribbon(aes(x=time, ymin=R.ll, ymax=R.uu, fill=timeframe), alpha=0.8) +
  scale_x_date(limits=c(as.Date('2020-01-25'), max(mrsi$date)))  +
  scale_fill_manual(name="Dataset:", values = c("WA w/o Yakima"=both_clades_col, "614D"=d_clade_col, "614G"=g_clade_col, "Yakima"=yakima_col))+
  scale_y_continuous(sec.axis = sec_axis(~ .*30 - 75, name="workplace mobility")) +
    scale_color_viridis_d()+
  theme_minimal() +
  geom_line(data=subset[which(subset$sub_region_2!="Snohomish County" & subset$sub_region_2!="Skagit County"),], aes(x=date, y=val/30+2.5, linetype=sub_region_2,group=sub_region_2))+
  scale_linetype(name="Mobility in County:")+ 
  xlab("") +
  ylab("effective reproduction number")+
  coord_cartesian(ylim=c(0,3)) +
    geom_vline(xintercept=vline1, alpha=0.5)+
  geom_vline(xintercept=vline2, alpha=0.5)+
  geom_vline(xintercept=vline3, alpha=0.5)
plot(p)
ggsave(plot=p + theme(legend.position = "none"), file=paste(path, 'figures/coal_R0_mobility_yakima.pdf', sep='/'), height=2.5,width=6)

