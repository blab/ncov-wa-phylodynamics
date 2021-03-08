library(ggplot2)
library("coda")
library("colorblindr") 
library(grid)
library(gridExtra)

path = "/Users/nmueller/Documents/github/hCoV-19_WA"

muts = read.table(paste(path, '/', 'results/subs_and_muts.tsv', sep=''), sep='\t', header = T)
clusters = read.table(paste(path, '/', 'results/cluster_assignment.tsv', sep=''), sep='\t', header = T)

clusters$date = c();
clusters$aa = c();
clusters$mutations = c();

for (i in seq(1,length(clusters$cluster))){
  ind = which(muts$Sequence==as.character(clusters[i, "strain"]))
  if (length(ind==1)){
    clusters[i, "date"] = muts[ind, "Date"]
    clusters[i, "aa"] = muts[ind, "AA"]
    clusters[i, "mutations"] = muts[ind, "Mutations"]
  }
}

clusters$date = as.Date(clusters$date)

uni_cl = unique(clusters$cluster)
dat = data.frame()
for (i in seq(1, length(uni_cl))){
  ind = which(clusters$cluster==uni_cl[[i]])
  min_ind = which.min(clusters[ind, "date"])
  dat = rbind(dat, data.frame(cl=uni_cl[[i]], 
                              size = length(ind), 
                              first = min(clusters[ind, "date"]), 
                              last = max(clusters[ind, "date"]),
                              mAA = mean(clusters[ind, "aa"]),
                              mMutations = mean(clusters[ind, "mutations"]),
                              mAA_min = clusters[ind[min_ind], "aa"],
                              mMutations_min = clusters[ind[min_ind], "mutations"]
                              ))
}


dat = dat[which(dat$first<as.Date("2020-06-01") & dat$first>as.Date("2020-03-01")),]

dat$date_num = as.numeric(dat$first)

sum.vals = data.frame()
for (i in seq(1,1)){
  dat$isSuccess = 1
  dat[which(dat$size <= i), "isSuccess"]=0
  dat$ratio = dat$mAA_min/dat$mMutations_min
  dat$syn = dat$mMutations_min - dat$mAA_min
  
  s = glm(isSuccess ~ mAA_min + syn + date_num + 1, family=binomial, dat)
  # s = glm(size ~ date_num + syn + mAA_min + 1, family=binomial, dat)
  sum.dat = summary(s)
  sum.vals = rbind(sum.vals, data.frame(cutoff = i,
                                        intercept=round(sum.dat$coefficients[1,4],3),
                                        date=round(sum.dat$coefficients[2,4],3),
                                        synonymous=round(sum.dat$coefficients[3,4],3),
                                        non_synonymous=round(sum.dat$coefficients[4,4],3)))
}


p.t=grid.table(sum.vals)


p1=ggplot(dat,aes(x=syn, y=isSuccess)) + 
  geom_count() + 
  geom_smooth(method = lm, color="black") +
  geom_smooth(color="black",linetype=2) +
  xlab("synonymous substitutions first sample") +
  ylab("introduction led to\ndetected local transmission") + 
  coord_cartesian(ylim=c(-0.1,1.1))+
  theme_minimal()

p2=ggplot(dat,aes(x=mAA_min, y=isSuccess)) + 
  geom_count() + 
  geom_smooth(method = lm, color="black") +
  geom_smooth(color="black",linetype=2) +
  xlab("non-synonymous substitutions first sample") +
  ylab("introduction led to\ndetected local transmission") + 
  coord_cartesian(ylim=c(-0.1,1.1))+
  theme_minimal()


p3=ggplot(dat,aes(x=first, y=isSuccess)) + 
  geom_count() + 
  geom_smooth(method = lm, color="black") +
  geom_smooth(color="black",linetype=2) +
  xlab("time first sample") +
  ylab("introduction led to\ndetected local transmission") + 
  coord_cartesian(ylim=c(-0.1,1.1))+
  theme_minimal()

rownames(sum.dat$coefficients)[2] = "non-synonymous\nsubstitutions"
rownames(sum.dat$coefficients)[3] = "synonymous\nsubstitutions"
rownames(sum.dat$coefficients)[4] = "time first sample"
  
library("gridExtra")
library("ggpubr")
p = ggarrange(p2, p1, p3, tableGrob(round(sum.dat$coefficients,4)),
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
plot(p)
ggsave(p, file=paste(path, "/figures/syn_non_syn.png", sep=""), width=10, height=6)

