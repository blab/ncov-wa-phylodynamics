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
  clusters[i, "date"] = muts[ind, "Date"]
  clusters[i, "aa"] = muts[ind, "AA"]
  clusters[i, "mutations"] = muts[ind, "Mutations"]
}

clusters$date = as.Date(clusters$date)

uni_cl = unique(clusters$cluster)
dat = data.frame()
for (i in seq(1, length(uni_cl))){
  ind = which(clusters$cluster==uni_cl[[i]])
  dat = rbind(dat, data.frame(cl=uni_cl[[i]], 
                              size = length(ind), 
                              first = min(clusters[ind, "date"]), 
                              last = max(clusters[ind, "date"]),
                              mAA = mean(clusters[ind, "aa"]),
                              mMutations = mean(clusters[ind, "mutations"])))
}


ggplot(dat, aes(x=mAA, y=size)) + 
  geom_point() + 
  scale_y_log10() + 
  geom_smooth(method='lm', formula= y~x)

ggplot(dat,aes(x=mAA, y=as.numeric(last-first))) + 
  geom_point() + 
  geom_smooth(method='lm', formula= y~x)


ggplot(dat) + 
  geom_point(aes(x=first, y=last-first)) 



