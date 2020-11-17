library(ggplot2)
library("coda")
library("colorblindr") 
library(grid)
library(gridExtra)
library(ggtree)

path = "/Users/nmueller/Documents/github/hCoV-19_WA"

lf = list.files(paste(path, "/results/", sep=""), pattern="cluster_assignment.*.tsv")
dat = data.frame()
for (i in seq(1,length(lf))){
  t = read.table(paste(path, "/results/", lf[[i]], sep=""), header=T, sep="\t")
  tmp = strsplit(lf[[i]], split="_")[[1]]
  val = gsub("assignment", "100", gsub(".tsv", "", tmp[[length(tmp)]]))
  
  uni_cl = unique(t$cluster)
  cl_size = rep(0,0)
  for (j in seq(1, length(uni_cl))){
    cl_size[j] = sum(uni_cl[[j]]==t$cluster)
  }
  dat = rbind(dat, data.frame(sub=as.numeric(val)/100, nr = length(unique(t$cluster)), m = mean(cl_size), u = quantile(cl_size, 0.975), median=quantile(cl_size, 0.5)))
}

p = ggplot(dat, aes(x=sub, y=nr)) +
  geom_smooth() +
  geom_point() +
  coord_cartesian(ylim=c(0,max(dat$nr))) +
  theme_minimal()+
  xlab("subsampling proportion of background sequences") +
  ylab("number of clusters")
plot(p)

ggsave(plot=p, filename=paste(path, "/figures/subsampling.png", sep=""), height=3, width=5)


p = ggplot(dat, aes(x=sub, y=m)) +
  geom_point() +
  coord_cartesian(ylim=c(0,100)) +
  theme_minimal()+
  xlab("subsampling proportion of background sequences") +
  ylab("mean cluster size")
plot(p)

ggsave(plot=p, filename=paste(path, "/figures/subsampling_mean.png", sep=""), height=3, width=5)
