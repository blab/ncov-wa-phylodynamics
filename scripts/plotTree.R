library(graphics)
library(ape)
library("gridExtra")
library(ggplot2)
# library(treeio)
library(ggtree)
# library("OutbreakTools")

# clear workspace
rm(list = ls())
# clock rate used for time tree
clock.rate = 0.0008

d_clade_col = "#6DC1B3"
g_clade_col = "#F6C445"


# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

tr <- read.tree("../../ncov/results/wa_state/tree.nwk")
tr <- ladderize(tr, right = TRUE)
# read the meta data
meta.wa = read.csv("../data/combined_meta.tsv", header=T, sep="\t")
meta = read.csv("../../ncov/data/example_metadata.tsv", header=T, sep="\t")
clade = read.csv("../../ncov-severity/across-states/strain_to_clade.tsv", header=T, sep="\t")
cluster = read.csv("../results/cluster_assignment.tsv", header=T, sep="\t")

tr$edge.length = tr$edge.length/clock.rate
edges <- tr$edge
node_labs <- tr$node.label
tip_labs <- tr$tip.label

use=FALSE
ind.cluster = c()
for (i in seq(1, length(tr$tip.label))){
  ind = which(meta.wa$strain==tr$tip.label[i])
  if (length(ind)==1){
    ind.cluster = which(cluster$strain==tr$tip.label[i])
    if (length(ind.cluster)==1){
      cl_size = sum(cluster$cluster==cluster[ind.cluster, "cluster"])
      new.tip_labels <- data.frame(id=tr$tip.label[i], study="this study", region="Washington State", county=meta.wa[ind,"location"], date=as.Date(meta.wa[ind,"date"]), cluster=cluster[ind.cluster, "cluster"], cl_size=cl_size)
    }else{
      new.tip_labels <- data.frame(id=tr$tip.label[i], study="this study", region=meta.wa[ind,"region"], county=meta.wa[ind,"location"], date=as.Date(meta.wa[ind,"date"]), cluster=0,cl_size=0)
    }
  }else{
    ind = which(meta$strain==tr$tip.label[i])
    new.tip_labels <- data.frame(id=tr$tip.label[i], study="other study", region=meta[ind,"region"], county=NA, date=as.Date(meta[ind,"date"]), cluster=0,cl_size=0)
  }
  
  ind.clade = which(clade$strain==tr$tip.label[i])
  if (length(ind.cluster)==1){
    new.clade_label = data.frame(id=tr$tip.label[i], clade=clade[ind.clade, "clade"])
  }else{
    new.clade_label = data.frame(id=tr$tip.label[i], clade=clade[ind.clade, "clade"])
  }
  if (i==1){
    tip_labels = new.tip_labels
    clade_labels = new.clade_label
  }else{
    tip_labels = rbind(tip_labels, new.tip_labels)
    clade_labels = rbind(clade_labels, new.clade_label)
  }
}

# in here but not sure if needed
rownames(tip_labels) <- NULL
rownames(clade_labels) <- clade_labels$id
clade_labels$id = c()

# root height in true time, has to be done by hand
root_height = max(tip_labels$date)-max(node.depth.edgelength(tr))*366

labels = c(as.Date("2020-01-01"),
           as.Date("2020-02-01"),
           as.Date("2020-03-01"),
           as.Date("2020-04-01"),
           as.Date("2020-05-01"),
           as.Date("2020-06-01"),
           as.Date("2020-07-01"))
labels_str = c("Jan",
               "Feb",
               "Mar",
               "Apr",
               "May",
               "Jun",
               "Jul")

label_vals = as.numeric(labels-root_height)/366

colors = c("Washington State" ="#268457",
           "Europe"="#14309A",
           "North America"="#BF0B30",
           "Asia"="#3CB3DF",
           "Oceania"="#030303",
           "South America"="#FBE000",
           "Africa"="#F8A93D",
           "grey"="#A5A5A5",
           "D"=d_clade_col,
           "G"=g_clade_col)

cluster_cols = c()

p <- ggtree(tr) %<+% tip_labels + geom_tippoint(aes(color=region), shape=16)+ 
  scale_color_manual(values=colors) + scale_size_manual(values=c(0.5,1)) +aes(color="grey")

p <- p + theme_minimal() +
  theme(panel.grid.major.x=element_line(color="grey20", size=0.3),
        panel.grid.major.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(breaks=label_vals,labels = labels_str)

plot(p)

p=p +theme(legend.position="none")
p2 = gheatmap(p, clade_labels,  color=NA, offset = 0, width=0.025, colnames = F ) + scale_fill_manual(values=colors) + theme(legend.position = "none")
print(p2)

ggsave(plot=p2,"../Figures/Tree.pdf",width=8, height=8)





library(RColorBrewer)
cluster_cols.names = c()
uni_sizes = c(0,unique(tip_labels[!is.na(tip_labels$cl_size),"cl_size"]))
n = length(uni_sizes)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=col_vector[seq(1,length(col_vector),length(col_vector)/n)])

greys = col_vector[seq(1,length(col_vector),length(col_vector)/n)]

for (i in seq(1,length(uni_sizes))){
  cluster_cols[[i+1]] = greys[[i %% (n-1) +1]]
  cluster_cols.names[[i+1]] = as.character(uni_sizes[[i]])
}
names(cluster_cols) = cluster_cols.names


tip_labels$size = tip_labels$cl_size>0

tip_labels$cluster_str = as.character(tip_labels$cl_size)

p <- ggtree(tr,color="grey") %<+% tip_labels + geom_tippoint(aes(color=cluster_str, size=size), shape=16)+ 
  scale_color_manual(values =cluster_cols ) + scale_size_manual(values=c(0,1)) 

p <- p + theme_minimal() +
  theme(panel.grid.major.x=element_line(color="grey20", size=0.3),
        panel.grid.major.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(breaks=label_vals,labels = labels_str)

print(p + theme(legend.position="none"))


ggsave(plot=p,"../Figures/Cluster_tree.pdf",width=12, height=8)


tr_raw <- read.tree("../../ncov/results/wa_state/tree_raw.nwk")

use=FALSE
ind.cluster = c()
for (i in seq(1, length(tr_raw$tip.label))){
  ind.cluster = which(cluster$strain==tr_raw$tip.label[i])
  if (length(ind.cluster)==1){
    cl_size = sum(cluster$cluster==cluster[ind.cluster, "cluster"])
    new.tip_labels <- data.frame(id=tr_raw$tip.label[i], cluster=cluster[ind.cluster, "cluster"], cl_size = cl_size)
  }else{
    new.tip_labels <- data.frame(id=tr_raw$tip.label[i], cluster=0,cl_size = cl_size)
  }
  if (i==1){
    tip_labels_raw = new.tip_labels
  }else{
    tip_labels_raw = rbind(tip_labels_raw, new.tip_labels)
  }
}
tip_labels_raw$size = tip_labels_raw$cluster>0

tip_labels_raw$cluster_str = as.character(tip_labels_raw$cl_size)

p <- ggtree(tr_raw,color="grey", layout='fan', size=0.5) %<+% tip_labels_raw + geom_tippoint(aes(color=cluster_str, size=size), shape=16)+ 
  scale_color_manual(values = cluster_cols) + scale_size_manual(values=c(0,1)) 

print(p + theme(legend.position="none"))
ggsave(plot=p,"../Figures/Cluster_tree_raw.pdf",width=8, height=8)


# plot a figure for the legend
dt=data.frame(color=c("Washington State",
                            "Europe",
                            "North America",
                            "Asia",
                            "Oceania",
                            "South America",
                            "Africa"))
                            
dt$x=1
dt$y=1
p <- ggplot(dt, aes(x=x,y=y,color=color), shape=16, size=10) +
  geom_point()+
  scale_color_manual(values=colors) + theme_minimal()
ggsave(plot=p,"../Figures/tree_legend.pdf",width=2, height=2)

