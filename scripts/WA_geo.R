
library(ggplot2)
library("coda")
library("colorblindr") 
library(grid)
library(gridExtra)
library(ggtree)

path = "/Users/nmueller/Documents/github/hCoV-19_WA"

# define the colors for bd and coal methods
both_clades_col = "#268457"
d_clade_col = "#0072B2"
g_clade_col = "#D55E00"

end_plot = as.Date("2020-07-01")

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



# get the WA sequences by county


# read in the UW VIrology county data
UW.county = read.delim(paste(path, "../ncov-severity/data/WA_df.tsv", sep='/'), header=T, sep="\t")
# read in metadata
cases = read.delim(paste(path, "../ncov/data/example_metadata.tsv", sep='/'), header=T, sep="\t")

# add the location info for the UW sequences
for (i in seq(1,length(UW.county$strain))){
  ind  = which(cases$strain==gsub("hCoV-19/","",UW.county[i,"strain"]))
  if (!is.na(UW.county[i,"county"])){
    cases[ind, "location"] = paste(UW.county[i,"county"], "County")
  }
}

write.table(cases[which(cases$division=='Washington'),], file=paste(path, "data/combined_meta.tsv", sep='/'), sep="\t", quote = FALSE,row.names=FALSE)

cases_tot = cases

cases = cases[which(cases$division=='Washington' & cases$location!=""),]
cases$region=""

for (i in seq(1,length(cases$location))){
  ind = which(region_map$county==as.character(cases[i,"location"]))
  if (length(ind)==1){
    cases[i, "region"] = as.character(region_map[ind, "region"])
  }
}



# get the corresponding strains
strains = read.delim(paste(path, "../ncov-severity/across-states/strain_to_clade.tsv", sep='/'), header=T, sep="\t")
cases$clade = "other"
# map to clade
for (i in seq(1, length(cases$strain))){
  cl = strains[which(strains$strain==as.character(cases[i,"strain"])), "clade"]
  if (length(cl)==1){
    cases[i,"clade"] = as.character(cl)
  }
}
cases$time = as.Date(cases$date)


# Get the number of cases by county


onset_cases = read.table(paste(path,'/data/PUBLIC_Tests_by_Specimen_Collection.csv',sep=""), sep=",", header=T)
onset_cases$Day = as.Date(onset_cases$Day)
testing = data.frame(date=onset_cases$Day, count=onset_cases$Positive, county=onset_cases$County)
testing = testing[-which(testing$county=="Unassigned"),]

for (i in seq(1,length(testing$county))){
  testing[i,"region"] = region_map[which(region_map$county==testing[i,"county"]), "region"]
}

testing$region2 = "WA w/o Yakima"
testing[which(testing$county=="Yakima County"), "region2"] = "Yakima"
testing$region3 = "WA w/o Yakima and King County"
testing[which(testing$county=="Yakima County"), "region3"] = "Yakima and King County"
testing[which(testing$county=="King County"), "region3"] = "Yakima and King County"
testing$region4 = "WA"


for (i in seq(1,length(cases$location))){
  ind = which(region_map$county==as.character(cases[i,"location"]))
  if (length(ind)==1){
    cases[i, "region"] = as.character(region_map[ind, "region"])
  }
}

# combine counts from different counties, but the same region
uni_dates=unique(testing$date)
uni_regions=unique(testing$region)
uni_regions2=unique(testing$region2)
uni_regions3=unique(testing$region3)
uni_regions4=unique(testing$region4)

for (i in seq(4,length(uni_dates)-3)){
  for (j in seq(1,length(uni_regions))){
    nr_cases = 0
    for (k in seq(i-3,i+3)){
      nr_cases = nr_cases+sum(testing[which(testing$date==uni_dates[[k]] & testing$region==uni_regions[[j]]), "count"])
    }
    nr_cases = nr_cases/7
    new.dat = data.frame(date=uni_dates[[i]], region=uni_regions[[j]], count=nr_cases)
    if (i==4 && j==1){
      combined.testing = new.dat
    }else{
      combined.testing = rbind(combined.testing, new.dat)
    }
  }
}


for (i in seq(4,length(uni_dates)-3)){
  for (j in seq(1,length(uni_regions2))){
    nr_cases = 0
    for (k in seq(i-3,i+3)){
      nr_cases = nr_cases+sum(testing[which(testing$date==uni_dates[[k]] & testing$region2==uni_regions2[[j]]), "count"])
    }
    nr_cases = nr_cases/7
    new.dat = data.frame(date=uni_dates[[i]], region=uni_regions2[[j]], count=nr_cases)
    if (i==4 && j==1){
      combined.testing2 = new.dat
    }else{
      combined.testing2 = rbind(combined.testing2, new.dat)
    }
  }
}

for (i in seq(4,length(uni_dates)-3)){
  for (j in seq(1,length(uni_regions3))){
    nr_cases = 0
    for (k in seq(i-3,i+3)){
      nr_cases = nr_cases+sum(testing[which(testing$date==uni_dates[[k]] & testing$region3==uni_regions3[[j]]), "count"])
    }
    nr_cases = nr_cases/7
    new.dat = data.frame(date=uni_dates[[i]], region=uni_regions3[[j]], count=nr_cases)
    if (i==4 && j==1){
      combined.testing3 = new.dat
    }else{
      combined.testing3 = rbind(combined.testing3, new.dat)
    }
  }
}

for (i in seq(4,length(uni_dates)-3)){
  for (j in seq(1,length(uni_regions4))){
    nr_cases = 0
    for (k in seq(i-3,i+3)){
      nr_cases = nr_cases+sum(testing[which(testing$date==uni_dates[[k]] & testing$region4==uni_regions4[[j]]), "count"])
    }
    nr_cases = nr_cases/7
    new.dat = data.frame(date=uni_dates[[i]], region=uni_regions3[[j]], count=nr_cases)
    if (i==4 && j==1){
      combined.testing4 = new.dat
    }else{
      combined.testing4 = rbind(combined.testing4, new.dat)
    }
  }
}




p <- ggplot(combined.testing) +
  geom_histogram(aes(x=date,y=count), stat="identity", position="stack") +
  facet_wrap(.~region,ncol=3, scales="free_y")+
  scale_x_date() + 
  theme_minimal()
  # scale_y_log10()

plot(p)



# Compare clades over time and county


transform=1
for (i in seq(1,length(uni_regions))){
  p <- ggplot(data=cases[which(cases$region==uni_regions[[i]]),])+
    geom_histogram(aes(x=time, group=clade, fill=clade), binwidth=7)+
    scale_fill_manual(values=c("D"=d_clade_col, "G"=g_clade_col, "other"=both_clades_col)) +
    scale_x_date(limits=c(min(cases$time), end_plot))+
    geom_line(data=combined.testing[which(combined.testing$region==uni_regions[[i]]),], aes(x=date, y=count/transform), stat="identity", position="stack", color=both_clades_col)+
    xlab("")+
scale_y_continuous(name="daily positive tests ", sec.axis = sec_axis(~ 1/transform*., name = "weekly nr sequences")) +
    theme( axis.title.y.right = element_text( angle = 90))+
coord_cartesian(ylim=c(0.01,500)) +

    theme_minimal() +
    theme(legend.position = "none")
  plot(p)
  ggsave(p, file=paste(path, "/figures/", uni_regions[[i]], '.pdf', sep=""), width=4, height=3)
}

transform=1

p <- ggplot(data=cases[which(cases$location!="Yakima County"),])+
  geom_histogram(aes(x=time, group=clade, fill=clade), binwidth=7)+
  scale_fill_manual(values=c("D"=d_clade_col, "G"=g_clade_col, "other"=both_clades_col)) +
  scale_x_date(limits=c(min(cases$time), end_plot))+
  geom_line(data=combined.testing2[which(combined.testing2$region!="Yakima"),], aes(x=date, y=count/transform), stat="identity", position="stack", color=both_clades_col)+
  xlab("")+
scale_y_continuous(name="daily positive tests ", sec.axis = sec_axis(~ 1/transform*., name = "weekly nr sequences")) +
    theme( axis.title.y.right = element_text( angle = 90))+
coord_cartesian(ylim=c(0.01,500)) +

  theme_minimal() +
  theme(legend.position = "none")
plot(p)
ggsave(p, file=paste(path, "/figures/Wa_wo_yakima.pdf", sep=""), width=4, height=3)

p <- ggplot(data=cases[which(cases$location=="Yakima County"),])+
  geom_histogram(aes(x=time, group=clade, fill=clade), binwidth=7)+
  scale_fill_manual(values=c("D"=d_clade_col, "G"=g_clade_col, "other"=both_clades_col)) +
  scale_x_date(limits=c(min(cases$time), end_plot))+
  geom_line(data=combined.testing2[which(combined.testing2$region=="Yakima"),], aes(x=date, y=count/transform), stat="identity", position="stack", color=both_clades_col)+
  xlab("")+
scale_y_continuous(name="daily positive tests ", sec.axis = sec_axis(~ 1/transform*., name = "weekly nr sequences")) +
    theme( axis.title.y.right = element_text( angle = 90))+
coord_cartesian(ylim=c(0.01,500)) +

  theme_minimal() +
  theme(legend.position = "none")
plot(p)
ggsave(p, file=paste(path, "/figures/yakima.pdf", sep=""), width=4, height=3)


p <- ggplot(data=cases)+
  geom_histogram(aes(x=time, group=clade, fill=clade), binwidth=7)+
  scale_fill_manual(values=c("D"=d_clade_col, "G"=g_clade_col, "other"=both_clades_col)) +
  scale_x_date(limits=c(min(cases$time), end_plot))+
  geom_line(data=combined.testing4, aes(x=date, y=count/transform), stat="identity", position="stack", color=both_clades_col)+
  xlab("")+
scale_y_continuous(name="daily positive tests ", sec.axis = sec_axis(~ 1/transform*., name = "weekly nr sequences")) +
    theme( axis.title.y.right = element_text( angle = 90))+
coord_cartesian(ylim=c(0.01,500)) +

  theme_minimal() +
  theme(legend.position = "none")
plot(p)
ggsave(p, file=paste(path, "/figures/Wa.pdf", sep=""), width=4, height=3)



require(rvest)
system(paste("curl -fsSL https://www.doh.wa.gov/Emergencies/NovelCoronavirusOutbreak2020COVID19/DataDashboard >> ", path, "/results/sampling.json", sep=""))
tmp = readLines(paste(path, "/results/sampling.json", sep="")) # line with a warning
first = T
for (i in seq(1,length(tmp))){
  if (startsWith(tmp[[i]],"<script type=\"application/json\"")){
    tmp2 = strsplit(tmp[[i]], split='\\[')[[1]]
    for (j in seq(1,length(tmp2))){
      if (grepl("Date:",tmp2[[j]])){
        tmp3 = strsplit(tmp2[[j]], split=",")[[1]]
        for (k in seq(1,length(tmp3))){
          tmp4 = strsplit(gsub("\"","",tmp3[[k]]), split="<br />")[[1]]
          if (length(tmp4)==3){
            tmp4 = strsplit(tmp4, split="\\:\\s+")
            new.dat = data.frame(date = as.Date(tmp4[[1]][[2]]), count = as.numeric(tmp4[[2]][[2]]), detection=gsub(']','',tmp4[[3]][[2]]))
            if (first){
              testing = new.dat
              first = F
            }else{
              testing = rbind(testing, new.dat)
            }
          }
        }
      }
    }
    break;
  }
}
library(usmap)
p<-plot_usmap("county",include = c("WA"))
ggsave(p, file=paste(path, "figures/WA_map.pdf", sep="/"),width=9,useDingbats=FALSE )

mobility = read.delim(paste(path,'/data/applemobilitytrends-2020-05-19.csv',sep=""), sep=",", header=T,na.strings=c("",".","NA"))


mobility = mobility[which(mobility$geo_type=="county" & mobility$sub.region=="Washington"),]

header.labs = labels(mobility)[[2]]

first=T
for (i in seq(1,length(mobility$geo_type))){
  if (grepl("County", mobility[i,"region"])){
    print(mobility[i,"region"])
    for (j in seq(1,length(header.labs))){
      if (startsWith(header.labs[[j]], "X2020")){
        new.dat = data.frame(time=as.Date(gsub("\\.","-",gsub("X","", header.labs[[j]]))), val=as.numeric(mobility[i,j]), mode=mobility[i,"transportation_type"], county=mobility[i,"region"])
        if(first){
          mob.dat = new.dat
          first=F
        }else{
          mob.dat = rbind(mob.dat, new.dat)
        }
      }
    }
  }
}
king_county = mob.dat[which(mob.dat$county=="King County"),]
king_county$county = c()

p = ggplot(mob.dat, aes(x=time, y=val))+
  geom_line()+
  geom_smooth(aes(group=mode))+
  geom_smooth(data=king_county, aes(x=time, y=val), color="red")+
  facet_wrap(.~county, ncol=4, scales = "free_y") +
  theme_minimal() +
  scale_x_date(limits=c(as.Date('2020-01-01'), as.Date('2020-05-19'))) +
  ylab("mobility (driving)") + 
  xlab("")
plot(p)
  
ggsave(p, file=paste(path, "/figures/mobility.pdf", sep=""), width=15, height=12)

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
subset = mobility[which(mobility$sub_region_2=="Skagit County" | mobility$region=="King" | mobility$region=="Pierce" | mobility$sub_region_2=="Snohomish County"),]

# define time period for which to look into correlation
start = as.Date('2020-02-29')
end = as.Date('2020-03-20')

mob.val = labels(subset)[[2]]

first = T

for (l in seq(10,length(mob.val)-2)){
  print(mob.val[[l]])
  subset$val = NA
  
  for (j in seq(4,length(subset$country_region_code)-3)){
    subset[j,"val"]=mean(subset[seq(j-3,j+3),mob.val[[l]]])
  }
  
  king_county = subset[which(subset$sub_region_2=="King County"),]
  king_county$sub_region_2 = c()

}
subset$bestlag=9
subset[which(subset$sub_region_2=="Snohomish County"),"bestlag"] = 6

p.trend = ggplot(subset[which(subset$sub_region_2!="King County" & subset$sub_region_2!="Yakima County"),], aes(x=date, y=val))+
  geom_line(aes(color="true values"))+
  geom_line(aes(x=date-bestlag, color="shifted values"))+
  # geom_smooth()+
  geom_line(data=king_county, aes(x=date, y=val, color="King County"))+
  facet_wrap(.~sub_region_2, ncol=3, scales = "free_y") +
  theme_minimal() +
  scale_color_discrete(name="")+
  scale_x_date(limits=c(start,end+14))+
  geom_text(data=subset[which(subset$sub_region_2!="King County" & subset$sub_region_2!="Yakima County"),], aes(x = end+5, y = -10,label=paste("shifted by\n",bestlag,"days")))+
  ylab("workplace mobility") +
  theme(legend.position = "top")+
  xlab("")
plot(p.trend)
ggsave(p.trend, file=paste(path, "/figures/mobility_trend.png", sep=""), width=9, height=4)


lats_longs = read.table(file=paste(path, "/data/lat_longs.tsv", sep=""), sep="\t")
lats_longs$V1 = as.character(lats_longs$V1)

if (require("maps")) {
  county <- map_data("county")
  states <- map_data("state")
}


connections = read.table(file=paste(path, "/results/connections.tsv", sep=""), sep="\t", header=T)
for (i in seq(1,length(connections$from_county))){
  ind = which(lats_longs$V1==connections[i,"from_county"])
  connections[i,"from_lat"] = lats_longs[ind,"V2"]
  connections[i,"from_long"] = lats_longs[ind,"V3"]
  ind = which(lats_longs$V1==connections[i,"to_county"])
  connections[i,"to_lat"] = lats_longs[ind,"V2"]
  connections[i,"to_long"] = lats_longs[ind,"V3"]
}


p.pairs = ggplot(connections[which(connections$pairwise==0),])+
  geom_polygon(data = county, aes(long, lat, group = group), fill="grey", color="black", alpha=0.1, size=0.15) +
  geom_polygon(data = states, aes(long, lat, group = group), fill=NA, color="black", size=0.5) +
  geom_segment(aes(x=from_long,xend=to_long,y=from_lat, yend=to_lat, size=pairs)) +
  scale_size(range=c(0,2), limits=c(0,max(connections[which(connections$prob>0.95 & connections$pairwise==0), "pairs"])))+
  coord_map("albers",  at0 = 45.5, lat1 = 29.5, xlim=c(-124,-117), ylim=c(46,49)) +
  theme_bw( ) +
  theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

plot(p.pairs)

ggsave(p.pairs, file=paste(path, "/figures/county_pairs.pdf", sep=""), width=8, height=4)


