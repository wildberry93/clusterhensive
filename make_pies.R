require(ggplot2)
require(RColorBrewer)
require(cowplot)

draw.barplot <- function(df){
  df$color[df$V1 == "Bacteria"] <- myColors[[4]]
  df$color[df$V1 == "Archaea"] <- myColors[[2]]
  df$color[df$V1 == "Eukaryota"] <- myColors[[3]]
  df$color[df$V1 == "Viruses"] <- myColors[[1]]
  bp <- ggplot(df, aes(x="", y=V2))+
    geom_bar(width = 1, stat = "identity", fill=df$color)+
    coord_flip()+
    theme_nothing()+
    scale_y_discrete(expand=c(0,0))
  return(bp)
}

myColors <- brewer.pal(4,"Set1")

filelist = list.files(pattern = ".*.txt")
names(filelist) <- filelist
datalist = lapply(filelist, function(x)read.table(x, header=F)) 

plots <- lapply(datalist, draw.barplot)

for(a in 1:length(plots)){
  name = paste(strsplit(names(plots[a]),"\\.")[[1]],".png",sep = "")[[1]]
  png(name,width = 100, height = 5)
  print(plots[a])
  
  dev.off()
}