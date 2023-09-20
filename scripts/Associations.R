setwd("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/2_Vocal interaction networks/Vocal interactions")

#install.packages(c("igraph", "bipartite", "asnipe", "assortnet", "ggplot2", "ggmap", "rnetcarto", "ecodist", "igraphdata", "statnet", "RColorBrewer", "tidyverse", "geosphere", "mapproj"))

library(tidyr) #or you can just load the whole tidyverse with library(tidyverse)
library(asnipe)
library(igraph)
library(RColorBrewer)
library(stringr)
library(sna)
library(ggplot2)

#Using https://dshizuka.github.io/networkanalysis/example_make_sparrow_network.html


#Input data
batdat21=read.csv('associations_2021.csv', header = TRUE, sep = ";", colClasses = c("numeric", "character", "Date", "character", "character", "character","character","character","character","character","character","character","character"))
batdat22=read.csv('associations_2022.csv', header = TRUE, sep = ";", colClasses = c("numeric", "character", "Date", "character", "character", "character","character","character","character","character","character"))

head(batdat21)
head(batdat22)


#Add attribute data
attrib_21=read.csv('supp_21.csv', header = TRUE, sep = ";", colClasses = c("character", "character"))
attrib_21
attrib_22=read.csv('supp_22.csv', header = TRUE, sep = ";", colClasses = c("character", "character"))
attrib_22


#Prepare association data 2021
batcols=grep("Bat",colnames(batdat21))
batcols
batdat21[,batcols]
gather(batdat21[,batcols])
bat.ids=unique(gather(batdat21[,batcols])$value)
bat.ids
bat.ids=bat.ids[is.na(bat.ids)==F]
bat.ids
bat.ids_df_21 <- data.frame(bat.ids)
csv_path  = 'bats_21.csv'
write.csv(bat.ids_df_21, csv_path, row.names = FALSE)

m1_21=apply(batdat21[,batcols], 1, function(x) match(bat.ids,x))
m1_21

m1_21[is.na(m1_21)]=0
m1_21[m1_21>0]=1
m1_21

rownames(m1_21)=bat.ids #rows are bat ids
colnames(m1_21)=paste('group', 1:ncol(m1_21), sep="_") #columns are group IDs (just "group_#")
m1_21

rowSums(m1_21)
set.seed(2)
png(file="C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/2_Vocal interaction networks/Vocal interactions/Histogram captures 2021.png",
    width=600, height=400)
par(oma=c(1,1,1,1)) # all sides have 3 lines of space
par(mar=c(5,4,4,2) + 0.8)
hist (rowSums(m1_21), main="", xlab="Number of observations", ylab="Frequency", cex.lab = 2.5, cex.axis = 1.5)
dev.off()

average_n_21 <- mean(rowSums(m1_21)) #average number of times a single bat was observed
average_n_21

sd_n_21 <- sd(rowSums(m1_21)) #sd number of times a single bat was observed
sd_n_21

average_gs_21 <- mean(colSums(m1_21)) #average group size
average_gs_21

sd_gs_21 <- sd(colSums(m1_21)) #sd group size
sd_gs_21


#Prepare association data 2022
batcols=grep("Bat",colnames(batdat22))
batcols
batdat22[,batcols]
gather(batdat22[,batcols])
bat.ids=unique(gather(batdat22[,batcols])$value)
bat.ids
bat.ids=bat.ids[is.na(bat.ids)==F]
bat.ids
bat.ids_df_22 <- data.frame(bat.ids)
csv_path  = 'bats_22.csv'
write.csv(bat.ids_df_22, csv_path, row.names = FALSE)

m1_22=apply(batdat22[,batcols], 1, function(x) match(bat.ids,x))
m1_22

m1_22[is.na(m1_22)]=0
m1_22[m1_22>0]=1
m1_22

rownames(m1_22)=bat.ids #rows are bat ids
colnames(m1_22)=paste('group', 1:ncol(m1_22), sep="_") #columns are group IDs (just "group_#")
m1_22

rowSums(m1_22)
set.seed(2)
png(file="C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/2_Vocal interaction networks/Vocal interactions/Histogram captures 2022.png",
    width=600, height=400)
par(oma=c(1,1,1,1)) # all sides have 3 lines of space
par(mar=c(5,4,4,2) + 0.8)
hist (rowSums(m1_22), main="", xlab="Number of observations", ylab="Frequency", cex.lab = 2.5, cex.axis = 1.5)
dev.off()

average_n_22 <- mean(rowSums(m1_22)) #average number of times a single bat was observed
average_n_22

sd_n_22 <- sd(rowSums(m1_22)) #sd number of times a single bat was observed
sd_n_22

average_gs_22 <- mean(colSums(m1_22)) #average group size
average_gs_22

sd_gs_22 <- sd(colSums(m1_22)) #average group size
sd_gs_22

#Create adjacency data for associations (2021)
adj_21=get_network(t(m1_21), data_format="GBI", association_index = "SRI") # the adjacency matrix
g_21=graph_from_adjacency_matrix(adj_21, "undirected", weighted=T) #the igraph object
V(g_21)$sex <- factor(attrib_21[match(V(g_21)$name, attrib_21$bat_id), "sex"])

mismatched_names <- setdiff(V(g_21)$name, attrib_21$bat_id) #check for issues


#Create adjacency data for associations (2022)
adj_22=get_network(t(m1_22), data_format="GBI", association_index = "SRI") # the adjacency matrix
g_22=graph_from_adjacency_matrix(adj_22, "undirected", weighted=T) #the igraph object
V(g_22)$sex <- factor(attrib_22[match(V(g_22)$name, attrib_22$bat_id), "sex"])

mismatched_names <- setdiff(V(g_22)$name, attrib_22$bat_id) #check for issues


#Select individuals for which more than 10 observations were made (not incorporated into adjacency data for associations, above)
#m2_21=m1_21[which(rowSums(m1_21)>10),]
#adj_21=get_network(t(m2_21), data_format="GBI","SRI")
#write.csv(adj_21, "assoc_21.csv", row.names=TRUE)
#rowSums(m2_21)

#m2_22=m1_22[which(rowSums(m1_22)>10),]
#adj_22=get_network(t(m2_22), data_format="GBI","SRI")
#rowSums(m2_22)
#write.csv(adj_22, "assoc_22.csv", row.names=TRUE)
#rowSums(m2_22)


#Create communities (2021)
com_21=fastgreedy.community(g_21) #community detection method
com_21[[11]]

# Assign community names
mapping <- c(
  11,  # Old membership 1 should be corrected to 11
  8,   # Old membership 2 should be corrected to 8
  1,   # Old membership 3 should be corrected to 1
  2,   # Old membership 4 should be corrected to 2
  6,   # Old membership 5 should be corrected to 6
  9,   # Old membership 6 should be corrected to 9
  4,   # Old membership 7 should be corrected to 4
  10,  # Old membership 8 should be corrected to 10
  7,   # Old membership 9 should be corrected to 7
  5,   # Old membership 10 should be corrected to 5
  3    # Old membership 11 should be corrected to 3
)

corrected_memberships <- mapping[com_21$membership]

com_21$membership <- corrected_memberships

community_assignments <- com_21$membership

color_mapping <- c(
  "female" = "slateblue",
  "male" = "gold",
  "unknown" = "gray"
)
node_sex <- V(g_21)$sex
node_colors <- sapply(node_sex, function(sex) {
  color_mapping[sex]
})
V(g_21)$color <- node_colors

#Prepare to save figure 
set.seed(2)
png(file="C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/2_Vocal interaction networks/Vocal interactions/Association networks 2021.png",
    width=600, height=600)

# Iterate through com_21 and assign community names to nodes in g_21
plot(
  g_21, 
  vertex.color = node_colors,  # Use the calculated node colors
  vertex.size = 10,  # Adjust the size of the nodes as needed
  vertex.label = community_assignments, # Use the extracted community assignments as labels
  edge.width=E(g_21)$weight*10
)

title(main = "2021", cex.main = 2)

legend("topleft", legend = names(color_mapping), fill = color_mapping, cex = 1.5)

dev.off()


#Create communities (2022) 
com_22=fastgreedy.community(g_22) #community detection method
com_22[[12]]

# Assign community names
mapping <- c(
  3,  
  4, 
  1,  
  9, 
  10, 
  11,  
  7, 
  5, 
  6, 
  2, 
  12,
  8
)

corrected_memberships <- mapping[com_22$membership]

com_22$membership <- corrected_memberships

community_assignments <- com_22$membership

color_mapping <- c(
  "female" = "slateblue",
  "male" = "gold",
  "unknown" = "gray"
)
node_sex <- V(g_22)$sex
node_colors <- sapply(node_sex, function(sex) {
  color_mapping[sex]
})
V(g_22)$color <- node_colors

#Prepare to save figure 
set.seed(2)
png(file="C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/2_Vocal interaction networks/Vocal interactions/Association networks 2022.png",
    width=600, height=600)

# Iterate through com_21 and assign community names to nodes in g_21
plot(
  g_22, 
  vertex.color = node_colors,  # Use the calculated node colors
  vertex.size = 10,  # Adjust the size of the nodes as needed
  vertex.label = community_assignments, # Use the extracted community assignments as labels
  edge.width=E(g_22)$weight*10
)

title(main = "2022", cex.main = 2)

dev.off()

#Estimate modulatity
modularity(com_21)
modularity(com_22)


#Estimate community measures
g.cluster=transitivity(g_21, "global")
g.cluster
l.cluster=transitivity(g_21, "local")
l.cluster
av.l.cluster=transitivity(g_21, "localaverage")
av.l.cluster

g.cluster=transitivity(g_22, "global")
g.cluster
l.cluster=transitivity(g_22, "local")
l.cluster
av.l.cluster=transitivity(g_22, "localaverage")
av.l.cluster

#Estimate community size
community_sizes_21 <- table(com_21$membership)
community_sizes_22 <- table(com_22$membership)

# Calculate the average community size and sd
average_community_size_21 <- mean(community_sizes_21)
average_community_size_22 <- mean(community_sizes_22)

average_community_size_21
average_community_size_22

sd_community_size_21 <- sd(community_sizes_21)
sd_community_size_22 <- sd(community_sizes_22)

sd_community_size_21
sd_community_size_22


#___________________________________________________________________________________________________#

#Vocal interactions (general trends)

setwd("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/2_Vocal interaction networks/Vocal interactions")

calls21=read.csv('vocal_interactions_2021.csv', header = TRUE, sep = ";", colClasses = c("numeric", "character", "character", "character", "character","character","numeric","numeric","numeric"))
calls22=read.csv('vocal_interactions_2022.csv', header = TRUE, sep = ";", colClasses = c("numeric", "character", "character", "character", "character","character","numeric","numeric","numeric"))

head(calls21)
head(calls22)

summary(calls21)
summary(calls22)

length(unique(calls21$Dyad))
length(unique(calls22$Dyad))

summary(subset(calls21, flight_time == 300))
summary(subset(calls22, flight_time == 300))

#Subset data for which we were able to record for the whole 300 s
calls21_300 <- subset(calls21, flight_time == 300)
calls22_300 <- subset(calls22, flight_time == 300)

sd(calls21_300$n_inquiry)
sd(calls21_300$n_response)

sd(calls22_300$n_inquiry, na.rm = TRUE)
sd(calls22_300$n_response, na.rm = TRUE)

