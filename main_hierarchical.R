#############################
#############################
#  needed additional functions:
#  This function checks if a set objset of objects (given as a 0-1 vector) is
#  a minimal generator in the sense that every proper subset of objset generates
#  a smaller extent than objset
#  Note that (assuming possible duplicates for all objects) for hierarchical
#  nominal scaling the minimal generators are identical to the ufg premises
#  which are exactly the one-element premises and
#  the two-element premises of objects with attributes that are not exactly
#  identical
#############################
#############################

check_if_is_minimal_generator <- function(objset,context){
  if(all(objset==0)){return(FALSE)}
  extent <- ddandrda:::operator_closure_obj_input(objset,context)
  indexs <- which(objset==1)
  for(k in indexs){
    reduced_objset <- objset
	  reduced_objset[k] <- 0
	if(all(ddandrda:::operator_closure_obj_input(
	  reduced_objset,context)==extent)){return(FALSE)}
  }
  return(TRUE)
}




# This function generates from hierarchical nominal data that are given in the
# format like the ISCO-08 format of the Allbus data set (ZA5280_v2-0-1.sav) the
# corresponding formal context
compute_hierarchical_scaling_vec <- function(values){
   # computes hierarchical scaling
   data_values <- rep("",length(values))
   for(k in (1:length(data_values))){
     if(nchar(values[k])==3){data_values[k] <- paste("0",as.character(values[k]),sep="") }
     if(nchar(values[k])==4){data_values[k] <- as.character(values[k]) }
     if(nchar(values[k]) <=2){print("ERROR");return(NULL)}
   }
   m <- length(data_values)
   values <- rep("",m*4)
   t <- 1
   for(l in c(1,2,3,4)){#c(1,10,100,1000)){
     for(k in (1:m)){
       values[t] <- substr(as.character(data_values[k]),1,l)#data_values[k]%/%l
       t <- t+1
     }
   }
   unique_values <- unique(values)
   print(unique_values)
   X <- NULL
   for(value in unique_values){
     temp <- rep(0,m)
     L <-nchar(value)
     for(k in (1:m)){
         temp[k] <- substr(as.character(data_values[k]),1,L)== value
     }
     X <- cbind(X,temp)
    }
    colnames(X) <- unique_values
return(X)}



#####
# This function computes the contributing values for the ufg depth given by one element ufg premises
ufg_1_depth <- function(context,weights=rep(1,nrow(context))){

  # computes all ufg premises of cardinality 1 and the corresponding depth values


  m <- nrow(context)
  ans <- rep(0,m)
  S <- 0
  for(k in (1:(m))){
      print(k)
      extent <- rep(0,m)
      extent[k] <- 1
      if(check_if_is_minimal_generator(extent,context)){
        extent <- ddandrda:::operator_closure_obj_input(extent,context)
        ans[which(extent==1)]=ans[which(extent==1)]+weights[k]
        S <- S+weights[k]
      }
    }

  return(list(depths = ans/S,number_of_ufgs=S))
}

#########
# This function computes the contributing values for the ufg depth given by two element ufg premises
ufg_2_depth <- function(context,weights=rep(1,nrow(context))){

  # computes all ufg premises of cardinality 2 and the corresponding depth values, the corresponding intents are also computed
  m <- nrow(context)
  ans <- rep(0,m)
  t <- 1
  S <- 0
  for(k in (1:(m-1))){
    print(k)
    for(l in((k+1):m)){
	    extent <- rep(0,m)
      extent[c(k,l)] <- 1
      if(check_if_is_minimal_generator(extent,context)){
        extent <- ddandrda:::operator_closure_obj_input(extent,context)
        ans[which(extent==1)]=ans[which(extent==1)]+weights[k]*weights[l]
	      S <- S+weights[k]*weights[l]
        t <- t+1
      }
    }
  }
return(list(depths=ans/S,number_of_ufgs=S) )
}







#######################################
## Actual Script
#######################################

# Allbus Analysis
library(foreign)
library(plyr)
library(phangorn)
library(ggplot2)
library(ggtree)
# read the Allbus 2021 data ( https://search.gesis.org/research_data/ZA5280 )
dat <- read.spss("ZA5280_v2-0-1.sav",use.value.labels=FALSE)

# select only persons, for which the ISCO-08 status (occupational status according
# to ISCO-08) is available (dat$isco08>0)
indexs <- which(dat$isco08>0)# & dat$sex %in% c(1,1) &dat$eastwest==1& dat$age%in% (20:30))# &  dat$eastwest==1 & dat$sex==2)
# data points coded with 4 digits from 0-9
x <- dat$isco08[indexs]
# different weights are used for east and west because of an oversampling of
# east germany
weights <- dat$wghthew[indexs]
context <- compute_hierarchical_scaling_vec(x)
rownames(context) <- x
# make a weighted data table (first column contains the names of the
# hierarchical categories, y_weighted are the reweighted counts)
weighted_repr <- oofos:::get_weighted_representation(x=cbind(x,context),y=weights)
weighted_context <- weighted_repr$x[,-1]
rownames(weighted_context) <- weighted_repr$x[,1]
y_weighted <- weighted_repr$y_weighted
# xw contains the categories for the weighted context
xw <- weighted_repr$x_weighted[,1]



# compute contributions of 1-element and 2-element ufg premises, as well as the
# final vector that contains ufg-depths
ufg_1 <- ufg_1_depth(weighted_context,weights=y_weighted)
ufg_2 <- ufg_2_depth(weighted_context,weights=y_weighted)
ufg   <- ufg_1$depths + ufg_2$depths








# data points with highest depths:
i1 <- which(ufg_1$depth==max(ufg_1$depths))
xw[i1]
# [1] 4110
#deepest point w.r.t. ufg_1: 4110: Allgemeine Buerrokraefte

 i2 <- which(ufg_2$depths==max(ufg_2$depths))
 xw[i2]
 # [1] 3341 3343 3342 3344
 # deepest points w.r.t. ufg_2 (all subcategories of category 334: Sekretariatsfachkraefte)
 #
 # 3341: Sekretariatsleiter
 # 3342: Sekretariatsfachkr?fte im juristischen Bereich
 # 3343: Sekretariatsfachkr?fte in Verwaltung und Gesch?ftsleitung
 # 3344: Sekretariatsfachkr?fte im Gesundheitswesen


 i <- which(ufg==max(ufg))
 xw[i]
 # [1] 3221
 #deepest point w.r.t. ufg depth: 3221: Nicht akademische Krankenpflegefachkraefte

 max(ufg)
 # maximal depth
 # [1] 0.9262753


 # Further analysis
 length(unique(ufg))
 # 285 unique depth values



length(unique(ddandrda:::compute_quasiconcave_hull(ufg,weighted_context)))
# [1] 3 contour sets sets of D^{qc}



# modus on the finest hierarchical level
(table(x))[which.max(table(x))]

# 4110
#  112
# Kategorialer Modus: 4110: Allgemeine Buerokraefte (of course identical to l. 164)

# Top down approach: Iterative modus starting from the most coarse level followed by successively looking at more fine-grained levels
j <- seq_len(nrow(context))
for(k in (1:4)){
  o <- order(colSums(context[j,]),decreasing=TRUE)
  j <- which(context[,o[k]]==1)
}
colnames(context)[o[4]]

# [1] "3343" Sekretariatsfachkraefte in Verwaltung und Geschaeftsleitung / Administrative and Executive Secretaries


# smallest depthe:

which(ufg==min(ufg))
# [1] 281
xw[281]
# [1] 6210
# datapoint with smallest depth: 6210: Forstarbeitskraefte und verwandte Berufe

# computation of the (corresponding intents of the) upper level sets
depth_values=sort(unique(ufg),decreasing=TRUE)

for( threshold in depth_values){
  indexs <- which(ufg>= threshold)
  extent <- rep(0,length(xw))
  extent[indexs] <- 1
  intent <- ddandrda:::calculate_psi(extent, weighted_context)
  print(colnames(weighted_context)[which(intent==1)])
  if(all(intent==0)){break}
  }

## generalized Tukey depth:

D_tukey <- ddandrda::compute_tukeys_depth(context,context)
table(D_tukey)
# D_tukey
# 0.708518518518517 0.758888888888886
# 1913               787
#######
# graphic still to be improved
library(ggsci)
col_indexs <- which(colnames(weighted_context)%in%c(3,4))
row_indexs <- which(weighted_context[,col_indexs[1]]==1 | weighted_context[,col_indexs[2]]==1)
Dist <- as.matrix(dist(weighted_context[row_indexs,],method="manhattan"))

tree <- upgma(Dist)
p <- ggtree(tree,size=1) +  geom_tiplab(align=TRUE, linesize=2)
d1 <- data.frame(id=tree$tip.label, val=y_weighted[row_indexs])
p2 <- facet_plot(p, panel='bar', data=d1,  geom=geom_segment,aes(x=0, xend=val, y=y, yend=y), size=3, color='blue4')
p2 + theme_tree2()

f <- function(x){log(1+0.1*x)}
### other graphic
x_min <- 3000
x_max <- 5000
indexs <- which(x>=x_min & x <=x_max)
df <- data.frame(x=x[indexs])
lwd1 <- 1
lwd2 <- 0.8
cols=(pal_bmj(palette = c("default"), alpha = 1))(n=6)
CC=1/1056
pdf("allbusplot1.pdf")
ggplot(df, aes(x=x))  + theme(axis.text.x = element_text(size = 10, angle = 90)) + scale_x_continuous(breaks = seq(3000, 4999, by = 100)) + geom_vline(alpha=0.5,xintercept=3343,col=cols[4],lwd=lwd2) +geom_vline(xintercept=3221,alpha=0.5,col=cols[5],lwd=lwd2) +geom_vline(xintercept=4110,col=cols[6],lwd=lwd2,alpha=0.5) + geom_histogram(aes(y = after_stat(f(count))),breaks=seq(x_min,x_max,1000),col="black",lwd=lwd1,alpha=0.3)  +labs(y="counts (log scale: y=ln(1+0.1*counts))",x="ISCO-08 classification of occupation") +
geom_histogram(aes(y = after_stat(f(count))),breaks=seq(x_min,x_max,100),alpha=0.2,col=cols[1],lwd=lwd1) +
geom_histogram(aes(y = after_stat(f(count))),breaks=seq(x_min,x_max,10),col=cols[2],alpha=0.5,lwd=lwd1) +
geom_histogram(aes(y = after_stat(f(count))),breaks=seq(x_min,x_max,1),col=cols[3],alpha=1,lwd=lwd1)
dev.off()

# further graphic

x_min <- 3200
x_max <- 3400
indexs <- which(x>=x_min & x <=x_max)
df <- data.frame(x=x[indexs])
lwd1 <- 1
lwd2 <- 0.8
cols=(pal_bmj(palette = c("default"), alpha = 1))(n=6)
CC=.1
pdf("allbusplot2.pdf")


ggplot(df, aes(x=x))  + theme(axis.text.x = element_text(size = 10, angle = 90)) + scale_x_continuous(breaks = seq(3000, 4999, by = 100)) + geom_vline(alpha=0.5,xintercept=3343,col=cols[4],lwd=lwd2) +geom_vline(xintercept=3221,alpha=0.5,col=cols[5],lwd=lwd2)    +labs(y="counts (log sacle: y=ln(1+0.1*counts))",x="ISCO-08 classification of occupation")  +
  geom_histogram(aes(y = after_stat(f(count))),breaks=seq(x_min,x_max,100),alpha=0.2,col=cols[1],lwd=lwd1) +
  geom_histogram(aes(y = after_stat(f(count))),breaks=seq(x_min,x_max,10),col=cols[2],alpha=0.5,lwd=lwd1) +
  geom_histogram(aes(y = after_stat(f(count))),breaks=seq(x_min,x_max,1),col=cols[3],alpha=1,lwd=lwd1) + scale_x_continuous(breaks = scales::pretty_breaks(n = 16))
dev.off()

#plt <- ggplot(df) + geom_bar(aes(x=x, y=counts), stat="identity")
#print(plt)


#ggplot(df, aes(x=x)) + geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,1000),col="black",alpha=0.3) + geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,100),alpha=0.2,col="blue") + geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,10),col="red",alpha=0.5) + geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,1),col="green",alpha=0.08)


#ggplot(df, aes(x=x)) + geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,10),col="black",alpha=0.3) + geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,1),alpha=0.2,col="blue") #+ geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,10),col="red",alpha=0.5) + geom_histogram(aes(y = after_stat(count / sum(count))),breaks=seq(0,10000,1),col="green",alpha=0.08)
##

