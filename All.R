rm(list=ls())

setwd('this is where the datafile directory would be specified')
data = read.csv("dataonefive.csv")

data$newstateid<-as.character(data$stateid)
kidnapdat <- data[!is.na(data$lpwomgirl),]

# define clustered standard errors function
cl <- function(dat,fm, cluster) {
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum))
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  return(list(coeftest(fm, vcovCL), n=round(length(fm$fitted.values),0), rsq=round(summary(fm)$r.squared,2)))
}

# fix data formatting issues and subset kidnapping data so that cl function will work
data$newstateid<-as.character(data$stateid)

#Creating a reservation year variable
data$year.res<-NULL
data$year.res[data$stateid == 'ANDHRA PRADESH']<-1995
data$year.res[data$stateid == 'ASSAM']<-2002
data$year.res[data$stateid == 'BIHAR']<-2001
data$year.res[data$stateid == 'GUJARAT']<-1995
data$year.res[data$stateid == 'HARYANA']<-1995
data$year.res[data$stateid == 'HIMACHAL PRADESH']<-1995
data$year.res[data$stateid == 'JAMMU & KASHMIR']<-2001
data$year.res[data$stateid == 'KARNATAKA']<-1987
data$year.res[data$stateid == 'KERALA']<-1991
data$year.res[data$stateid == 'MADHYA PRADESH']<-1994
data$year.res[data$stateid == 'MAHARASHTRA']<-1992
data$year.res[data$stateid == 'ORISSA']<-1992
data$year.res[data$stateid == 'PUNJAB']<-1994
data$year.res[data$stateid == 'RAJASTHAN']<-1995
data$year.res[data$stateid == 'TAMIL NADU']<-1996
data$year.res[data$stateid == 'UTTAR PRADESH']<-2006
data$year.res[data$stateid == 'WEST BENGAL']<-1993

##TABLE 2
##Duration model

install.packages("survival")
library(survival)
library(splines)
library(stargazer)

#Generating a new dataset with all variables at their 1985 values
newdata <- data[which(data$year==1985),]

Date.imp <-survreg(Surv(newdata$year.res)~ newdata$pcgsdp  + newdata$pfemale + newdata$plit + newdata$pwlit + newdata$prural, data=newdata, dist="weibull")
summary(Date.imp)
stargazer(Date.imp, covariate.labels = c("GDP per capita", "Female-male ratio", "Literacy", "Womens literacy", "Percent rural"), title = "Duration model predicting time to state level implementation of the national reservation policy",
          column.labels="Reservation implementation year")


##################################################################
##TABLE 1
#Analysis of the data in Iyer et al. on and after 1995 in states where the reservation policy was implemented in 1995 and after

#Subsetting data to years post-1994 and reservation post-1994
better.data <- data[which(data$year.res>=1995),]
better.data <- data[which(better.data$year>=1995),]
View(better.data)

#Running the initial model with post-1995 data, column 1
better.data$newstateid<-as.character(better.data$stateid)
new.column1<-lm(better.data$lpcr_womtot ~ better.data$postwres + factor(better.data$year) + factor(better.data$newstateid), subset=(better.data$majstate==1))
summary(new.column1)
column1.cl<- cl(better.data, new.column1, better.data$newstateid)
column1.cl

#Running the initial model with post-1995 data, column 6
better.data$newstateid<-as.character(better.data$stateid)
new.total.crimes<-lm(better.data$lpcr_womtot ~ better.data$postwres + better.data$newstateid:better.data$year + better.data$pfemale + better.data$prural + better.data$plit + better.data$pfarm + better.data$womancm + better.data$pcgsdp + better.data$ppol_strengt +
              factor(better.data$year) + factor(better.data$newstateid), subset=(better.data$year>=1985 & better.data$majstate==1))
summary(new.total.crimes)
column2.cl<- cl(better.data, new.total.crimes, better.data$newstateid)
column2.cl

#Table 3 column 1 code for total crimes against women in Iyer
model1w<-lm(data$lpcr_womtot ~ data$postwres +  
              factor(data$year) + factor(data$newstateid), subset=(data$year>=1985 & data$majstate==1))
summary(model1w)
model1cw<- cl(data, model1w, data$newstateid)
model1cw

#Table 3 column 6 code for total crimes against women in Iyer
model6w<-lm(data$lpcr_womtot ~ data$postwres + data$newstateid:data$year + data$pfemale + data$prural + data$plit + data$pfarm + data$womancm + data$pcgsdp + data$ppol_strengt +
              factor(data$year) + factor(data$newstateid), subset=(data$year>=1985 & data$majstate==1))
summary(model6w)
model6cw<- cl(data, model6w, data$newstateid)
model6cw

#Outputting the results from the model restricted to post-1995
table2.1<-as.matrix(c(round(column1.cl[[1]]["better.data$postwres",1],3),paste("[",round(column1.cl[[1]]["better.data$postwres",2],3),"]"),round(column1.cl[[3]],2), round(column1.cl[[2]],0) ))
table2.2<-as.matrix(c(round(column2.cl[[1]]["better.data$postwres",1],3),paste("[",round(column2.cl[[1]]["better.data$postwres",2],3),"]"),round(column2.cl[[3]],2), round(column2.cl[[2]],0) ))
table2.3<-as.matrix(c(round(model1cw[[1]]["data$postwres",1],3),paste("[",round(model1cw[[1]]["data$postwres",2],3),"]"),round(model1cw[[3]],2), round(model1cw[[2]],0) ))
table2.4<-as.matrix(c(round(model6cw[[1]]["data$postwres",1],3),paste("[",round(model6cw[[1]]["data$postwres",2],3),"]"),round(model6cw[[3]],2), round(model6cw[[2]],0) ))

table2_out <- matrix(NA,nrow=8,ncol=2)
table2_out[1:4,1] <- table2.1
table2_out[1:4,2] <- table2.2
table2_out[5:8,1] <- table2.3
table2_out[5:8,2] <- table2.4

table2_out
colnames(table2_out) <- c("No controls (1)", "Control for state-specific time trends + other controls (2)")
rownames(table2_out)<- c("Total crimes against women per 1,000 women", "", "R2", "Observations",
                         "Total crimes against women per 1,000 women", "", "R2", "Observations")

library(stargazer)
stargazer(table2_out,title="Women's Political Representation and Crimes against Women")


########################
##FIGURE 6
#Using Zelig to look at changes in female-male ratio
install.packages("Zelig")
library(Zelig)
remove.packages("Zelig")
install.packages("Zelig", repos = "http://r.iq.harvard.edu/archived/", type = "source")

Date.imp <-survreg(Surv(newdata$year.res)~ pcgsdp  + pfemale + plit + pwlit + prural, data=newdata, dist="weibull")
summary(Date.imp)

z.out <- zelig(Surv(year.res,) ~ pcgsdp  + pfemale + plit + pwlit + prural, model = "weibull", data=newdata)
summary(z.out)
z.out$coef

# Trying a female-male ratio 1.0

#Setting female-male ratio at 1.0 and all other variables at their median 
smalldata<- cbind(newdata$year.res, newdata$pcgsdp, newdata$pfemale, newdata$plit, newdata$pwlit, newdata$prural)
colnames(smalldata)<-c("year.res", "pcgsdp", "pfemale", "plit", "pwlit", "prural")
X<-smalldata[,2:6]
X2 <- as.matrix(cbind(1,X))
colnames(X2)[1] <- "Intercept"
Xc <- apply(X2, MARGIN = 2, FUN = median)
Xc[3]<- 1


#setting x values
x.evs <- setx(z.out, Intercept = 1, pcgsdp = 1.1879470, pfemale = 1.0, plit = 0.4032032, pwlit = 0.3091190, prural = 0.3091190)

#Simulating predicted values
set.seed(1234)
zelig.sim <- sim(z.out, x.evs)
summary(zelig.sim)
names(zelig.sim)
plot(zelig.sim)
summary(zelig.sim)
yl <- "Predicted Probability of Voting"

#Graphing predicted probabilities
graphics.off()
plot(density(zelig.sim$qi$pr), xlab = "Year of policy implementation", main= "Predicted year of policy implementation for an \n average state with a female-male ratio equal to 1")
mean(zelig.sim$qi$pr)

##################################
##Code for Figures 2-5
##################################

rm(list=ls())

setwd('Where the working directory would be specified')
require(foreign)
tables1to5<-as.data.frame(read.dta('tables1to5.dta'))

tables1to5$newstateid<-as.character(tables1to5$stateid)
tables1to5$numid<-as.numeric(factor(tables1to5$newstateid))

data<-tables1to5

data$year.res<-NULL
data$year.res[data$stateid == 'ANDHRA PRADESH']<-1995
data$year.res[data$stateid == 'ASSAM']<-2002
data$year.res[data$stateid == 'BIHAR']<-2001
data$year.res[data$stateid == 'GUJARAT']<-1995
data$year.res[data$stateid == 'HARYANA']<-1995
data$year.res[data$stateid == 'HIMACHAL PRADESH']<-1995
data$year.res[data$stateid == 'JAMMU & KASHMIR']<-2001
data$year.res[data$stateid == 'KARNATAKA']<-1987
data$year.res[data$stateid == 'KERALA']<-1991
data$year.res[data$stateid == 'MADHYA PRADESH']<-1994
data$year.res[data$stateid == 'MAHARASHTRA']<-1992
data$year.res[data$stateid == 'ORISSA']<-1992
data$year.res[data$stateid == 'PUNJAB']<-1994
data$year.res[data$stateid == 'RAJASTHAN']<-1995
data$year.res[data$stateid == 'TAMIL NADU']<-1996
data$year.res[data$stateid == 'UTTAR PRADESH']<-2006
data$year.res[data$stateid == 'WEST BENGAL']<-1993

data$year.diff<-data$year - data$year.res
data$view.year<-data$year
min(data$year.res)

data$duration<-data$year.res-1985



#install.packages('ggplot2')

require('ggplot2')

data$res<-"Post"
data$res[data$postwres==0]<-"Pre"
###############################################################
##Figure 2

p<-ggplot(data, aes(y=pcr_womtot, x=year, group=newstateid, colour = factor(res))) + 
  geom_point() +facet_wrap( ~ newstateid, ncol = 3) +geom_vline(xintercept = 1993) +
  xlab("Year")+ ylab("Total Crimes against women per 1000 women") +  
  scale_colour_discrete(name = "Reservation Implementation")+theme(legend.position="bottom")

print(p)

##Figure 3

p1<-ggplot(data, aes(y=pwomgirl, x=year, group=newstateid, colour = factor(res))) + 
  geom_point() +facet_wrap( ~ newstateid, ncol = 3) +geom_vline(xintercept = 1993) +
  xlab("Year")+ ylab("Kidnapping of women per 1000 women") +  
  scale_colour_discrete(name = "Reservation Implementation")+theme(legend.position="bottom")

print(p1)

###############################################################
#Figure 4

p2<-ggplot(data, aes(y=prape2, x=year, group=newstateid, colour = factor(res))) + 
  geom_point() +facet_wrap( ~ newstateid, ncol = 3) +geom_vline(xintercept = 1993) +
  xlab("Year")+ ylab("Rapes per 1000 women") +  
  scale_colour_discrete(name = "Reservation Implementation")+theme(legend.position="bottom")

print(p2)

##########################################################


#############################################################
rm(list=ls())
require('ggplot2')

new.data <- read.csv("This dataset is not publicly available: Obtained through personal communication with Professor Laxmi Iyer")
#sub<-new.data[new.data$Data<10000,]

sub<-new.data[new.data$Category!="ARSON",]
sub<-sub[sub$Category!="ATMTMURDER",]
sub<-sub[sub$Category!="CULP",]
sub<-sub[sub$Category!="CUSTODIAL RAPE",]
sub<-sub[sub$Category!="OTHER RAPE",]
sub<-sub[sub$Category!="OTHERS",]
sub<-sub[sub$Category!="DACOITY",]
sub<-sub[sub$Category!="PREP",]
sub<-sub[sub$Category!="BURGLARY",]
sub<-sub[sub$Category!="AUTOTHEFT",]
sub<-sub[sub$Category!="OTHER THEFT",]
sub<-sub[sub$Category!="BREACH TRUST",]
sub<-sub[sub$Category!="COUNTERFIET",]
sub<-sub[sub$Category!="HURT",]
sub<-sub[sub$Category!="THEFT",]

sub2<-na.omit(sub)

###############################################################
##Figure 5
p10<-ggplot(sub2, aes(y=Data, x=Year)) + 
  geom_point() +facet_wrap( ~ Category, ncol = 3) +geom_vline(xintercept = 1995) +
  xlab("Year")+ ylab("Number of Crimes") 

print(p10)


