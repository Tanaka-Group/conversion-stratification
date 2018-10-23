# Conversion oSCORAD to EASI with Isotonic regression

# Initialisation ----------------------------------------------------------

rm(list=ls()) # Clear Workspace

load("data")
rm(Y,data_healthy)

library(ggplot2)
library(corrplot)

source("utils.R")

score=processing_score(data_ad)
df=score[,c("Patient","Week","oSCORAD","EASI")]
df=na.omit(df)

MCID=6.6 # EASI
lmax=c(83,72)

set.seed(15793245) # Reproductibility

# Isotonic Regression EASI vs oSCORAD -----------------------------------------------------

tmp=df[!(df$Patient==120 & df$Week==0),] # Remove outlier
tmp=tmp[order(tmp$oSCORAD),] # Order

ir=with(tmp,isoreg(oSCORAD,EASI))
plot(ir,pch="+",main="",xlab="oSCORAD",ylab="EASI")

abc=abacus_isoreg(ir,lmax)

gp=with(ir,data.frame(oSCORAD=x,EASI=y,Fit=yf))
gp=rbind(gp,data.frame(oSCORAD=lmax[1],EASI=NA,Fit=lmax[2]))
cf=ggplot(data=gp,aes(x=oSCORAD))+
  geom_point(aes(y=EASI),size=1.5)+
  geom_smooth(aes(y=Fit),se=FALSE,colour="#E69F00",size=2,method = "loess")+
  theme_bw(base_size = 15)
cf

# Analysis Isotonic regression --------------------------------------------

# # Residuals
res=ir$y-ir$yf
qqnorm(res);qqline(res) # qqplot
(std_ir=sd(res)) # Std residuals (centered on zero)
lm1=lm(EASI~oSCORAD,data=df) # Comparison with linear model
(std_lm=summary(lm1)$sigma)

(exp_acc= 1-2*pnorm(-abs(MCID/std_ir))) # Expected accuracy

# Crossvalidation EASI vs oSCORAD -----------------------------------------

point_cv=data.frame(Patient=c(),Week=c(),Actual=c(),Prediction=c(),Iteration=c(),Fold=c())

n_cross=25
k=4

tmp=df[!(df$Patient==120 & df$Week==0),] # Remove outlier

pt=as.numeric(as.character(unique(tmp$Patient)))

for (i in 1:n_cross){
  id=sample(cut(pt,breaks=k,labels=FALSE)) # K-fold
  for (j in 1:k){
    test=tmp[grep(paste(pt[id==j],collapse = "|"),tmp$Patient),]
    train=tmp[grep(paste(pt[id!=j],collapse = "|"),tmp$Patient),]
    
    train=train[order(train$oSCORAD),] # Order
    fit=isoreg(train$oSCORAD,train$EASI)
    
    abc=abacus_isoreg(fit,lmax)
    
    p=rep(NA,dim(test)[1])
    for (l in 1:dim(test)[1]){p[l]=pred_abacus(abc,test$oSCORAD[l])}

    point_cv=rbind(point_cv,
                   data.frame(Patient=test$Patient,
                              Week=test$Week,
                              Actual=test$EASI,
                              Prediction=p,
                              Iteration=rep(i,nrow(test)),
                              Fold=rep(j,nrow(test))
                   ))
  }
}

# Analysis crossvalidation ------------------------------------------------

# RMSE, R? and acc
(m=metrics(point_cv,MCID))

# Distribution of RMSE
point_cv$res2=with(point_cv,(Actual-Prediction)^2)
rmse_dis=aggregate(res2~Iteration+Fold,point_cv,mean)
rmse_dis$rmse=sqrt(rmse_dis$res2)
mean(rmse_dis$rmse)
sd(rmse_dis$rmse)

# Prediction plot
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
easi_range=0:(ceiling(max(point_cv$Actual,point_cv$Prediction)))
acl=data.frame(Acc=easi_range,Acc_lower=(easi_range-MCID)*(easi_range>MCID),
               Acc_upper=(easi_range+MCID)*(easi_range+MCID<72)+72*(!easi_range+MCID<72))

cp=ggplot()+
  geom_point(data=point_cv,aes(x=Actual,y=Prediction))+
  geom_ribbon(data=acl,
              aes(x=Acc,ymin=Acc_lower,ymax=Acc_upper),alpha=0.3,fill=cbbPalette[2])+
  geom_line(data=acl,aes(x=Acc,y=Acc),colour=cbbPalette[7],lwd=2,lty=1)+
  labs(x="EASI",y="EASI predicted by oSCORAD")+
  theme_bw(base_size = 15)
cp

# Abacus: oSCORAD to EASI ------------------------------------------------------------------

ir=with(tmp,isoreg(oSCORAD,EASI))
abc=abacus_isoreg(ir,c(83,72))

os_to_easi=data.frame(oSCORAD=0:83,EASI=rep(NA,84))
for (i in 0:83){os_to_easi[i+1,"EASI"]=pred_abacus(abc,i)}

ggplot(data=os_to_easi,aes(x=oSCORAD,y=EASI))+geom_point()+theme_bw(base_size = 15)

# Abacus: EASI to oSCORAD -------------------------------------------------

ir=with(tmp,isoreg(EASI,oSCORAD))
abc=abacus_isoreg(ir,c(72,83))

easi_to_os=data.frame(EASI=0:72,oSCORAD=rep(NA,73))
for (i in 0:72){easi_to_os[i+1,"oSCORAD"]=pred_abacus(abc,i)}

ggplot(data=easi_to_os,aes(x=EASI,y=oSCORAD))+geom_point()+theme_bw(base_size = 15)

# Stratification ----------------------------------------------------------

stratification_heatmap()

ct=(0:4)+.1 # Value TPR
# ct=c(.8,1.2,2.8,3.2,4.8) # Custom one
xl=c(0,ct,5);xl=(xl[-1]+xl[-length(xl)])/2; # Between cutoffs

stratification_heatmap(ct=ct,xl=xl,xi=c(0,5))
