# Processing Score --------------------------------------------------------

processing_score <- function(data_ad){
  n=dim(data_ad)[1]
  Week=c(rep(0,n),rep(2,n),rep(4,n),rep(8,n),rep(12,n),rep(24,n))
  
  SCORAD=c(data_ad$SCORAD_TOT.2,data_ad$SCORAD_TOT.3,
           data_ad$SCORAD_TOT.4,data_ad$SCORAD_TOT.5,
           data_ad$SCORAD_TOT.6,data_ad$SCORAD_TOT.7)
  
  EASI=c(data_ad$EASI_TOTAL.2,data_ad$EASI_TOTAL.3,
         data_ad$EASI_TOTAL.4,data_ad$EASI_TOTAL.5,
         data_ad$EASI_TOTAL.6,data_ad$EASI_TOTAL.7)
  
  POEM=c(data_ad$POEM_TOTAL.2,data_ad$POEM_TOTAL.3,
         data_ad$POEM_TOTAL.4,data_ad$POEM_TOTAL.5,
         data_ad$POEM_TOTAL.6,data_ad$POEM_TOTAL.7)
  
  VAS_ITCH=c(data_ad$ITCH_VAS.2,data_ad$ITCH_VAS.3,
             data_ad$ITCH_VAS.4,data_ad$ITCH_VAS.5,
             data_ad$ITCH_VAS.6,data_ad$ITCH_VAS.7)
  
  VAS_SLEEP_LOSS=c(data_ad$SCSLEEPLOSS.2,data_ad$SCSLEEPLOSS.3,
                   data_ad$SCSLEEPLOSS.4,data_ad$SCSLEEPLOSS.5,
                   data_ad$SCSLEEPLOSS.6,data_ad$SCSLEEPLOSS.7)
  
  oSCORAD=SCORAD-VAS_ITCH-VAS_SLEEP_LOSS
  
  IGA=c(data_ad$IGA.2,data_ad$IGA.3,
        data_ad$IGA.4,data_ad$IGA.5,
        data_ad$IGA.6,data_ad$IGA.7)
  
  PGA=c(data_ad$PAGA.2,data_ad$PAGA.3,
        data_ad$PAGA.4,data_ad$PAGA.5,
        data_ad$PAGA.6,data_ad$PAGA.7)
  
  score=data.frame(Patient=data_ad$PT,Week,SCORAD,oSCORAD,EASI,POEM,VAS_ITCH,VAS_SLEEP_LOSS,IGA,PGA)
  
  return(score)
  
}

# Abacus Isotonic Regression ------------------------------------------------------------------

abacus_isoreg <- function(ir,lmax){
  # Process the isoreg model in order to find correspondence between x and y
  
  if (!ir$isOrd){
    o=order(ir$x)
    ir$x=ir$x[o]
    ir$y=ir$y[o]
  }

  smt=loess.smooth(ir$x,ir$yf,degree=2)
  
  # Add limit values
  smt$x=c(0,smt$x,lmax[1]);smt$y=c(0,smt$y,lmax[2])
  
  # Round value at 0.5
  smt$x=round(smt$x*2)/2
  smt$y=round(smt$y*2)/2
  
  return(smt)
}

# Prediction Isotonic Regression ------------------------------------------

pred_abacus <- function(abc,x0){
  # Prediction (weighted mean) of x0 given an abacus
  
  if (x0==0){p=0}else{
    inf=tail(which(abc$x<=x0),1)
    sup=head(which(abc$x>=x0),1)
    if (inf>=sup){p=abc$y[inf]}else{
      dinf=x0-abc$x[inf]
      dsup=abc$x[sup]-x0
      p=(abc$y[inf]*dsup+abc$y[sup]*dinf)/(dsup+dinf)
      p=round(p*2)/2
    }
  }
  
  return(p)
}

# Metrics -----------------------------------------------------------------

metrics <- function(point_cv,MCID){
  res=list()
  
  resid=with(point_cv,Actual-Prediction)
  resid0=with(point_cv,Actual-mean(Actual))
  
  res$rmse=sqrt(mean(resid^2))
  res$rmse0=sqrt(mean(resid0^2))
  
  res$r2=1-(res$rmse/res$rmse0)^2
  
  res$acc=mean(abs(resid)<MCID)
  res$acc0=mean(abs(resid0)<MCID)
  
  res$kappa=(res$acc-res$acc0)/(1-res$acc0)
  
  return(res)
}

# Heatmap stratification --------------------------------------------------

stratification_heatmap <- function(file="abacus.csv",ct=c(.5,1.5,2.5,3.5,4.5),xl=0:5,xi=0:5){
  # ct cutoffs, xl position of severity labels, xi IGA values
  library(ggplot2)
  library(readr)
  
  abacus <- read_delim(file, ";", escape_double = FALSE, trim_ws = TRUE)
  
  sl=c("Clear","Almost\nclear","Mild","Moderate","Severe","Very\nsevere") # Severity label
  cl=c(rep("black",5),"white") # Colour for severity labels
  cl=rep(c("black","white"),each=3) # Colour for severity labels
  
  pred <- function(abc,lY,x0){
    # lY: EASI or oSCORAD
    # x0: iga to predict
    inf=tail(which(abc$IGA<=x0),1)
    sup=head(which(abc$IGA>=x0),1)
    if (length(inf)==0){p=x0/abc$IGA[sup]*abc[sup,lY][[1]]}else{
      if (length(sup)==0){p=x0/abc$IGA[inf]*abc[inf,lY][[1]]}else{
        if (inf>=sup){p=abc[inf,lY][[1]]}else{
          dinf=x0-abc$IGA[inf]
          dsup=abc$IGA[sup]-x0
          p=(abc[inf,lY][[1]]*dsup+abc[sup,lY][[1]]*dinf)/(dsup+dinf)
        }
      }
    }
    p=round(p*2)/2
    return(p)
  }
  
  xi=sort(c(xi,ct)) # Position for IGA
  xe=0*xi # Position for EASI
  xo=xe # Position for oSCORAD
  
  for (i in 1:length(xi)){
    xe[i]=pred(abacus,"EASI",xi[i])
    xo[i]=pred(abacus,"oSCORAD",xi[i])
  }
  
  ds=data.frame(x=c(rep(-.5,3),rep(xi,3)),
                y=c(1.1,1.2,1.3,rep(c(1.1,1.2,1.3),each=length(xi))),
                l=c("IGA","EASI","oSCORAD",xi,xe,xo))
  
  ds$x[ds$x==0]=-.1;xl[1]=-.1 # Repositioning text for IGA=0
  
  p=ggplot()+
    geom_tile(data=data.frame(IGA=seq(0,5,.01)),aes(x=IGA,fill=IGA,y=.5))+ # x Colour as a function of IGA
    scale_fill_gradientn(colours=c("#FFFFFF",rev(heat.colors(3)),"#000000"))+ # Colour gradient
    geom_text(aes(x=xl,y=.5,label=sl),colour=cl,fontface="bold",size=8)+ # Severity labels
    geom_text(data=ds,aes(x=x,y=y,label=l),fontface="bold",size=7)+ # Scale
    geom_segment(aes(x=ct,xend=ct,y=0,yend=1),colour=head(cl,5),size=2,linetype="longdash")+ # Vertical line cutoff
    labs(x="",y="")+ # Remove axis labels
    ylim(c(0,2.5))+ # Window size (y)
    theme_classic(base_size = 15)+
    theme(legend.position="none",
          axis.text.y=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),
          axis.line.y=element_blank(),axis.line.x=element_blank()) # Remove graphical elements
  
  return(p)
}
