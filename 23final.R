cuttab<-function(fit,var,wdtmp,cutpoint=NULL,pdec=NULL,p.plrt=NULL) {
  UseMethod("cuttab")
}

cuttab.glm<-function(fit,var,wdtmp,cutpoint=NULL,pdec=NULL,p.plrt=NULL) {
  fit<-fit;var<-var;wdtmp<-as.data.frame(wdtmp);cutpoint<-cutpoint
  if (inherits(fit,  "glm") & fit[["family"]][["family"]]=="binomial") {
    out<-cuttab.glmbinomial(fit = fit,var = var,wdtmp = wdtmp,cutpoint=cutpoint,pdec<-pdec,p.plrt<-p.plrt)
  }
  if (inherits(fit,  "glm") & fit[["family"]][["family"]]=="gaussian") {
    out<-cuttab.glmgaussian(fit = fit,var = var,wdtmp = wdtmp,cutpoint=cutpoint,pdec<-pdec,p.plrt<-p.plrt)
  }
  out   
}

cuttab.lrm<-function(fit,var,wdtmp,pdec=NULL) {
  if (is.null(pdec)) {pdec<-3} else {pdec<-pdec}
  cut.p<-gettpval(fit=fit,var =var,wdtmp=wdtmp)
  dt<-as.data.frame(wdtmp)
  x <-wdtmp[,var]
  X1 <-(x<=cut.p)*(x-cut.p)
  X2 <-(x> cut.p)*(x-cut.p)
  dt <-cbind(dt,X1,X2)
  dd <- datadist(dt)
  options(datadist="dd")
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2];fml1<-paste(fml,"+X2");
  fml2<-transformfml(fit,var)
  fml2<-paste(fml2,"+X1+X2");
  fit0<-lrm(formula(fml),weights=dt$weights,data=dt)
  fit1<-lrm(formula(fml1),weights=dt$weights,data=dt)
  fit2<-lrm(formula(fml2),weights=dt$weights,data=dt)
  t1.1<-summary(fit0)
  t1.2<-round(exp(summary(fit0)[var,4]),pdec);t1.3<-round(anova(fit0)[var,3],pdec);
  t1.4<-round(exp(summary(fit0)[var,6]),pdec);t1.5<-round(exp(summary(fit0)[var,7]),pdec)
  t1<-paste0(t1.2,"(",t1.4,"-",t1.5,")",t1.3)
  t2<-coef(summary(fit1));t2.1<-summary(fit1)$conf.int;
  t2.2<-round(coef(summary(fit1))[var,2],pdec);t2.3<-round(coef(summary(fit1))[var,5],pdec)
  t2.4<-round(summary(fit1)$conf.int[var,3],pdec);t2.5<-round(summary(fit1)$conf.int[var,4],pdec)
  t2<-paste0(t2.2,"(",t2.4,"-",t2.5,")",t2.3)
  t3<-coef(summary(fit2));t3.1<-summary(fit2)$conf.int;
  t3.2<-round(t3[paste0("X2"),2],pdec);t3.3<-round(t3[paste0("X2"),5],pdec);
  t3.4<-round(t3.1[paste0("X2"),3],pdec);t3.5<-round(t3.1[paste0("X2"),4],pdec)
  t3<-paste0(t3.2,"(",t3.4,"-",t3.5,")",t3.3)
  i1<-paste0("<",cut.p)
  i2<-paste0(">",cut.p)
  plrt<-pvformat(1-pchisq(2*(logLik(fit2)[1]-logLik(fit0)[1]),1),3)
  d1<-c("Model 1 Fitting model by standard linear regression",t1)
  d2<-c("Model 2 Fitting model by two-piecewise linear regression","")
  d3<-c("Inflection point",cut.p)
  d4<-c(i1,t2)
  d5<-c(i2,t3)
  d6<-c("P for likelihood ratio test",plrt)
  dat<-rbind(d1,d2,d3,d4,d5,d6)
  colnames(dat)<-c("Outcome","the effect size, 95%CI, P value")
  dat
}

cuttab.cph<-function(fit,var,wdtmp,cutpoint=NULL,pdec=NULL,p.plrt=NULL) {
  if (is.null(pdec)) {pdec<-3} else {pdec<-pdec}
  if (is.null(p.plrt)) {p.plrt<-3} else {p.plrt<-p.plrt}
  cut.p<-gettpval(fit=fit,var =var,wdtmp=wdtmp)
  if (!is.null(cutpoint)) {cut.p<-cutpoint}
  dt<-as.data.frame(wdtmp)
  x <-wdtmp[,var]
  X1 <-(x<=cut.p)*(x-cut.p)
  X2 <-(x> cut.p)*(x-cut.p)
  dt <-cbind(dt,X1,X2)
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2];fml1<-paste(fml,"+X2");
  fml2<-transformfml(fit,var)
  fml2<-paste(fml2,"+X1+X2");
  fit0<-coxph(formula(fml),weights=wdtmp$weights,data=wdtmp)
  fit1<-coxph(formula(fml1),weights=wdtmp$weights,data=wdtmp)
  fit2<-coxph(formula(fml2),weights=wdtmp$weights,data=wdtmp)
  t1<-coef(summary(fit0));t1.1<-summary(fit0)$conf.int
  t1.2<-round(coef(summary(fit0))[var,2],pdec);t1.3<-round(coef(summary(fit0))[var,5],pdec);
  t1.4<-round(summary(fit0)$conf.int[var,3],pdec);t1.5<-round(summary(fit0)$conf.int[var,4],pdec)
  t1<-paste0(t1.2,"(",t1.4,"-",t1.5,")",t1.3)
  t2<-coef(summary(fit1));t2.1<-summary(fit1)$conf.int;
  t2.2<-round(coef(summary(fit1))[var,2],pdec);t2.3<-round(coef(summary(fit1))[var,5],pdec)
  t2.4<-round(summary(fit1)$conf.int[var,3],pdec);t2.5<-round(summary(fit1)$conf.int[var,4],pdec)
  t2<-paste0(t2.2,"(",t2.4,"-",t2.5,")",t2.3)
  t3<-coef(summary(fit2));t3.1<-summary(fit2)$conf.int;
  t3.2<-round(t3[paste0("X2"),2],pdec);t3.3<-round(t3[paste0("X2"),5],pdec);
  t3.4<-round(t3.1[paste0("X2"),3],pdec);t3.5<-round(t3.1[paste0("X2"),4],pdec)
  t3<-paste0(t3.2,"(",t3.4,"-",t3.5,")",t3.3)
  i1<-paste0("<",cut.p)
  i2<-paste0(">",cut.p)
  plrt<-pvformat(1-pchisq(2*(logLik(fit2)[1]-logLik(fit0)[1]),1),p.plrt)
  d1<-c("Model 1 Fitting model by standard linear regression",t1)
  d2<-c("Model 2 Fitting model by two-piecewise linear regression","")
  d3<-c("Inflection point",cut.p)
  d4<-c(i1,t2)
  d5<-c(i2,t3)
  d6<-c("P for likelihood ratio test",plrt)
  dat<-rbind(d1,d2,d3,d4,d5,d6)
  colnames(dat)<-c("Outcome","the effect size, 95%CI, P value")
  dat
}

cuttab.glmbinomial<-function(fit,var,wdtmp,cutpoint=NULL,pdec=NULL,p.plrt=NULL) {
  if (is.null(pdec)) {pdec<-3} else {pdec<-pdec}
  if (is.null(p.plrt)) {p.plrt<-3} else {p.plrt<-p.plrt}
  cut.p<-gettpval(fit=fit,var =var,wdtmp=wdtmp)
  if (!is.null(cutpoint)) {cut.p<-cutpoint}
  dt<-as.data.frame(wdtmp)
  x <-wdtmp[,var]
  X1 <-(x<=cut.p)*(x-cut.p)
  X2 <-(x> cut.p)*(x-cut.p)
  dt <-cbind(dt,X1,X2)
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2];fml1<-paste(fml,"+X2");
  fml2<-transformfml(fit,var)
  fml2<-paste(fml2,"+X1+X2");
  fit0<-glm(formula(fml),family = binomial("logit"),data=wdtmp)
  fit1<-glm(formula(fml1),family = binomial("logit"),data=wdtmp)
  fit2<-glm(formula(fml2),family = binomial("logit"),data=wdtmp)
  t1<-coef(summary(fit0));t1.1<-confint(fit0)
  t1.2<-round(exp(coef(summary(fit0))[var,1]),pdec)
  t1.3<-round(coef(summary(fit0))[var,4],pdec)
  t1.4<-round(exp(confint(fit0)[var,1]),pdec)
  t1.5<-round(exp(confint(fit0)[var,2]),pdec)
  t1<-paste0(t1.2,"(",t1.4,"-",t1.5,")",t1.3)
  t2<-coef(summary(fit1));t2.1<-confint(fit1)
  t2.2<-round(exp(coef(summary(fit1))[var,1]),pdec)
  t2.3<-round(coef(summary(fit1))[var,4],pdec)
  t2.4<-round(exp(confint(fit1)[var,1]),pdec)
  t2.5<-round(exp(confint(fit1)[var,2]),pdec)
  t2<-paste0(t2.2,"(",t2.4,"-",t2.5,")",t2.3)
  t3<-coef(summary(fit2));t3.1<-confint(fit2)
  t3.2<-round(exp(t3[paste0("X2"),1]),pdec);t3.3<-round(t3[paste0("X2"),4],pdec);
  t3.4<-round(exp(t3.1[paste0("X2"),1]),pdec);t3.5<-round(exp(t3.1[paste0("X2"),2]),pdec)
  t3<-paste0(t3.2,"(",t3.4,"-",t3.5,")",t3.3)
  i1<-paste0("<",cut.p)
  i2<-paste0(">",cut.p)
  plrt<-pvformat(1-pchisq(2*(logLik(fit2)[1]-logLik(fit0)[1]),1),p.plrt)
  d1<-c("Model 1 Fitting model by standard linear regression",t1)
  d2<-c("Model 2 Fitting model by two-piecewise linear regression","")
  d3<-c("Inflection point",cut.p)
  d4<-c(i1,t2)
  d5<-c(i2,t3)
  d6<-c("P for likelihood ratio test",plrt)
  dat<-rbind(d1,d2,d3,d4,d5,d6)
  colnames(dat)<-c("Outcome","the effect size, 95%CI, P value")
  dat
}

cuttab.glmgaussian<-function(fit,var,wdtmp,cutpoint=NULL,pdec=NULL,p.plrt=NULL) {
  if (is.null(pdec)) {pdec<-3} else {pdec<-pdec}
  if (is.null(p.plrt)) {p.plrt<-3} else {p.plrt<-p.plrt}
  cut.p<-gettpval(fit=fit,var =var,wdtmp=wdtmp)
  if (!is.null(cutpoint)) {cut.p<-cutpoint}
  dt<-as.data.frame(wdtmp)
  x <-wdtmp[,var]
  X1 <-(x<=cut.p)*(x-cut.p)
  X2 <-(x> cut.p)*(x-cut.p)
  dt <-cbind(dt,X1,X2)
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2];fml1<-paste(fml,"+X2");
  fml2<-transformfml(fit,var)
  fml2<-paste(fml2,"+X1+X2");
  fit0<-glm(formula(fml),family = gaussian(link = "identity"),data=wdtmp)
  fit1<-glm(formula(fml1),family = gaussian(link = "identity"),data=wdtmp)
  fit2<-glm(formula(fml2),family = gaussian(link = "identity"),data=wdtmp)
  t1<-coef(summary(fit0));t1.1<-confint(fit0)
  t1.2<-round(coef(summary(fit0))[var,1],pdec)
  t1.3<-round(coef(summary(fit0))[var,4],pdec)
  t1.4<-round(confint(fit0)[var,1],pdec)
  t1.5<-round(confint(fit0)[var,2],pdec)
  t1<-paste0(t1.2,"(",t1.4,"-",t1.5,")",t1.3)
  t2<-coef(summary(fit1));t2.1<-confint(fit1)
  t2.2<-round(coef(summary(fit1))[var,1],pdec)
  t2.3<-round(coef(summary(fit1))[var,4],pdec)
  t2.4<-round(confint(fit1)[var,1],pdec)
  t2.5<-round(confint(fit1)[var,2],pdec)
  t2<-paste0(t2.2,"(",t2.4,"-",t2.5,")",t2.3)
  t3<-coef(summary(fit2));t3.1<-confint(fit2)
  t3.2<-round(t3[paste0("X2"),1],pdec);t3.3<-round(t3[paste0("X2"),4],pdec);
  t3.4<-round(t3.1[paste0("X2"),1],pdec);t3.5<-round(t3.1[paste0("X2"),2],pdec)
  t3<-paste0(t3.2,"(",t3.4,"-",t3.5,")",t3.3)
  i1<-paste0("<",cut.p)
  i2<-paste0(">",cut.p)
  plrt<-pvformat(1-pchisq(2*(logLik(fit2)[1]-logLik(fit0)[1]),1),p.plrt)
  d1<-c("Model 1 Fitting model by standard linear regression",t1)
  d2<-c("Model 2 Fitting model by two-piecewise linear regression","")
  d3<-c("Inflection point",cut.p)
  d4<-c(i1,t2)
  d5<-c(i2,t3)
  d6<-c("P for likelihood ratio test",plrt)
  dat<-rbind(d1,d2,d3,d4,d5,d6)
  colnames(dat)<-c("Outcome","the effect size, 95%CI, P value")
  dat
}
############
gettpval<-function(fit,var,wdtmp, tppmin=NA, tppmax=NA,dec=NA) {
  UseMethod("gettpval")
}


gettpval.lrm<-function(fit,var,wdtmp, tppmin=NA, tppmax=NA,dec=NA) {
  if (missing(dec)) dec<-3
  wdtmp<-as.data.frame(wdtmp)
  xTMP <- wdtmp[,var]
  tmp.ss<-seq(0.05,0.95,0.05)
  tp<-quantile(xTMP,probs=tmp.ss,na.rm=TRUE) 
  tmp.llk<-rep(NA,length(tmp.ss))
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2]
  fml<-paste(fml,"+tmp.X")
  if (!is.na(tppmin) & !is.na(tppmax)) {
    tp2.min = tppmin; tp2.max = tppmax;
  } else {
    for (k in (1:length(tmp.ss))) {
      tmp.X<-(xTMP > tp[k])*(xTMP-tp[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- lrm(formula(fml),weights=wdtmp$weights,data=wdtmp)
      tmp.llk[k]<-logLik(tmp.mdl);
      rm(wdtmp1, tmp.X)
    }
    tp1<-tmp.ss[which.max(tmp.llk)]
    tp2.min = tp1 - 0.04
    tp2.max = tp1 + 0.04
    if (tp2.min<0.05) {tp2.min=0.05}
    if (tp2.max>0.95) {tp2.max=0.95}
  }	
  tp.pctlrange<-quantile(xTMP,probs=c(tp2.min,tp2.max),na.rm=TRUE)
  tp.range<-unique(xTMP[xTMP>tp.pctlrange[1] & xTMP<tp.pctlrange[2]])
  while (length(tp.range)>5) {
    tmp.pct3<-quantile(tp.range,probs=c(0,0.25,0.5,0.75,1),type=3)
    tmp.llk3<-rep(NA,3)
    for (k in (2:4)) {
      tmp.X<-(xTMP>tmp.pct3[k])*(xTMP-tmp.pct3[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- lrm(formula(fml),weights=wdtmp$weights,data=wdtmp)
      
      tmp.llk3[k-1]<-logLik(tmp.mdl);
      rm(wdtmp1, tmp.X)
    }
    tmp.min3<-which.max(tmp.llk3)
    tp.range<-tp.range[tp.range>=tmp.pct3[tmp.min3] & tp.range<=tmp.pct3[tmp.min3+2]]
  }
  if (length(tp.range)>0) {
    if (length(tp.range)==1) {tp.val=tp.range[1];} else {
      tmp.llk<-rep(NA,length(tp.range))
      for (k in (1:length(tp.range))) {
        tmp.X<-(xTMP>tp.range[k])*(xTMP-tp.range[k]); wdtmp1<-cbind(wdtmp,tmp.X)
        tmp.mdl<- lrm(formula(fml),weights=wdtmp$weights,data=wdtmp)
        
        tmp.llk[k]<-logLik(tmp.mdl);
        rm(wdtmp1, tmp.X)
      }
      tp.val<-tp.range[which.max(tmp.llk)]
    }
  } else { tp.val<-tp.pctlrange[1];}
  return(round(tp.val,dec));
}

gettpval.cph<-function(fit,var,wdtmp, tppmin=NA, tppmax=NA,dec=NA) {
  if (missing(dec)) dec<-3
  wdtmp<-as.data.frame(wdtmp)
  xTMP <- wdtmp[,var]
  tmp.ss<-seq(0.05,0.95,0.05)
  tp<-quantile(xTMP,probs=tmp.ss,na.rm=TRUE) 
  tmp.llk<-rep(NA,length(tmp.ss))
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2]
  fml<-paste(fml,"+tmp.X")
  if (!is.na(tppmin) & !is.na(tppmax)) {
    tp2.min = tppmin; tp2.max = tppmax;
  } else {
    for (k in (1:length(tmp.ss))) {
      tmp.X<-(xTMP > tp[k])*(xTMP-tp[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- coxph(formula(fml),weights=wdtmp$weights,data=wdtmp,na.action=na.omit,method="breslow")
      tmp.llk[k]<-tmp.mdl$loglik[2];
      rm(wdtmp1, tmp.X)
    }
    tp1<-tmp.ss[which.max(tmp.llk)]
    tp2.min = tp1 - 0.04
    tp2.max = tp1 + 0.04
    if (tp2.min<0.05) {tp2.min=0.05}
    if (tp2.max>0.95) {tp2.max=0.95}
  }	
  tp.pctlrange<-quantile(xTMP,probs=c(tp2.min,tp2.max),na.rm=TRUE)
  tp.range<-unique(xTMP[xTMP>tp.pctlrange[1] & xTMP<tp.pctlrange[2]])
  while (length(tp.range)>5) {
    tmp.pct3<-quantile(tp.range,probs=c(0,0.25,0.5,0.75,1),type=3)
    tmp.llk3<-rep(NA,3)
    for (k in (2:4)) {
      tmp.X<-(xTMP>tmp.pct3[k])*(xTMP-tmp.pct3[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- coxph(formula(fml),data=wdtmp,na.action=na.omit)
      
      tmp.llk3[k-1]<-tmp.mdl$loglik[2];
      rm(wdtmp1, tmp.X)
    }
    tmp.min3<-which.max(tmp.llk3)
    tp.range<-tp.range[tp.range>=tmp.pct3[tmp.min3] & tp.range<=tmp.pct3[tmp.min3+2]]
  }
  if (length(tp.range)>0) {
    if (length(tp.range)==1) {tp.val=tp.range[1];} else {
      tmp.llk<-rep(NA,length(tp.range))
      for (k in (1:length(tp.range))) {
        tmp.X<-(xTMP>tp.range[k])*(xTMP-tp.range[k]); wdtmp1<-cbind(wdtmp,tmp.X)
        tmp.mdl<- coxph(formula(fml),data=wdtmp,na.action=na.omit)
        
        tmp.llk[k]<-tmp.mdl$loglik[2];
        rm(wdtmp1, tmp.X)
      }
      tp.val<-tp.range[which.max(tmp.llk)]
    }
  } else { tp.val<-tp.pctlrange[1];}
  return(round(tp.val,dec));
}

gettpval.glmbinomial<-function(fit,var,wdtmp, tppmin=NA, tppmax=NA,dec=NA) {
  fit<-fit;var<-var;wdtmp<-as.data.frame(wdtmp)
  if (missing(dec)) dec<-3
  xTMP <- wdtmp[,var]
  tmp.ss<-seq(0.05,0.95,0.05)
  tp<-quantile(xTMP,probs=tmp.ss,na.rm=TRUE) 
  tmp.llk<-rep(NA,length(tmp.ss))
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2]
  fml<-paste(fml,"+tmp.X")
  if (!is.na(tppmin) & !is.na(tppmax)) {
    tp2.min = tppmin; tp2.max = tppmax;
  } else {
    for (k in (1:length(tmp.ss))) {
      tmp.X<-(xTMP > tp[k])*(xTMP-tp[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- glm(formula(fml),weights=wdtmp$weights,data=wdtmp,family = binomial("logit"))
      tmp.llk[k]<-logLik(tmp.mdl);
      rm(wdtmp1, tmp.X)
    }
    tp1<-tmp.ss[which.max(tmp.llk)]
    tp2.min = tp1 - 0.04
    tp2.max = tp1 + 0.04
    if (tp2.min<0.05) {tp2.min=0.05}
    if (tp2.max>0.95) {tp2.max=0.95}
  }	
  tp.pctlrange<-quantile(xTMP,probs=c(tp2.min,tp2.max),na.rm=TRUE)
  tp.range<-unique(xTMP[xTMP>tp.pctlrange[1] & xTMP<tp.pctlrange[2]])
  while (length(tp.range)>5) {
    tmp.pct3<-quantile(tp.range,probs=c(0,0.25,0.5,0.75,1),type=3)
    tmp.llk3<-rep(NA,3)
    for (k in (2:4)) {
      tmp.X<-(xTMP>tmp.pct3[k])*(xTMP-tmp.pct3[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- glm(formula(fml),weights=wdtmp$weights,data=wdtmp,family = binomial("logit"))
      
      tmp.llk3[k-1]<-logLik(tmp.mdl);
      rm(wdtmp1, tmp.X)
    }
    tmp.min3<-which.max(tmp.llk3)
    tp.range<-tp.range[tp.range>=tmp.pct3[tmp.min3] & tp.range<=tmp.pct3[tmp.min3+2]]
  }
  if (length(tp.range)>0) {
    if (length(tp.range)==1) {tp.val=tp.range[1];} else {
      tmp.llk<-rep(NA,length(tp.range))
      for (k in (1:length(tp.range))) {
        tmp.X<-(xTMP>tp.range[k])*(xTMP-tp.range[k]); wdtmp1<-cbind(wdtmp,tmp.X)
        tmp.mdl<- glm(formula(fml),weights=wdtmp$weights,data=wdtmp,family = binomial("logit"))
        
        tmp.llk[k]<-logLik(tmp.mdl);
        rm(wdtmp1, tmp.X)
      }
      tp.val<-tp.range[which.max(tmp.llk)]
    }
  } else { tp.val<-tp.pctlrange[1];}
  return(round(tp.val,dec));
}

gettpval.glmgaussian<-function(fit,var,wdtmp, tppmin=NA, tppmax=NA,dec=NA) {
  fit<-fit;var<-var;wdtmp<-as.data.frame(wdtmp)
  if (missing(dec)) dec<-3
  xTMP <- wdtmp[,var]
  tmp.ss<-seq(0.05,0.95,0.05)
  tp<-quantile(xTMP,probs=tmp.ss,na.rm=TRUE) 
  tmp.llk<-rep(NA,length(tmp.ss))
  call<-fit[["call"]]
  ff<- paste0(call)
  fml<-ff[2]
  fml<-paste(fml,"+tmp.X")
  if (!is.na(tppmin) & !is.na(tppmax)) {
    tp2.min = tppmin; tp2.max = tppmax;
  } else {
    for (k in (1:length(tmp.ss))) {
      tmp.X<-(xTMP > tp[k])*(xTMP-tp[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- glm(formula(fml),weights=wdtmp$weights,data=wdtmp,gaussian(link = "identity"))
      tmp.llk[k]<-logLik(tmp.mdl);
      rm(wdtmp1, tmp.X)
    }
    tp1<-tmp.ss[which.max(tmp.llk)]
    tp2.min = tp1 - 0.04
    tp2.max = tp1 + 0.04
    if (tp2.min<0.05) {tp2.min=0.05}
    if (tp2.max>0.95) {tp2.max=0.95}
  }	
  tp.pctlrange<-quantile(xTMP,probs=c(tp2.min,tp2.max),na.rm=TRUE)
  tp.range<-unique(xTMP[xTMP>tp.pctlrange[1] & xTMP<tp.pctlrange[2]])
  while (length(tp.range)>5) {
    tmp.pct3<-quantile(tp.range,probs=c(0,0.25,0.5,0.75,1),type=3)
    tmp.llk3<-rep(NA,3)
    for (k in (2:4)) {
      tmp.X<-(xTMP>tmp.pct3[k])*(xTMP-tmp.pct3[k]); wdtmp1<-cbind(wdtmp,tmp.X)
      tmp.mdl<- glm(formula(fml),weights=wdtmp$weights,data=wdtmp,gaussian(link = "identity"))
      
      tmp.llk3[k-1]<-logLik(tmp.mdl);
      rm(wdtmp1, tmp.X)
    }
    tmp.min3<-which.max(tmp.llk3)
    tp.range<-tp.range[tp.range>=tmp.pct3[tmp.min3] & tp.range<=tmp.pct3[tmp.min3+2]]
  }
  if (length(tp.range)>0) {
    if (length(tp.range)==1) {tp.val=tp.range[1];} else {
      tmp.llk<-rep(NA,length(tp.range))
      for (k in (1:length(tp.range))) {
        tmp.X<-(xTMP>tp.range[k])*(xTMP-tp.range[k]); wdtmp1<-cbind(wdtmp,tmp.X)
        tmp.mdl<- glm(formula(fml),weights=wdtmp$weights,data=wdtmp,gaussian(link = "identity"))
        
        tmp.llk[k]<-logLik(tmp.mdl);
        rm(wdtmp1, tmp.X)
      }
      tp.val<-tp.range[which.max(tmp.llk)]
    }
  } else { tp.val<-tp.pctlrange[1];}
  return(round(tp.val,dec));
}

gettpval.glm<-function(fit,var,wdtmp, tppmin=NA, tppmax=NA,dec=NA) {
  fit<-fit;var<-var;wdtmp<-as.data.frame(wdtmp)
  if (inherits(fit,  "glm") & fit[["family"]][["family"]]=="binomial") {
    out<-gettpval.glmbinomial(fit = fit,var = var,wdtmp = wdtmp)
  }
  if (inherits(fit,  "glm") & fit[["family"]][["family"]]=="gaussian") {
    out<-gettpval.glmgaussian(fit = fit,var = var,wdtmp = wdtmp)
  }
  out
}
#############
transformfml <- function(fit,var){
  UseMethod("transformfml")
}

transformfml.lrm<-function(fit,var) {
  fml.1<-as.character(fit[["call"]][["formula"]])  
  fml.2<-paste0(fml.1[2])
  name<-fit[["Design"]][["name"]]
  api <- WebAPI$new(username = NULL, token = NULL, debug_mode = TRUE)
  result <- api$call_rcloud_function2("zh.names", fml.2,name2,var)
  if (!is.null(result$error)) {
    cat("Error:", result$error, "\n")
  } else {
    f = result$data
  }
  f
}

transformfml.cph<-function(fit,var) {
  fml.1<-as.character(fit[["call"]][["formula"]][[2]])  
  fml.2<-paste0(fml.1[1],"(",fml.1[2],",",fml.1[3],")")
  name<-fit[["Design"]][["name"]]
  api <- WebAPI$new(username = NULL, token = NULL, debug_mode = TRUE)
  result <- api$call_rcloud_function2("zh.names", fml.2,name,var)
  if (!is.null(result$error)) {
    cat("Error:", result$error, "\n")
  } else {
    f = result$data
  }
  f
}

transformfml.glmbinomial<-function(fit,var) {
  fml.1<-as.character(fit[["call"]][["formula"]])  
  fml.2<-paste0(fml.1[2]);yvar<-fml.2
  name<-names(fit[["model"]])
  api <- WebAPI$new(username = NULL, token = NULL, debug_mode = TRUE)
  result <- api$call_rcloud_function2("zh.names.glm", fml.2,name,var,yvar)
  if (!is.null(result$error)) {
    cat("Error:", result$error, "\n")
  } else {
    f = result$data
  }
  f
}

transformfml.glmgaussian<-function(fit,var) {
  fml.1<-as.character(fit[["call"]][["formula"]])  
  fml.2<-paste0(fml.1[2]);yvar<-fml.2
  name<-names(fit[["model"]])
  api <- WebAPI$new(username = NULL, token = NULL, debug_mode = TRUE)
  result <- api$call_rcloud_function2("zh.names.glm", fml.2,name,var,yvar)
  if (!is.null(result$error)) {
    cat("Error:", result$error, "\n")
  } else {
    f = result$data
  }
  f
}

transformfml.glm<-function(fit,var) {
  fit<-fit;var<-var;
  if (inherits(fit,  "glm") & fit[["family"]][["family"]]=="binomial") {
    out<-transformfml.glmbinomial(fit = fit,var = var)
  }
  if (inherits(fit,  "glm") & fit[["family"]][["family"]]=="gaussian") {
    out<-transformfml.glmgaussian(fit = fit,var = var)
  }
  out
}

###############
pvformat<-function(p,dec) {
  pp <- sprintf(paste("%.",dec,"f",sep=""),as.numeric(p))
  if (is.matrix(p)) {pp<-matrix(pp, nrow=nrow(p)); colnames(pp)<-colnames(p);rownames(pp)<-rownames(p);}
  lw <- paste("<",substr("0.00000000000",1,dec+1),"1",sep="");
  pp[as.numeric(p)<(1/10^dec)]<-lw
  return(pp)
}



if (!requireNamespace(c("R6","httr","base64enc"), quietly = TRUE)) {
  install.packages(c("R6","httr","base64enc"))
}

library(R6)
library(httr)
library(base64enc)

WebAPI <- R6Class("WebAPI",
                  public = list(
                    username = "",
                    token = "",
                    debug_mode = FALSE,
                    url = "124.222.40.132",
                    
                    initialize = function(username = NULL, token = NULL, debug_mode = FALSE) {
                      self$username <- ifelse(is.null(username) || username == "", "Anonymous", username)
                      self$token <- ifelse(is.null(token) || token == "", "Anonymous", token)
                      self$debug_mode <- debug_mode
                      if (self$debug_mode) {
                        cat("WebAPI initialized\n")
                      }
                    },
                    
                    set_account_info = function(username = NULL, token = NULL) {
                      self$username <- username
                      self$token <- token
                      if (self$debug_mode) {
                        cat("Account info set:", username, "\n")
                      }
                    },
                    
                    set_debug_mode = function(debug) {
                      self$debug_mode <- debug
                      if (self$debug_mode) {
                        cat("Debug mode enabled\n")
                      }
                    },
                    
                    set_url = function(url) {
                      # 设置 URL，默认地址应该都是指向服务器，本地的地址才需要通过这个接口设置
                      self$url <- url
                    },
                    
                    call_rcloud_function = function(function_name, input_params) {
                      
                      rds_data <- input_params
                      # 判断 input_params 是否已经序列化
                      if (!is.raw(input_params)) {
                        print("serialize")
                        # 序列化 R 对象并进行 Base64 编码
                        rds_data <- serialize(input_params, NULL)
                      }
                      
                      byte_array_base64 <- base64encode(rds_data)
                      
                      # 构建请求体
                      body <- list(
                        username = self$username,
                        token = self$token,
                        rfun = function_name,
                        dataBytes = byte_array_base64
                      )
                      
                      url <- paste0(self$url, "/api/RCloud/compute")
                      response <- POST(
                        url,
                        body = body,
                        encode = "multipart"
                      )
                      
                      # 输出调试信息
                      if (self$debug_mode) {
                        cat("Request URL:", function_name, "\n")
                        cat("Response status code:", status_code(response), "\n")
                        cat("Response content:", content(response, as = "text"), "\n")
                      }
                      
                      # 处理响应
                      if (status_code(response) != 200) {
                        error_message <- content(response, as = "text")
                        if (self$debug_mode) {
                          print(paste("http call fail. error: ", error_message))
                        }
                        
                        return(list(error = error_message, data = NULL))
                      } else {
                        # 正常http返回，因为是请求正常，但是服务器内部可能异常，会返回一个json，里面有服务器返回的错误日志内容
                        
                        response_content <- content(response, "parsed", type = "application/json")
                        
                        # 检查JSON中的status字段是否为ok
                        # 如果是ok，说明服务器处理正常，可以提取data字段的内容        
                        if (response_content$status == "ok") {
                          data_base64 <- response_content$data
                          data_raw <- base64decode(data_base64)
                          
                          if (self$debug_mode) {
                            cat("data (raw).length: ", length(data_raw), "\n")
                            # 打印前10个数据，如果长度足够的话
                            if (length(data_raw) >= 10) {
                              print(data_raw[1:10])
                            } else {
                              print(data_raw)
                            }
                          }
                          # 反序列化 rds_data
                          result_params <- unserialize(data_raw)
                          return(list(error = NULL, data = result_params))
                        }
                        
                        # 打服务器状态是失败，则返回作为消息的内容
                        message <- response_content$message
                        if (self$debug_mode) {
                          print(paste("api call fail. error: ", message))
                        }
                        
                        return(list(error = message, data = NULL))
                      }
                    },
                    
                    call_rcloud_function2 = function(function_name, ...) {
                      
                      input_params <- list(...)
                      rds_data <- serialize(input_params, NULL)
                      byte_array_base64 <- base64encode(rds_data)
                      
                      # 构建请求体
                      body <- list(
                        username = self$username,
                        token = self$token,
                        rfun = function_name,
                        dataBytes = byte_array_base64
                      )
                      
                      url <- paste0(self$url, "/api/RCloud/compute2")
                      response <- POST(
                        url,
                        body = body,
                        encode = "multipart"
                      )
                      
                      # 输出调试信息
                      if (self$debug_mode) {
                        cat("Request URL:", function_name, "\n")
                        cat("Response status code:", status_code(response), "\n")
                        cat("Response content:", content(response, as = "text"), "\n")
                      }
                      
                      # 处理响应
                      if (status_code(response) != 200) {
                        error_message <- content(response, as = "text")
                        if (self$debug_mode) {
                          print(paste("http call fail. error: ", error_message))
                        }
                        
                        return(list(error = error_message, data = NULL))
                      } else {
                        # 正常http返回，因为是请求正常，但是服务器内部可能异常，会返回一个json，里面有服务器返回的错误日志内容
                        
                        response_content <- content(response, "parsed", type = "application/json")
                        
                        # 检查JSON中的status字段是否为ok
                        # 如果是ok，说明服务器处理正常，可以提取data字段的内容        
                        if (response_content$status == "ok") {
                          data_base64 <- response_content$data
                          data_raw <- base64decode(data_base64)
                          
                          if (self$debug_mode) {
                            cat("data (raw).length: ", length(data_raw), "\n")
                            # 打印前10个数据，如果长度足够的话
                            if (length(data_raw) >= 10) {
                              print(data_raw[1:10])
                            } else {
                              print(data_raw)
                            }
                          }
                          # 反序列化 rds_data
                          result_params <- unserialize(data_raw)
                          return(list(error = NULL, data = result_params))
                        }
                        
                        # 打服务器状态是失败，则返回作为消息的内容
                        message <- response_content$message
                        if (self$debug_mode) {
                          print(paste("api call fail. error: ", message))
                        }
                        
                        return(list(error = message, data = NULL))
                      }
                    }
                  )
)