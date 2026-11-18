# Default: reproduce results from the original posterior parameter draws.
# Optional argument: directory of parameter draws produced by refit.R.
source('code/common.R')
args<-commandArgs(trailingOnly=TRUE); folder<-if(length(args)) args[1] else 'data/posterior'
dir.create('output',showWarnings=FALSE)
summary_rows<-list();fit_rows<-list()
for(n in model_names) {
  z<-read_draws(n,folder);d<-model_data(n)
  pars<-setdiff(names(z),c('chain','iteration'))
  summary_rows[[n]]<-data.frame(model=n,parameter=pars,mean=sapply(z[pars],mean),sd=sapply(z[pars],sd),N=nrow(d),row.names=NULL)
  ll<-loglik_matrix(z,d)
  stopifnot(all(is.finite(ll)))
  lpd<-apply(ll,2,function(x){a<-max(x);a+log(mean(exp(x-a)))})
  pointwise<- -2*(lpd-apply(ll,2,var))
  row<-data.frame(model=n,N=nrow(d),WAIC=sum(pointwise),WAIC_SE=sqrt(ncol(ll)*var(pointwise)))
  if(requireNamespace('loo',quietly=TRUE)) {
    l<-loo::loo(ll)
    row$elpd_LOO<-l$estimates['elpd_loo','Estimate'];row$elpd_LOO_SE<-l$estimates['elpd_loo','SE']
    row$max_pareto_k<-max(l$diagnostics$pareto_k)
  }
  fit_rows[[n]]<-row
}
write.csv(do.call(rbind,summary_rows),'output/posterior_summary.csv',row.names=FALSE)
write.csv(do.call(rbind,fit_rows),'output/model_fit.csv',row.names=FALSE)
if(!requireNamespace('loo',quietly=TRUE)) message('LOO skipped: install.packages("loo") and rerun to obtain elpd_LOO.')
# Table 1: preserve the original complete-case filter; correct its duplicated sti_auth column.
desc<-list()
for(s in 1:2){
  d<-complete(get(paste0('study',s)),c('t0_gen','t1_gen','sti_auth','sti_lr','leftright'))
  vars<-intersect(c('t0_gen','t1_gen','t2_gen','leftright','sti_auth','sti_lr'),names(d))
  desc[[s]]<-do.call(rbind,lapply(vars,function(v){x<-d[[v]];x<-x[!is.na(x)];data.frame(study=s,variable=v,N=length(x),mean=mean(x),SD=sd(x),min=min(x),max=max(x))}))
}
write.csv(do.call(rbind,desc),'output/descriptive_statistics.csv',row.names=FALSE)
# Figures 3 and 5: same variables, all model-ready rows, with a loess smooth.
for(s in 1:2){
 d<-get(paste0('study',s));d$x<-abs(d$sti_auth-d$t0_gen);d$y<-abs(d$t1_gen-d$t0_gen)
 d<-complete(d,c('x','y'))
 png(paste0('output/figure',if(s==1) 3 else 5,'.png'),width=1000,height=900,res=140)
 plot(d$x,d$y,xlab='Distance to stimulus',ylab='Absolute opinion change',pch=ifelse(d$t1_gen%in%c(0,1),2,1))
 grid<-seq(min(d$x),max(d$x),length.out=200); fit<-loess(y~x,d,span=.75)
 lines(grid,predict(fit,newdata=data.frame(x=grid)),lwd=2);dev.off()
}
# Figure 4: posterior mean response at ideological endpoints and both moral frames.
p<-sapply(read_draws('m12',folder),mean);x<-seq(0,1,length.out=501)
png('output/figure4.png',width=1000,height=900,res=140)
plot(x,x,type='n',xlab='Pre-stimulus opinion',ylab='Predicted post-stimulus opinion',ylim=c(0,1));abline(0,1,col='grey')
k<-0
for(b1 in 0:1)for(b2 in 0:1){k<-k+1;d<-data.frame(t0_gen=x,sti_auth=0,b1=b1,b2=b2);lines(x,pmin(1,pmax(0,predict_mu(as.list(p),d))),lty=k,lwd=2)}
legend('topleft',c('Moral 0 / left','Moral 0 / right','Moral 1 / left','Moral 1 / right'),lty=1:4,bty='n');dev.off()
cat('Analysis outputs written to output/.\n')
