# Refit the 12 models underlying Tables 2-5. Requires JAGS 4 and runjags.
if(!requireNamespace('runjags',quietly=TRUE)) stop('Install JAGS 4 and the R package runjags first.')
library(runjags)
source('code/common.R')
dir.create('output/refit',recursive=TRUE,showWarnings=FALSE)
for(s in 1:2) {
  source(paste0('code/models_study',s,'.R'))
  for(n in model_names[substr(model_names,2,2)==as.character(s)]) {
    d<-model_data(n);k<-as.integer(substr(n,3,3))
    args<-list(y_cens=d$y,x1=d$t0_gen,x2=d$sti_auth)
    if(k==2) args<-c(args,if(s==1) list(b1=d$b1,b2=d$b2) else list(g1=d$g1))
    fit<-do.call(get(paste0('bayes_m',k)),args)
    draws<-do.call(rbind,lapply(seq_along(fit$mcmc),function(j){
      m<-as.matrix(fit$mcmc[[j]]); m<-m[,!grepl('^loglik',colnames(m)),drop=FALSE]
      data.frame(chain=j,iteration=seq_len(nrow(m)),m)
    }))
    con<-gzfile(paste0('output/refit/',n,'.csv.gz'),'wt');write.csv(draws,con,row.names=FALSE);close(con)
  }
}
cat('Refits saved. Run: Rscript code/analyze.R output/refit\n')
