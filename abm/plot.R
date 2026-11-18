# Dependency-free plotting of ABM output. Run from package root after abm/run.py.
files<-list.files('output',pattern='^abm_.*[.]csv$',full.names=TRUE)
for(f in files){
 d<-read.csv(f);legacy<-grepl('illustration',f)
 metrics<-c('spread','dispersion',if(legacy) 'coverage_legacy' else 'coverage_interval')
 png(sub('[.]csv$','.png',f),width=1500,height=500,res=140);par(mfrow=c(1,3),mar=c(4,4,1,1))
 for(v in metrics){
  xvar<-if(legacy) 'gamma' else 'bubble'
  plot(range(d[[xvar]]),c(0,1),type='n',xlab=if(legacy) 'Homophily (gamma)' else 'Bubble size',ylab=v)
  groups<-if(legacy) unique(d$bubble) else 0
  for(g in seq_along(groups)){
   z<-if(legacy) d[d$bubble==groups[g],] else d
   means<-aggregate(z[[v]],list(x=z[[xvar]]),mean)
   lines(means[[1]],means[[2]],lty=g,lwd=2)
  }
 }
 dev.off()
}
