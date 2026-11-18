# Run entry points from the package root. No files outside this directory are read.
study1 <- read.csv('data/study1.csv')
study2 <- read.csv('data/study2.csv')
complete <- function(x, vars) x[complete.cases(x[,vars,drop=FALSE]),,drop=FALSE]
model_data <- function(name) {
  if (substr(name,2,2)=='1') {
    d <- complete(study1,c('t0_gen','t1_gen','sti_auth'))
    if(grepl('ref|lag',name)) d<-complete(d,'t2_gen')
    d$b1<-d$stimulus_morality
    # Preserve the original reference/lagged calls, which used stimulus direction.
    d$b2<-if(grepl('ref|lag',name)) d$stimulus_direction else d$leftright/100
  } else {
    d<-complete(study2,c('t0_gen','t1_gen','sti_auth','sti_lr','leftright'))
    d<-d[ifelse(d$topic=='dev',d$ExpCondition!=0,d$ExpCondition!=1),]
    d$g1<-abs(d$leftright-d$sti_lr)/100
  }
  d$y<-if(grepl('lag',name)) d$t2_gen else d$t1_gen
  d
}
model_names <- c('m10','m11','m12','m10.ref','m11.ref','m12.ref',
                 'm10.lag','m11.lag','m12.lag','m20','m21','m22')
read_draws <- function(name, folder='data/posterior') read.csv(gzfile(file.path(folder,paste0(name,'.csv.gz'))))
# Exact predictor used by the archived JAGS functions. gamma is paper beta_0.
predict_mu <- function(p,d) {
  delta<-d$sti_auth-d$t0_gen
  mu<-d$t0_gen+p[['alpha']]*delta
  if('gamma'%in%names(p)) mu<-mu-p[['gamma']]*abs(delta)*delta
  if('b_moral'%in%names(p)) mu<-mu+(p[['b_moral']]*d$b1+p[['b_right']]*d$b2+p[['b_moralright']]*d$b1*d$b2)*delta
  if('zeta'%in%names(p)) mu<-mu-p[['zeta']]*d$g1*abs(delta)*delta
  mu
}
loglik_matrix <- function(draws,d) {
  t(vapply(seq_len(nrow(draws)),function(i) {
    p<-draws[i,];mu<-predict_mu(p,d);s<-abs(p$sigma)
    dnorm(d$y,mu,s,log=TRUE)-log(pnorm(1,mu,s)-pnorm(0,mu,s))
  },numeric(nrow(d))))
}
