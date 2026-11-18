library(ggplot2)

# setting the theme and creating labels ----

theme_set(theme_minimal() +
            theme(axis.text=element_text(size=11,color='black'),
                  axis.title=element_text(size=11),
                  plot.margin = unit(c(1, 1, 0.5, 0.5), "lines"),
                  panel.border = element_rect(fill=NA,color="black", size=1,linetype="solid"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
)

labs <- c(expression(paste(gamma,' = ',0.5)),
          expression(paste(gamma,' = ',1)),
          expression(paste(gamma,' = ',2)),
          expression(paste(gamma,' = ',4)))

# influence functions and weights ----

fw    <- function(x,y=0,gamma=1){ 1 - (gamma * abs(y - x)) }
fw_bc <- function(x,y=0,gamma=1){ (1 - abs(y - x)) ^ gamma }

cf    <- function(x, y=0, alpha=1, gamma=1){ 
  x2 = x + alpha * fw(x,y,gamma) * (y-x) 
  if(x2>1){ x2 = 1 }
  if(x2<0){ x2 = 0 }
  return(x2)
}
cf_bc <- function(x, y=0, alpha=1, gamma=1){ 
  x2 = x + alpha * fw_bc(x,y,gamma) * (y-x) 
  if(x2>1){ x2 = 1 }
  if(x2<0){ x2 = 0 }
  return(x2)
}

    # the function used in Bayesian models
ncf   <- function(x, y=0, alpha=1, gamma=1, beta=NA){ 
  if(is.na(beta)){ beta = gamma * alpha }  # maybe not true for alpha<1
  x2 = x + alpha * (y - x) - beta * abs(y - x) * (y - x)
  if(x2>1){ x2 = 1 }
  if(x2<0){ x2 = 0 }
  return(x2)
}

# plotting opinion change ----

  # WLI

df <- expand.grid( x1 = c(0:1000) / 1000 ,
                   gamma = c(0.5,1,2,4),
                   x2 = NA)
for(i in 1:nrow(df)) df$x2[i] <- cf(df$x1[i], gamma=df$gamma[i])
df$gamma <- as.factor(df$gamma)

(p.change.wli <- ggplot(df,aes(x=x1, y=x2, linetype=gamma, color=gamma)) + 
  geom_abline(intercept=0, slope=1, lty=1, color='grey50', size=.5) + 
  geom_line(size=.8) +
  scale_x_continuous(breaks=c(0:5)/5, limits=c(0,1), name=expression(o['i'])) + 
  scale_y_continuous(breaks=c(0:5)/5, limits=c(0,1), name=expression(o['i,t+1'])) + 
  scale_color_grey("",labels=labs) +
  scale_linetype_manual("",labels=labs,values=c(1,2,3,4)) +
  theme(
    legend.title = element_blank(),
    legend.position = 'none', #c(.2,.75),
    legend.text = element_text(size=11),#, margin = margin(b=-10, unit='mm')),
    legend.text.align = 0,
    legend.key.width = unit(9,'mm')
    ))



  # BC

df <- expand.grid( x1 = c(0:1000) / 1000 ,
                   gamma = c(0.5,1,2,4),
                   x2 = NA)
for(i in 1:nrow(df)) df$x2[i] <- cf_bc(df$x1[i], gamma=df$gamma[i])
df$gamma <- as.factor(df$gamma)

(p.change.bc <- ggplot(df,aes(x=x1, y=x2, linetype=gamma, color=gamma)) + 
  geom_abline(intercept=0, slope=1, lty=1, color='grey50', size=.5) + 
  geom_line(size=.8) +
  scale_x_continuous(breaks=c(0:5)/5, limits=c(0,1), name=expression(o['i'])) + 
  scale_y_continuous(breaks=c(0:5)/5, limits=c(0,1), name=expression(o['i,t+1'])) + 
  scale_color_grey("",labels=labs) +
  scale_linetype_manual("",labels=labs,values=c(1,2,3,4)) +
  theme(
    legend.title = element_blank(),
    legend.position = 'none', #c(.2,.75),
    legend.text = element_text(size=11),
    legend.text.align = 0,
    legend.key.width = unit(9,'mm')
  ))


## plotting the weights ----

  # WLI

df <- expand.grid( x1 = c(0:1000) / 1000 ,
                   gamma = c(0.5,1,2,4),
                   x2 = NA )
for(i in 1:nrow(df)) df$x2[i] <- fw(df$x1[i], gamma=df$gamma[i])
df$gamma <- as.factor(df$gamma)

(p.weight.wli <- ggplot(df,aes(x=x1, y=x2, linetype=gamma, color=gamma)) + 
  geom_abline(intercept=0, slope=0, lty=1, color='grey50', size=.5) + 
  geom_line(size=.6) +
  scale_x_continuous(breaks=c(0:5)/5, limits=c(0,1), name=expression(abs(o['j'] - o['i,t']))) + 
  scale_y_continuous(breaks=c(-3:1), limits=c(-3,1), name=expression(w['ijt'])) + 
  scale_color_grey("",labels=labs) +
  scale_linetype_manual("",labels=labs,values=c(1,2,3,4)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(.23,.27),
    legend.text = element_text(size=11),
    legend.text.align = 0,
    legend.key.width = unit(7,'mm')
  ))


  # BC

df <- expand.grid( x1 = c(0:1000) / 1000 ,
                   gamma = c(0.5,1,2,4),
                   x2 = NA )
for(i in 1:nrow(df)) df$x2[i] <- fw_bc(df$x1[i], gamma=df$gamma[i])
df$gamma <- as.factor(df$gamma)

(p.weight.bc <- ggplot(df,aes(x=x1, y=x2, linetype=gamma, color=gamma)) + 
    geom_abline(intercept=0, slope=0, lty=1, color='grey50', size=.5) + 
    geom_line(size=.6) +
    scale_x_continuous(breaks=c(0:5)/5, limits=c(0,1), name=expression(abs(o['j'] - o['i,t']))) + 
    scale_y_continuous(breaks=c(-3:1), limits=c(-3,1), name=expression(w['ijt'])) + 
    scale_color_grey("",labels=labs) +
    scale_linetype_manual("",labels=labs,values=c(1,2,3,4)) +
    theme(
      legend.title = element_blank(),
      legend.position = c(.23,.27),
      legend.text = element_text(size=11),
      legend.text.align = 0,
      legend.key.width = unit(7,'mm')
    ))




# saving the plots ----

ggsave('output/influence_predict_wli.eps',p.change.wli,width=35*2,height=32*2,units=c('mm'))
ggsave('output/influence_predict_bc.eps', p.change.bc, width=35*2,height=32*2,units=c('mm'))
ggsave('output/influence_weights_wli.eps',p.weight.wli,width=35*2,height=32*2,units=c('mm'))
ggsave('output/influence_weights_bc.eps', p.weight.bc, width=35*2,height=32*2,units=c('mm'))

