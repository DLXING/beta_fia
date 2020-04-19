load('NNeighb.RData')
load('data/dat.RData')
source('src/fns4plot.R',encoding = 'UTF-8')
n=999
## Fig 2 (beta ~ lat) ------
sites=sitesNN[[paste(100)]]
res=resNN[[paste(100)]]
sites$DbP=apply(res$bp,1,function(x)x[1])# Raw beta
sites$DbPses=apply(res$bp,1,function(x)(x[1]-mean(x[-1]))/sd(x[-1]))# SES beta
sites$agg=apply(res[[paste('Imor',3,sep='')]], 1, function(x)(x[1]-mean(x[-1]))/sd(x[-1])) # SES agg
beta_res=sites

win.metafile('beta_bp_lat_N100.wmf',4,6)
par(mfrow=c(2,1),mar=c(1.6,3.6,1.5,1),oma=c(2,0,0,0),mgp=c(2.3,.5,0))
# Raw beta
fm=summary(lm(DbP~LAT,beta_res))
r=round(fm$r.sq,2)
p=round(1-pf(fm$f[1],fm$f[2],fm$f[3]),3)

r=c(r,round(cor.test(beta_res$agg,beta_res$DbP)$estimate,2))
p=c(p,round(cor.test(beta_res$agg,beta_res$DbP)$p.value,3))

r=paste('=',format(r,nsmall=2),', ',sep='')
p0=p;p[p<.001]='<0.001';p[p0>=.001]=paste('=',format(p0[p0>=.001],nsmall=3),sep='');rm(p0)
par(xpd=NA)
plot(beta_res$LAT,beta_res$DbP,xlab='',ylab=expression(paste(beta, '-diversity',sep='')),las=1, cex=.5+2*(beta_res$agg-min(beta_res$agg))/(max(beta_res$agg)-min(beta_res$agg)),pch=21,bg='white',col='grey')
text(par('usr')[1], par('usr')[4], 'A',adj=c(3.5,-.1),font=2,cex=1.2)
par(xpd=F)
text(par('usr')[1], par('usr')[3], substitute(paste(italic(R^2),r,italic(P),p),list(r=r[1],p=p[1])),adj=c(-.1,-.4))

x=seq(min(beta_res$LAT),max(beta_res$LAT),by=.1)
y=fm$coef[1,1]+fm$coef[2,1]*x
lines(x,y,lwd=2)

# SES beta
fm=summary(lm(DbPses~LAT+I(LAT^2),beta_res))
r=round(fm$r.sq,2)
p=round(1-pf(fm$f[1],fm$f[2],fm$f[3]),3)

r=c(r,round(cor.test(beta_res$agg,beta_res$DbPses)$estimate,2))
p=c(p,round(cor.test(beta_res$agg,beta_res$DbPses)$p.value,3))


r=paste('=',format(r,nsmall=2),', ',sep='')
p0=p;p[p<.001]='<0.001';p[p0>=.001]=paste('=',format(p0[p0>=.001],nsmall=3),sep='');rm(p0)
par(xpd=NA)
plot(beta_res$LAT,beta_res$DbPses,xlab='Latitude (°)',ylab=expression(paste(beta, '-deviation',sep='')),las=1, cex=.5+2*(beta_res$agg-min(beta_res$agg))/(max(beta_res$agg)-min(beta_res$agg)),pch=21,bg='white',col='grey')
text(par('usr')[1], par('usr')[4], 'B',adj=c(3.5,-.1),font=2,cex=1.2)

par(xpd=F)
text(par('usr')[1], par('usr')[3], substitute(paste(italic(R^2),r,italic(P),p),list(r=r[1],p=p[1])),adj=c(-.1, -.4))
text(par('usr')[1]/2+par('usr')[2]/2, par('usr')[4], substitute(paste(italic(r)[agg],R,italic(P),p),list(R=r[2],p=p[2])),adj=c(.5,1.4))

x=seq(min(beta_res$LAT),max(beta_res$LAT),by=.1)
y=fm$coef[1,1]+fm$coef[2,1]*x+fm$coef[3,1]*x^2
lines(x,y,lwd=2)

dev.off()


## Fig 3 VP (beta_dev ~ lat + agg)---------
plot_vp_bp(sitesNN,resNN,NN=c(2,5,10,20,50,100),aggtype = 3)

## Fig 4 (range&SD env ~ lat) --------
win.metafile('Fig4.wmf',3.5,7)
par(mfrow=c(3,1),mar=c(1.6,3.6,1.5,3.6),oma=c(2,0,0,0),mgp=c(2.3,.5,0),xpd=NA)
plot(meanRangeSize/1e6~LAT,sitesNN[["100"]],las=1,xlab='',ylab=expression(paste('Mean range size (1e6  ',km^2,')',sep='')),type='n')
bdev=apply(resNN[["100"]]$bp,1,function(x)(x[1]-mean(x[-1]))/sd(x[-1]))
points(sitesNN[["100"]]$LAT,(bdev-min(bdev))/(max(bdev)-min(bdev))*(5-1.5)+1.5,pch=16,col=grey(.5))
points(meanRangeSize/1e6~LAT,sitesNN[["100"]])
axis(4,at=(seq(40,82,by=10)-min(bdev))/(max(bdev)-min(bdev))*(5-1.5)+1.5,labels = seq(40,82,by=10),las=1,col.ticks =grey(.3),col.axis=grey(.3))
text(par('usr')[1], par('usr')[4], 'A',adj=c(3.5,-.1),font=2,cex=1.2)
summary(lm(bdev~sitesNN[['100']]$meanRangeSize))
text(par('usr')[1]/2+par('usr')[2]/2, par('usr')[4], substitute(paste(italic(R)[range.size-beta[dev]]^2,r),list(r=' = 0.05')),adj=c(.5,1.5))
mtext(expression(paste(beta, '-deviation',sep='')),4,line=2.5,col=grey(.3),cex=.9)

b=sapply(plotsNN[["100"]], function(x){sd(all.plot[all.plot$plotID%in%x,'BIO15'],na.rm=T)});plot(sitesNN[["100"]]$LAT,b,xlab='',ylab='SD Precipitation Seasonality (BIO15)',las=1)
text(par('usr')[1], par('usr')[4], 'B',adj=c(3.5,-.1),font=2,cex=1.2)
summary(lm(bdev~b))
text(par('usr')[1]/2+par('usr')[2]/2, par('usr')[4], substitute(paste(italic(R)[SD.BIO15 - beta[dev]]^2,r),list(r=' = 0.14')),adj=c(.9,2))

s=sapply(plotsNN[["100"]], function(x){sd(all.plot[all.plot$plotID%in%x,'Topsoil_Sand_Fraction'],na.rm=T)});plot(sitesNN[["100"]]$LAT,s,xlab='Latitude (°)',ylab='SD Topsoil_Sand_Fraction',las=1)
text(par('usr')[1], par('usr')[4], 'C',adj=c(3.5,-.1),font=2,cex=1.2)
summary(lm(bdev~s))
text(par('usr')[1]/2+par('usr')[2]/2, par('usr')[4], substitute(paste(italic(R)[SD.soil.texture~-~beta[dev]]^2,r),list(r=' = 0.15')),adj=c(.8,2.0))

dev.off()

## Fig 5 VP (beta_dev ~ lat + SDenv)---------
plot_vp_env_sel(sitesNN,resNN,NN=c(2,5,10, 20,50,100), env = all.plot,var = colnames(all.plot)[5:40])

