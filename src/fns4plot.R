plot_vp_bp=function(sitesNN,resNN,NN,aggtype=NULL){
    resList=list()
    for(nn in NN){
        sites=sitesNN[[paste(nn)]]
        res=resNN[[paste(nn)]]
        Dorder=apply(res$bp,1,function(x)(x[1]-mean(x[-1]))/sd(x[-1]))
        for(i in 1:3){
            Dorder=cbind(Dorder,apply(res$bS2[i,,1:n],1,function(x)(x[1]-mean(x[-1]))/sd(x[-1])))
            Dorder=cbind(Dorder,apply(res$bJ2[i,,1:n],1,function(x)(x[1]-mean(x[-1]))/sd(x[-1])))
        }
        colnames(Dorder)=c('DbP','DS0','DJ0','DS1','DJ1','DS2','DJ2')
        sites=cbind(sites,Dorder)
        sites$agg=apply(res[[paste('Imor',aggtype,sep='')]], 1, function(x)(x[1]-mean(x[-1]))/sd(x[-1]))
        sites=sites[sites$DbP!=Inf,]
        resList[[paste(nn)]]=sites
    }
    require(vegan)
    win.metafile(paste('varpart_agg',aggtype,'_bp.wmf',sep=''),5,5)
    #layout(matrix(c(1,2,3,4,1,5,6,7),nrow=2,byrow=T), widths= c(1.7,1,1,1))
    par(mar=c(2,1,1.5,2),oma=c(3,3,.1,0),mgp=c(1.9,.5,0),xpd=NA)
    
    vp=sapply(resList,function(x){
        varpart(x[,'DbP'],~x$agg,~x$LAT+I(x$LAT^2))$part$indfract$Adj.R.squared
    })
    barplot(vp*100,las=1,legend.text = paste('[',letters[1:4],']',sep=''),args.legend = list(bg='white'),xpd = NA,lwd=.5)
    #title(expression(paste("Proportional species turnover (",beta['P'],")")))
    
    mtext('Spatial extent of regional community\n(Number of local plots)',1,line=1,outer=T)
    mtext('Proportion of variation (%)', side=2, line=1, outer=T)
    dev.off()
}

plot_vp_env_sel=function(sitesNN,resNN,NN,env=all.plot,var=c('Topsoil_Sand_Fraction','BIO7')){
    resList=list()
    for(nn in NN){
        sites=sitesNN[[paste(nn)]]
        res=resNN[[paste(nn)]]
        plots=plotsNN[[paste(nn)]]
        Dorder=apply(res$bp,1,function(x)(x[1]-mean(x[-1]))/sd(x[-1]))
        for(i in 1:3){
            Dorder=cbind(Dorder,apply(res$bS2[i,,1:n],1,function(x)(x[1]-mean(x[-1]))/sd(x[-1])))
            Dorder=cbind(Dorder,apply(res$bJ2[i,,1:n],1,function(x)(x[1]-mean(x[-1]))/sd(x[-1])))
        }
        colnames(Dorder)=c('DbP','DS0','DJ0','DS1','DJ1','DS2','DJ2')
        sites=cbind(sites,Dorder)
        ENV=sapply(plots, function(x){sd(env[env$plotID%in%x,var[1]],na.rm=T)})
        if(length(var)>1){
            for(i in 2:length(var)){
                ENV=cbind(ENV,sapply(plots, function(x){sd(env[env$plotID%in%x,var[i]],na.rm=T)}))
            }
        }
        colnames(ENV)=var
        sites=cbind(sites,ENV)
        sites=sites[sites$DbP!=Inf,]
        resList[[paste(nn)]]=sites
    }
    require(vegan)
    require(packfor)
    win.metafile('varpart_env_bp.wmf',5,5)
    par(mar=c(2,1,1.5,2),oma=c(3,3,.1,0),mgp=c(1.9,.5,0),xpd=NA)
    
    vp=sapply(resList,function(x){
        x=na.omit(x[,c('DbP',var,'LAT')])
        var.sel=forward.sel(x[,'DbP'],x[,var])$variables
        tmp=paste('x$',var.sel,sep='')
        envForm=paste('~',tmp[1],sep='')
        if(length(var.sel)>1){
            for(j in 2:length(var.sel)){
                envForm=paste(envForm,tmp[j],sep='+')
            }
        }
        
        varpart(x[,'DbP'],as.formula(envForm),~x$LAT+I(x$LAT^2))$part$indfract$Adj.R.squared
    })
    barplot(vp*100,las=1,legend.text = paste('[',letters[1:4],']',sep=''),args.legend = list(bg='white'),xpd = NA,lwd=.5)
    
    
    mtext('Spatial extent of regional community\n(Number of local plots)',1,line=1,outer=T)
    mtext('Proportion of variation (%)', side=2, line=1, outer=T)
    dev.off()
    #############
    win.metafile(paste('varpart_env_Hill.wmf',sep=''),6,4)
    par(mfrow=c(2,3),mar=c(2,1,1.5,.1),oma=c(3,3,.1,0),mgp=c(1.9,.5,0),xpd=NA)
    q=c('= 0','= 1','= 2','= 0','= 1','= 2')
    index=c('DS0','DS1','DS2','DJ0','DJ1','DJ2')
    for(i in 1:6){
        vp=sapply(resList,function(x){
            x=na.omit(x[,c(index[i],var,'LAT')])
            var.sel=forward.sel(x[,index[i]],x[,var])$variables
            tmp=paste('x$',var.sel,sep='')
            envForm=paste('~',tmp[1],sep='')
            if(length(var.sel)>1){
                for(j in 2:length(var.sel)){
                    envForm=paste(envForm,tmp[j],sep='+')
                }
            }
            varpart(x[,index[i]],as.formula(envForm),~x$LAT+I(x$LAT^2))$part$indfract$Adj.R.squared
        })
        if(i==1)
            barplot(vp*100,las=2,las=2, xaxt=if(i%in%c(4:6)) 's' else 'n', yaxt = if(i%in%c(1,4))'s'else 'n', xpd=NA,lwd=.5,legend.text = paste('[',letters[1:4],']',sep=''),args.legend = list(bg='white'))
        else
            barplot(vp*100,las=2,las=2, xaxt=if(i%in%c(4:6)) 's' else 'n', yaxt = if(i%in%c(1,4))'s'else 'n', xpd=NA,lwd=.5)
        if(i%in%c(1:3))title(main = substitute(paste('Sørensen, ',italic(q),Q),list(Q=q[i])))
        if(i%in%c(4:6))title(main = substitute(paste('Jaccard, ',italic(q),Q),list(Q=q[i])))
    }
    
    mtext('Spatial extent of regional community (Number of local plots)',1,line=1,outer=T)
    mtext('Proportion of variation (%)', side=2, line=1, outer=T)
    dev.off()
}
