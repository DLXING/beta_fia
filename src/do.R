load('data/FIA.RData')
range_size=read.csv('range_size_updated.csv')
## function to generate the metacommunities for a sample of locations
############################
NN_smpl=function(plot.dat, tree.dat, range.dat, NN=100){
    require(sp)
    require(e1071)
    sites0=plot.dat[,1:3]
    sites0$gridID=paste(sites0$LON%/%1,sites0$LAT%/%1)
    sites=sites0[0,]
    set.seed(123)
    # sample the focal plots (i.e. centres of the metacommunities) - one from each 1*1 grid
    for(id in unique(sites0$gridID)){
        tmp=subset(sites0,gridID==id)
        tmp=tmp[sample(nrow(tmp),1),]
        sites=rbind(sites,tmp)
    }
    tree.dat$subplot=sites0$plotID[match(tree.dat$subplot,sites0$plotID)]
    
    # generate a metacommunity that include the first NN neigbouring plots
    sites$nspp=NA
    sites$meanRangeSize=NA
    sites$meanLATextent=NA
    comp.dat=list()
    nn_plots=list()
    pairwise_dists=list()
    del=c()
    
    for(i in 1:nrow(sites)){
        loc_i=as.matrix(sites[i,c('LON','LAT')])
        neigb_i=sites0[which(abs(sites0[,'LON']-loc_i[1])<=5 & abs(sites0[,'LAT']-loc_i[2])<=5),]
        dist_i=spDistsN1(as.matrix(neigb_i[,c('LON','LAT')]), loc_i, longlat = T)
        names(dist_i)=neigb_i$plotID
        dist_i=dist_i[!duplicated(dist_i)]
        
        pneigb_i=names(dist_i[order(dist_i)][1:NN])
        nn_plots[[i]]=pneigb_i
        
        if(any(is.na(pneigb_i)))
            del=c(del,i)
        else{
            neigb_i=neigb_i[which(neigb_i$plotID%in%pneigb_i),]
            dist_i=spDists(as.matrix(neigb_i[,c('LON','LAT')]), longlat = T)
            x=dist_i[lower.tri(dist_i)]
            pairwise_dists[[i]]=c(min=min(x),max=max(x),mean=mean(x),sd=sd(x),sk=skewness(x),kt=kurtosis(x))
            ind_i=tree.dat[which(tree.dat[,'subplot']%in%pneigb_i),]
            abund_i=as.data.frame.matrix(with(ind_i,table(sp,subplot)))
            abund_i=abund_i[,as.character(pneigb_i)]
            sites$nspp[i]=nrow(abund_i)
            sites$meanRangeSize[i]=mean(range.dat$areaKM2[range.dat$Latin.Name%in%rownames(abund_i)],na.rm=T)
            sites$meanLATextent[i]=mean(range.dat$LATextent[range.dat$Latin.Name%in%rownames(abund_i)],na.rm = T)
            
            comp.dat[[i]]=as.matrix(abund_i)
        }
        if(i%%20==1)
            cat(i,'\t')
    }
    sel=1:nrow(sites)
    sel=sel[!sel%in%del]
    comp.dat=comp.dat[sel]
    nn_plots=nn_plots[sel]
    pairwise_dists=pairwise_dists[sel]
    sites=sites[sel,]
    
    return(list(comp.dat=comp.dat, sites=sites, nn_plots=nn_plots, pairwise_dists=pairwise_dists))
}
## function to convert the species X plots composition matrixs to individual level data frames that have a single tree in each row
############################
mat_to_ind=function(comp.dat){
    for(i in 1:length(comp.dat)){
        mat=comp.dat[[i]]
        spnames=rownames(mat);abund=rowSums(mat)
        sp=rep(spnames,abund)
        plotnames=colnames(mat)
        subplot=NULL
        for(j in 1:nrow(mat))
            subplot=c(subplot,rep(plotnames,mat[j,]))
        comp.dat[[i]]=data.frame(sp=sp,subplot=subplot)
    }
    comp.dat
}

############################
## run the analyses
library(Rcpp)
sourceCpp('src/fns.cpp')
sitesNN=list()
plotsNN=list()
pwdistsNN=list()
resNN=list()
NN=c(2,5,10,20,50,100)
for(nn in NN){
    cat('metacommunity size:',nn,' - generating the metacommunities...\n')
    dat=NN_smpl(all.plot,all.tree,range_size,nn)
    cat('do the analyses...\n')
    sitesNN[[paste(nn)]]=dat$sites
    plotsNN[[paste(nn)]]=dat$nn_plots
    pwdistsNN[[paste(nn)]]=dat$pairwise_dists
    ind.dat=mat_to_ind(dat$comp.dat)
    resNN[[paste(nn)]]=do_simu(dat$comp.dat,ind.dat,nsim=999)
}
save(NN,sitesNN,plotsNN,pwdistsNN,resNN,file = 'NNeighb.RData')
