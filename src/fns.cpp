// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace arma;

CharacterVector csample_char(const CharacterVector& x) {
    // same function as sample in R
    CharacterVector ret = 
        Rcpp::RcppArmadillo::sample(x, x.size(), false, NumericVector::create()) ;
    return ret ;
}

NumericMatrix table2(const CharacterVector& x, const CharacterVector& y){
    // same function as table in R
    CharacterVector xlevs =sort_unique(x);
    CharacterVector ylevs =sort_unique(y);
    IntegerVector xi = match(x, xlevs)-1;
    IntegerVector yi = match(y, ylevs)-1;
    NumericMatrix res(xlevs.size(),ylevs.size());
    for(int i=0; i<x.size(); i++)
        res(xi[i],yi[i]) +=1.0;
    rownames(res)=xlevs;
    colnames(res)=ylevs;
    return res;
}

// [[Rcpp::export]]
NumericMatrix commsimu(const DataFrame& ind){
    // the individual level null model of Kraft et al (2011). 
    //
    // ind is a dataframe with each row a individual tree. 
    // it has two columns: sp and subplot. 
    //
    // the funciton returns a species X plots composition matrix.

    CharacterVector swap = csample_char(ind("subplot"));
    NumericMatrix comp=table2(ind("sp"),swap);
    return comp;
}


// [[Rcpp::export]]
DataFrame Beta(const NumericMatrix& comp){
    // function to calculate beta diversity based on Hill numbers of 
    // orders 0,1, and 2 (Chao & Chiu 2016 MEE).
    // 
    double N=comp.ncol();
    NumericVector q=NumericVector::create(0,1,2);
    NumericVector C=clone(q),U=clone(q);
    
    NumericVector p = rowSums(comp)/sum(comp);
    double Dg, Da;
    p=p[p>0];

    NumericVector pa = comp/sum(comp);
    pa=pa[pa>0];

    // q=0
    Dg=p.size();
    Da=pa.size()/N;
    C[0]=(Dg/Da-N)/(1.0-N);
    U[0]=(N*Da/Dg-1.0)/(N-1.0);
    // q=1
    Dg=exp(-sum(p*log(p)));
    Da=exp(-sum(pa*log(pa))-log(N));
    C[1]=1-log(Dg/Da)/log(N);
    U[1]=C[1];
    // q=2
    Dg=1.0/sum(p*p);
    Da=1.0/sum(pa*pa)/N;
    C[2]=(N*Da/Dg-1.0)/(N-1.0);
    U[2]=(Dg/Da-N)/(1.0-N);

    DataFrame res=DataFrame::create(_("q")=q, _("DSorensen")= 1-C, _("DJaccard")=1-U);
    return res;
}

// [[Rcpp::export]]
NumericVector morisita(NumericMatrix mat){
    // function to calculate the mean morisita index for:
    // all non-singleton species (res1),
    // the most abundant species (res2),
    // common species (#ind >= #plots; res3),
    // and dominate species that accumulate the first 50% of total abundance (res4).
    int n=mat.ncol(), m=mat.nrow();
    arma::mat M(mat.begin(),m,n,false);//copy mat to M (same mem)
    vec abund = sum(M,1);
    double ab_tot=sum(abund);
    uvec Ind = find(abund>1);
    if(Ind.size()<1){
        return NumericVector::create(NA_REAL, NA_REAL, NA_REAL, NA_REAL);
    } else{
        abund=abund(Ind);
        arma::mat M2 = M.rows(Ind);
        M2 = M2 % M2;
        vec res = n*(sum(M2,1) - abund) / (abund % abund - abund);
        double res1=mean(res), res2=NA_REAL, res3=NA_REAL, res4=NA_REAL;
        
        Ind=abund.index_max();
        res2 = mean(res(Ind));
        
        Ind=find(abund>=n);
        if(Ind.size()>=1){
            res3=mean(res(Ind));
        }
        
        uvec ab_order=sort_index(abund);
        vec ab_sort = abund(ab_order);
        vec res_sort = res(ab_order);
        vec ab_cum = cumsum(ab_sort);
        Ind=find(ab_cum>=(0.5*ab_tot));
        if(Ind.size()>=1){
            res4=mean(res_sort(Ind));
        }
        return NumericVector::create(res1,res2,res3,res4);
    }
}



// [[Rcpp::export]]
double betaP(NumericMatrix comp){
    // the proportional species turnover
    int N = comp.ncol(), S = comp.nrow();
    arma::mat comp1(comp.begin(),S,N,true);//copy comp to comp1 (diff mem)
    comp1.elem( find(comp1 > 0) ).ones();
    rowvec alpha=sum(comp1,0);
    double b=(double)S;
    b/= mean(alpha);
    return 1-1/b;
}

// [[Rcpp::export]]
List do_simu(List comp_dat, List ind_dat, int nsim=999){
    // The main function for calculating beta diversity and aggregation
    // and their standardized effect sizes.
    //
    // comp_dat is the list of metacommunities, each of which is a 
    // species X plots matrix containing the composition data.
    //
    // ind_dat is the list of individual level data converted from the comp_dat.
    //
    int m = comp_dat.size();
    arma::mat Imor(m,nsim+1), Imor2(m,nsim+1), Imor3(m,nsim+1), 
    Imor4(m,nsim+1), bp(m,nsim+1);
    arma::cube bS2(4,m,nsim+1), bJ2(4,m,nsim+1);
    Imor.fill(0.0);
    Imor2.fill(0.0);
    Imor3.fill(0.0);
    Imor4.fill(0.0);
    bp.fill(0.0);
    bS2.fill(0.0);
    bJ2.fill(0.0);

    NumericVector q=NumericVector::create(0,1,2);
    
    for(int i=0; i<m; i++){
        NumericMatrix comp=comp_dat[i];
        DataFrame ind_i(ind_dat[i]);
        NumericVector mor=morisita(comp);
        Imor(i,0)+=mor[0];
        Imor2(i,0)+=mor[1];
        Imor3(i,0)+=mor[2];
        Imor4(i,0)+=mor[3];
        bp(i,0)+=betaP(comp);
        DataFrame btmp=Beta(comp);
        NumericVector DS =btmp("DSorensen");
        NumericVector DJ =btmp("DJaccard");
        for(int j=0; j<q.size(); j++){
            bS2(j,i,0)+=DS[j];
            bJ2(j,i,0)+=DJ[j];
            
        }
        // run the simulation
        for(int k=1; k<nsim+1; k++){
            comp=commsimu(ind_i);
            mor = morisita(comp);
            Imor(i,k)+=mor[0];
            Imor2(i,k)+=mor[1];
            Imor3(i,k)+=mor[2];
            Imor4(i,k)+=mor[3];
            bp(i,k)+=betaP(comp);
            btmp=Beta(comp);
            DS =btmp("DSorensen");
            DJ =btmp("DJaccard");
            for(int j=0; j<q.size(); j++){
                bS2(j,i,k)+=DS[j];
                bJ2(j,i,k)+=DJ[j];
                
            }
            
        }
        if(i%20==1)
            Rcout<<i<<"\t"; //printf("%i\t",  i); //
    }
    printf("\n");
    return List::create(_("Imor")=Imor, _("Imor2")=Imor2, 
                        _("Imor3")=Imor3, _("Imor4")=Imor4,_("bp")=bp, 
                        _("bS2")=bS2, _("bJ2")=bJ2);
}

