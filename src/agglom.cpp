/*
Agglomerative partitioning of LC-MS peaks
author: Martin Loos, mloos@envibee.ch
Copyright (c) 2019 envibee. All rights reserved.
*/

#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

using namespace Rcpp;

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RMATRIX2(m,i,j) (INTEGER(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RVECTOR(m,i) (REAL(m)[i])
#define RRow(m) (INTEGER(GET_DIM(m))[0])
#define RCol(m) (INTEGER(GET_DIM(m))[1])



/************************************************************************/
/* EIC clustering *******************************************************/
/************************************************************************/
// [[Rcpp::export]]
SEXP getEIC_new(      
	SEXP mz,
	SEXP RT,
	SEXP intens,
	SEXP orderedint,
	SEXP orderedret,
	SEXP dmzdens,
	SEXP ppm2,
	SEXP drtdens,
	SEXP merged2
){

           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(intens = AS_NUMERIC(intens));
           PROTECT(orderedint = AS_INTEGER(orderedint));
           PROTECT(orderedret = AS_INTEGER(orderedret));
           PROTECT(dmzdens = AS_NUMERIC(dmzdens));
           PROTECT(ppm2 = AS_INTEGER(ppm2));
           PROTECT(drtdens = AS_NUMERIC(drtdens));
           PROTECT(merged2 = AS_INTEGER(merged2));
           double *ret, *mass, *intensity;
           mass = NUMERIC_POINTER(mz);
           ret = NUMERIC_POINTER(RT);
           intensity = NUMERIC_POINTER(intens);
           int *ordint, *ordret;
           ordint = INTEGER_POINTER(orderedint);
           ordret = INTEGER_POINTER(orderedret);
           double dmzdens2 = NUMERIC_VALUE(dmzdens);
           int ppm3 = INTEGER_VALUE(ppm2);
           double drtdens2 = NUMERIC_VALUE(drtdens);
           int merged3 = INTEGER_VALUE(merged2);
           int leng = LENGTH(RT);
           int m,n,i,k=0,clustnumb,maxat=0,maxit=0;
           double delmz;
           SEXP clusters;
           PROTECT(clusters = Rf_allocMatrix(REALSXP, leng, 13));
           double *clus;
           clus = REAL(clusters);
           for(m=0;m<13;m++){
               for(n=0;n<leng;n++){
                   clus[(m*leng)+n]=0;
               }
           }
           int *at;
           at = new int[leng];

           /* initialize with most intense measurement ************************/
           clustnumb=1;
           if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint)-1)))/1e6);}else{delmz=dmzdens2;}
           clus[0]=(*(mass+(*(ordint)-1))-(2*delmz));       /* low mass boundary **************/
           clus[(1*leng)]=(*(mass+(*(ordint)-1))+(2*delmz));/* high mass boundary *************/
           clus[(2*leng)]=(*(ret+(*(ordint)-1))-drtdens2);  /* low RT boundary ****************/
           clus[(3*leng)]=(*(ret+(*(ordint)-1))+drtdens2);  /* high RT boundary ***************/
           clus[(4*leng)]=1;                                /* number of measurements *********/
           clus[(5*leng)]=*(mass+(*(ordint)-1));            /* mass sum ***********************/
           clus[(6*leng)]=clustnumb;                        /* cluster ID *********************/
           clus[(7*leng)]=0;                                /* merged (1) or not (0)? *********/
           clus[(8*leng)]=*(intensity+(*(ordint)-1));       /* maximum intensity in a cluster */
           clus[(9*leng)+(*(ordint)-1)]=clustnumb;          /* cluster ID for measurement *****/
           clus[(10*leng)]=0;                               /* variance, set later if merged **/
           clus[(11*leng)]=*(mass+(*(ordint)-1));           /* lowest mass in cluster *********/
           clus[(12*leng)]=*(mass+(*(ordint)-1));           /* highest mass in cluster ********/

           /* assign all other peaks ******************************************/
           for(n=1;n<leng;n++){
               /* check for possible fit to existing clusters *****************/
               maxat=0;
               for(m=0;m<clustnumb;m++){
                if(*(mass+(*(ordint+n)-1))<=clus[(1*leng)+m]){
                   if(*(mass+(*(ordint+n)-1))>=clus[(0*leng)+m]){
                           if(*(ret+(*(ordint+n)-1))>=clus[(2*leng)+m]){
                               if(*(ret+(*(ordint+n)-1))<=clus[(3*leng)+m]){
                                       at[maxat]=m;
                                       maxat++;
                               }
                           }
                       }
                   }
               }
               /* (a) not assignable - create new cluster *********************/
               if(maxat==0){
                   clustnumb++;
                   clus[(9*leng)+(*(ordint+n)-1)]=clustnumb;
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   clus[0+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))-(2*delmz));
                   clus[(1*leng)+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))+(2*delmz));
                   clus[(2*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))-drtdens2);
                   clus[(3*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))+drtdens2);
                   clus[(4*leng)+(clustnumb-1)]=1;
                   clus[(5*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(6*leng)+(clustnumb-1)]=clustnumb;
                   clus[(7*leng)+(clustnumb-1)]=0;
                   clus[(8*leng)+(clustnumb-1)]=*(intensity+(*(ordint+n)-1));
                   clus[(11*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(12*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   continue;
               }
               /* check for other peaks at same RT ****************************/
               maxit=maxat;
               for(i=0;i<leng;i++){ /* find RT-order index */
                   if(*(ordint+n)==*(ordret+i)){
                       k=i;
                       break;
                   }
               }

               if(k>0){ /* backward over RT-order */
                   for(i=(k-1);i>=0;i--){
                       if(  *(ret+(*(ordret+k)-1))==*(ret+(*(ordret+i)-1))  ){
                           if((clus[(9*leng)+(*(ordret+i)-1)])==0){
                               continue;
                           }else{
                               for(m=0;m<maxat;m++){
                                   if(clus[(9*leng)+(*(ordret+i)-1)]==(at[m]+1)){
                                       at[m]=-9999;
                                       maxit--;
                                       break;
                                   }
                               }
                           }
                       }else{
                          break;
                       }
                   }
               }

               if(k<(leng-1)){ /* forward over RT-order */
                   for(i=(k+1);i<leng;i++){
                       if(*(ret+(*(ordret+k)-1))==*(ret+(*(ordret+i)-1))){
                           if((clus[(9*leng)+(*(ordret+i)-1)])==0){
                               continue;
                           }else{
                               for(m=0;m<maxat;m++){
                                   if(clus[(9*leng)+(*(ordret+i)-1)]==(at[m]+1)){
                                       at[m]=-9999;
                                       maxit--;
                                       break;
                                   }
                               }
                           }
                       }else{
                           break;
                       }
                   }
               }

               if(maxit<maxat){
                   i=maxat;
                   k=0;
                   for(m=0;m<i;m++){
                       if(at[m]==-9999){
                           k++;
                           maxat--;
                       }else{
                           at[m-k]=at[m];
                       }
                   }
               }
               maxat=maxit;
               /* (a) not assignable - create new cluster *********************/
               if(maxat==0){
                   clustnumb++;
                   clus[(9*leng)+(*(ordint+n)-1)]=clustnumb;
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   clus[0+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))-(2*delmz));
                   clus[(1*leng)+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))+(2*delmz));
                   clus[(2*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))-drtdens2);
                   clus[(3*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))+drtdens2);
                   clus[(4*leng)+(clustnumb-1)]=1;
                   clus[(5*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(6*leng)+(clustnumb-1)]=clustnumb;
                   clus[(7*leng)+(clustnumb-1)]=0;
                   clus[(8*leng)+(clustnumb-1)]=*(intensity+(*(ordint+n)-1));
                   clus[(11*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(12*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   continue;
               }
               /* (b) fits exactly to one cluster *****************************/
               if(maxat==1){
                   clus[(9*leng)+(*(ordint+n)-1)]=(at[0]+1);
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   /* shrink lower & upper mass bounds */
                   if(clus[(0*leng)+(at[0])]<(*(mass+(*(ordint+n)-1))-(2*delmz))){clus[(0*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))-(2*delmz));}
                   if(clus[(1*leng)+(at[0])]>(*(mass+(*(ordint+n)-1))+(2*delmz))){clus[(1*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))+(2*delmz));}
                   /* adapt lower & upper RT bounds */
                   if(clus[(2*leng)+(at[0])]>(*(ret+(*(ordint+n)-1))-drtdens2)){clus[(2*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))-drtdens2);}
                   if(clus[(3*leng)+(at[0])]<(*(ret+(*(ordint+n)-1))+drtdens2)){clus[(3*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))+drtdens2);}
                   clus[(4*leng)+(at[0])]=clus[(4*leng)+(at[0])]+1;
                   clus[(5*leng)+(at[0])]=clus[(5*leng)+(at[0])]+*(mass+(*(ordint+n)-1));
                   if(*(mass+(*(ordint+n)-1))<clus[(11*leng)+(at[0])]){clus[(11*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   if(*(mass+(*(ordint+n)-1))>clus[(12*leng)+(at[0])]){clus[(12*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   continue;
               }
               /* (c) fits to several clusters - m/z density clustering *******/
               if(maxat>1){
                   delmz=fabs(*(mass+(*(ordint+n)-1))- ( clus[(5*leng)+(at[0])]/clus[(4*leng)+(at[0])] ));
                   for(m=0;m<maxat;m++){
                       if(
                           fabs(*(mass+(*(ordint+n)-1))- ( clus[(5*leng)+(at[m])]/clus[(4*leng)+(at[m])] ))<
                           delmz
                       ){
                           at[0]=at[m];
                           delmz=fabs(*(mass+(*(ordint+n)-1))- ( clus[(5*leng)+(at[0])]/clus[(4*leng)+(at[0])] ));
                       }
                   }
                   clus[(9*leng)+(*(ordint+n)-1)]=(at[0]+1);
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   /* shrink lower & upper mass bounds */
                   if(clus[(0*leng)+(at[0])]<(*(mass+(*(ordint+n)-1))-(2*delmz))){clus[(0*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))-(2*delmz));}
                   if(clus[(1*leng)+(at[0])]>(*(mass+(*(ordint+n)-1))+(2*delmz))){clus[(1*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))+(2*delmz));}
                   /* adapt lower & upper RT bounds */
                   if(clus[(2*leng)+(at[0])]>(*(ret+(*(ordint+n)-1))-drtdens2)){clus[(2*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))-drtdens2);}
                   if(clus[(3*leng)+(at[0])]<(*(ret+(*(ordint+n)-1))+drtdens2)){clus[(3*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))+drtdens2);}
                   clus[(4*leng)+(at[0])]=clus[(4*leng)+(at[0])]+1;
                   clus[(5*leng)+(at[0])]=clus[(5*leng)+(at[0])]+*(mass+(*(ordint+n)-1));
                   if(*(mass+(*(ordint+n)-1))<clus[(11*leng)+(at[0])]){clus[(11*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   if(*(mass+(*(ordint+n)-1))>clus[(12*leng)+(at[0])]){clus[(12*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   continue;
              }
           }

           /* *************************************************************** */
           /* merge isobaric mass cluster *************************************/

		   if((merged3==1)&(clustnumb>1)){

               /* set mean ****************************************************/
               for(n=0;n<clustnumb;n++){
                   clus[(5*leng)+(n)]=(clus[(5*leng)+(n)]/clus[(4*leng)+(n)]);
               }
               /* define exclusion-list of cluster with same RT ***************/
               std::vector<int> blackone;
               std::vector<int> blacktwo;
               int blacksize;
               for(n=1;n<leng;n++){
                   for(k=(n-1);k>=0;k--){
                       if(*(ret+(*(ordret+k)-1))==*(ret+(*(ordret+n)-1))){
                           if( /* confine to overlapping masses */
                                (clus[(0*leng)+(int(clus[(9*leng)+(*(ordret+n)-1)])-1)] < clus[(1*leng)+(int(clus[(9*leng)+(*(ordret+k)-1)])-1)]) &
                                (clus[(1*leng)+(int(clus[(9*leng)+(*(ordret+n)-1)])-1)] > clus[(0*leng)+(int(clus[(9*leng)+(*(ordret+k)-1)])-1)])
                           ){          /* smallest index always first */
                               if((clus[(9*leng)+(*(ordret+k)-1)])<(clus[(9*leng)+(*(ordret+n)-1)])){
                                       blackone.push_back(int(clus[(9*leng)+(*(ordret+k)-1)]));
                                       blacktwo.push_back(int(clus[(9*leng)+(*(ordret+n)-1)]));
                               }else{
                                       blackone.push_back(int(clus[(9*leng)+(*(ordret+n)-1)]));
                                       blacktwo.push_back(int(clus[(9*leng)+(*(ordret+k)-1)]));
                               }
                           }
                       }else{
                           break;
                       }
                   }
               }
               blacksize=blackone.size();
               /* find mergeable clusters: indices & mass differences *********/
               /* cross-check for overlaps in RT between cluster **************/
               std::vector<int> clusterone;
               std::vector<int> clustertwo;
               std::vector<double> clusterdiff;
               int doit;
               for(n=0;n<(clustnumb-1);n++){
                   for(m=(n+1);m<clustnumb;m++){
                       if(
                           (clus[(0*leng)+(n)]<clus[(1*leng)+(m)]) &
                           (clus[(1*leng)+(n)]>clus[(0*leng)+(m)])
                       ){ /* check for overlapping masses, ppm window */
                           if(
                               (
                                   (clus[(0*leng)+(n)]<=clus[(11*leng)+(m)]) &
                                   (clus[(1*leng)+(n)]>=clus[(12*leng)+(m)])
                                )||(
                                   (clus[(0*leng)+(m)]<=clus[(11*leng)+(n)]) &
                                   (clus[(1*leng)+(m)]>=clus[(12*leng)+(n)])
                                )
                           ){ /* check for overlapping masses, lowest/highest centroid masses vs. ppm window */
                               if(blacksize==0){
                                   clusterone.push_back(n+1);
                                   clustertwo.push_back(m+1);
                                   clusterdiff.push_back(fabs(clus[(5*leng)+(n)]-clus[(5*leng)+(m)])); // difference of mean mass
                               }else{
                                   doit=0;
                                   for(k=0;k<blacksize;k++){ // check exclusion list
                                       if(blackone[k]==(n+1)){
                                            if(blacktwo[k]==(m+1)){
                                                doit=1;
                                                break;
                                            }
                                       }
                                   }
                                   if(doit==0){
                                       clusterone.push_back(n+1);
                                       clustertwo.push_back(m+1);
                                       clusterdiff.push_back(fabs(clus[(5*leng)+(n)]-clus[(5*leng)+(m)]));
                                    }
                               }
                           }
                       }
                   }
                }
                /* merge cluster ***********************************************/
                int mergesize=clusterone.size(), stay, gone;
                std::vector<int> clusterase;
				if(mergesize>0){				
					while(mergesize>0){
                        /* find cluster pair closest in mean mass difference */
                        delmz=clusterdiff[0];
                        m=0;
						if(mergesize>1){
                            for(k=1;k<mergesize;k++){
                               if(clusterdiff[k]<delmz){
                                   m=k;
                                   delmz=clusterdiff[k];
                               }
                            }
                    	}
                        /* reassign measurements */
						for(n=0;n<leng;n++){
                            if(clus[(9*leng)+n]==clustertwo[m]){
                                clus[(9*leng)+n]=clusterone[m];
                            }
                        }
						/* update cluster entries */
                        clus[(4*leng)+(clusterone[m]-1)]=(clus[(4*leng)+(clusterone[m]-1)]+clus[(4*leng)+(clustertwo[m]-1)]); /* number of measurements */
                        clus[(6*leng)+(clustertwo[m]-1)]=clusterone[m]; /* cluster ID */
                        clus[(7*leng)+(clustertwo[m]-1)]=1;            /* merged? */
                        /* delete that pair & all links to second (=merged) cluster & check blacklist  */
                        gone=clustertwo[m];
                        stay=clusterone[m];
                        mergesize=clusterone.size();
                        n=0;
                        while(n<mergesize){ /* remove merges */
                            if( (clusterone[n]==gone) || (clustertwo[n]==gone) ){
                                clusterone.erase(clusterone.begin()+n);
                                clustertwo.erase(clustertwo.begin()+n);
                                clusterdiff.erase(clusterdiff.begin()+n);
                                mergesize=clusterone.size();
                            }else{
                                n++;
                            }
                        }
                        if((blacksize>0) & (mergesize>0)){
                            clusterase.clear(); /* check blacklist-entries */
                            for(k=0;k<blacksize;k++){
                                if(blackone[k]==gone){
                                    clusterase.push_back(blacktwo[k]);
                                }
                                 if(blacktwo[k]==gone){
                                    clusterase.push_back(blackone[k]);
                                }
                            }
                            if(clusterase.size()>0){ /* remove indirect blacklisted links */
                                for(m=0;(unsigned)m<clusterase.size();m++){
                                    mergesize=clusterone.size();
                                    n=0;
                                    while(n<mergesize){
                                        if( (clustertwo[n]==stay) || (clusterone[n]==clusterase[m])){
                                            clusterone.erase(clusterone.begin()+n);
                                            clustertwo.erase(clustertwo.begin()+n);
                                            clusterdiff.erase(clusterdiff.begin()+n);
                                            mergesize=clusterone.size();
                                        }else{
                                            n++;
                                        }
                                    }
                                }
                            }
                        }
                   }
               }
               /* output ******************************************************/
               delete[] at;
               UNPROTECT(10);
               return(clusters);
            }else{
               /* output ******************************************************/
               delete[] at;
               UNPROTECT(10);
               return(clusters);
           }
       }



/******************************************************************************/
/* calculate theta ************************************************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP series_relat(
    SEXP homol_peaks_relat,
    SEXP range_mz,
    SEXP range_RT
){

    PROTECT(homol_peaks_relat = AS_NUMERIC(homol_peaks_relat));
    PROTECT(range_mz = AS_NUMERIC(range_mz));
    PROTECT(range_RT = AS_NUMERIC(range_RT));

    int len1 = (RRow(homol_peaks_relat)-1);
    int i = 0, j = 0, n = 0, m = 0;
    double aa = 0, ab = 0, a2 = 0, ba = 0 , bb = 0, b2 = 0, c2 = 0, to_cos, wink;

    SEXP thetas = PROTECT(Rf_allocVector(REALSXP, (len1 + 1)));
    for(i = 0;i <= len1; i++){
        RVECTOR(thetas,i) = R_PosInf;
    }

    for(i = 0; i < len1; i++){
        for(j = i; j <= len1; j++){
            if( RMATRIX(homol_peaks_relat, j, 0) != RMATRIX(homol_peaks_relat, i, 0) ){
                j--;
                break;
            }
        }
        if((j-i)> 0 ){
            void R_CheckUserInterrupt(void);
            for(m = i; m < j; m++){
                aa = (RMATRIX(homol_peaks_relat,m,2) / RVECTOR(range_mz,0));
                ba = (RMATRIX(homol_peaks_relat,m,3) / RVECTOR(range_RT,0));
                for(n = (m + 1); n <= j; n++){
                    ab = (RMATRIX(homol_peaks_relat,n,2) / RVECTOR(range_mz,0));
                    bb = (RMATRIX(homol_peaks_relat,n,3) / RVECTOR(range_RT,0));
                    a2 = ((aa * aa) + (ab * ab));
                    b2 = ((ba * ba) + (bb * bb));
                    c2 = sqrt(a2 * b2);
                    to_cos = (  ((aa * ba) + (ab * bb))   /   c2  );
                    if(ISNA(to_cos) || ISNAN(to_cos)){
                        wink = 0;
                    }else{
                        wink = acos(to_cos);
                        wink = (wink * 180 / M_PI);
                    }
                    if(!R_FINITE(RVECTOR(thetas, n)) || wink<abs(RVECTOR(thetas, n))){
                        RVECTOR(thetas,n) = wink;
                    }
                }
            }
            i = j;
        }
    }

    UNPROTECT(4);
    return thetas;

}


/******************************************************************************/
/* calculate moving counts ****************************************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP moving_count(
    SEXP homol_peaks_relat,
    SEXP deldel
){

    PROTECT(homol_peaks_relat = AS_NUMERIC(homol_peaks_relat));
    PROTECT(deldel = AS_NUMERIC(deldel));
    int len = RRow(homol_peaks_relat);
    int i=0,n=0;

    SEXP counts = PROTECT(Rf_allocVector(REALSXP,len));
    for(i = 0; i < len; i++){
        RVECTOR(counts, i) = 1;
    }

    for(i = 0; i < (len - 1); i++){
        void R_CheckUserInterrupt(void);
        for(n = (i + 1); n < len; n++){
            if((RMATRIX(homol_peaks_relat, n, 2) - RMATRIX(homol_peaks_relat, i, 2)) <= RVECTOR(deldel, 0)){
                RVECTOR(counts, i)++;
                RVECTOR(counts, n)++;
            }else{
                break;
            }
        }
    }

    UNPROTECT(3);
    return counts;

}


/******************************************************************************/
/* compare row-wise ***********************************************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP Ccompare(
    SEXP a,
    SEXP b,
    SEXP results
){
        PROTECT(a = AS_NUMERIC(a));
        PROTECT(b = AS_NUMERIC(b));
        PROTECT(results = AS_NUMERIC(results));
        int dims, len1, len2, i = 0, m = 0, at = 0, k = 0;
		bool stopit;
        len1 = RRow(a);
        len2 = RRow(b);
        dims = RCol(a);

		for(i = 0; i < len1; i++){
			stopit = false;
			for(m = 0; m < dims; m++){
				while( RMATRIX(b, at, m) < RMATRIX(a, i, m) ){
					at++;
					if(at >= len2){
						stopit = true;
						break;
					}
					if(m > 0){
						for(k = 0; k < m; k++){
							if(RMATRIX(b, at, k) != RMATRIX(a, i, k)){
								stopit = true;
								break;							
							}
						}
						if(stopit){break;}
					}	
				}
				if(stopit){break;}
				if(RMATRIX(b, at, m) > RMATRIX(a, i, m)){break;}	
				if(m == (dims - 1)){RMATRIX(results, i, 0) = (at + 1);}
			}
			if(at >= len2){break;}
		}
					
        UNPROTECT(3);
        return(results);
}


/******************************************************************************/
/* find peak intensities & correct them ***************************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP correct_intens(
                        SEXP corfac,
                        SEXP sampleID,
                        SEXP intens,
                        SEXP sampleID_peak
){


            PROTECT(corfac = AS_NUMERIC(corfac));
            PROTECT(sampleID = AS_INTEGER(sampleID));
            PROTECT(intens = AS_NUMERIC(intens));
            PROTECT(sampleID_peak = AS_INTEGER(sampleID_peak));
            double *corfac2;
            corfac2 = NUMERIC_POINTER(corfac);
            int *sampleID2;
            sampleID2 = INTEGER_POINTER(sampleID);
            double *intens2;
            intens2 = NUMERIC_POINTER(intens);
            int *sampleID_peak2;
            sampleID_peak2 = INTEGER_POINTER(sampleID_peak);

            int leng_sampleID = LENGTH(sampleID);
            int leng_peak = LENGTH(intens);
            int n,m;

            SEXP new_intens;        /* stores corrected intensities*/
            PROTECT(new_intens = NEW_NUMERIC(leng_peak));
            double *atnew;
            atnew = NUMERIC_POINTER(new_intens);
            for(n=0;n<leng_peak;n++){*(atnew+n) = 0;}

            for(n=0;n<leng_peak;n++){
                for(m=0;m<leng_sampleID;m++){
                    if(*(sampleID_peak2+n)==*(sampleID2+m)){
                        *(atnew+n)=*(intens2+n)/ *(corfac2+m);
                        continue;
                    }
                }
            }

            UNPROTECT(5);
            return new_intens;
}


/******************************************************************************/
/* build groups from profile-profile relation *********************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP metagroup(
                            SEXP proffrom, /* must be sorted */
                            SEXP profto
                            ){

            PROTECT(proffrom = AS_INTEGER(proffrom));
            PROTECT(profto = AS_INTEGER(profto));
            int *proffrom2;
            proffrom2 = INTEGER_POINTER(proffrom);
            int *profto2;
            profto2 = INTEGER_POINTER(profto);
            int leng = LENGTH(proffrom);
            int n, m, g, k, atwhat, atto_n, atfrom_n, stay;

            SEXP group;
            PROTECT(group = NEW_INTEGER(leng));
            int *grouped;
            grouped = INTEGER_POINTER(group);
            for(n = 0; n < leng; n++){*(grouped+n) = 0;}
            SEXP from;
            PROTECT(from = NEW_INTEGER(leng));
            int *atfrom;    /* stores indices */
            atfrom = INTEGER_POINTER(from);
            for(n = 0; n < leng; n++){*(atfrom + n) = 0;}
            SEXP to = 0;        /* stores profile IDs */
            PROTECT(to = NEW_INTEGER(leng));
            int *atto;
            atto = INTEGER_POINTER(to);
            for(n = 0; n < leng; n++){*(atto+n) = 0;}

            g = 1;
            for(n = 0; n < leng; n++){
                if(*(grouped+n) == 0){
                    *(grouped+n) = g;
                    /* initialize  */
                    *(atto)=*(proffrom2+n); /* cause it could be there repeatedly */
                    *(atto+1)=*(profto2+n);
                    atto_n=2;
                    stay=1;
                    while(stay>0){
                        /* to -> from */
                        /* given profile IDs in to, search indices in from having these IDs */
                        atfrom_n=0;
                        atwhat=n;
                        for(m=0;m<atto_n;m++){
                           if(*(atto+m)>=*(proffrom2+atwhat)){
                               for(k=atwhat;k<leng;k++){
                                    if( *(proffrom2+k) > *(atto+m) ){
                                        atwhat=k;
                                        break;
                                    }
                                    if( *(proffrom2+k) == *(atto+m) ){
                                        if(*(grouped+k)==0){
                                            *(grouped+k)=g;
                                            *(atfrom+atfrom_n)=k;
                                            atfrom_n++;
                                            atwhat=k;
                                        }
                                    }
                                }
                            }else{
                               for(k=atwhat;k>n;k--){
                                    if( *(proffrom2+k) < *(atto+m) ){
                                        atwhat=k;
                                        break;
                                    }
                                    if( *(proffrom2+k) == *(atto+m) ){
                                        if(*(grouped+k)==0){
                                            *(grouped+k)=g;
                                            *(atfrom+atfrom_n)=k;
                                            atfrom_n++;
                                            atwhat=k;
                                        }
                                    }
                                }
                            }
                        }
                        /* from -> to */
                        /* write with indices in from profile IDs into to */
                        atto_n=0;
                        if(atfrom_n>0){
                            for(k=0;k<atfrom_n;k++){
                                *(atto+atto_n)=*(profto2+*(atfrom+k));
                                atto_n++;
                            }
                        }else{
                            stay=0;
                        }
                    }
                    g++;
                }
            }

            UNPROTECT(5);
            return group;

}


/******************************************************************************/
/* count peaks in neigbourhood of dmz and dRT *********************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP profpeakprof(
                            SEXP ProPeak_pro,
                            SEXP ProPeak_peak,
                            SEXP Peak_peak1,
                            SEXP Peak_peak2,
                            SEXP Peak_score,
                            SEXP PeakPro
                             ){

            PROTECT(ProPeak_pro = AS_INTEGER(ProPeak_pro));
            PROTECT(ProPeak_peak = AS_INTEGER(ProPeak_peak));
            PROTECT(Peak_peak1 = AS_INTEGER(Peak_peak1));
            PROTECT(Peak_peak2 = AS_INTEGER(Peak_peak2));
            PROTECT(Peak_score = AS_NUMERIC(Peak_score));
            PROTECT(PeakPro = AS_INTEGER(PeakPro));
            int leng1 = LENGTH(ProPeak_pro);
            int leng2 = LENGTH(Peak_peak1);
            int n,m,atpeak;

            int *ProPeak_pro2;
            ProPeak_pro2 = INTEGER_POINTER(ProPeak_pro);
            int *ProPeak_peak2;
            ProPeak_peak2 = INTEGER_POINTER(ProPeak_peak);
            int *Peak_peak12;
            Peak_peak12 = INTEGER_POINTER(Peak_peak1);
            int *Peak_peak22;
            Peak_peak22 = INTEGER_POINTER(Peak_peak2);
            double *Peak_score2;
            Peak_score2 = NUMERIC_POINTER(Peak_score);
            int *PeakPro2;
            PeakPro2 = INTEGER_POINTER(PeakPro);

            SEXP relations;
            PROTECT(relations = Rf_allocMatrix(REALSXP, leng2, 3));
            double *relat;
            relat = REAL(relations);
            for(m=0;m<3;m++){
               for(n=0;n<leng2;n++){
                   relat[(m*leng2)+n]=0;
               }
            }

            atpeak=0;
            for(m=0;m<leng1;m++){
                for(n=atpeak;n<leng2;n++){
                    if(*(Peak_peak12+n)>*(ProPeak_peak2+m)){
                        atpeak=n;
                        break;
                    }else{
                        if(*(Peak_peak12+n)==*(ProPeak_peak2+m)){
                            relat[n]=*(ProPeak_pro2+m);
                            relat[(1*leng2)+n]=*(PeakPro2+(*(Peak_peak22+n)-1));
                            relat[(2*leng2)+n]=*(Peak_score2+n);
                        }
                    }
                }
            }


            UNPROTECT(7);
            return relations;

}


/******************************************************************************/
/* merge profiles **********************************************************/
/******************************************************************************/
// [[Rcpp::export]]
      SEXP mergeProfiles(
                             SEXP mz_lower,
                             SEXP mz_upper,
                             SEXP RT_lower,
                             SEXP RT_upper,
                             SEXP intens,
                             SEXP sam,
                             SEXP orderedint,
                             SEXP orderedsam,
                             SEXP supress
                       ){

           PROTECT(mz_lower = AS_NUMERIC(mz_lower));
           PROTECT(mz_upper = AS_NUMERIC(mz_upper));
           PROTECT(RT_lower = AS_NUMERIC(RT_lower));
           PROTECT(RT_upper = AS_NUMERIC(RT_upper));
           PROTECT(intens = AS_NUMERIC(intens));
           PROTECT(sam = AS_INTEGER(sam));
           PROTECT(orderedint = AS_INTEGER(orderedint));
           PROTECT(orderedsam = AS_INTEGER(orderedsam));
           PROTECT(supress = AS_INTEGER(supress));

           double *ret_lower, *ret_upper, *mass_lower, *mass_upper, *intensity;
           mass_lower = NUMERIC_POINTER(mz_lower);
           mass_upper = NUMERIC_POINTER(mz_upper);
           ret_lower = NUMERIC_POINTER(RT_lower);
           ret_upper = NUMERIC_POINTER(RT_upper);
           intensity = NUMERIC_POINTER(intens);
           int *samp, *ordint, *ordsam;
           samp = INTEGER_POINTER(sam);
           ordint = INTEGER_POINTER(orderedint);
           ordsam = INTEGER_POINTER(orderedsam);
           int supr = INTEGER_VALUE(supress);
           int leng = LENGTH(RT_lower);
           int m,n,i,k=0,clustnumb,maxat=0,maxit=0;

           SEXP clusters;
           PROTECT(clusters = Rf_allocMatrix(REALSXP, leng, 9));
           double *clus;
           clus = REAL(clusters);
           for(m=0;m<13;m++){
               for(n=0;n<leng;n++){
                   clus[(m*leng)+n]=0;
               }
           }
           int *at;
           at = new int[leng];

           /* initialize with most intense measurement ************************/
           clustnumb=1;
           clus[0]=(*(mass_lower+(*(ordint)-1)));       /* low mass boundary **************/
           clus[(1*leng)]=(*(mass_upper+(*(ordint)-1)));/* high mass boundary *************/
           clus[(2*leng)]=(*(ret_lower+(*(ordint)-1)));  /* low RT boundary ****************/
           clus[(3*leng)]=(*(ret_upper+(*(ordint)-1)));  /* high RT boundary ***************/
           clus[(4*leng)]=1;                                /* number of profiles *********/
           clus[(5*leng)]=clustnumb;                        /* cluster ID *********************/
           clus[(6*leng)]=0;                                /* merged (1) or not (0)? *********/
           clus[(7*leng)]=*(intensity+(*(ordint)-1));       /* maximum intensity in a cluster */
           clus[(8*leng)+(*(ordint)-1)]=clustnumb;          /* cluster ID for measurement *****/

           /* assign all other peaks ******************************************/
           for(n=1;n<leng;n++){
               /* check for possible fit to existing clusters *****************/
               maxat=0;
               for(m=0;m<clustnumb;m++){
                if(*(mass_upper+(*(ordint+n)-1))>=clus[(0*leng)+m]){
                   if(*(mass_lower+(*(ordint+n)-1))<=clus[(1*leng)+m]){
                           if(*(ret_lower+(*(ordint+n)-1))>=clus[(2*leng)+m]){
                               if(*(ret_lower+(*(ordint+n)-1))<=clus[(3*leng)+m]){
                                    at[maxat]=m;
                                    maxat++;
                               }
                           }
                       }
                   }
               }
               /* (a) not assignable - create new cluster *********************/
               if(maxat==0){
                   clustnumb++;
                   clus[(8*leng)+(*(ordint+n)-1)]=clustnumb;
                   clus[0+(clustnumb-1)]=(*(mass_lower+(*(ordint+n)-1)));
                   clus[(1*leng)+(clustnumb-1)]=(*(mass_upper+(*(ordint+n)-1)));
                   clus[(2*leng)+(clustnumb-1)]=(*(ret_lower+(*(ordint+n)-1)));
                   clus[(3*leng)+(clustnumb-1)]=(*(ret_upper+(*(ordint+n)-1)));
                   clus[(4*leng)+(clustnumb-1)]=1;
                   clus[(5*leng)+(clustnumb-1)]=clustnumb;
                   clus[(6*leng)+(clustnumb-1)]=0;
                   clus[(7*leng)+(clustnumb-1)]=*(intensity+(*(ordint+n)-1));
                   continue;
               }
               /* check for other peaks in same sample ****************************/
               maxit=maxat;
               for(i=0;i<leng;i++){ /* find sample-order index */
                   if(*(ordint+n)==*(ordsam+i)){
                       k=i;
                       break;
                   }
               }
               if(k>0){ /* backward over sample-order */
                   for(i=(k-1);i>=0;i--){
                       if(  *(samp+(*(ordsam+k)-1))==*(samp+(*(ordsam+i)-1))  ){
                           if((clus[(8*leng)+(*(ordsam+i)-1)])==0){
                               continue;
                           }else{
                               for(m=0;m<maxat;m++){
                                   if(clus[(8*leng)+(*(ordsam+i)-1)]==(at[m]+1)){
                                       at[m]=-9999;
                                       maxit--;
                                       break;
                                   }
                               }
                           }
                       }else{
                           break;
                       }
                   }
               }
               if(k<(leng-1)){ /* forward over sample-order */
                   for(i=(k+1);i<leng;i++){
                       if(*(samp+(*(ordsam+k)-1))==*(samp+(*(ordsam+i)-1))){
                           if((clus[(8*leng)+(*(ordsam+i)-1)])==0){
                               continue;
                           }else{
                               for(m=0;m<maxat;m++){
                                   if(clus[(8*leng)+(*(ordsam+i)-1)]==(at[m]+1)){
                                       at[m]=-9999;
                                       maxit--;
                                       break;
                                   }
                               }
                           }
                       }else{
                           break;
                       }
                   }
               }
               if(maxit<maxat){
                   i=maxat;
                   k=0;
                   for(m=0;m<i;m++){
                       if(at[m]==-9999){
                           k++;
                           maxat--;
                       }else{
                           at[m-k]=at[m];
                       }
                   }
               }
               /* (a) not assignable - create new cluster *********************/
               if( (maxit==0) & (supr!=1)){
                   clustnumb++;
                   clus[(8*leng)+(*(ordint+n)-1)]=clustnumb;
                   clus[0+(clustnumb-1)]=(*(mass_lower+(*(ordint+n)-1)));
                   clus[(1*leng)+(clustnumb-1)]=(*(mass_upper+(*(ordint+n)-1)));
                   clus[(2*leng)+(clustnumb-1)]=(*(ret_lower+(*(ordint+n)-1)));
                   clus[(3*leng)+(clustnumb-1)]=(*(ret_upper+(*(ordint+n)-1)));
                   clus[(4*leng)+(clustnumb-1)]=1;
                   clus[(5*leng)+(clustnumb-1)]=clustnumb;
                   clus[(6*leng)+(clustnumb-1)]=0;
                   clus[(7*leng)+(clustnumb-1)]=*(intensity+(*(ordint+n)-1));
                   continue;
               }
               /* (b) fits to one or several cluster **************************/
               /* always assign to the most intense cluster at[0] *************/
               if(maxat>0){
                    clus[(8*leng)+(*(ordint+n)-1)]=(at[0]+1);
                    if(clus[(0*leng)+(at[0])]<(*(mass_lower+(*(ordint+n)-1)))){clus[(0*leng)+(at[0])]=(*(mass_lower+(*(ordint+n)-1)));}
                    if(clus[(1*leng)+(at[0])]>(*(mass_upper+(*(ordint+n)-1)))){clus[(1*leng)+(at[0])]=(*(mass_upper+(*(ordint+n)-1)));}
                    if(clus[(2*leng)+(at[0])]<(*(ret_lower+(*(ordint+n)-1)))){clus[(2*leng)+(at[0])]=(*(ret_lower+(*(ordint+n)-1)));}
                    if(clus[(3*leng)+(at[0])]>(*(ret_upper+(*(ordint+n)-1)))){clus[(3*leng)+(at[0])]=(*(ret_upper+(*(ordint+n)-1)));}
                    clus[(5*leng)+(at[0])]=clus[(5*leng)+(at[0])]+1;
                    clus[(6*leng)+(clustnumb-1)]=1;
               }
           }
           /* output **********************************************************/
           delete[] at;
           UNPROTECT(10);
           return(clusters);

       }



/******************************************************************************/
/* count peaks in neigbourhood of dmz and dRT *********************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP neighbour(        SEXP mz, /* must be sorted */
                       SEXP rt,
                       SEXP sample,
                       SEXP maxsample,
                       SEXP ppm,
                       SEXP dmz,
                       SEXP drt
                      ){


            PROTECT(mz = AS_NUMERIC(mz));
            PROTECT(rt = AS_NUMERIC(rt));
            PROTECT(sample = AS_INTEGER(sample));
            PROTECT(maxsample = AS_INTEGER(maxsample));
            PROTECT(ppm = AS_INTEGER(ppm));
            PROTECT(dmz = AS_NUMERIC(dmz));
            PROTECT(drt = AS_NUMERIC(drt));
            double *mass;
            mass = NUMERIC_POINTER(mz);
            double *ret;
            ret = NUMERIC_POINTER(rt);
            int *sam;
            sam = INTEGER_POINTER(sample);
            int ppm2 = INTEGER_VALUE(ppm);
            int maxsample2 = INTEGER_VALUE(maxsample);
            double dmass = NUMERIC_VALUE(dmz);
            double dret = NUMERIC_VALUE(drt);
            int n,p,m;
            double lowmass, lowret, highret;
            int leng = LENGTH(rt);
            SEXP clusters;
            PROTECT(clusters = Rf_allocMatrix(INTSXP, maxsample2, maxsample2));
            int *clus;
            clus = INTEGER(clusters);
            for(m=0;m<maxsample2;m++){
               for(n=0;n<maxsample2;n++){
                   clus[(m*maxsample2)+n]=0;
               }
            }

            for(n=1;n<leng;n++){
                if(ppm2 == 1){
                    lowmass=(*(mass+n)-((*(mass+n)*dmass)/1E6));
                }else{
                    lowmass=(*(mass+n)-(dmass/1000));
                }
                lowret=(*(ret+n)-dret);
                highret=(*(ret+n)+dret);
                for(p=(n-1);p>=0;p--){ /* only towards lower mass */
                    if((*(mass+p))<=(lowmass)){
                        break;
                    }else{
                        if((*(ret+p)>=lowret)&&(*(ret+p)<=highret)){
                            clus[((*(sam+n)-1)*maxsample2)+(*(sam+p)-1)]=(clus[((*(sam+n)-1)*maxsample2)+(*(sam+p)-1)]+1);
                            clus[((*(sam+p)-1)*maxsample2)+(*(sam+n)-1)]=(clus[((*(sam+p)-1)*maxsample2)+(*(sam+n)-1)]+1);
                        }
                    }
                }
            }

            UNPROTECT(8);
            return clusters;

}

/******************************************************************************/
/* agglomerative partitioning of all data based on dmz and dRT ****************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP Cagglom(          SEXP mz, /* must be sorted */
                       SEXP rt,
                       SEXP sample,
                       SEXP ppm,
                       SEXP dmz,
                       SEXP drt
                      ){


            PROTECT(mz = AS_NUMERIC(mz));
            PROTECT(rt = AS_NUMERIC(rt));
            PROTECT(sample = AS_INTEGER(sample));
            PROTECT(ppm = AS_INTEGER(ppm));
            PROTECT(dmz = AS_NUMERIC(dmz));
            PROTECT(drt = AS_NUMERIC(drt));
            double *mass;
            mass = NUMERIC_POINTER(mz);
            double *ret;
            ret = NUMERIC_POINTER(rt);
            int *sam;
            sam = INTEGER_POINTER(sample);
            int ppm2 = INTEGER_VALUE(ppm);
            double dmass = NUMERIC_VALUE(dmz);
            double dret = NUMERIC_VALUE(drt);
            int n,m,p,k=0,found=0;
            double lowmass, highmass, lowret, highret;
            int leng = LENGTH(rt);
            SEXP outit;
            PROTECT(outit = NEW_INTEGER(leng));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n=0;n<leng;n++){*(at+n) = 0;}
            SETLENGTH(outit,leng);
            int *these;
            these = new int[leng];
            int *those;
            those = new int[leng];
            int untilthese=0,untilthose=0,atthese;

            for(n=0;n<leng;n++){
                if(*(at+n) == 0){ /* unassigned entry ? */
                    k++;
                    *(at+n)=k;
                    found=1;
                    these[0]=n;
                    untilthese=1;
                    atthese=1;
                    while(found == 1){
                        found=0;
                        if(atthese==1){
                            untilthose=0;
                            for(m=0;m<untilthese;m++){
                                if(ppm2 == 1){
                                   lowmass=(*(mass+these[m])-((*(mass+these[m])*dmass)/1E6));
                                   highmass=(*(mass+these[m])+((*(mass+these[m])*dmass)/1E6));
                                }else{
                                   lowmass=(*(mass+these[m])-(dmass/1000));
                                   highmass=(*(mass+these[m])+(dmass/1000));
                                }
                                lowret=(*(ret+these[m])-dret);
                                highret=(*(ret+these[m])+dret);
                                if(these[m]>0){
                                    for(p=(these[m]-1);p>=0;p--){ /* towards lower mass */
                                        if((*(mass+p))<=(lowmass)){
                                            break;
                                        }else{
                                            if(*(at+p)==0){
                                                if(*(sam+these[m])!=*(sam+p)){
                                                    if((*(ret+p)>=lowret)&&(*(ret+p)<=highret)){
                                                        those[untilthose]=p;
                                                        untilthose++;
                                                        *(at+p)=k;
                                                        found=1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if(these[m]<(leng-1)){
                                    for(p=(these[m]+1);p<leng;p++){ /* towards higher mass */
                                        if((*(mass+p))>=(highmass)){
                                            break;
                                        }else{
                                            if(*(at+p)==0){
                                                if(*(sam+these[m])!=*(sam+p)){
                                                    if((*(ret+p)>=lowret)&&(*(ret+p)<=highret)){
                                                        those[untilthose]=p;
                                                        untilthose++;
                                                        *(at+p)=k;
                                                        found=1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            atthese=0;
                        }else{
                            untilthese=0;
                            for(m=0;m<untilthose;m++){
                                if(ppm2 == 1){
                                   lowmass=(*(mass+those[m])-((*(mass+those[m])*dmass)/1E6));
                                   highmass=(*(mass+those[m])+((*(mass+those[m])*dmass)/1E6));
                                }else{
                                   lowmass=(*(mass+those[m])-(dmass/1000));
                                   highmass=(*(mass+those[m])+(dmass/1000));
                                }
                                lowret=(*(ret+those[m])-dret);
                                highret=(*(ret+those[m])+dret);
                                if(those[m]>0){
                                    for(p=(those[m]-1);p>=0;p--){
                                        if((*(mass+p))<=(lowmass)){
                                            break;
                                        }else{
                                            if(*(at+p)==0){
                                                if(*(sam+those[m])!=*(sam+p)){
                                                    if((*(ret+p)>=lowret)&(*(ret+p)<=highret)){
                                                        these[untilthese]=p;
                                                        untilthese++;
                                                        *(at+p)=k;
                                                        found=1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                if(those[m]<(leng-1)){
                                    for(p=(those[m]+1);p<leng;p++){
                                       if((*(mass+p))>=(highmass)){
                                           break;
                                        }else{
                                             if(*(at+p)==0){
                                                if(*(sam+those[m])!=*(sam+p)){
                                                    if((*(ret+p)>=lowret)&(*(ret+p)<=highret)){
                                                        these[untilthese]=p;
                                                        untilthese++;
                                                        *(at+p)=k;
                                                        found=1;
                                                    }
                                                }
                                             }
                                        }
                                    }
                                }
                            }
                            atthese=1;
                        }
                    }
                }
            }

            delete[] these;
            delete[] those;
            SETLENGTH(outit,leng);
            UNPROTECT(7);
            return outit;

}

/******************************************************************************/
/* assemble indices for return value of Cagglom *******************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP indexed(      SEXP index, /* must be sorted */
                   SEXP maxindex,
                   SEXP many
            ){

            PROTECT(index = AS_NUMERIC(index));
            PROTECT(maxindex = AS_NUMERIC(maxindex));
            PROTECT(many = AS_INTEGER(many));
            double *ind;
            ind = NUMERIC_POINTER(index);
            int maxind = INTEGER_VALUE(maxindex);
            int many2 = INTEGER_VALUE(many);
			int leng = LENGTH(index);
            int n,from,to=0,countit,atind;
            SEXP outit;
            PROTECT(outit = Rf_allocMatrix(INTSXP, maxind, many2));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n=0;n<(maxind*many2);n++){
                *(at+n) = 0;
            }

            from=1;
            to=1;
            countit=1;
            atind=0;
            for(n=1;n<leng;n++){
                if(*(ind+n)!=*(ind+n-1)){
                	if( (*(ind+n-1)!=0) ){
                   		*(at+atind)=from;
                        *(at+((maxind*1)+atind))=to;
                   		*(at+((maxind*2)+atind))=countit;
                   		atind++;
                    };
                	from=(n+1);
                    countit=1;
                    to=from;
                }else{
                    to++;
                    countit++;
                }
            }

            n--;
            if((*(ind+n)!=0)  ){
                *(at+atind)=from;
                *(at+((maxind*1)+atind))=to;
                *(at+((maxind*2)+atind))=countit;
            }

            SETLENGTH(outit,maxind*many2);
            UNPROTECT(4);
            return outit;

}

/******************************************************************************/
/* assemble indices for return value of Cagglom *******************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP fill_timeset(      SEXP timeset,
                        SEXP sampleID,
                        SEXP intensity,
                        SEXP lengtimeset
                      ){

           PROTECT(timeset = AS_NUMERIC(timeset));
           PROTECT(sampleID = AS_NUMERIC(sampleID));
           PROTECT(intensity = AS_NUMERIC(intensity));
           PROTECT(lengtimeset = AS_INTEGER(lengtimeset));
           double *tset, *id, *inte;
           tset = NUMERIC_POINTER(timeset);
           id = NUMERIC_POINTER(sampleID);
           inte = NUMERIC_POINTER(intensity);
           int leng1 = LENGTH(intensity);
           int leng2 = INTEGER_VALUE(lengtimeset);
           int m,n;
           SEXP vec;
           PROTECT(vec = Rf_allocMatrix(REALSXP, leng2, 2));
           double *at;
           at = REAL(vec);
           for(m=0;m<2;m++){
               for(n=0;n<leng2;n++){
                   at[(m*leng2)+n]=0;
               }
           }

           for(m=0;m<leng1;m++){
               for(n=0;n<leng2;n++){
                   if(tset[(1*leng2)+n]!=0){
                       if(tset[(1*leng2)+n]==*(id+m)){
                           at[n]=*(inte+m);
                           continue;
                       }
                   }
                   if(tset[(2*leng2)+n]!=0){
                       if(tset[(2*leng2)+n]==*(id+m)){
                           at[(leng2)+n]=*(inte+m);
                       }
                   }
               }
           }

           UNPROTECT(5);
           return vec;

}

/******************************************************************************/
/* get maximum difference in an time ordered series of intensities ************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP meandel(      SEXP timeset,
                   SEXP subit,
                   SEXP subrat,
                   SEXP numtime,
                   SEXP getwhat,
                   SEXP lags,
                   SEXP threshold,
                   SEXP notrend
           ){

           PROTECT(timeset = AS_NUMERIC(timeset));
           PROTECT(subit = AS_INTEGER(subit));
           PROTECT(subrat = AS_NUMERIC(subrat));
           PROTECT(numtime = AS_NUMERIC(numtime));
           PROTECT(getwhat = AS_NUMERIC(getwhat));
           PROTECT(lags = AS_NUMERIC(lags));
           PROTECT(threshold = AS_NUMERIC(threshold));
           PROTECT(notrend = AS_INTEGER(notrend));
           double *inte,*lagit, *tim;
           inte = NUMERIC_POINTER(timeset);
           lagit = NUMERIC_POINTER(lags);
           tim = NUMERIC_POINTER(numtime);
           int leng = LENGTH(numtime);
           int lagnumb = LENGTH(lags);
           int subtractit = INTEGER_VALUE(subit);
           double subratio = NUMERIC_VALUE(subrat);
           double thres = NUMERIC_VALUE(threshold);
           int getit = INTEGER_VALUE(getwhat);
           int notrend2 = INTEGER_VALUE(notrend);
           double intdifto,blindint,atint,maxint,intsum,intcount,meanint,meanint_all,varint,varint_all,nowint;
           int m,n,k,maxat=0,minat,nowat,from,to=0,tosam,doit;
           SEXP scored;
           PROTECT(scored = Rf_allocMatrix(REALSXP,6,lagnumb));
           double *at;
           at = REAL(scored);

           /* run blind subtraction *******************************************/
           if(subtractit==1){
               /* any blind available? */
               minat=-1;
               for(n=0;n<leng;n++){
                   if(*(inte+(leng*2)+n)!=0){
                       if(minat==-1){
                           minat=n;
                       }
                       maxat=n;
                   }
               }
               /* interpolate blind values ***********************************/
               if(minat!=-1){
                   for(n=0;n<minat;n++){
                       *(inte+(leng*4)+n)=*(inte+(leng*4)+minat);
                   }
                   for(n=maxat;n<leng;n++){
                       *(inte+(leng*4)+n)=*(inte+(leng*4)+maxat);
                   }
                   if(minat<(maxat-1)){
                       from=minat;
                       while(from<maxat){
                           for(n=(from+1);n<leng;n++){
                               if(*(inte+(leng*2)+n)!=0){
                                   to=n;
                                   break;
                               }
                           }
                           if(from<(to-1)){
                               for(n=(from+1);n<=(to-1);n++){
                                  *(inte+(leng*4)+n) = (  *(inte+(leng*4)+from) + ( (*(inte+(leng*4)+to)-*(inte+(leng*4)+from)) * fabs((*(tim+n)-*(tim+from))/(*(tim+to)-*(tim+from))) )  );
                               }
                           }
                           from=to;
                       }
                   }
               }
               /* filter: above threshold subrat? ******************************/
               doit=0;
               for(n=0;n<(leng);n++){
                   if(*(inte+leng+n)!=0){
                       if(*(inte+(leng*3)+n)>(*(inte+(leng*4)+n)*subratio)){
                           *(inte+n)=2;
                           doit=1;
                       }else{
                           *(inte+n)=1;
                       }
                   }
               }
           }else{
               doit=1;
               for(n=0;n<(leng);n++){
                   *(inte+n)=2;
               }
           }

           /* detect peaks & store maximum value over all lags ****************/
           if(doit==1){
                for(m=0;m<lagnumb;m++){ // for each lag
                    doit=0;
                    k=0;
                    to=-1;
                    intsum=0;
                    intcount=0;
                    maxint=0;
                    maxat=0;
                    blindint=0;
                    /* initialize: get first nonblind = starting point */
                    while( (to+1)<leng ){
                        to++;
                        if(*(inte+(leng)+to)!=0){
                            intsum=(intsum+*(inte+(leng*3)+to));
                            intcount++;
                            tosam=to;
                            maxint=*(inte+(leng*3)+tosam);
                            blindint=*(inte+(leng*4)+tosam);
                            maxat=tosam;
                            *(inte+(leng*5)+(leng*m)+tosam)=*(inte+(leng*3)+tosam);
                            break;
                        }
                    }
                    from=to;
                    /* iterate over all non-blind time points */
                    while((to+1)<leng){
                        to++;
                        if(*(inte+(leng)+to)!=0){
                            /* deal with to */
                            intsum=(intsum+*(inte+(leng*3)+to));
                            intcount++;
                            tosam=to;
                            /* deal with from */
                            while(  (fabs(*(tim+to)-*(tim+from))>*(lagit+m)) && ((from+1)<leng)  ){
                                if(*(inte+(leng*3)+from)!=0){
                                    intsum=(intsum-*(inte+(leng*3)+from));
                                    intcount--;
                                }
                                from++;
                            }
                            /* check result */
                            if(intcount>0){
                                intdifto=(intsum/intcount);
                                if(*(inte+(leng*3)+tosam)>intdifto){
                                    *(inte+(leng*5)+(leng*m)+tosam)=intdifto;
                                    if(*(inte+(leng*3)+tosam)>maxint){
                                        maxint=*(inte+(leng*3)+tosam);
                                        blindint=*(inte+(leng*4)+tosam);
                                        maxat=tosam;
                                    }
                                    doit=1;
                                }else{
                                    if(doit==1){
                                        *(inte+(leng*5)+(leng*lagnumb)+(leng*m)+k)=maxint;
                                        *(inte+(leng*5)+(leng*lagnumb*2)+(leng*m)+k)=*(tim+maxat);
                                        *(inte+(leng*5)+(leng*lagnumb*3)+(leng*m)+k)=blindint;
                                        k++;
                                        doit=0;
                                    }
                                    *(inte+(leng*5)+(leng*m)+tosam)=*(inte+(leng*3)+tosam);
                                    intsum=*(inte+(leng*3)+tosam);
                                    intcount=1;
                                    maxint=*(inte+(leng*3)+tosam);
                                    blindint=*(inte+(leng*4)+tosam);
                                    maxat=tosam;
                                    from=tosam;
                               }
                            }
                        }
                    }

                    /* get absolute deviation & mean **************************/
                    /* have maximum peak excluded for meanint & varint ********/
                    /* subtract blind from maximum ****************************/
                    intsum=0;
                    for(to=0;to<k;to++){
                        intsum=(intsum+*(inte+(leng*5)+(leng*lagnumb)+(leng*m)+to));
                    }
                    maxint=0;
                    meanint_all=0;
                    varint_all=0;
                    maxat=-1;
                    if(k>1){ // more than one trend per lag detected?
                        for(to=0;to<k;to++){
                            atint=*(inte+(leng*5)+(leng*lagnumb)+(leng*m)+to);
                            blindint=*(inte+(leng*5)+(leng*lagnumb*3)+(leng*m)+to);
                            meanint=((intsum-*(inte+(leng*5)+(leng*lagnumb)+(leng*m)+to))/(k-1));
                            varint=0;
                            for(from=0;from<k;from++){
                                if(from!=to){
                                    varint=(varint+fabs(*(inte+(leng*5)+(leng*lagnumb)+(leng*m)+from)-meanint));
                                }
                            }
                            varint=(varint/(k-1));
                            if( ((atint-blindint)>maxint) && (atint>(blindint*subratio)) && (atint>(meanint+(thres*varint))) ){
                                maxint=(atint-blindint);
                                meanint_all=meanint;
                                varint_all=varint;
                                maxat=to;
                            }
                        }
                    }else{
                        maxint=intsum;
                    }
                    if(maxat>-1){
                        *(at+(m*6)+1)=*(inte+(leng*5)+(leng*lagnumb)+(leng*m)+maxat);
                    }else{
                        *(at+(m*6)+1)=0;
                    }
                    *(at+(m*6)+2)=varint_all;
                    *(at+(m*6)+4)=maxint;
                    *(at+(m*6)+5)=meanint_all;
                    /* instead of global trend, report maximum sample-blind intensity */
                    if(notrend2==1){
                        maxint=0;
                        for(n=0;n<(leng);n++){
                            if((*(inte+(leng*3)+n)-*(inte+(leng*4)+n))>maxint){
                                maxint=(*(inte+(leng*3)+n)-*(inte+(leng*4)+n));
                            }
                        }
                        *(at+(m*6)+4)=maxint; // overwrite above value
                    }
                    /* get absolute deviation & mean **************************/
                    /* have current peak excluded for meanint & varint ********/
                    /* subtract blind from newest *****************************/
                    /* get latest=newest value */
                    nowint=0;
                    nowat=0;
                    for(to=(leng-1);to>=0;to--){
                        if(*(inte+leng+to)!=0){
                            if(*(inte+to)==2){
                                nowint=*(inte+(leng*3)+to);
                                nowat=to;
                                break;
                            }else{
                                nowint=0;
                                nowat=to;
                                break;
                            }
                        }
                    }
                    if(nowint>0){
                        meanint=(intsum/k);
                        varint=0;
                        for(to=0;to<k;to++){
                            varint=(varint+fabs(*(inte+(leng*5)+(leng*lagnumb)+(leng*m)+to)-meanint));
                        }
                        if(k>0){
                            varint=(varint/k);
                        }
                        if(varint!=0){
                            *(at+(m*6))=nowint;
                            if( (nowint>(meanint+(thres*varint))) ){
                                *(at+(m*6)+3)=(nowint-*(inte+(leng*4)+nowat));
                            }else{
                                *(at+(m*6)+3)=0;
                            }
                        }else{
                            *(at+(m*6))=nowint;
                            if(nowint>0){
                                *(at+(m*6)+3)=(nowint-*(inte+(leng*4)+nowat));
                            }else{
                                *(at+(m*6)+3)=0;
                            }
                        }
                    }else{
                        *(at+(m*6))=0;
                        *(at+(m*6)+3)=0;
                    }
                }
           }else{
                for(m=0;m<lagnumb;m++){
                    *(at+(m*6))=0;      /* nowint */
                    *(at+(m*6)+1)=0;    /* maxint */
                    *(at+(m*6)+2)=0;    /* varint */
                    *(at+(m*6)+3)=0;    /* nowint-blind*/
                    *(at+(m*6)+4)=0;    /* maxint-blind */
                    *(at+(m*6)+5)=0;    /* meanint */
                }
           }

           if(getit==1){
               UNPROTECT(9);
               return scored;
           }else{
               UNPROTECT(9);
               return timeset;
           }

}

/******************************************************************************/
/* get maximum difference in an orderes series of intensities *****************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP intdiff(      SEXP timeset,
                   SEXP subit,
                   SEXP subrat,
                   SEXP numtime,
                   SEXP getwhat
           ){

           PROTECT(timeset = AS_NUMERIC(timeset));
           PROTECT(subit = AS_INTEGER(subit));
           PROTECT(subrat = AS_NUMERIC(subrat));
           PROTECT(numtime = AS_NUMERIC(numtime));
           PROTECT(getwhat = AS_NUMERIC(getwhat));
           double *inte;
           inte = NUMERIC_POINTER(timeset);
           int leng = LENGTH(timeset)/5;
           int subtractit = INTEGER_VALUE(subit);
           double subratio = NUMERIC_VALUE(subrat);
           int getit = INTEGER_VALUE(getwhat);
           double *tim;
           tim = NUMERIC_POINTER(numtime);
           double intdiffor=0, intdifback=0, minint,maxint;
           int m,n,maxat=0,minat,from,to=0,doit;
           SEXP intensidif;
           PROTECT(intensidif = Rf_allocMatrix(REALSXP, 1, 5));
           double *at;
           at = REAL(intensidif);

           /* run blind subtraction *******************************************/
           if(subtractit==1){
               /* any blind available? */
               minat=-1;
               for(n=0;n<leng;n++){
                   if(*(inte+(leng*4)+n)!=0){
                       if(minat==-1){
                           minat=n;
                       }
                       maxat=n;
                   }
               }
               /* interpolate blind values ***********************************/
               if(minat!=-1){
                   for(n=0;n<minat;n++){
                       *(inte+(leng*4)+n)=*(inte+(leng*4)+minat);
                   }
                   for(n=maxat;n<leng;n++){
                       *(inte+(leng*4)+n)=*(inte+(leng*4)+maxat);
                   }
                   if(minat<(maxat-1)){
                       from=minat;
                       while(from<maxat){
                           for(n=(from+1);n<leng;n++){
                               if(*(inte+(leng*4)+n)!=0){
                                   to=n;
                                   break;
                               }
                           }
                           if(from<(to-1)){
                               for(n=(from+1);n<=(to-1);n++){
                                  *(inte+(leng*4)+n) = (  *(inte+(leng*4)+from) + ( (*(inte+(leng*4)+to)-*(inte+(leng*4)+from)) * fabs((*(tim+n)-*(tim+from))/(*(tim+to)-*(tim+from))) )  );
                               }
                           }
                           from=to;
                       }
                   }
               }
               /* filter: above threshold subrat? subtract ********************/
               doit=0;
               for(n=0;n<(leng);n++){
                   if(*(inte+leng+n)!=0){
                       if(*(inte+(leng*3)+n)>(*(inte+(leng*4)+n)*subratio)){
                           *(inte+n)=2;
                           doit=1;
                       }else{
                           *(inte+n)=1;
                           if( *(inte+(leng*3)+n)<*(inte+(leng*4)+n) ){
                               *(inte+(leng*3)+n)=*(inte+(leng*4)+n);
                           }
                       }
                   }
               }
           }else{
               doit=1;
               for(n=0;n<(leng);n++){
                   *(inte+n)=2;
               }
           }

           /* forward globally ************************************************/
           if(doit==1){
               maxat=0;
               minat=0;
               maxint=0;
               intdiffor=0;
               for(n=0;n<(leng-1);n++){
                       minint=*(inte+(leng*3)+n);
                       if(maxat<=n){
                           for(m=(n+1);m<leng;m++){
                               if(*(inte+(leng*3)+m)>maxint){
                                   if(*(inte+m)==2){
                                       maxat=m;
                                       maxint=*(inte+(leng*3)+m);
                                   }
                               }
                           }
                       }
                       if((maxint-minint)>intdiffor){
                           minat=n;
                           intdiffor=(maxint-minint);
                       }
               }
               /* backward latest *************************************************/
               intdifback=0;
               if( (*(inte+(leng-1))>0) & (*(inte+(leng-1))==2)){
                   minint=*(inte+(leng-1));
                   for(m=(leng-2);m>=0;m--){
                       if(*(inte+m)>*(inte+(leng-1))){
                           break;
                       }else{
                           if(*(inte+(leng*3)+m)<minint){
                               minint=*(inte+(leng*3)+m);
                           }
                       }
                   }
                   intdifback=(*(inte+(leng*3)+(leng-1))-minint);
               }
           }

           if(getit==1){
               *(at)=intdiffor;
               *(at+1)=intdifback;
               *(at+2)=*(inte+maxat);
               *(at+3)=minat;
               *(at+4)=maxat;
               UNPROTECT(6);
               return intensidif;
           }else{
                UNPROTECT(6);
               return timeset;
           }

}

/******************************************************************************/
/* assemble data for plotting *************************************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP Cplot_prof(       SEXP RTlim_low,
                        SEXP RTlim_up,
                        SEXP mzlim_low,
                        SEXP mzlim_up,
                        SEXP mz,
                        SEXP RT,
                        SEXP intensity,
                        SEXP sampleID,
                        SEXP color1,
                        SEXP color2,
                        SEXP whatcolor
                        ){

           PROTECT(RTlim_low = AS_NUMERIC(RTlim_low));
           PROTECT(RTlim_up = AS_NUMERIC(RTlim_up));
           PROTECT(mzlim_low = AS_NUMERIC(mzlim_low));
           PROTECT(mzlim_up = AS_NUMERIC(mzlim_up));
           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(intensity = AS_NUMERIC(intensity));
           PROTECT(sampleID = AS_NUMERIC(sampleID));
           PROTECT(color1 = AS_NUMERIC(color1));
           PROTECT(color2 = AS_NUMERIC(color2));
           PROTECT(whatcolor = AS_INTEGER(whatcolor));
           double RTlim_low2 = NUMERIC_VALUE(RTlim_low);
           double RTlim_up2 = NUMERIC_VALUE(RTlim_up);
           double mzlim_low2 = NUMERIC_VALUE(mzlim_low);
           double mzlim_up2 = NUMERIC_VALUE(mzlim_up);
           int whatcol = INTEGER_VALUE(whatcolor);
           double *ret, *mass, *inte, *sam, maxint=0;
           ret = NUMERIC_POINTER(RT);
           mass = NUMERIC_POINTER(mz);
           inte = NUMERIC_POINTER(intensity);
           sam = NUMERIC_POINTER(sampleID);
           double *colorit1,*colorit2;
           colorit1 = NUMERIC_POINTER(color1);
           colorit2 = NUMERIC_POINTER(color2);
           int leng2 = LENGTH(RT);
           int n,counter1=0,counter2=0;
           int *does;
           does = new int[leng2];


           for(n=0;n<leng2;n++){
               if((*(mass+n)>=mzlim_low2) & (*(mass+n)<=mzlim_up2)){
                   if((*(ret+n)>=RTlim_low2) & (*(ret+n)<=RTlim_up2)){
                       does[n]=1;
                       if(*(inte+n)>maxint){
                           maxint=*(inte+n);
                       }
                       counter1++;
                   }else{does[n]=0;};
               }else{does[n]=0;};
           }

           if(counter1>0){
               SEXP ans;
               PROTECT(ans = Rf_allocMatrix(REALSXP, counter1, 6));
               double *rans;
               rans = REAL(ans);
               if(whatcol==1){
                   for(n=0;n<leng2;n++){
                       if(does[n]==1){
                           rans[counter2]=*(mass+n);
                           rans[counter2+counter1]=*(inte+n);
                           rans[counter2+(2*counter1)]=*(ret+n);
                           rans[counter2+(3*counter1)]=*(sam+n);
                           rans[counter2+(4*counter1)]=1;
                           rans[counter2+(5*counter1)]=0;
                           counter2++;
                       }
                   }
               }
               if(whatcol==2){
                   for(n=0;n<leng2;n++){
                       if(does[n]==1){
                           rans[counter2]=*(mass+n);
                           rans[counter2+counter1]=*(inte+n);
                           rans[counter2+(2*counter1)]=*(ret+n);
                           rans[counter2+(3*counter1)]=*(sam+n);
                           rans[counter2+(4*counter1)]=*(colorit1+n)+1;
                           rans[counter2+(5*counter1)]=0;
                           counter2++;
                       }
                   }
               }
               if(whatcol==3){
                   for(n=0;n<leng2;n++){
                       if(does[n]==1){
                           rans[counter2]=*(mass+n);
                           rans[counter2+counter1]=*(inte+n);
                           rans[counter2+(2*counter1)]=*(ret+n);
                           rans[counter2+(3*counter1)]=*(sam+n);
                           rans[counter2+(4*counter1)]=*(colorit2+n)+1;
                           rans[counter2+(5*counter1)]=0;
                           counter2++;
                       }
                   }
               }
               rans[(5*counter1)]=maxint;
               rans[(5*counter1)+1]=counter1;
               SETLENGTH(ans, counter1*6);
               delete[] does;
			   UNPROTECT(12);
               return ans;
           }else{
               delete[] does;
			   UNPROTECT(11);
               return R_NilValue;
           }

}

/******************************************************************************/
/* RT-binning of data for plotting ********************************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP binRT_prof(   SEXP RT,
                   SEXP intensity,
                   SEXP binRT,
                   SEXP colorit,
                   SEXP what
                   ){

           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(intensity = AS_NUMERIC(intensity));
           PROTECT(binRT = AS_NUMERIC(binRT));
           PROTECT(colorit = AS_NUMERIC(colorit));
           PROTECT(what = AS_INTEGER(what));
           int leng2 = LENGTH(RT);
           int leng3 = LENGTH(binRT);
           int what2 = INTEGER_VALUE(what);
           double *ret, *intens, *binit, *col;
           ret = NUMERIC_POINTER(RT);
           intens = NUMERIC_POINTER(intensity);
           binit = NUMERIC_POINTER(binRT);
           col = NUMERIC_POINTER(colorit);
           int n,m,l,counter;
           SEXP ans;
           PROTECT(ans = Rf_allocMatrix(REALSXP, leng3-1, 2));
           double *rans;
           rans = REAL(ans);
           for(m=0;m<(leng3-1);m++){
               rans[m]=0;
               rans[m+leng3-1]=0;
           };


           if(what2==1){
		       l=0;
               counter=0;
               for(n=0;n<leng2;n++){
                   if(*(col+n)>0){
                       for(m=l;m<(leng3-1);m++){
                           if((*(ret+n)>=*(binit+m))&(*(ret+n)<*(binit+m+1))){
                               if(rans[m]==0){counter++;};
                                    rans[m]=(rans[m]+*(intens+n));
                                    rans[m+leng3-1]=((*(binit+m)+*(binit+m+1))/2);
                               l=m;
                               break;
                           }
                       }
                   }
               }
           }else{
		       l=0;
               counter=0;
               for(n=0;n<leng2;n++){
                   if(*(col+n)>0){
                       for(m=l;m<(leng3-1);m++){
                           if((*(ret+n)>=*(binit+m))&(*(ret+n)<*(binit+m+1))){
                               if(rans[m]==0){counter++;};
                                   if(rans[m]<*(intens+n)){
                                       rans[m]=*(intens+n);
                                       rans[m+leng3-1]=((*(binit+m)+*(binit+m+1))/2);
                                   }
                               l=m;
                               break;
                           }
                       }
                   }
               }
           }


           /* omit 0-entries in ans2, i.e., resize to within mzlimit set in R */
           SEXP ans2;
           PROTECT(ans2 = Rf_allocMatrix(REALSXP, counter, 2));
           double *rans2;
           rans2 = REAL(ans2);
           n=0;
           for(m=0;m<(leng3-1);m++){
               if(rans[m]!=0){
                   rans2[n]=rans[m];
                   rans2[counter+n]=rans[m+leng3-1];
                   n++;
               }
           }
           UNPROTECT(7);
           return ans2;
}

/******************************************************************************/
/* RT-binning of data for plotting ********************************************/
/******************************************************************************/
// [[Rcpp::export]]
SEXP binmz_prof(   SEXP mz,
                   SEXP intensity,
                   SEXP binmzs,
                   SEXP colorit
                   ){

           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(intensity = AS_NUMERIC(intensity));
           PROTECT(binmzs = AS_NUMERIC(binmzs));
           PROTECT(colorit = AS_NUMERIC(colorit));
           int leng2 = LENGTH(mz);
           int leng3 = LENGTH(binmzs);
           double *mass, *intens, *binit, *col;
           mass = NUMERIC_POINTER(mz);
           intens = NUMERIC_POINTER(intensity);
           binit = NUMERIC_POINTER(binmzs);
           col = NUMERIC_POINTER(colorit);
           int n,m,l,counter;
           SEXP ans;
           PROTECT(ans = Rf_allocMatrix(REALSXP, leng3-1, 2));
           double *rans;
           rans = REAL(ans);
           for(m=0;m<(leng3-1);m++){
               rans[m]=0;
               rans[m+leng3-1]=0;
           };

           l=0;
           counter=0;
           for(n=0;n<leng2;n++){
               if(*(col+n)>0){
                   for(m=l;m<(leng3-1);m++){
                       if((*(mass+n)>=*(binit+m))&(*(mass+n)<*(binit+m+1))){
                           if(rans[m]==0){counter++;};
                               if(rans[m]<*(intens+n)){
                                   rans[m]=*(intens+n);
                                   rans[m+leng3-1]=((*(binit+m)+*(binit+m+1))/2);
                               }
                           l=m;
                           break;
                       }
                   }
               }
           }

           /* omit 0-entries in ans2, i.e., resize to within mzlimit set in R */
           SEXP ans2;
           PROTECT(ans2 = Rf_allocMatrix(REALSXP, counter, 2));
           double *rans2;
           rans2 = REAL(ans2);
           n=0;
           for(m=0;m<(leng3-1);m++){
               if(rans[m]!=0){
                   rans2[n]=rans[m];
                   rans2[counter+n]=rans[m+leng3-1];
                   n++;
               }
           }
           UNPROTECT(6);
           return ans2;

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


