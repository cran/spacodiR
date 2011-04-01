#include <stdio.h>
#include <R.h> 
#include "Xatools.h"
#define NR_END 1
#define FREE_ARG char* 
#ifndef ERRORFILE
	#define ERRORFILE "error.txt"
#endif

/*
np: #plots
ns: #species
sp_plot: #species-plot abundances
abundtype: type of abundance (0=presence/absence, 1=relative frequency, 2=number of indiv)
distmat: dist matrix btw species
stat: output (diversity partition through space and/or time)
*/



void spacodi(int *np, int *ns, double *sp_plot, double *distmat, int *abundtype,
			 int *Ndclass, double *dclass, double *Ist_out, double *Pst_out, double *Bst_out, double *PIst_out,
			 double *pairwiseIst_out, double *pairwisePst_out, double *pairwiseBst_out, double *pairwisePIst_out) // Modified by jme 01-12-10
{ 
	int NP,NS,abundt,Ndivc;
	int s1,s2,p1,p2,c, **sps,counter; //Modified by jme 01-12-10
	double **fps, **dist, divc[102];
	double phylod,maxdist,F1,F2;
	float ***DivIijc,***DivPijc,***DivPIijc,***DivSijc;
	double **Ist,**Pst,**Bst,**PIst;
	double divIb[102],divIw[102],divPb[102],divPw[102],divBb[102],divBw[102],divPIb[102],divPIw[102],sumIb,sumIw1,sumIw2,sumPb,sumPw1,sumPw2,sumSb,sumSw1,sumSw2,sumPIb,sumPIw1,sumPIw2;
	double Istc[102],Pstc[102],Bstc[102],PIstc[102];

	NP=(*np);				//# plots
	NS=(*ns);				//# species
	abundt=(*abundtype);	//type of abundance (0, 1 or 2)
	Ndivc=(*Ndclass);		//# of classes of divergence intervals

	fps=dmatrix(0,NP,0,NS);	//frequency per plot and species
	dist=dmatrix(0,NS,0,NS);//divergence matrix between species	

	for(s1=1;s1<=NS;s1++)for(s2=1;s2<=NS;s2++) dist[s1][s2]=distmat[(s1-1)+NS*(s2-1)];
	for(p1=1;p1<=NP;p1++)for(s1=1;s1<=NS;s1++) fps[p1][s1]=sp_plot[(s1-1)+NS*(p1-1)];
	for(c=1;c<=Ndivc;c++) divc[c]=dclass[c-1];

	//if dist classes are given, check if last class is larger or equal to max dist, 
	// otherwise add a class
	maxdist=0.;
	for(s1=1;s1<NS;s1++)for(s2=s1+1;s2<=NS;s2++) if(maxdist<dist[s1][s2]) maxdist=dist[s1][s2];
	if(Ndivc){
		if(divc[Ndivc]<maxdist){
			Ndivc++;
			divc[Ndivc]=maxdist;
			//(*Ndclass)++;
			//dclass[Ndivc-1]=maxdist;
		}
	}
	else{
		Ndivc=1;
		divc[1]=maxdist;
	}

	//identify species present in each plot and attribute a new numerotation (to speed loops)
	// sps[plot][0]=number of species in plot
	// sps[plot][new sp number]=absolute sp number 
	// where 'new sp number' ranges from 1 to number of species in plot, 
	// and 'absolute sp number' ranges from 1 to NS
	sps=imatrix(0,NP,0,NS);
	for(p1=0;p1<=NP;p1++){
		sps[p1][0]=0;
		for(s1=1;s1<=NS;s1++) if(fps[p1][s1]){
			sps[p1][0]++;
			sps[p1][sps[p1][0]]=s1;
		}
	}

	//transform abundances in relative frequencies per plot
	//sum of abundances per plot stored in fps[plot][0]
	for(p1=1;p1<=NP;p1++){
		fps[p1][0]=0.;
		for(s1=1;s1<=sps[p1][0];s1++) fps[p1][0]+=fps[p1][sps[p1][s1]];
		for(s1=1;s1<=sps[p1][0];s1++) fps[p1][sps[p1][s1]]/=fps[p1][0];
	}


	//create 3-dim arrays to store values per pair of plots and per divergence classes
	// c=0 for pairs intra-sp, c=-1 for sums over all classes
	DivIijc=f3tensor(0,NP,0,NP,-1,Ndivc);	//freq of pairs of ind 
	DivPijc=f3tensor(0,NP,0,NP,-1,Ndivc);	//mean dist between ind
	DivSijc=f3tensor(0,NP,0,NP,-1,Ndivc);	//freq of pairs of species
	DivPIijc=f3tensor(0,NP,0,NP,-1,Ndivc);	//mean dist between species	


	//compute diversity within and between plots 
	for(p1=1;p1<=NP;p1++)for(p2=p1;p2<=NP;p2++) {
		for(c=-1;c<=Ndivc;c++) DivIijc[p1][p2][c]=DivPijc[p1][p2][c]=DivSijc[p1][p2][c]=DivPIijc[p1][p2][c]=0.f;

		for(s1=1;s1<=sps[p1][0];s1++){
			F1=fps[p1][sps[p1][s1]];
			for(s2=1;s2<=sps[p2][0];s2++){
				F2=fps[p2][sps[p2][s2]];
				if(p1==p2 && s1==s2 && abundt==2) F2=((F2*fps[p1][0])-1.)/(fps[p1][0]-1.); //sample size correction when abundnaces are individuals counts

				if(sps[p1][s1]==sps[p2][s2]){ 
					c=0;
					phylod=0.;
				}
				else{
					phylod=dist[sps[p1][s1]][sps[p2][s2]];	//phyletic distance between species (0 for a species with itself)
					c=1;
					while(divc[c]<phylod) c++;
				}

				DivIijc[p1][p2][c]+=(float)(F1*F2);				//prob identity		
				DivPijc[p1][p2][c]+=(float)(F1*F2*phylod);		//mean dist btw ind
				DivSijc[p1][p2][c]++;							//# pairs of sp
				DivPIijc[p1][p2][c]+=(float)(phylod);			//mean dist btw sp
			}  //end of loop s2
		}  //end loop s1
		
		//convert into mean dist btw ind or sp per class
		for(c=0;c<=Ndivc;c++) DivPijc[p1][p2][c]/=DivIijc[p1][p2][c];
		for(c=1;c<=Ndivc;c++) DivPIijc[p1][p2][c]/=DivSijc[p1][p2][c];

		//sums
		for(c=0;c<=Ndivc;c++) DivIijc[p1][p2][-1]+=DivIijc[p1][p2][c];
		for(c=1;c<=Ndivc;c++) DivSijc[p1][p2][-1]+=DivSijc[p1][p2][c];

		//proportions
		for(c=0;c<=Ndivc;c++) DivIijc[p1][p2][c]/=DivIijc[p1][p2][-1];
		for(c=1;c<=Ndivc;c++) DivSijc[p1][p2][c]/=DivSijc[p1][p2][-1];

		//mean dist btw ind or sp cumulated over classes
		for(c=1;c<=Ndivc;c++){
			DivPijc[p1][p2][-1]+=DivPijc[p1][p2][c]*DivIijc[p1][p2][c];
			DivPIijc[p1][p2][-1]+=DivPIijc[p1][p2][c]*DivSijc[p1][p2][c];
		}

	}
	
	Ist=dmatrix(0,NP,0,NP);
	Pst=dmatrix(0,NP,0,NP);
	Bst=dmatrix(0,NP,0,NP);
	PIst=dmatrix(0,NP,0,NP);

	for(c=0;c<=Ndivc;c++) divIb[c]=divIw[c]=divPb[c]=divPw[c]=divBb[c]=divBw[c]=divPIb[c]=divPIw[c]=0.;
	//pairwise Ist & cie
	counter=0; //Modified by jme 01-12-10
	for(p1=1;p1<=NP;p1++)for(p2=p1+1;p2<=NP;p2++) {
		pairwiseIst_out[counter]=Ist[p1][p2]=Ist[p2][p1]=1.- (((1.-DivIijc[p1][p1][0])+(1.-DivIijc[p2][p2][0]))/2.) / (1.-DivIijc[p1][p2][0]); //Modified by jme 01-12-10
		pairwisePst_out[counter]=Pst[p1][p2]=Pst[p2][p1]=1.- ((DivPijc[p1][p1][-1]+DivPijc[p2][p2][-1])/2.) / DivPijc[p1][p2][-1]; //Modified by jme 01-12-10
		pairwiseBst_out[counter]=Bst[p1][p2]=Bst[p2][p1]=1.- ((DivPijc[p1][p1][-1]/(1.-DivIijc[p1][p1][0])+DivPijc[p2][p2][-1]/(1.-DivIijc[p2][p2][0]))/2.) / (DivPijc[p1][p2][-1]/(1.-DivIijc[p1][p2][0])); //Modified by jme 01-12-10
		pairwisePIst_out[counter]=PIst[p1][p2]=PIst[p2][p1]=1.- ((DivPIijc[p1][p1][-1]+DivPIijc[p2][p2][-1])/2.) / DivPIijc[p1][p2][-1]; //Modified by jme 01-12-10
		counter=counter+1; //Modified by jme 01-12-10

		//sums for global stat
		divIb[0]+=(1.-DivIijc[p1][p2][0]);
		divIw[0]+=( (1.-DivIijc[p1][p1][0]) + (1.-DivIijc[p2][p2][0]) )/2.;
		divPb[0]+=DivPijc[p1][p2][-1];
		divPw[0]+=(DivPijc[p1][p1][-1]+DivPijc[p2][p2][-1])/2.;
		divBb[0]+=DivPijc[p1][p2][-1]/(1.-DivIijc[p1][p2][0]);
		divBw[0]+=( DivPijc[p1][p1][-1]/(1.-DivIijc[p1][p1][0]) + DivPijc[p2][p2][-1]/(1.-DivIijc[p2][p2][0]) )/2.;
		divPIb[0]+=DivPIijc[p1][p2][-1];
		divPIw[0]+=(DivPIijc[p1][p1][-1]+DivPIijc[p2][p2][-1])/2.;
	}
	Istc[0]=1.-divIw[0]/divIb[0];
	Pstc[0]=1.-divPw[0]/divPb[0];
	Bstc[0]=1.-divBw[0]/divBb[0];
	PIstc[0]=1.-divPIw[0]/divPIb[0];

	//global Ist by dist classes
	for(p1=1;p1<=NP;p1++)for(p2=p1+1;p2<=NP;p2++) {
		sumIb=sumIw1=sumIw2=sumPb=sumPw1=sumPw2=sumSb=sumSw1=sumSw2=sumPIb=sumPIw1=sumPIw2=0.;
		for(c=1;c<=Ndivc;c++){
			sumIb+=DivIijc[p1][p2][c];
			sumIw1+=DivIijc[p1][p1][c];
			sumIw2+=DivIijc[p2][p2][c];
			sumPb+=DivPijc[p1][p2][c]*DivIijc[p1][p2][c];
			sumPw1+=DivPijc[p1][p1][c]*DivIijc[p1][p1][c];
			sumPw2+=DivPijc[p2][p2][c]*DivIijc[p2][p2][c];
			divBb[c]+=sumPb/sumIb;
			divBw[c]+=((sumPw1/sumIw1)+(sumPw2/sumIw2))/2.;
			divPb[c]+=sumPb/(sumIb+DivIijc[p1][p2][0]);
			divPw[c]+=((sumPw1/(sumIw1+DivIijc[p1][p1][0]))+(sumPw2/(sumIw2+DivIijc[p2][p2][0])))/2.;

			sumSb+=DivSijc[p1][p2][c];
			sumSw1+=DivSijc[p1][p1][c];
			sumSw2+=DivSijc[p2][p2][c];
			sumPIb+=DivPIijc[p1][p2][c]*DivSijc[p1][p2][c];
			sumPIw1+=DivPIijc[p1][p1][c]*DivSijc[p1][p1][c];
			sumPIw2+=DivPIijc[p2][p2][c]*DivSijc[p2][p2][c];
			divPIb[c]+=sumPIb/sumSb;
			divPIw[c]+=((sumPIw1/sumSw1)+(sumPIw2/sumSw2))/2.;		
		}
	}
	for(c=1;c<=Ndivc;c++){
		Pstc[c]=1.-divPw[c]/divPb[c];
		Bstc[c]=1.-divBw[c]/divBb[c];
		PIstc[c]=1.-divPIw[c]/divPIb[c];
	}



	(*Ist_out)=Istc[0]; // Modified by cetp 17-08-09
	(*Pst_out)=Pstc[0]; // Modified by cetp 17-08-09
	(*Bst_out)=Bstc[0]; // Modified by cetp 17-08-09
	(*PIst_out)=PIstc[0];// Modified by cetp 17-08-09

	//free memory
	free_dmatrix(fps,0,NP,0,NS);
	free_dmatrix(dist,0,NS,0,NS);	
	free_imatrix(sps,0,NP,0,NS);
	free_f3tensor(DivIijc,0,NP,0,NP,-1,Ndivc);
	free_f3tensor(DivPijc,0,NP,0,NP,-1,Ndivc);
	free_f3tensor(DivSijc,0,NP,0,NP,-1,Ndivc);	
	free_f3tensor(DivPIijc,0,NP,0,NP,-1,Ndivc);	
	free_dmatrix(Ist,0,NP,0,NP);
	free_dmatrix(Pst,0,NP,0,NP);
	free_dmatrix(Bst,0,NP,0,NP);
	free_dmatrix(PIst,0,NP,0,NP);


} 




/********************************************************************/

double *dvector(long int nl, long int nh)
/*allocate a double vector with subscript range[nl..nh]*/
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in dvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_dvector(double *v, long int nl, long int nh)
/*free a double vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/****************************************************************************/

double **dmatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	/*allocate pointers to rows*/
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double *)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_dmatrix(double **m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/
int **imatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	/*allocate pointers to rows*/
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in imatrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in imatrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_imatrix(int **m, long nrl,long nrh,long ncl, long nch)
/*free an int matrix allocated by imatrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

float ***f3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a float 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float **)));
	if (!t) {
		write(ERRORFILE,"\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float *)));
	if (!t[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float )));
	if (!t[nrl][ncl]) {
		write(ERRORFILE,"\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_f3tensor(float ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an float f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	


/****************************************************************************/

