/*
This code implements the 3D B0 corrected CEST map, positive and negative images at ppm of interest.
function [pimgout,nimgout] = getB0correctedCEST3Dimgfit(ppm,ppmlist,posimages,negimages,refimg,B0map,mask,b0ppmstep)
*/

/* You can include any C libraries that you normally use */

#include "mex.h"   /*This C library is required*/
#include <math.h>

#define PI 3.141592653589793

/* Underlying Math routines */

#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}
static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)


void nrerror(char *error_text)
{
 	// mexPrintf("%s\n",error_text);
}

double *vector(int nl,int nh)
{
	double *v;

	v=(double *)calloc((unsigned) (nh-nl+1), sizeof(double));
	if (!v)
        nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(int nl,int nh)
{
	int *v;

	v=(int *)calloc((unsigned) (nh-nl+1), sizeof(int));
	if (!v)
        nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(int nl,int nh)
{
	double *v;

	v=(double *)calloc((unsigned) (nh-nl+1), sizeof(double));
	if (!v)
        nrerror("allocation failure in dvector()");
	return v-nl;
}

double **matrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m = (double **) calloc((unsigned) (nrh-nrl+1), sizeof(double*));
	if (!m)
        nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++)
    {
		m[i] = (double *) calloc((unsigned) (nch-ncl+1), sizeof(double));
		if (!m[i])
            nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) calloc((unsigned) (nrh-nrl+1), sizeof(double*));
	if (!m)
        nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++)
    {
		m[i]=(double *) calloc((unsigned) (nch-ncl+1), sizeof(double));
		if (!m[i])
            nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(int nrl,int nrh,int ncl,int nch)
{
	int i,**m;

	m=(int **)calloc((unsigned) (nrh-nrl+1), sizeof(int*));
	if (!m)
        nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++)
    {
		m[i]=(int *)calloc((unsigned) (nch-ncl+1), sizeof(int));
		if (!m[i])
            nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

double **submatrix(double **a, int oldrl, int oldrh,int oldcl,int oldch,int newrl,int newcl)
{
	int i,j;
	double **m;

	m=(double **) calloc((unsigned) (oldrh-oldrl+1), sizeof(double*));
	if (!m)
        nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++)
        m[j]=a[i]+oldcl-newcl;

	return m;
}

void free_vector(double *v,int nl,int nh)
{
	free(v+nl);
}

void free_ivector(int *v,int nl,int nh)
{
	free(v+nl);
}

void free_dvector(double *v,int nl,int nh)
{
	free(v+nl);
}

void free_matrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--)
        free((m[i]+ncl));
	free( (m+nrl));
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--)
        free((m[i]+ncl));
	free( (m+nrl));
}


void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--)
        free(m[i]+ncl);
	free(m+nrl);
}

void free_submatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--)
        free((m[i]+ncl));
	free( (m+nrl));
}

double **convert_matrix(double *a, int nrl, int nrh, int ncl,int nch)
{
	int i,j,nrow,ncol;
	double **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (double **) calloc((unsigned) (nrow), sizeof(double*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++)
        m[j]=a+ncol*i-ncl;
	return m;
}

void free_convert_matrix(double **b, int nrl, int nrh, int ncl, int nch)
{
	free(b+nrl);
}

void gaussj(double **a, int n,double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++)
        ipiv[j]=0;
	for (i=1;i<=n;i++)
    {
		big=0.0;
		for (j=1;j<=n;j++)
        {
			if (ipiv[j] != 1)
            {
				for (k=1;k<=n;k++)
                {
					if (ipiv[k] == 0)
                    {
						if (fabs(a[j][k]) >= big)
                        {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
                    else if (ipiv[k] > 1)
                        nrerror("GAUSSJ: Singular Matrix-1");
				}
            }
        }
		++(ipiv[icol]);
		if (irow != icol)
        {
			for (l=1;l<=n;l++)
                SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++)
                SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0)
            nrerror("GAUSSJ: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++)
            a[icol][l] *= pivinv;
		for (l=1;l<=m;l++)
            b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
        {
			if (ll != icol)
            {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
        }
	}
	for (l=n;l>=1;l--)
    {
		if (indxr[l] != indxc[l])
        {
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

void covsrt(double **covar, int ma,int *lista, int mfit)
{
	int i,j;
	double swap;

	for (j=1;j<ma;j++)
    {
		for (i=j+1;i<=ma;i++)
            covar[i][j]=0.0;
    }
	for (i=1;i<mfit;i++)
    {
		for (j=i+1;j<=mfit;j++)
        {
			if (lista[j] > lista[i])
				covar[lista[j]][lista[i]]=covar[i][j];
			else
				covar[lista[i]][lista[j]]=covar[i][j];
		}
    }
	swap=covar[1][1];
	for (j=1;j<=ma;j++)
    {
		covar[1][j]=covar[j][j];
		covar[j][j]=0.0;
	}
	covar[lista[1]][lista[1]]=swap;
	for (j=2;j<=mfit;j++)
        covar[lista[j]][lista[j]]=covar[1][j];
	for (j=2;j<=ma;j++)
    {
		for (i=1;i<=j-1;i++)
            covar[i][j]=covar[j][i];
    }
}

void polyfit(double *x, double *y, double *sig, int ndata, double *a, int ma, int *lista, int mfit,
          double **covar,double *chisq, void (*funcs)(double,double *,int) )
{
	int k,kk,j,ihit,i;
	double ym,wt,sum,sig2i,**beta,*afunc;

	beta=matrix(1,ma,1,1);
	afunc=vector(1,ma);
	kk=mfit+1;
	for (j=1;j<=ma;j++)
    {
		ihit=0;
		for (k=1;k<=mfit;k++)
        {
			if (lista[k] == j)
                ihit++;
        }
		if (ihit == 0)
			lista[kk++]=j;
		else if (ihit > 1)
            nrerror("Bad LISTA permutation in LFIT-1");
	}
	if (kk != (ma+1))
        nrerror("Bad LISTA permutation in LFIT-2");
	for (j=1;j<=mfit;j++)
    {
		for (k=1;k<=mfit;k++)
            covar[j][k]=0.0;
		beta[j][1]=0.0;
	}
	for (i=1;i<=ndata;i++)
    {
		(*funcs)(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma)
        {
			for (j=(mfit+1);j<=ma;j++)
				ym -= a[lista[j]]*afunc[lista[j]];
        }
		sig2i=1.0/SQR(sig[i]);
		for (j=1;j<=mfit;j++)
        {
			wt=afunc[lista[j]]*sig2i;
			for (k=1;k<=j;k++)
				covar[j][k] += wt*afunc[lista[k]];
			beta[j][1] += ym*wt;
		}
	}
	if (mfit > 1)
    {
		for (j=2;j<=mfit;j++)
        {
			for (k=1;k<=j-1;k++)
				covar[k][j]=covar[j][k];
        }
    }
	gaussj(covar,mfit,beta,1);
	for (j=1;j<=mfit;j++)
        a[lista[j]]=beta[j][1];
	*chisq=0.0;
	for (i=1;i<=ndata;i++)
    {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++)
            sum += a[j]*afunc[j];
		*chisq += SQR((y[i]-sum)/sig[i]);
	}
	covsrt(covar,ma,lista,mfit);
	free_vector(afunc,1,ma);
	free_matrix(beta,1,ma,1,1);
}

#undef SWAP
#undef SQR


void xpoly(double x,double *afunc,int mma)
{
	int i;

	afunc[1]=1.0;
	for (i = 2; i <= mma;i++)
        afunc[i]=x*afunc[i-1];
}


void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;

    u = calloc( n, sizeof(double) );
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else
    {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++)
    {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
        qn=un=0.0;
	else
    {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free(u);
}
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
	int klo,khi,k;
	double h,b,a;

	klo=0;
	khi=n-1;
	while (khi-klo > 1)
    {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double * calc_b0corrected_cest(double ppmval, double *ppmlist, double *pimg,
                               double *nimg, double *B0map, long *mask, int nppm,int nelem)


   // nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Inside mexFunction---*/
    /*Declarations*/

    double  *ppmlistasc; //, *pimgout, *nimgout; /* Pointers to ppmlist,W,H,pos and negative images, and mask arrays */
    double *out;
		// double ppmval, *ppmlist, *pimg, *nimg, *refimg, *B0map, *mask
    double maxppm, minppm, stepppm, B0val, temp;  // ,refval;
    double zsp2p, zsp2n, *xp, *xn, *yp, *yn, ppmpos, ppmneg;
    // mwSize ndim;
    // mwSize nintppm, ndim, ndimimg;
    // mwSize nppm, nelem, nelemimg, ielem, ippm1, ippm2;
    int ielem, ippm1, ippm2, intfactor, i, j, k, NDEG, nintppm;
    const mwSize *dim_array_mask, *dim_array_image;
    // mwSize *psort;
    int *psort;

    // mwSize intfactor, i, j, k, NDEG;
    double *xint,*ypint,*ynint, *y2;


    int i1, i2;
    int *listap, *listan, NTERM, NPTS;
    double *xzsp, *afunc;
	double *yp1, chisqp,*ap,*sigmap,**covarp;
	double *yn1, chisqn,*an,*sigman,**covarn;
    double temp1,temp2;

    void (*fp) (double, double *, int);
    fp = xpoly;

    NDEG = 2;
    NTERM = NDEG+1;

    // ppmval = mxGetScalar(prhs[0]);
    // mexPrintf("ppmval = %f\n", ppmval);

    // ppmlist = mxGetPr(prhs[1]);
    // nppm=mxGetNumberOfElements(prhs[1]);

    intfactor = 100;
    nintppm = 1 + (nppm-1)*intfactor;
    NPTS = nintppm;

    // for (ippm1 = 0; ippm1 < nppm; ippm1++)
        // mexPrintf(" ppmlist[%i] = %f \n", ippm1, ppmlist[ippm1]);

    // pimg = mxGetPr(prhs[2]);
    // nimg = mxGetPr(prhs[3]);
    // nelemimg = mxGetNumberOfElements(prhs[2]);
    // ndimimg=mxGetNumberOfDimensions(prhs[2]);

    // refimg = mxGetPr(prhs[4]);
    // B0map = mxGetPr(prhs[5]);
    // mask = mxGetPr(prhs[6]);

    /* Get the dimensions in the mask array */
    // nelem = mxGetNumberOfElements(prhs[6]);
    // ndim=mxGetNumberOfDimensions(prhs[6]);
    // TODO
    // dim_array_mask=mxGetDimensions(prhs[6]);

    // TODO
    // plhs[0] = mxCreateNumericArray(ndim,dim_array_mask,mxDOUBLE_CLASS,mxREAL);
    // pimgout = mxGetPr(plhs[0]);
    // plhs[1] = mxCreateNumericArray(ndim,dim_array_mask,mxDOUBLE_CLASS,mxREAL);
    // nimgout = mxGetPr(plhs[1]);
		out = (double *)malloc(nelem * 2 * sizeof(double));
		// pimgout = (double *)malloc(nelem * sizeof(double));
		// nimgout = (double *)malloc(nelem * sizeof(double));


    // out = (double **)malloc(2 * sizeof(double *))


    xp = (double *)malloc(nppm * sizeof(double)); /* positive ppm values for fit */
    xn = (double *)malloc(nppm * sizeof(double)); /* negative ppm values for fit */
    yp = (double *)malloc(nppm * sizeof(double)); /* positive yvalues for fit */
    yn = (double *)malloc(nppm * sizeof(double)); /* negative yvalues for fit */
    y2 = (double *)malloc( nppm * sizeof(double) );  /* storage for second derivative */

    xint  = (double *)malloc(nintppm * sizeof(double)); /* Scaled to mean and std of xin */
    ypint = (double *)malloc(nintppm * sizeof(double)); /* Scaled to mean and std of xin */
    ynint = (double *)malloc(nintppm * sizeof(double)); /* Scaled to mean and std of xin */

    ppmlistasc = malloc(nppm *  sizeof(double));
    psort      = malloc(nppm *  sizeof(int));
    maxppm = -10000.0;
    minppm = 10000.0;
    for (ippm1 =0; ippm1 < nppm; ippm1++)
    {
        if (ppmlist[ippm1] > maxppm)
            maxppm = ppmlist[ippm1];
        if (ppmlist[ippm1] < minppm)
            minppm = ppmlist[ippm1];
    }
    stepppm = (maxppm - minppm)/(nppm-1);
    for (ippm1 =0; ippm1 < nppm; ippm1++)
    {
        ppmlistasc[ippm1] = minppm + ippm1*stepppm;
    }
    for (ippm1 = 0; ippm1 < nppm; ippm1++)
    {
        psort[ippm1] = ippm1;
        for (ippm2 = 0; ippm2 < nppm; ippm2++)
        {
            if ( fabs(ppmlistasc[ippm2] - ppmlist[ippm1]) < 0.001 )
                psort[ippm1] = ippm2;
        }
    }

    for (ippm1 = 0; ippm1 < nintppm; ippm1++)
    {
        xint[ippm1] = minppm + ippm1 * (stepppm/intfactor);
    }

    afunc  = vector(1,NTERM);
	xzsp   = vector( 1, NPTS );
    for (i1 =0; i1 < NPTS; i1++)
    {
        xzsp[i1+1] = xint[i1];
    }

	listap = ivector( 1,NTERM);
	ap     = vector ( 1, NTERM );
	covarp = matrix( 1, NTERM, 1, NTERM );

	listan = ivector( 1,NTERM);
    for (i2 = 1; i2 <= NTERM; i2++)
    {
        listap[i2] = i2;
        listan[i2] = i2;
    }

	an     = vector ( 1, NTERM );
	covarn = matrix( 1, NTERM, 1, NTERM );

	yp1     = vector( 1 , NPTS );
	sigmap = vector( 1, NPTS );
	yn1     = vector( 1 , NPTS );
	sigman = vector( 1, NPTS );

    for (ielem = 0; ielem < nelem; ielem++)
    {
        // pimgout[ielem] = 0.0;
				out[ielem] = 0.0;
				out[ielem + nelem] = 0.0;
        // nimgout[ielem] = 0.0;
				// out[1][ielem] = 0.0;
        B0val = B0map[ielem];
        // refval = refimg[ielem];

        if ( (mask[ielem] > 0) )
        {
            ppmpos = ppmval;
            ppmneg = ppmval;


            for  (ippm1 = 0; ippm1 < nppm; ippm1++)
            {
                i1 = psort[ippm1]*nelem + ielem;
                xp[ippm1] = ppmlistasc[ippm1]-B0val;
                xn[ippm1] = ppmlistasc[ippm1]+B0val;
//                 yp[ippm1] = 100.0*pimg[i1]/(refval+0.1);
//                 yn[ippm1] = 100.0*nimg[i1]/(refval+0.1);
                yp[ippm1] = pimg[i1];
                yn[ippm1] = nimg[i1];
            }

            if (nppm > 2)   /* Use cubic spline interpolation */
            {

                /* Call spline to get second derivatives and splint to interpolate Positive side */
                spline(xp, yp, nppm, yp[0], yp[nppm-1], y2);
                for (i2 =0; i2< nintppm; i2++)
                {
                    splint(xp, yp, y2, nppm, xint[i2], &ypint[i2]);
                }

                /* Call spline to get second derivatives and splint to interpolate Negative side */
                spline(xn, yn, nppm, yn[0], yn[nppm-1], y2);
                for (i2 =0; i2< nintppm; i2++)
                {
                    splint(xn, yn, y2, nppm, xint[i2], &ynint[i2]);
                }
            }
            else /* Use liner interpolation */
            {
                ypint[0] = yp[0];
                ypint[nintppm-1] = yp[nppm-1];
                ynint[0] = yn[0];
                ynint[nintppm-1] = yn[nppm-1];

                for (i2 = 1; i2 < nintppm-1 ; i2++)
                {
                    temp1 = xint[i2]-xint[0];
                    // i = (mwSize)temp1/stepppm;
                    i = (int)temp1/stepppm;
                    j = (fabs(temp1 - i*stepppm) < 1.0e-6) ? i : i+1;
                    if (i == j)
                    {
                        ypint[i2] = yp[i];
                        ynint[i2] = yn[i];
                    }
                    else
                    {
                        temp2 = (xint[i2] - xp[i])/stepppm;
                        ypint[i2] = yp[i] + ( yp[j] - yp[i] ) * temp2;
                        temp2 = (xint[i2] - xn[i])/stepppm;
                        ynint[i2] = yn[i] + ( yn[j] - yn[i] ) * temp2;
                    }
                }
            }

            for (i2 = 0; i2 < NPTS; i2++)
            {
                yp1[i2+1] = ypint[i2];
                yn1[i2+1] = ynint[i2];
                sigmap[i2+1] = 1.0;
                sigman[i2+1] = 1.0;
            }
            for (i2 = 1; i2 <= NTERM; i2++)
            {
                listap[i2] = i2;
                listan[i2] = i2;
            }
            /* Now call polyfit routine for positive side to get polynomial coeffs */
            polyfit(xp,yp,sigmap,nppm,ap,NTERM,listap,NTERM,covarp,&chisqp,fp);
            xpoly(ppmpos,afunc,NTERM);
            zsp2p = 0.0;
            for (i2=1;i2<=NTERM;i2++)
                zsp2p += (afunc[i2]*ap[i2]);

            /* Now call polyfit routine for negative side to get polynomial coeffs */
            polyfit(xn,yn,sigman,nppm,an,NTERM,listan,NTERM,covarn,&chisqn,fp);
            xpoly(ppmneg,afunc,NTERM);
            zsp2n = 0.0;
            for (i2=1;i2<=NTERM;i2++)
                zsp2n += (afunc[i2]*an[i2]);

            /* Load B0 corrected  CEST images and map */
//             pimgout[ielem] = zsp2p*refval/100.0;
//             nimgout[ielem] = zsp2n*refval/100.0;
            // pimgout[ielem] = zsp2p;
						out[ielem] = zsp2p;
						out[ielem + nelem] = zsp2n;
            // nimgout[ielem] = zsp2n;
        }
    }

    free(psort);
    free(ppmlistasc);
    free(xint);
    free(ypint);
    free(ynint);
    free(yp);
    free(yn);

		free_vector(xzsp,1,NPTS);
		free_vector(afunc,1,NTERM);

    free_matrix(covarp,1,NTERM,1,NTERM);
		free_vector(sigmap,1,NPTS);
		free_vector(yp1,1,NPTS);
		free_vector(ap,1,NTERM);
		free_ivector(listap,1,NTERM);

    free_matrix(covarn,1,NTERM,1,NTERM);
		free_vector(sigman,1,NPTS);
		free_vector(yn1,1,NPTS);
		free_vector(an,1,NTERM);
		free_ivector(listan,1,NTERM);

    return out;
}
