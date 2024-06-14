/*You can include any C libraries that you normally use*/
#include "math.h"
#include "mex.h"   /*This C library is required*/

//cc -fPIC -shared -o calc_b0_map.so calc_b0_map.c

/*
This code implements the maximum asymmetry WASSR B0map algorithm
function [freqWMap] = wassr_b0map(ppmlist, posimage, negimage, mask4all
 */

#define PI 3.141592653589793

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
	float h,b,a;

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

double symWASSR(double xval, double *xsp1, int nsp1, double *zsp1, double *xsp2, int nsp2, double *zsp2, double xstep2)
{
    double SSE,*SSE_1, freqWJ,cf;
    int count,j,index;
    double A1, A2, A3;

    cf = 2*xval;

    SSE_1 = calloc( nsp1, sizeof(double) );
    count = 0;

#ifdef SKIP
/* Matlab code from Anup */
for i = 1:W
    for j = 1:H
        if(mask4all(i,j)>0)
            zSpectW = squeeze(double(ImgfreqW1Smooth(i,j, :)));
            Ref_zSpectW = zSpectW(nufreqW:-1:1);  %%%% reflect zSpect about centeral freq or ref freq
            IntSp_RefzSpectW = spline(tW, Ref_zSpectW(:),tsW);
            [minVal, minIndex] = min(zSpectW);
            start_c = -freqW(minIndex);  %%careful (-ve sign)--initial guess
            c_bnd = (start_c-0.2): IntStep : (start_c +0.2); % in ppm
            SSE_W = zeros(1, length(c_bnd));

            for k = 1: length(c_bnd)
                SSE_W(k) = symWASSR(c_bnd(k), freqW, nufreqW, zSpectW, IntSp_RefzSpectW, tsNew, f_N);
            end
            [minSSE, minSSEIndex] =  min(SSE_W);
            EstiB0shift(i,j) = -c_bnd(minSSEIndex);
        else
            EstiB0shift(i,j) = 0;
        end
    end
end

function SSE = symWASSR(c_bnd, freqW, nufreqW, vect1, ImgIntSp, tsNew, f_N)
    SSE = 1000000.0;
    cW =2*c_bnd;
    SSE_1 = zeros(1, nufreqW);
    count = 0;
    if(cW<0)
         for j = 1: nufreqW
             freqWJ = (freqW(j) + cW);
            if(freqWJ >freqW(1) && freqWJ <(f_N + cW))
               count = j;
               index = tsNew == int32(10000*freqWJ);
               SSE_1(j) = (vect1(j) - ImgIntSp(index)).^2;
            end
         end
    else
        for j = 1: nufreqW
            freqWJ = (freqW(j) + cW);
            if(freqWJ >(freqW(1)+ cW) && freqWJ <f_N)
                count = j;
                index = tsNew == int32(10000*freqWJ);
                SSE_1(j) = (vect1(j) - ImgIntSp(index)).^2;
            end
        end
    end
    if (count>0)
        SSE = sum(SSE_1);
    end
end

#endif

/*  Begin Kim Paper - Magnetic Resonance in Medicine 61:1441-1450 (2009)*/

    for (j = 0; j < nsp1; j++)
    {
        A1 = xsp1[j] + cf;
        if (cf < 0)
        {
            A2 = xsp1[0];
            A3 = cf + xsp1[nsp1-1];
        }
        else
        {
            A2 = cf + xsp1[0];
            A3 = xsp1[nsp1-1];
        }
        if ( A1 > A2 && A1 < A3 )
        {
            index = (int)( (A1 - xsp2[0])/xstep2 );
            SSE_1[count] = pow( (zsp1[j] - zsp2[index]),2 );
            count++;
        }
    }

/*  End  Kim WASSR Paper- Magnetic Resonance in Medicine 61:1441?1450 (2009) */


    if (count > 0)
    {
        SSE = 0.0;
        for (j = 0; j < count; j++)
            SSE += SSE_1[j];
    }
    else
        SSE = 1.0e38;
    free (SSE_1);
    return(SSE);
}

double * calc_b0_map(double *ppmlist, double *pimg, double *nimg, long *mask,
                     float b0ppmstep, int nppm, int nelemimg, int nelem)
{

    int intfactor, ielem, ippm1, ippm2;
    const mwSize *dim_array_mask, *dim_array_image;
    int i1, i2,nzsp1, nzsp2, rangezsp2, minzsp2index, minSSEindex, W, H;
    double maxppm, minppm, stepppm, stepppm2, stepppm2inv, ftemp1, ftemp2, maxpimg, maxnimg;
    double zspymin, zspxmin; // *B0map,
    double temp1,temp2;


    maxnimg = -1.0E38;
    maxpimg = -1.0E38;
    for ( i1 = 0; i1 < nelemimg; i1++)
    {
        if (pimg[i1] > maxpimg)
            maxpimg = pimg[i1];
        if (nimg[i1] > maxnimg)
            maxnimg = nimg[i1];
    }


    double *B0map =  (double *)malloc(nelem * sizeof(double));

    double *ppmlistasc = (double *)malloc(nppm * sizeof(double));
    int *psort = (int *)malloc(nppm * sizeof(int));
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


    intfactor = floor(stepppm/b0ppmstep);
    nzsp1 = 2*nppm-1;
    nzsp2 = intfactor*(nzsp1-1)+1;
    stepppm2 = (stepppm/intfactor);
    stepppm2inv = (intfactor/stepppm);
    rangezsp2 = floor( 0.2*stepppm2inv );

    double *zspp1 = (double *)malloc(nppm * sizeof(double));
    double *zspn1 = (double *)malloc(nppm * sizeof(double));
    double *zsp1 = (double *)malloc(nzsp1 * sizeof(double));  /* raw z spectrum sorted */
    double * xzsp = (double *)malloc(nzsp1 * sizeof(double));  /* raw z spectrum  cordinate */
    double * y2 = (double *)malloc(nzsp1 * sizeof(double));  /* storage for second derivative */

    double * xzsp2 = (double *)malloc(nzsp2 * sizeof(double));  /* Interpolated Zspetrum ppm coordinates */
    double * xzsp2r = (double *)malloc(nzsp2* sizeof(double));  /* Interpolated reflected Zspetrum ppm coordinates */
    double * zsp2 = (double *)malloc(nzsp2 * sizeof(double));  /* Interpolated Zspetrum */
    double * zsp2r = (double *)malloc(nzsp2 * sizeof(double));  /* Interpolated reflected Zspetrum */


    for (ippm1 = 0; ippm1 < nppm; ippm1++)
    {
        ippm2 = ippm1 + nppm -1;
        xzsp[ippm1] = -ppmlist[nppm-1-psort[ippm1]];
        xzsp[ippm2] = ppmlist[psort[ippm1]];
    }

    for (i2 =0; i2< nzsp2; i2++)
    {
        xzsp2[i2] = -maxppm + i2 * stepppm2;
        xzsp2r[nzsp2-1-i2] = xzsp2[i2];
    }

    for (ielem = 0; ielem < nelem; ielem++)
    {
        B0map[ielem] = 0.0;
        if (mask[ielem] > 0) {

            /* Create z spectrum from -maxppm to +maxppm; */
            for  (ippm1 = 0; ippm1 < nppm; ippm1++)
            {
                i1 = ippm1*nelem + ielem;
                zspp1[ippm1] = pimg[i1];
                zspn1[ippm1] = nimg[i1];
            }
            for (ippm1 = 0; ippm1 < nppm; ippm1++)
            {
                ippm2 = ippm1 + nppm -1;
                zsp1[ippm1] = zspn1[nppm-1-psort[ippm1]];
                zsp1[ippm2] = zspp1[psort[ippm1]];
            }

            /* Find zspectrum minimum position and value */
            zspymin = 1.0e38;
            for (ippm1 = 0; ippm1 < nzsp1; ippm1++)
            {
                if (zsp1[ippm1] < zspymin)
                {
                    zspymin     = zsp1[ippm1];
                    zspxmin     = xzsp[ippm1];
                }
            }

            /* Call spline to get second derivatives */
            spline(xzsp, zsp1, nzsp1, zsp1[0], zsp1[nzsp1-1], y2);
            /* Call splint  for interpolations and find index for zsp2 minimum */
            zspymin = 1.0e38;
            for (i2 =0; i2< nzsp2; i2++)
            {
                splint(xzsp, zsp1, y2, nzsp1, xzsp2[i2], &zsp2[i2]);
                if (zsp2[i2] < zspymin)
                {
                    zspymin     = zsp2[i2];
                    zspxmin     = xzsp2[i2];
                }

                zsp2r[i2] = zsp2[nzsp2-1-i2];
            }
            /* Apply maximum symmetry algorithm to find B0 map position */
            temp1 = -zspxmin;
            ftemp2 = 1.0E30;
            minSSEindex = -rangezsp2;
            for (i1 = -rangezsp2; i1 <= rangezsp2; i1++)
            {
                temp2 = temp1 + (i1 * stepppm2);
                ftemp1 = symWASSR(temp2, xzsp, nzsp1, zsp1, xzsp2, nzsp2, zsp2r,stepppm2);
                if (ftemp1 < ftemp2)
                {
                    minSSEindex = i1;
                    ftemp2 = ftemp1;
                }
            }
            B0map[ielem] = -( temp1 + (minSSEindex * stepppm2) );
        }
    }

    free(psort);
    free(ppmlistasc);
    free(zspp1);
    free(zspn1);
    free(y2);
    free(zsp1);
    free(zsp2);
    free(zsp2r);
    free(xzsp);
    free(xzsp2);
    free(xzsp2r);


    return B0map;
}
