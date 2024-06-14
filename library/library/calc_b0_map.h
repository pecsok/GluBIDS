void spline(double *x, double *y, int n, double yp1,
                   double ypn, double *y2);

void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)

double symWASSR(double xval, double *xsp1, int nsp1, double *zsp1,
                       double *xsp2, int nsp2, double *zsp2, double xstep2);

double * calc_b0_map(double *ppmlist, double *pimg, double *nimg, long *mask,
					 float b0ppmstep, int nppm, int nelemimg, int nelem)