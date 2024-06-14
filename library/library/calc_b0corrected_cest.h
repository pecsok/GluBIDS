void nrerror(char *error_text)

double *vector(int nl,int nh)

int *ivector(int nl,int nh)

double *dvector(int nl,int nh)

double **matrix(int nrl,int nrh,int ncl,int nch)

double **dmatrix(int nrl,int nrh,int ncl,int nch)

int **imatrix(int nrl,int nrh,int ncl,int nch)

double **submatrix(double **a, int oldrl, int oldrh,int oldcl,int oldch,int newrl,int newcl)

void free_vector(double *v,int nl,int nh)

void free_ivector(int *v,int nl,int nh)

void free_matrix(double **m,int nrl,int nrh,int ncl,int nch)

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)

void free_dvector(double *v,int nl,int nh)

void free_submatrix(double **m,int nrl,int nrh,int ncl,int nch)

double **convert_matrix(double *a, int nrl, int nrh, int ncl,int nch)

void free_convert_matrix(double **b, int nrl, int nrh, int ncl, int nch)

void gaussj(double **a, int n,double **b, int m)

void covsrt(double **covar, int ma,int *lista, int mfit)

void polyfit(double *x, double *y, double *sig, int ndata, double *a, int ma, int *lista, int mfit,
          double **covar,double *chisq, void (*funcs)(double,double *,int) )

void xpoly(double x,double *afunc,int mma)

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)

void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)

double * calc_b0corrected_cest(double ppmval, double *ppmlist, double *pimg,
  double *nimg, double *B0map, long *mask, double *refimg, int nppm, int nelem)
