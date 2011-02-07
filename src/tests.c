#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/* ***********************************************************************
 *                                                                       *
 *                          Declaration of functions                     *
 *                                                                       *
 * ********************************************************************* */



/* Functions coming from the package ade4 */

void vecpermut (double *A, int *num, double *B);
double alea (void);
void aleapermutvec (double *a);
void trirapideintswap (int *v, int i, int j);
void trirapideint (int *x , int *num, int gauche, int droite);
void sqrvec (double *v1);
void getpermutation (int *numero, int repet);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAAtB (double **a, double **b);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
void taballoc (double ***tab, int l1, int c1);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void freetab (double **tab);
void freevec (double *vec);
void freeintvec (int *vec);
void matcentrage (double **A, double *poili, char *typ);
void matmodiffc (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcm (double **tab, double *poili);
void DiagobgComp (int n0, double **w, double *d, int *rang);





/* Functions from the package adehabitat */
void rpath(double **xp, double *rcx, double *rcy, double **asc, 
	   double **tabdist, double *dt, 
	   double *angles, double *xc, double *yc,
	   double *cs, int r);
void randpath(double *xpr, double *rcrx, double *rcry, double *ascr, 
	      double *xcr, double *ycr, double *csr,
	      double *tabdistr, double *dtr, double *anglesr, 
	      int *nlasc, int *ncasc, int *nltdr, int *nlocsr);
void dtmp(double x1, double x2, double y1, double y2, 
	  double *di);
void fptt(double *x, double *y, double *t, int pos, double radius, 
	  double *fptto, int nlo);
void fipati(double *x, double *y, double *t, 
	    int nlo, int nra, double *rad, 
	    double **fpt);
void fipatir(double *xr, double *yr, double *tr, 
	     int *nlocs, double *radius, int *nrad, 
             double *fptr);
void perclu(double **map, int nr, int nc, double *x, double *y,
	    int nmax, int *nreel, double *pm);
void perclur(double *mapr, int *nrm, int *ncm, double *probamr,
	     double *xr, double *yr, int *nmaxr, int *nreel);
void resolpol(double a, double b, double c, double *x1, double *x2, 
	      int *warn);
void discretraj(double *x, double *y, double *dat, double *xn, 
		double *yn, int n, int nn, double *datn, 
		double u, int *neff);
void discretrajr(double *xr, double *yr, double *datr, double *xnr, 
		 double *ynr, int *nr, int *nnr, double *datnr, 
		 double *xdeb, double *ydeb, double *ur, double *dat0, 
		 int *neff);
void permutR2n(double *xyr, int *nro, int *nrepr, 
	       double *R2nr, double *dtr, double *dtsimr);
void runsltr(int *xr, int *nr, double *res, int *nrepr);
void testindepangl (double *sim, double *ang, int *nang, int *debut, 
                    int *fin, int *ndeb, int *ni);
void testindepdist (double *sim, double *di, int *ndi, int *debut, 
		    int *fin, int *ndeb, int *ni);
void prepquart (double *dtur, int *nur, double *dtrr, double *resr, int *nrr, 
		int *ncr, double *rur);
void optcut (double **Pid, double **mk, int *maxk);
void optcutr (double *Pidr, double *mkr, int *maxk, int *lr, int *pr, 
	      double *baye);
void partraj(double **Pid, int *maxk, double **Mk, double **Mkd, 
	     double **res);
void partrajr(double *Pidr, double *curmar, int *curmodr, int *curlocr, 
	      int *lr, int *Dr, int *Kmr);
void acfdist (double *sim, double *di, int *ndi, double *diper, int *nbobs, 
	      int *indxnona, int *nsim, int *lag);
void acfangl (double *sim, double *ang, int *nang, double *angper, 
	      int *nbobs, int *indxnona, int *nsim, int *lag);

/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of ADE-4                    *****
 *********               --------------------                    *****
 *********************************************************************
 *********************************************************************
 */



/**************************/
double alea (void)
{
    double w;
    w = ((double) rand())/ (double)RAND_MAX;
    return (w);
}

/*************************/
void aleapermutvec (double *a)
{
    /* Randomly permutes the elements of a vector a
       Manly p. 42 The vector is modified
       from Knuth 1981 p. 139 */
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
	j=lig-i+1;
	k = (int) (j*alea()+1);
	/* k = (int) (j*genrand()+1); */
	if (k>j) k=j;
	z = a[j];
	a[j]=a[k];
	a[k] = z;
    }
}


/*******************/	
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
 * A is a vector with n elements
 * B is a vector with n elements
 * num is a random permutation of the n first integers
 * B contains in output the permuted elements of A
 * ---------------------------------------*/
    
    int lig, lig1, lig2, i, k;
    
    lig = A[0];
    lig1 = B[0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (lig!=lig2) ) {
	/* err_message ("Illegal parameters (vecpermut)");
	   closelisting(); */
    }
    
    for (i=1; i<=lig; i++) {
	k=num[i];
	B[i] = A[k];
    }
}

/********* Centring accrding to row weights poili **********/	
void matcentrage (double **A, double *poili, char *typ)
{
    
    if (strcmp (typ,"nc") == 0) {
	return;
    } else if (strcmp (typ,"cm") == 0) {
	matmodifcm (A, poili);
	return;
    } else if (strcmp (typ,"cn") == 0) {
	matmodifcn (A, poili);
	return;
    } else if (strcmp (typ,"cp") == 0) {
	matmodifcp (A, poili);
	return;
    } else if (strcmp (typ,"cs") == 0) {
	matmodifcs (A, poili);
	return;
    } else if (strcmp (typ,"fc") == 0) {
	matmodiffc (A, poili);
	return;
    } else if (strcmp (typ,"fl") == 0) {
	matmodifcm (A, poili);
	return;
    }
}

/*********************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a complete disjonctive table with n rows and m columns
 * poili is a vector with n components
 * The process returns tab centred by column
 * with weighting poili (sum=1)
 * centring type multple correspondances
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    for (i=1;i<=l1;i++) tab[i][j] = 0;
	} else {
	    
	    for (i=1;i<=l1;i++) {
		z = tab[i][j]/x - 1.0;
		tab[i][j] = z;
	    }
	}
    }
    freevec (poimoda);
}

/*********************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows and p columns
 * poili is a vector with n components
 * the function returns tab normed by column
 * with the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid, x, z, y, v2;
    int 			i, j, l1, c1;
    double		*moy, *var;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    
    vecalloc(&moy, c1);
    vecalloc(&var, c1);
    
    
/*--------------------------------------------------
 * centred and normed table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    for (i=1;i<=l1;i++) {
	poid=poili[i];
	for (j=1;j<=c1;j++) {
	    x = tab[i][j] - moy[j];
	    var[j] = var[j] + poid * x * x;
	}
    }
    
    for (j=1;j<=c1;j++) {
	v2 = var[j];
	if (v2<=0) v2 = 1;
	v2 = sqrt(v2);
	var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	y = var[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    z = z / y;
	    tab[j][i] = z;
	}
    }
    
    freevec(moy);
    freevec(var);
    
}

/*********************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows, p columns
 * poili is a vector with n components
 * The function returns tab standardised by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
	double		x,poid, z, y, v2;
	int 			i, j, l1, c1;
	double		*var;
	
	l1 = tab[0][0];
	c1 = tab[1][0];
	vecalloc(&var, c1);
	

/*--------------------------------------------------
 * calculation of the standardised table
 --------------------------------------------------*/
	
	for (i=1;i<=l1;i++) {
	    poid=poili[i];
	    for (j=1;j<=c1;j++) {
		x = tab[i][j];
		var[j] = var[j] + poid * x * x;
	    }
	}
	
	for (j=1;j<=c1;j++) {
	    v2 = var[j];
	    if (v2<=0) v2 = 1;
	    v2 = sqrt(v2);
	    var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
	    y = var[i];
	    for (j=1;j<=l1;j++) {
		z = tab[j][i];
		z = z / y;
		tab[j][i] = z;
	    }
	}
	freevec(var);
}


/**********/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and p colonnes
 * poili is a vector with n components
 * The function returns tab centred by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, c1;
    double		*moy, x, z;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    vecalloc(&moy, c1);
    
    
/*--------------------------------------------------
 * Centred table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    tab[j][i] = z;
	}
    }
    freevec(moy);
}

/*********************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and m columns
 * of number >=0
 * poili is a vector with n components
 * The function returns tab doubly centred
 * for the weighting poili (sum=1)
 * centring type simple correspondance analysis
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	x = 0;
	for (j=1;j<=m1;j++) {
	    x = x + tab[i][j];
	}
	if (x!=0) {
	    for (j=1;j<=m1;j++) {
		tab[i][j] = tab[i][j]/x;
	    }
	}	
    }
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    /* err_message("column has a nul weight (matmodiffc)"); */
	}
	
	for (i=1;i<=l1;i++) {
	    z = tab[i][j]/x - 1.0;
	    tab[i][j] = z;
	}
    }
    freevec (poimoda);
}









/*****************/
void getpermutation (int *numero, int repet)
/*----------------------
 * affects a random permutation of the first n integers
 * in an integer vector of length n
 * First vecintalloc is needed
 * *numero is a vector of integer
 * repet is an integer which can take any arbitrary value
 * used in the seed of the pseudo-random number generation process
 * if it is increased in repeated calls (e.g. simulation), it is ensured that
 * two calls returns different results (seed=clock+repet)
 ------------------------*/
{
    int i, n, seed;
    int *alea;
    
    n=numero[0];
    vecintalloc (&alea,n);
    
    /*-------------
     * numbering in numero
     -----------*/
    for (i=1;i<=n;i++) {
	numero[i]=i;
    }
    
    /*-------------
     * affects random numbers in alea
     ----------------*/
    seed = clock();
    seed = seed + repet;
    srand(seed);
    for (i=1;i<=n;i++) {
	alea[i]=rand();
    }
    
    trirapideint (alea , numero, 1, n);
    freeintvec (alea);
}

/*****************************************/
/* Sorting: used in getpermutation */

void trirapideint (int *x , int *num, int gauche, int droite)
{
    int j, dernier, milieu, t;
    
    if ( (droite-gauche)<=0) return;
    
    milieu = (gauche+droite)/2;
    trirapideintswap (x, gauche, milieu);
    trirapideintswap (num, gauche, milieu);
    
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
	if (x[j] < t) {
	    dernier = dernier + 1;
	    trirapideintswap (x, dernier, j);	
	    trirapideintswap (num, dernier, j);
	}
    }
    trirapideintswap (x, gauche, dernier);
    trirapideintswap (num, gauche, dernier);
    
    trirapideint (x, num, gauche, dernier-1);
    trirapideint (x, num, dernier+1, droite);
    
}

/**************************************/
/* Sorting: used in trirapideint */

void trirapideintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}

/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
 * Square root of the elements of a vector
 --------------------------------------------------*/
{
    int i, c1;
    double v2;
    
    c1 = v1[0];
    
    for (i=1;i<=c1;i++) {
	v2 = v1[i];
	/* if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)"); */
	v2 = sqrt(v2);
	v1[i] = v2;
    }
}

/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
 * Eigenstructure of a matrix. See
 * T. FOUCART Analyse factorielle de tableaux multiples,
 * Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
 * de LEBART et coll.
 --------------------------------------------------*/
{
    double			*s;
    double			a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double			dble;
    int				ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    ni = 100;
    if (n0 == 1) {
	d[1] = w[1][1];
	w[1][1] = 1.0;
	*rang = 1;
	freevec (s);
	return;
    }
    
    for (i2=2;i2<=n0;i2++) {
	
	b=0.0;
	c=0.0;
	i=n0-i2+2;
	k=i-1;
	if (k < 2) goto Et1;
	for (l=1;l<=k;l++) {
	    c = c + fabs((double) w[i][l]);
	}
	if (c != 0.0) goto Et2;
	
    Et1:	s[i] = w[i][k];
	goto Etc;
	
    Et2:	for (l=1;l<=k;l++) {
	x = w[i][l] / c;
	w[i][l] = x;
	b = b + x * x;
    }
	xp = w[i][k];
	ix = 1;
	if (xp < 0.0) ix = -1;
		
/*		q = -sqrt(b) * ix; */
	dble = b;
	dble = -sqrt(dble);
	q = dble * ix;
	
	s[i] = c * q;
	b = b - xp * q;
	w[i][k] = xp - q;
	xp = 0;
	for (m=1;m<=k;m++) {
	    w[m][i] = w[i][m] / b / c;
	    q = 0;
	    for (l=1;l<=m;l++) {
		q = q + w[m][l] * w[i][l];
	    }
	    m1 = m + 1;
	    if (k < m1) goto Et3;
	    for (l=m1;l<=k;l++) {
		q = q + w[l][m] * w[i][l];
	    }
	    
	Et3:		s[m] = q / b;
	    xp = xp + s[m] * w[i][m];
	}
	bp = xp * 0.5 / b;
	for (m=1;m<=k;m++) {
	    xp = w[i][m];
	    q = s[m] - bp * xp;
	    s[m] = q;
	    for (l=1;l<=m;l++) {
		w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
	    }
	}
	for (l=1;l<=k;l++) {
	    w[i][l] = c * w[i][l];
	}
	
    Etc:	d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
	
	k = i - 1;
	if (d[i] == 0.0) goto Et4;
	for (m=1;m<=k;m++) {
	    q = 0.0;
	    for (l=1;l<=k;l++) {
		q = q + w[i][l] * w[l][m];
	    }
	    for (l=1;l<=k;l++) {
		w[l][m] = w[l][m] - q * w[l][i];
	    }
	}
	
    Et4:	d[i] = w[i][i];
	w[i][i] = 1.0;
	if (k < 1) goto Et5;
	for (m=1;m<=k;m++) {
	    w[i][m] = 0.0;
	    w[m][i] = 0.0;
	}
	
    Et5:;
    }
    
    for (i=2;i<=n0;i++) {
	s[i-1] = s[i];
    }
    s[n0] = 0.0;
    
    for (k=1;k<=n0;k++) {
	
	m = 0;
	
    Et6: 	for (j=k;j<=n0;j++) {
	if (j == n0) goto Et7;
	ab = fabs((double) s[j]);
	ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
	if (ab < ep) goto Et7;
    }
	
    Et7: 	isnou = 1;
	h = d[k];
	if (j == k) goto Eta;
	if (m < ni) goto Etd;
	
	/* err_message("Error: can't compute matrix eigenvalues"); */
	
    Etd:	m = m + 1;
	q = (d[k+1]-h) * 0.5 / s[k];
	
/*		t = sqrt(q * q + 1.0); */
	dble = q * q + 1.0;
	dble = sqrt(dble);
	t = dble;
	
	if (q < 0.0) isnou = -1;
	q = d[j] - h + s[k] / (q + t * isnou);
	u = 1.0;
	v = 1.0;
	h = 0.0;
	jk = j-k;
	for (ijk=1;ijk<=jk;ijk++) {
	    i = j - ijk;
	    xp = u * s[i];
	    b = v * s[i];
	    if (fabs((double) xp) < fabs((double) q)) goto Et8;
	    u = xp / q;
	    
/*			t = sqrt(u * u + 1); */
	    dble = u * u + 1.0;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = q * t;
	    v = 1 / t;
	    u = u * v;
	    goto Et9;
	    
	Et8:		v = q / xp;
	    
/*			t = sqrt(1 + v * v); */
	    dble = 1.0 + v * v;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = t * xp;
	    u = 1 / t;
	    v = v * u;
	    
	Et9:
	    q = d[i+1] - h;
	    t = (d[i] - q) * u + 2.0 * v * b;
	    h = u * t;
	    d[i+1] = q + h;
	    q = v * t - b;
	    for (l=1;l<=n0;l++) {
		xp = w[l][i+1];
		w[l][i+1] = u * w[l][i] + v * xp;
		w[l][i] = v * w[l][i] - u * xp;
	    }
	}
	d[k] = d[k] - h;
	s[k] = q;
	s[j] = 0.0;
	
	goto Et6;
	
    Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
	
	i = ij - 1;
	l = i;
	h = d[i];
	for (m=ij;m<=n0;m++) {
	    if (d[m] >= h) {
		l = m;
		h = d[m];
	    }
	}
	if (l == i) {
	    goto Etb;
	} else {
	    d[l] = d[i];
	    d[i] = h;
	}
	for (m=1;m<=n0;m++) {
	    h = w[m][i];
	    w[m][i] = w[m][l];
	    w[m][l] = h;
	}
	
    Etb:;
    } /* for (ij=2;ij<=n0;ij++) */
    
    /* final:; */
    *rang = 0;
    for (i=1;i<=n0;i++) {
	/*
	  if (d[i] / d[1] < 0.00001) d[i] = 0.0;
	  if (d[i] != 0.0) *rang = *rang + 1;
	*/
	if (d[i] > 0.0) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoCompbg */







/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Matrix product AB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (i=1;i<=lig;i++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (j=1;j<=col;j++) {
		s = s + a[i][j] * b[j][k];
	    }
	    c[i][k] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Matrix product AtA
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=j;k<=col;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][k] * a[i][j];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
 * Matrix product AtB
 --------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][j] * b[i][k];
	    }
	    c[j][k] = s;
	}		
    }
}


/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
 * Matrix product B = AAt
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=lig;j++) {
	for (k=j;k<=lig;k++) {
	    s = 0;
	    for (i=1;i<=col;i++) {
		s = s + a[j][i] * a[k][i];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/*******************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
 * Produit matriciel AtB
 * les lignes de B sont permutÚes par la permutation permut
 --------------------------------------------------*/
{
    int j, k, i, i0, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		i0 = permut[i];
		s = s + a[i][j] * b[i0][k];
	    }
	    c[j][k] = s;
	}		
    }
}

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
 * Dynamic Memory Allocation for a table (l1, c1)
 --------------------------------------------------*/
{
    int i, j;
    
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
	for (i=0;i<=l1;i++) {
	    if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
		return;
		for (j=0;j<i;j++) {
		    free(*(*tab+j));
		}
	    }
	}
    }
    
    **(*tab) = l1;
    **(*tab+1) = c1;
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
 * Memory Allocation for a vector of length n
 --------------------------------------------------*/
{
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/*****************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
 * Memory allocation for an integer vector of length  n
 --------------------------------------------------*/
{
    if ( (*vec = (int *) calloc(n+1, sizeof(int))) != NULL) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
 * Free memory for a table
 --------------------------------------------------*/
{
    int 	i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
	free((char *) *(tab+i) );
    }
    free((char *) tab);
}

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
 * Free memory for a vector
 --------------------------------------------------*/
{
    free((char *) vec);	
}

/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* Free memory for an integer  vector
--------------------------------------------------*/
{
    
    free((char *) vec);
    
}














/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of adehabitat               *****
 *********               -------------------------               *****
 *********************************************************************
 *********************************************************************
 */



void dedans(double *pts, double *xc, double *yc, double *na,
	    double cs, double **asc)
{
    int nl, nc, i, ligne, colo;
    double x, y;
  
    x = pts[1];
    y = pts[2];
    
    nl = xc[0];
    nc = yc[0];
    
    ligne = 0;
    colo = 0;
    
    for (i = 1; i <= nl; i++) {
	if (((xc[i] - cs/2) <= x) && ((xc[i] + cs/2) > x))
	    ligne = i;
    }
    
    for (i = 1; i <= nc; i++) {
	if (((yc[i] - cs/2) <= y) && ((yc[i] + cs/2) > y))
	    colo = i;
    }
    *na = asc[ligne][colo]; 
}





/* ****************************************************************
   *                                                              *
   *   Randomization of a traject, based on the distances between *
   *   successive relocations and the time lag, as well as the    *
   *   angles between successive steps                            *
   *                                                              *
   **************************************************************** */

void rpath(double **xp, double *rcx, double *rcy, double **asc, 
	   double **tabdist, double *dt, 
	   double *angles, double *xc, double *yc,
	   double *cs, int r)
{
    /* Declaration */
    int i, j, k, l, m, nsam, nltd, *index, *indangles, nlocs;
    double *pts, na, interv, *dobs, dech, anglech, ang;
    
    /* Memory allocation */
    vecalloc(&pts, 2);
    nltd = tabdist[0][0];
    nlocs = xp[0][0];
    vecintalloc(&indangles, (nlocs-2));
    
    /* fills the vector indangle */
    for (i = 1; i <= (nlocs - 2); i++) {
	indangles[i] = i;
    }
    
    /* 1. First loc of the traject */
    k = 0;
    ang = 0;
    
    while (k==0) {
	
	/* Random draw of relocations coordinates */
	xp[1][1] = (alea())*(rcx[2]-rcx[1]) + rcx[1];
	xp[1][2] = (alea())*(rcy[2]-rcy[1]) + rcy[1];
	
	pts[1] = xp[1][1];
	pts[2] = xp[1][2];
	
	/* Verifies that the loc is inside the study area */
	dedans(pts, xc, yc, &na, *cs, asc);
	if (fabs(na + 9999) > 0.000000001)
	    k = 1;
	
    }
    
    
    /* loop for the following relocations */
    for (i = 1; i <= (nlocs-1); i++) {
	interv = dt[i];
	
	/* How many distances for the observed time lag ? */
	nsam = 0;
	for (j = 1; j <= nltd; j++) {
	    if (fabs(tabdist[j][1] - interv) < 0.000000001)
		nsam++;
	}
	
	/* Table of distances */
	vecalloc(&dobs, nsam);
	
	/* the vector index will be used to draw a random relocation */
	vecintalloc(&index, nsam);
	for (l = 1; l <= nsam; l++) {
	    index[l] = l;
	}
	
	/* In a first time, gets the corresponding distances */
	m = 1;
	for (j = 1; j <= nltd; j++) {
	    if (fabs(tabdist[j][1] - interv) < 0.000000001) {
		dobs[m] = tabdist[j][2];
		m++;
	    }
	}
	
	k = 0;
	while (k == 0) {
	    /* Sampled Distance */
	    r = (int) (alea() * 100);
	    getpermutation(index, j * r);
	    dech = dobs[index[1]];
	    
	    /* Sampled Angles */
	    getpermutation(indangles, j * r);
	    anglech = angles[indangles[1]];
	    
	    /* update the angles */
	    ang = (ang + anglech);
	    
	    /* The new coordinates */
	    xp[i+1][1] = xp[i][1] + dech * cos(ang);
	    xp[i+1][2] = xp[i][2] + dech * sin(ang);
	    
	    pts[1] = xp[i+1][1];
	    pts[2] = xp[i+1][2];
	    
	    dedans(pts, xc, yc, &na, *cs, asc);
	    if (fabs(na + 9999) > 0.000000001)
		k = 1;
	}
	freevec(dobs);
	freeintvec(index);
    }
    
    /* Free memory */
    freeintvec(indangles);
    freevec(pts);
}







/* Test of rpath with R */

void randpath(double *xpr, double *rcrx, double *rcry, double *ascr, 
	      double *xcr, double *ycr, double *csr,
	      double *tabdistr, double *dtr, double *anglesr, 
	      int *nlasc, int *ncasc, int *nltdr, int *nlocsr)
{
    /* declaration of the variables */
    int i, j, k, r, nlocs, nltd;
    double **xp, *rcx, *rcy, **asc, **tabdist, *dt, *angles;
    double *xc,*yc, cs;
    
    /* Memory allocation */
    nlocs = *nlocsr;
    nltd = *nltdr;
    cs = *csr;
    
    taballoc(&xp, nlocs, 2);
    vecalloc(&rcx, 2);
    vecalloc(&rcy, 2);
    taballoc(&asc, *nlasc, *ncasc);
    vecalloc(&xc, *nlasc);
    vecalloc(&yc, *ncasc);
    taballoc(&tabdist, nltd, 2);
    vecalloc(&dt, nlocs-1);
    vecalloc(&angles, nlocs-2);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xp[i][j] = xpr[k];
	    k++;
	}
    }
    rcx[1] = rcrx[0];
    rcx[2] = rcrx[1];
    rcy[1] = rcry[0];
    rcy[2] = rcry[1];
    
    k = 0;
    for (i = 1; i <= *nlasc; i++) {
	for (j = 1; j <= *ncasc; j++) {
	    asc[i][j] = ascr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= nltd; i++) {
	for (j = 1; j <= 2; j++) {
	    tabdist[i][j] = tabdistr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= (nlocs - 1); i++) {
	dt[i] = dtr[i-1];
    }
    for (i = 1; i <= *nlasc; i++) {
	xc[i] = xcr[i-1];
    }
    for (i = 1; i <= *ncasc; i++) {
	yc[i] = ycr[i-1];
    }
    for (i = 1; i <= (nlocs - 2); i++) {
	angles[i] = anglesr[i-1];
    }
    
    /* test of randpath */
    r = 1;
    rpath(xp, rcx, rcy, asc, tabdist, dt, angles,
	  xc, yc, &cs, r);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xpr[k] = xp[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(xp);
    freevec(rcx);
    freevec(rcy);
    freetab(asc);
    freetab(tabdist);
    freevec(dt);
    freevec(angles);
    freevec(xc);
    freevec(yc);
}




/* *********************************************************************
 *                                                                     *
 *                   First passage time                                *
 *                                                                     *
 ***********************************************************************/

/* compute the distance cetween two points */
void dtmp(double x1, double x2, double y1, double y2, 
	  double *di)
{
  *di = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}



/* compute the FPT for ONE relocation */
void fptt(double *x, double *y, double *t, int pos, double radius, double *fptto, int nlo)
{
    /* Declaration */
    int ok, pos2, naar, naav, na;
    double di, dt, dt2, di2, fptar, fptav;
    
    ok = 0;
    di = 0;
    di2 = 0;
    dt = 0;
    dt2 = 0;
    naar = 1;
    naav = 1;
    fptar = 0;
    fptav = 0;
  
  
    /* Search of the first loc outside the circle (before) */
    pos2 = pos;
    while (ok == 0) {
	pos2 = pos2 - 1;
	if (pos2 > 0) {
	    dtmp(x[pos2], x[pos], y[pos2], y[pos], &di);
	    if (di >= radius)
		ok = 1;
	} else {
	    ok = 1;
	    naar = 0;
	}
    }
    
    /* computes the linear approximation */
    if (naar > 0) {
	dt = abs(t[pos] - t[pos2]);
	dt2 = abs(t[pos] - t[(pos2+1)]);
	dtmp(x[(pos2+1)], x[pos], y[(pos2+1)], y[pos], &di2);
	fptar = dt2 + ( (dt - dt2) * (radius - di2) / (di - di2) );
    }
    
    
    /* Search of the first loc outside the circle (after) */
    pos2 = pos;
    ok = 0;
    while (ok == 0) {
	pos2 = pos2 + 1;
	if (pos2 <= nlo) {
	    dtmp(x[pos2], x[pos], y[pos2], y[pos], &di);
	    if (di >= radius)
		ok = 1;
	} else {
	    ok = 1;
	    naav = 0;
	}
    }
    
    /* Computes linear approximation */
    if (naav > 0) {
	dt = abs(t[pos2] - t[pos]);
	dt2 = abs(t[(pos2-1)] - t[pos]);
	dtmp(x[(pos2-1)], x[pos], y[(pos2-1)], y[pos], &di2);
	fptav = dt2 + ( (dt - dt2) * (radius - di2) / (di - di2) );
    }
    
    na = naar * naav;
    if (na > 0) {
	*fptto = fptar + fptav;
    } else {
	*fptto = -1;
    }
    
}



/* Computes the FPT for all relocations */
void fipati(double *x, double *y, double *t, 
	    int nlo, int nra, double *rad, 
	    double **fpt)
{
    /* Declaration */
    int i, j;
    double val;
    
    /* Computes the FPT */
    for (i=1; i<=nra; i++) {
	for (j=1; j<=nlo; j++) {
	    fptt(x, y, t, j, rad[i], &val, nlo);
	    fpt[j][i] = val;
	}
    }
}

/* for external call from within R */
void fipatir(double *xr, double *yr, double *tr, 
	     int *nlocs, double *radius, int *nrad, 
             double *fptr)
{
    /* Declaration */
    int i, j, k, nlo, nra;
    double *x, *y, *t, *rad, **fpt;
    
    /* Memory allocation */
    nlo = *nlocs;
    nra = *nrad;
    
    vecalloc(&x, nlo);
    vecalloc(&y, nlo);
    vecalloc(&t, nlo);
    vecalloc(&rad, nra);
    taballoc(&fpt, nlo, nra);
  
    /* R to C */
    for (i = 1; i <= nlo; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
	t[i] = tr[i-1];
    }
    
    for (i = 1; i <= nra; i++) {
	rad[i] = radius[i-1];
    }
    
    /* main function */
    fipati(x,y,t, nlo, nra, rad, fpt);
    
    /* C to R */
    k = 0;
    for (i=1; i<= nlo; i++) {
	for (j = 1; j<=nra; j++) {
	    fptr[k]=fpt[i][j];
	    k++;
	}
    }
    
    /* free memory */
    freetab(fpt);
    freevec(x);
    freevec(y);
    freevec(t);
    freevec(rad);
}






/* *********************************************************************
 *                                                                     *
 *                   Percolation cluster                               *
 *                                                                     *
 ***********************************************************************/


void perclu(double **map, int nr, int nc, double *x, double *y,
	    int nmax, int *nreel, double *pm)
{
    /* Declaration */
    int i,j, encore, k, l, dir, *vois, *rvois, *cvois, choix, xt, yt, cons, *reord, len;
    double **rr, **cc, *cs, ptir;
    
    /* Memory allocation */
    xt = (int) x[1];
    yt = (int) y[1];
    len = nr;
    if (nc < nr)
	len = nc;
    encore = 1;
    i = 1;
    j = 1;
    l = 0;
    k = 2;
    ptir = 0;
    dir = 1;
    choix = 1;
    cons = 0;
    
    taballoc(&rr, nr, nc);
    taballoc(&cc, nr, nc);
    vecintalloc(&vois, 4);
    vecintalloc(&reord, 4);
    vecintalloc(&rvois, 4);
    vecintalloc(&cvois, 4);
    vecalloc(&cs, 4);
    
    
    /* Rows and columns matrices */
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    rr[i][j] = (double) i;
	    cc[i][j] = (double) j;
	}
    }
    
    
    cs[1] = pm[1];
    for (i = 2; i <= 4; i++) {
	cs[i] = cs[i-1] + pm[i];
    }
    
    /* Beginning of the loop */
    while (encore == 1) {
    
	/* Storage of the neighbouring */
	vois[1] = (int) map[xt][yt+1];
	vois[2] = (int) map[xt+1][yt];
	vois[3] = (int) map[xt][yt-1];
	vois[4] = (int) map[xt-1][yt];
	
	rvois[1] = (int) rr[xt][yt+1];
	rvois[2] = (int) rr[xt+1][yt];
	rvois[3] = (int) rr[xt][yt-1];
	rvois[4] = (int) rr[xt-1][yt];
	
	cvois[1] = (int) cc[xt][yt+1];
	cvois[2] = (int) cc[xt+1][yt];
	cvois[3] = (int) cc[xt][yt-1];
	cvois[4] = (int) cc[xt-1][yt];
	
	/* Re-order of the neighbouring according to the direction */
	l = 1;
	for (i = dir; i <= 4; i++) {
	    reord[l] = i;
	    l++;
	}
	i = 1;
	while (l != 4) {
	    reord[l] = i;
	    l++;
	    i++;
	}
	
	/* random draw of a direction */
	ptir = alea();
	choix = 4;
	if (ptir <= pm[1]) {
	    choix = 1;
	}
	if ((ptir > pm[1])&&(ptir <= pm[2])) {
	    choix = 2;
	}
	if ((ptir > pm[2])&&(ptir <= pm[3])) {
	    choix = 3;
	}
	
	/* And again, until the direction lead us into an accessible area */
	cons = reord[choix];
	
	while (vois[cons] == 1) {
	    ptir = alea();
	    choix = 4;
	    if (ptir <= pm[1]) {
		choix = 1;
	    }
	    if ((ptir > pm[1])&&(ptir <= pm[2])) {
		choix = 2;
	    }
	    if ((ptir > pm[2])&&(ptir <= pm[3])) {
		choix = 3;
	    }
	    cons = reord[choix];
	}
	
	/* Stores all information */
	xt = (int) (rvois[choix]);
	yt = (int) (cvois[choix]);
	
	x[k] = (double) xt;
	y[k] = (double) yt;
	
	if ((xt==1)|(xt==len)|(yt==1)|(yt==len))
	    encore = 0;
	if (k==nmax)
	    encore = 0;
	
	for (i = 1; i <= 4; i++) {
	    if ( ((int) rvois[i]) == xt) {
		if ( ((int) cvois[i]) == yt) {
		    dir = i;
		}
	    }
	}
	
	*nreel = k;
	k++;
    }
    
    /* free memory */
    freeintvec(vois);
    freeintvec(rvois);
    freeintvec(cvois);
    freeintvec(reord);
    freevec(cs);
    freetab(rr);
    freetab(cc);
}





/* For external call from within R */
void perclur(double *mapr, int *nrm, int *ncm, double *probamr,
	     double *xr, double *yr, int *nmaxr, int *nreel)
{
    /* Declaration */
    double **map, *pm, *x, *y;
    int i, j, k, nr, nc, nmax;
    
    /* Memory allocation */
    nr = *nrm;
    nc = *ncm;
    nmax = *nmaxr;
    
    taballoc(&map, nr, nc);
    vecalloc(&x, nmax);
    vecalloc(&y, nmax);
    vecalloc(&pm, 4);
    
    /* R to C */
    x[1] = xr[0];
    y[1] = yr[0];
    
    k = 0;
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    map[i][j] = mapr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= 4; i++) {
	pm[i] = probamr[i-1];
    }
    
    /* Main function */
    perclu(map, nr, nc, x, y, nmax, nreel, pm);
    
    /* C to R */
    for (i = 1; i <= *nreel; i++) {
	xr[i-1] = x[i];
	yr[i-1] = y[i];
    }
    
    /* free memory */
    freevec(x);
    freevec(y);
    freevec(pm);
    freetab(map);
}



/* *********************************************************************
 *                                                                     *
 *         Rediscretization algorithm for a traject                    *
 *                                                                     *
 ***********************************************************************/



/* Resolves a quadratic equation of the type
 */

void resolpol(double a, double b, double c, double *x1, double *x2, int *warn)
{
    double delta;
    delta = (b * b) - 4 * a * c;
    *warn = 0;
    if (delta > 0) {
	*x1 = (-b - sqrt(delta)) / (2 * a);
	*x2 = (-b + sqrt(delta)) / (2 * a);
    } else {
	*warn = 1;
    }
}




void discretraj(double *x, double *y, double *dat, double *xn, 
		double *yn, int n, int nn, double *datn, 
		double u, int *neff)
{
    /* Declaration */
    double R, xt, yt, a, b, c, pente, ori, x1, x2, di1, di2;
    int k, m, p, fini, ok, warn, *dedans, lo, new, pp;
    
    /* memory allocation */
    fini = 0;
    k = 1;
    p = 2;
    m = 1;
    ok = 0;
    a = 0;
    b = 0;
    c = 0;
    pente = 0;
    ori = 0;
    x1 = 0;
    x2 = 0;
    lo = 0;
    di1 = 0;
    di2 = 0;
    *neff = 0;
    new = 0;
    pp = 1;
    
    
    vecintalloc(&dedans,2);
    
    /* Main algorithm */
    while (fini == 0) {
	
	dedans[1] = 0;
	dedans[2] = 0;
	ok = 0;
	xt = xn[k];
	yt = yn[k];
	k++;
	new = 0;
	
	/* Determines the "upper" point */
	while (ok == 0) {
	    if (new == 1)
		p++;
	    R = sqrt((((x[p] - xt) * (x[p] - xt)) + ((y[p] - yt) * (y[p] - yt))));
	    if (R > u) {
		ok = 1;
	    } else {
		if (p == n) {
		    fini = 1;
		    ok = 1;
		}
	    }
	    new = 1;
	}
	m = p-1;
    
    if (fini == 0) {
	/* Does the difference between x[p] and x[m] = 0? */
	if ((abs(x[p] - x[m]) > 0.000000000001)) {
	    /* Computes the slope between m and p */
	    pente = (y[p] - y[m]) / (x[p] - x[m]); /* when diff(x) == 0 ? */
	    /* The intercept */
	    ori = y[p] - (pente * x[p]);
	    /* The parameters of the polynomial equation */
	    a = 1 + (pente * pente);
	    b = (-2 * xt) + (2 * pente * ori) - (2 * pente * yt);
	    c = (xt * xt) + (yt * yt) + (ori * ori) - (2 * ori * yt) - (u * u);
	    resolpol(a, b, c, &x1, &x2, &warn);
	    /* 
	       A line cuts a circle with radius u at two points. One has 
	       (i) to identify the point the closest from m,n and 
	       (ii) to keep the one on the segment m-p
	    */
	    
	    
	    /* Which one are in the interval ? */
	    if (x1 >= x[m]) {
		if (x1 < x[p]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x1 >= x[p]) {
		if (x1 < x[m]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x2 >= x[m]) {
		if (x2 < x[p]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    if (x2 >= x[p]) {
		if (x2 < x[m]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    
	    /* What is the minimum distance to m ? */
	    if ((dedans[1] + dedans[2]) > 1) {
		di1 = fabs((double) (x[p] - x1));
		di2 = fabs((double) (x[p] - x2));
		
		/* verify that xk-1 is not in the same interval. Otherwise one increase of 1 */
		if (di1 < di2) {
		    lo = 2;
		} else {
		    lo = 1;
		}
		if (pp == p) {
		    if (di1 < di2) {
			lo = 1;
		    } 
		    if (di2 < di1) {
			lo = 2;
		    } 
		}
	    }
	    
	    /* storage of the coordinates */
	    if (lo == 1) {
		xn[k] = x1;
		yn[k] = (pente * x1) + ori;
	    }
	    if (lo == 2) {
		xn[k] = x2;
		yn[k] = (pente * x2) + ori;
	    }

	} else { /* We change x and y coordinates */
	    
	    /* Computes the slope between m and p */
	    pente =  (x[p] - x[m]) / (y[p] - y[m]);
	    /* The intercept */
	    ori = x[p] - (pente * y[p]);
	    /* The parameters of the polynomial equation */
	    a = 1 + (pente * pente);
	    b = (-2 * yt) + (2 * pente * ori) - (2 * pente * xt);
	    c = (xt * xt) + (yt * yt) + (ori * ori) - (2 * ori * xt) - (u * u);
	    resolpol(a, b, c, &x1, &x2, &warn);
	    /* 
	       A line cuts a circle with radius u at two points. One has 
	       (i) to identify the point the closest from m,n and 
	       (ii) to keep the one on the segment m-p
	    */
	    
	    
	    /* Which one are in the interval ? */
	    if (x1 >= y[m]) {
		if (x1 < y[p]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x1 >= y[p]) {
		if (x1 < y[m]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x2 >= y[m]) {
		if (x2 < y[p]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    if (x2 >= y[p]) {
		if (x2 < y[m]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    
	    /* What is the minimum distance to m ? */
	    if ((dedans[1] + dedans[2]) > 1) {
		di1 = fabs((double) (y[p] - x1));
		di2 = fabs((double) (y[p] - x2));
		
		/* verify that yk-1 is not in the same interval. Otherwise one increase of 1 */
		if (di1 < di2) {
		    lo = 2;
		} else {
		    lo = 1;
		}
		if (pp == p) {
		    if (di1 < di2) {
			lo = 1;
		    } 
		    if (di2 < di1) {
			lo = 2;
		    } 
		}
	    }
	    
	    /* storage of the coordinates */
	    if (lo == 1) {
		yn[k] = x1;
		xn[k] = (pente * x1) + ori;
	    }
	    if (lo == 2) {
		yn[k] = x2;
		xn[k] = (pente * x2) + ori;
	    }
	}
	
	/* Computes the nnew date (linear approximation) */
	di1 = sqrt((((xn[k] - x[m]) * (xn[k] - x[m])) + ((yn[k] - y[m]) * (yn[k] - y[m]))));
	R = sqrt((((x[p] - x[m]) * (x[p] - x[m])) + ((y[p] - y[m]) * (y[p] - y[m]))));
	di2 = dat[p] - dat[m];
	datn[k] = dat[m] + (di1 * di2 / R);
    }
    if (k == nn) {
	fini = 1;
    }
    pp = p;
    }
    
    /* Free memory */
    *neff = k;
    freeintvec(dedans);
}



/* For external Call from within R */

void discretrajr(double *xr, double *yr, double *datr, double *xnr, 
		 double *ynr, int *nr, int *nnr, double *datnr, 
		 double *xdeb, double *ydeb, double *ur, double *dat0, int *neff)
{
    /* Declaration */
    int i, n, nn;
    double *x, *y, *xn, *yn, *dat, *datn, u;
    
    /* Memory allocation */
    n = *nr;
    nn = *nnr;
    u = *ur;
    
    vecalloc(&x, n);
    vecalloc(&y, n);
    vecalloc(&xn, nn);
    vecalloc(&yn, nn);
    vecalloc(&dat, n);
    vecalloc(&datn, nn);
    
    /* R to C */
    for (i = 1; i <= n; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
	dat[i] = datr[i-1];
    }
    
    xn[1] = *xdeb;
    yn[1] = *ydeb;
    datn[1] = *dat0;
    
    /* Main function  */
    discretraj(x, y, dat, xn, yn, n, nn, datn, u, neff);
    
    /* C to R */
    for (i = 1; i <= nn; i++) {
	xnr[i-1] = xn[i];
	ynr[i-1] = yn[i];
	datnr[i-1] = datn[i];
    }
    
    /* Free memory */
    freevec(x);
    freevec(y);
    freevec(xn);
    freevec(yn);
    freevec(dat);
    freevec(datn);
}








void permutR2n(double *xyr, int *nro, int *nrepr, 
	       double *R2nr, double *dtr, double *dtsimr)
{
  double **xy, **R2n, *xp, *dt, **dtsim, tt;
  int n, i, j, k, *index, nr;
  
  n = *nro;
  nr = *nrepr;
  tt = 0;
  vecalloc(&xp, 2);
  vecalloc(&dt, n);
  taballoc(&R2n, n, nr);
  taballoc(&dtsim, n, nr);
  vecintalloc(&index, n);
  taballoc(&xy, n, 2);
  
  k = 0;
  for (i = 1; i <= n; i++) {
    dt[i] = dtr[i-1];
    for (j = 1; j <= 2; j++) {
      xy[i][j] = xyr[k];
      k++;
    }
  }
  
  for (k = 1; k <= nr; k++) {
    
    getpermutation(index, k);
    j = index[1];
    xp[1] = 0;
    xp[2] = 0;
    tt = 0;

    for (i = 1; i <= n; i++) {
      j = index[i];
      xp[1] = xp[1] + xy[j][1];
      xp[2] = xp[2] + xy[j][2];
      R2n[i][k] = (xp[2] * xp[2]) + (xp[1] * xp[1]);
      tt = tt + dt[j];
      dtsim[i][k] = tt;
    }
  }
  
  
  k = 0;
  for (i = 1; i <= n; i++) {
    for (j = 1; j<= nr; j++) {
      R2nr[k] = R2n[i][j];
      dtsimr[k] = dtsim[i][j];
      k++;
    }
  }
  
  
  freevec(xp);
  freevec(dt);
  freeintvec(index);
  freetab(xy);
  freetab(dtsim);
  freetab(R2n);
}




void runsltr(int *xr, int *nr, double *res, int *nrepr)
{
  int i, j, n, *x, *xb, nrep, nbsui, nz, nu, *numero;
  double m, s, n1, n2;
  
  n = *nr;
  nrep = *nrepr;
  
  vecintalloc(&x, n);
  vecintalloc(&xb, n);
  vecintalloc(&numero, n);

  
  for (i = 1; i <= n; i++) {
    x[i] = xr[i];
    numero[i] = i;
  }
  
  nbsui = 1;
  nz = 0;
  nu = 1;
  if (x[1] == 0)
    nz = 1;
  if (x[1] == 1)
    nu = 1;
  
  for (i = 2; i <= n; i++) {
    if (x[i] == 0)
      nz++;
    if (x[i] == 1)
      nu++;
    if (x[i-1] != x[i])
      nbsui++;
  }
  
  n1 = ((double) nz);
  n2 = ((double) nu);
  
  m = 1 + 2 * n1 * n2 / (n1 + n2);
  s = sqrt(2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)/((n1 + n2) * (n1 + n2) * (n1 + n2 - 1)));
  
  res[0] = (((double) nbsui) - m) / s;
  
  
  for (j = 1; j <= nrep; j++) {
    nbsui = 1;
    getpermutation(numero, j);
    trirapideint(numero , x, 1, n);
    for (i = 2; i <= n; i++) {
      if (x[i-1] != x[i])
	nbsui++;
    }
    res[j] = (((double) nbsui) - m) / s;
  }
  
  freeintvec(x);
  freeintvec(xb);
  freeintvec(numero);

}




void testindepangl (double *sim, double *ang, int *nang, int *debut, int *fin, int *ndeb, int *ni){
  
  int i,j,k;
  double *angle;
  vecalloc(&angle, *nang);
  for (i=1; i<=*nang; i++) {
    angle[i] = ang[i-1];
  }
  for (k=0;k<=(*ndeb-1);k++){
    for (i=debut[k]; i<=(fin[k]-1); i++) {
      sim[0]=sim[0]+(1.0-cos(angle[i+1] - angle[i]));
    }
  }
  for (j=1; j<=*ni; j++) {
    aleapermutvec(angle);
    for (k=0;k<=(*ndeb-1);k++){
      for (i=debut[k]; i<=(fin[k]-1); i++) {
	sim[j]=sim[j]+(1.0-cos(angle[i+1] - angle[i]));
      }
    }
  }
  for (j=0; j<=*ni+1; j++) {
    sim[j]=2*sim[j];

  }
  freevec(angle);
}

void testindepdist (double *sim, double *di, int *ndi, int *debut, int *fin, int *ndeb, int *ni){
  
  int i,j,k;
  double *dist;
  vecalloc(&dist, *ndi);
  for (i=1; i<=*ndi; i++) {
    dist[i] = di[i-1];
  }
  for (k=0;k<=(*ndeb-1);k++){
    for (i=(debut[k]); i<=(fin[k]-1); i++) {
      sim[0]=sim[0]+pow(dist[i+1] - dist[i],2);
    }
  }
  for (j=1; j<=*ni; j++) {
    aleapermutvec(dist);
    for (k=0;k<=(*ndeb-1);k++){
      for (i=(debut[k]); i<=(fin[k]-1); i++) {
	sim[j]=sim[j]+pow(dist[i+1] - dist[i],2);
      }
    }

  }
  freevec(dist);
}


void prepquart (double *dtur, int *nur, double *dtrr, double *resr, int *nrr, 
		int *ncr, double *rur)
{
  int i,j,k, nu, nr, nc, *rt;
  double *dtu, **dtr, **res, **ru, tmp1;
  
  nu = *nur;
  nr = *nrr;
  nc = *ncr;
  
  vecalloc(&dtu, nu);
  vecintalloc(&rt, nc);
  taballoc(&dtr, nr, nc);
  taballoc(&res, nr, nc);
  taballoc(&ru, nu, nc);
  
  for (i = 1; i <= nu; i++) {
    dtu[i] = dtur[i-1];
  }
  
  k = 0;
  for (i = 1; i <= nr; i++) {
    for (j = 1; j<= nc; j++) {
      dtr[i][j] = dtrr[k];
      res[i][j] = resr[k];
      k++;
    }
  }
  
  for (i = 1; i <= nu; i++) {
    for (k = 1; k <= nc; k++) {
      rt[k] = 0;
    }
    tmp1 = dtu[i];
    for (j = 1; j <= nr; j++) {
      for (k = 1; k <= nc; k++) {
	if ((fabs((double) (dtr[j][k] - tmp1)))< 0.0000000001)
	  rt[k] = j;
      }
    }
    for (k = 1; k <= nc; k++) {
      if ((fabs((double) rt[k] )) < 0.00000000001) {
	ru[i][k] = -1;
      } else {
	j = rt[k];
	ru[i][k] = res[j][k];
      }
    }
  }
  
  
  k = 0;
  for (i = 1; i <= nu; i++) {
    for (j = 1; j<= nc; j++) {
      rur[k] = ru[i][j];
      k++;
    }
  }
  
  freevec(dtu);
  freeintvec(rt);
  freetab(dtr);
  freetab(res);
  freetab(ru);
}





/* *********************************************************************
 *                                                                     *
 *                   partition d'un trajet                             *
 *                                                                     *
 ***********************************************************************/

void optcut (double **Pid, double **mk, int *maxk) 
{
    /* Declaration of variables */
    int l, p, i, ii, Km, k, j;
    double **mkd, **mkdn, tmp, mi1, mi2, mi3;
    
    
    
    /* memory allocation */
    l = Pid[0][0];   /* The length of the sequence */
    p = Pid[1][0];   /* The number of models */
    Km = *maxk;     /* the max number of partition to compute */
    i = 0;     /* the sequence index used in the paper (from 0 to l-1) */
    ii = 0;    /* the position index used in the program (from 1 to l) */
    j = 0;    /* the position index for the model (from 1 to p) */
    k = 0;    /* the index for the partition */
    taballoc(&mkd, l, p); /* The table for the log-probabilities of a 
			     partition given a model 
			     for a k partition
			     (to be re-used for each k) */
    taballoc(&mkdn, l, p); /* The table for the log-probabilities 
			      of a partition given a model
			     for a k+1 partition
			     (to be re-used for each k) */
    mi1 = 0;
    mi2 = 0;
    mi3 = 0;
    
    
    /* 1. First computes the probability of one partition*/
    
    /* 1.1 Fills the table mkd   */
    
    for (j = 1; j <= p; j++) {
	mkd[1][j] = log(Pid[1][j]);
    }
    for (j = 1; j <= p; j++) {
	for (ii = 2; ii <= l; ii++) {
	    mkd[ii][j] = mkd[ii-1][j] + log(Pid[ii][j]);
	}
    }
    
    /* 1.2 Fills the table mk */
    for (ii = 1; ii <= l; ii++) {
	
        /* computes the max */
	tmp = mkd[ii][1];
	for (j = 1; j <= p; j++) {
	    if (mkd[ii][j] - tmp > 0.0000000000001) {
		tmp = mkd[ii][j];
	    }
	}
	
	/* removes the max */
	for (j = 1; j <= p; j++) {
	    mkd[ii][j] = mkd[ii][j] - tmp;
	}
	
	/* computes the mean */
	mk[ii][1] = exp(mkd[ii][1]) / ( (double) p ) ;
	for (j = 2; j <= p; j++) {
	    mk[ii][1] = mk[ii][1] + (exp(mkd[ii][j]) / ( (double) p ) );
	}
	
        /* takes the log and adds the max again */
	mk[ii][1] = log(mk[ii][1]) + tmp;
	
	/* adds the max again for mkd */
	for (j = 1; j <= p; j++) {
	    mkd[ii][j] = mkd[ii][j] + tmp;
	}
    }
    
    
    
    
    /* 2. computes mk for each possible r-partitions: a loop */
    for (k = 2; k <= Km; k++) {
	
	
        /* First start at the first step */
	i = k-1;
	ii = k;
	
	/* Computes mkdn for each d and ii */
	for (j = 1; j <= p; j++) {
	    
	    /* Removes the max */
	    tmp = mk[ii-1][k-1];
	    if (mk[ii-1][k-1] - mkd[ii-1][j] < -0.00000000001)
		tmp = mkd[ii-1][j];
	    mi1 = mk[ii-1][k-1] - tmp;
	    mi2 = mkd[ii-1][j] - tmp;
	    
	    /* Computes for i = k-1 */
	    mkdn[ii][j] = log(Pid[ii][j]) + 
		log((((double) k) - 1) / (((double) i) * (((double) p) - 1))) +
		log(((double) p) * exp(mi1) - exp(mi2));
	    
	    /* adds the max again */
	    mkdn[ii][j] = mkdn[ii][j] + tmp;
	}
	


	/* fills mk */
	/* computes the max */
	tmp = mkdn[ii][1];
	for (j = 1; j <= p; j++) {
	    if (mkdn[ii][j] - tmp > 0.00000000000000001) {
		tmp = mkdn[ii][j];
	    }
	}
	
	/* removes the max */
	for (j = 1; j <= p; j++) {
	    mkdn[ii][j] = mkdn[ii][j] - tmp;
	}
	
	/* computes the mean */
	mk[ii][k] = exp(mkdn[ii][1]) / ( (double) p ) ;
	for (j = 2; j <= p; j++) {
	    mk[ii][k] = mk[ii][k] + (exp(mkdn[ii][j]) / ( (double) p ) );
	}
	
	/* takes the log and adds the max again */
	mk[ii][k] = log(mk[ii][k]) + tmp;
	
	/* adds the max again for mkd */
	for (j = 1; j <= p; j++) {
	    mkdn[ii][j] = mkdn[ii][j] + tmp;
	}
	
	
	/* Then increase the sequence */
	for (ii = (k + 1) ; ii <= l; ii++) {
	    
	    i = ii - 1;
	    
	    /* again computes mkdn for each d */
	    for (j = 1; j <= p; j++) {
		
		tmp = mkdn[ii-1][j];
		if (tmp - mk[ii-1][k-1] < -0.00000000000001) {
		    tmp = mk[ii-1][k-1];
		}
		if (tmp - mkd[ii-1][j] < -0.000000000000000001) {
		    tmp = mkd[ii-1][j];
		} 
		mi1 = mkdn[ii-1][j] - tmp;
		mi2 = mk[ii-1][k-1] - tmp;
		mi3 = mkd[ii-1][j] - tmp;
		
		mkdn[ii][j] = log(Pid[ii][j]) +
		    log((( ((double) (i - k + 1)) / ((double) i)) * 
			exp(mi1)) + (((((double) k) - 1) / 
				    (((double) i) * 
				     (((double) p) - 1))) * 
			((((double) p) * exp(mi2)) - 
			 exp(mi3)))) + tmp;
	    }


	    
	    /* ... and again fills mk */
	    /* computes the max */
	    tmp = mkdn[ii][1];
	    for (j = 1; j <= p; j++) {
		if (mkdn[ii][j] - tmp > 0.000000000000001) {
		    tmp = mkdn[ii][j];
		}
	    }
	    
	    /* removes the max */
	    for (j = 1; j <= p; j++) {
		mkdn[ii][j] = mkdn[ii][j] - tmp;
	    }
	    
	    /* computes the mean */
	    mk[ii][k] = exp(mkdn[ii][1]) / ( (double) p ) ;
	    for (j = 2; j <= p; j++) {
		mk[ii][k] = mk[ii][k] + (exp(mkdn[ii][j]) / ( (double) p ) );
	    }
	    
	    /* takes the log and adds the max again */
	    mk[ii][k] = log(mk[ii][k]) + tmp;
	    
	    /* adds the max again for mkd */
	    for (j = 1; j <= p; j++) {
		mkdn[ii][j] = mkdn[ii][j] + tmp;
	    }
	    
	}
	
	/* and finally, replace mkd with mkdn for next k */    
	for (ii = 1; ii <= l; ii++) {
	    for (j = 1; j <= p; j++) {
		mkd[ii][j] = mkdn[ii][j];
	    }
	}
    }
    
    
    /* free memory */
    freetab(mkd);
    freetab(mkdn);

}



/* Now, the R version */

void optcutr (double *Pidr, double *mkr, int *maxk, int *lr, int *pr, 
	      double *baye) 
{
    /* Declaration */
    int i, j, k, l, p, Km;
    double **Pid, **mk, tmp, *mik, *msk, *kk;
    
    /* Memory allocation */
    l = *lr;
    Km = *maxk;
    p = *pr;
    taballoc(&Pid, l, p);
    taballoc(&mk, l, Km);
    vecalloc(&mik, Km);
    vecalloc(&msk, Km);
    vecalloc(&kk, Km);

    
    /* Fills local variables */
    k = 0;
    for (i = 1; i <= l; i++) {
	for (j = 1; j <= p; j++) {
	    Pid[i][j] = Pidr[k];
	    k++;
	}
    }
    
    optcut(Pid, mk, maxk);
    
    /* keeps the last row */
    for (k = 1; k <= Km; k++) {
	kk[k] = mk[l][k];
    }
    /* computes the max */
    tmp = kk[1];
    for (k = 1; k <= Km; k++) {
	if (tmp - kk[k] < -0.0000000001)
	    tmp = kk[k];
    }
    
    /* removes the max */
    for (k = 1; k <= Km; k++) {
	kk[k] = kk[k] - tmp;
    }
    
    /* Computes mik */
    mik[1] = exp(kk[1]);
    for (k = 2; k <= Km; k++) {
	mik[k] = mik[k-1] + (exp(kk[k]) / ((double) k ));
    }
    
    /* adds the max again */
    for (k = 1; k <= Km; k++) {
	mik[k] = log(mik[k]) + tmp;
    }
    
    /* Computes msk */
    msk[Km] = exp(kk[Km]) / ((double) (l - Km + 1)) ;
    for (k = Km-1; k >= 1; k--) {
	msk[k] = msk[k+1] + exp(kk[k]) / ((double) (l - k + 1)  );
    }

    /* adds the max again */
    for (k = 1; k <= Km; k++) {
	msk[k] = log(msk[k]) + tmp;
    }
    
    
    /* The bayes factor */
    for (k = 2; k <= (Km - 1); k++) {
	baye[k-2] = msk[k] - mik[k-1];
    }

    /* ... and back to R */
    k = 0;
    for (j = 1; j<=Km; j++) {
	mkr[k] = mk[l][j];
	k++;
    }
    


    /* Free memory */
    freevec(msk);
    freevec(mik);
    freevec(kk);
    freetab(Pid);
    freetab(mk);
}



/* *********************************************************************
 *                                                                     *
 *                   partitioning of a trajectory                      *
 *                                                                     *
 ***********************************************************************/

void partraj(double **Pid, int *maxk, double **Mk, double **Mkd, 
	     double **res)
{
    /* declaration of variables */
    int i, j, k, m, l, D, Km;
    double **Mkk, **cumPid, tmp;
    
    /* Memory allocation */
    l = Pid[0][0]; /* length of the sequence */
    D = Pid[1][0]; /* Number of models */
    Km = *maxk;    /* Partition size */
    tmp = 0;
    m = 0;
    
    taballoc(&Mkk, Km, D); /* Contains mkd for i = k, for all models */
    taballoc(&cumPid, l, D); /* For the probability of the sequences */
    
    
    /* First Compute the prediction of the sequences */
    for (j = 1; j <= D; j++) {
	cumPid[1][j] = Pid[1][j];
    }
    for (i = 2; i <= l; i++) {
	for (j = 1; j <= D; j++) {
	    cumPid[i][j] = cumPid[i-1][j] + Pid[i][j];
	    /* fills res */
	    res[i][j] = cumPid[i][j];
	}
    }
    
    
    /* Computes M1 */    
    for (i = 1; i <= l; i++) {
	Mk[i][1] = cumPid[i][1];
	for (j = 2; j <= D; j++) {
	    if (Mk[i][1] < cumPid[i][j]) {
		Mk[i][1] = cumPid[i][j];
	    }
	}
    }
    
    /* Then, computes Mkk for all k <= Km */
    for (j = 1; j <= D; j++) {
	Mkk[1][j] = Pid[1][j];
    }
    
    for (k = 2; k <= Km; k++) {
	/* Computes Mkk */
	for (j = 1; j <= D; j++) {
	    Mkk[k][j] = Pid[k][j] + Mk[k-1][k-1];
	}
	/* Update Mk */
	Mk[k][k] = Mkk[k][1];
	for (j = 2; j <= D; j++) {
	    if (Mkk[k][j] - Mk[k][k] < 0.000000000001)
		Mk[k][k] = Mkk[k][j];
	}
    }
    
    
    /* Now, compute Mkd */
    for (k = 2; k <= Km; k++) {
	
	/* Deletes for i < k */
	for (i = 1; i < k; i++) {
	    for (j = 1; j <= D; j++) {
		Mkd[i][j] = 0;
	    }
	}
		
	/* first fill the first line of Mkd */
	for (j = 1; j <= D; j++) {
	    Mkd[k][j] = Mkk[k][j];
	}

	/* Computes Mkd */
	for (i = (k+1); i <= l; i++) {
	    for (j = 1; j <= D; j++) {
		tmp = Mk[i-1][k-1];
		if (Mkd[i-1][j] - tmp > 0.000000000000000001) {
		    tmp = Mkd[i-1][j];
		}
		Mkd[i][j] = Pid[i][j] + tmp;
	    }
	    
	    /* Update Mk */
	    Mk[i][k] = Mkd[i][1];
	    for (j = 1; j <= D; j++) {
		if (Mkd[i][j] - Mk[i][k] > 0.000000000000000001) {
		    Mk[i][k] = Mkd[i][j];
		}
	    }
	    
	}
	
	/* Fills res */
	for (i = 1; i <= l; i++) {
	    for (j = 1; j <= D; j++) {
		res[((k-1) * l)+i][j] = Mkd[i][j];
	    }
	}
    }

    /* free memory */
    freetab(Mkk);
    freetab(cumPid);
}



void partrajr(double *Pidr, double *curmar, int *curmodr, int *curlocr, 
	      int *lr, int *Dr, int *Kmr)
{
    /* Variable declaration */
    int l, D, Km, i, j, k, m, n, new, *curloc, *curmod;
    double **Mk, **Mkd, **Pid, **res, tmp, *curma, **grap;
    
    /* Memory allocation */
    l = *lr;
    D = *Dr;
    Km = *Kmr;
    m = 0;
    n = 0;
    tmp = 0;
    new = 0;
    taballoc(&Mk, l, Km);
    taballoc(&Mkd, l, D);
    taballoc(&Pid, l, D);
    taballoc(&res, (l * Km), D);
    taballoc(&grap, l, D);
    vecalloc(&curma, Km);
    vecintalloc(&curmod, Km);
    vecintalloc(&curloc, (Km+1));
    
    
    /* R to C */
    k = 0;
    for (i = 1; i <= l; i++) {
	for (j = 1; j <= D; j++) {
	    Pid[i][j] = log(Pidr[k]);
	    k++;
	}
    }
    
    /* The main algorithm */
    partraj(Pid, Kmr, Mk, Mkd, res);
    
    
    /* Backtracking */
    curloc[1] = l;
    for (k = Km; k>=1; k--) {
	if (k > 1) {
	    
	    m = Km - k + 2;
	    
            /* Stores the graph */
	    for (i = 1; i <= l; i++) {
		for (j = 1; j <= D; j++) {
		    if (res[(l * (k-1)) + i][j] > Mk[i][k-1]) {
			grap[i][j] = 1;
		    } else {
			grap[i][j] = 0;
		    }
		}
	    }
	    
	    /* The best model for the last step */
	    n = (l * (k-1)) + curloc[m-1];
	    curma[m-1] = res[n][1];
	    curmod[m-1] = 1;
	    for (j = 2; j <= D; j++) {
		if (res[n][j] > curma[m-1]) {
		    curma[m-1] = res[n][j];
		    curmod[m-1] = j;
		}
	    }
	    
	    /* Keep this model until ? */
	    n = curloc[m-1];
	    j = curmod[m-1];
	    while (grap[n][j] > 0.0000000001)
		n--;
	    curloc[m] = n;
	} else {
	    m = Km+1;
	    curloc[m] = 1;
	    curma[m-1] = res[curloc[m-1]][1];
	    curmod[m-1] = 1;
	    for (j = 1; j <= D; j++) {
		if (res[curloc[m-1]][j] > curma[m-1]) {
		    curma[m-1] = res[curloc[m-1]][j];
		    curmod[m-1] = j;
		}
	    }
	}
	
	
    }
    
    
    
    /* C to R */
    for (i = 1; i <= Km; i++) {
	curmodr[i-1] = curmod[i];
	curlocr[i-1] = curloc[i];
	curmar[i-1] = curma[i];
    }
    curlocr[Km] = curloc[Km+1];
    
    
    /* free memory */
    freetab(Mk);
    freetab(Mkd);
    freetab(res);
    freetab(grap);
    freevec(curma);
    freeintvec(curmod);
    freeintvec(curloc);
}






/* acfdist and acfang */


void acfdist (double *sim, double *di, int *ndi, double *diper, int *nbobs, int *indxnona, int *nsim, int *lag){
  
  int i,j, *nona, *numero;
  double nobslag;
  vecintalloc(&nona, *nbobs);
  vecintalloc(&numero, *nbobs);
  for (i=1; i<=*nbobs; i++) {
    nona[i] = indxnona[i-1]-1;
  }
  nobslag=0;
  /* compute observed statistic */
  for (i=1;i<=(*nbobs);i++){
    if((nona[i]+*lag)<=(*ndi-1)){
      if(!ISNAN(di[nona[i] + *lag])) {
	sim[0]=sim[0]+pow(di[nona[i]+*lag] - di[nona[i]],2);
	nobslag=nobslag+1; /* number of observed differences */

      }
    }
  }
 
  
  /* permute only observed values */
  for (j=1; j<=*nsim; j++) {
    getpermutation(numero,j);
    for (i=1;i<=(*nbobs);i++){
      diper[nona[i]]=di[nona[numero[i]]];
    }
    for (i=1;i<=(*nbobs);i++){
      if((nona[i]+*lag)<=(*ndi-1)){
	if(!ISNAN(di[nona[i] + *lag])) {
	  sim[j]=sim[j]+pow(diper[nona[i]+*lag] - diper[nona[i]],2);
	}
      }
    }
  }

  for (j=0; j<=*nsim; j++) {
    sim[j]=sim[j]/nobslag; /* rescale the stat by dividing by the number of observed differences */
  }
  freeintvec(nona);
  freeintvec(numero);
}

void acfangl (double *sim, double *ang, int *nang, double *angper, int *nbobs, int *indxnona, int *nsim, int *lag){
  
  int i,j, *nona, *numero;
  double nobslag;  
  vecintalloc(&nona, *nbobs);
  vecintalloc(&numero, *nbobs);
  for (i=1; i<=*nbobs; i++) {
    nona[i] = indxnona[i-1]-1;
  }
  nobslag=0;
  /* compute observed statistic */
  for (i=1;i<=(*nbobs);i++){
    if((nona[i]+*lag)<=(*nang-1)){
      if(!ISNAN(ang[nona[i] + *lag])) {
	sim[0]=sim[0]+2.0*(1.0-cos(ang[nona[i]+*lag] - ang[nona[i]]));
	nobslag=nobslag+1; /* number of observed differences */
      }
    }
  }
  /* permute only observed values */
  for (j=1; j<=*nsim; j++) {
    getpermutation(numero,j);
    for (i=1;i<=(*nbobs);i++){
      angper[nona[i]]=ang[nona[numero[i]]];
    }
    for (i=1;i<=(*nbobs);i++){
      if((nona[i]+*lag)<=(*nang-1)){
	if(!ISNAN(ang[nona[i] + *lag])) {
	  sim[j]=sim[j]+2.0*(1.0-cos(angper[nona[i]+*lag] - angper[nona[i]]));
	}
      }
    }
  }
  for (j=0; j<=*nsim; j++) {
    sim[j]=sim[j]/nobslag; /* rescale the stat by dividing by the number of observed differences */
  }
  freeintvec(nona);
  freeintvec(numero);

}






/* *********************************************************************
 *                                                                     *
 *                       Rasterizing a trajectory                      *
 *                                                                     *
 ***********************************************************************/




SEXP RasterPas(SEXP df, SEXP xllr, SEXP yllr, SEXP cs, SEXP type1)
{
    int npas, i, j, nso, k;
    SEXP xl, yl, resu, so, xso, yso, dfso;
    double x1, y1, x2, y2, 
	xc, yc, dist, xt, yt, csi, xll, yll;
    
    npas = length(VECTOR_ELT(df,0)) - 1;
    csi = REAL(cs)[0];
    nso = INTEGER(type1)[0];
    xll = REAL(xllr)[0];
    yll = REAL(yllr)[0];
    
    PROTECT(xl = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(yl = coerceVector(VECTOR_ELT(df,1), REALSXP));
    
    
    if (nso > 0) {
	PROTECT(xso = allocVector(REALSXP, nso));
	PROTECT(yso = allocVector(REALSXP, nso));
	PROTECT(so = allocVector(REALSXP, nso));
    }
    
    nso = 0;
    
    for (i = 0; i < npas; i++) {
	
	/* premier filtre */
	x1 = REAL(xl)[i];
	x2 = REAL(xl)[i+1];
	y1 = REAL(yl)[i];
	y2 = REAL(yl)[i+1];
	
	/* identification des coordonnées initiales et finales */	
	dist = pythag(x2-x1, y2-y1);
	k = (int) round(50.0*dist/csi);
	if (k == 0)
	    k++;
	xc = -230876.0;
	yc = -230876.0;
	
	for (j = 0; j<k; j++) {
	    yt = y1 + ((double) j)*(y2-y1)/((double) k);
	    xt = x1 + ((double) j)*(x2-x1)/((double) k);
	    xt = xll + round((xt - xll)/csi)*csi;
	    yt = yll + round((yt - yll)/csi)*csi;
	    if (pythag(yt - yc, xt - xc) > 0.00000000001) {
		xc = xt;
		yc = yt;
		if (INTEGER(type1)[0] != 0) {
		    REAL(xso)[nso] = xc;
		    REAL(yso)[nso] = yc;
		    REAL(so)[nso] = ((double) (i+1));
		} 
		nso++;
	    }
	}
    }
    
    if (INTEGER(type1)[0] != 0) {
	PROTECT(dfso = allocVector(VECSXP, 3));
    	SET_VECTOR_ELT(dfso, 0, xso);
	SET_VECTOR_ELT(dfso, 1, yso);
	SET_VECTOR_ELT(dfso, 2, so);
	UNPROTECT(6);
	
	return(dfso);
    } else {
	PROTECT(resu = allocVector(INTSXP, 1));
	INTEGER(resu)[0] = nso;
	PROTECT(dfso = RasterPas(df, xllr, yllr, cs, resu));
	UNPROTECT(4);
	return(dfso);
    }

}
