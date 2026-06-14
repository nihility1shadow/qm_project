#ifndef _YY_CPLUSPLUS_H_
#define _YY_CPLUSPLUS_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <cmath>

#if 0
#ifndef _cplusplus
  #include <complex.h>
	typedef _Complex double         dcomplex;
	typedef _Complex float          fcomplex;
	typedef _Complex long double    lcomplex;
#else
  #include <complex>
  typedef complex<float>       fcomplex;
  typedef complex<double>      dcomplex;
  typedef complex<long double> lcomplex;
#endif
#endif

#include <yyy_complex.h>

template<class T> inline T   MIN(const T &a, const T &b) { return (a) > (b) ? (b) : (a);}
template<class T> inline T   MAX(const T &a, const T &b) { return (a) < (b) ? (b) : (a);}
template<class T> inline T   SQR(const T &a) {return (a*a);}
template<class T> inline T   TRI(const T &a) {return (a*a*a);}
template<class T> inline T   QUD(const T &a) {T b = a*a; return (b*b);}
template<class T> inline int SIGN(const T &b) { if(b==0) return 0; return  (b >0 ? 1 : -1); }

template<class T> inline 
int __abscompare_(const void *a, const void *b) {
  T *x = (T *) a,
    *y = (T *) b;

  return (abs(*x) > abs(*y) ? -1 : 1);
}

template<class T> inline void SWAP(T *a, T *b) {
  T tmp;

  tmp = *a;
  *a  = *b;
  *b  = tmp;

  return;
}

template<class T> inline void free1d(T *&a) {
  if(a) {
    free(a);
    a = 0;
  }
  return;
}

template<class T> inline void free2d(T **&a) {
  if(a) {
    free(a[0]);
    a[0] = 0;
    free(a);
    a = 0;
  }

  return;
}

template<class T> inline void free3d(T ***&a) {
  if(a) {
    free2d<T>(*a);
    free(a);
    a = 0;
  }

  return;
}

template<class T> inline void free4d(T ****&a) {
  if(a) {
    free3d<T>(*a);
    free(a);
    a = 0;
  }

  return;
}

template<class T> inline T *array1d(const size_t m) {
  T *a;

#ifdef __ENABLE_CHECK_
  if(m==0) {
    fprintf(stderr, "Warning: array1d --- size of the array is zero.\n");
    abort();
  }
#endif

  a  = (T *)calloc(m, sizeof(T));
  if(a == NULL) {
    perror("failed to calloc memory.\n");
    abort();
  }
  bzero(a, sizeof(T)*m);

  return a;
}

template<class T> inline T *realloc1d(T *arr, const size_t m, const size_t mold) {
  if(m == mold) return arr;

  T *ptr = (T *)malloc(ptr, sizeof(T)*m);
  if(!ptr) {
    fprintf(stderr, "Error: realloc1d --- realloc for 1d array fails.\n");
    return (NULL);
  }

  size_t min = m < mold ? m : mold;
  memcpy(ptr, arr, sizeof(T)*min);
  free(arr);

  return ptr;
}

template<class T> inline T **array2d(const size_t n, const size_t m) {
  size_t i;
  T **a;

#ifdef __ENABLE_CHECK_
  if(m==0 || n==0) {
    fprintf(stderr, "Warning: array2d --- size of the array is zero.\n");
    abort ();
  }
#endif

  a  = (T **)calloc(n, sizeof(T *));
  *a = (T *)calloc(n*m, sizeof(T));
  if(*a == NULL || a == NULL) {
    fprintf(stderr, "Error: array2d --- failed to allocate memory.\n");
    abort ();
  }
  for(i=1; i<n; i++)  *(a+i) = *(a+i-1) + m;

  i = sizeof(T)*m*n;
  bzero(*a, i);

  return a;
}

template<class T> inline T **realloc2d(T **arr, const size_t n, const size_t m, 
    const size_t nold, const size_t mold) {
  if(n == nold && m == mold) return arr;

  T **ptr;
  size_t min = n*m;

#ifdef  __ENABLE_CHECK_
  if(m==0 || n==0) {
    fprintf(stderr, "Warning: realloc2d --- size of the array is zero.\n");
    abort ();
  }
#endif

  ptr =  array2d<T>(n, m);
  if(mold*nold<min) min = mold*nold;

  memcpy(*ptr, *arr, sizeof(T)*min);
  free2d(arr);
  return ptr;
}

template<class T> inline T **IRarray2d(const size_t n, const size_t *m){
    size_t i, ntot=0;
    T **a;

    for(i=0; i<n; i++) ntot += m[i];

    a = (T **)malloc(sizeof(T *)*n);
    if(a)
        *a = (T *) malloc(sizeof(T)*ntot);
    else{
        printf("malloca failed.\n");
        exit (errno);
    }

    if(*a == NULL){
           printf("malloca failed.\n");
         exit (errno);
    }

    for(i=1; i<n; i++)  *(a+i) = *(a+i-1) + m[i-1];

    bzero(*a, sizeof(T)*ntot);
    return a;
}

template <class T> inline T ***array3d(const size_t m, const size_t n, const size_t l){
    size_t i, j;
    T ***a;

    a = (T ***)malloc(sizeof(T **)*m);
   *a = (T  **)malloc(sizeof(T  *)*m*n);
  **a = (T   *)malloc(sizeof(T   )*m*n*l);

    for(i=1; i<m; i++) a[i] = a[i-1] + n;
    for(i=0; i<m; i++) {
      if(i) a[i][0] = a[i-1][0] + n*l;
      for(j=1; j<n; j++) a[i][j] = a[i][j-1] + l;
    }

    bzero(**a, sizeof(T)*n*m*l);
    return a;
}

template <class T> inline T ****array4d(const size_t m, const size_t n, const size_t l, const size_t k){
    size_t i, j, q;
    T ****a;

     a = (T ****)malloc(sizeof(T ***)*m);
    *a = (T  ***)malloc(sizeof(T  **)*m*n);
   **a = (T   **)malloc(sizeof(T   *)*m*n*l);
  ***a = (T    *)malloc(sizeof(T    )*m*n*l*k);

    for(i=1; i<m; i++) a[i] = a[i-1] + n;
    for(i=0; i<m; i++) {
      if(i) a[i][0] = a[i-1][0] + n*l;
      for(j=1; j<n; j++) a[i][j] = a[i][j-1] + l;
      for(j=0; j<n; j++) for(q=0; q<l; q++) {
        a[i][j][q] = a[0][0][0] + i*n*l*k + j*l*k + q*k;
      }
    }

    bzero(***a, sizeof(T)*n*m*l*k);
    return a;
}

inline void strprintf(FILE *FL, const char *str, const int n) {
 for(int j=0; j<n; j++) {
   fprintf(FL, "%2d", static_cast<int>(str[j]));
  }
  fprintf(FL, "\n");

  return;
}

/*
 * transpose a square matrix
 */
template<class T> inline
void transpose_sqr_matrix(const size_t m, T **mtx) {
  T f;
  for(int j=0; j<m; j++) for(int k=0; k<j; k++) {
    f = mtx[j][k];
    mtx[j][k] = mtx[k][j];
    mtx[k][j] = f;
  }

  return;
}

template<class T> inline
void matrix_print(FILE *fp, const size_t row, const size_t col, const T *arr,
		const char *pfmt, const char *delim, 
		const char *begining,     const char *ending, 
		const char *begining_row, const char *ending_row){
	size_t i, j, len;
	FILE *FL;
	char *fmt;
	len = strlen(pfmt);
	fmt = new char [len+3];
  strcpy(fmt, pfmt);
  strcat(fmt, "%s");

	if(fp == NULL) FL = stdout;
	else           FL = fp;

	if(begining) fprintf(FL, "%s", begining);
	for(i=0; i<row; i++) {
		if(begining_row) fprintf(FL, "%s", begining_row);
		for(j=0; j<col-1; j++) {
			if(delim) fprintf(FL, fmt, arr[i*col+j], delim);
			else      fprintf(FL,pfmt, arr[i*col+j]);
		}
		if(ending_row) fprintf(FL,  fmt, arr[i*col+col-1], ending_row);
		else           fprintf(FL, pfmt, arr[i*col+col-1]);
		if(i==row-1) {
				if(ending) fprintf(FL, "%s\n", ending);
				else fprintf(FL, "\n");
		}
		else {
			if(delim) fprintf(FL, "%s\n", delim);
			else      fprintf(FL, "\n");
		}
	}

	delete [] fmt;
	return;
}

#endif
