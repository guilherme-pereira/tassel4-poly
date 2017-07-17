//
//  net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix.c
//  
//
//  Created by Peter Bradbury on 7/5/13.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vecLib/cblas.h>
#include <limits.h>
#include "net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix.h"

extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double A[], int *lda, double S[], double U[],
                    int *ldu, double VT[], int *ldvt, double work[], int *lwork, int *info);

extern void dsyevr_(char *jobz, char *range, char *uplo, int *n, double A[], int *lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double w[], double z[], int *ldz, int isuppz[], double work[], int *lwork, int iwork[], int *liwork, int *info);

extern void dgelsd_(int *M, int *N, int *NRHS, double A[], int *LDA, double B[], int *LDB, double S[], double *RCOND, int *RANK, double WORK[], int *LWORK, int IWORK[], int *INFO);

extern void dgelsy_(int *M, int *N, int *NRHS, double A[], int *LDA, double B[], int *LDB, int JPVT[], double *RCOND, int *RANK, double WORK[], int *LWORK, int *INFO);

extern void dgesdd_(char *JOBZ, int *M, int *N, double A[], int *LDA, double S[], double U[], int *LDU, double *VT, int *LDVT, double WORK[], int *LWORK, int IWORK[], int *INFO);

JNIEXPORT void JNICALL Java_net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix_multMatrices
(JNIEnv *env, jclass myclass, jdoubleArray A, jint nrowsA, jint ncolsA, jdoubleArray B, jint nrowsB, jint ncolsB, jdoubleArray C, jdouble alpha, jdouble beta, jboolean transA, jboolean transB) {
    
    int M, N, K, LDA, LDB, LDC;
    double alph, bet;
    jboolean acopy = 0, bcopy = 0, ccopy = 0;
    
    jdouble *aptr = (*env)->GetDoubleArrayElements(env, A, &acopy);
    jdouble *bptr = (*env)->GetDoubleArrayElements(env, B, &bcopy);
    jdouble *cptr = (*env)->GetDoubleArrayElements(env, C, &ccopy);
    enum CBLAS_TRANSPOSE trA, trB;
    
    
    if (transA) {
        trA = CblasTrans;
        M = ncolsA;
        K = nrowsA;
    }
    else {
        trA = CblasNoTrans;
        M = nrowsA;
        K = ncolsA;
    }
    
    if (transB) {
        trB = CblasTrans;
        N = nrowsB;
    } else {
        trB = CblasNoTrans;
        N = ncolsB;
    }
    
    LDA = nrowsA;
    LDB = nrowsB;
    LDC = M;
    alph = alpha;
    bet = beta;
    
    cblas_dgemm(CblasColMajor, trA, trB, M, N, K, alph, aptr, LDA, bptr, LDB, bet, cptr, LDC);
    
    
    (*env)->ReleaseDoubleArrayElements(env, A, aptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, B, bptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, C, cptr, 0);
    return;
}

JNIEXPORT jint JNICALL Java_net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix_solveLSdgelsd
(JNIEnv *env, jclass myclass, jdoubleArray A, jint Arows, jint Acols, jdoubleArray B, jint Bcols, jdouble rcond, jintArray rank) {
    int m = Arows;
    int n = Acols;
    int nrhs = Bcols;
    jboolean acopy = 0, bcopy = 0, rcopy = 0;
    jdouble *aptr = (*env)->GetDoubleArrayElements(env, A, &acopy);
    int lda = m;
    jdouble *bptr = (*env)->GetDoubleArrayElements(env, B, &bcopy);
    int ldb = m;
    int sizeS = m*n;
    double *S;
    S = malloc(sizeS * sizeof(double));
    int RANK = 0;
    double worksize[1];
    int lwork = -1;
    int smlsiz = 25;
    int minmn = m > n ? n : m;
    int nlvl = ((int)log2(minmn/(smlsiz + 1))) + 1;
    if (nlvl < 0) nlvl = 0;
    int liwork = 3 * minmn * nlvl + 11 * minmn + 1000;
    
    int *iwork;
    iwork = malloc(liwork * sizeof(int));
    int info = 0;
    
    dgelsd_(&m, &n, &nrhs, aptr, &lda, bptr, &ldb, S, &rcond, &RANK, worksize, &lwork, iwork, &info);
    
    lwork = worksize[0];
    double *work;
    work = malloc(lwork * sizeof(double));
    
    dgelsd_(&m, &n, &nrhs, aptr, &lda, bptr, &ldb, S, &rcond, &RANK, work, &lwork, iwork, &info);
    
    free(S);
    free(work);
    free(iwork);
    
    (*env)->ReleaseDoubleArrayElements(env, A, aptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, B, bptr, 0);
    jint *rptr = (*env)->GetIntArrayElements(env, rank, &rcopy);
    rptr[0] = RANK;
    (*env)->ReleaseIntArrayElements(env, rank, rptr, 0);
    return (jint) info;
}

JNIEXPORT jint JNICALL Java_net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix_solveLSdgelsy
(JNIEnv *env, jclass myclass, jdoubleArray A, jint Arows, jint Acols, jdoubleArray B, jint Bcols, jdouble rcond, jintArray rank) {
    int m = Arows;
    int n = Acols;
    int nrhs = Bcols;
    jboolean acopy = 0, bcopy = 0, rcopy = 0;
    jdouble *aptr = (*env)->GetDoubleArrayElements(env, A, &acopy);
    int lda = m;
    jdouble *bptr = (*env)->GetDoubleArrayElements(env, B, &bcopy);
    int ldb = m;
    int RANK = 0;
    double worksize[1];
    int lwork = -1;
    int info = 0;
    int *JPVT;
    JPVT = malloc(n * sizeof(int));
    memset(JPVT, 0, n);
    
    dgelsy_(&m, &n, &nrhs, aptr, &lda, bptr, &ldb, JPVT, &rcond, &RANK, worksize, &lwork, &info);
    
    lwork = worksize[0];
    double *work;
    work = malloc(lwork * sizeof(double));
    
    dgelsy_(&m, &n, &nrhs, aptr, &lda, bptr, &ldb, JPVT, &rcond, &RANK, work, &lwork, &info);
    
    free(JPVT);
    free(work);
    
    (*env)->ReleaseDoubleArrayElements(env, A, aptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, B, bptr, 0);
    jint *rptr = (*env)->GetIntArrayElements(env, rank, &rcopy);
    rptr[0] = RANK;
    (*env)->ReleaseIntArrayElements(env, rank, rptr, 0);
    return (jint) info;
    
}

JNIEXPORT jint JNICALL Java_net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix_singularValueDecompositionDgesdd
(JNIEnv *env, jclass myclass, jchar jobz, jint m, jint n, jdoubleArray A, jint lda, jdoubleArray S, jdoubleArray U, jint ldu, jdoubleArray VT, jint ldvt) {
    char JOBZ;
    int M, N, LDA, LDU, LDVT;
    jboolean acopy = 0, scopy = 0, ucopy = 0, vtcopy = 0;
    
    JOBZ = jobz;
    M = m;
    N = n;
    LDA = lda;
    LDU = ldu;
    LDVT = ldvt;
    
    jdouble *aptr = (*env)->GetDoubleArrayElements(env, A, &acopy);
    jdouble *sptr = (*env)->GetDoubleArrayElements(env, S, &scopy);
    jdouble *uptr = (*env)->GetDoubleArrayElements(env, U, &ucopy);
    jdouble *vtptr = (*env)->GetDoubleArrayElements(env, VT, &vtcopy);
    
    double worksize[1];
    int lwork = -1;
    int info = 0;
    int liwork;
    if (M < N) liwork = 8*M;
    else liwork = 8*N;
    int* iwork;
    iwork = malloc(liwork * sizeof(int));
    
    dgesdd_(&JOBZ, &M, &N, aptr, &LDA, sptr, uptr, &LDU, vtptr, &LDVT, worksize, &lwork, iwork, &info);
    
    lwork = worksize[0];
    double* work;
    work = malloc(lwork * sizeof(double));
    
    dgesdd_(&JOBZ, &M, &N, aptr, &LDA, sptr, uptr, &LDU, vtptr, &LDVT, work, &lwork, iwork, &info);
    
    free(iwork);
    free(work);
    
    (*env)->ReleaseDoubleArrayElements(env, A, aptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, S, sptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, U, uptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, VT, vtptr, 0);
    
    return (jint) info;
}

JNIEXPORT jint JNICALL Java_net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix_singularValueDecompositionDgesvd
(JNIEnv *env, jclass myclass, jchar jobu, jchar jobvt, jint m, jint n, jdoubleArray A, jint lda, jdoubleArray S, jdoubleArray U, jint ldu, jdoubleArray VT, jint ldvt){
    char JOBU, JOBVT;
    int M, N, LDA, LDU, LDVT;
    jboolean acopy = 0, scopy = 0, ucopy = 0, vtcopy = 0;
    
    M = m;
    N = n;
    LDA = lda;
    LDU = ldu;
    JOBU = jobu;
    JOBVT = jobvt;

    jdouble *aptr = (*env)->GetDoubleArrayElements(env, A, &acopy);
    jdouble *sptr = (*env)->GetDoubleArrayElements(env, S, &scopy);
    jdouble *uptr = (*env)->GetDoubleArrayElements(env, U, &ucopy);
    jdouble *vtptr = (*env)->GetDoubleArrayElements(env, VT, &vtcopy);
    
    double worksize[1];
    int lwork = -1;
    int info = 0;
    
    dgesvd_(&JOBU, &JOBVT, &M, &N, aptr, &LDA, sptr, uptr, &LDU, vtptr, &LDVT, worksize, &lwork, &info);
    
    lwork = worksize[0];
    double *work;
    work = malloc(lwork * sizeof(double));
    dgesvd_(&JOBU, &JOBVT, &M, &N, aptr, &LDA, sptr, uptr, &LDU, vtptr, &LDVT, work, &lwork, &info);
    
    free(work);
    
    (*env)->ReleaseDoubleArrayElements(env, A, aptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, S, sptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, U, uptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, VT, vtptr, 0);
    
    return (jint) info;
}

JNIEXPORT jint JNICALL Java_net_maizegenetics_matrixalgebra_Matrix_BlasDoubleMatrix_eigenValueSymmetricDecomposition
(JNIEnv *env, jclass myclass, jint order, jdoubleArray A, jdoubleArray eigenvalues, jdoubleArray eigenvectors){
    char jobz = 'V';
    char range = 'A';
    char uplo = 'U';
    int N = order;
    jboolean acopy = 0;
    double *Aptr = (*env)->GetDoubleArrayElements(env, A, &acopy);
    int lda = N;
    double vl = 0;
    double vu = 0;
    int il = 0;
    int iu = 0;
    double abstol = 0;
    int M = 0;
    jboolean wcopy = 0;
    double *Wptr = (*env)->GetDoubleArrayElements(env, eigenvalues, &wcopy);
    jboolean zcopy = 0;
    double *Zptr = (*env)->GetDoubleArrayElements(env, eigenvectors, &zcopy);
    int ldz = N;
    int sizeisuppz = 2 * N;
    int *isuppz;
    isuppz = malloc(sizeisuppz * sizeof(int));
    double worksize[1];
    int lwork = -1;
    int iworksize[1];
    int liwork = -1;
    int info = 0;
    
    dsyevr_(&jobz, &range, &uplo, &N, Aptr, &lda, &vl, &vu, &il, &iu, &abstol, &M, Wptr, Zptr, &ldz, isuppz, worksize, &lwork, iworksize, &liwork, &info);
    
    lwork = worksize[0];
    liwork = iworksize[0];
    double *work;
    work = malloc(lwork * sizeof(double));
    int *iwork;
    iwork = malloc(liwork * sizeof(int));
    
    dsyevr_(&jobz, &range, &uplo, &N, Aptr, &lda, &vl, &vu, &il, &iu, &abstol, &M, Wptr, Zptr, &ldz, isuppz,
            work, &lwork, iwork, &liwork, &info);
    
    free(isuppz);
    free(iwork);
    free(work);
    
    (*env)->ReleaseDoubleArrayElements(env, A, Aptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, eigenvalues, Wptr, 0);
    (*env)->ReleaseDoubleArrayElements(env, eigenvectors, Zptr, 0);
    
    return (jint) info;
}
