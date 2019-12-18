#include <iostream> //cout
#include <cmath> //log, exp
#include <algorithm> //min
#include "func.h" //logsum, sum, first_geq, ...

using namespace std;

void DP(double* output, int* I1_l, unsigned int* I2_start,
        int* I2_l, unsigned int* I2_end, int n, int m, unsigned int g,
        int mode=0)
{   
    if (I2_end[m]>g)
    {
        cout << "Wrong input! Intervals are out of range" << endl;
        return;
    }
    int L = 0;
    for (int i=1; i<=n; i++)
    {
        L = max(L, I1_l[i]);
    }
    L++;
    double N = -1;
    double **** N1 = new double *** [L];
    for (int hh=0; hh<L; hh++)
    {
        N1[hh] = new double ** [n+1];   
        for (int ii=0; ii<=n; ii++)
        {
            N1[hh][ii] = new double * [m+1];
            for (int kk=0; kk<=m; kk++)
            {
                N1[hh][ii][kk] = new double [2];
                for (int aa=0; aa<=1; aa++)
                {
                    N1[hh][ii][kk][aa] = -1;
                }
            }
        }
    }
    int cval = 0;
    int diff = 0;
    int f1 = 0;
    int f2 = 0;
    int* Lsum = new int [n+1];
    Lsum[1] = I1_l[1];
    for (int i=2; i<=n; i++)
    {
       Lsum[i] = Lsum[i-1] + I1_l[i]; 
    }
    for (int h=1; h<=g; h++)
    {
        for (int i=1; i<=n; i++)
        {
            cval =  c(I1_l[i], h, m, I2_start, I2_end);
            f1 = f(h-I1_l[i], m, I2_start, I2_end);
            f2 = f(h-1, m, I2_start, I2_end);
            for (int k=0; k<=m; k++)
            {
                for (int a=0; a<=1; a++)
                {
                    diff = cval-a;
                    if ((i==1)&&(h>=I1_l[1])&&(k==diff))
                    {
                        N = 0;
                    }
                    if ((i>1)&&(h>=Lsum[i])&&(k-diff>=0)&&(diff>=0))
                    {

                        N = N1[(h-I1_l[i])%L][i-1][k-diff][f1];
                    }
                    if (h==1)
                    {
                        N1[h%L][i][k][a] = N;
                    }
                    else
                    {
                        N1[h%L][i][k][a] = logsum(N, N1[(h-1)%L][i][k][min(a,f2)]);
                    }
                    N = -1;
                }
            }
        }
    }
    double Sum1 = -1;
    for (int k=0; k<=m; k++)
    {
        Sum1 = logsum(Sum1, N1[g%L][n][k][0]);
    }
    int nn = g-Lsum[n]+n;
    double ratio = lognchoosek(nn, n) - Sum1;
    // cout << exp(ratio) << endl;

    if (mode==1)
    {
        for (int k=0; k<=m; k++)
        {
            output[k] = N1[g%L][n][k][0];
        }
    }
    else
    {
        double* CumSum = new double [m+1];
        for (int k=0; k<=m; k++)
        {
            if (k==0)
            {
                CumSum[m] = N1[g%L][n][m][0];
            }
            else 
                CumSum[m-k] = logsum(CumSum[m-k+1], N1[g%L][n][m-k][0]);
        }
        for (int k=0; k<=m; k++)
        {
            if (CumSum[k]==-1)
                output[k] = 0;
            else 
                output[k] = exp(CumSum[k] - Sum1);
        }
        delete [] CumSum;
    }
    

    for (int hh=0; hh<L; hh++)
    {  
        for (int ii=0; ii<=n; ii++)
        {
            for (int kk=0; kk<=m; kk++)
            {
                delete [] N1[hh][ii][kk];
            }
            delete [] N1[hh][ii];
        }
        delete [] N1[hh];
    }
    delete [] N1;
    delete [] Lsum;
}