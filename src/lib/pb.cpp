#include <iostream>
#include <cmath> //ceil
#include <algorithm> //min

using namespace std;

void poisbinom(double* p_value_poisbinom, int* I1_l, int* I2_l, int n, int m, unsigned int g)
{
    double** q = new double* [n+1];
    for(int i=0; i<n+1; i++)
        q[i] = new double [m+1];
    double prod;
    double* p = new double [m+1];
    for(int j=1; j<=m; j++){
        int z = ceil(double(I2_l[j])/2);
        z = 1;
        prod = 1;
        for(int i=1; i<=n; i++){
            q[i][j] = 0;
            if (z<=min(I1_l[i], I2_l[j]))
                q[i][j] = I1_l[i] + I2_l[j] - 2*z + 1;
            q[i][j] = q[i][j]/g;
            prod *= (1-q[i][j]);
        }
        p[j] = 1-prod;
    }

    double** x = new double* [m+1];
    for(int k=0; k<=m; k++)
       x[k] = new double [m+1];
    for(int k=0; k<=m; k++)
    {
       for(int j=k; j<=m; j++)
       {
            if(k==0 && j==0)
                x[k][j] = 1;
            else if(k==0)
                x[k][j] = x[k][j-1] * (1-p[j]);
            else if(k==j)
                x[k][j] = x[k-1][j-1] * (p[j]);
            else
            {
                x[k][j] = x[k][j-1] * (1-p[j]) + x[k-1][j-1] * (p[j]);
            }

       }
    }
    double* CumSum = new double [m+1];
    for(int k=0; k<=m; k++)
    {
        if (k==0)
        {
            CumSum[m] = x[m][m];
        }
        else CumSum[m-k] = CumSum[m-k+1] + x[m-k][m];
    }
    for(int k=0; k<=m; k++)
    {
        p_value_poisbinom[k] = CumSum[k];
    }

    for(int k=0; k<=m; k++)
    {
        delete [] x[k];
    }
    for(int i=0; i<n+1; i++)
        delete [] q[i];
    delete [] q;
    delete [] p;
    delete [] x;
    delete [] CumSum;
}

void poisbinom(double* p_value_poisbinom, double* p, int m)
{
    double** x = new double* [m+1];
    for(int k=0; k<=m; k++)
       x[k] = new double [m+1];
    for(int k=0; k<=m; k++)
    {
       for(int j=k; j<=m; j++)
       {
            if(k==0 && j==0)
                x[k][j] = 1;
            else if(k==0)
                x[k][j] = x[k][j-1] * (1-p[j]);
            else if(k==j)
                x[k][j] = x[k-1][j-1] * (p[j]);
            else
            {
                x[k][j] = x[k][j-1] * (1-p[j]) + x[k-1][j-1] * (p[j]);
            }

       }
    }
    double* CumSum = new double [m+1];
    for(int k=0; k<=m; k++)
    {
        if (k==0)
        {
            CumSum[m] = x[m][m];
        }
        else CumSum[m-k] = CumSum[m-k+1] + x[m-k][m];
    }
    for(int k=0; k<=m; k++)
    {
        p_value_poisbinom[k] = CumSum[k];
    }

    for(int k=0; k<=m; k++)
    {
        delete [] x[k];
    }
    delete [] x;
    delete [] CumSum;
}