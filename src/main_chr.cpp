#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>       // stoi, stoul     
#include <sstream>
#include <math.h>
#include <algorithm>
#include "dp.h"
#include "pb.h"
#include "prep.h"
#include "func.h"
#include <thread>

using namespace std;

static const char *const HEADER = "\nISTAT v1.0.0 © The Regents of the University of California. All Rights Reserved.\n\n";
static const char *const USAGE = "Usage:\n\tistat <QUERY_INTERVALS_FILE> <REFERENCE_INTERVALS_FILE> <CHROMOSOME_SIZES_FILE>\n\t<SCALING_FACTOR> <OUTPUT_FILE> <THREADING_MODE>\n\nDescription:\n\tComputes the p-value of overlap between query and reference intervals.\n";

void chrom(int* n_chr, int* m_chr, double* NM_chr, int* NI, double* eta_chr, double* Pj_chr, 
             string query_path, string reference_path, string chrom_size_path, int chr_id, string id, 
             int scale)
{
    int n;
    int m;
    unsigned int g = 0;

    n = scan_input(query_path, id);

    unsigned int* I1_start_0 = new unsigned int [n+1];
    unsigned int* I1_end_0 = new unsigned int [n+1];

    read_input(I1_start_0, I1_end_0, query_path, id);

    m = scan_input(reference_path, id);

    unsigned int* I2_start_0 = new unsigned int [m+1];
    unsigned int* I2_end_0 = new unsigned int [m+1];
    
    read_input(I2_start_0, I2_end_0, reference_path, id);

    ifstream yourfile (chrom_size_path);
    if (yourfile.is_open())
    {
        string line;
        string token;
        while(getline(yourfile, line))
        {
            stringstream iss;
            iss << line;
            getline(iss, token, '\t');
            if(token==id)
            {
                getline(iss, token);
                g = stoul(token);
                g = ceil(double(g)/scale);
            }
        }
        yourfile.close();
    }
    else cout << "Unable to open file" << endl;

    n = merge_adj(I1_start_0, I1_end_0, n);
    m = merge_adj(I2_start_0, I2_end_0, m);

    // for (int i=0; i<=m; i++)
    //     cout << I2_start_0[i] << "\t" << I2_end_0[i] << endl;

    int* I1_l = new int [n+1];
    int* I1_bar_l = new int [n+1];
    int* I2_l = new int [m+1];
    int* I2_bar_l = new int [m+1];

    compute_length(I1_l, I1_bar_l, I1_start_0, I1_end_0, n);
    compute_length(I2_l, I2_bar_l, I2_start_0, I2_end_0, m);

    unsigned int* I1_start = new unsigned int [n+1];
    unsigned int* I1_end = new unsigned int [n+1];
    unsigned int* I2_start = new unsigned int [m+1];
    unsigned int* I2_end = new unsigned int [m+1];

    scale_interval(I1_start, I1_end, I1_l, I1_bar_l, n, scale);
    scale_interval(I2_start, I2_end, I2_l, I2_bar_l, m, scale);

    // cout << g << " ";
    g = max(g, max(I1_end[n], I2_end[m])); //preventing ceiling artifact
    n_chr[chr_id] = n;
    m_chr[chr_id] = m;
    eta_chr[chr_id] = 1;

    if (m>1)
    {
        int mid = (m-1)/2;
        nth_element(I2_bar_l+2, I2_bar_l+2+mid, I2_bar_l+m+1);
        double median = I2_bar_l[mid+2];
        eta_chr[chr_id] = n*median/g;
    }
    // cout << n << " " << m << " " << g << " ";

    if(n==0)
    {
        NI[chr_id] = 0;
        eta_chr[chr_id] = 1;
        NM_chr[0] = 0;
    }
    else if(m==0)
    {
        NI[chr_id] = 0;
        //NM_chr[0] = 0; //or lognchoosek(g-sumI1l+n,n)
        int nn = g-sum(I1_l, n)+n;
        NM_chr[0] = lognchoosek(nn, n);
    }

    if (n==0 || m==0)
    {
        delete [] I1_start_0;
        delete [] I1_end_0;
        delete [] I2_start_0;
        delete [] I2_end_0;
        delete [] I1_l;
        delete [] I1_bar_l;
        delete [] I2_l;
        delete [] I2_bar_l;
        delete [] I1_start;
        delete [] I1_end;
        delete [] I2_start;
        delete [] I2_end;
        // cout << endl;
        return;
    }

    int N_intersect = 0;
    for(int i=1; i<=m; i++)
    {
        N_intersect += Indicator(I2_l[i], I2_end[i], n, I1_start, I1_end);
    }
    NI[chr_id] = N_intersect;

    DP(NM_chr, I1_l, I2_start, I2_l, I2_end, n, m, g, 1);
    double Sum1 = -1;
    for(int k=0; k<=m; k++) {
        Sum1 = logsum(Sum1, NM_chr[k]);
    }
    int nn = g-sum(I1_l, n)+n;
    double ratio = lognchoosek(nn, n) - Sum1;
    // cout << exp(ratio) << " ";

    double* CumSum = new double [m+1];
    for(int k=0; k<=m; k++) {
        if (k==0) {
            CumSum[m] = NM_chr[m];
        }
        else CumSum[m-k] = logsum(CumSum[m-k+1], NM_chr[m-k]);
    }

    double** q = new double* [n+1];
    for(int i=0; i<n+1; i++)
        q[i] = new double [m+1];
    double prod;
    // double* p = new double [m+1];
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
        Pj_chr[j] = 1-prod;
    }
    // poisbinom(Pj_chr[1], p, m);

    delete [] I1_start_0;
    delete [] I1_end_0;
    delete [] I2_start_0;
    delete [] I2_end_0;
    delete [] I1_l;
    delete [] I1_bar_l;
    delete [] I2_l;
    delete [] I2_bar_l;
    delete [] I1_start;
    delete [] I1_end;
    delete [] I2_start;
    delete [] I2_end;
    delete [] CumSum;
    for(int i=0; i<n+1; i++)
        delete [] q[i];
    delete [] q;
    // delete [] p;

    return;
}

void call_from_thread(int* n_chr, int* m_chr, double* NM_chr, int* NI, double* eta_chr, double* Pj_chr, 
                      string query_path, string reference_path, string chrom_size_path, int chr_id,
                      string id, int scale)
{
    chrom(n_chr, m_chr, NM_chr, NI, eta_chr, Pj_chr, query_path, reference_path, chrom_size_path,
                      chr_id, id, scale);
}

void combine_chr(double* p_value_g, int m, int num_chr, int* m_chr, double** NM)
{
    double** M;
    M = new double* [num_chr+1];
    for(int q=0; q<=num_chr; q++)
    {
        M[q] = new double [m+1];
        for(int k=0; k<=m; k++)
            M[q][k] = -1;
    }
    for(int q=1; q<=num_chr; q++)
    {
        for(int k=0; k<=m; k++)
        {
            if (k==0)
            {
                M[q][0] = 0;
                int flag = 0;
                for(int j=1; j<=q; j++){
                    M[q][0] += NM[j][0];
                    if (NM[j][0]==-1)
                        flag = 1;
                }
                if (flag==1)
                    M[q][0] = -1;
            }
            else if (q==1)
                M[1][k] = NM[1][k]; //could lower max_m_chr by taking care of this
            else
            {
                double temp = -1;
                for(int i=0; i<=min(k,m_chr[q]); i++)
                {
                    temp = M[q-1][k-i]+NM[q][i];
                    if (NM[q][i]==-1 || M[q-1][k-i]==-1)
                        temp = -1;
                    M[q][k] = logsum(M[q][k], temp);
                }
            }
        }
    }

    double Sum1 = -1;
    for(int k=0; k<=m; k++) {
        Sum1 = logsum(Sum1, M[num_chr][k]);
    }

    double* CumSum = new double [m+1];
    for(int k=0; k<=m; k++) {
        if (k==0) {
            CumSum[m] = M[num_chr][m];
        }
        else CumSum[m-k] = logsum(CumSum[m-k+1], M[num_chr][m-k]);
    }

    for(int k=0; k<=m; k++) {
        p_value_g[k] = exp(CumSum[k] - Sum1);
    }

    for(int i=0; i<=num_chr; i++)
    {
        delete [] M[i];
    }
    delete [] M;
    delete [] CumSum;
    return;
}

int main(int argc, char** argv)
{
    cout << HEADER;

    if (argc < 7)
    {
        cout << USAGE;
        return 1;
    }
    string query_path = argv[1];
    string reference_path = argv[2];
    string chrom_size_path = argv[3];
    int scale = stoi(argv[4]);
    string resultpath = argv[5];
    string thread_mode = argv[6];
    int num_chr = 24;
    cout << "scale is " << scale << endl;
    int max_m_chr = 100000;
    double** Pj;
    Pj = new double* [num_chr+1];
    for(int i=0; i<=num_chr; i++)
    {
        Pj[i] = new double [max_m_chr];
        for(int k=0; k<max_m_chr; k++)
        {
            Pj[i][k] = 0;
        }
    }

    double** NM; //log scale
    NM = new double* [num_chr+1];
    for(int i=0; i<=num_chr; i++)
    {
        NM[i] = new double [max_m_chr];
        for(int k=0; k<max_m_chr; k++)
            {
                NM[i][k] = -1;
            }
    }

    int* n_chr = new int [num_chr+1];
    int* m_chr = new int [num_chr+1];
    int* NI = new int [num_chr+1];
    double* eta_chr = new double [num_chr+1];

    thread *th = new thread [num_chr];

    if (thread_mode=="p")
    {//Mltithreading
        for(int i=0; i<num_chr; i++)
        {
            string id;
            if (i+1==23)
                id = "chrX";
            else if (i+1==24)
                id = "chrY";
            else
                id = "chr" + to_string(i+1);
            th[i] = thread(call_from_thread, n_chr, m_chr, NM[i+1], NI, eta_chr, Pj[i+1],
                           query_path, reference_path, chrom_size_path, i+1, id, scale);
        }
        for(int i=0; i<num_chr; i++)
        {
            th[i].join();
        }
    }
    else
    {// Single thread
        for(int i=1; i<=num_chr; i++)
        {
            string id;
            if (i==23)
                id = "chrX";
            else if (i==24)
                id = "chrY";
            else
                id = "chr" + to_string(i);
            chrom(n_chr, m_chr, NM[i], NI, eta_chr, Pj[i], query_path, reference_path, chrom_size_path,
                          i, id, scale);
        }
    }
    // int job_id[8] = {0, 1, 2, 3, 4, 5, 6, num_chr-1};
    // for (int j=0; j<num_chr/8; j++)
    // {
    //     for(auto i: job_id)
    //     {
    //         string id;
    //         if (i+1==23)
    //             id = "chrX";
    //         else if (i+1==24)
    //             id = "chrY";
    //         else
    //             id = "chr" + to_string(i+1);
    //         th[i] = thread(call_from_thread, n_chr, m_chr, NM[i+1], NI, Pj[i+1],
    //                        query_path, reference_path, chrom_size_path, i+1, id, scale);
    //     }
    //     for(auto i: job_id)
    //     {
    //         th[i].join();
    //     }
    //     job_id[0] += 7;
    //     job_id[1] += 7;
    //     job_id[2] += 7;
    //     job_id[3] += 7;
    //     job_id[4] += 7;
    //     job_id[5] += 7;
    //     job_id[6] += 7;
    //     job_id[7]--;
    // }

    int n = sum(n_chr, num_chr);
    int m = sum(m_chr, num_chr);
    int N_intersect = sum(NI, num_chr);
    cout << "Number of query intervals (n) is " << n << endl;
    cout << "Number of reference intervals (m) is " << m << endl;
    cout << "Overlap is " << N_intersect << endl;

    // int mid = num_chr/2;
    // nth_element(eta_chr+1, eta_chr+1+mid, eta_chr+num_chr+1);
    // double eta = eta_chr[1+mid];
    double eta = *min_element(eta_chr+1, eta_chr+num_chr+1);
    cout << "eta is " << eta << endl;

    double* p_value_dp = new double [m+1];
    combine_chr(p_value_dp, m, num_chr, m_chr, NM);
    cout << "DP p-value is " << p_value_dp[N_intersect] << endl;

    double* p = new double [m+1];
    int idx = 0;
    for(int i=1; i<=num_chr; i++)
    {
      for(int k=1; k<=m_chr[i]; k++)
      {
        p[idx+k] = Pj[i][k];
      }
      idx += m_chr[i];
    }
    double* p_value_poisbinom = new double [m+1];
    poisbinom(p_value_poisbinom, p, m);
    cout << "PB p-value is " << p_value_poisbinom[N_intersect] << endl;

    // write2file(resultpath, p_value_dp, p_value_poisbinom, m);
    write2file(resultpath, p_value_dp, p_value_poisbinom, m);

    for(int i=0; i<=num_chr; i++)
    {
        delete [] Pj[i];
    }
    delete [] Pj;
    for(int i=0; i<=num_chr; i++)
    {
        delete [] NM[i];
    }
    delete [] NM;
    delete [] n_chr;
    delete [] m_chr;
    delete [] NI;
    delete [] eta_chr;
    delete [] th;
    delete [] p_value_dp;
    delete [] p;
    delete [] p_value_poisbinom;
    return 0;
}
