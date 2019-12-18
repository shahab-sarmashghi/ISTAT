#ifndef PREP_H
#define PREP_H

int scan_input(std::string path, std::string id);
// scans the input file to count the number of intervals 

void read_input(unsigned int* I_start, unsigned int* I_end, std::string path, std::string id);
// reads the input and stores the start and end of intervals

int merge_adj(unsigned int* I_start, unsigned int* I_end, int n_int);
// merges adjacent intervals

void compute_length(int* I_l, int* I_bar_l, unsigned int* I_start, unsigned int* I_end, int n_int);
// computes the length of intervals and gaps based on start and end of intervals

void scale_interval(unsigned int* I_start, unsigned int* I_end, int* I_l, int* I_bar_l, int n_int, int scale);
// scales intervals and gaps through dividing the length by scale

#endif