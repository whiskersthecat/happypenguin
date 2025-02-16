#ifndef happypenguin
#define happypenguin

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <set>
#include <map>
#include <sstream>
#include <cstring>
#include <chrono>
#include <experimental/filesystem>
#include <sys/types.h>
#include <sys/stat.h>
#include <iomanip>
#include <chrono>
using namespace std;

void printProgramInfo();
void printCategorizationResults(double time);
void makeDirectories(string base_dir);
void readSNPs(string polymorphism_file);
void readChromosomes(string alignment_genome_Afile, string alignment_genome_Bfile);
void analyzeReads();
void memoryCleanup();
void parseCIGAR(stringstream& read, vector<int>& mm_locs, int& len, unsigned long& loc, string& chrom_name, bool& dir);
void categorizeRead(int categorization, string& read_name, string& sam_line_1, string& sam_line_2, string& textoutput, bool silenced, int len);
void createSNPfile(string polymorphism_file);
void generatePlots(string base_dir);

#endif
