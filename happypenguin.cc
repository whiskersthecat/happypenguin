// happypenguin v2.0.0
// Peter Reifenstein
// February 25, 2025

#include "happypenguin.h"
const string VERSION = "2.0.1";
#define NUM_CATEGORIES 9
#define NUM_RECOMBINANT_CATEGORIES 3
#define CHOMP_CATEGORIES_INDEX 4

// version 2.0.0. added consensus matrix
// version 2.0.1. changed 1 and X to combine to x (instead of 1)
// version 2.0.1. do not flush buffer every iteration. Write '\n' token instead of endl, buffer will flush (write to file) automatically when it is full

unsigned int*** SNP_coordinate;
char*** SNP_nucleotide;

int num_chromosomes;

ofstream*** output_SAM_files;
ofstream output_TXT_files[NUM_CATEGORIES + 1];
ofstream summary_file;
ofstream log_file;
ofstream recombinantinfo_files[NUM_RECOMBINANT_CATEGORIES];
ofstream readlength_files[NUM_CATEGORIES + 1];
ofstream error_file;

vector<tuple<long long,string>> recombinantinfo_vector[NUM_RECOMBINANT_CATEGORIES];

long genome1size = 0;

ifstream sam_files[2];

unsigned long long totalErrors = 0;
unsigned long long totalErrorBP = 0;

map<string, int> chromosome_numbers;
map<int, string> chromosome_names[2];

static string read_SNP_symbols[]      = {"1","2","O","X","-"};
// static string consensus_SNP_symbols[] = {"1","2","O","X","-","!","?","$"};


static char read_to_concensus[5][5] = {{'1','!','?','x','-'},
                                          {'!','2','?','x','-'},
                                          {'?','?','O','$','-'},
                                          {'x','x','$','X','-'},
                                          {'-','-','-','-','-'}
                                        };

// static string SNP_meanings[] = {"Genotype1","Genotype2","BothMismatch","NeitherMismatch"};
static char ALIGN_DIRECTIONS[] = {'F','R'};
static map<char, short> GT_ID = {{'1' , 0} , {'2' , 1}};

static string CATEGORY_NAMES[] = {"Genotype1", "Genotype2", "LowQualityCrossover", "MediumQualityCrossover", "HighQualityCrossover","ChompBadAlign", "ChompNoSNPs", "ChompBadGT", "ChompMultimapRec"};
int num_reads_in_category[NUM_CATEGORIES];

int main(int argc, char* argv[]) {
    if(argc < 4) {
        printProgramInfo();
    }
    string dir_string = "happypenguin_" + VERSION + "_output";
    if(argc > 4) {
        string append = argv[4];
        dir_string += "_on_" + append;
    }

    makeDirectories(dir_string);
    readChromosomes(argv[1], argv[2]);
    readSNPs(argv[3]);

    auto begin = chrono::high_resolution_clock::now();
    analyzeReads();
    cout << "[happypenguin] Done analyzing reads" << endl;

    auto end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::microseconds>(end - begin).count();

    printCategorizationResults((double)elapsed / (double) 1000000);

    memoryCleanup();

    generatePlots(dir_string);
}

void categorizeRead(int categorization, string& read_name, string& sam_line_A, string& sam_line_B, string& textoutput, bool silenced, int len) {

    // 1. print textoutput to standard output, if not silenced
    textoutput.append(" ***** Read is " + CATEGORY_NAMES[categorization] + "\n");
    if(!silenced) cout << textoutput << endl;

    // 2. increment counter based on categorization
    num_reads_in_category[categorization]++;

    // 3. write the SAM files with the read
    *output_SAM_files[0][categorization] << read_name << sam_line_A << '\n';
    *output_SAM_files[1][categorization] << read_name << sam_line_B << '\n';

    // 4. write the text file with the read and the all reads text file
    output_TXT_files[categorization] << textoutput << '\n';
    output_TXT_files[NUM_CATEGORIES] << textoutput << '\n';

    // 5. write the length of the read to its corresponding file
    readlength_files[categorization] << len << '\n';
    readlength_files[NUM_CATEGORIES] << len << '\n';

    return;
}


string reconstructAlignedSequence(stringstream& read, int& len, unsigned long& loc, string& chrom_name, bool& dir) {
    // Lsat_Sal__Chrom2_(7016910..7020531)	2064	SERH_U_Chr_2	35850229	1	4=1X82=1X79=1X60=1X111=1X20=1X3=1X97=1X194=1X12=1X28=1X12=1X21=1X15=1X3=1I49=1X10=1X1=1X144=1D76=1X9=1X1=1X5=2565H
    // readname                             flag    chrom_aligned   position    qual    CIGARstring

    int num, flag, mapq;
    string cigar_str;
    string seq;

    string str;

    read >> flag;
    read >> chrom_name;
    read >> loc;
    read >> mapq;

    read >> cigar_str; // field 6
    read >> str; 
    read >> num;
    read >> num;
    read >> seq;   // field 10

    dir = flag & 16; // check if 16 present in flag, indicating aligned reverse direction

    stringstream cigar(cigar_str);
    char token;
    int val;
    string align_seq;
    int seq_pos = 0;
    
    if(!(cigar.str() == "*")) {
        while (cigar) {
            cigar >> val;
            cigar >> token;
            if (token == 'M' || token == '=' || token == 'X') {
                align_seq += seq.substr(seq_pos, val);
                seq_pos += val;
            } 
            else if (token == 'I') {
                // ignore and consume seq
                seq_pos += val;
            }
            else if (token == 'D' || token == 'N') {
                // seq deleted from query
                for (int i = 0; i < val; i++) {
                    align_seq += "X";
                }
            }
            else if (token == 'S') {
                // consume seq
                seq_pos += val;
            }
            else if (token == 'H' || token == 'P') {
                // skip token
            }
            else {
                throw("[cigar_parse] Invalid cigar token: " + token);
            }
        }
    }
    len = align_seq.length();

    return align_seq;

}

void analyzeReads() {
    
    // the header fields are already consumed when finding the chromosome sizes

    cout << "[happypenguin] Analyzing reads" << endl;

    string current_read_name = "";
    string next_read_name = "";

    sam_files[0] >> next_read_name;
    sam_files[1] >> next_read_name;

    int read_number = 0;
    
    while(sam_files[0]) {

        string textoutput;

        bool silenced = true;
        if(read_number % 1000 == 0)
            silenced = false;
        read_number++;

        current_read_name = next_read_name;
        textoutput.append("[" + to_string(read_number) + "] Processing Read " + current_read_name + '\n');


        // 1. reconstruct aligned sequences for both alignments

        string align_seq[2];
        string line;
        int align_len[2];
        unsigned long align_loc[2];
        unsigned long align_end_loc[2];
        int align_chrom[2];
        string align_chrom_name[2];
        int n_alternate_alignments[2];
        string SAM_line[2];
        bool align_dir[2];

        for(int i = 0; i < 2; i++) {
            getline(sam_files[i], SAM_line[i]);
            stringstream read(SAM_line[i]);

            align_seq[i] = reconstructAlignedSequence(read, align_len[i], align_loc[i], align_chrom_name[i], align_dir[i]);

            // parseCIGAR(read, (mmlocations[i]), align_len[i], align_loc[i], align_chrom_name[i], align_dir[i]);
            align_end_loc[i] = align_loc[i] + align_len[i];
            align_chrom[i] = chromosome_numbers[align_chrom_name[i]];

            n_alternate_alignments[i] = 0;
            sam_files[i] >> next_read_name;
            while(next_read_name == current_read_name && sam_files[i]) { // stop looking for alternates if a new read name comes up or end of file reached
                // alternate alignment
                n_alternate_alignments[i]++;
                getline(sam_files[i], line);
                sam_files[i] >> next_read_name;
            }
        }

        // 2. look up SNPs on genome 1 alignment, load into memory

        vector<int> variantLocations[2];
        vector<char> variantNucleotides[2];

        for(unsigned int pos = align_loc[0]; pos < align_end_loc[0]; ++pos) {
            if(SNP_coordinate[0][align_chrom[0]][pos] != 0) {
                // find coordinates of variant site along both alignments
                unsigned int relative_variant_pos[2];
                relative_variant_pos[0] = pos - align_loc[0];
                relative_variant_pos[1] = SNP_coordinate[0][align_chrom[0]][pos] - align_loc[1];

                // record nucleotides on both genomes
                variantLocations[0].push_back(relative_variant_pos[0]);
                variantLocations[1].push_back(relative_variant_pos[1]);
                variantNucleotides[0].push_back(SNP_nucleotide[0][align_chrom[1]][pos]);
                variantNucleotides[1].push_back(SNP_nucleotide[1][align_chrom[1]][pos]);
            }
        }

        for (int i = 0; i < 2; i++) {
            textoutput.append("Mapped to chromosome " + to_string(align_chrom[i] + 1)  + " from position " + to_string(align_loc[i]) + " to " + to_string(align_end_loc[i]) + " (Total Length " + to_string(align_len[i]) + ")");
            if(n_alternate_alignments[i] > 0)
                textoutput.append(" with " + to_string(n_alternate_alignments[i]) + " alternate mappings(s)");
            textoutput.append("\n");
        }

        int numSNPs = variantLocations[0].size();

        textoutput.append("SNP Locations (Total " + to_string(numSNPs) + "): ");
        for(int i = 0; i < numSNPs; i++) {
            textoutput.append(to_string(variantLocations[0][i]) + "/" + to_string(variantLocations[1][i]) + ":" + (variantNucleotides[0][i]) + "/" + (variantNucleotides[1][i]) + "\t");
        }
        textoutput.append("\n");

        // 3. extract nucleotides at SNP locations for both alignments (A, C, G, T, X, or -), and call their genotype
        char readNucleotides[2][numSNPs];
        short readSNPs[2][numSNPs];
        for(int n = 0; n < numSNPs; n++) {
            for(int i = 0; i < 2; i++) {
                if (i > 0 && (variantLocations[i][n] < 0 || variantLocations[i][n] > align_len[i])) {
                    // this SNP position is out of bounds
                    // indicates read did not align in similar location that genomic SNP calling
                    readNucleotides[i][n] = '-';
                    readSNPs[i][n] = 4;
                }
                else {
                    readNucleotides[i][n] = align_seq[i][variantLocations[i][n]];
                    if(readNucleotides[i][n] == 'X') readSNPs[i][n] = 3;
                    else if(readNucleotides[i][n] == variantNucleotides[0][n]) readSNPs[i][n] = 0;
                    else if(readNucleotides[i][n] == variantNucleotides[1][n]) readSNPs[i][n] = 1;
                    else readSNPs[i][n] = 2;
                }
            }
        }

        for (int i = 0; i < 2; i++) {
            textoutput.append("  Genotype Matrix:  ");
            for (int n = 0; n < numSNPs; n++) {
                textoutput.append(read_SNP_symbols[readSNPs[i][n]]);
            }
            textoutput.append("\n");
        }

        // 4. generate concensus matrix
        char readConsensus[numSNPs];
        short numGTswitches = -1;
        short breakpoint = -1;
        char lastSNP = 'N';
        char newSNP = 'N';
        char firstSNP = 'N';
        
        int numUsefulSNPs = 0;
        int numIgnoredSNPs = 0;
        for (int n = 0; n < numSNPs; n++) {
            newSNP = read_to_concensus[readSNPs[0][n]][readSNPs[1][n]];
            readConsensus[n] = newSNP;
            if(newSNP == '1' || newSNP == '2') {
                if(firstSNP == 'N')
                    firstSNP = newSNP;
                numUsefulSNPs ++;
                if(newSNP != lastSNP) {
                    lastSNP = newSNP;
                    numGTswitches ++;
                    breakpoint = numUsefulSNPs;
                }
            } else
                numIgnoredSNPs ++;
        }

        textoutput.append(" Consensus Matrix:  ");
        for (int n = 0; n < numSNPs; n++) {
            textoutput.append(1, readConsensus[n]);
        }
        textoutput.append("\n");

        // 5. categorization
        // #CHOMP CATEGORY 1
        if(align_chrom[0] != align_chrom[1]) {
            categorizeRead(CHOMP_CATEGORIES_INDEX + 1, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }
        // #CHOMP CATEGORY 2
        if(numUsefulSNPs < 1) {
            categorizeRead(CHOMP_CATEGORIES_INDEX + 2, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }
        // #CHOMP CATEGORY 3
        if(numGTswitches > 1) {
            categorizeRead(CHOMP_CATEGORIES_INDEX + 3, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }

        // #GENOTYPE 1 or GENOTYPE 2
        if(numGTswitches == 0) {
            // 1. Count Errors
            // TO IMPLEMENT
            // totalErrors += mmlocations[1].size() - expectedSNPlocations[1].size();
            // totalErrorBP += align_len[1];

            // 2. Categorize
            categorizeRead(GT_ID[firstSNP], current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }
        
        // ### RECOMBINANT

        // Check for alternate mappings
        bool alternate_mappings = false;
        if(n_alternate_alignments[0] > 0 || n_alternate_alignments[1] > 0) {
            alternate_mappings = true;
        }

        // #CHOMP CATEGORY 4
        if(alternate_mappings) {
            categorizeRead(CHOMP_CATEGORIES_INDEX + 4, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);
            continue;
        }

        int quality = 0;

        // #LOW QUALITY RECOMBINANT
        if(breakpoint == 2 || breakpoint == numUsefulSNPs) quality = 0;

        // #MEDIUM QUALITY RECOMBINANT
        else if(breakpoint == 3 || breakpoint == numUsefulSNPs - 1) quality = 1;

        // #HIGH QUALITY RECOMBINANT
        else quality = 2;

        // write to recombinant info

        // "name\tchromosome\tcoordinate1\tcoordinate2\tfirstSNP\tnumUsefulSNPs\tnumIgnoredSNPs\tbreakpoint\talignDir(A:B)\tnumAlternateAlignments(A:B)" << endl;

        stringstream rinfo;
        rinfo << current_read_name << '\t' << (align_chrom[0] + 1) << '\t' << setw(10) << align_loc[0] << '\t' << setw(10) << align_loc[1] << '\t' << firstSNP << '\t' << numUsefulSNPs << '\t' << numIgnoredSNPs << '\t' << breakpoint << '\t' << ALIGN_DIRECTIONS[align_dir[0]] << ':' << ALIGN_DIRECTIONS[align_dir[1]] << '\t' << n_alternate_alignments[0] << ':' << n_alternate_alignments[1];
        if(align_dir[0] != align_dir[1]) log_file << "[anomolous case] Crossover read " << current_read_name << " has alignments in different directions" << endl;
        categorizeRead(2 + quality, current_read_name, SAM_line[0], SAM_line[1], textoutput, silenced, align_len[1]);

        long long coordinatevalue = genome1size * align_chrom[0] + align_loc[0];
        recombinantinfo_vector[quality].push_back(make_tuple(coordinatevalue,rinfo.str()));

        if(!silenced) cout << textoutput << endl;

    }

    for(int i = 0; i < NUM_RECOMBINANT_CATEGORIES; i++) {
        sort(recombinantinfo_vector[i].begin(), recombinantinfo_vector[i].end());
        for(auto it = recombinantinfo_vector[i].begin(); it !=recombinantinfo_vector[i].end(); it++ ) {
            recombinantinfo_files[i] << get<1>(*it) << endl;
        }
        
    }
    
    

    error_file << "Total Error Rate: " << 100 * ((double) totalErrors / (double) totalErrorBP) << "%" << endl;

    return;
}

void printProgramInfo() {
    const char* info = "\nhappypenguin version 2.1\nPeter Reifenstein Feb 2025\n            _____        ____\n           / (^) \\______/ (^) \\ \n          /      ________      \\ \n/\\__/\\____|__o___\\______/__o___|____/\\/\\/\\ \n\nA tool to classify Nanopore or HiFi DNA reads as belonging to Genome A, Genome B, or a Crossover\n\ncompile:\ng++ happypenguin.cc -std=c++11 -O3 -o happypenguin\n\nusage: \n./happypenguin alignment_A.sam alignment_B.sam variations.tsv [output_name]\n\n-- alignment_A.sam\nincludes all reads aligned to genome A, with match/mismatch information in the CIGAR string\nfor example, to align using the minimap2 aligner (https://github.com/lh3/minimap2):\n./minimap2 --eqx -a --secondary=no genomeA.fasta reads.fasta\n\n-- variations.tsv\na tab seperated value file with the format specified by SYRI (https://schneebergerlab.github.io/syri/fileformat.html):\n  genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ \nregion_type beginning with \"SYN\" and variation_type = \"SNP\" indicates a SNP in a syntenic region\ngenomeA_chrom and genomeA_pos give the coordinates of the SNP in genome A\nfields marked with ~ are not used and can be any string\n";
    cout << info << endl;
    exit(1);
}

void printCategorizationResults(double time) {
    stringstream summary;
    long total_reads = 0;
    long total_discarded_reads = 0;

    // Get total reads
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        total_reads += num_reads_in_category[i];
    }

    summary << left;
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        double percentage = 100 * (double) num_reads_in_category[i] / (double) total_reads;
        if(i > CHOMP_CATEGORIES_INDEX) {
            if(i == CHOMP_CATEGORIES_INDEX + 1)
                summary << endl;
            total_discarded_reads += num_reads_in_category[i];
            summary << "Total " << setw(26) << (CATEGORY_NAMES[i] + ":") << setw(12) << num_reads_in_category[i] << fixed << setprecision(3) << percentage << "%"  << endl;
        } else {
            summary << "Total " << setw(26) << (CATEGORY_NAMES[i] + ":") << setw(12) << num_reads_in_category[i] << fixed << setprecision(3) << percentage  << "%"  << endl;
        }
    }

    double attrition_rate = (double)total_discarded_reads / (double)total_reads;

    summary << endl << "Processed [" << total_reads << "] reads, " << 100 * attrition_rate << "% were thrown out (chomped)" << endl;

    
    summary << "Elapsed time: " << time << " seconds" << endl;

    string summary_str = summary.str();
    cout << summary_str << endl;
    summary_file << summary_str << endl;

    return;
}

void makeDirectories(string base_dir) {
    mkdir((base_dir).c_str(), 0777);

    string alignments_dir = base_dir + "/sorted_alignments";
	mkdir((alignments_dir).c_str(), 0777);

    // open ofstreams for writing SAM files
    output_SAM_files = (ofstream***) malloc(2 * sizeof(ofstream**));

    for(int i = 0; i < 2; i++) {
        string SAM_names_prefix = "toGenome";
        SAM_names_prefix += read_SNP_symbols[i];
        output_SAM_files[i] = (ofstream**) malloc(NUM_CATEGORIES * sizeof(ofstream*));
        string gt_dir = alignments_dir + "/" + SAM_names_prefix;
        mkdir((gt_dir).c_str(), 0777);
        for(int j = 0; j < NUM_CATEGORIES; j++) {
            output_SAM_files[i][j] = new ofstream(gt_dir + "/" + SAM_names_prefix + "_" + CATEGORY_NAMES[j] + ".SAM");
        }
    }

    // open ofstreams for writing to text files
    string txt_dir = base_dir + "/text_files";
	mkdir((txt_dir).c_str(), 0777);
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        output_TXT_files[i].open(txt_dir + "/" + CATEGORY_NAMES[i] + ".txt");
    }
    output_TXT_files[NUM_CATEGORIES].open(txt_dir + "/" + "AllReads" + ".txt");

    // open ofstream for writing summary
    summary_file.open(base_dir + "/summary.txt");

    // open ofstream for writing log
    log_file.open(base_dir + "/log.txt");

    // open ofstreams for writing recombinant info
    string info_dir = base_dir + "/crossover_info";
	mkdir((info_dir).c_str(), 0777);
    for(int i = 0; i < NUM_RECOMBINANT_CATEGORIES; i++) {
        recombinantinfo_files[i].open(info_dir + "/" + CATEGORY_NAMES[i + 2] + "Info.tsv");
        // recombinantinfo_files[i] << "name alignmentSNPs(A:B)  numSNPs(A:B)    relativeBreakpoint  numIgnoredSNPs  numAlternateAlignments(A:B)   firstSNP    alignDir(A:B)   chromosome(A:B)   len(A:B)    alignloc(A:B)" << endl;
    
        recombinantinfo_files[i] << "name\tchromosome\tcoordinate1\tcoordinate2\tfirstSNP\tnumSNPS\tnumIgnoredSNPs\tbreakpoint\talignDir(A:B)\tnumAlternateAlignments(A:B)" << endl;

    }

    string misc_dir = base_dir + "/misc";
    mkdir((misc_dir).c_str(), 0777);

    // open ofstreams for writing lengths of reads
    string len_dir = misc_dir + "/read_lengths";
    mkdir((len_dir).c_str(), 0777);
    for(int i = 0; i < NUM_CATEGORIES; i++) {
        readlength_files[i].open(len_dir + "/" + CATEGORY_NAMES[i] + "_lengths" + ".tsv");
    }
    readlength_files[NUM_CATEGORIES].open(len_dir + "/" + "AllReads_lengths" + ".tsv");

    // open ofstream for writing error rate
    error_file.open(misc_dir + "/Error_Rate.txt");

    return;
}

void createSNPfile(string polymorphism_file) {
    ifstream p_file(polymorphism_file);
    ofstream p_SNP_file(polymorphism_file + ".happypenguin.SNPs");
    ofstream p_info_file(polymorphism_file + ".happypenguin.info");
    ofstream p_VCF_files[2];
    for(int i = 0; i < 2; i++) {
        p_VCF_files[i].open(polymorphism_file + ".happypenguin.SNPs" + ".genome" + read_SNP_symbols[i] + ".VCF");
        p_VCF_files[i] << "##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" << endl;
    }
    p_SNP_file << "chromosome_num\tgenome_1_location\tgenome_2_location\tgenome_1_nucleotide\tgenome_2_nucleotide" << endl;

    if(!p_file) { cout << "[ERROR] cannot open " << polymorphism_file << endl; exit(1); }
    string chrBloc_str, unique_id, parent_id, type, prev_parent_id, str;
    string seq[2];
    string chrName[2];
    //int chrAloc, chrBloc, chrAnum, chrBnum, chrAname, chrBname;
    int chrLoc[2];
    int chrNum[2];
    unsigned long totalcount = 0;
    int num_SNPs_on_chrom [num_chromosomes] = {0};
    int num_syntenic_regions_on_chrom [num_chromosomes] = {0};
    while(p_file >> chrName[0]) {
        p_file >> chrLoc[0]; p_file >> chrLoc[0]; p_file >> seq[0]; p_file >> seq[1];
        p_file >> chrName[1]; p_file >> chrBloc_str;  p_file >> chrBloc_str;
        p_file >> unique_id; p_file >> parent_id; p_file >> type; p_file >> str;

        // SALH_C_Chr_1	1	26851	-	-	-	-	-	NOTAL1	-	NOTAL	-
        if(parent_id.substr(0, 3) == "SYN" && type == "SNP") {
            // found SNP in syntenic region
            chrLoc[1] = stoi(chrBloc_str);
            totalcount++;
            chrNum[0] = chromosome_numbers[chrName[0]];
            chrNum[1] = chromosome_numbers[chrName[1]];
            num_SNPs_on_chrom[chrNum[0]] ++;
            if(chrNum[0] != chrNum[1]) {
                log_file << "[createSNPfile] [warning] Different chromosomes in 'syntenic region'" << endl;
            }

            // write to SNP file
            p_SNP_file << chrNum[0] + 1 << '\t' << chrLoc[0] << '\t' << chrLoc[1] << '\t' << seq[0] << '\t' << seq[1] << endl;

            // COMMENT BACK!
            // // write to both VCF files
            // for(int i = 0; i < 2; i++) {
            //     // CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            //     p_VCF_files[i] << (chromosome_names[i][chrNum[i]] + "\t" + to_string(chrLoc[i]) + "\t" + unique_id + "\t" + seq[i] + "\t" + seq[i ^ true] + "\t.\tPASS\t." ) << endl;
            // }

            if(parent_id != prev_parent_id) {
                num_syntenic_regions_on_chrom[chrNum[0]]++;
            }
            prev_parent_id = parent_id;
            //if(totalcount % 400000 == 0)
            //    cout << "[SNP parse] Found SNP in syntenic region on chromosome " << chrAnum << " A coordinate: " << chrAloc << " B coordinate: " << chrBloc << endl;
        }
    }
    for(int i = 0; i < num_chromosomes; i++) {
        p_info_file << "Chromosome " << i + 1 << " has total " << num_SNPs_on_chrom[i] << " SNPs among " << num_syntenic_regions_on_chrom[i] << " syntenic regions" << endl;
    }
    
}

void readSNPs(string polymorphism_file) {

    // Lsat_PI251246_v14_Chr1  3435181 3435181 C       T       Lsat_PI251246_v14_Chr4  360187218       360187218       SNP13644        INVTR7135       SNP     -
    // ref_chrom               ref_pos ref_pos ref_seq qry_seq qry_chrom               qry_pos         qry_pos         unique_id       parent_id       type    copy-status
    //genomeA_chrom  genomeA_pos  ~  ~  ~  genomeB_chrom  genomeB_pos  ~  ~  region_type  variation_type  ~ \nwhere the region_type beginning with "SYN" indicates a syntenic region and variation_type = "SNP" indicates a SNP\nand genomeA_chrom and genomeA_pos give the coordinates of the SNP in genome A\n(the other fields marked with ~ are not used and can be any string)
    //
    //
    // (the other fields marked with ~ are not used and can be any string)
    // look for parent_id = SNYx (i.e. the variation occurs in the xth syntenic alignment region) and type = SNP   

    cout << "[happypenguin] Reading in Polymorphism file" << endl;

    ifstream p_SNP_file_test(polymorphism_file + ".happypenguin.SNPs");
    if(!p_SNP_file_test) {
        cout << "[SNP parse] No SNP file created yet from polymorphisms file, creating SNP file and VCF files..." << endl;
        createSNPfile(polymorphism_file);
    }

    cout << "[happypenguin] Reading in SNP file" << endl;
    ifstream p_SNP_file(polymorphism_file + ".happypenguin.SNPs");
    string header;
    getline(p_SNP_file, header);
    int chr1loc, chr2loc;
    char nuc1, nuc2;
    int chromosome;
    int count = 0;
    // all three arrays are indexed by location on the first genome
    while(p_SNP_file >> chromosome >> chr1loc >> chr2loc >> nuc1 >> nuc2) {
        SNP_coordinate[0][chromosome - 1][chr1loc] = chr2loc;
        SNP_nucleotide[0][chromosome - 1][chr1loc] = nuc1;
        SNP_nucleotide[1][chromosome - 1][chr1loc] = nuc2;
        // SNP_coordinates[1][chromosome - 1][chrBloc] = 1;
        count ++;
    }
}

void readChromosomes(string alignment_genome_Afile, string alignment_genome_Bfile) {
    sam_files[0].open(alignment_genome_Afile);
    sam_files[1].open(alignment_genome_Bfile);

    if(!sam_files[0] || !sam_files[1]) { cout << "[ERROR] cannot open one or more SAM files" << endl; exit(1); }

    cout << "[happypenguin] Reading chromosomes from SAM files" << endl;
    vector<int> genome_sizes[2];
    string chrom_name;

    char c;
    string str;
    int len;
    int n = 0;
    for(int i = 0; i < 2; i++) {
        n = 0;
        while(true) {
            sam_files[i] >> c;
            if(c!= '@') {
                sam_files[i].unget();
                break;
            }
            
            getline(sam_files[i], str);
            str = "@" + str;
            stringstream header_line(str);
            header_line << str;
            for(int k = 0; k < NUM_CATEGORIES; k++) {
                *output_SAM_files[i][k] << str << endl;
            }

            header_line >> str;
            if(str == "@SQ") { // genomic alignment reference sequence (i.e. chromosome)
                header_line >> str;
                str = str.substr(3);
                header_line >> c;header_line >> c;header_line >> c;
                header_line >> len;

                log_file << "[chromosome] " << str << " len: " << len << endl;
                chromosome_numbers[str] = n;
                chromosome_names[i][n] = str;
                genome_sizes[i].push_back(len);
                n++;
            } else { // consume the rest of the field
                getline(header_line, str); 
            }
            
        }
    }
    num_chromosomes = genome_sizes[1].size();

    genome1size = accumulate(genome_sizes[0].begin(), genome_sizes[0].end(), 0);

    cout << "[happypenguin] Allocating Memory for SNP locations and SNP nucleotides" << endl;

    SNP_coordinate    = (unsigned int***) malloc(2 * sizeof(unsigned int**));
    SNP_coordinate[0] = (unsigned int**)  malloc(num_chromosomes * sizeof(unsigned int*));
    SNP_coordinate[1] = (unsigned int**)  malloc(num_chromosomes * sizeof(unsigned int*));

    SNP_nucleotide    = (char***) malloc(2 * sizeof(char**));
    SNP_nucleotide[0] = (char**)  malloc(num_chromosomes * sizeof(char*));
    SNP_nucleotide[1] = (char**)  malloc(num_chromosomes * sizeof(char*));

    for(int i = 0; i < num_chromosomes; i++) {
        // memory allocated with calloc will be zero initialized
        SNP_coordinate[0][i] = (unsigned int*) calloc((genome_sizes[0][i] + 1), sizeof(unsigned int));
        SNP_nucleotide[0][i] = (char*) calloc((genome_sizes[0][i] + 1), sizeof(char));
        SNP_nucleotide[1][i] = (char*) calloc((genome_sizes[0][i] + 1), sizeof(char));
        log_file << "[calloc] Allocated memory for chromosome " << i + 1 << " SNPs (total " << genome_sizes[0][i] + 1 << " bytes)" << endl;
    }
    return;
}

void memoryCleanup() {
    log_file << "[happypenguin] Deallocating heap memory" << endl;
    for(int i = 0; i < num_chromosomes; i++) {
        free(SNP_coordinate[0][i]);
        free(SNP_nucleotide[0][i]);
        free(SNP_nucleotide[1][i]);
        // free(SNP_coordinate[1][i]);
    }
    free(SNP_coordinate[0]);
    free(SNP_coordinate[1]);
    free(SNP_coordinate);

    free(SNP_nucleotide[0]);
    free(SNP_nucleotide[1]);
    free(SNP_nucleotide);

    for(int i = 0; i < NUM_CATEGORIES; i++) {
        free(output_SAM_files[0][i]);
        free(output_SAM_files[1][i]);
    }
    free(output_SAM_files[0]);
    free(output_SAM_files[1]);
    free(output_SAM_files);

    log_file << "[happypenguin] Deallocated heap memory" << endl;
    return;
}

void generatePlots(string base_dir) {
    string command = "python3 generatePlots.py " + base_dir;
    cout << "[happypenguin] Generating plots" << endl;
    int result = system(command.c_str());
    if(result == 1) {
        cout << "[happypenguin] Error generating plots" << endl;
    } else {
        cout << "[happypenguin] Plots generated" << endl;
    }
}