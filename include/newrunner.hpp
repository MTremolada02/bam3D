#ifndef newrunner_hpp
#define newrunner_hpp

#include <htslib/sam.h>
#include <map>

struct UserInputBam3D : UserInput { // additional input
	bool hist_none        =false;
	bool hist_global      = false;
	bool hist_by_chrom    = false;

	uint8_t decompression_threads = 4;
	uint8_t histogram_mode = hist_none; // histogram mode

};

struct ReadStats {
    uint64_t readN   = 0;
    uint64_t qc_fail = 0;
    uint64_t unmapped = 0;
    uint64_t secondary = 0;
    uint64_t supplementary = 0;
    uint64_t primary = 0;
    uint64_t mapQ0 = 0;

	long double mean_insert=0;
	long double quadratic_mean=0;

	double error_rate=0;
};

struct PairStats {
    uint64_t pairN        = 0;
    uint64_t proper_pairs = 0;
    uint64_t good_pairs   = 0;
	bool good_read1=false; //dovrebbe essere locale
	bool good_read2=false; //dovrebbe essere locale

    uint64_t UMone_sided  = 0;
    uint64_t UMtwo_sided  = 0;
    uint64_t duplicated   = 0;
    uint64_t sameCr       = 0;  // cis

    uint64_t read1 = 0;
    uint64_t read2 = 0;
};

enum class Maptype :uint8_t {N=0, U=1, M=2, R=3};

struct QnameStats { 
	Maptype type1;
	Maptype type2;

	uint64_t UU=0;
	uint64_t RU=0;
	uint64_t UR=0;
	uint64_t WW=0;
	uint64_t DD=0;
	uint64_t MU=0;
	uint64_t MM=0;
	uint64_t MR=0;
	uint64_t NM=0;
	uint64_t NU=0;
	uint64_t NR=0;
	uint64_t NN=0;
};

class Runner {
    
    UserInputBam3D userInput;
	ReadStats readStats;
	PairStats pairStats;
	QnameStats qnameStats;
    
public:
    
    void loadInput(UserInputBam3D userInput);
	long double update_mean_tlen(long double,uint64_t, bam1_t*);
	long double update_quadratic_mean_tlen(long double,uint64_t, bam1_t*);
	double error_rate(uint64_t,uint64_t);
	uint16_t Alignstarts(const bam1_t*);
	void qname_stats(std::vector<bam1_t*> &,int);
	void flag_inspector(bam1_t*);
	void histo_global_distance(std::unordered_map<uint64_t, uint64_t>&);
	void histo_chrom_distance(std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>>&); 
    int data_vector(std::vector<std::vector<bam1_t*>> &, bam1_t* &, bool &, std::vector<int> &, samFile *, bam_hdr_t *);
	void processReads(bam_hdr_t*, std::vector<std::vector<bam1_t*>> &, int , std::vector<int>&);
	void output();
	void run();
    
};

#endif /* newrunner_hpp */