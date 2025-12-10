#ifndef runner_hpp
#define runner_hpp

#include <htslib/sam.h>

struct UserInputBam3D : UserInput { // additional input
	uint8_t decompression_threads = 4;
};

struct GraphStats{
	std::vector<long double> edge; //Histo
    std::vector<uint64_t> counts; //Histo
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

class Runner {
    
    UserInputBam3D userInput;
	ReadStats readStats;
	PairStats pairStats;
	GraphStats graphStats;
    
public:
    
    void loadInput(UserInputBam3D userInput);
	void bam_sorting();
	long double update_mean_tlen(long double,uint64_t, bam1_t*);
	long double update_quadratic_mean_tlen(long double,uint64_t, bam1_t*);
	double error_rate(uint64_t,uint64_t);
	void qname_group(bam1_t*,std::string&,std::vector<bam1_t*> &);
	void qname_stats(std::string&,std::vector<bam1_t*> &);
	void flag_inspector(bam1_t*);
	void histo_Pdistance();
	void update_histodata(bam1_t*); 
	void processHeader(bam_hdr_t*); 
	void processReads(samFile* , bam_hdr_t* , bam1_t*);
	void output();
	void run();
    
};

#endif /* runner_hpp */


/*  
DOMANDE:
-problema numero 1: se il mate di una coppia per un errore non è stato segnato la coppia è invalida ma non sapremo mai quale (forse con grandi numeri è inutile) 
-ma il pos e il pos del mate non hanno di mezzo la sequenza del frammento stesso? non dovrebbe essere pos dell'ultima base mappata del primo frammento e inizio del mate?

DOMANDE A YING:
-le mapped/e non sono delle primary o in generale di tutte? //risposto samtools
-mean insert è tra le coppie buone o tutti record in generale?
-i duplicati sono dei singoli read non delle paia, o l'averli calcolati con pairtools cambia la loro situa?
-error rate: è tra tutte le letture o solo le mapped? e si calcola con NM o con cigar?(NM non è sempre presente) //risposto samtools
-distanza minima sensata tra i read per p(s)
*/