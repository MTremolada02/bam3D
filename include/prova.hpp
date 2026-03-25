#ifndef prova_hpp
#define prova_hpp

#include <htslib/sam.h>
#include <map>
#include <vector>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

struct UserInputBam3D : UserInput { // additional input
	bool hist_none        =false;
	bool hist_global      = false;
	bool hist_by_chrom    = false;
	bool single_read_stats = false;
	bool pair_read_stats = false;

	uint8_t decompression_threads = 4;
	//uint8_t histogram_mode = hist_none;

};

struct ReadStats {
    uint64_t readN  = 0;
    uint64_t qc_fail = 0;
    uint64_t unmapped = 0;
    uint64_t secondary = 0;
    uint64_t supplementary = 0;
    uint64_t primary = 0;
    uint64_t mapQ0 = 0;

    uint64_t mismatched_bases=0;
    uint64_t total_mapped_base=0;
	std::size_t av_counter=0;
	std::size_t avf_counter=0;
//	std::map<uint32_t, uint64_t> insert_hist; //bin , counter
//	const uint32_t bin_size = 100;   // 100 bp per bin
//	uint64_t sd_insert=0;
//	uint64_t mean_insert=0;

	long double mean_insert=0;
	long double quadratic_mean=0;
	long double mean_insert_filtr=0;
	long double quadratic_mean_filtr=0;

	double error_rate=0;
};

struct PairStats {
	std::size_t n_group;
	std::size_t counter1=0;
	std::size_t counter2=0;
    uint64_t pairN        = 0;
    uint64_t proper_pairs = 0;
    uint64_t good_pairs   = 0;
	bool good_read1=false; //dovrebbe essere locale
	bool good_read2=false; //dovrebbe essere locale

    uint64_t UMone_sided  = 0;
    uint64_t UMtwo_sided  = 0;
    uint64_t duplicated   = 0;
    uint64_t UNmapped     = 0;
    uint64_t sameCr       = 0;  // cis

    uint64_t read1 = 0;
    uint64_t read2 = 0;
};

enum class Maptype :uint8_t {N=0, U=1, M=2, R=3};

struct QnameStats { 
	uint64_t LOST;
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

class Bam_record_vector {
public:
    explicit Bam_record_vector(std::size_t initial_capacity);
    ~Bam_record_vector();

    // non-copyable
    Bam_record_vector(const Bam_record_vector&) = delete;
    Bam_record_vector& operator=(const Bam_record_vector&) = delete;

    // movable
    Bam_record_vector(Bam_record_vector&& other) noexcept;
    Bam_record_vector& operator=(Bam_record_vector&& other) noexcept;

    // main API
    bam1_t* push_back(const bam1_t* src);

    // container-like helpers
    void clear() noexcept;
    std::size_t size() const noexcept;
    std::size_t capacity() const noexcept;
	std::size_t get_size_wanted() const noexcept;
	bool is_file_end() const noexcept;

	bool add_record(samFile *fp_in,bam_hdr_t *bamHdr);

    bam1_t* operator[](std::size_t i) noexcept;
    const bam1_t* operator[](std::size_t i) const noexcept;

private:
    void expand(std::size_t new_capacity);

    std::vector<bam1_t*> slots;
	std::size_t size_wanted=0;
    std::size_t used = 0; 
    int hiwater_data = 0;
	bool file_end=false; 
};

class Runner {

    UserInputBam3D userInput;
	ReadStats readStats;
	PairStats pairStats;
	QnameStats qnameStats;

private:

private:

    std::unordered_map<uint64_t,uint64_t> global_dist_count;
    std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>> chrom_dist_count;
    
public:
    void loadInput(UserInputBam3D userInput);
	uint64_t cigar_mapped_bases(const bam1_t*);
	long double update_mean_tlen(long double,uint64_t, bam1_t*);
	long double update_quadratic_mean_tlen(long double,uint64_t, bam1_t*);
//	void estimate_insert_stats();
	double error_rate(uint64_t,uint64_t);
	uint16_t Alignstarts(const bam1_t*);
	void qname_stats(Bam_record_vector &);
	void flag_inspector(bam1_t*);
	void histo_global_distance(std::unordered_map<uint64_t, uint64_t>&);
	void histo_chrom_distance(std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>>&); 
        void data_vector(Bam_record_vector &,samFile *,bam_hdr_t *);
	void data_vector(Bam_record_vector &, bam1_t *,bool &, samFile *, bam_hdr_t *);
	void processReads(Bam_record_vector &);
	void output();
	void run();
    
};

#endif /* prova_hpp */
