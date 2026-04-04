#ifndef def_hpp
#define def_hpp

#include <htslib/sam.h>
#include <map>
#include <vector>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

//AGGIUNTA PER GRAFICI
struct TlenClassRow {
    std::string run;
    std::string length_class;
    std::string ref_name;
    uint64_t contig_len = 0;
    uint64_t abs_tlen = 0;
    double log10_abs_tlen = 0.0;
    uint8_t mapq1 = 0;
    uint8_t mapq2 = 0;
};

struct ReadErrorRow {
    std::string run;
    std::string ref_name;
    uint64_t contig_len = 0;
    uint64_t read_len1 = 0;
    uint64_t read_len2 = 0;
    uint64_t nm1 = 0;
    uint64_t nm2 = 0;
    double nm_rate1 = 0.0;
    double nm_rate2 = 0.0;
    uint8_t mapq1 = 0;
    uint8_t mapq2 = 0;
};


//////////////////////////////////////////////////////////////////////////////

struct UserInputBam3D : UserInput { // additional input
	bool hist_none        =false;
	bool hist_global      = false;
	bool hist_by_chrom    = false;
	bool single_read_stats = false;
	bool pair_read_stats = false;

	uint8_t decompression_threads = 4;
	//uint8_t histogram_mode = hist_none;

};

struct Graph {
	std::unordered_map<uint32_t, uint64_t> graph_intra_count;
	std::map<std::pair<uint32_t,uint32_t>, uint64_t> graph_inter_count;

	std::unordered_map<uint64_t, uint64_t> Ps_binned_dist_count;
	std::unordered_map<uint64_t, uint64_t> ff_binned_dist_count; //forse un po' grande uint_64 
	std::unordered_map<uint64_t, uint64_t> fr_binned_dist_count;
	std::unordered_map<uint64_t, uint64_t> rf_binned_dist_count;
	std::unordered_map<uint64_t, uint64_t> rr_binned_dist_count; 

	std::unordered_map<uint32_t, uint64_t> isize_hist;
	double log_bin_factor = std::pow(10.0, 1.0 / 10.0); // 10 bin per decade: pow(10.0, 1.0 / 10.0)
	double inv_log_bin_factor = 1.0 / std::log(log_bin_factor);
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

    uint64_t one_side  = 0;
    uint64_t two_side_mapped  = 0;
    uint64_t duplicated   = 0;
    uint64_t dupl=0;
    uint64_t UNmapped     = 0;
    uint64_t cis       = 0;  // cis

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

uint64_t dbg_unresolved = 0;
uint64_t dbg_not_2plus1 = 0;
uint64_t dbg_outer_bad = 0;
uint64_t dbg_other_no_index = 0;
uint64_t dbg_other_type_N = 0;
uint64_t dbg_other_type_M = 0;
uint64_t dbg_ignore_inner = 0;
uint64_t dbg_geom_pass = 0;
uint64_t dbg_geom_fail = 0;
uint64_t dbg_gap_large = 0;
uint64_t dbg_outer_bad_otherN = 0;
uint64_t dbg_outer_bad_otherM = 0;
uint64_t dbg_outer_bad_otherU = 0;
uint64_t fallback = 0;
uint64_t dbg_not3 =0;
uint64_t dbg_unmapped =0;
uint64_t dbg_outer_noindex = 0;
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
	Graph graph;

private:

private:

	std::vector<TlenClassRow> tlen_class_rows;
	std::vector<ReadErrorRow> read_error_rows;

    std::unordered_map<uint64_t,uint64_t> global_dist_count;
    std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>> chrom_dist_count;
    
public:

	void collect_tlen_for_ecdf_violin(const bam1_t*, const bam1_t*, bam_hdr_t*, const std::string&);
	void collect_read_error_data(const bam1_t*, const bam1_t*, bam_hdr_t*, const std::string&);
	uint64_t Runner::get_nm_mismatches(const bam1_t*);
	void write_tlen_class_section(std::ofstream&);
	void write_read_error_section(std::ofstream&);

    void loadInput(UserInputBam3D userInput);
	void write_section_header(std::ofstream&, const std::string&,const std::string&);
	void write_all_stats_file(const std::string&, bam_hdr_t*);
	void write_binned_map(std::ofstream&, const std::string&, const std::unordered_map<uint64_t, uint64_t>&);
	void write_pair_types_section(std::ofstream&);
	void write_reference_graph(std::ofstream&, bam_hdr_t*);
	void update_log_binned_distance(uint64_t, std::unordered_map<uint64_t, uint64_t>&);
	void update_pair_plots_from_records(const bam1_t*, const bam1_t*);
	void update_reference_graph(const bam1_t*, const bam1_t*);
	uint64_t cigar_mapped_bases(const bam1_t*);
	double percentage(std::size_t, double);
	long double update_mean_tlen(long double,uint64_t, bam1_t*);
	long double update_quadratic_mean_tlen(long double,uint64_t, bam1_t*);
//	void estimate_insert_stats();
	double error_rate(uint64_t,uint64_t);
	uint16_t Alignstarts(const bam1_t*);
	void qname_stats(Bam_record_vector &, bam_hdr_t*);
	//void qname_stats(Bam_record_vector &);
	void flag_inspector(bam1_t*);
	void histo_global_distance(std::unordered_map<uint64_t, uint64_t>&);
	void histo_chrom_distance(std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>>&); 
    void data_vector(Bam_record_vector &,samFile *,bam_hdr_t *);
	void data_vector(Bam_record_vector &, bam1_t *,bool &, samFile *, bam_hdr_t *);
	void processReads(Bam_record_vector &, bam_hdr_t*);
	//void processReads(Bam_record_vector &);
	int Alignend(const bam1_t* b);	
	int inter_align_gap_on_query(const bam1_t* left_seg, const bam1_t* right_seg);
	void output();
	void write_analytics_file(const std::string&);
	void run();
    
};



#endif /* def_hpp */
