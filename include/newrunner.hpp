#ifndef newrunner_hpp
#define newrunner_hpp

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
    const char* current_qname() const noexcept;

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

    std::vector<bam1_t*> slots; //o meglio allocarla con new?
	std::size_t size_wanted;
    std::size_t used = 0; //il primo libero
    int hiwater_data = 0;
	bool file_end=false; 
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
    void data_vector(Bam_record_vector &,samFile *,bam_hdr_t *);
	void data_vector(Bam_record_vector &, bam1_t *,bool &, samFile *, bam_hdr_t *);
	void processReads(Bam_record_vector &);
	void output();
	void run();
    
};

#endif /* newrunner_hpp */

/*
class bam_record_vector {
private:
	bam1_t** reads;      // specie di array di puntatori a *bam1_t (scatola da cui partono tutti i puntatori a bam1_t)
    uint8_t capacity;        // Quanti puntatori bam1_t allocati ho
    uint8_t size;            // a quante read sono arrivata per questo group (counter) (è la prima libera perchè sempre a i+1)

public: 
    qnameGroup(){ //costruttore di default
		size=0;
		capacity=4;
		reads = (bam1_t**)malloc(capacity * sizeof(bam1_t*)); //tengo da parte della memoria(grande tanto quando un bam1_t) per i puntatori , ALLOCATI CON MALLOC SONO CONTIGUI!
        for (uint8_t i = 0; i < capacity; ++i) {//alloca le 4 strutture vuote
            reads[i] = bam_init1();
        }
    }

	~qnameGroup() {
        if (reads) {
            for (int i = 0; i < capacity; ++i) {
                bam_destroy1(reads[i]);
            }
            free(reads); // free a block allocated by malloc(def)
        }
    }
	
	void expand(){ //push_back se size==capacity
		uint8_t old_capacity = capacity;
		capacity *= 2; // Raddoppiamo lo spazio
		reads = (bam1_t**)realloc(reads, capacity * sizeof(bam1_t*)); //se c'è spazio lo aumenta oltre se stesso (se è finita la memoria reads punerà a NULL e crasha tutto)
		for (uint8_t i = old_capacity; i < capacity; ++i) {//aggiungiamo costruttori vuoti in coda
     		reads[i] = bam_init1();
    	}

	}

	void clear(){size=0;} //contatore uguale a 0 (capacity rimane uguale)

	inline uint8_t get_size(){return size;}
	
	inline uint8_t get_capacity(){return capacity;}
	
	void add_size() {
		if(capacity==size){
			expand();
		}
		++size;
	}
	
	inline bam1_t* get_read(uint8_t index) const {
        return reads[index];
    }

	void add_read(bam1_t* source) {
    	if (size == capacity) {
       	 expand(); 
   	 	}
   		bam_copy1(reads[size], source);
    	size++;
	}

	const char* get_current_qname() {  // se la scatola è vuota mi ridarà null e crasha
		return bam_get_qname(get_read[size-1]);
	}

};
*/