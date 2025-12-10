#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "functions.h"
#include "global.h"
#include "runner.hpp"

void Runner::loadInput(UserInputBam3D userInput) {
    this->userInput = userInput;
}

void Runner::qname_group(bam1_t *bamdata,std::string &qname,std::vector<bam1_t*> &group){

	std::string current_qname=bam_get_qname(bamdata);

	if(qname==current_qname){ //stesso qname
		group.push_back(bam_dup1(bamdata)); //duplica il bam1_t puntato da bamdata e ne salva il puntatore nel vettore

	} else {
		qname_stats(qname,group);
			
		for (bam1_t* rec : group) { //per ogni puntatore in group
    		bam_destroy1(rec); // libera il bam1_t a cui il puntatore rec si riferisce
		}
	group.clear();

	qname=current_qname;
	group.push_back(bam_dup1(bamdata)); 
	} 
}

void Runner::qname_stats(std::string &qname, std::vector<bam1_t*> &group){
//fare statistiche
}

long double Runner::update_mean_tlen(long double prev_mean,std::uint64_t k, bam1_t* bamdata){  //<x>
    long double xk = std::abs((long double)bamdata->core.isize);  // TLEN del record
    return (xk / k) + ((k - 1) / (long double)k) * prev_mean;
}

long double Runner::update_quadratic_mean_tlen(long double prev_mean,std::uint64_t k, bam1_t* bamdata){ //<x^2> FORSE SBAGLIATA
	long double xk = std::abs((long double)bamdata->core.isize);  // TLEN del record
    return (pow(xk,2) / k) + ((k - 1) / (long double)k) * pow(prev_mean,2);
}

double Runner::error_rate(uint64_t mismatched_bases,uint64_t total_base){
	return (total_base==0) ? 0.0 : (long double)mismatched_bases/total_base;
}

void Runner::histo_Pdistance (){
	std::fstream myfile;
	myfile.open("Pair_by_distance.txt",std::ios::out); //agginugere il path (i vecchi dati vengono cancellati e sovrascritti)

	if(!myfile.is_open()){
		std::cout<<"pair_by_distance not open"<<std::endl;
		return;
	}

	for(long double i=0.0L; i<graphStats.edge.size()-1; ++i) {
		myfile<<(sqrt(graphStats.edge[i]*graphStats.edge[i+1]))<<"\t"<<graphStats.counts[i]<<std::endl;
	}

	myfile.close();  
}

void Runner::flag_inspector (bam1_t* bamdata) {
	uint16_t flag= bamdata-> core.flag;

	//sui singoli record

	if (flag & BAM_FQCFAIL){++readStats.qc_fail;} //return?
	if (flag & BAM_FUNMAP) {++readStats.unmapped;} //così i mapped sono di tutti come in samtools ma ho la richiesta una volta sola e non dentro e fuori dall'else
	if (flag & BAM_FPROPER_PAIR) {++pairStats.proper_pairs;}

	if(flag & BAM_FSUPPLEMENTARY) {
		++readStats.supplementary;
		return;  // o si tolgono i return e mettere le cose di coppia nell'else come fa samtools flagstats
	} else if (flag & BAM_FSECONDARY) {
		++readStats.secondary;
		//std::cout<<"secondary=="<<stats.secondary<<std::endl;
		return;
	} else {readStats.primary++;}

	//sulle coppie, faccio tutte e due insieme così ho già filtrato dall calcolo le coppie supplementary e secondary
	
	if(flag & BAM_FPAIRED || flag & BAM_FPROPER_PAIR) {
		++pairStats.pairN;

		if(flag & BAM_FDUP) {++pairStats.duplicated;}
		else if (flag & BAM_FREAD1) { // così ne prendo solo una e non due non so se ha senso
			++pairStats.read1;

			if ((flag & BAM_FUNMAP && !(flag & BAM_FMUNMAP))^(flag & BAM_FMUNMAP && !(flag & BAM_FUNMAP))) {++pairStats.UMone_sided;} // statistica fatta sul singolo se no è doppia
			else if (flag & BAM_FUNMAP && (flag & BAM_FMUNMAP)) {++pairStats.UMtwo_sided;}
	
			if(!(flag & BAM_FUNMAP) && !(flag & BAM_FMUNMAP)) {
				pairStats.good_read1=true;
				++pairStats.good_pairs;
			}
		}
		else if (flag & BAM_FREAD2) {
			++pairStats.read2;
			if(!(flag & BAM_FUNMAP) && !(flag & BAM_FMUNMAP)) {
				pairStats.good_read2=true;
			}

		}	
	} 
}
	
void Runner::processHeader(bam_hdr_t *bamHdr) {
	long double l_max=0;
	long double l_min=1.0; //QUANTE??
	int n_bins=0;
	int bins_per_decade = 20;
	double step=1.0/bins_per_decade; // sempre supponendo 20 bin per decade
	
	for(int i=0; i<bamHdr->n_targets; ++i) {
		l_max=std::max<long double>(l_max,bamHdr->target_len[i]);
	}

    n_bins=(int)bins_per_decade*(log10(l_max)-log10(l_min)); //se è vero che ci sono 20 bin per decade:(10^) ->lo step in scala logaritmica è 1/20=0.05
	graphStats.edge.assign(n_bins + 1, 0.0L);
    graphStats.counts.assign(n_bins, 0L);

	for(int i=0; i<=n_bins; ++i){
		graphStats.edge[i] = pow(10.0L, log10(l_min) + i * step); //edge[i]=l_min​⋅10^(i⋅Δ) ogni 20 bins si ha fatto una decade
	}
}

void Runner::update_histodata(bam1_t* bamdata){
    long long sep = std::llabs((long long)bamdata->core.pos - (long long)bamdata->core.mpos);

    // fuori range? salta
    if (sep < graphStats.edge.front() || sep > graphStats.edge.back()) return;

    // trova il bin
    int bin = std::upper_bound(graphStats.edge.begin(), graphStats.edge.end(), (long double)sep)- graphStats.edge.begin() - 1; // upper_bound mi ritorna il bordo dell'edge subito più grande del mio valore di sep es115 edge{100,110,150...}, upper_bound mi ritorna 150, upperbound[]-begin[0]= indice della casella dell'edge max:[3] ma [2-3) corrisponde al bin 2 quindi il -1

    if (bin >= 0 && bin < (int)graphStats.counts.size()) {graphStats.counts[bin]++;} //controllo più aumento del count del bin corrisponente (così "esistono i bin")
}

void Runner::processReads(samFile *fp_in, bam_hdr_t *bamHdr, bam1_t *bamdata) { 
	bool qname_sorted =(std::string(bamHdr->text, bamHdr->l_text).find("SO:queryname") != std::string::npos); // perchè la funzione string.find() ritorna npos;
	std::vector<bam1_t*> group;
	std::string qname;
	bool first = true;

	uint64_t av_counter=0;
	
	uint64_t mismatched_bases=0;
	uint64_t total_base=0;

	processHeader(bamHdr);
	
	while(sam_read1(fp_in, bamHdr, bamdata)>0) {
		pairStats.good_read1=false;
		pairStats.good_read2=false;
		++readStats.readN;
		
		if (bamdata->core.qual==0) {++readStats.mapQ0;}
		flag_inspector(bamdata);

		if (pairStats.good_read1 && bamdata->core.tid == bamdata->core.mtid) {
			++av_counter;
			++pairStats.sameCr;

			readStats.mean_insert = update_mean_tlen(readStats.mean_insert, av_counter, bamdata);    //agigungere controllo sull'orientamento delle letture di r1 e r2!
			readStats.quadratic_mean=update_quadratic_mean_tlen(readStats.mean_insert,av_counter, bamdata);
			update_histodata(bamdata);
		}

		if (pairStats.good_read1 || pairStats.good_read2) { 
			uint8_t* nm_ptr = bam_aux_get(bamdata, "NM");
    		uint64_t nm = nm_ptr ? bam_aux2i(nm_ptr) : 0;

    		mismatched_bases += nm;  
    		uint64_t aligned = bam_cigar2rlen(bamdata->core.n_cigar, bam_get_cigar(bamdata)); //bam_cigar2rlen(int n_cigar, const uint32_t *cigar):This function returns the sum of the lengths of the M, I, S, = and X operations in @p cigar (these are the operations that "consume" query bases
    		total_base += aligned;
		} 

		if(qname_sorted){
 			if (first) {  // primo record
            	qname = bam_get_qname(bamdata);
            	group.push_back(bam_dup1(bamdata));
            	first = false;
        	} else {
            qname_group(bamdata, qname, group);
        	}
		}

	}

	if (qname_sorted && !group.empty()) {
    	qname_stats(qname,group);
    	for (bam1_t *rec : group) bam_destroy1(rec);
    		group.clear();
	}


	pairStats.pairN=pairStats.pairN/2; //conto pairN solo con i read1? if(stats.pair1==stats.pair2) {lo divido} else {errore qualcosa non va nel file :/} //controllo sulle paia?
	readStats.error_rate=error_rate(mismatched_bases,total_base);
}

void Runner::output(){ //si può scrivere meglio con una funzione di perfentu
	std::cout<<"Qc_fail:"<<readStats.qc_fail<<std::endl;
	std::cout<<"Tot_record:"<<readStats.readN<<std::endl;
	std::cout<<"Non_primary:"<<readStats.secondary+readStats.supplementary<<std::endl; 
	std::cout<<"Reads_mapped:"<<readStats.readN-readStats.unmapped<<std::endl;       
	std::cout<<"%mapped:"<< 100-((readStats.unmapped*100)/(long double)readStats.readN)<<"%"<<std::endl;    //non mi escono percentuali se non forzo il double
	std::cout<<"Proper_pairs:"<<((pairStats.proper_pairs*100)/(long double)readStats.readN)<<"%"<<std::endl;
	std::cout<<"%MapQ0:"<< ((readStats.mapQ0*100)/(long double)readStats.readN)<<"%"<<std::endl;
	std::cout<<"Pairs:"<<pairStats.pairN<<std::endl; // quanto debole come condizione il /2? controllo migliore con messaggio di errore per una frazione senza resto 0?
	std::cout<<"Read1:"<<pairStats.read1<<std::endl; // non esplicitamente chiesti ma di controllo sul /2 
	std::cout<<"Read2:"<<pairStats.read2<<std::endl;// non esplicitamente chiesti ma di controllo sul /2 
	std::cout<<"%One_sided:"<< ((pairStats.UMone_sided*100)/(long double)pairStats.pairN)<<"%"<<std::endl;  
	std::cout<<"%Two_sided:"<< ((pairStats.UMtwo_sided*100)/(long double)pairStats.pairN)<<"%"<<std::endl; 
	std::cout<<"%Duplicated:"<< ((pairStats.duplicated*100)/(long double)pairStats.pairN)<<"%"<<std::endl;
	std::cout<<"%CIS:"<< ((pairStats.sameCr*100)/(long double)pairStats.pairN)<<"%"<<std::endl;  
	std::cout<<"mean_insert:"<<readStats.mean_insert<<std::endl; 
	std::cout<<"insert SD:"<<sqrt(readStats.quadratic_mean-pow(readStats.mean_insert,2))<<std::endl; //rad(<x^2>-<x>^2)
	std::cout<<"error_rate:"<<readStats.error_rate<<std::endl;
}

void Runner::run() {
	
	std::size_t numFiles = userInput.inFiles.size();
	lg.verbose("Processing " + std::to_string(numFiles) + " files");
	
	for (uint32_t i = 0; i < numFiles; ++i) {
		
		//uint64_t readCounter = 0;
		std::string file = userInput.file('r', i);
		std::string ext = getFileExt(file);
		
		samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); //open bam file
		bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
		bam1_t *bamdata = bam_init1(); //initialize an alignment

		//std::vector<bam1_t*> bamBatch;
		
		htsThreadPool tpool_read = {NULL, 0};
		tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
		if (tpool_read.pool) {	hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
		} else { lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");}

		processReads(fp_in,bamHdr, bamdata);
		output();
		histo_Pdistance();

		bam_hdr_destroy(bamHdr);
		bam_destroy1(bamdata);
		sam_close(fp_in);
		if (tpool_read.pool)
			hts_tpool_destroy(tpool_read.pool);
	}
}




/*
	corpo del process read
	for (bam1_t *bamdata : record) {// per processare un vettore con i record vanno salvati tutti, non basta quello che si era fatto l'altra volta

		flag_inspector(record);

		uint32_t len = bamdata->core.l_qseq; // length of the read.
		uint8_t *seq = bam_get_seq(bamdata); // seq string
		std::unique_ptr<std::string> inSequenceQuality;
		
		std::unique_ptr<std::string> inSequence = std::make_unique<std::string>();
		inSequence->resize(len);
		for(uint32_t i=0; i<len; ++i)
			inSequence->at(i) = seq_nt16_str[bam_seqi(seq,i)]; //gets nucleotide id and converts them into IUPAC id.
		
		uint8_t* qual = bam_get_qual(bamdata);
		if (qual && (len > 0) && qual[0] != 0xFF) {
			inSequenceQuality = std::make_unique<std::string>(len, '\0');
			for (uint32_t i = 0; i < len; ++i)
				(*inSequenceQuality)[i] = static_cast<char>(qual[i] + 33);
		} else {
			inSequenceQuality = std::make_unique<std::string>(len, '!'); // No per-base qualities; synthesize minimal qualities
		}
		//std::cout<<*inSequence<<"\t"<<*inSequenceQuality<<std::endl;
	} */

		/*while(sam_read1(fp_in,bamHdr,bamdata) > 0) { //immaginiamo siano tutte sorted per qname
			
			bamBatch.push_back(bamdata);
			++readCounter;
			if (readCounter % 50 == 0) {
				processReads(bamBatch);
				bamBatch.clear();
			}
		}*/