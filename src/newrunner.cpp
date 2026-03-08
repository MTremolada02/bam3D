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
#include <unordered_map>
#include <map>
#include <algorithm>

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "functions.h"
#include "global.h"
#include "newrunner.hpp"

void Runner::loadInput(UserInputBam3D userInput) {
    this->userInput = userInput;
}

//se le read partono uguali ma finiscono diverse lo considero un walk in ogni caso, anche pairtools dovrebbe fare così, se è un miglioramento si può pensare ad un'implementazione
uint16_t Runner::Alignstarts(const bam1_t* b){//legge il cigar e riporta le basi segnate nel read a sinistra prima delle basi mappate
	const uint32_t* cigar = bam_get_cigar(b);
	uint16_t bases=0;

	for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        int8_t op  = bam_cigar_op(cigar[i]);//cigar operator character
        int8_t len = bam_cigar_oplen(cigar[i]);//quante basi per lettera (es.50S)

		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {break;}//inizia l'allineamento con M,X,=
		if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {bases+=len;}//sommo quante basi S e I ci sono prima che la read si allinei al riferimento
	}
	return bases;
}

void Runner::qname_stats(Bam_record_vector &group) {//sono contanta di essere riuscita ad adattarlo! ora vedremo
	std::size_t begin=0;//fist_read
	std::size_t end=0;
	const char* qname;

	std::vector <int> r1_side;
	r1_side.reserve(5);//inventato
	std::vector <int> r2_side;
	r2_side.reserve(5);

	std::size_t inner=0;
	std::size_t outer=0;
	std::size_t other=0;

	bool rescued=false;
	bool R1_chim=false;
	bool R2_chim=false;

	bool cis=false;
	bool facing=false;
	bool distance=false;

	uint8_t mapped_count1=0;
	uint8_t mapped_count2=0;

	Maptype a=qnameStats.type1=Maptype::N;;
	Maptype b=qnameStats.type2=Maptype::N;;

	while(begin<group.size()){ 
		qname=bam_get_qname(group[begin]);
		end=begin+1;
		
    	while(end < group.size() && strcmp(bam_get_qname(group[end]), qname) == 0){++end;}////////////////////group [begin,end) perchè end=begin+1 che è il primo diverso

	/*	if(begin >= group.size()){
    std::cout << "BAD BEGIN " << begin << " size=" << group.size() << std::endl;
    break;
}
		std::cout << "\nGROUP begin=" << begin
          << " end=" << end
          << " size=" << (end - begin)
          << " qnamei=" << bam_get_qname(group[begin])
		  << " qnamef=" << bam_get_qname(group[end-1])<<std::endl;*/

		r1_side.clear();
		r2_side.clear();

	    inner=0;
	    outer=0;
	    other=0;

		rescued=false;
		R1_chim=false;
		R2_chim=false;

		cis=false;
		facing=false;
		distance=false;

		mapped_count1=0;
		mapped_count2=0;

		a=qnameStats.type1=Maptype::N;
		b=qnameStats.type2=Maptype::N;

//////walked e rescued. Una molecola candidata WW viene rescued: se ha una sola ligazione reale, ma appare come walk per effetti geometrici/tecnici.
		for(std::size_t j=begin; j<end; ++j){
			if (group[j]->core.flag & BAM_FSECONDARY) continue; //rimangono solo i supplementary e i primary mappati
    		if (group[j]->core.flag & BAM_FUNMAP) continue;

    		if (group[j]->core.flag & BAM_FREAD1) {
				r1_side.push_back(j);//unico dai
				std::cout<<"r1_side"<<j<<std::endl;
			}
			
    		if (group[j]->core.flag & BAM_FREAD2) { 
				r2_side.push_back(j);
				std::cout<<"r2_side"<<j<<std::endl;
			}
		}

		if(r1_side.size()>=2 || r2_side.size()>=2){
			if((r1_side.size()==2 && r2_side.size()==1) || (r2_side.size()==2 && r1_side.size()==1)){ 
				if	(r1_side.size() == 2) {
					inner = r1_side.at(0);
					outer = r1_side.at(1);
					other = r2_side.at(0);
					R1_chim=true;

				}else if (r2_side.size() == 2){
					inner = r2_side.at(0);
					outer = r2_side.at(1);
					other = r1_side.at(0);
					R2_chim=true;
				}
			} else {
				std::cout<<"ww1"<<std::endl;
				++qnameStats.WW;//return;
				begin=end;
				continue; 
			}

			auto start_i = Alignstarts(group[inner]);
			auto start_ou = Alignstarts(group[outer]);

			if (start_i > start_ou){std::swap(inner, outer);}//se partono dello stesso posto più è lunga la parte non allineata più mi avvicino alla ligazione
//RIGUARDARE SENSO BIOLOGICO

			bool rev_i = group[inner]->core.flag & BAM_FREVERSE;
			bool rev_ot = group[other]->core.flag & BAM_FREVERSE;
			//controlli 
			if((group[inner])->core.tid==(group[other])->core.tid) {cis=true;} 
			if((!rev_i &&  rev_ot && group[inner]->core.pos <= group[other]->core.pos) || ( rev_i && !rev_ot && group[other]->core.pos <= group[inner]->core.pos)) {facing=true;}
			if(llabs(group[inner]->core.pos - group[other]->core.pos) <= 2000) {distance=true;}////di default 2000 pb

			if(cis && facing && distance) {
				rescued=true;
			}else{
				++qnameStats.WW; 
				std::cout<<"WW2"<<std::endl;
				begin=end;
				continue;  //rerturn;
			}
		}

//////pair stats
		//devo rianalizzare lo stesso gruppo nel caso non fosse stato un WW
		for(std::size_t i=begin; i<end; ++i){
			if(group[i]->core.flag & BAM_FUNMAP) {continue;}
			if(group[i]->core.flag & BAM_FSUPPLEMENTARY) {continue;} //elimino i supplementari (devo farlo?)
			if(!(group[i]->core.flag & BAM_FPAIRED)) {break;}// se non è una coppia è inutile fare la statistica

			if(!(group[i]->core.flag & BAM_FSECONDARY) && group[i]->core.flag & BAM_FDUP) { //per ora uso i duplicati marcati nel bamfile, si può fare un mappa in cui si salvano tutte le coppie e si controllano realmente i duplicati
				++qnameStats.DD;
				begin=end;
				continue;//return;
			}
			if(group[i]->core.flag & BAM_FREAD1) {
				if(group[i]->core.flag & BAM_FSECONDARY) {mapped_count1++;
				}else{mapped_count1++;} //dovrebbe essere il primary
			}
			//if(group[tmp]->core.flag & BAM_FREAD1) {++mapped_count1;}
			if(group[i]->core.flag & BAM_FREAD2) {++mapped_count2;} //conta sia secondary che primary della R2
		}

		a=(mapped_count1==0 ? Maptype::N :(mapped_count1==1 ? Maptype::U : Maptype::M));
		b=(mapped_count2==0 ? Maptype::N :(mapped_count2==1 ? Maptype::U : Maptype::M));

		if(rescued && R1_chim) {a=Maptype::R;}
		if(rescued && R2_chim) {b=Maptype::R;}

		if(a==Maptype::U && b==Maptype::R) {
			++qnameStats.UR;
			begin=end;
			continue;//return;
		} //unico di cui importa l'ordine

		if (static_cast<uint8_t>(a) < static_cast<uint8_t>(b)) { std::swap(a, b);} //raggruppo UM e MU perchè a sarà sempre M rispetto a U

		if(a==Maptype::U && b==Maptype::U) ++qnameStats.UU;
		if(a==Maptype::M && b==Maptype::M) ++qnameStats.MM;
		if(a==Maptype::N && b==Maptype::N) ++qnameStats.NN;

		if(a==Maptype::M && b==Maptype::U) ++qnameStats.MU;
		if(a==Maptype::U && b==Maptype::N) ++qnameStats.NU;
		if(a==Maptype::M && b==Maptype::N) ++qnameStats.NM;

		if(a==Maptype::R && b==Maptype::U) ++qnameStats.RU;
		if(a==Maptype::R && b==Maptype::M) ++qnameStats.MR;
		if(a==Maptype::R && b==Maptype::N) ++qnameStats.NR;

		begin=end;
		
	}
}
	
long double Runner::update_mean_tlen(long double prev_mean,std::uint64_t k, bam1_t* bamdata){  //<x>
    long double xk = std::abs((long double)bamdata->core.isize);  // TEN dLel record
    return (xk / k) + ((k - 1) / (long double)k) * prev_mean;									
}

long double Runner::update_quadratic_mean_tlen(long double prev_mean,std::uint64_t k, bam1_t* bamdata){ //<x^2> FORSE SBAGLIATA
	long double xk = std::abs((long double)bamdata->core.isize);  // TLEN del record
	long double xk2 = xk * xk;
 //   return (pow(xk,2) / k) + ((k - 1) / (long double)k) * pow(prev_mean,2);
	return prev_mean + (xk2 - prev_mean) / (long double)k;
}

double Runner::error_rate(uint64_t mismatched_bases,uint64_t total_base){
	return (total_base==0) ? 0.0 : (long double)mismatched_bases/total_base;
}

void Runner::histo_global_distance (std::unordered_map<uint64_t,uint64_t>& global_dist_count){
	std::fstream myfile;
	myfile.open("Pair_by_global_distance.txt",std::ios::out); //agginugere il path (i vecchi dati vengono cancellati e sovrascritti) 

	if(!myfile.is_open()){
		std::cout<<"pair_by_global_distance not open"<<std::endl;
		return;
	}

	myfile << "distance" << "\t" << "count" << "\n";
	for (auto i = global_dist_count.begin(); i != global_dist_count.end(); ++i) {
		myfile << i->first << "\t" << i->second<< "\n";
	}
	myfile.close();  
}

void Runner::histo_chrom_distance (std::map<uint32_t,std::unordered_map<uint64_t,uint64_t>>& chrom_dist_count) { 
	std::fstream myfile;
	myfile.open("Pair_chromosome_by_distance.txt",std::ios::out); //aggiungere il path (i vecchi dati vengono cancellati e sovrascritti)

	if(!myfile.is_open()){
		std::cout<<"pair_chromosome_by_distance not open"<<std::endl;
		return;
	}
	
	myfile << "\n# chromosome" <<"\t"<<"distance"<<"\t"<<"counter"<< "\n";


	for (auto i = chrom_dist_count.begin(); i != chrom_dist_count.end(); ++i) {

    	uint32_t chrom = i->first;
    	const auto& dist_map = i->second;

    	for (auto j = dist_map.begin(); j != dist_map.end(); ++j) {
       	 myfile <<i->first<< "\t" << j->first << "\t" << j->second << "\n";
    	}
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
		return;
	} else {readStats.primary++;}

	//sulle coppie, faccio tutte e due insieme così ho già filtrato dall calcolo le coppie supplementary e secondary
	
	if(flag & BAM_FPAIRED || flag & BAM_FPROPER_PAIR) {
		++pairStats.pairN;

		if(flag & BAM_FDUP) {++pairStats.duplicated;} //WARNING! se c'è solo una read1 ma segnata come duplicato così non la conto nelle rad (basta toglier l'else if)
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
	

void Runner::processReads(Bam_record_vector &vectorbox) {
	if(userInput.single_read_stats || userInput.hist_global || userInput.hist_by_chrom){
		uint64_t av_counter=0;
		
		uint64_t mismatched_bases=0;
		uint64_t total_base=0;

		uint32_t chrom=0;
		uint64_t dist=0;

		if(userInput.single_read_stats){ 
			for(int i=0;i<vectorbox.size();++i){
					pairStats.good_read1=false;
					pairStats.good_read2=false;
					++readStats.readN;
					////FLAG STATS;		
					if (!(vectorbox[i]->core.flag & BAM_FUNMAP) && vectorbox[i]->core.qual==0) {++readStats.mapQ0;}
					flag_inspector(vectorbox[i]);
 
					if (pairStats.good_read1 && vectorbox[i]->core.tid == vectorbox[i]->core.mtid) {
						++pairStats.sameCr;

						if(((vectorbox[i]->core.flag & BAM_FREVERSE) != (vectorbox[i]->core.flag & BAM_FMREVERSE)) && std::abs((long double)vectorbox[i]->core.isize)>0){//così hanno sempre orientamenti opposti
							++av_counter;			
							readStats.mean_insert = update_mean_tlen(readStats.mean_insert, av_counter, vectorbox[i]);   
							//	readStats.quadratic_mean=update_quadratic_mean_tlen(readStats.mean_insert,av_counter, bamdata);
							readStats.quadratic_mean=update_quadratic_mean_tlen(readStats.quadratic_mean,av_counter, vectorbox[i]);

						}
						if(userInput.hist_global){ //HISTO_GLOBAL_DATA
							dist=llabs(vectorbox[i]->core.pos - vectorbox[i]->core.mpos); // dovrebbero essere degli uint64_t quindi non serve forzare il double ne arrotondare
							++global_dist_count[dist]; 
						}
					
						if(userInput.hist_by_chrom){	 //HISTO_CHROM_DATA
							chrom=vectorbox[i]->core.tid;
							dist=llabs(vectorbox[i]->core.pos - vectorbox[i]->core.mpos);
							++chrom_dist_count[chrom][dist];
						}
					}

					if (pairStats.good_read1 || pairStats.good_read2) { 
						uint8_t* nm_ptr = bam_aux_get(vectorbox[i], "NM");//diff tra a read e il riferimento
						uint64_t nm = nm_ptr ? bam_aux2i(nm_ptr) : 0;

						mismatched_bases += nm;  
						uint64_t aligned = bam_cigar2rlen(vectorbox[i]->core.n_cigar, bam_get_cigar(vectorbox[i])); //bam_cigar2rlen(int n_cigar, const uint32_t *cigar):This function returns the sum of the lengths of the M, I, S, = and X operations in @p cigar (these are the operations that "consume" query bases
						total_base += aligned;
					} 
			}

			if ((readStats.readN % 2) == 0) {
				pairStats.pairN=readStats.readN/2;
			} else {pairStats.pairN += (readStats.readN % 2) * 2 >= 2 ? 1 : 0;} //arromtonda all'intero più vicino 
			readStats.error_rate=error_rate(mismatched_bases,total_base);
		}
	}
	if(userInput.pair_read_stats){ 
		qname_stats(vectorbox);
	}
}


void Runner::output(){
	
		std::cout<<"Qc_fail:"<<readStats.qc_fail<<std::endl;
		std::cout<<"Tot_record:"<<readStats.readN<<std::endl;
		std::cout<<"Non_primary:"<<readStats.secondary+readStats.supplementary<<std::endl; 
		std::cout<<"Reads_mapped:"<<readStats.readN-readStats.unmapped<<std::endl;       
		std::cout<<"%mapped:"<< 100-((readStats.unmapped*100)/(long double)readStats.readN)<<"%"<<std::endl;    
		std::cout<<"Proper_pairs:"<<((pairStats.proper_pairs*100)/(long double)readStats.readN)<<"%"<<std::endl;
		std::cout<<"%MapQ0:"<< ((readStats.mapQ0*100)/(long double)readStats.readN)<<"%"<<std::endl;
		std::cout<<"MapQ0:"<<readStats.mapQ0<<std::endl;
		std::cout<<"Pairs:"<<pairStats.pairN<<std::endl; 
		std::cout<<"Read1:"<<pairStats.read1<<std::endl; 
		std::cout<<"Read2:"<<pairStats.read2<<std::endl;
		std::cout<<"unmapped:"<<readStats.unmapped<<std::endl;
		std::cout<<"%One_sided:"<< ((pairStats.UMone_sided*100)/(long double)pairStats.pairN)<<"%"<<std::endl;  
		std::cout<<"%Two_sided:"<< ((pairStats.UMtwo_sided*100)/(long double)pairStats.pairN)<<"%"<<std::endl; 
		std::cout<<"%Duplicated:"<< ((pairStats.duplicated*100)/(long double)pairStats.pairN)<<"%"<<std::endl;
		std::cout<<"%CIS:"<< ((pairStats.sameCr*100)/(long double)pairStats.pairN)<<"%"<<std::endl;  
		std::cout<<"mean_insert:"<<readStats.mean_insert<<std::endl; 
		std::cout<<"insert SD:"<<sqrt(readStats.quadratic_mean-pow(readStats.mean_insert,2))<<std::endl; //rad(<x^2>-<x>^2)
		std::cout<<"error_rate:"<<readStats.error_rate<<std::endl;
		std::cout<<"|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<std::endl;
		std::cout<<"Samtools_Stats_Tot_record(Tot_record-Non_Primary):"<<readStats.readN-(readStats.secondary+readStats.supplementary)<<std::endl;
		std::cout<<"%mapped:"<< 100-((readStats.unmapped*100)/(long double)readStats.readN-(readStats.secondary+readStats.supplementary))<<"%"<<std::endl; 
		std::cout<<"Reads_mapped:"<<(readStats.readN-(readStats.secondary+readStats.supplementary))-readStats.unmapped<<std::endl;
		std::cout<<"%MapQ0:"<< ((readStats.mapQ0*100)/(long double)(readStats.readN-(readStats.secondary+readStats.supplementary)))<<"%"<<std::endl;
		std::cout<<"MapQ0:"<<readStats.mapQ0<<std::endl;
		std::cout<<"Pairs:"<<(readStats.readN-(readStats.secondary+readStats.supplementary))/2<<std::endl;
		std::cout<<"Pairs_tot_veri?"<<(readStats.readN)/2<<std::endl;
		std::cout<<"%One_sided:"<< ((pairStats.UMone_sided*100)/(long double)readStats.unmapped)<<"%"<<std::endl;  
		std::cout<<"%Two_sided:"<< ((pairStats.UMtwo_sided*100)/(long double)readStats.unmapped)<<"%"<<std::endl; 
		std::cout<<"%Duplicated:"<< pairStats.duplicated<<std::endl;
		std::cout<<"%Duplicated:"<< ((pairStats.duplicated*100)/(long double)((readStats.readN-(readStats.secondary+readStats.supplementary))/2))<<"%"<<std::endl;
		std::cout<<"%CIS:"<< ((pairStats.sameCr*100)/(long double)((readStats.readN-(readStats.secondary+readStats.supplementary))/2))<<"%"<<std::endl;  
		
		std::cout<<"UU"<<":"<<"MM"<<":"<<"NN"<<":"<<"UM"<<":"<<"UN"<<":"<<"NM"<<"\t"<<qnameStats.UU<<":"<<qnameStats.MM<<":"<<qnameStats.NN<<":"<<qnameStats.MU<<":"<<qnameStats.NU<<":"<<qnameStats.NM<<std::endl;
		std::cout<<"DD"<<":"<<"WW"<<":"<<"UR"<<":"<<"RU"<<":"<<"RN"<<":"<<"RM"<<"\t"<<qnameStats.DD<<":"<<qnameStats.WW<<":"<<qnameStats.UR<<":"<<qnameStats.RU<<":"<<qnameStats.NR<<":"<<qnameStats.MR<<std::endl;
}

void Runner::data_vector(Bam_record_vector &vectorbox,samFile *fp_in,bam_hdr_t *bamHdr){
	vectorbox.clear();
	for (int i=0; i<vectorbox.get_size_wanted();++i){ 
		if(!vectorbox.add_record(fp_in,bamHdr))
			break;
	}
}

void Runner::data_vector(Bam_record_vector &vectorbox, bam1_t *bridge_read,bool &first, samFile *fp_in, bam_hdr_t *bamHdr){
	const char* qname;
	const char* current_qname;
	bool bridge=true;
	vectorbox.clear();

	for (int i=0;i<vectorbox.get_size_wanted();++i){
		if(bridge){
			if(first){//:( come faccio
				vectorbox.add_record(fp_in,bamHdr);
				first =false;
			}else{
				vectorbox.push_back(bridge_read);
			}
			bridge=false;
			continue;
		}
		vectorbox.add_record(fp_in,bamHdr);
	}

	qname=bam_get_qname(vectorbox[vectorbox.size()-1]);
	for(;;){
		if(sam_read1(fp_in, bamHdr, bridge_read)>=0){ 
			current_qname=bam_get_qname(bridge_read);
			if(strcmp(current_qname, qname) != 0){
				break;}
			vectorbox.push_back(bridge_read); 
		}else{break;}                 
	}
}



void Runner::run() {

	std::size_t numFiles = userInput.inFiles.size();
	lg.verbose("Processing " + std::to_string(numFiles) + " files");
	
	for (uint32_t i = 0; i < numFiles; ++i) {

		global_dist_count.clear();//svuotare le mappe prima di ogni file o si mescoleranno (se non è l'obiettivo)s
		chrom_dist_count.clear();

		std::string file = userInput.file('r', i);
		std::string ext = getFileExt(file);
		
		samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); 
		if (!fp_in) {std::cout<<"hts_open has failed"<<std::endl;}
		bam_hdr_t *bamHdr = sam_hdr_read(fp_in); 
		if (!bamHdr) {std::cout<<"sam_hdr_read has failed"<<std::endl;}
		
		htsThreadPool tpool_read = {NULL, 0};
		tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
		if (tpool_read.pool) {	hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
		} else { lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");}

		bool qname_sorted =(std::string(bamHdr->text, bamHdr->l_text).find("SO:queryname") != std::string::npos); // perchè la funzione string.find() ritorna npos;
		if (!userInput.single_read_stats && !userInput.pair_read_stats) {//default
			userInput.single_read_stats = true;
    		if (qname_sorted) {
        	userInput.pair_read_stats = true;
    		} else { std::cout<<"Warning: input BAM file is not qname sorted, pair read statistics will not be computed."<<std::endl;}
		}
		if(userInput.pair_read_stats && !qname_sorted){
			userInput.pair_read_stats=false;
			std::cout<<"Error: to compute pair read statistics the input BAM file must be qname sorted."<<std::endl;
			exit(1);
		}
		std::size_t j=10;//set real capacity
		bool first=true;
 
		Bam_record_vector records_vector(j); 
		bam1_t *bridge_read=bam_init1(); //È UN PUNTATORE

		while(!(records_vector.is_file_end())){ 
			if(userInput.pair_read_stats){ 
				data_vector(records_vector, bridge_read,first, fp_in, bamHdr);
			}else if(userInput.single_read_stats){
				data_vector(records_vector,fp_in,bamHdr);
			}
			processReads(records_vector);
                }

		if(userInput.hist_global){histo_global_distance(global_dist_count);}

		if(userInput.hist_by_chrom){histo_chrom_distance(chrom_dist_count);}
    		

		output();
		
		bam_hdr_destroy(bamHdr);
		bam_destroy1(bridge_read);
		sam_close(fp_in);
		if (tpool_read.pool)
		hts_tpool_destroy(tpool_read.pool);
	}
}


//////////////////////////////////////////////////////////////////////////////////////////class functions definition

Bam_record_vector::Bam_record_vector(std::size_t initial_capacity)
    : used(0), hiwater_data(0), size_wanted(0), file_end(false)
{
    slots.reserve(initial_capacity);
	size_wanted=initial_capacity;
    for (std::size_t i = 0; i < initial_capacity; ++i){ 
        slots.push_back(bam_init1());
	}
}

Bam_record_vector::~Bam_record_vector() {
   for (auto* b : slots) bam_destroy1(b);
}

Bam_record_vector::Bam_record_vector(Bam_record_vector&& other) noexcept //move contructor
    : slots(std::move(other.slots)),
    used(other.used),
    hiwater_data(other.hiwater_data),
	size_wanted(other.size_wanted),
    file_end(other.file_end)
{
	other.slots.clear();
    other.used= 0;
    other.hiwater_data= 0;
	other.size_wanted = 0;
    other.file_end = false;
}
Bam_record_vector& Bam_record_vector::operator=(Bam_record_vector&& other) noexcept{ //omve assignment operator
    if (this == &other) return *this;
    	for (auto* b : slots) bam_destroy1(b);

    	slots= std::move(other.slots);
    	used= other.used;
   		hiwater_data= other.hiwater_data;
		size_wanted = other.size_wanted;
    	file_end = other.file_end;

		other.slots.clear();
    	other.used= 0;
    	other.hiwater_data= 0;
		other.size_wanted = 0;
        other.file_end = false;
    return *this;
}

void Bam_record_vector::clear() noexcept { used=0;}
std::size_t Bam_record_vector::size() const noexcept { return used; }
std::size_t Bam_record_vector::capacity() const noexcept { return slots.size(); }
std::size_t Bam_record_vector::get_size_wanted() const noexcept {return size_wanted;}
bool Bam_record_vector::is_file_end() const noexcept {return file_end;}

bam1_t* Bam_record_vector::operator[](std::size_t i) noexcept { return slots.at(i); }
const bam1_t* Bam_record_vector::operator[](std::size_t i) const noexcept { return slots[i]; }

bool Bam_record_vector::add_record(samFile *fp_in,bam_hdr_t *bamHdr){
	if(used==slots.size()){expand(slots.empty() ? 10 : slots.size() * 2);}
	if(sam_read1(fp_in, bamHdr, slots[used])>=0){
		++used;
		return true;
	} else {
		file_end=true;
		return false;
	}
}

bam1_t* Bam_record_vector::push_back(const bam1_t* src) { // da sorgente al primo slot libero del vectorbox.
    if (used== slots.size())
        expand(slots.empty() ? 10 : slots.size() * 2);//se p empty per qualche motivo a cosa lo metto aiuto ahahha

    bam1_t* dst = slots[used];

	if (src->l_data > dst->m_data) {
    	bam1_t* new_dst = bam_dup1(src);
	 	if (!new_dst)
        throw std::bad_alloc();

        bam_destroy1(dst);
        slots[used] = new_dst;
        dst = new_dst;
    } else {
        bam_copy1(dst, src);
    }
    ++used;

    hiwater_data= std::max<int>(hiwater_data, dst->l_data);

	return dst;
}

void Bam_record_vector::expand(std::size_t new_capacity) {
    if (new_capacity <= slots.size()) return;
    slots.reserve(new_capacity);
    while (slots.size() < new_capacity) {
        bam1_t* b = bam_init1();
        if (!b) throw std::bad_alloc();

		if (hiwater_data > 0) {
   			b->data = (uint8_t*)malloc(hiwater_data);
   	 		if (!b->data) {
       			bam_destroy1(b);
       			throw std::bad_alloc();
   			}
   			b->m_data = hiwater_data;
    		b->l_data = 0;
		}
		slots.push_back(b);
    }
}
