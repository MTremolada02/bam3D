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

void Runner::qname_stats(std::vector<bam1_t*> &group, int group_size){
	std::vector <bam1_t*> r1_side, r2_side;
	std::vector <bam1_t*> *chim = nullptr;
	std::vector <bam1_t*> *other= nullptr; //potrebbe essere solo un bam e non un vettore di bam? tanto è solo uno !!!!!!!!!!!!!!!
	bool rescued =false;
	bool R1=false;
	bool R2=false;

	bool cis=false;
	bool facing=false;
	bool distance=false;

	uint8_t mapped_count1=0;
	uint8_t mapped_count2=0;
	Maptype a=qnameStats.type1=Maptype::N;
	Maptype b=qnameStats.type2=Maptype::N; //così se esiste solo uno dei due read (anche se con flag paired attiva), lo mette a null;
//////walked e rescued. Una molecola candidata WW viene rescued: se ha una sola ligazione reale, ma appare come walk per effetti geometrici/tecnici.
	for(int i=0; i<group_size; ++i){
		if (group[i]->core.flag & BAM_FSECONDARY) continue; //rimangono solo i supplementary e i primary mappati
    	if (group[i]->core.flag & BAM_FUNMAP) continue;

    	if (group[i]->core.flag & BAM_FREAD1) r1_side.push_back(group[i]);
    	if (group[i]->core.flag & BAM_FREAD2) r2_side.push_back(group[i]);
	}

	if(r1_side.size()>=2 || r2_side.size()>=2){
		if((r1_side.size()==2 && r2_side.size()==1) || (r2_side.size()==2 && r1_side.size()==1)){ 
	if	(r1_side.size() == 2) {
				chim = &r1_side;
				other = &r2_side;
				R1=true;
			}else if (r2_side.size() == 2){
				chim = &r2_side;
				other = &r1_side;
				R2=true;
			}
		} else {
		++qnameStats.WW;
		return;
		}

		if (Alignstarts((*chim)[0]) > Alignstarts((*chim)[1])) { //se partono dello stesso posto più è lunga la parte non allineata più mi avvicino alla ligazione
			std::swap((*chim)[0], (*chim)[1]);///////non tanto chiaro
		}//ora chim[0] è l'outer, [1]è l'inner e l'altro è automaticamente l'other
		bool rev_i = (*chim)[1]->core.flag & BAM_FREVERSE;
		bool rev_o = (*other)[0]->core.flag & BAM_FREVERSE;
		//controlli 
		if(((*chim)[1])->core.tid==((*other)[0])->core.tid) {cis=true;} 
		if((!rev_i &&  rev_o && (*chim)[1]->core.pos <= (*other)[0]->core.pos) || ( rev_i && !rev_o && (*other)[0]->core.pos <= (*chim)[1]->core.pos)) {facing=true;}
		if(llabs((*chim)[1]->core.pos - (*other)[0]->core.pos) <= 2000) {distance=true;}////di default 2000 pb

		if(cis && facing && distance) {
			rescued=true;
		}else{
			++qnameStats.WW; 
			return;  
		}
	}

//////pair stats
	for(int i=0; i<group_size; ++i){
		if(group[i]->core.flag & BAM_FUNMAP) {continue;}
		if(group[i]->core.flag & BAM_FSUPPLEMENTARY) {continue;} //elimino i supplementari (devo farlo?)
		if(!(group[i]->core.flag & BAM_FPAIRED)) {break;}// se non è una coppiaè inutile fare la statistica

		if(!(group[i]->core.flag & BAM_FSECONDARY) && group[i]->core.flag & BAM_FDUP) { //per ora uso i duplicati marcati nel bamfile, si può fare un mappa in cui si salvano tutte le coppie e si controllano realmente i duplicati
			++qnameStats.DD;
			return;
		}
		if(group[i]->core.flag & BAM_FREAD1) {
			if(group[i]->core.flag & BAM_FSECONDARY) {mapped_count1++;}
			else{mapped_count1++;} //dovrebbe essere il primary
		}
		if(group[i]->core.flag & BAM_FREAD1) {++mapped_count1;}
		if(group[i]->core.flag & BAM_FREAD2) {++mapped_count2;} //conta sia secondary che primary
	}

	a=(mapped_count1==0 ? Maptype::N :(mapped_count1==1 ? Maptype::U : Maptype::M));
	b=(mapped_count2==0 ? Maptype::N :(mapped_count2==1 ? Maptype::U : Maptype::M));

	if(rescued && R1) {a=Maptype::R;}
	if(rescued && R2) {b=Maptype::R;}

	if(a==Maptype::U && b==Maptype::R) {
		++qnameStats.UR;
		return;
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
	

void Runner::processReads(bam_hdr_t* bamHdr , std::vector<std::vector<bam1_t*>> &vectorbox, int boxes_filled, std::vector<int> &group_counter) {
	bool qname_sorted =(std::string(bamHdr->text, bamHdr->l_text).find("SO:queryname") != std::string::npos); // perchè la funzione string.find() ritorna npos;!!!!!! posso anche passarla ocme argomento
	uint64_t av_counter=0;
	
	uint64_t mismatched_bases=0;
	uint64_t total_base=0;

	std::unordered_map <uint64_t,uint64_t> global_dist_count; //per ora la tengo così
	std::map <uint32_t,std::unordered_map<uint64_t,uint64_t>> chrom_dist_count;
	uint32_t chrom=0;
	uint64_t dist=0;
	
	for(int i=0;i<boxes_filled;++i){ //sono sempre lo stesso numero aparte l'ultimo giro che può essere incompleto
		for(int k=0; k<group_counter[i]; ++k){   //j
			pairStats.good_read1=false;
			pairStats.good_read2=false;
			++readStats.readN;
			////FLAG STATS;		
			if (!(vectorbox[i][k]->core.flag & BAM_FUNMAP) && vectorbox[i][k]->core.qual==0) {++readStats.mapQ0;}
			flag_inspector(vectorbox[i][k]);

			if (pairStats.good_read1 && vectorbox[i][k]->core.tid == vectorbox[i][k]->core.mtid) {
				++pairStats.sameCr;

				if(((vectorbox[i][k]->core.flag & BAM_FREVERSE) != (vectorbox[i][k]->core.flag & BAM_FMREVERSE)) && std::abs((long double)vectorbox[i][k]->core.isize)>0){//così hanno sempre orientamenti opposti
					++av_counter;			
					readStats.mean_insert = update_mean_tlen(readStats.mean_insert, av_counter, vectorbox[i][k]);   
					//	readStats.quadratic_mean=update_quadratic_mean_tlen(readStats.mean_insert,av_counter, bamdata);
 					readStats.quadratic_mean=update_quadratic_mean_tlen(readStats.quadratic_mean,av_counter, vectorbox[i][k]);

				}
				if(userInput.hist_global){ //HISTO_GLOBAL_DATA
					dist=llabs(vectorbox[i][k]->core.pos - vectorbox[i][k]->core.mpos); // dovrebbero essere degli uint64_t quindi non serve forzare il double ne arrotondare
					++global_dist_count[dist]; 
				}
			
				if(userInput.hist_by_chrom){	 //HISTO_CHROM_DATA
					chrom=vectorbox[i][k]->core.tid;
					dist=llabs(vectorbox[i][k]->core.pos - vectorbox[i][k]->core.mpos);
					++chrom_dist_count[chrom][dist];
				}
			}

			if (pairStats.good_read1 || pairStats.good_read2) { 
				uint8_t* nm_ptr = bam_aux_get(vectorbox[i][k], "NM");//diff tra a read e il riferimento
    			uint64_t nm = nm_ptr ? bam_aux2i(nm_ptr) : 0;

    			mismatched_bases += nm;  
    			uint64_t aligned = bam_cigar2rlen(vectorbox[i][k]->core.n_cigar, bam_get_cigar(vectorbox[i][k])); //bam_cigar2rlen(int n_cigar, const uint32_t *cigar):This function returns the sum of the lengths of the M, I, S, = and X operations in @p cigar (these are the operations that "consume" query bases
    			total_base += aligned;
			} 
		}
		if(qname_sorted){ 
			qname_stats(vectorbox[i], group_counter[i]); //fa le statistiche sui qname per ogni gruppo
		}
	}	

	if(userInput.hist_global){histo_global_distance(global_dist_count);}
	if (userInput.hist_by_chrom){histo_chrom_distance(chrom_dist_count);}

	if ((readStats.readN % 2) == 0) {
		pairStats.pairN=readStats.readN/2;
	} else {pairStats.pairN += (readStats.readN % 2) * 2 >= 2 ? 1 : 0;} //arromtonda all'intero più vicino 
	readStats.error_rate=error_rate(mismatched_bases,total_base);
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
	std::cout<<"UU"<<":"<<"MM"<<":"<<"NN"<<":"<<"UM"<<":"<<"UN"<<":"<<"NM"<<"\t"<<qnameStats.UU<<":"<<qnameStats.MM<<":"<<qnameStats.NN<<":"<<qnameStats.MU<<":"<<qnameStats.NU<<":"<<qnameStats.NM<<std::endl;
	std::cout<<"DD"<<":"<<"WW"<<":"<<"UR"<<":"<<"RU"<<":"<<"RN"<<":"<<"RM"<<"\t"<<qnameStats.DD<<":"<<qnameStats.WW<<":"<<qnameStats.UR<<":"<<qnameStats.RU<<":"<<qnameStats.NR<<":"<<qnameStats.MR<<std::endl;
	
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
	

}
 

int Runner::data_vector(std::vector<std::vector<bam1_t*>> &vectorbox, bam1_t* &bridge_read, bool &bridge, std::vector<int> &group_counter, samFile *fp_in, bam_hdr_t *bamHdr){
	bool qname_sorted =(std::string(bamHdr->text, bamHdr->l_text).find("SO:queryname") != std::string::npos); // perchè la funzione string.find() ritorna npos;
	std::string qname;
	std::string current_qname;
	long unsigned int j=0;//indice dei vettori di record interni a vectorbox:ultima posizione piena zero compreso (tipo di dato per vectorbox[i-1].size()<=j+1 che evidentemente è un long unsigned int)
	int boxes_filled=0;
	if(!bridge){bool first = true;}//fist reading

	if(!qname_sorted){
		for (int i=0;i<10;++i){
			if(sam_read1(fp_in, bamHdr, vectorbox[i][0]) <0){
				return i;
			}else{
				boxes_filled=i+1;
				group_counter[i]=1;//ultima posizione piena del vettore è la zeresima
			}
		}
		return boxes_filled;
	}else{
		for(int i=0; i<11; ++i){ //ora vector box è nella posizione [i][0] piena del record i-esimo corrente (first e bridge a parte)

			if (bridge){//primo record dopo ogni cambio di vettore (quindi i,j=1->[0][0])
				bam_copy1(vectorbox[0][0],bridge_read);
				qname = bam_get_qname(vectorbox[0][0]);
				group_counter[0]=1;
			    boxes_filled=1;
				bridge=false; 
				continue;
			}else{
				if(sam_read1(fp_in, bamHdr, vectorbox[i][0]) <0){
					return i+1;
				}
			}

			if (first){ // per il primo in assoluto quindi [0][0] del primissimo thread di tutto il file
				qname = bam_get_qname(vectorbox[i][0]);
				group_counter[0]=1;
				boxes_filled=1;
				first = false;
				continue;
			}
				
			current_qname=bam_get_qname(vectorbox[i][0]);

			if(strcmp(current_qname.c_str(), qname.c_str())==0){//devo salvarlo nella posizione i precedente ma j successiva al j attuale del sottovettore del vectorbox
				if(vectorbox[i-1].size()<=j+1){//controllo se esiste già uno spazio da sovrascrivere nella pozione i-1
					vectorbox[i-1].push_back(bam_init1());
					bam_copy1(vectorbox[i-1][j+1], vectorbox[i][0]);
				}else {//lo spazio c'è già:saremo al secondo giro. potrei usare swap ma così uso puntatori, posso solo sovrascrivere
					bam_copy1(vectorbox[i-1][j+1], vectorbox[i][0]);
				}
				++j;//al secondo giro j vale 1
				group_counter[i-1]=j+1;
				--i;
			} else {//quindi qname diverso
				if(i==10){ //ultimo giro
					bridge=true;
					bam_copy1(bridge_read,vectorbox[10][0]);
				} else {
					qname=bam_get_qname(vectorbox[i][0]);
					j=0; //primo record del nuovo gruppo con j=0
					group_counter[i]=1;
					boxes_filled=i+1;
				}
			}
		}
		return boxes_filled;
	}
} 

void Runner::run() {
	std::size_t numFiles = userInput.inFiles.size();
	lg.verbose("Processing " + std::to_string(numFiles) + " files");
	
	for (uint32_t i = 0; i < numFiles; ++i) {

		std::string file = userInput.file('r', i);
		std::string ext = getFileExt(file);
		
		samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); //open bam file
		bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
		
		htsThreadPool tpool_read = {NULL, 0};
		tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
		if (tpool_read.pool) {	hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
		} else { lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");}

		std::vector<std::vector<bam1_t*>> vectorbox; //vettori di vettori, i quali conterranno almeno un record al loro interno
		std::vector<int> group_counter(10); //vettore di indici
		bool bridge=false;
		bam1_t *bridge_read= bam_init1(); 
		
		for (int i=0; i<11; ++i){vectorbox[i].push_back(bam_init1());} // lo riempio e rimarrà con le stract vuote (11 perchè l'ultimo mi serve per il bridge) 
		//ne ho bisogno uno per ogni thread
//decompression e reading thread
		
		while(true){
			int boxes_filled=data_vector(vectorbox, bridge_read, bridge, group_counter,fp_in, bamHdr); //riempie vectorbox con i record del file in input
			processReads(bamHdr, vectorbox, boxes_filled, group_counter); //statistiche tipo samtools
			if (!bridge) break; // Se non c'è un bridge, il file è finito
		}

		output();
		//devo farne uno per thread

	for (long unsigned int i = 0; i < vectorbox.size(); ++i) { //warning: comparison of integer expressions of different signedness: ‘int’ and ‘std::vector<std::vector<bam1_t*> >::size_type’ {aka ‘long unsigned int
    	for (bam1_t* b : vectorbox[i]) {
      	  bam_destroy1(b);
    	}
    	vectorbox[i].clear(); // Opzionale: azzera la size per sicurezza
	}
		bam_hdr_destroy(bamHdr);
		bam_destroy1(bridge_read);
		sam_close(fp_in);
		if (tpool_read.pool)
		hts_tpool_destroy(tpool_read.pool);
	}
}


