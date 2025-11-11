#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include "functions.h"
#include "global.h"
#include "runner.hpp"

void Runner::loadInput(UserInputBam3D userInput) {
    this->userInput = userInput;
}

void Runner::run() {
	
	std::size_t numFiles = userInput.inFiles.size();
	lg.verbose("Processing " + std::to_string(numFiles) + " files");
	
	for (uint32_t i = 0; i < numFiles; ++i) {
		
		std::string file = userInput.file('r', i);
		std::string ext = getFileExt(file);
		
		//Sequences* readBatch = new Sequences;
		samFile *fp_in = hts_open(userInput.file('r', i).c_str(),"r"); //open bam file
		bam_hdr_t *bamHdr = sam_hdr_read(fp_in); //read header
		bam1_t *bamdata = bam_init1(); //initialize an alignment
		
		htsThreadPool tpool_read = {NULL, 0};
		tpool_read.pool = hts_tpool_init(userInput.decompression_threads);
		if (tpool_read.pool)
			hts_set_opt(fp_in, HTS_OPT_THREAD_POOL, &tpool_read);
		else
			lg.verbose("Failed to generate decompression threadpool with " + std::to_string(userInput.decompression_threads) + " threads. Continuing single-threaded");
		
		while(sam_read1(fp_in,bamHdr,bamdata) > 0) {
			
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
			std::cout<<*inSequence<<"\t"<<*inSequenceQuality<<std::endl;
		}
		bam_destroy1(bamdata);
		sam_close(fp_in);
		if (tpool_read.pool)
			hts_tpool_destroy(tpool_read.pool);
	}
}
