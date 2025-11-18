#ifndef runner_hpp
#define runner_hpp

#include <htslib/sam.h>

struct UserInputBam3D : UserInput { // additional input
	uint8_t decompression_threads = 4;
};

class Runner {
    
    UserInputBam3D userInput;
	uint64_t readN = 0;
        
public:
    
    void loadInput(UserInputBam3D userInput);
	void processReads(std::vector<bam1_t*> &readBatch);
    void run();
    
};

#endif /* runner_hpp */
