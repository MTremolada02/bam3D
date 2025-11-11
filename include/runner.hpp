#ifndef runner_hpp
#define runner_hpp

struct UserInputBam3D : UserInput { // additional input
	uint8_t decompression_threads = 4;
};

class Runner {
    
    UserInputBam3D userInput;
        
public:
    
    void loadInput(UserInputBam3D userInput);
    void run();
    
};

#endif /* runner_hpp */
