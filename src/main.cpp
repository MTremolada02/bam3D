#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <getopt.h>
#include <unordered_map>

#include "log.h"
#include "global.h"
#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "threadpool.h"

#include <runner.hpp>
#include <main.hpp>

std::string version = "0.0.1";
std::string toolName = "bam3D";

std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int verbose_flag;
int maxThreads = 0;
int cmd_flag;
int tabular_flag;
UserInputBam3D userInput;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;
std::vector<Log> logs;

void printHelp() {
	printf("%s [command]\n", toolName.c_str());
	printf("\nOptions:\n");
	printf("-v --version software version.\n");
	printf("--cmd print $0 to stdout.\n");
	exit(0);
}

std::string getArgs(char* optarg, unsigned int argc, char **argv) {
    
    std::string cmd;
    bool record = false;

    for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
        
        if (optarg != argv[arg_counter] && !record) {
            continue;
        }else{
            record = true;
            if(optarg != argv[arg_counter]){
                cmd += ' ';
                cmd += argv[arg_counter];
            }
        }
    }
    
    return cmd;
    
}

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    
    std::string cmd;
    
    bool isPipe = false; // to check if input is from pipe
    
    if (argc == 1) { // tool name with no arguments
        printf("%s [command]\n-h for additional help.\n", toolName.c_str());
        exit(0);
    }
    
    static struct option long_options[] = { // struct mapping long options
		{"input-bam", required_argument, 0, 'b'},
		
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"hist-global", no_argument, 0, 'g'},
        {"hist-chrom", no_argument, 0, 'm'},
        
        {0, 0, 0, 0}
    };
    
//    const static std::unordered_map<std::string,int> tools{
//        {"tool1",1},
//        {"tool2",2}
//    };
    
    const static std::unordered_map<std::string,int> modes{
        {"parseBam",1},
        {"mode2",2}
    };
    
    userInput.mode = modes.count(argv[1]) ? modes.at(argv[1]) : 0;
    if (userInput.mode == 0) {
        fprintf(stderr, "Unrecognized mode (%s).\n", argv[1]);
        return EXIT_FAILURE;
    }
        
    switch (userInput.mode) {
        default:
        case 0: { // mode1
            
            while (arguments) { // loop through argv
                
                int option_index = 0;
                
                c = getopt_long(argc, argv, "-:b:vhgm",
                                long_options, &option_index);
                
                if (c == -1) { // exit the loop if run out of options
                    break;
                }
                
                switch (c) {
                    case ':': // handle options without arguments
                        switch (optopt) { // the command line option last matched
                            case 'f':
                                break;
                                
                            default:
                                fprintf(stderr, "option -%c is missing a required argument.\n", optopt);
                                return EXIT_FAILURE;
                        }
                        break;
                    default: // handle positional arguments
                        
                        //                switch (tools.count(optarg) ? tools.at(optarg) : 0) {
                        //                    case 1:
                        //                        cmd = "tool1/build/bin/mytool1" + getArgs(optarg, argc, argv);;
                        //                        arguments = false;
                        //                        break;
                        //                    case 2:
                        //                        cmd = "tool2/build/bin/tool2" + getArgs(optarg, argc, argv);;
                        //                        arguments = false;
                        //                        break;
                        //                }
                    case 0: // case for long options without short options
                        
                        //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                        //                  splitLength = atoi(optarg);
                        break;
					case 'b': // read input bam files
						if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
							userInput.pipeType = 'r'; // pipe input is a sequence
						}else{ // input is a regular file
							optind--;
							for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
								ifFileExists(argv[optind]);
								userInput.inFiles.push_back(argv[optind]);
							}
						}
						break;
                    case 'g': //histogram options
                    case 'm':
                        if(c=='g'){userInput.hist_global=true;}
                        if(c=='m'){userInput.hist_by_chrom=true;}
                        break;

                    case '?': // unrecognized option
                        if (optopt != 0)
                            fprintf(stderr, "Unrecognized option: %c\n", optopt);
                        else
                            fprintf(stderr, "Unrecognized option: %s\n", argv[optind-1]);
                        printHelp();
                        break;
                    case 'v': // software version
                        printf("%s v%s\n", toolName.c_str(), version.c_str());
                        printf("Giulio Formenti giulio.formenti@gmail.com\n");
						printf("Margherita Tremolada margheritatremolada@gmail.com\n");
                        exit(0);
                    case 'h': // help
						printHelp();
                }
                
                if    (argc == 2 || // handle various cases in which the output should include summary stats
                       (argc == 3 && pos_op == 2) ||
                       (argc == 4 && pos_op == 3)) {
                }
            }
        }
    }
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
//    std::cout<<"Invoking: "<<cmd<<std::endl;
//    std::system(cmd.c_str());
    
    Runner runner;
    threadPool.init(maxThreads); // initialize threadpool
    maxMem = (userInput.maxMem == 0 ? get_mem_total(3) * 0.9 : userInput.maxMem); // set memory limit
    runner.loadInput(userInput); // load user input
    lg.verbose("User input loaded");
    runner.run(); // run algorithms
    threadPool.join(); // join threads
    exit(EXIT_SUCCESS);
}
