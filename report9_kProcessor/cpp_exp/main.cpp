#include <iostream>


/*

fasta_file="/home/mabuelanin/Desktop/kexpression_experiment/symbolic/data/drtamer_data/reps_unitigs_SRR11015356_k75.fa"
names_file = fasta_file + ".names"

kSize = 25
Q = 30
mode = 1 # kmers
chunk_size = int(1e3)

kf = kp.kDataFrameMQF(kSize, Q, mode)

print("Indexing ...")

ckf = kp.index(kf, {"mode" : mode}, fasta_file , chunk_size, names_file)

print("Serializing the index ...")

ckf.save("idx_k25_SRR11015356_k75.unitigs")

*/

#include <iostream>
#include "kDataFrame.hpp"
#include <string>
#include <vector>
#include <stdint.h>
#include <gqf.h>
#include "algorithms.hpp"
#include "Utils/kmer.h"
#include "Utils/utils.hpp"
#include <cmath>
#include <map>

using namespace std;


int main(){

    
    string fasta_file = "/home/mabuelanin/Desktop/kexpression_experiment/symbolic/data/drtamer_data/reps_unitigs_SRR11015356_k75.fa";
    string names_file = fasta_file + ".names";

    int kSize = 25;
    int Q = 27;
    int mode = 1 ;
    int chunk_size = 100000;

    map<string,int> params = {{"mode", mode}};

    kDataFrame * kf = new kDataFrameMQF(kSize, Q, mode);

    // index(kDataFrame *frame, std::map<std::string, int> parse_params, string filename, int chunks, string names_fileName);
    cerr << "Indexing ... " << endl;
    colored_kDataFrame * ckf = kProcessor::index(kf, params, fasta_file , chunk_size, names_file);
    cerr << "Serializing ... " << endl;
    ckf->save("../idx_reps_unitigs_SRR11015356_k75");


    return 0;
}