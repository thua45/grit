//
// Created by Tinghua Huang on 19/12/30.
//

#include <getopt.h>
#include <map>
#include <string>
#include <iostream>

#include "grit.h"

//#include "gerea.h"

char const short_options[] = "hm:i:b:n:z:s:t:p:o:d:";
struct option long_options[] =
        {
            {"help", 0, NULL, 'h'},
            {"motif", 1, NULL, 'm'},
            {"homoseq", 1, NULL, 'i'},
            {"background", 1, NULL, 'b'},
            {"sampling_times", 1, NULL, 'n'},
            {"sampling_size", 1, NULL, 'z'},
            {"lratio", 1, NULL, 's'},
            {"pvalue", 1, NULL, 't'},
            {"score", 1, NULL, 'p'},
            {"output", 1, NULL, 'o'},
            {"dir", 1, NULL, 'd'},
            {0, 0, 0, 0}
        };

int main(int argc, char *argv[]) {

    if (argc == 1)
    {
        printf("grit -m motif -i homoseq -b bgseq [-n sampling_times] [-z sampling_size] [-s rs] [-t pvalue] [-p pscore] -o output\n");
        exit(1);

    }

    string session_name = "my_session";
    unsigned long permutation_n = 1000;
    string mat_file = "motif.txt";
    string bgseq_file = "bgseq.txt";
    string seq_file = "homoseq.txt";
    string output_file = "output.txt";
    unsigned int sampling_size = 200;
    unsigned int sampling_times = 10;
    unsigned int lratio = 1;  //1 for averare lratio, 2 for maximum lratio
    double p_score = 0.0;
    double threshold = 0.05;
    string output_dir = "./";

    int c;
    while((c=getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
    {
        switch (c)
        {
            case 'h':
                printf("grit -m motif -i homoseq -b bgseq [-n sampling_times] [-z sampling_size] [-s rs] [-t pvalue] [-p pscore] -o output\n");
                exit(0);
                break;
            case 'm':
                mat_file = optarg;
                break;
            case 'i':
                seq_file = optarg;
                break;
            case 'b':
                bgseq_file = optarg;
                break;
                //cout<<optarg<<endl;
                //target_n = atoi(optarg);
                //break;
            case 'n':
                sampling_times = atoi(optarg);
                break;
            case 'z':
                sampling_size = atoi(optarg);
                break;
            case 's':
                lratio = atoi(optarg);
                break;
            case 't':
                threshold = atof(optarg);
                break;
            case 'p':
                //cout<<optarg<<endl;
                p_score = atof(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            default :
                cout << "invalid option: " << optarg << endl;
                exit(1);
        }
    } 

    raw_score raw_score_instance;
    raw_score_instance.threshold = threshold;
    raw_score_instance.output_file = output_file;
    raw_score_instance.hit_thresh = p_score;
    raw_score_instance.sampling_times = sampling_times;
    raw_score_instance.sampling_times = sampling_times;
    raw_score_instance.lratio = lratio;
    //raw_score_instance.permutation_n = permutation_n;
    //raw_score_instance.dir = output_dir;
    //raw_score_instance.session = session_name;

    raw_score_instance.mat_file = mat_file;
    raw_score_instance.get_motifs();
    //cout << "Read motifs finished" << endl;

    raw_score_instance.bgseq_file = bgseq_file;
    //raw_score_instance.get_seqs();
    raw_score_instance.get_bgseqs();
    raw_score_instance.get_rdm_bsgroups();
    raw_score_instance.get_mutibgbase_probs();
    //raw_score_instance.get_bgsbase_probs();
    //cout << "Read background sequences finished" << endl;

    raw_score_instance.seq_file = seq_file;
    //raw_score_instance.get_seqs();
    raw_score_instance.get_mutiseqs();
    //cout << "Read sequences finished" << endl;

    raw_score_instance.get_mutibase_probs();
    //raw_score_instance.analysis_motifs_wks();
    //cout<<"test 03 finished."<<endl;
    raw_score_instance.analysis_motifs_st();

    //cout << "test07" << endl;
    //cout << "Analysis sequences finished" << endl;

}

