/***************************************
*
* GEREA
* Gene Expression Reulator Enrhchment Analysis
* This program was designed to
* search the enriched gene expresion
* regulators
* Mod Data: 10/20/2019
* Author: Tinghua Huang
*
***************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>  // accumulate
#include <string>
#include <vector>
#include <cmath>  // log
#include <libgen.h>  // basename
#include <assert.h>
#include<sstream>

#include "matrix.h"  // my 2D matrix class
//#include "wks.h"

using namespace std;
using mcf::matrix;
using namespace mcf;
using mcf::mcf_die;

class raw_score {
typedef unsigned uint;

template<typename _Tp>
struct greater : public binary_function<_Tp, _Tp, bool>
{
    bool operator()(const _Tp& __x, const _Tp& __y) const
    { return __x > __y; }
};

template<typename _Tp>
struct less : public binary_function<_Tp, _Tp, bool>
{
    bool operator()(const _Tp& __x, const _Tp& __y) const
    { return __x < __y; }
};

public:
	//using mcf::markov;
	//using args::alphsize;
	string mat_file = "motif.txt";
	string seq_file = "seq.txt";
	string bgseq_file = "bgseq.txt";
	string stat_file = "stat.txt";
	string score_file = "score.txt";
	string output_file = "output.txt";
	const uint alphsize = 4;
	double pseudocount = 0.375;
	bool mask_lower = false;
	double hit_thresh = 6.0;
	string session = "my_work";
	string dir = "./";
	unsigned long permutation_n = 1000;
	double threshold = 0.05;
	double ifnv = -1.79e308;
	double lowest_logsocre = -1.0e10;
	unsigned int sampling_times = 10;
	unsigned int sampling_size = 200;
	unsigned int lratio = 1;

	//typedef unsigned uint;
	typedef vector<matrix<double> > matvec;

	struct seq_set_info {  // misc info about a sequence set
		size_t num;  // number of seqs
		size_t len;  // total bp
		double gc;   // GC composition

		//seq_set_info(const vector<vector<uint> > & seqs);
	};

	struct result {
		//uint motif;  // index of motif
		double raw_score;
		vector<double> pvalues;
		vector<double> seq_scores;  // score for each sequence

		bool operator > (const result & r) const
		{ return raw_score > r.raw_score; }
		bool operator < (const result & r) const
		{ return raw_score < r.raw_score; }
	};

	struct hit {
		bool found;
		uint motif;  // index of motif
		uint strand;  // 0 or 1: palindromes??
		uint location;
		double score;
		string site;
		bool operator< (const hit & h) const
		{ return location < h.location; }
	};

	struct stat {
		//uint motif;  // index of motif
		string motif_name;
		string homoseq_name;
		vector<pair<string, double> > gene_list;
		double pvalue, fdr;
		double ccv, sd;
		double aglib_pvalue;
		double nc_pvalue;
		double diff, uc_diff;
		double score;
		hit bst_hit;
		bool operator > (const stat & r) const
		{ return pvalue > r.pvalue; }
		bool operator < (const stat & r) const
		{ return pvalue < r.pvalue; }

	};

	struct homo_seq {
		vector<vector<uint> > seqs;
		vector<string> seq_names;
		string homo_name;

	};

	vector<matrix<float> > ss_motifs;  // single stranded motifs
	vector<matvec> ds_motifs;  // double strand motifs
	vector<string> motif_names;
	vector<vector<uint> > bgseqs;
	vector<vector<unsigned int> > rdm_bsgroups;
	vector<homo_seq> mutiseqs;

	//vector<vector<vector<uint> > > mutiseqs;
	//vector<string> homo_names;
	//vector<vector<string> > mutiseq_names;
	vector<string> bgseq_names;
	//vector<vector<double> > bg_base_probs;  // residue abundances per sequence
	vector<vector<vector<double> > > mutibg_base_probs;
	vector<vector<vector<double> > > muti_base_probs;  // residue abundances per sequence

	//vector<double> dp;  // dynamic programming matrix used in combine_scores
	vector<double> bgscores;
	vector<result> results;  // motif raw scores & pvalues
	vector<vector<hit> > hits;  // motif locations in each sequence
	//seq_set_info seq_info;
	vector<seq_set_info> bgseq_infos;

	unsigned (*translator)(char);

	  // (I think uint(-1) == ~0u == UINT_MAX)
	//double scan_seq(const vector<uint> & seq, const matvec & motif, const vector<double> & base_probs, uint seqnum = ~0u, uint motnum = ~0u);
	void scan_seq(const vector<uint> & seq, const matvec & motif, const vector<double> & base_probs, unsigned int &lratio, double &rs, bool &success);
	void scan_hits(const vector<uint> & seq, const matvec & motif, const vector<double> & base_probs, hit &bst_hit, bool &found_hit);
	  // void init_dp(uint num_seqs);
	  // double combine_scores(const vector<double> & scores);
	  // void get_hits();
	//void get_raw_scores();

	//void get_raw_scores(matvec &ds_motif, vector<vector<uint> > &seqs, vector<vector<double> > &base_probs, result &myresult);
	//void get_raw_scores(matvec &ds_motif, vector<vector<uint> > &seqs, vector<vector<double> > &base_probs, vector<double> &raw_scores);
	//void get_bg_raw_scores();

	  // void rand_test(const vector<vector<uint> > & myseqs, const vector<vector<double> > & b_probs, vector<matvec> & motifs, vector<uint> & losses);
	  // void copy_masks(const vector<uint> & source, vector<uint> & dest);
	void get_base_probs(const vector<uint> & seq, vector<double> & probs);
	void get_mutibase_probs();
	void get_mutibgbase_probs();
	  // double get_gc_composition(const vector<vector<uint> > & seqs);
	  // void bg_fragment(const vector<vector<uint> > & bg_seqs, vector<uint> & frag, uint len, const vector<uint> & frag_num, uint frag_tot);
	  // void shuffle_bgseq(const vector<vector<uint> > & bg_seqs);
	  // void shuffle_mono();
	  // void shuffle_di();
	  // void shuffle_mat();
	void get_ss_motifs();
	void fifth_column(const matrix<float> & m1, matrix<double> & m2);
	void get_ds_motifs(const vector<matrix<float> > & ss_mots, vector<matvec> & ds_mots);
	void get_motifs();
	void get_bgseqs();
	void get_rdm_bsgroups();
	bool get_fasta(std::istream & strm, std::vector<unsigned> & seq, std::string & title, unsigned (*translator)(char));
	void get_mutiseqs();
	bool get_homoseq(std::istream & strm, homo_seq &homoseq);
	
	void get_fastaseqs(std::istream & strm, homo_seq &homoseq);

	  // bool is_significant(const result & r);
	  // void print_hits(uint wn, uint ws);
	  // void print_per_sequence(uint wn);
	  // void print_results();
	istream & get_simple_pssm(std::istream & strm, matrix<float> & mat, string & title, unsigned alphsize);  // not allowed to repeat the default argument?
	void normalize_pssm(matrix<float> & pssm, const std::vector<float> & pseudos);
	void count_residues(const std::vector<unsigned> & seq, std::vector<size_t> & counts, unsigned alphsize);
	//void set_seq_info(const vector<vector<uint> > & seqs);
	void set_seq_info(const vector<vector<uint> > & seqs, vector<vector<unsigned int> > &groups);

	//void analysis_motifs();
	//void analysis_motifs_wks();
	void analysis_motifs_st();
	void print_score(ofstream &writer, string &sname, vector<pair<string, double> > &gene_score);
	void print_score(ofstream &writer, string &sname, vector<double> &gene_score);
	//void wks_analysis();

};

string doubleToString(const double &dbNum, unsigned long dot_n);
string long_to_string(unsigned long val);
bool stat_less_p(raw_score::stat a, raw_score::stat b);
bool stat_less_fdr(raw_score::stat a, raw_score::stat b);
bool hit_larger_score(raw_score::hit a, raw_score::hit b);
double mixed_studentt_avg(double one, vector<double> &obs, vector<vector<double> > &bkg, double &ccv, double &sd);
double average(vector<vector<double> > &bkg);

