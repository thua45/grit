/***************************************
*
* RAW_SCORE
* This program was designed to
* caculate the raw score of the transcription
* factor binding site in promoter sequence
* Mod Data: 12/20/2019
* Author: Tinghua Huang
*
***************************************/

#include "cdflib.hpp"
#include "grit.h"
//#include "statistics.h"
//using namespace alglib_impl;

inline unsigned nolower_translator(char c);
inline unsigned DNA_to_number(char c);
inline char number_to_DNA(unsigned b);
double average(vector<double> &data_vec);

double average(vector<double> &data_vec) {
  double total = 0;
  if (data_vec.size() == 0) {
    return 0.0;
  }
  for (unsigned int di = 0; di < data_vec.size(); ++di) {
    total += data_vec[di];
  }
  return total / data_vec.size();
}

// read motifs from a file
void raw_score::get_motifs()
{
  get_ss_motifs();
  get_ds_motifs(ss_motifs, ds_motifs);
}

void raw_score::get_ss_motifs()
{
  ifstream file(mat_file.c_str());
  if (!file) mcf_die("Sorry, couldn't open " + mat_file);

  const vector<float> pseudos(alphsize, pseudocount);
  matrix<float> matf;
  string title;

  while (get_simple_pssm(file, matf, title, alphsize)) {
    if (matf.rows() == 0) mcf_die("Empty matrix not allowed: " + title);
    motif_names.push_back(title);
    normalize_pssm(matf, pseudos);
    //    transform(matf[0], matf[1], matf[0], bind2nd(multiplies<float>(), 0.5f));
    ss_motifs.push_back(matf);
  }

  if (!file.eof())  // catches some but not all errors
    mcf_die("Sorry, couldn't understand the matrix file " + mat_file);
}

istream & raw_score::get_simple_pssm(std::istream & strm, matrix<float> & mat, string & title, unsigned alphsize) {
std::string t;
  matrix<float> m(0, alphsize);
  char c = 0;
  bool titflag = false;  // have we read a title line yet?

  while (strm >> c) {
    if (c == '>') {
      if (titflag || m.rows() != 0) {
	      strm.unget();
	      break;
      } else {
	      std::getline(strm, t);
	      titflag = true;
      }
    } else if (c == '#') {  // skip comments
      std::string junk;
      std::getline(strm, junk);
    } else {
      strm.unget();
      std::vector<double> v;
      for (unsigned i = 0; i < alphsize; ++i) {
	      double d;
	      strm >> d;
	      v.push_back(d);
      }
      if (!strm)
	      return strm;  // failed to read alphsize doubles
        m.push_row(v.begin());
    }
  }

  // if reached EOF but read something, clear the stream state:
  if (strm.eof() && (titflag || m.rows() != 0))
    strm.clear();

  if (strm) {
    title = t;
    mat = m;
  }
  return strm;

}

void raw_score::normalize_pssm(matrix<float> & pssm, const std::vector<float> & pseudos) {
  const unsigned r = pssm.rows();
  const unsigned c = pssm.cols();
  assert(pseudos.size() == c);

  for (std::vector<float>::iterator k = pssm[0]; k < pssm[r]; k += c) {
    for (unsigned i = 0; i < c; ++i)
      k[i] += pseudos[i];
    normalize(k, k+c);
  }
}

// finish preprocessing matrices
// 3. add an extra column of all zeros (for 'n' bases)
// 4. make reverse complements
void raw_score::get_ds_motifs(const vector<matrix<float> > & ss_mots, vector<matvec> & ds_mots)
{
  for (uint m = 0; m < ss_mots.size(); ++m) {
    matrix<float> matf(ss_mots[m]);
    matrix<double> matd(matf.rows(), alphsize+1);
    ds_mots.push_back(matvec());
    fifth_column(matf, matd);
    ds_mots.back().push_back(matd);
    //continue;
    /*
    if (args::strand == 2 &&
	    !matf.is_rotate180()) {  // if the motif isn't palindromic
      matf.rotate180();  // reverse complement
      fifth_column(matf, matd);
      ds_mots.back().push_back(matd);
    }
    */
  }
}

// add a column to a matrix, with all cells = zero
void raw_score::fifth_column(const matrix<float> & m1, matrix<double> & m2)
{
  for (uint r = 0; r < m1.rows(); ++r) {
    copy(m1[r], m1[r+1], m2[r]);
    m2[r][alphsize] = 0;
  }
}

void raw_score::get_mutiseqs()
{
  ifstream file(seq_file.c_str());
  if (!file) mcf_die("Sorry, couldn't open file " + seq_file);
  translator = (mask_lower ? nolower_translator : DNA_to_number);

  while (1 == 1) {
    homo_seq one_homoseq;
    mutiseqs.push_back(one_homoseq);
    //string n;
    if (!get_homoseq(file, mutiseqs.back()))
      break;
    //seq_names.push_back(n);
  }

  mutiseqs.pop_back();
  //set_seq_info(seqs);

}

// Industrial strength fasta-format sequence reader
// appends the sequence to anything already in seq
bool raw_score::get_homoseq(std::istream & strm, homo_seq &homoseq)
{
  std::string ht;
  //unsigned old_size = homoseq.seqs.size();
  char c;
  char junk_c;
  bool titflag = false;  // have we read a title line yet?

  while (strm >> c) {
    if (c == '/') {
      if (titflag) {
	      strm.unget();
	      break;
      } else {
        strm >> junk_c;
	      getline(strm, ht);
        homoseq.homo_name = ht;
	      titflag = true;
        get_fastaseqs(strm, homoseq);
        break;
        //titflag = false;
      }
    } else if (c == '#') {  // skip comments
      string junk;
      getline(strm, junk);
    }
  }

  // if reached EOF but read something, clear the stream state:
  if (strm.eof()) {
    strm.clear();
    return false;
  }
  else
    return true;

}

void raw_score::get_fastaseqs(std::istream & strm, homo_seq &homoseq)
{
  while (1 == 1) {
    homoseq.seqs.push_back(vector<uint>(0));  // push back an empty vector
    string n;
    if (!get_fasta(strm, homoseq.seqs.back(), n, translator))
    {
      homoseq.seq_names.push_back(n);
      break;
    }
    homoseq.seq_names.push_back(n);
  }

  //homoseq.seqs.pop_back();
  //set_seq_info(seqs);

}

void raw_score::get_bgseqs()
{
  ifstream file(bgseq_file.c_str());
  if (!file) mcf_die("Sorry, couldn't open file " + bgseq_file);
  unsigned (*translator)(char) = (mask_lower ? nolower_translator : DNA_to_number);

  char c;
  while (file >> c) {
    if (c == '>') {
      file.unget();
      bgseqs.push_back(vector<uint>(0));  // push back an empty vector
      string n;
      if (!get_fasta(file, bgseqs.back(), n, translator))
      {
        bgseq_names.push_back(n);
        break;
      }
      bgseq_names.push_back(n);
    }
  }
  //set_seq_info(bgseqs);
  //set_seq_info(bgseqs, rdm_bsgroups);
  /*
  while (1 == 1) {
    bgseqs.push_back(vector<uint>(0));  // push back an empty vector
    string n;
    if (!get_fasta(file, bgseqs.back(), n, translator))
      break;
    bgseq_names.push_back(n);
  }

  bgseqs.pop_back();
  set_seq_info(bgseqs);
  */

}

void raw_score::get_rdm_bsgroups()
{
  srand((unsigned int)(time(0))); //important!

  list<unsigned int> index_pool;
  list<unsigned int>::iterator itor;
  for (unsigned int i = 0; i < bgseqs.size(); i++) {
    index_pool.push_back(i);
  }
  //srand((unsigned int)time(NULL));
  vector<unsigned int> tmp_samples;
  for (unsigned int stime = 0; stime < sampling_times; stime++) {
    tmp_samples.clear();
    for (unsigned int ssize = 0; ssize < sampling_size; ssize++) {
      if (index_pool.size() == 0) {
        cout<<"not enough number of bgseqs for randomization."<<endl;
        exit(1);
      }
      int rdm_i = rand() % index_pool.size();
      unsigned int tmp_i = 0;
      itor = index_pool.begin();
      while(itor != index_pool.end()) {
        if (rdm_i == tmp_i) {
          tmp_samples.push_back(*itor);
          index_pool.erase(itor);
          break;
        }
        itor++;
        tmp_i++;
      }
    }
    rdm_bsgroups.push_back(tmp_samples);
  }

  set_seq_info(bgseqs, rdm_bsgroups);

}

// DNA-to-number translator that only recognizes uppercase
inline unsigned nolower_translator(char c)
{
  switch(c) {
  case 'A':
    return 0u;
  case 'C':
    return 1u;
  case 'G':
    return 2u;
  case 'T':
    return 3u;
  default:
    return 4u;
  }
}

inline char number_to_DNA(unsigned b)
{
  static const char lookup[] = "acgtn";
  assert(b < 5);
  return lookup[b];
}

// should use static array instead of switch - ask Jianhua how
inline unsigned DNA_to_number(char c)
{
  switch(c) {
  case 'a':
  case 'A':
    return 0u;
  case 'c':
  case 'C':
    return 1u;
  case 'g':
  case 'G':
    return 2u;
  case 't':
  case 'T':
    return 3u;
  default:
    return 4u;
  }
}

// Industrial strength fasta-format sequence reader
// appends the sequence to anything already in seq
bool raw_score::get_fasta(
  std::istream & strm,
  std::vector<unsigned> & seq,
  std::string & title,
  unsigned (*translator)(char))  // not allowed to repeat the default argument?
{
  std::string t;
  unsigned old_size = seq.size();
  char c;
  bool titflag = false;  // have we read a title line yet?
  bool homoflag = false;

  while (strm >> c) {
    if (c == '>') {
      if (titflag || seq.size() != old_size) {
	      strm.unget();
	      break;
      } else {
	      std::getline(strm, t);
	      titflag = true;
        title = t;
      }
    } else if (c == '/') {
	    strm.unget();
      homoflag = true;
	    break;
    } else if (c == '#') {  // skip comments
      std::string junk;
      std::getline(strm, junk);
    } else if (isalpha(c)) {
      seq.push_back(translator(c));
    }
  }

  // if reached EOF but read something, clear the stream state:
  if (strm.eof()) {
    strm.clear();
    return false;
  } else if (homoflag) {
    return false;
  } else 
    return true;
    
}

void raw_score::set_seq_info(const vector<vector<uint> > & seqs, vector<vector<unsigned int> > &groups) {
  for (unsigned int gi = 0; gi < groups.size(); gi++) {
    seq_set_info seq_info;
    seq_info.num = groups[gi].size();
    seq_info.len = 0;
    for (uint si = 0; si < groups[gi].size(); ++si)
      seq_info.len += seqs[groups[gi][si]].size();

    vector<size_t> counts;  // no Laplace's rule
    for (size_t si = 0; si < groups[gi].size(); ++si)
      count_residues(seqs[groups[gi][si]], counts, alphsize);
    const size_t zero = 0;
    seq_info.gc = double(counts[1] + counts[2]) / accumulate(counts.begin(), counts.end(), zero);
    bgseq_infos.push_back(seq_info);
  }
}

/*
void raw_score::set_seq_info(const vector<vector<uint> > & seqs) {
  seq_info.num = seqs.size();
  seq_info.len = 0;
  for (uint s = 0; s < seqs.size(); ++s)
    seq_info.len += seqs[s].size();

  vector<size_t> counts;  // no Laplace's rule
  for (size_t s = 0; s < seqs.size(); ++s)
    count_residues(seqs[s], counts, alphsize);
  const size_t zero = 0;
  seq_info.gc = double(counts[1] + counts[2]) / accumulate(counts.begin(), counts.end(), zero);

}
*/

void raw_score::count_residues(const std::vector<unsigned> & seq, std::vector<size_t> & counts, unsigned alphsize)
{
  if (counts.size() < alphsize)
    counts.resize(alphsize);

  for (std::vector<unsigned>::const_iterator n = seq.begin(); n < seq.end(); ++n) {
    unsigned i = *n;
    if (i < alphsize)
      ++counts[i];
  }
}

/*
void raw_score::get_raw_scores(matvec &ds_motif, vector<vector<uint> > &seqs, vector<vector<double> > &base_probs, vector<double> &raw_scores) {
  for (uint s = 0; s < seqs.size(); ++s)
    raw_scores.push_back(scan_seq(seqs[s], ds_motif, base_probs[s], lratio));
    //results[m].raw_score = combine_scores(results[m].seq_scores);
  


}
*/

/*
void raw_score::get_raw_scores(matvec &ds_motif, vector<vector<uint> > &seqs, vector<vector<double> > &base_probs, result &myresult)
{


  for (uint s = 0; s < seqs.size(); ++s)
    myresult.seq_scores.push_back(scan_seq(seqs[s], ds_motif, base_probs[s], lratio));
    //results[m].raw_score = combine_scores(results[m].seq_scores);
  
  for (uint s = 0; s < seqs.size(); ++s) {
    cout << myresult.seq_scores[s] << endl;
  }
  cout << endl;

}
*/

/*
void raw_score::get_bg_raw_scores()
{
  cerr << "Getting raw scores..." << endl;

  bg_base_probs.resize(bgseqs.size());
  for (unsigned s = 0; s < bgseqs.size(); ++s)
    get_base_probs(bgseqs[s], bg_base_probs[s]);

  bgscores.resize(bgseqs.size());
  
  for (uint m = 0; m < ds_motifs.size(); ++m) {
    bgscores.clear();
    for (uint s = 0; s < bgseqs.size(); ++s)
      bgscores.push_back(scan_seq(bgseqs[s], ds_motifs[m], bg_base_probs[s]));
    //results[m].raw_score = combine_scores(results[m].seq_scores);
    cout << motif_names[m] << endl;
    for (uint s = 0; s < bgseqs.size(); ++s) {
      cout << bgscores[s] << endl;
    }
    cout << endl;
  }

}
*/

// get background abundances of residues in the sequence
void raw_score::get_mutibase_probs()
{
  for (unsigned int hi = 0; hi < mutiseqs.size(); ++hi) {
    vector<vector<double> > base_prob;
    muti_base_probs.push_back(base_prob);
    muti_base_probs.back().resize(mutiseqs[hi].seqs.size());
    for (unsigned s = 0; s < mutiseqs[hi].seqs.size(); ++s)
      get_base_probs(mutiseqs[hi].seqs[s], muti_base_probs.back()[s]);
  }
  
}

// get background abundances of residues in the sequence
void raw_score::get_mutibgbase_probs()
{
  for (unsigned int gi = 0; gi < rdm_bsgroups.size(); gi++) {
    vector<vector<double> > base_prob;
    mutibg_base_probs.push_back(base_prob);
    mutibg_base_probs.back().resize(rdm_bsgroups[gi].size());
    for (unsigned si = 0; si < rdm_bsgroups[gi].size(); ++si)
      get_base_probs(bgseqs[rdm_bsgroups[gi][si]], mutibg_base_probs.back()[si]);
  }

}

// get background abundances of residues in the sequence
void raw_score::get_base_probs(const vector<uint> &seq, vector<double> &probs)
{
  probs.clear();  // just in case
  vector<size_t> counts;  // no Laplace's rule
  count_residues(seq, counts, alphsize);
  const size_t zero = 0;
  const double tot = accumulate(counts.begin(), counts.end(), zero);
  for (uint i = 0; i < alphsize; ++i)
    probs.push_back(counts[i] / tot);
}

// get the score for one motif in one sequence:
// prob(sequence | motif) / prob(sequence | null model)
// default arguments shouldn't be repeated, apparently
void raw_score::scan_seq(const vector<uint> & seq, const matvec & motif, const vector<double> & b_probs, unsigned int &lratio, double &rs, bool &success)
{
  //cout<<"test-04.2"<<endl;
  for (uint c = 0; c < alphsize; ++c) {
	  if (b_probs[c] == 0.0) {
      //return 0.0;
      success = false;
      return;
    }
  }
  //cout<<"test-04.3"<<endl;

  //double rs = 0;
  unsigned int ms = motif.size();
  //cout << "motif_size: " << ms << endl;

  // iterate over motif orientations:
  for (uint m = 0; m < ms; ++m) {
    matrix<double> pssm(motif[m]);
    if (pssm.rows() > seq.size())
      continue;

    unsigned int pssm_len = pssm.rows();
    //cout << "pssm_len: " << pssm_len << endl;

    // incorporate base probs into PSSM:
    for (uint r = 0; r < pssm.rows(); ++r)
      for (uint c = 0; c < alphsize; ++c)
	      pssm[r][c] /= b_probs[c];

    // finally, scan the PSSM against the sequence:
    if (lratio == 1) {
      double score = 0;
      double score_avg = 0;
      const uint posns = seq.size() - pssm.rows() + 1;
      for (vector<uint>::const_iterator n = seq.begin(); n != seq.begin() + posns; ++n) {
        double s = 1;
        for (uint k = 0; k < pssm.rows(); ++k)
          s *= pssm[k][*(n+k)];
        if (s > 0.0) {
          //score += log(s) / pssm_len;
          score += s / pssm_len;
          //score += s;
        }
      }
      score_avg = score / posns;
      if (score_avg > 0) {
        rs = log(score_avg);
        success = true;
      }
      else {
        success = false;
      }
      
    } else if (lratio == 2) {
      //double score = 0;
      double max_score = 0;
      const uint posns = seq.size() - pssm.rows() + 1;
      for (vector<uint>::const_iterator n = seq.begin(); n != seq.begin() + posns; ++n) {
        double s = 1;
        for (uint k = 0; k < pssm.rows(); ++k)
          s *= pssm[k][*(n+k)];
        if (s > 0) {
          double llr_score = s / pssm_len;
          if (llr_score > max_score) {
            max_score = llr_score;
          }
        }
      }
      if (max_score > 0) {
        rs = log(max_score);
        success = true;
      }
      else {
        success = false;
      }
      
    }
    else {
      cout << "lratio = " << lratio << " unsupported." << endl;
      exit(1);
    }
  }

  //cout << "tot_score: " << tot_score << endl;

  /*
  if (tot_score > 0.0) {
    return log(tot_score);
  }
  else {
    return lowest_logsocre;
  }
  */
  //return rs;

}

void raw_score::scan_hits(const vector<uint> & seq, const matvec & motif, const vector<double> & b_probs, hit &bst_hit, bool &found_hit)
{
  double tot_score = 0;

  // iterate over motif orientations:
  for (uint m = 0; m < motif.size(); ++m) {
    matrix<double> pssm(motif[m]);
    if (pssm.rows() > seq.size())
      continue;
    
    unsigned int pssm_len = pssm.rows();
    // incorporate base probs into PSSM:
    for (uint r = 0; r < pssm_len; ++r)
      for (uint c = 0; c < alphsize; ++c)
	      pssm[r][c] /= b_probs[c];

    // finally, scan the PSSM against the sequence:
    //double score = 0;
    vector<hit> seq_hits;
    const uint posns = seq.size() - pssm_len + 1;
    double pscore;
    for (vector<uint>::const_iterator n = seq.begin(); n != seq.begin() + posns; ++n) {
      double s = 1;
      for (uint k = 0; k < pssm_len; ++k)
	      s *= pssm[k][*(n+k)];
      //score += s;
      //cout << "s: " << s << endl;
      //cout << "s: " << s << ", log s: " << log(s) << endl;
      if (s > 0.0) {
        pscore = log(s) / pssm_len;
        //pscore = log(s);
      }
      if (pscore >= hit_thresh) {
        hit a_hit;
        //a_hit.motif = ???;
        a_hit.strand = m;
        a_hit.location = n - seq.begin();
        a_hit.score = pscore;
        const uint motif_width = motif[m].rows();
        string site;
        for (unsigned int k = 0; k < motif_width; ++k)
	        site += number_to_DNA(seq[a_hit.location + k]);
        a_hit.site = site;
	      seq_hits.push_back(a_hit);
      }
    }

    sort(seq_hits.begin(), seq_hits.end(), hit_larger_score);
    if (seq_hits.size() > 0) {
      found_hit = true;
      bst_hit = seq_hits[0];
    }
    else {
      found_hit = false;
    }
    //tot_score += score / posns / motif.size();
  }

  //return tot_score;
}

/*
void raw_score::scan_hits(const vector<uint> & seq, const matvec & motif, const vector<double> & b_probs, hit &bst_hit, bool &found_hit)
{
  double tot_score = 0;

  // iterate over motif orientations:
  for (uint m = 0; m < motif.size(); ++m) {
    matrix<double> pssm(motif[m]);
    if (pssm.rows() > seq.size())
      continue;

    // incorporate base probs into PSSM:
    for (uint r = 0; r < pssm.rows(); ++r)
      for (uint c = 0; c < alphsize; ++c)
	      pssm[r][c] /= b_probs[c];

    // finally, scan the PSSM against the sequence:
    //double score = 0;
    vector<hit> seq_hits;
    const uint posns = seq.size() - pssm.rows() + 1;
    for (vector<uint>::const_iterator n = seq.begin(); n != seq.begin() + posns; ++n) {
      double s = 1;
      for (uint k = 0; k < pssm.rows(); ++k)
	      s *= pssm[k][*(n+k)];
      //score += s;
      //cout << "s: " << s << endl;
      //cout << "s: " << s << ", log s: " << log(s) << endl;
      if (s > 0.0  && log(s) >= hit_thresh) {
        hit a_hit;
        //a_hit.motif = ???;
        a_hit.strand = m;
        a_hit.location = n - seq.begin();
        a_hit.score = s;
        const uint motif_width = motif[m].rows();
        string site;
        for (unsigned int k = 0; k < motif_width; ++k)
	        site += number_to_DNA(seq[a_hit.location + k]);
        a_hit.site = site;
	      seq_hits.push_back(a_hit);
      }
    }

    sort(seq_hits.begin(), seq_hits.end(), hit_larger_score);
    if (seq_hits.size() > 0) {
      found_hit = true;
      bst_hit = seq_hits[0];
    }
    else {
      found_hit = false;
    }
    //tot_score += score / posns / motif.size();
  }

  //return tot_score;
}
*/

bool stat_less_p(raw_score::stat a, raw_score::stat b)
{
  return a.pvalue < b.pvalue; //升序排列，如果改为return a>b，则为降序
}

bool stat_less_fdr(raw_score::stat a, raw_score::stat b)
{
  return a.fdr < b.fdr; //升序排列，如果改为return a>b，则为降序
}

bool hit_larger_score(raw_score::hit a, raw_score::hit b)
{
  return a.score > b.score; //升序排列，如果改为return a>b，则为降序
}


void raw_score::analysis_motifs_st() {
  srand((unsigned int)(time(0)));  //important!

  //ofstream score_writer;
  //string score_file = dir + "/" + session + ".score.txt";
  //string score_file = "test.score.txt";
  //score_writer.open(score_file.c_str(), ios::out);
  //if (!score_writer) mcf_die("Sorry, couldn't open " + score_file);
  ofstream stat_writer;
  //string stat_file = dir + "/" + session + ".stat.txt";
  string stat_file = output_file;
  stat_writer.open(stat_file.c_str(), ios::out);
  if (!stat_writer) mcf_die("Sorry, couldn't open " + stat_file);
  stat_writer<< "# motif_name\thomoseq_name\tbinding_site\tlocation\tscore\tccv\tsd\tpvalue\tfdr" << endl;

  vector<pair<string, double> > gene_score;
  vector<pair<string, double> > gene_list;
  //double bg_min, bg_max, bg_total_score;

  //bg_base_probs.resize(bgseqs.size());
  //for (unsigned s = 0; s < bgseqs.size(); ++s)
  //  get_base_probs(bgseqs[s], bg_base_probs[s]);

  double rs;
  bool rs_success;

  for (unsigned int mi = 0; mi < ds_motifs.size(); ++mi) {
    //cout << motif_names[mi] << endl;
    //score_writer << "//" << motif_names[mi] <<"\n";

    vector<vector<double> > bgseq_scores;
    vector<double> obs_scores;
    double score_one;
    //double pvalue = two_sample_ttest(exp_a, exp_b);
    bgseq_scores.clear();
    //double tmp_score;
    vector<double> tmp_score_vec;
    //cout<<"test 04.0 finished."<<endl;
    for (unsigned int gi = 0; gi < rdm_bsgroups.size(); gi++) {
      tmp_score_vec.clear();
      for (unsigned si = 0; si < rdm_bsgroups[gi].size(); ++si) {
      //for (uint s = 0; s < bgseqs.size(); ++s) {
        //add group info
        //cout<<"test 04.1 finished."<<endl;
        scan_seq(bgseqs[rdm_bsgroups[gi][si]], ds_motifs[mi], mutibg_base_probs[gi][si], lratio, rs, rs_success);
        //cout << tmp_score << endl;
        if (rs_success) {
          tmp_score_vec.push_back(rs);
        }
        
        //if (tmp_score > 0.0) {
        //  tmp_score_vec.push_back(tmp_score);
        //}
        //exp_a.push_back(scan_seq(bgseqs[s], ds_motifs[mi], bg_base_probs[s]));
      }
      bgseq_scores.push_back(tmp_score_vec);
      //string sname = "background";
      //print_score(score_writer, sname, exp_a);
    }
    //cout<<"test 04 finished."<<endl;

    //double average_bgseq_score = average(bgseq_scores);

    vector<stat> homo_stats;
    homo_stats.clear();
    unsigned int mutiseq_size = mutiseqs.size();
    //cout << mutiseq_size << endl;
    for (unsigned int hi = 0; hi < mutiseqs.size(); ++hi) {
      //cout << mutiseqs[hi].seqs.size() << endl;
      if (mutiseqs[hi].seqs.size() < 2)
        continue;
      scan_seq(mutiseqs[hi].seqs[0], ds_motifs[mi], muti_base_probs[hi][0], lratio, rs, rs_success);
      //cout << "score2" << endl;
      //cout << tmp_score << endl;
      if (rs_success) {
        score_one = rs;
      }
      else {
        continue;
      }
      
      //if (tmp_score > 0.0)
      //  score_one = tmp_score;
      //else
      //  continue;
        //score_one = -1.79e308;

      obs_scores.clear();
      for (uint s = 1; s < mutiseqs[hi].seqs.size(); ++s) {
        scan_seq(mutiseqs[hi].seqs[s], ds_motifs[mi], muti_base_probs[hi][s], lratio, rs, rs_success);
        if (rs_success) {
          obs_scores.push_back(rs);
        }
        //if (tmp_score > 0.0) {
        //  obs_scores.push_back(tmp_score);
        //}
        //exp_b.push_back(scan_seq(mutiseqs[hi].seqs[s], ds_motifs[mi], muti_base_probs[hi][s]));
      }
      //cout << mutiseqs[hi].seq_names[0] << endl;

      //string sname = "observes";
      //print_score(score_writer, sname, exp_b);

      stat a_stat;
      a_stat.motif_name = motif_names[mi];
      a_stat.homoseq_name = mutiseqs[hi].homo_name;
      //a_stat.aglib_pvalue = two_sample_ttest(exp_a, exp_b);
      if (obs_scores.size() == 0) {
        a_stat.pvalue = 1.0;
        //cout << "test01" << endl;
      }
      else {
        a_stat.pvalue = mixed_studentt_avg(score_one, obs_scores, bgseq_scores, a_stat.ccv, a_stat.sd);
        //a_stat.pvalue = studentt_one(score_one, obs_scores, bgseq_scores[0]);
        //a_stat.pvalue = studentt_one_sample(score_one, exp_a);
        //a_stat.pvalue = studentt_one_sample(score_one, exp_b);
        //a_stat.pvalue = studentt_two_sample(exp_b, exp_a);
        //cout << "test02" << endl;
        //a_stat.pvalue = studentt_one(exp_b, exp_a);
        //cout << "test03" << endl;
      }
      //cout << "pvalue: " << a_stat.pvalue << endl;
      //a_stat.nc_pvalue = 1.0 - studentt_one(exp_b[0], exp_b);
      //a_stat.diff = average(exp_b) - average_a;
      //!!!a_stat.diff = score_one - average_bgseq_score;
      //cout << "score_one - average_a: " << score_one << ", " << average_a << endl;
      //a_stat.score = exp_b[0];
      //a_stat.uc_diff = exp_b[0] - average(exp_b);

      bool found_hit;
      scan_hits(mutiseqs[hi].seqs[0], ds_motifs[mi], muti_base_probs[hi][0], a_stat.bst_hit, found_hit);
      if (found_hit) {
        a_stat.bst_hit.found = true;
      } else {
        a_stat.bst_hit.found = false;
      }
      if (a_stat.pvalue <= threshold) {
        //a_stat.gene_list = gene_list;
        //print_score(score_writer, mutiseqs[hi].homo_name, gene_list);
      }
      homo_stats.push_back(a_stat);
    }

    //cout << "test04" << endl;

    sort(homo_stats.begin(), homo_stats.end(), stat_less_p);
    double fdr;
    unsigned int homo_stats_size = homo_stats.size();
    for (unsigned int hi = 0; hi < homo_stats.size(); ++hi) {
      fdr = homo_stats[hi].pvalue * homo_stats_size / (hi + 1.0);
      if (fdr <= 1.0) {
        homo_stats[hi].fdr = fdr;
      }
      else {
        homo_stats[hi].fdr = 1.0;
      }
      
    }

    //cout << "test05" << endl;

    sort(homo_stats.begin(), homo_stats.end(), stat_less_fdr);
    for (unsigned int hi = 0; hi < homo_stats.size(); ++hi) {
      if (homo_stats[0].fdr <= threshold) {
        string sname = "BACKGROUND";
        //print_score(score_writer, sname, gene_score);
      }
      if (homo_stats[hi].fdr <= threshold) {
        //scan_hits(mutiseqs[hi].seqs[0], ds_motifs[mi], muti_base_probs[hi][0], bst_hit, found_hit);
        //cout<< homo_stats[hi].motif_name << "\t" << homo_stats[hi].homoseq_name << "\t" << homo_stats[hi].pvalue << endl;
        stat_writer<< homo_stats[hi].motif_name << "\t" << homo_stats[hi].homoseq_name << "\t";
        //stat_writer<< homo_stats[hi].score << "\t" << homo_stats[hi].uc_diff << "\t" << homo_stats[hi].nc_pvalue << "\t";
        if (homo_stats[hi].bst_hit.found) {
          stat_writer<< homo_stats[hi].bst_hit.site << "\t" << homo_stats[hi].bst_hit.location << "\t" << homo_stats[hi].bst_hit.score << "\t";
        }
        else {
          stat_writer<< "---" << "\t" << "---" << "\t" << "---" << "\t";
        }
        //stat_writer<< homo_stats[hi].aglib_pvalue << "\t";
        stat_writer<< homo_stats[hi].ccv << "\t" << homo_stats[hi].sd << "\t" << homo_stats[hi].pvalue << "\t" << homo_stats[hi].fdr << endl;
        //print_score(score_writer, homo_stats[hi].homoseq_name, homo_stats[hi].gene_list);
      }
    }


  }
  //cout << "test06" << endl;
  //score_writer.close();
  stat_writer.close();

}

double mixed_studentt_avg(double one, vector<double> &obs, vector<vector<double> > &bkg, double &ccv, double &sd) {
  if (bkg.size() == 0) {
    return 1.0;
    ccv = 0.0;
    sd = 0.0;
  }
  double total_tvalue = 0.0;
  double total_dfvalue = 0.0;
  double total_ccv = 0.0;
  double total_sd = 0.0;
  double tmp_pvalue, tmp_tvalue, tmp_df;
  double tmp_ccv, tmp_sd;
  for (unsigned int gi = 0; gi < bkg.size(); gi++) {
      mixed_ttest(one, obs, bkg[gi], tmp_pvalue, tmp_tvalue, tmp_df, tmp_ccv, tmp_sd);
      total_tvalue += tmp_tvalue;
      total_dfvalue += tmp_df;
      total_ccv += tmp_ccv;
      total_sd += tmp_sd;
  }
  double pvalue, tvalue, dfvalue;
  tvalue = total_tvalue / (double)bkg.size();
  dfvalue = total_dfvalue / (double)bkg.size();

  ccv = total_ccv / (double)bkg.size();
  sd = total_sd / (double)bkg.size();

  double bound;
  double df;
  double p;
  double q;
  int status;
  double t;
  int which = 1;

  cdft ( &which, &p, &q, &tvalue, &dfvalue, &status, &bound );

  if ( status != 0 ) {
    return 1.0;
  } else {
    return p;
  }

}

double average(vector<vector<double> > &bkg) {
  unsigned int total_n = 0;
  double total = 0.0;
  for (unsigned int gi = 0; gi < bkg.size(); gi++) {
    for (unsigned int si = 0; si < bkg[gi].size(); si++) {
      total += bkg[gi][si];
      total_n++;
    }
  }
  if (total_n == 0) {
    return 0.0;
  } else {
    return total / (double)total_n;
  }

}

/*
void raw_score::analysis_motifs_st() {
  srand(int(time(0)));  //important!

  //ofstream score_writer;
  //string score_file = dir + "/" + session + ".score.txt";
  //string score_file = "test.score.txt";
  //score_writer.open(score_file.c_str(), ios::out);
  //if (!score_writer) mcf_die("Sorry, couldn't open " + score_file);
  ofstream stat_writer;
  //string stat_file = dir + "/" + session + ".stat.txt";
  string stat_file = output_file;
  stat_writer.open(stat_file.c_str(), ios::out);
  if (!stat_writer) mcf_die("Sorry, couldn't open " + stat_file);
  stat_writer<< "# motif_name\thomoseq_name\tbinding_site\tlocation\tscore\tdiff\tpvalue\tfdr" << endl;

  vector<pair<string, double> > gene_score;
  vector<pair<string, double> > gene_list;
  //double bg_min, bg_max, bg_total_score;

  bg_base_probs.resize(bgseqs.size());
  for (unsigned s = 0; s < bgseqs.size(); ++s)
    get_base_probs(bgseqs[s], bg_base_probs[s]);

  for (unsigned int mi = 0; mi < ds_motifs.size(); ++mi) {
    //cout << motif_names[mi] << endl;
    //score_writer << "//" << motif_names[mi] <<"\n";

    vector<double> exp_a;
    vector<double> exp_b;
    double score_one;
    //double pvalue = two_sample_ttest(exp_a, exp_b);
    exp_a.clear();
    double tmp_score;
    for (uint s = 0; s < bgseqs.size(); ++s) {
      tmp_score = scan_seq(bgseqs[s], ds_motifs[mi], bg_base_probs[s]);
      if (tmp_score > 0.0) {
        exp_a.push_back(log(tmp_score));
      }
      //exp_a.push_back(scan_seq(bgseqs[s], ds_motifs[mi], bg_base_probs[s]));
    }
    //string sname = "background";
    //print_score(score_writer, sname, exp_a);

    double average_a = average(exp_a);

    vector<stat> homo_stats;
    homo_stats.clear();
    unsigned int mutiseq_size = mutiseqs.size();
    for (unsigned int hi = 0; hi < mutiseqs.size(); ++hi) {
      if (mutiseqs[hi].seqs.size() < 2)
        continue;
      tmp_score = scan_seq(mutiseqs[hi].seqs[0], ds_motifs[mi], muti_base_probs[hi][0]);
      if (tmp_score > 0.0)
        score_one = log(tmp_score);
      else
        continue;
        //score_one = -1.79e308;

      exp_b.clear();
      for (uint s = 1; s < mutiseqs[hi].seqs.size(); ++s) {
        tmp_score = scan_seq(mutiseqs[hi].seqs[s], ds_motifs[mi], muti_base_probs[hi][s]);
        if (tmp_score > 0.0) {
          exp_b.push_back(log(tmp_score));
        }
        //exp_b.push_back(scan_seq(mutiseqs[hi].seqs[s], ds_motifs[mi], muti_base_probs[hi][s]));
      }
      //cout << mutiseqs[hi].seq_names[0] << endl;

      //string sname = "observes";
      //print_score(score_writer, sname, exp_b);

      stat a_stat;
      a_stat.motif_name = motif_names[mi];
      a_stat.homoseq_name = mutiseqs[hi].homo_name;
      //a_stat.aglib_pvalue = two_sample_ttest(exp_a, exp_b);
      if (exp_b.size() == 0) {
        a_stat.pvalue = 1.0;
        //cout << "test01" << endl;
      }
      else {
        a_stat.pvalue = studentt_one(score_one, exp_b, exp_a);
        //a_stat.pvalue = studentt_one_sample(score_one, exp_a);
        //a_stat.pvalue = studentt_one_sample(score_one, exp_b);
        //a_stat.pvalue = studentt_two_sample(exp_b, exp_a);
        //cout << "test02" << endl;
        //a_stat.pvalue = studentt_one(exp_b, exp_a);
        //cout << "test03" << endl;
      }
      //a_stat.nc_pvalue = 1.0 - studentt_one(exp_b[0], exp_b);
      //a_stat.diff = average(exp_b) - average_a;
      a_stat.diff = score_one - average_a;
      //cout << "score_one - average_a: " << score_one << ", " << average_a << endl;
      //a_stat.score = exp_b[0];
      //a_stat.uc_diff = exp_b[0] - average(exp_b);

      bool found_hit;
      scan_hits(mutiseqs[hi].seqs[0], ds_motifs[mi], muti_base_probs[hi][0], a_stat.bst_hit, found_hit);
      if (found_hit) {
        a_stat.bst_hit.found = true;
      } else {
        a_stat.bst_hit.found = false;
      }
      if (a_stat.pvalue <= threshold) {
        //a_stat.gene_list = gene_list;
        //print_score(score_writer, mutiseqs[hi].homo_name, gene_list);
      }
      homo_stats.push_back(a_stat);
    }

    //cout << "test04" << endl;

    sort(homo_stats.begin(), homo_stats.end(), stat_less_p);
    double fdr;
    for (unsigned int hi = 0; hi < mutiseqs.size(); ++hi) {
      fdr = homo_stats[hi].pvalue * mutiseq_size / (hi + 1.0);
      if (fdr <= 1.0) {
        homo_stats[hi].fdr = homo_stats[hi].pvalue * mutiseq_size / (hi + 1.0);
      }
      else {
        homo_stats[hi].fdr = 1.0;
      }
      
    }

    //cout << "test05" << endl;

    sort(homo_stats.begin(), homo_stats.end(), stat_less_fdr);
    for (unsigned int hi = 0; hi < mutiseqs.size(); ++hi) {
      if (homo_stats[0].fdr <= threshold) {
        string sname = "BACKGROUND";
        //print_score(score_writer, sname, gene_score);
      }
      if (homo_stats[hi].fdr <= threshold) {
        //scan_hits(mutiseqs[hi].seqs[0], ds_motifs[mi], muti_base_probs[hi][0], bst_hit, found_hit);
        //cout<< homo_stats[hi].motif_name << "\t" << homo_stats[hi].homoseq_name << "\t" << homo_stats[hi].pvalue << endl;
        stat_writer<< homo_stats[hi].motif_name << "\t" << homo_stats[hi].homoseq_name << "\t";
        //stat_writer<< homo_stats[hi].score << "\t" << homo_stats[hi].uc_diff << "\t" << homo_stats[hi].nc_pvalue << "\t";
        if (homo_stats[hi].bst_hit.found) {
          stat_writer<< homo_stats[hi].bst_hit.site << "\t" << homo_stats[hi].bst_hit.location << "\t" << homo_stats[hi].bst_hit.score << "\t";
        }
        else {
          stat_writer<< "---" << "\t" << "---" << "\t" << "---" << "\t";
        }
        //stat_writer<< homo_stats[hi].aglib_pvalue << "\t";
        stat_writer<< homo_stats[hi].diff << "\t" << homo_stats[hi].pvalue << "\t" << homo_stats[hi].fdr << endl;
        //print_score(score_writer, homo_stats[hi].homoseq_name, homo_stats[hi].gene_list);
      }
    }


  }
  //cout << "test06" << endl;
  //score_writer.close();
  stat_writer.close();

}
*/

void raw_score::print_score(ofstream &writer, string &sname, vector<double> &gene_score) {
	//ofstream writer;
	/*writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file: " << file_name << "\n";
		exit(1);
	}*/
	writer << sname << "\t";
	string one_score = "";
	for (unsigned long gi = 0; gi < gene_score.size(); ++gi) {
    //one_score = gene_score[gi].first + "," + doubleToString(1.0 - gene_score[gi].second, 6) + "; ";
    one_score = doubleToString(gene_score[gi], 6) + ";";
		writer << one_score;
  }
  writer << endl;
	//writer.close();

}

void raw_score::print_score(ofstream &writer, string &sname, vector<pair<string, double> > &gene_score) {
	//ofstream writer;
	/*writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file: " << file_name << "\n";
		exit(1);
	}*/
	writer << sname << "\n";
	string one_score = "";
	for (unsigned long gi = 0; gi < gene_score.size(); ++gi) {
    //one_score = gene_score[gi].first + "," + doubleToString(1.0 - gene_score[gi].second, 6) + "; ";
    one_score = gene_score[gi].first + "," + doubleToString(gene_score[gi].second, 6) + "\n";
		writer << one_score;
  }
  writer << endl;
	//writer.close();

}

/*
void raw_score::analysis_motifs() {
  result myresult;
  for (unsigned int mi = 0; mi < ds_motifs.size(); ++mi) {
    cout << motif_names[mi] << endl;
    for (unsigned int hi = 0; hi < mutiseqs.size(); ++hi) {
      cout << mutiseqs[hi].homo_name << endl;
      myresult.seq_scores.clear();
      get_raw_scores(ds_motifs[mi], mutiseqs[hi].seqs, muti_base_probs[hi], myresult);
    }
    cout << endl;
  }
}
*/

string doubleToString(const double &dbNum, unsigned long dot_n) {
    char *chCode;
    chCode = new(std::nothrow)char[20];
	string code_str = "%." + long_to_string(dot_n) + "lf";  // .6 是控制输出精度的，6位小数
    sprintf(chCode, code_str.c_str(), dbNum);
    string strCode(chCode);
    delete []chCode;
    return strCode;

}

string long_to_string(unsigned long val)
{
    std::stringstream ss;
    ss << val;
    return ss.str();
}
