#include <iostream>
#include <string.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <time.h>
#include <iomanip>
#include <iterator>
#include <cctype>
#include <algorithm>
using namespace std;

/******************************************************************
 * Please noted that we exchange the query and template when user 
 * inputted and print the results 
 * 
 * Please use the following CMD to compile this program.
 * 1. conda create -n gcc-g++ python=3.10
 * 2. conda install -c conda-forge gcc gxx
 * 3. conda activate gcc-g++ && g++ -static -O3 -ffast-math -lm -o TMscoreC TMscoreC.cpp
 ******************************************************************/ 

const char* VERSION = "20260115";

typedef enum {
	PDB,
	CIF,
	AUTO
}FILE_FORMAT;

FILE_FORMAT g_file_format = AUTO;

typedef enum {
	PROTEIN, // 0
	DNA,  // 1
	RNA,  // 2
	LIGAND  //3
} MOLTYPE;

typedef enum {
	ROUGH,
	DETAIL
}PRINT_RESULT_TYPE;

typedef enum {
	FAST,
	NORMAL,
	SLOW,
	SUPERFAST
}ALIGN_TYPE;

typedef struct {
	string qchain;
	string tchain;
}CHAIN_PAIR;

typedef struct {
	string qchain;
	string tchain;
	string qaa;
	string taa;
	int qind;
	int tind;
	int qoind;
	int toind;
	char qoindsuf;
	char toindsuf;
	double dis2;
}ALIGN_PAIR; 

typedef struct {
	string chain1;
	string chain2;
	int ind1;
	int ind2;
}INTERFACE_PAIR; 

bool g_lig_match_in_order_or_content = false; // false means with content

double g_interact_dis_cut_pow2 = 64.;
bool g_fully_cared_ligatom_match = true;
bool g_use_rtmscore_to_search_best = false; // whether use rtmscore to search best superimpose.
bool g_use_chain_order_or_not = false;
bool g_use_res_orig_index_or_not = false; // record whether directly use the residue/nucletide index in pdb.
bool g_is_load_H2O = false; // record whether load the H2O (water) molecule.
double g_seqid_cutoff = 0.7;
PRINT_RESULT_TYPE g_print_result_type = DETAIL; // [DO NOT CHANGE IT] required to be DETAIL
ALIGN_TYPE g_ali_type = NORMAL;
bool g_go_detail = false;

double g_user_given_d0 = -1.; 
double g_user_given_d02 = -1.; 

bool g_is_sup_protein = true;
bool g_is_sup_dna     = true;
bool g_is_sup_rna     = true;
bool g_is_sup_ligand  = true;

bool g_is_inputed_chain_alignment = false;
string g_inputed_chain_alignment;

int g_maxmum_number_of_iterations = 20; // !maximum number of iterations for TM-score-monomer
int g_maxmum_number_of_L_init = 6; // !maximum number of L_init

bool g_is_output_in_detail = false;

string g_atomtype_nuc = " C3'";
string g_atomtype_res = " CA ";

string TOOL_EXE = "TM-score-Complex";

double g_eps = 1e-9;

class CSort{
	public:
		// need release the return pointer
		static int* quickDescendSortIndex(const int& n, const double* arr, const int& left, const int& right){
			int* indexes = new int[n];
			for (int i = 0; i < n; i++)
				indexes[i] = i;
			
			vector<int> stack;
			stack.push_back(left);
			stack.push_back(right);
			
			while (0 != stack.size()) {
				int end = stack[stack.size()-1];
				stack.pop_back();
				int begin = stack[stack.size()-1];
				stack.pop_back();
				
				int keyi = partition_of_quickDescendSortIndex(arr, indexes, begin, end);
				
				if (keyi+1 < end) {
					stack.push_back(keyi+1);
					stack.push_back(end);
				}
				if (begin < keyi-1) {
					stack.push_back(begin);
					stack.push_back(keyi-1);
				}
			}
			
			return indexes;
		}
		
		static int find_max_index(const int& n, const double* arr){
			int max_ind = 0;
			for (int i = 1; i < n; i++){
				if (arr[i] > arr[max_ind]){
					max_ind = i;
				}
			}
			
			return max_ind;
		}
	
	private:
		static int partition_of_quickDescendSortIndex(const double* arr, int* indexes, const int& __left, const int& __right) {
			int left = __left;
			int right = __right;
			int pivot = indexes[left];
	        while (left < right) {
	            while (left < right && arr[indexes[right]] < arr[pivot]) {  
	                right--;  
	            }  
	            indexes[left] = indexes[right];  
	            while (left < right && arr[indexes[left]] >= arr[pivot]) {  
	                left++; 
	            }  
	            indexes[right] = indexes[left];  
	        }
	        
	        indexes[left] = pivot;  
	        return left;  
		}
};


class CBaseFunc {
	public:
		
		static double gdt_ts(const vector<double>& dis2_vec);
		static double gdt_ha(const vector<double>& dis2_vec);
		static double rmsd(const vector<double>& dis2_vec);
		
		static int        sum(const vector<int>& vec){
			int i = 0, n = vec.size();
			int ans = 0;
			for (i = 0; i <n; i++){
				ans += vec[i];
			}
			
			return ans;
		}
		static void       print_logo();
		static void       print_help(const char* arg); 
		static string     stringTrim(const string& str);
		static void       toUpperString(string &str);
		static string     eraseAll(const string &str, char ch);
		static string     eraseAll(const string &str, const char* arr, int len);
		CBaseFunc();
		virtual    ~CBaseFunc();
		
		static char       nucletideAbbr3WordsTo1(const string& nuc);
		static string     aminoAcidAbbr1WordTo3(const char& acid);
		static char       aminoAcidAbbr3WordsTo1(const string& acid);
		static double**   new2Darr(const int& r, const int& c);
		static void       delete2Darr(double** pMtx, const int& r);
		static bool       isInit2Darr(double** mtx, const int& r, const int& c);
		static int**      new2DIntArr(const int& row, const int& col);
		static void       delete2DIntArr(const int& n, int ** Arr);
		static bool*      new1Dbool(const int& row);
		static bool**     new2DBoolArr(const int& row, const int& col);
		static void       delete2DBoolArr(const int& n, bool ** Arr);
		static int***     new3DIntArr(int row, int col, int thd);
		static void       delete3DIntArr(int row, int col, int*** Arr);
		static string     charVec2String(const vector<char>& vec);
		static vector<string> string2stringVec(const string& str);
		static vector<string> stringSplit(const string& str, const char& spit);
		static vector<string> stringSplit(const string& str, const char& spit1, const char& spit2);
		
		static double d0_of_tmscore(const int& Lnorm);
		static double d0_of_lsscore(const int& Lnorm);
		static double d0_of_tmscore_c3prime(const int& Lnorm);
		
		static double greedySearch(double** scoMtx, const int& rownum, const int& colnum, int* alivec, int* transpose_alivec);
		static double __2merExchangedSearch(double** scoMtx, const int& rownum, const int& colnum, int* alivec, int* transpose_alivec, const double& prev_score, const int& iter_max);
		static double __3merExchangedSearch(double** scoMtx, const int& rownum, const int& colnum, int* alivec, int* transpose_alivec, const double& prev_score, const int& iter_max); 
		
		static double cal_rot_tran_from_query_to_templ__(const vector<double*>& query, const vector<double*>& templ, double** out_u, double* out_t, const double& user_d0, const bool& fast);
		static double cal_rot_tran_from_query_to_templ__II(const vector<double*>& query, const vector<double*>& templ, double** out_u, double* out_t, const double& user_d0, const bool& fast);
		static double cal_rot_tran_from_query_to_templ__(const vector<double*>& query, const vector<double*>& templ, double** out_u, double* out_t, const double& user_d0, const int* q2t, const bool& fast);
		static double cal_rot_tran_from_query_to_templ__(const vector<double*>& query, const vector<double*>& templ, const vector<MOLTYPE>& moltypes, double** out_u, double* out_t, const double& user_d0_pro, const double& user_d0_dna,  const double& user_d0_rna, const double& user_d0_lig, const bool& fast);
		static double cal_rot_tran_from_query_to_templ__II(const vector<double*>& query, const vector<double*>& templ, const vector<MOLTYPE>& moltypes, double** out_u, double* out_t, const double& user_d0_pro, const double& user_d0_dna, const double& user_d0_rna, const double& user_d0_lig, const bool& fast); 
		static void score_fun(double* xt, double* yt, double* zt, double* xb, double* yb, double* zb, const double& d2, int& n_cut, int& n_ali, int* i_ali, const double& d02, const int& nseq, double& score);
		static void score_fun(double* xt, double* yt, double* zt, double* xb, double* yb, double* zb, MOLTYPE* mt, 
								const double& d_pro2, const double& d_dna2, const double& d_rna2, const double& d_lig2, 
								int& n_cut, int& n_ali, int* i_ali, 
								const double& d0_pro2, const double& d0_dna2, const double& d0_rna2, const double& d0_lig2, 
								const int& nseq, double& score);
		static double score_fun_once(const vector<double*>& rotted_qxyz, 
								const vector<double*>& txyz, 
								const vector<MOLTYPE>& mts,
								const double& d0_pro2, const double& d0_dna2, const double& d0_rna2, const double& d0_lig2,
								const int& total_res_num);
		
		static double cal_rot_tran_from_query_to_templ__for_rTMscore(
			const int& aligned_chain_num, 
			const int& chain_num, 
			const vector<double*>& query, const vector<double*>& templ, 
			map<int, int>& chain_index_corr_to_query__aa_num,
			map<int, MOLTYPE>& chain_index_corr_to_query__moltype,
            map<int, double>& chain_index_corr_to_query__d0,
            map<int, double>& chain_index_corr_to_query__d02,
			const vector<int>& chain_index_corr_to_query,
			double** out_u, double* out_t, const bool& fast);
		static void score_fun_rtmsco(
			const int& aligned_chain_num,
			const int& chain_num,
			double* xt, double* yt, double* zt, 
            map<int, int>& chain_index_corr_to_query__aa_num,
            map<int, double>& chain_index_corr_to_query__d2,
            map<int, double>& chain_index_corr_to_query__d02,
			int* chain_index_corr_to_query, 
			double* xb, double* yb, double* zb,
			int& n_cut, int& n_ali, int* i_ali, double& score);
		static double score_fun_rtmsco_once(
			const int& aligned_chain_num,
			const int& chain_num,
			const vector<double*>& rotted_qxyz, 
			const vector<double*>& txyz, 
			const vector<int>& chain_index_corr_to_query, 
			map<int, int>& chain_index_corr_to_query__aa_num,
			map<int, double>& chain_index_corr_to_query__d02);
		
		static double u3b_func(const vector<double*>& query, const vector<double*>& templ, const vector<MOLTYPE>& mts, const double& d02_pro, const double& d02_dna, const double& d02_rna, const double& d02_lig);
		static void u3b(double** x, double** y, const int& n, const int& mode, double** u, double* t);		
		static void u3b(const vector<double*>& query, const vector<double*>& templ, double** out_u, double* out_t);
		
		static double* rotateAndTrans(double* xyz, double** u, double* t);
		
		static bool is_same(const vector<char>& a, const vector<char>& b);
		static double rough_score(const MOLTYPE& mt, const vector<double*>& q, const vector<double*>& t, const int* q2t);
		
		static double distance2(const double* a, const double* b);
		static double distance2(const double* q, const double* t, double** qu, double* qt);
		static double fabs(const double& det);
		
		
		static void inverseRotAndTransMtx(double** input_rot_mtx, double* input_trans, double** output_rot_mtx, double* output_trans);
		
		static const vector<CHAIN_PAIR>* parse_matched_chain_pairs(const string& matched_chain_pair_path);
	private:
		static double yuzishiForInverse33Mtx(double** A, int i, int j); 
		static bool inverse33Mtx(double** A, double** invA);
};

class CVector{
	public:
		static bool is_same_contents(const int& alen, int* avec, const int& blen, int* bvec){
			if (alen == blen){
				double ans = true;
				bool* is_used = CBaseFunc::new1Dbool(blen);
				for (int i = 0; i < alen; i++){
					bool is_matched = false; 
					for (int j = 0; j < blen; j++) {
						if (is_used[j]) continue;
						if (avec[i] == bvec[j]){
							is_matched = true;
							is_used[j] = true;
							break;
						}
					}
					
					if (!is_matched){
						ans = false;
						break;
					}
				}
				
				delete[] is_used;
				return ans;
			}else return false;
		}
		
		static bool is_same_contents(const vector<string>& avec, const vector<string>& bvec){
			int alen = avec.size();
			int blen = bvec.size();
			if (alen == blen){
				double ans = true;
				bool* is_used = CBaseFunc::new1Dbool(blen);
				for (int i = 0; i < alen; i++){
					bool is_matched = false; 
					for (int j = 0; j < blen; j++) {
						if (is_used[j]) continue;
						if (avec[i] == bvec[j]){
							is_matched = true;
							is_used[j] = true;
							break;
						}
					}
					
					if (!is_matched){
						ans = false;
						break;
					}
				}
				
				delete[] is_used;
				return ans;
			}else return false;
		}
		
		static bool is_same_order(const vector<string>& avec, const vector<string>& bvec){
			int alen = avec.size();
			int blen = bvec.size();
			if (alen == blen){
				double ans = true;
				for (int i = 0; i < alen; i++){
					if (avec[i] != bvec[i]){
						ans = false;
						break;
					}
				}
				
				return ans;
			}else return false;
		}
		
		static bool is_same_contents(const vector<vector<string> >& avec, const vector<vector<string> >& bvec, const bool& consider_inner_order=false){
			int alen = avec.size();
			int blen = bvec.size();
			if (alen == blen){
				double ans = true;
				bool* is_used = CBaseFunc::new1Dbool(blen);
				for (int i = 0; i < alen; i++){
					bool is_matched = false; 
					for (int j = 0; j < blen; j++) {
						if (is_used[j]) continue;
						if (consider_inner_order? is_same_order(avec[i], bvec[j]) :is_same_contents(avec[i], bvec[j])){
							is_matched = true;
							is_used[j] = true;
							break;
						}
					}
					
					if (!is_matched){
						ans = false;
						break;
					}
				}
				
				delete[] is_used;
				return ans;
			}else return false;
		}
		
		static bool is_same_contents(const vector<vector<string>* >& avec, const vector<vector<string>* >& bvec, const bool& consider_inner_order=false){
			int alen = avec.size();
			int blen = bvec.size();
			if (alen == blen){
				double ans = true;
				bool* is_used = CBaseFunc::new1Dbool(blen);
				for (int i = 0; i < alen; i++){
					bool is_matched = false; 
					for (int j = 0; j < blen; j++) {
						if (is_used[j]) continue;
						if (consider_inner_order? is_same_order(*avec[i], *bvec[j]) :is_same_contents(*avec[i], *bvec[j])){
							is_matched = true;
							is_used[j] = true;
							break;
						}
					}
					
					if (!is_matched){
						ans = false;
						break;
					}
				}
				
				delete[] is_used;
				return ans;
			}else return false;
		}
};

class AtomVarDerWaalRadius{
	private:
		map<string, double> vdwRMap;
	public:
		AtomVarDerWaalRadius(){
			// the var der waal radius information is come from "Van der Waals Radii of Elements" (S. S. Batsanov)
			string varDerWaalRadiusStr = "H   1.2 \nHe 1.4 \nLi 2.2 \nBe 1.9 \nB 1.8 \nC 1.7 \nN 1.6 \nO 1.55 \nF 1.5 \nNe 1.5 \nNa 2.4 \nMg 2.2 \nAl 2.1 \nSi 2.1 \nP 1.95 \nS 1.8 \nCl 1.8 \nAr 1.9 \nK 2.8 \nCa 2.4 \nSc 2.3 \nTi 2.15 \nV 2.05 \nCr 2.05 \nMn 2.05 \nFe 2.05 \nCo 2.0 \nNi 2.0 \nCu 2.0 \nZn 2.1 \nGa 2.1 \nGe 2.1 \nAs 2.05 \nSe 1.9 \nBr 1.9 \nKr 2.0 \nRb 2.9 \nSr 2.55 \nY 2.4 \nZr 2.3 \nNb 2.152 \nMo 2.1 \nTc 2.05 \nRu 2.05 \nRh 2.0 \nPd 2.05 \nAg 2.1 \nCd 2.2 \nIn 2.2 \nSn 2.25 \nSb 2.2 \nTe 2.1 \nI 2.1 \nXe 2.2 \nCs 3.0 \nBa 2.7 \nLa 2.5 \nCe 2.35 \nPr 2.39 \nNd 2.29 \nPm 2.36 \nSm 2.29 \nEu 2.33 \nGd 2.37 \nTb 2.21 \nDy 2.29 \nHo 2.16 \nEr 2.35 \nTm 2.27 \nYb 2.42 \nLu 2.21 \nHf 2.25 \nTa 2.2 \nW 2.1 \nRe 2.05 \nOs 2.0 \nIr 2.0 \nPt 2.05 \nAu 2.1 \nHg 2.05 \nTl 2.2 \nPb 2.3 \nBi 2.3 \nPo 2.29 \nAt 2.36 \nRn 2.43 \nFr 2.56 \nRe 2.43 \nRa 2.43 \nAc 2.6 \nTh 2.4 \nPa 2.43 \nU 2.3 \nNp 2.21 \nPu 2.56 \nAm 2.56 \nCm 2.56 \nBk 2.56 \nCf 2.56 \nEs 2.56 \nFm 2.56 \nMd \nNo \nLr \nRf \nDb \nSg \nBh \nHs \nMt \nDs \nRg \nCn \nUut \nUuq \nUup \nUuh \nUus \nUuo ";
	
			vector<string> atomInfos = CBaseFunc::stringSplit(varDerWaalRadiusStr, '\n');
			for (int i = 0; i < atomInfos.size(); i++){
				vector<string> atomRadii = CBaseFunc::stringSplit(atomInfos[i], ' ');
				if (2 == atomRadii.size()){
					string key(atomRadii[0]);
					CBaseFunc::toUpperString(key);
					
					double value = 0.0;
					sscanf(atomRadii[1].c_str(), "%lf", &value);
					
					vdwRMap[key] = value;
				}
			}
		}
		
		double operator[](const string& atomtype){
			string at = CBaseFunc::eraseAll(atomtype, ' ');
			CBaseFunc::toUpperString(at);
			return vdwRMap[at];
		}
};

class CNWalign {
	private:
		void load_blosum62_for_protein();
		void load_blosumn_for_nucletide();
		void run_needleman_wunsch();
	public:
		static double nwalign(const vector<char>& aseq, const vector<int>& aoind, const vector<char>& bseq, const vector<int>& boind, const MOLTYPE& mol_type, int* out_ali, int& identity_ali_num);
	private:
		CNWalign(const string& aseq, const string& bseq, const MOLTYPE& mol_type);
		virtual ~CNWalign();
		
		const int& get_identity_num();
		const double identity_div_min_len();
		const int& operator [](const int& i); // return ali[i]
		const int* get_identity_ali();
	private:
		const string& aseq;
		const string& bseq;
		const MOLTYPE& mol_type;
		int alen;
		int blen;
		string imut_seq;  // corresponding to imut
		int imut_seq_len;
		int imut[23][23]; // e.g., blosum62
		int gap_open;
		int gap_extn;
		int identity_num;
		int* identity_ali;  // save the aligned information // index start from 0 // it is normal alignment not identity alignment like US-align's option "-TM-score 7"
};

class Molecule {
	private:
		MOLTYPE m_moltype;
		
		vector<string> m_all_info_vec;       // for saving
		vector<double*> m_all_xyz_vec;       // for saving, need release memory 
		vector<double*> m_cared_xyz_vec;     // for calculating, do not need release memory
		vector<int> m_orig_index_vec;        // for calculating
		vector<char> m_char_following_orig_index_vec;
		vector<char> m_seq_vec;				 // for protein, dna, rna, not ligand
		string m_seq_str;
		vector<string> m_cared_atomtype_vec;       // only for ligand
		vector<string> m_cared_atomsimpletype_vec; // only for ligand
		
	public:
		Molecule(const MOLTYPE& moltype, const vector<string>* contents);
		virtual ~Molecule();
		
		const MOLTYPE& get_moltype();
		const bool is_same_molytpe(const Molecule& mol);
		const double* operator [] (const int& i);
		const int size();
		const vector<string> to_str(double** u, double* t);
		const vector<string>& to_str();
		const string get_seq_str();
		const vector<char>& get_seq_vec();
		const vector<double*>& get_cared_xyz_vec();
		const int& get_ith_orig_index(const int& i);
		const char& get_ith_char_following_orig_index_vec(const int& i);
		const vector<int>& get_orig_index_vec();
		
		const string& get_cared_atomsimpletype_in_lig(const int& i);  // only for ligand
		const vector<string>& get_cared_atomtype_vec_in_lig(); //only for ligand
	private:
		void load_one_protein(const vector<string>* contents);
		void load_one_nucletide(const vector<string>* contents);
		void load_one_ligand(const vector<string>* contents);
};

class ShortestRoad{
private:
	double** adjmap;
	int* p_node_num;
	int*** allRoad8Floyd;
	int** allRoadLen8Floyd;

	vector<int> longestShortestRoad8Floyd;
public:
	ShortestRoad(double** adjmap, int node_num){
		this->adjmap = CBaseFunc::new2Darr(node_num, node_num);
		for (int i = 0; i < node_num; i++){
			for (int j = 0; j < node_num; j++){
				this->adjmap[i][j] = adjmap[i][j];
			}
		}
		
		this->p_node_num = new int[1];
		p_node_num[0] = node_num;
		
		allRoad8Floyd = CBaseFunc::new3DIntArr(node_num, node_num, node_num);
		allRoadLen8Floyd = CBaseFunc::new2DIntArr(node_num, node_num);
		
		floyd4LongestLengthRoadWithShortestDis();
	}
	
	int getNodeNum(){
		return *p_node_num;
	}
	
	vector<int> getLongestShortestRoad8Floyd(){
		return longestShortestRoad8Floyd;
	}

	virtual ~ShortestRoad(){
		if (NULL != adjmap)           CBaseFunc::delete2Darr(adjmap, *p_node_num);
		if (NULL != allRoad8Floyd)    CBaseFunc::delete3DIntArr(*p_node_num, *p_node_num, allRoad8Floyd);
		if (NULL != allRoadLen8Floyd) CBaseFunc::delete2DIntArr(*p_node_num, allRoadLen8Floyd);
		delete[] p_node_num;
	}

	int*** getAllRoad8Floyd(){
		return this->allRoad8Floyd;
	}

private: 	
	void floyd4LongestLengthRoadWithShortestDis(){
		int i = 0;
		int maxRoadLen = 0;
		double corrRoadDis = 0;

		double** disMtx = CBaseFunc::new2Darr(*p_node_num, *p_node_num);
		int** spot = CBaseFunc::new2DIntArr(*p_node_num, *p_node_num);
		int* maxLenRoadInfo = new int[*p_node_num]; // road info is "maxLenRoadInfo[0] -> ... -> maxLenRoadInfo[i] -> maxLenRoadInfo[i+1] -> ... -> maxLenRoadInfo[N]"
		
		for (i = 0; i < *p_node_num; i++) {
			for (int j = 0; j < *p_node_num; j++) {
				spot[i][j] = -1;
				disMtx[i][j] = adjmap[i][j];
				
				if (disMtx[i][j] <= 0)
					disMtx[i][j] = 999999999999999;
				if (i == j)
					disMtx[i][j] = 0;
			}
		}

		for (int k = 0; k < *p_node_num; k++){
			for (int i = 0; i < *p_node_num; i++){
				for (int j = 0; j < *p_node_num; j++){
					if (disMtx[i][j] > disMtx[i][k] + disMtx[k][j]) {
						disMtx[i][j] = disMtx[i][k] + disMtx[k][j];
						spot[i][j] = k;
					}
				}
			}
		}
	
		// search the longest length road with the longest shortest distance
		for (i = 0; i < *p_node_num; i++){
			for (int j = 0; j < *p_node_num; j++){
				int* roadLen = new int[1];
				int* tmpPath = new int[*p_node_num];
				
				for (int k = 0; k < *p_node_num; k++)
					tmpPath[k] = -1;
				roadLen[0] = 0;
				tmpPath[roadLen[0]++] = i;
				recursionSearchPath8FloydAlgorithm(spot, i, j, tmpPath, roadLen);
				double dis = disMtx[i][j];
				if (roadLen[0] > maxRoadLen){
					maxRoadLen = roadLen[0];
					for (int k = 0; k < *p_node_num; k++)
						maxLenRoadInfo[k] = tmpPath[k];
					corrRoadDis = dis;
				}else if (roadLen[0] == maxRoadLen){
					if (dis > corrRoadDis){
						maxRoadLen = roadLen[0];
						for (int k = 0; k < *p_node_num; k++)
							maxLenRoadInfo[k] = tmpPath[k];
						corrRoadDis = dis;
					}
				}
				for (int k = 0; k < *p_node_num; k++)
					allRoad8Floyd[i][j][k] = tmpPath[k];
				allRoadLen8Floyd[i][j] = roadLen[0];
				
				delete[] tmpPath;
				delete[] roadLen;
			}
		}
	
		for (i = 0; i < *p_node_num; i++){
			if (-1 == maxLenRoadInfo[i])
				break;
			longestShortestRoad8Floyd.push_back(maxLenRoadInfo[i]);
		}
	
		delete[] maxLenRoadInfo;
		CBaseFunc::delete2DIntArr(*p_node_num, spot);
		CBaseFunc::delete2Darr(disMtx, *p_node_num);
	}
	
	/****************************************************************************
	 * @param spot : the i to j throw node index, floyd algorithm spot matrix
	 * @param i : the 1th position
	 * @param j : the 2nd position
	 * @param onePath : record the shortest road from i, j
	 * @param roadLen : roadLen[0] record the road length
	 ****************************************************************************/
	void recursionSearchPath8FloydAlgorithm(int** spot, int i, int j, int* path, int* roadLen) {
		if (i == j) return;
		if (spot[i][j] == -1)
			path[roadLen[0]++] = j;
		else {
			recursionSearchPath8FloydAlgorithm(spot, i, spot[i][j], path, roadLen);
			recursionSearchPath8FloydAlgorithm(spot, spot[i][j], j, path, roadLen);
		}
	}
	
};

class LigAtomMatcher {
	private:
		int atnum;
		Molecule& lig;
		bool** is_same_imp_mtx;
		vector<vector<int>* > atgs; // atom groups
		vector<vector<vector<string>* >* > attypes_in_roads;
	public:
		LigAtomMatcher(Molecule* p_lig);
		virtual ~LigAtomMatcher();
		const int& size();
		const bool& is_same_important(const int& i, const int& j);
		const bool is_same_important(const int& i, const LigAtomMatcher& other, const int& other_i);
		
		const bool match(LigAtomMatcher& other);
		static double quick_identity_atom_align(Molecule* query, Molecule* templ, double** in_u, double* in_t, int* out_q2t, const double& d02);
		static double quick_identity_atom_align(Molecule* query, Molecule* templ, int* out_q2t, double** out_u, double* out_t, const double& d02);
		static double quick_identity_atom_align(LigAtomMatcher& querym, LigAtomMatcher& templm, double** in_u, double* in_t, int* out_q2t, const double& d02);
		static double quick_identity_atom_align(LigAtomMatcher& querym, LigAtomMatcher& templm, int* out_q2t, double** out_u, double* out_t, const double& d02);
	private:
		bool is_same_roads(const vector<vector<string>* >& ar, const vector<vector<string>* >& br);
		vector<vector<string>* >* extract_ith_roads(int** ar, const string* attypes);
		void extract_atgs();
};


inline const bool LigAtomMatcher::match(LigAtomMatcher& other){
	vector<vector<int>* >& qatgs = this->atgs;
	vector<vector<int>* >& tatgs = other.atgs;
	
	int mx = qatgs.size();
	int nx = tatgs.size();
	
	bool ans = true;
	if (mx == nx){
		bool* is_used = CBaseFunc::new1Dbool(nx);
		for (int i = 0; i < mx; i++){
			vector<int>& ith = *qatgs[i];
			int matched_j = -1;
			
			for (int j = 0; j < nx; j++){
				vector<int>& jth = *tatgs[j];
				if (is_used[j]) continue;
				if (is_same_important(ith[0], other, jth[0])){
					is_used[j] = true;
					matched_j = j;
					break;
				}
			}
			if (matched_j == -1 || qatgs[i]->size() != tatgs[matched_j]->size()){
				ans = false;
				break;
			}
		}
		delete[] is_used;
	}else{
		ans = false;
	}
	
	return ans;
}

inline double LigAtomMatcher::quick_identity_atom_align(Molecule* query, Molecule* templ, double** in_u, double* in_t, int* out_q2t, const double& d02){
	LigAtomMatcher qobj(query);
	LigAtomMatcher tobj(templ);
	
	vector<vector<int>* >& qatgs = qobj.atgs;
	vector<vector<int>* >& tatgs = tobj.atgs;
	
	int mx = qatgs.size();
	int nx = tatgs.size();
	
	int m = query->size();
	int n = templ->size();
	
	bool** is_matched = CBaseFunc::new2DBoolArr(m, n);
	bool* is_used = CBaseFunc::new1Dbool(nx);
	for (int i = 0; i < mx; i++){
		vector<int>& ith = *qatgs[i];
		int nn = ith.size();
		int matched_j = -1;
		for (int j = 0; j < nx; j++){
			if (is_used[j]) continue;
			vector<int>& jth = *tatgs[j];
			if (qobj.is_same_important(ith[0], tobj, jth[0])){
				int mm = jth.size();
				for (int k = 0; k < nn; k++){
					for (int l = 0; l < mm; l++){
						is_matched[ith[k]][jth[l]] = true;
					}
				}
				is_used[j] = true;
				matched_j = j;
				break;
			}
		}
		if (matched_j == -1){
			cout << "* WARNING : there are two ligands with same name but different structure topologies." << endl;
		}
	}
	
	double** scomtx = CBaseFunc::new2Darr(m, n);
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			double dis2 = CBaseFunc::distance2((*query)[i], (*templ)[j], in_u, in_t);
			if (is_matched[i][j])
				scomtx[i][j] = 1./ (1. + dis2/d02);
			else scomtx[i][j] = 0.;
		}
	}
	
	int* t2q = new int[n];
	double gs_score = CBaseFunc::greedySearch(scomtx, m, n, out_q2t, t2q);
	gs_score = CBaseFunc::__2merExchangedSearch(scomtx, m, n, out_q2t, t2q, gs_score, m>n?m:n);
	
	delete[] t2q;
	CBaseFunc::delete2DBoolArr(m, is_matched);
	CBaseFunc::delete2Darr(scomtx, m);
	delete[] is_used;
	
	return gs_score;
}

inline double LigAtomMatcher::quick_identity_atom_align(Molecule* query, Molecule* templ, int* out_q2t, double** out_u, double* out_t, const double& d02){
 	LigAtomMatcher qobj(query);
	LigAtomMatcher tobj(templ);
	
	vector<double*> qxyz_atgs;
	vector<double*> txyz_atgs;
	
	vector<vector<int>* >& qatgs = qobj.atgs;
	vector<vector<int>* >& tatgs = tobj.atgs;
	
	int mx = qatgs.size();
	int nx = tatgs.size();
	
	for (int i = 0; i < mx; i++){
		double* xyz = new double[3];
		xyz[0] = 0.;
		xyz[1] = 0.;
		xyz[1] = 0.;
		
		vector<int>& atg = *qatgs[i];
		int mm = atg.size();
		for (int j = 0; j < mm; j++){
			const double* qxyz = (*query)[atg[j]];
			xyz[0] += qxyz[0];
			xyz[1] += qxyz[1];
			xyz[2] += qxyz[2];
		}
		if (mm > 1){
			xyz[0] /= mm;
			xyz[1] /= mm;
			xyz[2] /= mm;
		}
		
		qxyz_atgs.push_back(xyz);	
	}
	
	for (int i = 0; i < nx; i++){
		double* xyz = new double[3];
		xyz[0] = 0.;
		xyz[1] = 0.;
		xyz[1] = 0.;
		
		vector<int>& atg = *tatgs[i];
		int mm = atg.size();
		for (int j = 0; j < mm; j++){
			const double* txyz = (*templ)[atg[j]];
			xyz[0] += txyz[0];
			xyz[1] += txyz[1];
			xyz[2] += txyz[2];
		}
		if (mm > 1){
			xyz[0] /= mm;
			xyz[1] /= mm;
			xyz[2] /= mm;
		}
		
		txyz_atgs.push_back(xyz);
	}
	
	int m = query->size();
	int n = templ->size();
	
	bool** is_matched = CBaseFunc::new2DBoolArr(m, n);
	vector<double*> matched_txyz_atgs;
	vector<double*> corr_qxyz_atgs;
	bool* is_used = CBaseFunc::new1Dbool(nx);
	for (int i = 0; i < mx; i++){
		vector<int>& ith = *qatgs[i];
		int nn = ith.size();
		int matched_j = -1;
		for (int j = 0; j < nx; j++){
			if (is_used[j]) continue;
			vector<int>& jth = *tatgs[j];
			if (qobj.is_same_important(ith[0], tobj, jth[0])){
				int mm = jth.size();
				for (int k = 0; k < nn; k++){
					for (int l = 0; l < mm; l++){
						is_matched[ith[k]][jth[l]] = true;
					}
				}
				is_used[j] = true;
				matched_j = j;
				break;
			}
		}
		if (matched_j == -1){
			cout << "* WARNING : there are two ligands with same name but different structure topologies." << endl;
		}else{
			corr_qxyz_atgs.push_back(qxyz_atgs[i]);
			matched_txyz_atgs.push_back(txyz_atgs[matched_j]);
		}
	}
	
	CBaseFunc::u3b(corr_qxyz_atgs, matched_txyz_atgs, out_u, out_t);
	
	double** scomtx = CBaseFunc::new2Darr(m, n);
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			double dis2 = CBaseFunc::distance2((*query)[i], (*templ)[j], out_u, out_t);
			if (is_matched[i][j])
				scomtx[i][j] = 1./ (1. + dis2/d02);
			else scomtx[i][j] = 0.;
		}
	}
	
	int* t2q = new int[n];
	double gs_score = CBaseFunc::greedySearch(scomtx, m, n, out_q2t, t2q);
	gs_score = CBaseFunc::__2merExchangedSearch(scomtx, m, n, out_q2t, t2q, gs_score, m>n?m:n);
	
	delete[] t2q;
	CBaseFunc::delete2DBoolArr(m, is_matched);
	CBaseFunc::delete2Darr(scomtx, m);
	delete[] is_used;
	for (int i = 0; i < mx; i++)
		delete[] qxyz_atgs[i];
	for (int i = 0; i < nx; i++)
		delete[] txyz_atgs[i];
	
	return gs_score;
}

inline double LigAtomMatcher::quick_identity_atom_align(LigAtomMatcher& querym, LigAtomMatcher& templm, double** in_u, double* in_t, int* out_q2t, const double& d02){
	LigAtomMatcher& qobj = querym;
	LigAtomMatcher& tobj = templm;
	Molecule* query = &querym.lig;
	Molecule* templ = &templm.lig;
	
	vector<vector<int>* >& qatgs = qobj.atgs;
	vector<vector<int>* >& tatgs = tobj.atgs;
	
	int mx = qatgs.size();
	int nx = tatgs.size();
	
	int m = query->size();
	int n = templ->size();
	
	bool** is_matched = CBaseFunc::new2DBoolArr(m, n);
	bool* is_used = CBaseFunc::new1Dbool(nx);
	for (int i = 0; i < mx; i++){
		vector<int>& ith = *qatgs[i];
		int nn = ith.size();
		int matched_j = -1;
		for (int j = 0; j < nx; j++){
			if (is_used[j]) continue;
			vector<int>& jth = *tatgs[j];
			if (qobj.is_same_important(ith[0], tobj, jth[0])){
				int mm = jth.size();
				for (int k = 0; k < nn; k++){
					for (int l = 0; l < mm; l++){
						is_matched[ith[k]][jth[l]] = true;
					}
				}
				is_used[j] = true;
				matched_j = j;
				break;
			}
		}
		if (matched_j == -1){
			cout << "* WARNING : there are two ligands with same name but different structure topologies." << endl;
		}
	}
	
	double** scomtx = CBaseFunc::new2Darr(m, n);
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			double dis2 = CBaseFunc::distance2((*query)[i], (*templ)[j], in_u, in_t);
			if (is_matched[i][j])
				scomtx[i][j] = 1./ (1. + dis2/d02);
			else scomtx[i][j] = 0.;
		}
	}
	
	int* t2q = new int[n];
	double gs_score = CBaseFunc::greedySearch(scomtx, m, n, out_q2t, t2q);
	gs_score = CBaseFunc::__2merExchangedSearch(scomtx, m, n, out_q2t, t2q, gs_score, m>n?m:n);
	
	delete[] t2q;
	CBaseFunc::delete2DBoolArr(m, is_matched);
	CBaseFunc::delete2Darr(scomtx, m);
	delete[] is_used;
	
	return gs_score;
}

inline double LigAtomMatcher::quick_identity_atom_align(LigAtomMatcher& querym, LigAtomMatcher& templm, int* out_q2t, double** out_u, double* out_t, const double& d02){
 	LigAtomMatcher& qobj = querym;
	LigAtomMatcher& tobj = templm;
	Molecule* query = &querym.lig;
	Molecule* templ = &templm.lig;
	
	vector<double*> qxyz_atgs;
	vector<double*> txyz_atgs;
	
	vector<vector<int>* >& qatgs = qobj.atgs;
	vector<vector<int>* >& tatgs = tobj.atgs;
	
	int mx = qatgs.size();
	int nx = tatgs.size();
	
	for (int i = 0; i < mx; i++){
		double* xyz = new double[3];
		xyz[0] = 0.;
		xyz[1] = 0.;
		xyz[2] = 0.;
		
		vector<int>& atg = *qatgs[i];
		int mm = atg.size();
		for (int j = 0; j < mm; j++){
			const double* qxyz = (*query)[atg[j]];
			xyz[0] += qxyz[0];
			xyz[1] += qxyz[1];
			xyz[2] += qxyz[2];
		}
		if (mm > 1){
			xyz[0] /= mm;
			xyz[1] /= mm;
			xyz[2] /= mm;
		}
		
		qxyz_atgs.push_back(xyz);	
	}
	
	for (int i = 0; i < nx; i++){
		double* xyz = new double[3];
		xyz[0] = 0.;
		xyz[1] = 0.;
		xyz[2] = 0.;
		
		vector<int>& atg = *tatgs[i];
		int mm = atg.size();
		for (int j = 0; j < mm; j++){
			const double* txyz = (*templ)[atg[j]];
			xyz[0] += txyz[0];
			xyz[1] += txyz[1];
			xyz[2] += txyz[2];
		}
		if (mm > 1){
			xyz[0] /= mm;
			xyz[1] /= mm;
			xyz[2] /= mm;
		}
		
		txyz_atgs.push_back(xyz);
	}
	
	int m = query->size();
	int n = templ->size();
	
	bool** is_matched = CBaseFunc::new2DBoolArr(m, n);
	vector<double*> matched_txyz_atgs;
	vector<double*> corr_qxyz_atgs;
	bool* is_used = CBaseFunc::new1Dbool(nx);
	for (int i = 0; i < mx; i++){
		vector<int>& ith = *qatgs[i];
		int nn = ith.size();
		int matched_j = -1;
		for (int j = 0; j < nx; j++){
			if (is_used[j]) continue;
			vector<int>& jth = *tatgs[j];
			if (qobj.is_same_important(ith[0], tobj, jth[0])){
				int mm = jth.size();
				for (int k = 0; k < nn; k++){
					for (int l = 0; l < mm; l++){
						is_matched[ith[k]][jth[l]] = true;
					}
				}
				is_used[j] = true;
				matched_j = j;
				break;
			}
		}
		if (matched_j == -1){
			cout << "* WARNING : there are two ligands with same name but different structure topologies." << endl;
		}else{
			corr_qxyz_atgs.push_back(qxyz_atgs[i]);
			matched_txyz_atgs.push_back(txyz_atgs[matched_j]);
		}
	}
	
	CBaseFunc::u3b(corr_qxyz_atgs, matched_txyz_atgs, out_u, out_t);
	
	double** scomtx = CBaseFunc::new2Darr(m, n);
	for (int i = 0; i < m; i++){
		for (int j = 0; j < n; j++){
			double dis2 = CBaseFunc::distance2((*query)[i], (*templ)[j], out_u, out_t);
			if (is_matched[i][j])
				scomtx[i][j] = 1./ (1. + dis2/d02);
			else scomtx[i][j] = 0.;
		}
	}
	
	int* t2q = new int[n];
	double gs_score = CBaseFunc::greedySearch(scomtx, m, n, out_q2t, t2q);
	gs_score = CBaseFunc::__2merExchangedSearch(scomtx, m, n, out_q2t, t2q, gs_score, m>n?m:n);
	
	delete[] t2q;
	CBaseFunc::delete2DBoolArr(m, is_matched);
	CBaseFunc::delete2Darr(scomtx, m);
	delete[] is_used;
	for (int i = 0; i < mx; i++)
		delete[] qxyz_atgs[i];
	for (int i = 0; i < nx; i++)
		delete[] txyz_atgs[i];
	
	return gs_score;
}


inline const int& LigAtomMatcher::size(){
	return atnum;
}

class Complex {
	private:
		vector<Molecule*> m_mols;
		vector<string> m_chains;
		vector<int> m_avaid_inds;
	public:
		Complex(const string& pdb); 
		const int size();
		const int all_size();
		Molecule* operator [](const int& i);
		Molecule* get_ith_mol_in_all(const int& i);
		Molecule* get_mol_base_on_chain(const string& chain);
		virtual ~Complex();
		
		const string& get_chain(const int& i);
		const vector<string> to_str();
		const vector<string> to_str(double** u, double* t);
		
		void save(const string& path);
		void save(const string& path, double** u, double* t);
		
		// for interface tm-score
		vector<INTERFACE_PAIR> parse_interface_pair(const double& interact_dis_cut_pow2);
	private:
		void load(const string& path);
		bool is_cif_func(const string& path);
		
		void load_pdb(const string& pdb);
		const MOLTYPE parse_atom_mol_type(const string& line);
		
		void load_cif(const string& cif);
		string transfer_ATOM_line_from_cif_to_pdb(const string &cifLine);
		string transfer_HETATM_line_from_cif_to_pdb(const string &cifLine);
};

inline vector<INTERFACE_PAIR> Complex::parse_interface_pair(const double& interact_dis_cut_pow2){
	vector<INTERFACE_PAIR> ans;
	int size = this->size();
	for (int i = 0; i < size; i++){
		Molecule* imol = (*this)[i];
		const string ichain = this->get_chain(i);
		int isize = imol->size();
		for (int j = i+1; j < size; j++){
			Molecule* jmol = (*this)[j];	
			const string jchain = this->get_chain(j);
			int jsize = jmol->size();
			
			for (int k = 0; k < isize; k++){
				const double* kxyz = (*imol)[k];
				for (int l = 0; l < jsize; l++){
					const double* lxyz = (*jmol)[l];
					double dis2 = CBaseFunc::distance2(kxyz, lxyz);
					if (dis2 <= interact_dis_cut_pow2){
						INTERFACE_PAIR ip;
						ip.chain1 = ichain;
						ip.chain2 = jchain;
						ip.ind1 = k;
						ip.ind2 = l;
						ans.push_back(ip);
					}
				}
			}	
		}
	}
	
	return ans;
}

class CTMscoreComplex {
	private:
		Complex* query;
		Complex* templ;
		int qsize;
		int tsize;
		
		double** individual_tmscore_mtx;
		
		double tmscore;
		double rtmscore;  // just a whole level, it is not ok for the cases with miss matched native chains 
		double use_seconds;
		double** u;
		double* t;
		double** inv_u;
		double* inv_t; 
		int* obj_level_ali;
		vector<ALIGN_PAIR> aa_level_ali;
		
		int*** qt_match_mtx;
		LigAtomMatcher** q_ligAtomMatch_obj_vec;
		LigAtomMatcher** t_ligAtomMatch_obj_vec;
		
		map<int, int> chain_index_corr_to_query__aa_num;
		map<int, int> chain_index_corr_to_templ__aa_num;
	    map<int, MOLTYPE> chain_index_corr_to_query__moltype;
	    map<int, double> chain_index_corr_to_query__d0;
	    map<int, double> chain_index_corr_to_query__d02;
		
		int total_res_in_pro_num;
		int total_nuc_in_dna_num;
		int total_nuc_in_rna_num;
		int total_atm_in_lig_num;
		int total_res_num;
		double d0_pro;
		double d0_dna;
		double d0_rna;
		double d0_lig;
		double d02_pro;
		double d02_dna;
		double d02_rna;
		double d02_lig;
	public:
		CTMscoreComplex(const string& qpdb, const string& tpdb, const ALIGN_TYPE& ali_type, const vector<CHAIN_PAIR>* user_chain_pair);
		void set_obj_level_ali(const vector<CHAIN_PAIR>& chainpair);
		const double& get_use_seconds();
		const double& get_tmscore();
		double** get_u();
		const double* get_t();
		const int* get_obj_level_ali();
		void print_result();
		void save_invUt(const string& savepath);
		void save_superposition_ditances(const string& savepath);
		void save_roted_query(const string& savepath);
		void save_roted_templ(const string& savepath);
		void save_sup_pdb_for_web(const string& savepath);
		virtual ~CTMscoreComplex();
	private:
		void align_monomer(const bool& fast);
		
		void align_multimer_fast_buf_inaccuracy_using_nwalign_and_greadsearch(const bool& fast);
		void align_multimer_normal_using_nwalign_and_greadsearch(const bool& fast);
		void align_multimer_slow_but_accuracy_using_nwalign_and_greadsearch(const bool& fast);
		void align_multimer_normal_using_nwalign_and_not_greadsearch(const bool& fast);
		void align_multimer_normal_using_not_nwalign_and_greadsearch(const bool& fast);
		void align_multimer_normal_using_not_nwalign_and_not_greadsearch(const bool& fast);
		
		void align_multimer_fast_buf_inaccuracy_using_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast);
		void align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast);
		void align_multimer_slow_but_accuracy_using_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast);
		void align_multimer_normal_using_nwalign_and_not_greadsearch_with_max_rtmscore(const bool& fast);
		void align_multimer_normal_using_not_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast);
		void align_multimer_normal_using_not_nwalign_and_not_greadsearch_with_max_rtmscore(const bool& fast);
		
		void align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand(const bool& fast);
		void align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand_II(const bool& fast);
		void align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand_III(const bool& fast);
		void align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore_fullycared_ligand(const bool& fast);
		void align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore_fullycared_ligand_III(const bool& fast);
		
		vector<string> formatRotMtx();
		vector<string> formatInvRotMtx();
		void invUt();
		
		void prepare_for_monomer(); 
		void prepare_for_multimer(const bool& use_atom_or_residue_index_order_or_not, const bool& use_chain_order_or_not);
		void generate_atom_alignment_of_each_ligand_pair_using_greedysearch(const bool& use_chain_order_or_not);
		void generate_atom_alignment_of_each_ligand_pair_using_index_order(const bool& use_chain_order_or_not);
		void generate_residue_alignment_of_each_molecule_pair_using_nwalign(const bool& use_chain_order_or_not);
		void generate_residue_alignment_of_each_molecule_pair_resdiue_index(const bool& use_chain_order_or_not);
		int* get_ij_qt_match_mtx(const int& i, const int& j);
		
		// for interface tm-score
		double calcluate_itmscore(const double& interact_dis_cut_pow2);
};

inline double CTMscoreComplex::calcluate_itmscore(const double& interact_dis_cut_pow2){
	vector<INTERFACE_PAIR> interface_pair_in_query = query->parse_interface_pair(interact_dis_cut_pow2);
	if (0 == interface_pair_in_query.size())
		return -1;
	
	map<string, double> chain_d02;
	for (int i = 0; i < qsize; i++){
		if (-1 != obj_level_ali[i]){
			double d02 = this->chain_index_corr_to_query__d02[i];
			chain_d02[query->get_chain(i)] = d02;
		}
	}
	
	// chain1-chain2, iTM-score
	map<string, double> iTM;
	// chain1-chain2, number of interacted residues
	map<string, int> interface_size;
	
	int n = this->aa_level_ali.size();
	int m = interface_pair_in_query.size();
	
	for (int i = 0; i < m; i++){
		INTERFACE_PAIR ip = interface_pair_in_query[i];
		string qchain1 = ip.chain1;
		string qchain2 = ip.chain2;
		int qind1 = ip.ind1;
		int qind2 = ip.ind2;
		
		string key = qchain2+"-"+qchain1;
		if (iTM.end() == iTM.find(key))
			key = qchain1+"-"+qchain2;
		
		if (interface_size.end() == interface_size.find(key))
			interface_size[key] = 1;
		else interface_size[key] = interface_size[key] + 1;
		
		double d02_1 = chain_d02[qchain1];
		double d02_2 = chain_d02[qchain2];
		
		bool is_dis2_1_ok = false;
		double dis2_1 = 0;
		string tchain_1;
		int tind_1 = -1;
		bool is_dis2_2_ok = false;
		double dis2_2 = 0;
		string tchain_2;
		int tind_2 = -1;
		
		for (int j = 0; j < n; j++){
			if (is_dis2_1_ok && is_dis2_2_ok)
				break;
			
			ALIGN_PAIR ap = this->aa_level_ali[j];
			string qchain = ap.qchain;
			int qind = ap.qind;
			double dis2 = ap.dis2;
		
			if (!is_dis2_1_ok && qchain == qchain1 && qind == qind1){
				dis2_1 = dis2;
				is_dis2_1_ok = true;
				tchain_1 = ap.tchain;
				tind_1 = ap.tind;
			}
					
			if (!is_dis2_2_ok && qchain == qchain2 && qind == qind2){
				dis2_2 = dis2;
				is_dis2_2_ok = true;
				tchain_2 = ap.tchain;
				tind_2 = ap.tind;
			}
		}
		
		if (is_dis2_1_ok && is_dis2_2_ok){
			Molecule* tchain1 = templ->get_mol_base_on_chain(tchain_1);
			Molecule* tchain2 = templ->get_mol_base_on_chain(tchain_2);
			double dis2_in_templ = CBaseFunc::distance2((*tchain1)[tind_1], (*tchain2)[tind_2]);
			if (dis2_in_templ <= interact_dis_cut_pow2){
				double tm1 = 1./(1. + dis2_1/d02_1);
				double tm2 = 1./(1. + dis2_2/d02_2);
				double tm = 2.*tm1*tm2/(tm1+tm2);
				
				if (iTM.end() == iTM.find(key))
					iTM[key] = tm;
				else iTM[key] = iTM[key] + tm;
			}
		}
	}
	
//	int aligned_ninterface = 0;
//	int ninterface = 0;
//	double riTM = 0;
//	for (map<string, int>::iterator iter = interface_size.begin(); iter != interface_size.end(); iter++){
//		string key = iter->first;
//		int isize = iter->second;
//		ninterface++;
//		
//		map<string, double>::iterator __iter = iTM.find(key); 
//		if (iTM.end() != __iter){
//			double iTMval = __iter->second;
//			iTMval /= isize;
//			
//			riTM += 1. / iTMval;
//			aligned_ninterface++;
//		}
//	}
//	
//	return aligned_ninterface * aligned_ninterface / (ninterface * riTM); 

	int ninterface = 0;
	double iTM_final = 0;
	for (map<string, int>::iterator iter = interface_size.begin(); iter != interface_size.end(); iter++){
		string key = iter->first;
		int isize = iter->second;
		ninterface++;
		
		map<string, double>::iterator __iter = iTM.find(key); 
		if (iTM.end() != __iter){
			double iTMval = __iter->second;
			iTMval /= isize;
			
			iTM_final += iTMval;
		}
	}
	
	return iTM_final / ninterface; 
}

int main(int argc, char** args) {
	TOOL_EXE = args[0];
	
	int i, j, pdbFileInd = 0;
	
	string query_pdb;  // Structure 2
	string templ_pdb;  // Structure 1
	bool isQueryProvide = false;
	bool isTemplProvide = false;
  	bool isSaveRotTemplInPDB = false;
  	string saveRotTemplInPDBPath;
  	bool isSaveUt = false;
	string saveUtPath;
	bool isSaveSuperpostion = false;
    string saveSuperpostionPath;
    bool isSaveWebNeedSuperposition = false;
	string saveWebNeedSuperpositionPath; 
  	
  	//--------------------------------------------------------------------------//
	//----------------          Load Input Parameters      ---------------------//
	//--------------------------------------------------------------------------//
	for (i = 1; i < argc; i++){
		if (0 == strcmp(args[i], "-mol") && i < argc-1){
			if (0 == strcmp(args[i+1], "prt")){
			    g_is_sup_protein = true;
				g_is_sup_dna = false;
				g_is_sup_rna = false;
				g_is_sup_ligand = false;
			}else if (0 == strcmp(args[i+1], "dna")){
				g_is_sup_protein = false;
				g_is_sup_dna = true;
				g_is_sup_rna = false;
				g_is_sup_ligand = false;
			}else if (0 == strcmp(args[i+1], "rna")){
				g_is_sup_protein = false;
				g_is_sup_dna = false;
				g_is_sup_rna = true;
				g_is_sup_ligand = false;
			}else if (0 == strcmp(args[i+1], "lig")){
				g_is_sup_protein = false;
				g_is_sup_dna = false;
				g_is_sup_rna = false;
				g_is_sup_ligand = true;
			}else if (0 == strcmp(args[i+1], "p+d")){
				g_is_sup_protein = true;
				g_is_sup_dna = true;
				g_is_sup_rna = false;
				g_is_sup_ligand = false;
			}else if (0 == strcmp(args[i+1], "p+r")){
				g_is_sup_protein = true;
				g_is_sup_dna = false;
				g_is_sup_rna = true;
				g_is_sup_ligand = false;
			}else if (0 == strcmp(args[i+1], "p+l")){
				g_is_sup_protein = true;
				g_is_sup_dna = false;
				g_is_sup_rna = false;
				g_is_sup_ligand = true;
			}else if (0 == strcmp(args[i+1], "d+r")){
				g_is_sup_protein = false;
				g_is_sup_dna = true;
				g_is_sup_rna = true;
				g_is_sup_ligand = false;
			}else if (0 == strcmp(args[i+1], "d+l")){
				g_is_sup_protein = false;
				g_is_sup_dna = true;
				g_is_sup_rna = false;
				g_is_sup_ligand = true;
			}else if (0 == strcmp(args[i+1], "r+l")){
				g_is_sup_protein = false;
				g_is_sup_dna = false;
				g_is_sup_rna = true;
				g_is_sup_ligand = true;
			}else if (0 == strcmp(args[i+1], "pdr")){
				g_is_sup_protein = true;
				g_is_sup_dna = true;
				g_is_sup_rna = true;
				g_is_sup_ligand = false;
			}else if (0 == strcmp(args[i+1], "pdl")){
				g_is_sup_protein = true;
				g_is_sup_dna = true;
				g_is_sup_rna = false;
				g_is_sup_ligand = true;
			}else if (0 == strcmp(args[i+1], "prl")){
				g_is_sup_protein = true;
				g_is_sup_dna = false;
				g_is_sup_rna = true;
				g_is_sup_ligand = true;
			}else if (0 == strcmp(args[i+1], "drl")){
				g_is_sup_protein = false;
				g_is_sup_dna = true;
				g_is_sup_rna = true;
				g_is_sup_ligand = true;
			}else{ // all
				g_is_sup_protein = true;
				g_is_sup_dna = true;
				g_is_sup_rna = true;
				g_is_sup_ligand = true;
			}
			i++;
		}else if (0 == strcmp(args[i], "-ffm") && i < argc-1){
			if (0 == strcmp(args[i+1], "pdb")){
			    g_file_format = PDB;
			}else if (0 == strcmp(args[i+1], "mmcif")){
				g_file_format = CIF;
			}
		}else if (0 == strcmp(args[i], "-s") && i < argc-1){
			if (0 == strcmp(args[i+1], "t")){
				g_use_rtmscore_to_search_best = false;
			}else{
				g_use_rtmscore_to_search_best = true;
			}
			i++;
		}else if (0 == strcmp(args[i], "-clig") && i < argc-1){
			if (0 == strcmp(args[i+1], "n")){
				g_fully_cared_ligatom_match = false; 
			}else{
				g_fully_cared_ligatom_match = true;
			}
			i++;
		}else if (0 == strcmp(args[i], "-d0") && i < argc-1){
			g_user_given_d0 = atof(args[i+1]);
			g_user_given_d02 = g_user_given_d0*g_user_given_d0;
			i++;
		}else if (0 == strcmp(args[i], "-da") && i < argc-1){
			if (0 == strcmp(args[i+1], "y")){
				g_use_chain_order_or_not = true; 
			}else{
				g_use_chain_order_or_not = false;
			}
			i++;
		}else if (0 == strcmp(args[i], "-ia") && i < argc-1){
			g_is_inputed_chain_alignment = true;
      		g_inputed_chain_alignment = args[i+1];
			i++;
		}else if (0 == strcmp(args[i], "-ri") && i < argc-1){
			if (0 == strcmp(args[i+1], "y")){
				g_use_res_orig_index_or_not = true;
			}else{
				g_use_res_orig_index_or_not = false;
			}
			i++;
		}else if (0 == strcmp(args[i], "-wt") && i < argc-1){
			if (0 == strcmp(args[i+1], "y")){
				g_is_load_H2O = true; 
			}else{
				g_is_load_H2O = false;
			}
			i++;
		}else if (0 == strcmp(args[i], "-atom-nuc") && i < argc-1){
			g_atomtype_nuc = args[i+1];
			i++;
		}else if (0 == strcmp(args[i], "-atom-res") && i < argc-1){
			g_atomtype_res = args[i+1];
			i++;
		}else if (0 == strcmp(args[i], "-odis") && i < argc-1){
			if (0 == strcmp(args[i+1], "y")){
				g_is_output_in_detail = true;
			}else {
				g_is_output_in_detail = false; 
			}
			i++;
		}else if (0 == strcmp(args[i], "-srm") && i < argc-1){
			isSaveUt = true;
      		saveUtPath = args[i+1];
      		i++;
		}else if (0 == strcmp(args[i], "-ssp") && i < argc-1){
			isSaveSuperpostion = true;
      		saveSuperpostionPath = args[i+1];
      		i++;
		}else if (0 == strcmp(args[i], "-o") && i < argc-1){
			isSaveRotTemplInPDB = true;
      		saveRotTemplInPDBPath = args[i+1];
      		i++;
      	}else if (0 == strcmp(args[i], "-oweb") && i < argc-1){
			isSaveWebNeedSuperposition = true;
			saveWebNeedSuperpositionPath = args[i+1];
      		i++;
		}else if (0 == strcmp(args[i], "-sid") && i < argc-1){
			g_seqid_cutoff = atof(args[i+1]);
			if (g_seqid_cutoff < 0.0) g_seqid_cutoff = 0.0;
			else if (g_seqid_cutoff > 1.0) g_seqid_cutoff = 1.0;
      		i++;
		}else if (0 == strcmp(args[i], "-nit") && i < argc-1){
			g_maxmum_number_of_iterations = atoi(args[i+1]);
			if (g_maxmum_number_of_iterations < 20) g_maxmum_number_of_iterations = 20;
			i++;			
		}else if (0 == strcmp(args[i], "-nLinit") && i < argc-1){
			g_maxmum_number_of_L_init = atoi(args[i+1]);
			if (g_maxmum_number_of_L_init < 6) g_maxmum_number_of_L_init = 6;
			i++;
		}else if (0 == strcmp(args[i], "-mode") && i < argc-1){
			if (0 == strcmp(args[i+1], "suf")){
				g_ali_type = SUPERFAST;
			}else if (0 == strcmp(args[i+1], "fast")){
				g_ali_type = FAST;
			}else if (0 == strcmp(args[i+1], "slw")){
				g_ali_type = SLOW;
			}else {
				if (0 == strcmp(args[i+1], "acc"))
					g_go_detail = true;
				g_ali_type = NORMAL;
			} 
			i++;
		}else if (0 == strcmp(args[i], "-h")){
			CBaseFunc::print_help(args[0]);
		}else if (0 == strcmp(args[i], "-v")){
			CBaseFunc::print_logo();
			exit(1);
		}else{
//			if (pdbFileInd == 0){
//				isQueryProvide = true;
//				query_pdb = args[i];
//			}else if (pdbFileInd == 1){
//				isTemplProvide = true;
//				templ_pdb = args[i];
//			}
			if (pdbFileInd == 0){
				// Structure 1
				isTemplProvide = true;
				templ_pdb = args[i];
			}else if (pdbFileInd == 1){				
			    // Structure 2
				isQueryProvide = true;
				query_pdb = args[i];
			}
			pdbFileInd++;
		}
	}
	
	if (!isQueryProvide){
		cout << "PLEASE PROVIDE PDB1.pdb FILE!!!" << endl << endl;
		CBaseFunc::print_help(args[0]);
	}else{
		fstream _file;
		_file.open(query_pdb.c_str(), ios::in);
		if(!_file){
			cout << query_pdb << " is not existed!" << endl;
			cout << "Please check it and input a correct PDB file path!" << endl;
			exit(1);
		}
	} 
	if (!isTemplProvide){
		cout << "PLEASE PROVIDE PDB2.pdb FILE!!!" << endl << endl;
		CBaseFunc::print_help(args[0]);
	}else{
		fstream _file;
		_file.open(templ_pdb.c_str(), ios::in);
		if(!_file){
			cout << templ_pdb << " is not existed!" << endl;
			cout << "Please check it and input a correct PDB file path!" << endl;
			exit(1);
		}
	} 
	
	if ((g_is_inputed_chain_alignment && g_use_chain_order_or_not) 
		 || (g_is_inputed_chain_alignment && g_ali_type == SLOW) 
		 || (g_use_chain_order_or_not && g_ali_type == SLOW)){
		cout << "WRONG USAGE: Note that, only one of the options of '-ia', '-da y', and '-pts slw' can be applied at the same time." << endl;
		cout << "But you apply two of options of '-ia', '-da y', and '-pts slw' at least. Please check and re-run it." << endl;
		exit(1);
	}
	
	const vector<CHAIN_PAIR>* chain_pair = NULL;
	if (g_is_inputed_chain_alignment)
		chain_pair = CBaseFunc::parse_matched_chain_pairs(g_inputed_chain_alignment);
	
	CTMscoreComplex* mtm = new CTMscoreComplex(query_pdb/*Structure 2*/, templ_pdb/*Structure 1*/, g_ali_type, chain_pair);
	mtm->print_result();
	
	if (isSaveRotTemplInPDB)
		mtm->save_roted_templ(saveRotTemplInPDBPath);
	
	if (isSaveUt)
		mtm->save_invUt(saveUtPath);
		
	if (isSaveSuperpostion)
		mtm->save_superposition_ditances(saveSuperpostionPath);
	
	if (isSaveWebNeedSuperposition)
		mtm->save_sup_pdb_for_web(saveWebNeedSuperpositionPath);
	
	if (NULL != chain_pair)
		delete chain_pair;
	delete mtm;
	
	return 0;
}

//int main(){
//	g_is_sup_protein = true;
//	g_is_sup_dna = false;
//	g_is_sup_rna = false;
//	g_is_sup_ligand = false;
//	
//	g_use_rtmscore_to_search_best = false;
//	CTMscoreComplex* mtm = new CTMscoreComplex("D:/install/Apache24/htdocs/codes/cpp/mTMscore/H1114/H1114.pdb", 
//			"D:/install/Apache24/htdocs/codes/cpp/mTMscore/H1114/H1114TS125_4.pdb", g_ali_type, NULL);
//	mtm->print_result();
//} 

inline void CTMscoreComplex::align_multimer_fast_buf_inaccuracy_using_nwalign_and_greadsearch(const bool& fast){
	int i, j, k, l, m, n, iL;
	char buf[2];
	string key;
	
	if (NULL == this->obj_level_ali){
		double gs_score, best_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(this->qsize, this->tsize);
		for (i = 0; i < this->qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < this->tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				for (k = 0; k < this->qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE& kmt = kmol->get_moltype();
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < this->tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[this->qsize];
				int* transpose_ali = new int[this->tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, this->qsize, this->tsize, ali, transpose_ali);
				if (gs_score > best_gs_score){
					bool is_ok = true;
					for (k = 0; k < this->qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						best_gs_score = gs_score;
						if (NULL != this->obj_level_ali)
							delete[] this->obj_level_ali;
						this->obj_level_ali = ali;
						ali = NULL;
					}		
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, this->qsize);
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query; 
	vector<ALIGN_PAIR> __a2b__aa__;
	vector<MOLTYPE> mts;
	string achain, bchain; 
	vector<string> aseq_vec, bseq_vec; 
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < this->qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}	
			}			
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]];
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		if (g_user_given_d0 <= 0)
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
		else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
		tmscore = tmscore * seqali_res_num / total_res_num;
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate rTMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
		this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
}


inline void CTMscoreComplex::align_multimer_fast_buf_inaccuracy_using_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast){
	int i, j, k, l, m, n, iL;
	char buf[2];
	string key;
	
	if (NULL == this->obj_level_ali){
		double gs_score, best_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(this->qsize, this->tsize);
		for (i = 0; i < this->qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < this->tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				for (k = 0; k < this->qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE& kmt = kmol->get_moltype();
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < this->tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[this->qsize];
				int* transpose_ali = new int[this->tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, this->qsize, this->tsize, ali, transpose_ali);
				if (gs_score > best_gs_score){
					bool is_ok = true;
					for (k = 0; k < this->qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							} 
						}
					}
					
					if (is_ok){
						best_gs_score = gs_score;
						if (NULL != this->obj_level_ali)
							delete[] this->obj_level_ali;
						this->obj_level_ali = ali;
						ali = NULL;
					}		
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
			}
		}
		CBaseFunc::delete2Darr(scomtx, this->qsize);
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query; 
	vector<ALIGN_PAIR> __a2b__aa__;
	vector<MOLTYPE> mts;
	string achain, bchain; 
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < this->qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);	
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}			
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		rtmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate TMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
			
		this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
}

inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_greadsearch(const bool& fast){
	int i, j, k, l, m, n, o, p, iL;
	char buf[2];
	string key;
	
	if (NULL == this->obj_level_ali){ 
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		for (i = 0; i < qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				for (k = 0; k < qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE& kmt = kmol->get_moltype();
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[qsize];
				int* transpose_ali = new int[tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
				gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				
				if (gs_score >= corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < qsize; k++){
							if (-1 == ali[k]) continue;
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const string& aseq = amol->get_seq_str();
							const string& bseq = bmol->get_seq_str();
							
							const MOLTYPE& amt = amol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, qsize); 
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]); 
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		if (g_user_given_d0 <= 0)
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
		else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
		tmscore = tmscore * seqali_res_num / total_res_num;
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate rTMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
		this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
}


inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand(const bool& fast){
	int i, j, k, l, m, n, o, p, iL;
	char buf[2];
	string key;
	
	int*** best_kl_lig_macther = NULL;
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		
		for (i = 0; i < qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				
				int*** kl_lig_macther = new int**[qsize];
				for (k = 0; k < qsize; k++)
					kl_lig_macther[k] = NULL;
				
				for (k = 0; k < qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE& kmt = kmol->get_moltype();
					if (LIGAND == kmt){
						kl_lig_macther[k] = new int*[tsize];
						for (l = 0; l < tsize; l++)
							kl_lig_macther[k][l] = NULL;
					}
					
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							if (kmt == LIGAND){
								k2l = new int[kxyzs.size()];
								LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
								kl_lig_macther[k][l] = k2l;
							}
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[qsize];
				int* transpose_ali = new int[tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
				gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				
//				cout << "DEBUG START" << endl;
//				for (int dei = 0; dei < qsize; dei++){
//					for (int dej = 0; dej < tsize; dej++){
//						cout << scomtx[dei][dej] << ' ';
//					}
//					cout << endl;
//				}
//				
//				cout << "========================= " << endl;
//				for (int dei = 0; dei < qsize; dei++)
//					cout << ali[dei] << ' ';
//				
//				cout << endl;
//				for (int dei = 0; dei < qsize; dei++)
//					cout << this->query->get_chain(dei) << ' ';
//				cout << endl;
//				for (int dei = 0; dei < qsize; dei++)
//					cout << this->templ->get_chain(ali[dei]) << ' ';
//				cout << endl;
//				
//				cout << "----------------- " << endl;
//				cout << endl << gs_score << endl;
//				
//				
//				
//				cout << "====== " << endl;
//				for (int dei = 0; dei < qsize; dei++)
//					cout << ali[dei] << ' ';
//				
//				cout << endl;
//				for (int dei = 0; dei < qsize; dei++)
//					cout << this->query->get_chain(dei) << ' ';
//				cout << endl;
//				for (int dei = 0; dei < qsize; dei++)
//					cout << this->templ->get_chain(ali[dei]) << ' ';
//				cout << endl;
//				
//				cout << "===== " << endl;
//				cout << endl << gs_score << endl;
//				
//				cout << "DEBUG END" << endl;
				
				if (gs_score >= corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < qsize; k++){
							if (-1 == ali[k]) continue;
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const string& aseq = amol->get_seq_str();
							const string& bseq = bmol->get_seq_str();
							
							const MOLTYPE& amt = amol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
							if (LIGAND == amt){
								a2b = kl_lig_macther[k][ali[k]];
							}
							
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
							
							if (NULL != best_kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != best_kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != best_kl_lig_macther[k][l])
												delete[] best_kl_lig_macther[k][l];
										}
										delete[] best_kl_lig_macther[k];
									}
								}
								delete[] best_kl_lig_macther;
							}
							best_kl_lig_macther = kl_lig_macther;
							kl_lig_macther = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
				
				// release kl_lig_macther
				if (NULL != kl_lig_macther){
					for (k = 0; k < qsize; k++){
						if (NULL != kl_lig_macther[k]){
							for (l = 0; l < tsize; l++){
								if (NULL != kl_lig_macther[k][l])
									delete[] kl_lig_macther[k][l];
							}
							delete[] kl_lig_macther[k];	
						}
					}
					delete[] kl_lig_macther;
				}
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, qsize); 
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]); 
			if (LIGAND == amt){
				if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
					a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
			}
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		if (g_user_given_d0 <= 0)
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
		else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
		tmscore = tmscore * seqali_res_num / total_res_num;
		
//		{
//			//=============================================
//			// [start] do chain mapping again
//			//=============================================
//			double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
//		
//			int*** kl_lig_macther = new int**[qsize];
//			for (k = 0; k < qsize; k++)
//				kl_lig_macther[k] = NULL;
//			
//			for (k = 0; k < qsize; k++){
//				Molecule* kmol = (*(this->query))[k];
//				const MOLTYPE& kmt = kmol->get_moltype();
//				if (LIGAND == kmt){
//					kl_lig_macther[k] = new int*[tsize];
//					for (l = 0; l < tsize; l++)
//						kl_lig_macther[k][l] = NULL;
//				}
//				
//				const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
//				vector<double*> roted_kxyzs;
//				for (l = 0; l < kxyzs.size(); l++)
//					roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
//				
//				for (l = 0; l < tsize; l++){
//					Molecule* lmol = (*(this->templ))[l];
//					const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
//					
//					int* k2l = this->get_ij_qt_match_mtx(k, l);
//					if (NULL == k2l){
//						scomtx[k][l] = 0.;
//					}else{
//						if (kmt == LIGAND){
//							k2l = new int[kxyzs.size()];
//							LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
//							kl_lig_macther[k][l] = k2l;
//						}
//						scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
//					}
//				}
//				
//				for (l = 0; l < roted_kxyzs.size(); l++)
//					delete[] roted_kxyzs[l];
//			}
//			
//			int* ali = new int[qsize];
//			int* transpose_ali = new int[tsize];
//			double gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
//			gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
//			gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
//			
//			cout << "DEBUG START" << endl;
//			for (int dei = 0; dei < qsize; dei++){
//				for (int dej = 0; dej < tsize; dej++){
//					cout << scomtx[dei][dej] << ' ';
//				}
//				cout << endl;
//			}
//			
//			cout << "========================= " << endl;
//			for (int dei = 0; dei < qsize; dei++)
//				cout << ali[dei] << ' ';
//			
//			cout << endl;
//			for (int dei = 0; dei < qsize; dei++)
//				cout << this->query->get_chain(dei) << ' ';
//			cout << endl;
//			for (int dei = 0; dei < qsize; dei++)
//				cout << this->templ->get_chain(ali[dei]) << ' ';
//			cout << endl;
//			
//			cout << "----------------- " << endl;
//			cout << endl << gs_score << endl;
//			cout << "DEBUG END" << endl; 
//			
//			bool is_ok = true;
//			for (k = 0; k < qsize; k++){
//				if (-1 != ali[k]){
//					if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
//						is_ok = false;
//						break;
//					}
//				}
//			}
//			
//			if (is_ok){
//				vector<double*> __rotted_aress;
//				vector<double*> __bress;
//				vector<MOLTYPE> __mts;
//				for (k = 0; k < qsize; k++){
//					if (-1 == ali[k]) continue;
//					Molecule* amol = (*(this->query))[k];
//					Molecule* bmol = (*(this->templ))[ali[k]];
//					
//					const string& aseq = amol->get_seq_str();
//					const string& bseq = bmol->get_seq_str();
//					
//					const MOLTYPE& amt = amol->get_moltype();
//					
//					const vector<double*> axyzs = amol->get_cared_xyz_vec();
//					const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
//					
//					int alen = axyzs.size();
//					int blen = bxyzs.size();
//					
//					int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
//					if (LIGAND == amt){
//						a2b = kl_lig_macther[k][ali[k]];
//					}
//					
//					for (l = 0; l < alen; l++){
//						if (-1 != a2b[l]){
//							__rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
//							__bress.push_back(bxyzs[a2b[l]]);
//							__mts.push_back(amt);
//						}
//					}
//				}
//				
//				double score = CBaseFunc::score_fun_once(__rotted_aress, __bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num); xxx notice g_user_given_d0
//				cout << "xxx : " <<  score << endl; 
//				if (score - tmscore > 1e-9){
//					tmscore = score;
//					if (NULL != this->obj_level_ali)
//						delete[] this->obj_level_ali;
//					this->obj_level_ali = ali;
//					ali = NULL;
//					
//					// reinitial information
//					aligned_chain_num = 0;
//					seqali_res_num = 0;
//					aress.clear();
//					bress.clear();
//					chain_index_corr_to_query.clear();
//					mts.clear();
//					__a2b__aa__.clear();
//					aseq_vec.clear();
//					bseq_vec.clear();
//					
//					for (i = 0; i < qsize; i++){
//						if (-1 == this->obj_level_ali[i]) continue;
//						
//						aligned_chain_num++;
//						
//						Molecule* amol = (*(this->query))[i];
//						Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
//						
//						const string& aseq = amol->get_seq_str();
//						const string& bseq = bmol->get_seq_str();
//						
//						const MOLTYPE& amt = amol->get_moltype();
//						
//						const vector<double*> axyzs = amol->get_cared_xyz_vec();
//						const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
//						
//						int alen = axyzs.size();
//						int blen = bxyzs.size();
//						
//						if (DETAIL == g_print_result_type){
//							achain = this->query->get_chain(i);
//							bchain = this->templ->get_chain(this->obj_level_ali[i]);
//							
//							if (LIGAND == amt){
//								aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
//								bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
//							}
//						}
//						
//						int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]); 
//						if (LIGAND == amt){
//							if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
//								a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
//						}
//						for (j = 0; j < alen; j++){
//							if (-1 != a2b[j]){
//								aress.push_back(axyzs[j]);
//								bress.push_back(bxyzs[a2b[j]]);
//								chain_index_corr_to_query.push_back(i);
//								mts.push_back(amt);
//								
//								if (DETAIL == g_print_result_type){
//									ALIGN_PAIR ap;
//									ap.qchain = achain;
//									ap.tchain = bchain;
//									ap.qind = j;
//									ap.tind = a2b[j];
//									ap.qoind = amol->get_ith_orig_index(j);
//									ap.toind = bmol->get_ith_orig_index(a2b[j]);
//									ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
//									ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
//									
//									if (LIGAND == amt){
//										ap.qaa = aseq_vec[j];
//										ap.taa = bseq_vec[a2b[j]]; 
//									}else{
//										ap.qaa = aseq[j];
//										ap.taa = bseq[a2b[j]]; 	
//									}
//									
//									__a2b__aa__.push_back(ap);
//								}
//								
//								seqali_res_num++;
//							}
//						}
//					}
//				}
//				
//				int nn = __rotted_aress.size();
//				for (k = 0; k < nn; k++)
//					delete[] __rotted_aress[k];
//			}
//			
//			if (NULL != ali)
//				delete[] ali;
//			delete[] transpose_ali;
//			
//			// release kl_lig_macther
//			if (NULL != kl_lig_macther){
//				for (k = 0; k < qsize; k++){
//					if (NULL != kl_lig_macther[k]){
//						for (l = 0; l < tsize; l++){
//							if (NULL != kl_lig_macther[k][l])
//								delete[] kl_lig_macther[k][l];
//						}
//						delete[] kl_lig_macther[k];	
//					}
//				}
//				delete[] kl_lig_macther;
//			}
//			
//			CBaseFunc::delete2Darr(scomtx, qsize); 
//			
//			//=============================================
//			// [end] do chain mapping again
//			//=============================================
//		}
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate rTMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
		this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
	
	// release best_kl_lig_macther
	if (NULL != best_kl_lig_macther){
		for (k = 0; k < qsize; k++){
			if (NULL != best_kl_lig_macther[k]){
				for (l = 0; l < tsize; l++){
					if (NULL != best_kl_lig_macther[k][l])
						delete[] best_kl_lig_macther[k][l];
				}
				delete[] best_kl_lig_macther[k];	
			}
		}
		delete[] best_kl_lig_macther;
	}
}


inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand_II(const bool& fast){
	int i, i2, j, j2, k, l, m, n, o, p, iL;
	char buf[2];
	string key;
	
	int*** best_kl_lig_macther = NULL;
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		
		for (i = 0; i < qsize; i++){
			Molecule* imol1 = (*(this->query))[i];
			const vector<double*> ixyzs1 = imol1->get_cared_xyz_vec();
			int ixyzs1_size = ixyzs1.size();
			for (i2 = i+1; i2 < qsize; i2++){
				Molecule* imol2 = (*(this->query))[i2];
				const vector<double*> ixyzs2 = imol2->get_cared_xyz_vec();
				int ixyzs2_size = ixyzs2.size();
				
				for (j = 0; j < tsize; j++){
					Molecule* jmol1 = (*(this->templ))[j];
					const vector<double*> jxyzs1 = jmol1->get_cared_xyz_vec();
					for (j2 = j+1; j2 < tsize; j2++){
						Molecule* jmol2 = (*(this->templ))[j2];
						const vector<double*> jxyzs2 = jmol2->get_cared_xyz_vec();
						
//						cout << "DEBUG " << i << ' ' << i2 << ' ' << j << ' ' << j2 << ' ' << best_sco << endl;
						{
							//--------------------------------------
							// [START] AB:AB
							//--------------------------------------
							vector<double*> ixyzs;
							vector<double*> jxyzs;
							
							int* i2j = this->get_ij_qt_match_mtx(i, j);
							if (NULL != i2j){
								for (int hi = 0; hi < ixyzs1_size; hi++){
									if (-1 != i2j[hi]){
										ixyzs.push_back(ixyzs1[hi]);
										jxyzs.push_back(jxyzs1[i2j[hi]]);
									}
								}
							}
							int* i22j2 = this->get_ij_qt_match_mtx(i2, j2);
							if (NULL != i22j2){
								for (int hi = 0; hi < ixyzs2_size; hi++){
									if (-1 != i22j2[hi]){
										ixyzs.push_back(ixyzs2[hi]);
										jxyzs.push_back(jxyzs2[i22j2[hi]]);
									}
								}
							}
							
							if (ixyzs.size() < 4)
								continue;
							
							// calculate the u and t
							CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, 8.0, true);
							
							int*** kl_lig_macther = new int**[qsize];
							for (k = 0; k < qsize; k++)
								kl_lig_macther[k] = NULL;
							
							for (k = 0; k < qsize; k++){
								Molecule* kmol = (*(this->query))[k];
								const MOLTYPE& kmt = kmol->get_moltype();
								if (LIGAND == kmt){
									kl_lig_macther[k] = new int*[tsize];
									for (l = 0; l < tsize; l++)
										kl_lig_macther[k][l] = NULL;
								}
								
								const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
								vector<double*> roted_kxyzs;
								for (l = 0; l < kxyzs.size(); l++)
									roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
								
								for (l = 0; l < tsize; l++){
									Molecule* lmol = (*(this->templ))[l];
									const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
									
									int* k2l = this->get_ij_qt_match_mtx(k, l);
									if (NULL == k2l){
										scomtx[k][l] = 0.;
									}else{
										if (kmt == LIGAND){
											k2l = new int[kxyzs.size()];
											LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
											kl_lig_macther[k][l] = k2l;
										}
										scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
									}
								}
								
								for (l = 0; l < roted_kxyzs.size(); l++)
									delete[] roted_kxyzs[l];
							}
							
							int* ali = new int[qsize];
							int* transpose_ali = new int[tsize];
							gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
							gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							if (gs_score >= corr_gs_score - 0.5){
								bool is_ok = true;
								for (k = 0; k < qsize; k++){
									if (-1 != ali[k]){
										if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
											is_ok = false;
											break;
										}
									}
								}
								
								if (is_ok){
									vector<double*> rotted_aress;
									vector<double*> bress;
									vector<MOLTYPE> mts;
									for (k = 0; k < qsize; k++){
										if (-1 == ali[k]) continue;
										Molecule* amol = (*(this->query))[k];
										Molecule* bmol = (*(this->templ))[ali[k]];
										
										const string& aseq = amol->get_seq_str();
										const string& bseq = bmol->get_seq_str();
										
										const MOLTYPE& amt = amol->get_moltype();
										
										const vector<double*> axyzs = amol->get_cared_xyz_vec();
										const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
										
										int alen = axyzs.size();
										int blen = bxyzs.size();
										
										int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
										if (LIGAND == amt){
											a2b = kl_lig_macther[k][ali[k]];
										}
										
										for (l = 0; l < alen; l++){
											if (-1 != a2b[l]){
												rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
												bress.push_back(bxyzs[a2b[l]]);
												mts.push_back(amt);
											}
										}
									}
									
									double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
									double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
									double score = score1>score2 ? score1 : score2;
									if (score > best_sco){
										best_sco = score;
										corr_gs_score = gs_score;
										if (NULL != this->obj_level_ali)
											delete[] this->obj_level_ali;
										this->obj_level_ali = ali;
										ali = NULL;
										
										if (NULL != best_kl_lig_macther){
											for (k = 0; k < qsize; k++){
												if (NULL != best_kl_lig_macther[k]){
													for (l = 0; l < tsize; l++){
														if (NULL != best_kl_lig_macther[k][l])
															delete[] best_kl_lig_macther[k][l];
													}
													delete[] best_kl_lig_macther[k];
												}
											}
											delete[] best_kl_lig_macther;
										}
										best_kl_lig_macther = kl_lig_macther;
										kl_lig_macther = NULL;
									}
									
									int nn = rotted_aress.size();
									for (k = 0; k < nn; k++)
										delete[] rotted_aress[k];
								}				
							}
							
							if (NULL != ali)
								delete[] ali;
							delete[] transpose_ali;
							
							// release kl_lig_macther
							if (NULL != kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != kl_lig_macther[k][l])
												delete[] kl_lig_macther[k][l];
										}
										delete[] kl_lig_macther[k];	
									}
								}
								delete[] kl_lig_macther;
							}
							
							//--------------------------------------
							// [END] AB:AB
							//--------------------------------------
						}
						{
							//--------------------------------------
							// [START] AB:BA
							//--------------------------------------
							vector<double*> ixyzs;
							vector<double*> jxyzs;
							
							int* i2j2 = this->get_ij_qt_match_mtx(i, j2);
							if (NULL != i2j2){
								for (int hi = 0; hi < ixyzs1_size; hi++){
									if (-1 != i2j2[hi]){
										ixyzs.push_back(ixyzs1[hi]);
										jxyzs.push_back(jxyzs2[i2j2[hi]]);
									}
								}
							}
							int* i22j = this->get_ij_qt_match_mtx(i2, j);
							if (NULL != i22j){
								for (int hi = 0; hi < ixyzs2_size; hi++){
									if (-1 != i22j[hi]){
										ixyzs.push_back(ixyzs2[hi]);
										jxyzs.push_back(jxyzs[i22j[hi]]);
									}
								}
							}
							
							if (ixyzs.size() < 4)
								continue;
							
							// calculate the u and t
							CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, 8.0, true);
							
							int*** kl_lig_macther = new int**[qsize];
							for (k = 0; k < qsize; k++)
								kl_lig_macther[k] = NULL;
							
							for (k = 0; k < qsize; k++){
								Molecule* kmol = (*(this->query))[k];
								const MOLTYPE& kmt = kmol->get_moltype();
								if (LIGAND == kmt){
									kl_lig_macther[k] = new int*[tsize];
									for (l = 0; l < tsize; l++)
										kl_lig_macther[k][l] = NULL;
								}
								
								const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
								vector<double*> roted_kxyzs;
								for (l = 0; l < kxyzs.size(); l++)
									roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
								
								for (l = 0; l < tsize; l++){
									Molecule* lmol = (*(this->templ))[l];
									const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
									
									int* k2l = this->get_ij_qt_match_mtx(k, l);
									if (NULL == k2l){
										scomtx[k][l] = 0.;
									}else{
										if (kmt == LIGAND){
											k2l = new int[kxyzs.size()];
											LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
											kl_lig_macther[k][l] = k2l;
										}
										scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
									}
								}
								
								for (l = 0; l < roted_kxyzs.size(); l++)
									delete[] roted_kxyzs[l];
							}
							
							int* ali = new int[qsize];
							int* transpose_ali = new int[tsize];
							gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
							gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							if (gs_score >= corr_gs_score - 0.5){
								bool is_ok = true;
								for (k = 0; k < qsize; k++){
									if (-1 != ali[k]){
										if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
											is_ok = false;
											break;
										}
									}
								}
								
								if (is_ok){
									vector<double*> rotted_aress;
									vector<double*> bress;
									vector<MOLTYPE> mts;
									for (k = 0; k < qsize; k++){
										if (-1 == ali[k]) continue;
										Molecule* amol = (*(this->query))[k];
										Molecule* bmol = (*(this->templ))[ali[k]];
										
										const string& aseq = amol->get_seq_str();
										const string& bseq = bmol->get_seq_str();
										
										const MOLTYPE& amt = amol->get_moltype();
										
										const vector<double*> axyzs = amol->get_cared_xyz_vec();
										const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
										
										int alen = axyzs.size();
										int blen = bxyzs.size();
										
										int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
										if (LIGAND == amt){
											a2b = kl_lig_macther[k][ali[k]];
										}
										
										for (l = 0; l < alen; l++){
											if (-1 != a2b[l]){
												rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
												bress.push_back(bxyzs[a2b[l]]);
												mts.push_back(amt);
											}
										}
									}
									
									double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
									double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
									double score = score1>score2 ? score1 : score2;
									if (score > best_sco){
										best_sco = score;
										corr_gs_score = gs_score;
										if (NULL != this->obj_level_ali)
											delete[] this->obj_level_ali;
										this->obj_level_ali = ali;
										ali = NULL;
										
										if (NULL != best_kl_lig_macther){
											for (k = 0; k < qsize; k++){
												if (NULL != best_kl_lig_macther[k]){
													for (l = 0; l < tsize; l++){
														if (NULL != best_kl_lig_macther[k][l])
															delete[] best_kl_lig_macther[k][l];
													}
													delete[] best_kl_lig_macther[k];
												}
											}
											delete[] best_kl_lig_macther;
										}
										best_kl_lig_macther = kl_lig_macther;
										kl_lig_macther = NULL;
									}
									
									int nn = rotted_aress.size();
									for (k = 0; k < nn; k++)
										delete[] rotted_aress[k];
								}				
							}
							
							if (NULL != ali)
								delete[] ali;
							delete[] transpose_ali;
							
							// release kl_lig_macther
							if (NULL != kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != kl_lig_macther[k][l])
												delete[] kl_lig_macther[k][l];
										}
										delete[] kl_lig_macther[k];	
									}
								}
								delete[] kl_lig_macther;
							}
							
							//--------------------------------------
							// [END] AB:BA
							//--------------------------------------
						}
					}
				}
		
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, qsize); 
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]); 
			if (LIGAND == amt){
				if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
					a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
			}
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		if (g_user_given_d0 <= 0)
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
		else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
		tmscore = tmscore * seqali_res_num / total_res_num;
		
		{
			//=============================================
			// [start] do chain mapping again
			//=============================================
			double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		
			int*** kl_lig_macther = new int**[qsize];
			for (k = 0; k < qsize; k++)
				kl_lig_macther[k] = NULL;
			
			for (k = 0; k < qsize; k++){
				Molecule* kmol = (*(this->query))[k];
				const MOLTYPE& kmt = kmol->get_moltype();
				if (LIGAND == kmt){
					kl_lig_macther[k] = new int*[tsize];
					for (l = 0; l < tsize; l++)
						kl_lig_macther[k][l] = NULL;
				}
				
				const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
				vector<double*> roted_kxyzs;
				for (l = 0; l < kxyzs.size(); l++)
					roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
				
				for (l = 0; l < tsize; l++){
					Molecule* lmol = (*(this->templ))[l];
					const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
					
					int* k2l = this->get_ij_qt_match_mtx(k, l);
					if (NULL == k2l){
						scomtx[k][l] = 0.;
					}else{
						if (kmt == LIGAND){
							k2l = new int[kxyzs.size()];
							LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
							kl_lig_macther[k][l] = k2l;
						}
						scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
					}
				}
				
				for (l = 0; l < roted_kxyzs.size(); l++)
					delete[] roted_kxyzs[l];
			}
			
			int* ali = new int[qsize];
			int* transpose_ali = new int[tsize];
			double gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
			gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			
//			cout << "DEBUG START" << endl;
//			for (int dei = 0; dei < qsize; dei++){
//				for (int dej = 0; dej < tsize; dej++){
//					cout << scomtx[dei][dej] << ' ';
//				}
//				cout << endl;
//			}
//			
//			cout << "========================= " << endl;
//			for (int dei = 0; dei < qsize; dei++)
//				cout << ali[dei] << ' ';
//			
//			cout << endl;
//			for (int dei = 0; dei < qsize; dei++)
//				cout << this->query->get_chain(dei) << ' ';
//			cout << endl;
//			for (int dei = 0; dei < qsize; dei++)
//				cout << this->templ->get_chain(ali[dei]) << ' ';
//			cout << endl;
//			
//			cout << "----------------- " << endl;
//			cout << endl << gs_score << endl;
//			cout << "DEBUG END" << endl; 
			
			bool is_ok = true;
			for (k = 0; k < qsize; k++){
				if (-1 != ali[k]){
					if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
						is_ok = false;
						break;
					}
				}
			}
			
			if (is_ok){
				vector<double*> __rotted_aress;
				vector<double*> __bress;
				vector<MOLTYPE> __mts;
				for (k = 0; k < qsize; k++){
					if (-1 == ali[k]) continue;
					Molecule* amol = (*(this->query))[k];
					Molecule* bmol = (*(this->templ))[ali[k]];
					
					const string& aseq = amol->get_seq_str();
					const string& bseq = bmol->get_seq_str();
					
					const MOLTYPE& amt = amol->get_moltype();
					
					const vector<double*> axyzs = amol->get_cared_xyz_vec();
					const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
					
					int alen = axyzs.size();
					int blen = bxyzs.size();
					
					int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
					if (LIGAND == amt){
						a2b = kl_lig_macther[k][ali[k]];
					}
					
					for (l = 0; l < alen; l++){
						if (-1 != a2b[l]){
							__rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
							__bress.push_back(bxyzs[a2b[l]]);
							__mts.push_back(amt);
						}
					}
				}
				
				double score = 0.; 
				if (g_user_given_d0 <= 0)
					score = CBaseFunc::score_fun_once(__rotted_aress, __bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
				else score = CBaseFunc::score_fun_once(__rotted_aress, __bress, mts, g_user_given_d02, g_user_given_d02, g_user_given_d02, g_user_given_d02, total_res_num);
				
				if (score - tmscore > 1e-9){
					tmscore = score;
					if (NULL != this->obj_level_ali)
						delete[] this->obj_level_ali;
					this->obj_level_ali = ali;
					ali = NULL;
					
					// reinitial information
					aligned_chain_num = 0;
					seqali_res_num = 0;
					aress.clear();
					bress.clear();
					chain_index_corr_to_query.clear();
					mts.clear();
					__a2b__aa__.clear();
					aseq_vec.clear();
					bseq_vec.clear();
					
					for (i = 0; i < qsize; i++){
						if (-1 == this->obj_level_ali[i]) continue;
						
						aligned_chain_num++;
						
						Molecule* amol = (*(this->query))[i];
						Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
						
						const string& aseq = amol->get_seq_str();
						const string& bseq = bmol->get_seq_str();
						
						const MOLTYPE& amt = amol->get_moltype();
						
						const vector<double*> axyzs = amol->get_cared_xyz_vec();
						const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
						
						int alen = axyzs.size();
						int blen = bxyzs.size();
						
						if (DETAIL == g_print_result_type){
							achain = this->query->get_chain(i);
							bchain = this->templ->get_chain(this->obj_level_ali[i]);
							
							if (LIGAND == amt){
								aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
								bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
							}
						}
						
						int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]); 
						if (LIGAND == amt){
							if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
								a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
						}
						for (j = 0; j < alen; j++){
							if (-1 != a2b[j]){
								aress.push_back(axyzs[j]);
								bress.push_back(bxyzs[a2b[j]]);
								chain_index_corr_to_query.push_back(i);
								mts.push_back(amt);
								
								if (DETAIL == g_print_result_type){
									ALIGN_PAIR ap;
									ap.qchain = achain;
									ap.tchain = bchain;
									ap.qind = j;
									ap.tind = a2b[j];
									ap.qoind = amol->get_ith_orig_index(j);
									ap.toind = bmol->get_ith_orig_index(a2b[j]);
									ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
									ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
									
									if (LIGAND == amt){
										ap.qaa = aseq_vec[j];
										ap.taa = bseq_vec[a2b[j]]; 
									}else{
										ap.qaa = aseq[j];
										ap.taa = bseq[a2b[j]]; 	
									}
									
									__a2b__aa__.push_back(ap);
								}
								
								seqali_res_num++;
							}
						}
					}
				}
				
				int nn = __rotted_aress.size();
				for (k = 0; k < nn; k++)
					delete[] __rotted_aress[k];
			}
			
			if (NULL != ali)
				delete[] ali;
			delete[] transpose_ali;
			
			// release kl_lig_macther
			if (NULL != kl_lig_macther){
				for (k = 0; k < qsize; k++){
					if (NULL != kl_lig_macther[k]){
						for (l = 0; l < tsize; l++){
							if (NULL != kl_lig_macther[k][l])
								delete[] kl_lig_macther[k][l];
						}
						delete[] kl_lig_macther[k];	
					}
				}
				delete[] kl_lig_macther;
			}
			
			CBaseFunc::delete2Darr(scomtx, qsize); 
			
			//=============================================
			// [end] do chain mapping again
			//=============================================
		}
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate rTMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
		this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
	
	// release best_kl_lig_macther
	if (NULL != best_kl_lig_macther){
		for (k = 0; k < qsize; k++){
			if (NULL != best_kl_lig_macther[k]){
				for (l = 0; l < tsize; l++){
					if (NULL != best_kl_lig_macther[k][l])
						delete[] best_kl_lig_macther[k][l];
				}
				delete[] best_kl_lig_macther[k];	
			}
		}
		delete[] best_kl_lig_macther;
	}
}

inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand_III(const bool& fast){
	int i, i2, j, j2, k, l, m, n, o, p, iL;
	char buf[2];
	string key;
	
	int*** best_kl_lig_macther = NULL;
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		double** all_single_scomtx = CBaseFunc::new2Darr(qsize, tsize);
		for (i = 0; i < qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				all_single_scomtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				individual_tmscore_mtx[i][j] = all_single_scomtx[i][j]; 
				
				int*** kl_lig_macther = new int**[qsize];
				for (k = 0; k < qsize; k++)
					kl_lig_macther[k] = NULL;
				
				for (k = 0; k < qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE& kmt = kmol->get_moltype();
					if (LIGAND == kmt){
						kl_lig_macther[k] = new int*[tsize];
						for (l = 0; l < tsize; l++)
							kl_lig_macther[k][l] = NULL;
					}
					
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							if (kmt == LIGAND){
								k2l = new int[kxyzs.size()];
								LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
								kl_lig_macther[k][l] = k2l;
							}
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[qsize];
				int* transpose_ali = new int[tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
				gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				
				if (gs_score >= corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < qsize; k++){
							if (-1 == ali[k]) continue;
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const string& aseq = amol->get_seq_str();
							const string& bseq = bmol->get_seq_str();
							
							const MOLTYPE& amt = amol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
							if (LIGAND == amt){
								a2b = kl_lig_macther[k][ali[k]];
							}
							
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
							
							if (NULL != best_kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != best_kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != best_kl_lig_macther[k][l])
												delete[] best_kl_lig_macther[k][l];
										}
										delete[] best_kl_lig_macther[k];
									}
								}
								delete[] best_kl_lig_macther;
							}
							best_kl_lig_macther = kl_lig_macther;
							kl_lig_macther = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
				
				// release kl_lig_macther
				if (NULL != kl_lig_macther){
					for (k = 0; k < qsize; k++){
						if (NULL != kl_lig_macther[k]){
							for (l = 0; l < tsize; l++){
								if (NULL != kl_lig_macther[k][l])
									delete[] kl_lig_macther[k][l];
							}
							delete[] kl_lig_macther[k];	
						}
					}
					delete[] kl_lig_macther;
				}
			}
		}
		
		if (g_go_detail){
			//=====================================================
			// [START] USING LARGE DISTANCE TO DO again 
			//=====================================================
			double larger_d0 = 30.;
			int* ali = new int[qsize];
			int* transpose_ali = new int[tsize];
			
			double __stm_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
			__stm_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, __stm_score, qsize);
			if (NULL != ali)
				delete[] ali;
			delete[] transpose_ali;
			
			if (__stm_score > 0.5 && best_sco < 0.5) {
				for (i = 0; i < qsize; i++){
					Molecule* imol = (*(this->query))[i];
					const vector<double*> ixyzs = imol->get_cared_xyz_vec();
					
					for (j = 0; j < tsize; j++){
						Molecule* jmol = (*(this->templ))[j];
						const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
						
						int* i2j = this->get_ij_qt_match_mtx(i, j);
						if (NULL == i2j)
							continue;
						
						// calculate the u and t
						CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, larger_d0, i2j, true);
						
						int*** kl_lig_macther = new int**[qsize];
						for (k = 0; k < qsize; k++)
							kl_lig_macther[k] = NULL;
						
						for (k = 0; k < qsize; k++){
							Molecule* kmol = (*(this->query))[k];
							const MOLTYPE& kmt = kmol->get_moltype();
							if (LIGAND == kmt){
								kl_lig_macther[k] = new int*[tsize];
								for (l = 0; l < tsize; l++)
									kl_lig_macther[k][l] = NULL;
							}
							
							const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
							vector<double*> roted_kxyzs;
							for (l = 0; l < kxyzs.size(); l++)
								roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
							
							for (l = 0; l < tsize; l++){
								Molecule* lmol = (*(this->templ))[l];
								const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
								
								int* k2l = this->get_ij_qt_match_mtx(k, l);
								if (NULL == k2l){
									scomtx[k][l] = 0.;
								}else{
									if (kmt == LIGAND){
										k2l = new int[kxyzs.size()];
										LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
										kl_lig_macther[k][l] = k2l;
									}
									scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
								}
							}
							
							for (l = 0; l < roted_kxyzs.size(); l++)
								delete[] roted_kxyzs[l];
						}
						
						int* ali = new int[qsize];
						int* transpose_ali = new int[tsize];
						gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
						gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
						gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
						
						if (gs_score >= corr_gs_score - 0.5){
							bool is_ok = true;
							for (k = 0; k < qsize; k++){
								if (-1 != ali[k]){
									if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
										is_ok = false;
										break;
									}
								}
							}
							
							if (is_ok){
								vector<double*> rotted_aress;
								vector<double*> bress;
								vector<MOLTYPE> mts;
								for (k = 0; k < qsize; k++){
									if (-1 == ali[k]) continue;
									Molecule* amol = (*(this->query))[k];
									Molecule* bmol = (*(this->templ))[ali[k]];
									
									const string& aseq = amol->get_seq_str();
									const string& bseq = bmol->get_seq_str();
									
									const MOLTYPE& amt = amol->get_moltype();
									
									const vector<double*> axyzs = amol->get_cared_xyz_vec();
									const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
									
									int alen = axyzs.size();
									int blen = bxyzs.size();
									
									int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
									if (LIGAND == amt){
										a2b = kl_lig_macther[k][ali[k]];
									}
									
									for (l = 0; l < alen; l++){
										if (-1 != a2b[l]){
											rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
											bress.push_back(bxyzs[a2b[l]]);
											mts.push_back(amt);
										}
									}
								}
								
								double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
								double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
								double score = score1>score2 ? score1 : score2;
								if (score > best_sco){
									best_sco = score;
									corr_gs_score = gs_score;
									if (NULL != this->obj_level_ali)
										delete[] this->obj_level_ali;
									this->obj_level_ali = ali;
									ali = NULL;
									
									if (NULL != best_kl_lig_macther){
										for (k = 0; k < qsize; k++){
											if (NULL != best_kl_lig_macther[k]){
												for (l = 0; l < tsize; l++){
													if (NULL != best_kl_lig_macther[k][l])
														delete[] best_kl_lig_macther[k][l];
												}
												delete[] best_kl_lig_macther[k];
											}
										}
										delete[] best_kl_lig_macther;
									}
									best_kl_lig_macther = kl_lig_macther;
									kl_lig_macther = NULL;
								}
								
								int nn = rotted_aress.size();
								for (k = 0; k < nn; k++)
									delete[] rotted_aress[k];
							}				
						}
						
						if (NULL != ali)
							delete[] ali;
						delete[] transpose_ali;
						
						// release kl_lig_macther
						if (NULL != kl_lig_macther){
							for (k = 0; k < qsize; k++){
								if (NULL != kl_lig_macther[k]){
									for (l = 0; l < tsize; l++){
										if (NULL != kl_lig_macther[k][l])
											delete[] kl_lig_macther[k][l];
									}
									delete[] kl_lig_macther[k];	
								}
							}
							delete[] kl_lig_macther;
						}
					}
				}
			}
			
			//=====================================================
			// [END] USING LARGE DISTANCE TO DO again 
			//=====================================================
		}
		
		//==========================================
		// [START] refinement chain mapping
		//==========================================
		if (NULL != this->obj_level_ali) {
			int i_num_which_ali_two_j_at_least = 0;
			for (i = 0; i < qsize; i++){
				int ali_j_num = 0;
				for (j = 0; j < tsize; j++){
					if (NULL != this->get_ij_qt_match_mtx(i, j)){
						ali_j_num++;
						if (ali_j_num >= 2){
							break;
						}
					}
				}
				if (ali_j_num >= 2)
					i_num_which_ali_two_j_at_least++; 
			}
			
			if (i_num_which_ali_two_j_at_least >= 2){
				int* refine_obj_level_ali = new int[qsize];
				for (i = 0; i < qsize; i++)
					refine_obj_level_ali[i] = this->obj_level_ali[i];
				
				for (i = 0; i < qsize; i++){
					j = refine_obj_level_ali[i];
					if (-1 == j) continue;
					
					Molecule* imol1 = (*(this->query))[i];
					const vector<double*> ixyzs1 = imol1->get_cared_xyz_vec();
					int ixyzs1_size = ixyzs1.size();
					
					Molecule* jmol1 = (*(this->templ))[j];
					const vector<double*> jxyzs1 = jmol1->get_cared_xyz_vec();
					
					for (i2 = i+1; i2 < qsize; i2++){
						j2 = refine_obj_level_ali[i2];
						if (-1 == j2) continue;
					
						Molecule* imol2 = (*(this->query))[i2];
						const vector<double*> ixyzs2 = imol2->get_cared_xyz_vec();
						int ixyzs2_size = ixyzs2.size();
						
						Molecule* jmol2 = (*(this->templ))[j2];
						const vector<double*> jxyzs2 = jmol2->get_cared_xyz_vec();
						
						{
							//--------------------------------------
							// [START] AB:AB
							//--------------------------------------
							vector<double*> ixyzs;
							vector<double*> jxyzs;
							
							int* i2j = this->get_ij_qt_match_mtx(i, j);
							if (NULL != i2j){
								for (int hi = 0; hi < ixyzs1_size; hi++){
									if (-1 != i2j[hi]){
										ixyzs.push_back(ixyzs1[hi]);
										jxyzs.push_back(jxyzs1[i2j[hi]]);
									}
								}
							}
							int* i22j2 = this->get_ij_qt_match_mtx(i2, j2);
							if (NULL != i22j2){
								for (int hi = 0; hi < ixyzs2_size; hi++){
									if (-1 != i22j2[hi]){
										ixyzs.push_back(ixyzs2[hi]);
										jxyzs.push_back(jxyzs2[i22j2[hi]]);
									}
								}
							}
							
							if (ixyzs.size() < 4)
								continue;
							
							// calculate the u and t
							CBaseFunc::cal_rot_tran_from_query_to_templ__II(ixyzs, jxyzs, u, t, 8.0, true);
							
							int*** kl_lig_macther = new int**[qsize];
							for (k = 0; k < qsize; k++)
								kl_lig_macther[k] = NULL;
							
							for (k = 0; k < qsize; k++){
								Molecule* kmol = (*(this->query))[k];
								const MOLTYPE& kmt = kmol->get_moltype();
								if (LIGAND == kmt){
									kl_lig_macther[k] = new int*[tsize];
									for (l = 0; l < tsize; l++)
										kl_lig_macther[k][l] = NULL;
								}
								
								const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
								vector<double*> roted_kxyzs;
								for (l = 0; l < kxyzs.size(); l++)
									roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
								
								for (l = 0; l < tsize; l++){
									Molecule* lmol = (*(this->templ))[l];
									const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
									
									int* k2l = this->get_ij_qt_match_mtx(k, l);
									if (NULL == k2l){
										scomtx[k][l] = 0.;
									}else{
										if (kmt == LIGAND){
											k2l = new int[kxyzs.size()];
											LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
											kl_lig_macther[k][l] = k2l;
										}
										scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
									}
								}
								
								for (l = 0; l < roted_kxyzs.size(); l++)
									delete[] roted_kxyzs[l];
							}
							
							int* ali = new int[qsize];
							int* transpose_ali = new int[tsize];
							gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
							gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							if (gs_score >= corr_gs_score - 0.5){
								bool is_ok = true;
								for (k = 0; k < qsize; k++){
									if (-1 != ali[k]){
										if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
											is_ok = false;
											break;
										}
									}
								}
								
								if (is_ok){
									vector<double*> rotted_aress;
									vector<double*> bress;
									vector<MOLTYPE> mts;
									for (k = 0; k < qsize; k++){
										if (-1 == ali[k]) continue;
										Molecule* amol = (*(this->query))[k];
										Molecule* bmol = (*(this->templ))[ali[k]];
										
										const string& aseq = amol->get_seq_str();
										const string& bseq = bmol->get_seq_str();
										
										const MOLTYPE& amt = amol->get_moltype();
										
										const vector<double*> axyzs = amol->get_cared_xyz_vec();
										const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
										
										int alen = axyzs.size();
										int blen = bxyzs.size();
										
										int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
										if (LIGAND == amt){
											a2b = kl_lig_macther[k][ali[k]];
										}
										
										for (l = 0; l < alen; l++){
											if (-1 != a2b[l]){
												rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
												bress.push_back(bxyzs[a2b[l]]);
												mts.push_back(amt);
											}
										}
									}
									
									double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
									double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
									double score = score1>score2 ? score1 : score2;
									if (score > best_sco){
										best_sco = score;
										corr_gs_score = gs_score;
										if (NULL != this->obj_level_ali)
											delete[] this->obj_level_ali;
										this->obj_level_ali = ali;
										ali = NULL;
										
										if (NULL != best_kl_lig_macther){
											for (k = 0; k < qsize; k++){
												if (NULL != best_kl_lig_macther[k]){
													for (l = 0; l < tsize; l++){
														if (NULL != best_kl_lig_macther[k][l])
															delete[] best_kl_lig_macther[k][l];
													}
													delete[] best_kl_lig_macther[k];
												}
											}
											delete[] best_kl_lig_macther;
										}
										best_kl_lig_macther = kl_lig_macther;
										kl_lig_macther = NULL;
									}
									
									int nn = rotted_aress.size();
									for (k = 0; k < nn; k++)
										delete[] rotted_aress[k];
								}				
							}
							
							if (NULL != ali)
								delete[] ali;
							delete[] transpose_ali;
							
							// release kl_lig_macther
							if (NULL != kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != kl_lig_macther[k][l])
												delete[] kl_lig_macther[k][l];
										}
										delete[] kl_lig_macther[k];	
									}
								}
								delete[] kl_lig_macther;
							}
							
							//--------------------------------------
							// [END] AB:AB
							//--------------------------------------
						}
					}
				}
				delete[] refine_obj_level_ali;	
			}
		}
		//==========================================
		// [END] refinement chain mapping
		//==========================================
		
		CBaseFunc::delete2Darr(all_single_scomtx, qsize);
		CBaseFunc::delete2Darr(scomtx, qsize); 
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]); 
			if (LIGAND == amt){
				if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
					a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
			}
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		if (g_user_given_d0 <= 0)
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
		else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
		tmscore = tmscore * seqali_res_num / total_res_num;
		
		{
			//=============================================
			// [start] do chain mapping again
			//=============================================
			double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		
			int*** kl_lig_macther = new int**[qsize];
			for (k = 0; k < qsize; k++)
				kl_lig_macther[k] = NULL;
			
			for (k = 0; k < qsize; k++){
				Molecule* kmol = (*(this->query))[k];
				const MOLTYPE& kmt = kmol->get_moltype();
				if (LIGAND == kmt){
					kl_lig_macther[k] = new int*[tsize];
					for (l = 0; l < tsize; l++)
						kl_lig_macther[k][l] = NULL;
				}
				
				const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
				vector<double*> roted_kxyzs;
				for (l = 0; l < kxyzs.size(); l++)
					roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
				
				for (l = 0; l < tsize; l++){
					Molecule* lmol = (*(this->templ))[l];
					const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
					
					int* k2l = this->get_ij_qt_match_mtx(k, l);
					if (NULL == k2l){
						scomtx[k][l] = 0.;
					}else{
						if (kmt == LIGAND){
							k2l = new int[kxyzs.size()];
							LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
							kl_lig_macther[k][l] = k2l;
						}
						scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
					}
				}
				
				for (l = 0; l < roted_kxyzs.size(); l++)
					delete[] roted_kxyzs[l];
			}
			
			int* ali = new int[qsize];
			int* transpose_ali = new int[tsize];
			double gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
			gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			
			bool is_ok = true;
			for (k = 0; k < qsize; k++){
				if (-1 != ali[k]){
					if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
						is_ok = false;
						break;
					}
				}
			}
			
			if (is_ok){
				vector<double*> __rotted_aress;
				vector<double*> __bress;
				vector<MOLTYPE> __mts;
				for (k = 0; k < qsize; k++){
					if (-1 == ali[k]) continue;
					Molecule* amol = (*(this->query))[k];
					Molecule* bmol = (*(this->templ))[ali[k]];
					
					const string& aseq = amol->get_seq_str();
					const string& bseq = bmol->get_seq_str();
					
					const MOLTYPE& amt = amol->get_moltype();
					
					const vector<double*> axyzs = amol->get_cared_xyz_vec();
					const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
					
					int alen = axyzs.size();
					int blen = bxyzs.size();
					
					int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
					if (LIGAND == amt){
						a2b = kl_lig_macther[k][ali[k]];
					}
					
					for (l = 0; l < alen; l++){
						if (-1 != a2b[l]){
							__rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
							__bress.push_back(bxyzs[a2b[l]]);
							__mts.push_back(amt);
						}
					}
				}
				
				double score = 0.; 
				if (g_user_given_d0 <= 0)
					score = CBaseFunc::score_fun_once(__rotted_aress, __bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
				else score = CBaseFunc::score_fun_once(__rotted_aress, __bress, mts, g_user_given_d02, g_user_given_d02, g_user_given_d02, g_user_given_d02, total_res_num);
				
				if (score - tmscore > 1e-9){
					tmscore = score;
					if (NULL != this->obj_level_ali)
						delete[] this->obj_level_ali;
					this->obj_level_ali = ali;
					ali = NULL;
					
					// reinitial information
					aligned_chain_num = 0;
					seqali_res_num = 0;
					aress.clear();
					bress.clear();
					chain_index_corr_to_query.clear();
					mts.clear();
					__a2b__aa__.clear();
					aseq_vec.clear();
					bseq_vec.clear();
					
					for (i = 0; i < qsize; i++){
						if (-1 == this->obj_level_ali[i]) continue;
						
						aligned_chain_num++;
						
						Molecule* amol = (*(this->query))[i];
						Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
						
						const string& aseq = amol->get_seq_str();
						const string& bseq = bmol->get_seq_str();
						
						const MOLTYPE& amt = amol->get_moltype();
						
						const vector<double*> axyzs = amol->get_cared_xyz_vec();
						const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
						
						int alen = axyzs.size();
						int blen = bxyzs.size();
						
						if (DETAIL == g_print_result_type){
							achain = this->query->get_chain(i);
							bchain = this->templ->get_chain(this->obj_level_ali[i]);
							
							if (LIGAND == amt){
								aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
								bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
							}
						}
						
						int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]); 
						if (LIGAND == amt){
							if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
								a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
						}
						for (j = 0; j < alen; j++){
							if (-1 != a2b[j]){
								aress.push_back(axyzs[j]);
								bress.push_back(bxyzs[a2b[j]]);
								chain_index_corr_to_query.push_back(i);
								mts.push_back(amt);
								
								if (DETAIL == g_print_result_type){
									ALIGN_PAIR ap;
									ap.qchain = achain;
									ap.tchain = bchain;
									ap.qind = j;
									ap.tind = a2b[j];
									ap.qoind = amol->get_ith_orig_index(j);
									ap.toind = bmol->get_ith_orig_index(a2b[j]);
									ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
									ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
									
									if (LIGAND == amt){
										ap.qaa = aseq_vec[j];
										ap.taa = bseq_vec[a2b[j]]; 
									}else{
										ap.qaa = aseq[j];
										ap.taa = bseq[a2b[j]]; 	
									}
									
									__a2b__aa__.push_back(ap);
								}
								
								seqali_res_num++;
							}
						}
					}
				}
				
				int nn = __rotted_aress.size();
				for (k = 0; k < nn; k++)
					delete[] __rotted_aress[k];
			}
			
			if (NULL != ali)
				delete[] ali;
			delete[] transpose_ali;
			
			// release kl_lig_macther
			if (NULL != kl_lig_macther){
				for (k = 0; k < qsize; k++){
					if (NULL != kl_lig_macther[k]){
						for (l = 0; l < tsize; l++){
							if (NULL != kl_lig_macther[k][l])
								delete[] kl_lig_macther[k][l];
						}
						delete[] kl_lig_macther[k];	
					}
				}
				delete[] kl_lig_macther;
			}
			
			CBaseFunc::delete2Darr(scomtx, qsize); 
			
			//=============================================
			// [end] do chain mapping again
			//=============================================
		}
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate rTMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
		this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
	
	// release best_kl_lig_macther
	if (NULL != best_kl_lig_macther){
		for (k = 0; k < qsize; k++){
			if (NULL != best_kl_lig_macther[k]){
				for (l = 0; l < tsize; l++){
					if (NULL != best_kl_lig_macther[k][l])
						delete[] best_kl_lig_macther[k][l];
				}
				delete[] best_kl_lig_macther[k];	
			}
		}
		delete[] best_kl_lig_macther;
	}
}


inline void CTMscoreComplex::prepare_for_monomer(){
	int i;
	double d0_tmp;
	for (i = 0; i < this->qsize; i++){
		Molecule* ai = (*(this->query))[i];
		const MOLTYPE& amt = ai->get_moltype();
		
		if (PROTEIN == amt){
			d0_tmp = CBaseFunc::d0_of_tmscore(ai->size());
		} else if (DNA == amt) {
			d0_tmp = CBaseFunc::d0_of_tmscore_c3prime(ai->size());
		} else if (RNA == amt) {
			d0_tmp = CBaseFunc::d0_of_tmscore_c3prime(ai->size());
		}else {
			d0_tmp = CBaseFunc::d0_of_lsscore(ai->size());
		}
		
		chain_index_corr_to_query__aa_num[i] = ai->size();
		chain_index_corr_to_query__moltype[i] = amt;
		chain_index_corr_to_query__d0[i] = d0_tmp;
		chain_index_corr_to_query__d02[i] = d0_tmp * d0_tmp;
	}
}

inline void CTMscoreComplex::prepare_for_multimer(const bool& use_atom_or_residue_index_order_or_not, const bool& use_chain_order_or_not){
	int i, j, iL;
	char buf[2];
	string key;
	double d0_tmp;
	
	total_res_in_pro_num = 0;
	total_nuc_in_dna_num = 0;
	total_nuc_in_rna_num = 0;
	total_atm_in_lig_num = 0;
	
	for (i = 0; i < this->qsize; i++){
		Molecule* ai = (*(this->query))[i];
		const MOLTYPE& amt = ai->get_moltype();
		
		if (PROTEIN == amt){
			total_res_in_pro_num += ai->size();
			d0_tmp = CBaseFunc::d0_of_tmscore(ai->size());
		} else if (DNA == amt) {
			total_nuc_in_dna_num += ai->size();
			d0_tmp = CBaseFunc::d0_of_tmscore_c3prime(ai->size());
		} else if (RNA == amt) {
			total_nuc_in_rna_num += ai->size();
			d0_tmp = CBaseFunc::d0_of_tmscore_c3prime(ai->size());
		}else {
			total_atm_in_lig_num += ai->size();
			d0_tmp = CBaseFunc::d0_of_lsscore(ai->size());
		}
		
		chain_index_corr_to_query__aa_num[i] = ai->size();
		chain_index_corr_to_query__moltype[i] = amt;
		chain_index_corr_to_query__d0[i] = d0_tmp;
		chain_index_corr_to_query__d02[i] = d0_tmp * d0_tmp;
	}
	
	for (i = 0; i < tsize; i++){
		Molecule* bi = (*(this->templ))[i];
		chain_index_corr_to_templ__aa_num[i] = bi->size();
	} 
	
	total_res_num = total_res_in_pro_num + total_nuc_in_dna_num + total_nuc_in_rna_num + total_atm_in_lig_num;
//	d0_pro = CBaseFunc::d0_of_tmscore(total_res_in_pro_num);
//	d0_dna = CBaseFunc::d0_of_tmscore_c3prime(total_nuc_in_dna_num);
//	d0_rna = CBaseFunc::d0_of_tmscore_c3prime(total_nuc_in_rna_num);
//	d0_lig = CBaseFunc::d0_of_lsscore(total_atm_in_lig_num);
	d0_pro = CBaseFunc::d0_of_tmscore(total_res_num);
	d0_dna = CBaseFunc::d0_of_tmscore_c3prime(total_res_num);
	d0_rna = CBaseFunc::d0_of_tmscore_c3prime(total_res_num);
	d0_lig = CBaseFunc::d0_of_lsscore(total_res_num);
	d02_pro = d0_pro*d0_pro;
	d02_dna = d0_dna*d0_dna;
	d02_rna = d0_rna*d0_rna;
	d02_lig = d0_lig*d0_lig;
	
	if (use_atom_or_residue_index_order_or_not){
		generate_atom_alignment_of_each_ligand_pair_using_index_order(use_chain_order_or_not);
		generate_residue_alignment_of_each_molecule_pair_resdiue_index(use_chain_order_or_not);
	}else{
		generate_atom_alignment_of_each_ligand_pair_using_greedysearch(use_chain_order_or_not);
		generate_residue_alignment_of_each_molecule_pair_using_nwalign(use_chain_order_or_not);
	}
	
}

inline void CTMscoreComplex::generate_residue_alignment_of_each_molecule_pair_resdiue_index(const bool& use_chain_order_or_not){
	if (use_chain_order_or_not){
		int i, j, ii, k, iali, jj, l; 
		char c1, c2;
		for (i = 0; i < this->qsize; i++){
			if (i >= tsize) break;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[i];
			
			const MOLTYPE& amt = amol->get_moltype();
			const MOLTYPE& bmt = bmol->get_moltype();
			if (LIGAND == amt) continue;
			
			if (NULL == qt_match_mtx[i]){
				qt_match_mtx[i] = new int*[tsize];
				for (j = 0; j < tsize; j++)
					qt_match_mtx[i][j] = NULL;
			}
			
			if (amt == bmt){  // the first charactor is the flag in {PROTEIN, DNA, RNA, LIGAND}
				int aL = amol->size();
				int bL = bmol->size();
				
				qt_match_mtx[i][i] = new int[aL];
				
				int identity_num = 0;
				bool* is_used = CBaseFunc::new1Dbool(bL);
				for (ii = 0; ii < aL; ii++){
					k = amol->get_ith_orig_index(ii);
					c1 = amol->get_ith_char_following_orig_index_vec(ii);
					iali = -1;
					for (jj = 0; jj < bL; jj++){
						if (is_used[jj]) continue;
						l = bmol->get_ith_orig_index(jj);
						c2 = bmol->get_ith_char_following_orig_index_vec(jj); 
						if (k == l && c1 == c2){
							iali = jj;
							identity_num++;
							break;
						}
					}
					qt_match_mtx[i][i][ii] = iali;
				}
				
				if (identity_num < 4){
					cout << "WRONG USAGE: the number of matched residue indexes between the two "<< (i+1) <<"-th molecules in the two input pdbs is less than 4." << endl;
					cout << "|-- It is hard to superimpose them. IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
					cout << "|---->"<<TOOL_EXE<<" pdb1.pdb pdb2.pdb" << endl;
					exit(1);
				}
				
				delete[] is_used;
			}else{
				cout << "WRONG USAGE: the molecule type is not matched between the two "<< (i+1) <<"-th molecules in the two input pdbs, respectively." << endl;
				cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
				cout << "|---->"<<TOOL_EXE<<" pdb1.pdb pdb2.pdb" << endl;
				exit(1);
			}
		}
	}else{ // use_chain_order_or_not == false
		int i, j, k, l, aL, bL, m, n;
		char c1, c2;
		for (i = 0; i < this->qsize; i++){
			aL = chain_index_corr_to_query__aa_num[i];
			Molecule* amol = (*(this->query))[i];
			const MOLTYPE& amt = amol->get_moltype(); 
			const string& aseq = amol->get_seq_str();
			
			if (LIGAND == amt) continue;
			
			if (NULL == qt_match_mtx[i]){
				qt_match_mtx[i] = new int*[tsize];
				for (j = 0; j < tsize; j++)
					qt_match_mtx[i][j] = NULL;
			}
			
			for (j = 0; j < this->tsize; j++){
				bL = chain_index_corr_to_templ__aa_num[j];
				Molecule* bmol = (*(this->templ))[j];
				const MOLTYPE& bmt = bmol->get_moltype(); 
				const string& bseq = bmol->get_seq_str();
				
				if (amt == bmt){  // the first charactor is the flag in {PROTEIN, DNA, RNA, LIGAND}
					int identity_num = 0;
					bool* is_used = CBaseFunc::new1Dbool(bL);
					qt_match_mtx[i][j] = new int[aL];
					for (m = 0; m < aL; m++){
						k = amol->get_ith_orig_index(m);
						c1 = amol->get_ith_char_following_orig_index_vec(m);
						int iali = -1;
						for (n = 0; n < bL; n++){
							l = bmol->get_ith_orig_index(n);
							c2 = bmol->get_ith_char_following_orig_index_vec(n);
							if (k == l && c1 == c2){
								iali = n;
								break;
							}
						}
						qt_match_mtx[i][j][m] = iali;
						if (-1 != iali) identity_num++;
					}
					
					if (identity_num < 4){
						delete[] qt_match_mtx[i][j];
						qt_match_mtx[i][j] = NULL;
					}
				}
			}
		}
	}
}

inline void CTMscoreComplex::generate_residue_alignment_of_each_molecule_pair_using_nwalign(const bool& use_chain_order_or_not){
	if (use_chain_order_or_not){
		int i, j, k;
		char buf[2];
		
		map<string, int*> is_similar_q2t_for_macromol;   // key = amt+aseq+'$'+bseq
		map<string, vector<char> > corr_q_vec_seq;
		map<string, vector<char> > corr_t_vec_seq;
		map<string, vector<int> > corr_q_vec_oind;
		map<string, vector<int> > corr_t_vec_oind;
		 
		map<string, int*>::iterator it;
		vector<string*> keys;
		
		for (i = 0; i < qsize; i++){
			if (i >= tsize) break;
			
			Molecule* ai = (*(this->query))[i];
			Molecule* bi = (*(this->templ))[i];
			const string& aseq = ai->get_seq_str();
			const string& bseq = bi->get_seq_str();
			const MOLTYPE& amt = ai->get_moltype();
			const MOLTYPE& bmt = bi->get_moltype();
			
			if (LIGAND == amt) continue;
			if (amt != bmt){
				cout << "ERROR: The types of the " << (i+1) << "-th molecules in the 1st and 2nd complexes are different." << endl;
				cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
				cout << "|----->"<< TOOL_EXE << " pdb1.pdb pdb2.pdb" << endl;
				exit(1);
			}
			
			const vector<char>& aseq_vec = ai->get_seq_vec();
			const vector<char>& bseq_vec = bi->get_seq_vec();
			const vector<int>& aoind = ai->get_orig_index_vec();
			const vector<int>& boind = bi->get_orig_index_vec();
			
			sprintf(buf, "%d", amt);
			string* key = new string(buf + aseq + "$" + bseq);
			if (is_similar_q2t_for_macromol.end() == is_similar_q2t_for_macromol.find(*key)){
				is_similar_q2t_for_macromol[*key] = NULL;
				corr_q_vec_seq[*key] = aseq_vec;
				corr_t_vec_seq[*key] = bseq_vec;
				corr_q_vec_oind[*key] = aoind;
				corr_t_vec_oind[*key] = boind;
			}
				
			keys.push_back(key);
		}
		
		for (it = is_similar_q2t_for_macromol.begin(); it != is_similar_q2t_for_macromol.end(); it++){
			string mt_aseq_bseq = it->first;
			char mt = mt_aseq_bseq[0];
			string aseq_bseq = mt_aseq_bseq.substr(1);
			vector<string> lc = CBaseFunc::stringSplit(aseq_bseq, '$');
			string& aseq = lc[0];
			string& bseq = lc[1];
			
			const vector<char>& aseq_vec = corr_q_vec_seq[mt_aseq_bseq];
			const vector<char>& bseq_vec = corr_t_vec_seq[mt_aseq_bseq];
			const vector<int>& aoind = corr_q_vec_oind[mt_aseq_bseq];
			const vector<int>& boind = corr_t_vec_oind[mt_aseq_bseq];
			
			int iL = aseq.size();
			int* i2j = NULL;
			if (aseq == bseq){
				i2j = new int[iL];
				for (i = 0; i < iL; i++)
					i2j[i] = i;
			}else{
				int* __i2j = new int[iL];
				
				int identical_ali_num = 0;
				double seqid = 0.;
				if (mt == '0' /*PROTEIN*/)
					seqid = CNWalign::nwalign(aseq_vec, aoind, bseq_vec, boind, PROTEIN, __i2j, identical_ali_num);
				else if (mt == '1' /*DNA*/)
					seqid = CNWalign::nwalign(aseq_vec, aoind, bseq_vec, boind, DNA, __i2j, identical_ali_num);
				else if (mt == '2' /*RNA*/)
					seqid = CNWalign::nwalign(aseq_vec, aoind, bseq_vec, boind, RNA, __i2j, identical_ali_num);
				
				if (seqid > g_seqid_cutoff && identical_ali_num > 4){
					i2j = new int[iL];
					for (i = 0; i < iL; i++)
						i2j[i] = __i2j[i];
				}
				
				delete[] __i2j;				
			}
			if (NULL != i2j)
				is_similar_q2t_for_macromol[mt_aseq_bseq] = i2j;
			else {
				cout << "ERROR: There is one pair molecule sequence is not homologous at least under the current sequence identity cutoff (" << g_seqid_cutoff << ")." << endl;
				cout << "|-- You can use the option of \"-sid\" to reset a smaller sequence identity cutoff and re-run " << TOOL_EXE << "." << endl;
				cout << "|-- Or, if you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
				cout << "|----->"<< TOOL_EXE << " pdb1.pdb pdb2.pdb" << endl;
				exit(1);
			}
		}
		
		for (i = 0; i < qsize; i++){
			if (i >= tsize) break;
			
			string& key = *keys[i];
			if (NULL == qt_match_mtx[i]){
				qt_match_mtx[i] = new int*[tsize];
				for (j = 0; j < tsize; j++)
					qt_match_mtx[i][j] = NULL;
			}
			
			int iL = this->chain_index_corr_to_query__aa_num[i];
			int* i2j = is_similar_q2t_for_macromol[key];
			if (NULL != i2j){
				qt_match_mtx[i][i] = new int[iL];
				for (k = 0; k < iL; k++)
					qt_match_mtx[i][i][k] = i2j[k];
			}
		}
		
		for (i = 0; i < keys.size(); i++){
			if (NULL != keys[i])
				delete keys[i];
		}
	}else{ // use_chain_order_or_not == false
		int i, j, k, iL;
		char buf[20];
	
		map<string, char> as;  // mt + aseq
		map<string, vector<char> > corr_a_vec_seq;
		map<string, vector<int> > corr_a_vec_oind;
		
		map<string, char> bs;  // mt + bseq
		map<string, vector<char> > corr_b_vec_seq;
		map<string, vector<int> > corr_b_vec_oind;
		
		map<string, char>::iterator ait;
		map<string, char>::iterator bit;
		map<string, int*> is_similar_q2t_for_macromol;   // key = amt+aseq+'$'+bmt+bseq
		
		vector<string*> akeys;
		vector<string*> bkeys;
		
		for (i = 0; i < this->qsize; i++){
			Molecule* ai = (*(this->query))[i];
			const string& aseq = ai->get_seq_str();
			const MOLTYPE& amt = ai->get_moltype();
			
			if (LIGAND != amt){
				sprintf(buf, "%d", amt);
				string* key = new string(buf + aseq);
				if (as.end() == as.find(*key)){
					as[*key] = 'X';
					
					const vector<char>& vec_seq = ai->get_seq_vec();
					const vector<int>& oind = ai->get_orig_index_vec();
					corr_a_vec_seq[*key] = vec_seq; 
					corr_a_vec_oind[*key] = oind;
				}
				akeys.push_back(key);
			}else akeys.push_back(NULL);
		}
		
		for (i = 0; i < this->tsize; i++){
			Molecule* bi = (*(this->templ))[i];
			const string& bseq = bi->get_seq_str();
			const MOLTYPE& bmt = bi->get_moltype();
			
			if (LIGAND != bmt){
				sprintf(buf, "%d", bmt);
				string* key = new string(buf + bseq);
				if (bs.end() == bs.find(*key)){
					bs[*key] = 'X';	
						
					const vector<char>& vec_seq = bi->get_seq_vec();
					const vector<int>& oind = bi->get_orig_index_vec();
					corr_b_vec_seq[*key] = vec_seq; 
					corr_b_vec_oind[*key] = oind;
				}
				bkeys.push_back(key);
			}else bkeys.push_back(NULL);
		}
		
		for (ait = as.begin(); ait != as.end(); ait++){
			string i_mt_seq = ait->first;
			char imt = i_mt_seq[0];
			iL = i_mt_seq.size() - 1;
			string iseq = i_mt_seq.substr(1);
			
			const vector<char>& iseq_vec = corr_a_vec_seq[i_mt_seq];
			const vector<int>& ioind = corr_a_vec_oind[i_mt_seq];
			
			for (bit = bs.begin(); bit != bs.end(); bit++){
				string j_mt_seq = bit->first;
				char jmt = j_mt_seq[0];
				string jseq = j_mt_seq.substr(1);
				
				const vector<char>& jseq_vec = corr_b_vec_seq[j_mt_seq];
				const vector<int>& joind = corr_b_vec_oind[j_mt_seq];
				
				int* i2j = NULL;
				if (imt == jmt){  // the first charactor is the flag in {PROTEIN, DNA, RNA, LIGAND}
					if (i_mt_seq == j_mt_seq){
						i2j = new int[iL];
						for (i = 0; i < iL; i++)
							i2j[i] = i;
					}else{
						int* __i2j = new int[iL];
				
						int identical_ali_num = 0;
						double seqid = 0.;
						if (imt == '0' /*PROTEIN*/)
							seqid = CNWalign::nwalign(iseq_vec, ioind, jseq_vec, joind, PROTEIN, __i2j, identical_ali_num);
						else if (imt == '1' /*DNA*/)
							seqid = CNWalign::nwalign(iseq_vec, ioind, jseq_vec, joind, DNA, __i2j, identical_ali_num);
						else if (imt == '2' /*RNA*/)
							seqid = CNWalign::nwalign(iseq_vec, ioind, jseq_vec, joind, RNA, __i2j, identical_ali_num);
						
						if (seqid > g_seqid_cutoff && identical_ali_num > 4){
							i2j = new int[iL];
							for (i = 0; i < iL; i++)
								i2j[i] = __i2j[i];
						}
						
						delete[] __i2j;			
					}
				}
				
				if (NULL != i2j) is_similar_q2t_for_macromol[i_mt_seq+"$"+j_mt_seq] = i2j;
			}
		}
		
		if (is_similar_q2t_for_macromol.size() <= 0){
			cout << "ERROR: There is no molecule sequence pair are homologous under the current sequence identity cutoff (" << g_seqid_cutoff << ")." << endl;
			cout << "|-- You can use the option of \"-sid\" to reset a smaller sequence identity cutoff and re-run " << TOOL_EXE << "." << endl;
			cout << "|-- Or, if you make sure the residue/nucletide index in two inputs are matched, you could try option of \"-ri y\"" << endl;
			cout << "|----->"<< TOOL_EXE << " -ri y pdb1.pdb pdb2.pdb" << endl;
			exit(1);
		} 
		
		for (i = 0; i < qsize; i++){
			if (NULL == akeys[i]) continue;
			
			string& ikey = *akeys[i];
			if (NULL == qt_match_mtx[i]){
				qt_match_mtx[i] = new int*[tsize];
				for (j = 0; j < tsize; j++)
					qt_match_mtx[i][j] = NULL;
			}
			
			iL = this->chain_index_corr_to_query__aa_num[i];
			for (j = 0; j < tsize; j++){
				if (NULL == bkeys[j]) continue;
				
				string& jkey = *bkeys[j];
				int* i2j = is_similar_q2t_for_macromol[ikey+"$"+jkey];
				if (NULL != i2j){
					qt_match_mtx[i][j] = new int[iL];
					for (k = 0; k < iL; k++)
						qt_match_mtx[i][j][k] = i2j[k];
				}
			}
		}
		
		for (i = 0; i < akeys.size(); i++){
			if (NULL != akeys[i])
				delete akeys[i];
		}
		
		for (i = 0; i < bkeys.size(); i++){
			if (NULL != bkeys[i])
				delete bkeys[i];
		}
			
		for (map<string, int*>::iterator sit = is_similar_q2t_for_macromol.begin(); sit != is_similar_q2t_for_macromol.end(); sit++){
			if (NULL != sit->second)
				delete[] sit->second;
		}
		map<string, int*>().swap(is_similar_q2t_for_macromol);	
	}
}

inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast){
	int i, j, k, l, iL;
	
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		for (i = 0; i < qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				for (k = 0; k < qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE kmt = kmol->get_moltype();
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[qsize];
				int* transpose_ali = new int[tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
				gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				
				if (gs_score >= corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < qsize; k++){
							if (-1 == ali[k]) continue;
							
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const MOLTYPE& amt = amol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, qsize);
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		rtmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate TMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
			
		this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
}


inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore_fullycared_ligand(const bool& fast){
	int i, j, k, l, iL;
	
	int*** best_kl_lig_macther = NULL;
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		for (i = 0; i < qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				
				int*** kl_lig_macther = new int**[qsize];
				for (k = 0; k < qsize; k++)
					kl_lig_macther[k] = NULL;
				
				for (k = 0; k < qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE kmt = kmol->get_moltype();
					if (LIGAND == kmt){
						kl_lig_macther[k] = new int*[tsize];
						for (l = 0; l < tsize; l++)
							kl_lig_macther[k][l] = NULL;
					}
					
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							if (LIGAND == kmt){
								k2l = new int[kxyzs.size()];
								LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
								kl_lig_macther[k][l] = k2l;
							}
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[qsize];
				int* transpose_ali = new int[tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
				gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				
				if (gs_score >= corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < qsize; k++){
							if (-1 == ali[k]) continue;
							
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const MOLTYPE& amt = amol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
							if (LIGAND == amt){
								a2b = kl_lig_macther[k][ali[k]];
							}
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
							
							if (NULL != best_kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != best_kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != best_kl_lig_macther[k][l])
												delete[] best_kl_lig_macther[k][l];
										}
										delete[] best_kl_lig_macther[k];
									}
								}
								delete[] best_kl_lig_macther;
							}
							best_kl_lig_macther = kl_lig_macther;
							kl_lig_macther = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
				
				// release kl_lig_macther
				if (NULL != kl_lig_macther){
					for (k = 0; k < qsize; k++){
						if (NULL != kl_lig_macther[k]){
							for (l = 0; l < tsize; l++){
								if (NULL != kl_lig_macther[k][l])
									delete[] kl_lig_macther[k][l];
							}
							delete[] kl_lig_macther[k];	
						}
					}
					delete[] kl_lig_macther;
				}
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, qsize);
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			if (LIGAND == amt){
				if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
					a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
			}
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		rtmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate TMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
			
		this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
	
	// release best_kl_lig_macther
	if (NULL != best_kl_lig_macther){
		for (k = 0; k < qsize; k++){
			if (NULL != best_kl_lig_macther[k]){
				for (l = 0; l < tsize; l++){
					if (NULL != best_kl_lig_macther[k][l])
						delete[] best_kl_lig_macther[k][l];
				}
				delete[] best_kl_lig_macther[k];	
			}
		}
		delete[] best_kl_lig_macther;
	}
}


inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore_fullycared_ligand_III(const bool& fast){
	int i, i2, j, j2, k, l, iL;
	
	int*** best_kl_lig_macther = NULL;
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
		double** all_single_scomtx = CBaseFunc::new2Darr(qsize, tsize);
		for (i = 0; i < qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				all_single_scomtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				individual_tmscore_mtx[i][j] = all_single_scomtx[i][j];
				
				int*** kl_lig_macther = new int**[qsize];
				for (k = 0; k < qsize; k++)
					kl_lig_macther[k] = NULL;
				
				for (k = 0; k < qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE kmt = kmol->get_moltype();
					if (LIGAND == kmt){
						kl_lig_macther[k] = new int*[tsize];
						for (l = 0; l < tsize; l++)
							kl_lig_macther[k][l] = NULL;
					}
					
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							if (LIGAND == kmt){
								k2l = new int[kxyzs.size()];
								LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
								kl_lig_macther[k][l] = k2l;
							}
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[qsize];
				int* transpose_ali = new int[tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
				gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
				
				if (gs_score >= corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < qsize; k++){
							if (-1 == ali[k]) continue;
							
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const MOLTYPE& amt = amol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
							if (LIGAND == amt){
								a2b = kl_lig_macther[k][ali[k]];
							}
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
							
							if (NULL != best_kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != best_kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != best_kl_lig_macther[k][l])
												delete[] best_kl_lig_macther[k][l];
										}
										delete[] best_kl_lig_macther[k];
									}
								}
								delete[] best_kl_lig_macther;
							}
							best_kl_lig_macther = kl_lig_macther;
							kl_lig_macther = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
				
				// release kl_lig_macther
				if (NULL != kl_lig_macther){
					for (k = 0; k < qsize; k++){
						if (NULL != kl_lig_macther[k]){
							for (l = 0; l < tsize; l++){
								if (NULL != kl_lig_macther[k][l])
									delete[] kl_lig_macther[k][l];
							}
							delete[] kl_lig_macther[k];	
						}
					}
					delete[] kl_lig_macther;
				}
			}
		}
		
		if (g_go_detail){
			//=====================================================
			// [START] USING LARGE DISTANCE TO DO again 
			//=====================================================
			double larger_d0 = 30.;
			int* ali = new int[qsize];
			int* transpose_ali = new int[tsize];
			
			double __stm_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
			__stm_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, __stm_score, qsize);
			if (NULL != ali)
				delete[] ali;
			delete[] transpose_ali;
			
			if (__stm_score > 0.5 && best_sco < 0.5) {
				for (i = 0; i < qsize; i++){
					Molecule* imol = (*(this->query))[i];
					const vector<double*> ixyzs = imol->get_cared_xyz_vec();
					
					for (j = 0; j < tsize; j++){
						Molecule* jmol = (*(this->templ))[j];
						const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
						
						int* i2j = this->get_ij_qt_match_mtx(i, j);
						if (NULL == i2j)
							continue;
						
						// calculate the u and t
						CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, larger_d0, i2j, true);
						
						int*** kl_lig_macther = new int**[qsize];
						for (k = 0; k < qsize; k++)
							kl_lig_macther[k] = NULL;
						
						for (k = 0; k < qsize; k++){
							Molecule* kmol = (*(this->query))[k];
							const MOLTYPE& kmt = kmol->get_moltype();
							if (LIGAND == kmt){
								kl_lig_macther[k] = new int*[tsize];
								for (l = 0; l < tsize; l++)
									kl_lig_macther[k][l] = NULL;
							}
							
							const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
							vector<double*> roted_kxyzs;
							for (l = 0; l < kxyzs.size(); l++)
								roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
							
							for (l = 0; l < tsize; l++){
								Molecule* lmol = (*(this->templ))[l];
								const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
								
								int* k2l = this->get_ij_qt_match_mtx(k, l);
								if (NULL == k2l){
									scomtx[k][l] = 0.;
								}else{
									if (kmt == LIGAND){
										k2l = new int[kxyzs.size()];
										LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
										kl_lig_macther[k][l] = k2l;
									}
									scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
								}
							}
							
							for (l = 0; l < roted_kxyzs.size(); l++)
								delete[] roted_kxyzs[l];
						}
						
						int* ali = new int[qsize];
						int* transpose_ali = new int[tsize];
						gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
						gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
						gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
						
						if (gs_score >= corr_gs_score - 0.5){
							bool is_ok = true;
							for (k = 0; k < qsize; k++){
								if (-1 != ali[k]){
									if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
										is_ok = false;
										break;
									}
								}
							}
							
							if (is_ok){
								vector<double*> rotted_aress;
								vector<double*> bress;
								vector<MOLTYPE> mts;
								for (k = 0; k < qsize; k++){
									if (-1 == ali[k]) continue;
									Molecule* amol = (*(this->query))[k];
									Molecule* bmol = (*(this->templ))[ali[k]];
									
									const string& aseq = amol->get_seq_str();
									const string& bseq = bmol->get_seq_str();
									
									const MOLTYPE& amt = amol->get_moltype();
									
									const vector<double*> axyzs = amol->get_cared_xyz_vec();
									const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
									
									int alen = axyzs.size();
									int blen = bxyzs.size();
									
									int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
									if (LIGAND == amt){
										a2b = kl_lig_macther[k][ali[k]];
									}
									
									for (l = 0; l < alen; l++){
										if (-1 != a2b[l]){
											rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
											bress.push_back(bxyzs[a2b[l]]);
											mts.push_back(amt);
										}
									}
								}
								
								double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
								double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
								double score = score1>score2 ? score1 : score2;
								if (score > best_sco){
									best_sco = score;
									corr_gs_score = gs_score;
									if (NULL != this->obj_level_ali)
										delete[] this->obj_level_ali;
									this->obj_level_ali = ali;
									ali = NULL;
									
									if (NULL != best_kl_lig_macther){
										for (k = 0; k < qsize; k++){
											if (NULL != best_kl_lig_macther[k]){
												for (l = 0; l < tsize; l++){
													if (NULL != best_kl_lig_macther[k][l])
														delete[] best_kl_lig_macther[k][l];
												}
												delete[] best_kl_lig_macther[k];
											}
										}
										delete[] best_kl_lig_macther;
									}
									best_kl_lig_macther = kl_lig_macther;
									kl_lig_macther = NULL;
								}
								
								int nn = rotted_aress.size();
								for (k = 0; k < nn; k++)
									delete[] rotted_aress[k];
							}				
						}
						
						if (NULL != ali)
							delete[] ali;
						delete[] transpose_ali;
						
						// release kl_lig_macther
						if (NULL != kl_lig_macther){
							for (k = 0; k < qsize; k++){
								if (NULL != kl_lig_macther[k]){
									for (l = 0; l < tsize; l++){
										if (NULL != kl_lig_macther[k][l])
											delete[] kl_lig_macther[k][l];
									}
									delete[] kl_lig_macther[k];	
								}
							}
							delete[] kl_lig_macther;
						}
					}
				}
			}
			
			//=====================================================
			// [END] USING LARGE DISTANCE TO DO again 
			//=====================================================
		}
		
		//==========================================
		// [START] refinement chain mapping
		//==========================================
		if (NULL != this->obj_level_ali) {
			int i_num_which_ali_two_j_at_least = 0;
			for (i = 0; i < qsize; i++){
				int ali_j_num = 0;
				for (j = 0; j < tsize; j++){
					if (NULL != this->get_ij_qt_match_mtx(i, j)){
						ali_j_num++;
						if (ali_j_num >= 2){
							break;
						}
					}
				}
				if (ali_j_num >= 2)
					i_num_which_ali_two_j_at_least++; 
			}
			
			if (i_num_which_ali_two_j_at_least >= 2){
				int* refine_obj_level_ali = new int[qsize];
				for (i = 0; i < qsize; i++)
					refine_obj_level_ali[i] = this->obj_level_ali[i];
				
				for (i = 0; i < qsize; i++){
					j = refine_obj_level_ali[i];
					if (-1 == j) continue;
					
					Molecule* imol1 = (*(this->query))[i];
					const vector<double*> ixyzs1 = imol1->get_cared_xyz_vec();
					int ixyzs1_size = ixyzs1.size();
					
					Molecule* jmol1 = (*(this->templ))[j];
					const vector<double*> jxyzs1 = jmol1->get_cared_xyz_vec();
					
					for (i2 = i+1; i2 < qsize; i2++){
						j2 = refine_obj_level_ali[i2];
						if (-1 == j2) continue;
					
						Molecule* imol2 = (*(this->query))[i2];
						const vector<double*> ixyzs2 = imol2->get_cared_xyz_vec();
						int ixyzs2_size = ixyzs2.size();
						
						Molecule* jmol2 = (*(this->templ))[j2];
						const vector<double*> jxyzs2 = jmol2->get_cared_xyz_vec();
						
						{
							//--------------------------------------
							// [START] AB:AB
							//--------------------------------------
							vector<double*> ixyzs;
							vector<double*> jxyzs;
							
							int* i2j = this->get_ij_qt_match_mtx(i, j);
							if (NULL != i2j){
								for (int hi = 0; hi < ixyzs1_size; hi++){
									if (-1 != i2j[hi]){
										ixyzs.push_back(ixyzs1[hi]);
										jxyzs.push_back(jxyzs1[i2j[hi]]);
									}
								}
							}
							int* i22j2 = this->get_ij_qt_match_mtx(i2, j2);
							if (NULL != i22j2){
								for (int hi = 0; hi < ixyzs2_size; hi++){
									if (-1 != i22j2[hi]){
										ixyzs.push_back(ixyzs2[hi]);
										jxyzs.push_back(jxyzs2[i22j2[hi]]);
									}
								}
							}
							
							if (ixyzs.size() < 4)
								continue;
							
							// calculate the u and t
							CBaseFunc::cal_rot_tran_from_query_to_templ__II(ixyzs, jxyzs, u, t, 8.0, true);
							
							int*** kl_lig_macther = new int**[qsize];
							for (k = 0; k < qsize; k++)
								kl_lig_macther[k] = NULL;
							
							for (k = 0; k < qsize; k++){
								Molecule* kmol = (*(this->query))[k];
								const MOLTYPE& kmt = kmol->get_moltype();
								if (LIGAND == kmt){
									kl_lig_macther[k] = new int*[tsize];
									for (l = 0; l < tsize; l++)
										kl_lig_macther[k][l] = NULL;
								}
								
								const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
								vector<double*> roted_kxyzs;
								for (l = 0; l < kxyzs.size(); l++)
									roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
								
								for (l = 0; l < tsize; l++){
									Molecule* lmol = (*(this->templ))[l];
									const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
									
									int* k2l = this->get_ij_qt_match_mtx(k, l);
									if (NULL == k2l){
										scomtx[k][l] = 0.;
									}else{
										if (kmt == LIGAND){
											k2l = new int[kxyzs.size()];
											LigAtomMatcher::quick_identity_atom_align(*(q_ligAtomMatch_obj_vec[k]), *(t_ligAtomMatch_obj_vec[l]), u, t, k2l, 4.0);
											kl_lig_macther[k][l] = k2l;
										}
										scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l);
									}
								}
								
								for (l = 0; l < roted_kxyzs.size(); l++)
									delete[] roted_kxyzs[l];
							}
							
							int* ali = new int[qsize];
							int* transpose_ali = new int[tsize];
							gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
							gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
							if (gs_score >= corr_gs_score - 0.5){
								bool is_ok = true;
								for (k = 0; k < qsize; k++){
									if (-1 != ali[k]){
										if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
											is_ok = false;
											break;
										}
									}
								}
								
								if (is_ok){
									vector<double*> rotted_aress;
									vector<double*> bress;
									vector<MOLTYPE> mts;
									for (k = 0; k < qsize; k++){
										if (-1 == ali[k]) continue;
										Molecule* amol = (*(this->query))[k];
										Molecule* bmol = (*(this->templ))[ali[k]];
										
										const string& aseq = amol->get_seq_str();
										const string& bseq = bmol->get_seq_str();
										
										const MOLTYPE& amt = amol->get_moltype();
										
										const vector<double*> axyzs = amol->get_cared_xyz_vec();
										const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
										
										int alen = axyzs.size();
										int blen = bxyzs.size();
										
										int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
										if (LIGAND == amt){
											a2b = kl_lig_macther[k][ali[k]];
										}
										
										for (l = 0; l < alen; l++){
											if (-1 != a2b[l]){
												rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
												bress.push_back(bxyzs[a2b[l]]);
												mts.push_back(amt);
											}
										}
									}
									
									double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
									double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
									double score = score1>score2 ? score1 : score2;
									if (score > best_sco){
										best_sco = score;
										corr_gs_score = gs_score;
										if (NULL != this->obj_level_ali)
											delete[] this->obj_level_ali;
										this->obj_level_ali = ali;
										ali = NULL;
										
										if (NULL != best_kl_lig_macther){
											for (k = 0; k < qsize; k++){
												if (NULL != best_kl_lig_macther[k]){
													for (l = 0; l < tsize; l++){
														if (NULL != best_kl_lig_macther[k][l])
															delete[] best_kl_lig_macther[k][l];
													}
													delete[] best_kl_lig_macther[k];
												}
											}
											delete[] best_kl_lig_macther;
										}
										best_kl_lig_macther = kl_lig_macther;
										kl_lig_macther = NULL;
									}
									
									int nn = rotted_aress.size();
									for (k = 0; k < nn; k++)
										delete[] rotted_aress[k];
								}				
							}
							
							if (NULL != ali)
								delete[] ali;
							delete[] transpose_ali;
							
							// release kl_lig_macther
							if (NULL != kl_lig_macther){
								for (k = 0; k < qsize; k++){
									if (NULL != kl_lig_macther[k]){
										for (l = 0; l < tsize; l++){
											if (NULL != kl_lig_macther[k][l])
												delete[] kl_lig_macther[k][l];
										}
										delete[] kl_lig_macther[k];	
									}
								}
								delete[] kl_lig_macther;
							}
							
							//--------------------------------------
							// [END] AB:AB
							//--------------------------------------
						}
					}
				}
				delete[] refine_obj_level_ali;	
			}
		}
		//==========================================
		// [END] refinement chain mapping
		//==========================================
		
		CBaseFunc::delete2Darr(all_single_scomtx, qsize);
		CBaseFunc::delete2Darr(scomtx, qsize);
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			if (LIGAND == amt){
				if (NULL != best_kl_lig_macther && NULL != best_kl_lig_macther[i])
					a2b = best_kl_lig_macther[i][this->obj_level_ali[i]];
			}
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		rtmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate TMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
			
		this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
	
	// release best_kl_lig_macther
	if (NULL != best_kl_lig_macther){
		for (k = 0; k < qsize; k++){
			if (NULL != best_kl_lig_macther[k]){
				for (l = 0; l < tsize; l++){
					if (NULL != best_kl_lig_macther[k][l])
						delete[] best_kl_lig_macther[k][l];
				}
				delete[] best_kl_lig_macther[k];	
			}
		}
		delete[] best_kl_lig_macther;
	}
}


inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_not_greadsearch(const bool& fast){
	if (this->qsize != this->tsize){
		cout << "WRONG USAGE: the molecule number between two input pdbs are different: " << this->qsize << ", " << this->tsize << endl;
		cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
		cout << "|----->"<<TOOL_EXE<<" pdb1.pdb pdb2.pdb" << endl;
	}
	
	int i, j, k, l, iL;
	char buf[2];
	string key, akey, bkey;
	
	
	this->obj_level_ali = new int[this->qsize];
	for (i = 0; i < this->qsize; i++){
		if (i < tsize)
			obj_level_ali[i] = i;
		else obj_level_ali[i] = -1;
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	
	int aligned_chain_num = 0;
	int seqali_res_num = 0;
	for (i = 0; i < this->qsize; i++){
		if (-1 == obj_level_ali[i]) continue;
		aligned_chain_num++;
		
		Molecule* amol = (*(this->query))[i];
		Molecule* bmol = (*(this->templ))[obj_level_ali[i]];
		
		const string& aseq = amol->get_seq_str();
		const string& bseq = bmol->get_seq_str();
		
		const MOLTYPE& amt = amol->get_moltype();
		
		const vector<double*> axyzs = amol->get_cared_xyz_vec();
		const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
		
		int alen = axyzs.size();
		int blen = bxyzs.size();
		
		if (DETAIL == g_print_result_type){
			achain = this->query->get_chain(i);
			bchain = this->templ->get_chain(this->obj_level_ali[i]);
			
			if (LIGAND == amt){
				aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
				bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
			}
		}
		
		int* a2b = this->get_ij_qt_match_mtx(i, obj_level_ali[i]);
		individual_tmscore_mtx[i][this->obj_level_ali[i]] = CBaseFunc::cal_rot_tran_from_query_to_templ__(axyzs, bxyzs, u, t, this->chain_index_corr_to_query__d0[i], a2b, fast);
		for (j = 0; j < alen; j++){
			if (-1 != a2b[j]){
				aress.push_back(axyzs[j]);
				bress.push_back(bxyzs[a2b[j]]);
				chain_index_corr_to_query.push_back(i);
				mts.push_back(amt);
				
				if (DETAIL == g_print_result_type){
					ALIGN_PAIR ap;
					ap.qchain = achain;
					ap.tchain = bchain;
					ap.qind = j;
					ap.tind = a2b[j];
					ap.qoind = amol->get_ith_orig_index(j);
					ap.toind = bmol->get_ith_orig_index(a2b[j]);
					ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
					ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
					
					if (LIGAND == amt){
						ap.qaa = aseq_vec[j];
						ap.taa = bseq_vec[a2b[j]]; 
					}else{
						ap.qaa = aseq[j];
						ap.taa = bseq[a2b[j]]; 	
					}
					
					__a2b__aa__.push_back(ap);
				}
				
				seqali_res_num++;
			}
		}
	}
	
	if (g_user_given_d0 <= 0)
		tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
	else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
	tmscore = tmscore * seqali_res_num / total_res_num;
	
	if (DETAIL == g_print_result_type){
		int n = aress.size();
		for (i = 0; i < n; i++){
			__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
			this->aa_level_ali.push_back(__a2b__aa__[i]);
		}
	}
	
	// calculate rTMscore
	vector<double*> rotted_axyzs;
	for (l = 0; l < aress.size(); l++)
		rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
	
	this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
	
	for (l = 0; l < rotted_axyzs.size(); l++)
		delete[] rotted_axyzs[l];
	vector<double*>().swap(rotted_axyzs);
}

inline void CTMscoreComplex::align_multimer_normal_using_nwalign_and_not_greadsearch_with_max_rtmscore(const bool& fast){
	int i, j, k, l, iL;
	char buf[2];
	string key, akey, bkey;
	
	this->obj_level_ali = new int[this->qsize];
	for (i = 0; i < this->qsize; i++){
		if (i < tsize) 
			obj_level_ali[i] = i;
		else obj_level_ali[i] = -1;
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	int aligned_chain_num = 0;
	int seqali_res_num = 0;
	for (i = 0; i < this->qsize; i++){
		if (-1 != this->obj_level_ali[i]){
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			individual_tmscore_mtx[i][this->obj_level_ali[i]] = CBaseFunc::cal_rot_tran_from_query_to_templ__(axyzs, bxyzs, u, t, this->chain_index_corr_to_query__d0[i], a2b, fast);
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
	}
	
	rtmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
		
	if (DETAIL == g_print_result_type){
		int n = aress.size();
		for (i = 0; i < n; i++){
			__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
			this->aa_level_ali.push_back(__a2b__aa__[i]);
		}
	}
	
	// calculate TMscore
	vector<double*> rotted_axyzs;
	for (l = 0; l < aress.size(); l++)
		rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
	this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
	
	for (l = 0; l < rotted_axyzs.size(); l++)
		delete[] rotted_axyzs[l];
	vector<double*>().swap(rotted_axyzs);
}


inline void CTMscoreComplex::align_multimer_normal_using_not_nwalign_and_not_greadsearch(const bool& fast){
	int i, ii, j, jj, k, l, iali;
	char buf[2];
	string key, akey, bkey;
	
	this->obj_level_ali = new int[this->qsize];
	for (i = 0; i < this->qsize; i++){
		if (i < tsize) 
			obj_level_ali[i] = i;
		else obj_level_ali[i] = -1;
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	int aligned_chain_num = 0;
	int seqali_res_num = 0;
	for (i = 0; i < this->qsize; i++){
		if (-1 == this->obj_level_ali[i]) 
			continue;
		
		aligned_chain_num++;
		
		Molecule* amol = (*(this->query))[i];
		Molecule* bmol = (*(this->templ))[i];
		
		const string& aseq = amol->get_seq_str();
		const string& bseq = bmol->get_seq_str();
		
		const MOLTYPE& amt = amol->get_moltype(); 
		const MOLTYPE& bmt = bmol->get_moltype();
		
		const vector<double*> axyzs = amol->get_cared_xyz_vec();
		const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
		
		int alen = axyzs.size();
		int blen = bxyzs.size();
		
		if (DETAIL == g_print_result_type){
			achain = this->query->get_chain(i);
			bchain = this->templ->get_chain(this->obj_level_ali[i]);
			
			if (LIGAND == amt){
				aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
				bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
			}
		}
		
		int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
		individual_tmscore_mtx[i][this->obj_level_ali[i]] = CBaseFunc::cal_rot_tran_from_query_to_templ__(axyzs, bxyzs, u, t, this->chain_index_corr_to_query__d0[i], a2b, fast);
		for (j = 0; j < alen; j++){
			if (-1 != a2b[j]){
				aress.push_back(axyzs[j]);
				bress.push_back(bxyzs[a2b[j]]);
				chain_index_corr_to_query.push_back(i);
				mts.push_back(amt);
				
				if (DETAIL == g_print_result_type){
					ALIGN_PAIR ap;
					ap.qchain = achain;
					ap.tchain = bchain;
					ap.qind = j;
					ap.tind = a2b[j];
					ap.qoind = amol->get_ith_orig_index(j);
					ap.toind = bmol->get_ith_orig_index(a2b[j]);
					ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
					ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
					
					if (LIGAND == amt){
						ap.qaa = aseq_vec[j];
						ap.taa = bseq_vec[a2b[j]]; 
					}else{
						ap.qaa = aseq[j];
						ap.taa = bseq[a2b[j]]; 	
					}
					
					__a2b__aa__.push_back(ap);
				}
				
				seqali_res_num++;
			}
		}
	}
	
	if (g_user_given_d0 <= 0)
		tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
	else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
	tmscore = tmscore * seqali_res_num / total_res_num;
	
	if (DETAIL == g_print_result_type){
		int n = aress.size();
		for (i = 0; i < n; i++){
			__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
			this->aa_level_ali.push_back(__a2b__aa__[i]);
		}
	}
	
	// calculate rTMscore
	vector<double*> rotted_axyzs;
	for (l = 0; l < aress.size(); l++)
		rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
	
	this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
	
	for (l = 0; l < rotted_axyzs.size(); l++)
		delete[] rotted_axyzs[l];
	vector<double*>().swap(rotted_axyzs);
}


inline void CTMscoreComplex::align_multimer_normal_using_not_nwalign_and_not_greadsearch_with_max_rtmscore(const bool& fast){
	if (this->qsize != this->tsize){
		cout << endl;
		cout << "WRONG USAGE: the molecule number between two input pdbs are different: " << this->qsize << ", " << this->tsize << endl;
		cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
		cout << "|----->"<<TOOL_EXE<<" pdb1.pdb pdb2.pdb" << endl << endl;
	}
	
	int i, ii, j, jj, k, l, iali;
	char buf[2];
	string key, akey, bkey;
	
	this->obj_level_ali = new int[this->qsize];
	for (i = 0; i < this->qsize; i++){
		if (i < tsize) obj_level_ali[i] = i;
		else obj_level_ali[i] = -1;
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	int aligned_chain_num = 0;
	int seqali_res_num = 0;
	for (i = 0; i < this->qsize; i++){
		if (-1 != obj_level_ali[i]){
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[i];
			
			const string& aseq = amol->get_seq_str(); 
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			const MOLTYPE& bmt = bmol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			individual_tmscore_mtx[i][this->obj_level_ali[i]] = CBaseFunc::cal_rot_tran_from_query_to_templ__(axyzs, bxyzs, u, t, this->chain_index_corr_to_query__d0[i], a2b, fast);
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
							
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
	}
	
	rtmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
		
	if (DETAIL == g_print_result_type){
		int n = aress.size();
		for (i = 0; i < n; i++){
			__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
			this->aa_level_ali.push_back(__a2b__aa__[i]);
		}
	}
	
	// calculate TMscore
	vector<double*> rotted_axyzs;
	for (l = 0; l < aress.size(); l++)
		rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
	this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
	
	for (l = 0; l < rotted_axyzs.size(); l++)
		delete[] rotted_axyzs[l];
	vector<double*>().swap(rotted_axyzs);
}


inline void CTMscoreComplex::align_multimer_normal_using_not_nwalign_and_greadsearch(const bool& fast){
	int i, j, k, l, m, n, aL, bL;
	
	char buf[100];
	string key, akey, bkey;
	
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(this->qsize, this->tsize);
		for (i = 0; i < this->qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < this->tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				
				for (k = 0; k < this->qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE& kmt = kmol->get_moltype();
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < this->tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[this->qsize];
				int* transpose_ali = new int[this->tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, this->qsize, this->tsize, ali, transpose_ali);
				if (gs_score > corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < this->qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < this->qsize; k++){
							if (-1 == ali[k]) continue;
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const string& aseq = amol->get_seq_str();
							const string& bseq = bmol->get_seq_str();
							
							const MOLTYPE& amt = amol->get_moltype();
							const MOLTYPE& bmt = bmol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx( k, ali[k]);
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];	
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, this->qsize);
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < this->qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			const MOLTYPE& bmt = bmol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]]; 	
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		if (g_user_given_d0 <= 0)
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
		else tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, g_user_given_d0, g_user_given_d0, g_user_given_d0, g_user_given_d0, fast);
		tmscore = tmscore * seqali_res_num / total_res_num;
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate rTMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
		
		this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
}


inline void CTMscoreComplex::align_multimer_normal_using_not_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast){
	int i, j, k, l, m, n, aL, bL;
	
	char buf[100];
	string key, akey, bkey;
	
	if (NULL == this->obj_level_ali){
		double gs_score, best_sco = 0., corr_gs_score = 0.;
		double** scomtx = CBaseFunc::new2Darr(this->qsize, this->tsize);
		for (i = 0; i < this->qsize; i++){
			Molecule* imol = (*(this->query))[i];
			const vector<double*> ixyzs = imol->get_cared_xyz_vec();
			
			for (j = 0; j < this->tsize; j++){
				Molecule* jmol = (*(this->templ))[j];
				const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
				
				int* i2j = this->get_ij_qt_match_mtx(i, j);
				if (NULL == i2j)
					continue;
				
				// calculate the u and t
				individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
				
				for (k = 0; k < this->qsize; k++){
					Molecule* kmol = (*(this->query))[k];
					const MOLTYPE& kmt = kmol->get_moltype(); 
					const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
					vector<double*> roted_kxyzs;
					for (l = 0; l < kxyzs.size(); l++)
						roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
					
					for (l = 0; l < this->tsize; l++){
						Molecule* lmol = (*(this->templ))[l];
						const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
						
						int* k2l = this->get_ij_qt_match_mtx(k, l);
						if (NULL == k2l){
							scomtx[k][l] = 0.;
						}else{
							scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
						}
					}
					
					for (l = 0; l < roted_kxyzs.size(); l++)
						delete[] roted_kxyzs[l];
				}
				
				int* ali = new int[this->qsize];
				int* transpose_ali = new int[this->tsize];
				gs_score = CBaseFunc::greedySearch(scomtx, this->qsize, this->tsize, ali, transpose_ali);
				if (gs_score > corr_gs_score - 0.5){
					bool is_ok = true;
					for (k = 0; k < this->qsize; k++){
						if (-1 != ali[k]){
							if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
								is_ok = false;
								break;
							}
						}
					}
					
					if (is_ok){
						vector<double*> rotted_aress;
						vector<double*> bress;
						vector<MOLTYPE> mts;
						for (k = 0; k < this->qsize; k++){
							if (-1 == ali[k]) continue;
							Molecule* amol = (*(this->query))[k];
							Molecule* bmol = (*(this->templ))[ali[k]];
							
							const MOLTYPE& amt = amol->get_moltype();
							
							const vector<double*> axyzs = amol->get_cared_xyz_vec();
							const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
							
							int alen = axyzs.size();
							int blen = bxyzs.size();
							
							int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
							for (l = 0; l < alen; l++){
								if (-1 != a2b[l]){
									rotted_aress.push_back(CBaseFunc::rotateAndTrans(axyzs[l], u, t));
									bress.push_back(bxyzs[a2b[l]]);
									mts.push_back(amt);
								}
							}
						}
						
						double score1 = CBaseFunc::u3b_func(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig) / total_res_num;
						double score2 = CBaseFunc::score_fun_once(rotted_aress, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
						double score = score1>score2 ? score1 : score2;
						if (score > best_sco){
							best_sco = score;
							corr_gs_score = gs_score;
							if (NULL != this->obj_level_ali)
								delete[] this->obj_level_ali;
							this->obj_level_ali = ali;
							ali = NULL;
						}
						
						int nn = rotted_aress.size();
						for (k = 0; k < nn; k++)
							delete[] rotted_aress[k];					
					}				
				}
				
				if (NULL != ali)
					delete[] ali;
				delete[] transpose_ali;
			}
		}
		
		CBaseFunc::delete2Darr(scomtx, this->qsize);
	}
	
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query;
	vector<MOLTYPE> mts;
	vector<ALIGN_PAIR> __a2b__aa__;
	string achain, bchain;
	vector<string> aseq_vec, bseq_vec;
	if (NULL != this->obj_level_ali){
		int aligned_chain_num = 0;
		int seqali_res_num = 0;
		for (i = 0; i < this->qsize; i++){
			if (-1 == this->obj_level_ali[i]) continue;
			
			aligned_chain_num++;
			
			Molecule* amol = (*(this->query))[i];
			Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
			
			const string& aseq = amol->get_seq_str();
			const string& bseq = bmol->get_seq_str();
			
			const MOLTYPE& amt = amol->get_moltype();
			
			const vector<double*> axyzs = amol->get_cared_xyz_vec();
			const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
			
			int alen = axyzs.size();
			int blen = bxyzs.size();
			
			if (DETAIL == g_print_result_type){
				achain = this->query->get_chain(i);
				bchain = this->templ->get_chain(this->obj_level_ali[i]);
				
				if (LIGAND == amt){
					aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
					bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
				}
			}
			
			int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
			for (j = 0; j < alen; j++){
				if (-1 != a2b[j]){
					aress.push_back(axyzs[j]);
					bress.push_back(bxyzs[a2b[j]]);
					chain_index_corr_to_query.push_back(i);
					mts.push_back(amt);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = j;
						ap.tind = a2b[j];
						ap.qoind = amol->get_ith_orig_index(j);
						ap.toind = bmol->get_ith_orig_index(a2b[j]);
						ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(j);
						ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[j]);
						
						if (LIGAND == amt){
							ap.qaa = aseq_vec[j];
							ap.taa = bseq_vec[a2b[j]]; 
						}else{
							ap.qaa = aseq[j];
							ap.taa = bseq[a2b[j]];
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					seqali_res_num++;
				}
			}
		}
		
		rtmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
		
		if (DETAIL == g_print_result_type){
			int n = aress.size();
			for (i = 0; i < n; i++){
				__a2b__aa__[i].dis2 = CBaseFunc::distance2(aress[i], bress[i], u, t);
				this->aa_level_ali.push_back(__a2b__aa__[i]);
			}
		}
		
		// calculate TMscore
		vector<double*> rotted_axyzs;
		for (l = 0; l < aress.size(); l++)
			rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], u, t));
			
		this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
		
		for (l = 0; l < rotted_axyzs.size(); l++)
			delete[] rotted_axyzs[l];
		vector<double*>().swap(rotted_axyzs);
	}
}

inline void CTMscoreComplex::align_multimer_slow_but_accuracy_using_nwalign_and_greadsearch(const bool& fast){
	int i, j, k, l, m, n, iL;
	char buf[2];
	string key;
	
	double** u = CBaseFunc::new2Darr(3, 3);
	double* t = new double[3];
	
	this->tmscore = 0.;
	double gs_score, corr_gs_score = 0.;
	double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
	for (i = 0; i < qsize; i++){
		Molecule* imol = (*(this->query))[i];
		const vector<double*> ixyzs = imol->get_cared_xyz_vec();
		
		for (j = 0; j < tsize; j++){
			Molecule* jmol = (*(this->templ))[j];
			const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
			
			int* i2j = this->get_ij_qt_match_mtx(i, j);
			if (NULL == i2j)
				continue;
			
			individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
			
			for (k = 0; k < qsize; k++){
				Molecule* kmol = (*(this->query))[k];
				const MOLTYPE& kmt = kmol->get_moltype();
				const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
				vector<double*> roted_kxyzs;
				for (l = 0; l < kxyzs.size(); l++)
					roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
				
				for (l = 0; l < tsize; l++){
					Molecule* lmol = (*(this->templ))[l];
					const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
					
					int* k2l = this->get_ij_qt_match_mtx(k, l);
					if (NULL == k2l){
						scomtx[k][l] = 0.;
					}else{
						scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
					}
				}
				
				for (l = 0; l < roted_kxyzs.size(); l++)
					delete[] roted_kxyzs[l];
			}
			
			int* ali = new int[qsize];
			int* transpose_ali = new int[tsize];
			gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
			gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			
			if (gs_score >= corr_gs_score - 0.5){
				bool is_ok = true;
				for (k = 0; k < qsize; k++){
					if (-1 != ali[k]){
						if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
							is_ok = false;
							break;
						}
					}
				}
				
				if (is_ok){
					vector<double*> aress;
					vector<double*> bress;
					vector<MOLTYPE> mts;
					vector<ALIGN_PAIR> __a2b__aa__;
					string achain, bchain; 
					vector<string> aseq_vec, bseq_vec;
					int seqali_res_num = 0;
					for (k = 0; k < qsize; k++){
						if (-1 == ali[k]) continue;
						Molecule* amol = (*(this->query))[k];
						Molecule* bmol = (*(this->templ))[ali[k]];
						
						const string& aseq = amol->get_seq_str();
						const string& bseq = bmol->get_seq_str();
						
						const MOLTYPE& amt = amol->get_moltype();
						
						const vector<double*> axyzs = amol->get_cared_xyz_vec();
						const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
						
						int alen = axyzs.size();
						
						if (DETAIL == g_print_result_type){
							achain = this->query->get_chain(k);
							bchain = this->templ->get_chain(ali[k]);
							
							if (LIGAND == amt){
								aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
								bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
							}
						}
						
						int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
						for (l = 0; l < alen; l++){
							if (-1 != a2b[l]){
								aress.push_back(axyzs[l]);
								bress.push_back(bxyzs[a2b[l]]);
								mts.push_back(amt);
								
								if (DETAIL == g_print_result_type){
									ALIGN_PAIR ap;
									ap.qchain = achain;
									ap.tchain = bchain;
									ap.qind = l;
									ap.tind = a2b[l];
									ap.qoind = amol->get_ith_orig_index(l);
									ap.toind = bmol->get_ith_orig_index(a2b[l]);
									ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(l);
									ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[l]);
									
									if (LIGAND == amt){
										ap.qaa = aseq_vec[l];
										ap.taa = bseq_vec[a2b[l]]; 
									}else{
										ap.qaa = aseq[l];
										ap.taa = bseq[a2b[l]]; 	
									}
									
									__a2b__aa__.push_back(ap);
								}
								
								seqali_res_num++;
							}
						}
					}
					
					double score = CBaseFunc::cal_rot_tran_from_query_to_templ__(aress, bress, mts, u, t, d0_pro, d0_dna, d0_rna, d0_lig, fast);
					score = score * seqali_res_num / total_res_num;
					if (score > this->tmscore){
						this->tmscore = score;
						corr_gs_score = gs_score;
						if (NULL != this->obj_level_ali)
							delete[] this->obj_level_ali;
						this->obj_level_ali = ali;
						ali = NULL;
						
						if (NULL != this->u){
							CBaseFunc::delete2Darr(this->u, 3);
							delete[] this->t;
						}
						this->u = u;
						this->t = t;
						u = CBaseFunc::new2Darr(3, 3);
						t = new double[3];
						
						if (DETAIL == g_print_result_type){
							int n = aress.size();
							vector<ALIGN_PAIR>().swap(this->aa_level_ali);
							for (k = 0; k < n; k++){
								__a2b__aa__[k].dis2 = CBaseFunc::distance2(aress[k], bress[k], this->u, this->t);
								this->aa_level_ali.push_back(__a2b__aa__[k]);
							}
						}
					}
				}	
			}
			
			if (NULL != ali)
				delete[] ali;
			delete[] transpose_ali;
		}
	}
	
	// calculate rTMscore
	vector<double*> aress;
	vector<double*> bress;
	vector<int> chain_index_corr_to_query; 
	int aligned_chain_num = 0;
	int seqali_res_num = 0;
	for (i = 0; i < qsize; i++){
		if (-1 == this->obj_level_ali[i]) continue;
		
		aligned_chain_num++;
		
		Molecule* amol = (*(this->query))[i];
		Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
		
		const vector<double*> axyzs = amol->get_cared_xyz_vec();
		const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
		
		int alen = axyzs.size();
		int blen = bxyzs.size();
		
		int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
		for (j = 0; j < alen; j++){
			if (-1 != a2b[j]){
				aress.push_back(axyzs[j]);
				bress.push_back(bxyzs[a2b[j]]);
				chain_index_corr_to_query.push_back(i);
			}
		}
	}
	
	vector<double*> rotted_axyzs;
	for (l = 0; l < aress.size(); l++)
		rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], this->u, this->t));
	this->rtmscore = CBaseFunc::score_fun_rtmsco_once(aligned_chain_num, this->qsize, rotted_axyzs, bress, chain_index_corr_to_query, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d02);
	for (l = 0; l < rotted_axyzs.size(); l++)
		delete[] rotted_axyzs[l];
	vector<double*>().swap(rotted_axyzs);
	
	if (NULL != u){
		CBaseFunc::delete2Darr(u, 3);
		delete[] t;
	}
	CBaseFunc::delete2Darr(scomtx, qsize);
}


inline void CTMscoreComplex::align_multimer_slow_but_accuracy_using_nwalign_and_greadsearch_with_max_rtmscore(const bool& fast){
	int i, j, k, l, m, n, iL;
	char buf[2];
	string key;
	
	double** u = CBaseFunc::new2Darr(3, 3);
	double* t = new double[3];
	
	this->tmscore = 0.;
	double gs_score, corr_gs_score = 0.;
	double** scomtx = CBaseFunc::new2Darr(qsize, tsize);
	for (i = 0; i < qsize; i++){
		Molecule* imol = (*(this->query))[i];
		const vector<double*> ixyzs = imol->get_cared_xyz_vec();
		
		for (j = 0; j < tsize; j++){
			Molecule* jmol = (*(this->templ))[j];
			const vector<double*> jxyzs = jmol->get_cared_xyz_vec();
			
			int* i2j = this->get_ij_qt_match_mtx(i, j);
			if (NULL == i2j)
				continue;
		
			individual_tmscore_mtx[i][j] = CBaseFunc::cal_rot_tran_from_query_to_templ__(ixyzs, jxyzs, u, t, this->chain_index_corr_to_query__d0[i], i2j, fast);
			
			for (k = 0; k < qsize; k++){
				Molecule* kmol = (*(this->query))[k];
				const MOLTYPE& kmt = kmol->get_moltype();
				const vector<double*> kxyzs = kmol->get_cared_xyz_vec();
				vector<double*> roted_kxyzs;
				for (l = 0; l < kxyzs.size(); l++)
					roted_kxyzs.push_back(CBaseFunc::rotateAndTrans(kxyzs[l], u, t));
				
				for (l = 0; l < tsize; l++){
					Molecule* lmol = (*(this->templ))[l];
					const vector<double*> lxyzs = lmol->get_cared_xyz_vec();
					
					int* k2l = this->get_ij_qt_match_mtx(k, l);
					if (NULL == k2l){
						scomtx[k][l] = 0.;
					}else{
						scomtx[k][l] = CBaseFunc::rough_score(kmt, roted_kxyzs, lxyzs, k2l); 
					}
				}
				
				for (l = 0; l < roted_kxyzs.size(); l++)
					delete[] roted_kxyzs[l];
			}
			
			int* ali = new int[qsize];
			int* transpose_ali = new int[tsize];
			gs_score = CBaseFunc::greedySearch(scomtx, qsize, tsize, ali, transpose_ali);
			gs_score = CBaseFunc::__2merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			gs_score = CBaseFunc::__3merExchangedSearch(scomtx, qsize, tsize, ali, transpose_ali, gs_score, qsize);
			
			if (gs_score >= corr_gs_score - 0.5){
				bool is_ok = true;
				for (k = 0; k < qsize; k++){
					if (-1 != ali[k]){
						if (NULL == this->get_ij_qt_match_mtx(k, ali[k])){
							is_ok = false;
							break;
						} 
					}
				}
				
				if (is_ok){
					vector<double*> aress;
					vector<double*> bress;
					vector<int> chain_index_corr_to_query;
					vector<ALIGN_PAIR> __a2b__aa__;
					string achain, bchain; 
					vector<string> aseq_vec, bseq_vec;
					int aligned_chain_num = 0;
					int seqali_res_num = 0;
					for (k = 0; k < qsize; k++){
						if (-1 == ali[k]) continue;
						
						aligned_chain_num++;
						
						Molecule* amol = (*(this->query))[k];
						Molecule* bmol = (*(this->templ))[ali[k]];
						
						const string& aseq = amol->get_seq_str();
						const string& bseq = bmol->get_seq_str();
						
						const MOLTYPE& amt = amol->get_moltype();
						
						const vector<double*> axyzs = amol->get_cared_xyz_vec();
						const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
						
						int alen = axyzs.size();
						
						if (DETAIL == g_print_result_type){
							achain = this->query->get_chain(k);
							bchain = this->templ->get_chain(ali[k]);
							
							if (LIGAND == amt){
								aseq_vec = CBaseFunc::stringSplit(aseq, ' ');
								bseq_vec = CBaseFunc::stringSplit(bseq, ' ');
							} 
						}
						
						int* a2b = this->get_ij_qt_match_mtx(k, ali[k]);
						for (l = 0; l < alen; l++){
							if (-1 != a2b[l]){
								aress.push_back(axyzs[l]);
								bress.push_back(bxyzs[a2b[l]]);
								chain_index_corr_to_query.push_back(k); 
								
								if (DETAIL == g_print_result_type){
									ALIGN_PAIR ap;
									ap.qchain = achain;
									ap.tchain = bchain;
									ap.qind = l;
									ap.tind = a2b[l];
									ap.qoind = amol->get_ith_orig_index(l);
									ap.toind = bmol->get_ith_orig_index(a2b[l]);
									ap.qoindsuf = amol->get_ith_char_following_orig_index_vec(l);
									ap.toindsuf = bmol->get_ith_char_following_orig_index_vec(a2b[l]);
									
									if (LIGAND == amt){
										ap.qaa = aseq_vec[l];
										ap.taa = bseq_vec[a2b[l]]; 
									}else{
										ap.qaa = aseq[l];
										ap.taa = bseq[a2b[l]]; 	
									}
									
									__a2b__aa__.push_back(ap);
								}
								
								seqali_res_num++;
							}
						}
					}
					
					double score = CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
														aligned_chain_num, 
														this->qsize, 
														aress, 
														bress, 
														chain_index_corr_to_query__aa_num,
														chain_index_corr_to_query__moltype,
														chain_index_corr_to_query__d0,
														chain_index_corr_to_query__d02,
														chain_index_corr_to_query,
														u, t, fast);
					if (score > this->rtmscore){
						this->rtmscore = score;
						corr_gs_score = gs_score;
						if (NULL != this->obj_level_ali)
							delete[] this->obj_level_ali;
						this->obj_level_ali = ali;
						ali = NULL;
						
						if (NULL != this->u){
							CBaseFunc::delete2Darr(this->u, 3);
							delete[] this->t;
						}
						this->u = u;
						this->t = t;
						u = CBaseFunc::new2Darr(3, 3);
						t = new double[3];
						
						if (DETAIL == g_print_result_type){
							int n = aress.size();
							vector<ALIGN_PAIR>().swap(this->aa_level_ali);
							for (k = 0; k < n; k++){
								__a2b__aa__[k].dis2 = CBaseFunc::distance2(aress[k], bress[k], this->u, this->t);
								this->aa_level_ali.push_back(__a2b__aa__[k]);
							}
						}
					}
				}
			}
			
			if (NULL != ali)
				delete[] ali;
			delete[] transpose_ali;
		}
	}
	
	// calculate TMscore
	vector<double*> aress;
	vector<double*> bress;
	vector<MOLTYPE> mts;
	int aligned_chain_num = 0;
	int seqali_res_num = 0;
	for (i = 0; i < qsize; i++){
		if (-1 == this->obj_level_ali[i]) continue;
		
		aligned_chain_num++;
		
		Molecule* amol = (*(this->query))[i];
		Molecule* bmol = (*(this->templ))[this->obj_level_ali[i]];
		
		const MOLTYPE& amt = amol->get_moltype();
		
		const string& aseq = amol->get_seq_str();
		const string& bseq = bmol->get_seq_str();
		
		const vector<double*> axyzs = amol->get_cared_xyz_vec();
		const vector<double*> bxyzs = bmol->get_cared_xyz_vec();
		
		int alen = axyzs.size();
		int blen = bxyzs.size();
		
		int* a2b = this->get_ij_qt_match_mtx(i, this->obj_level_ali[i]);
		for (j = 0; j < alen; j++){
			if (-1 != a2b[j]){
				aress.push_back(axyzs[j]);
				bress.push_back(bxyzs[a2b[j]]);
				mts.push_back(amt);
			}
		}
	}
	
	vector<double*> rotted_axyzs;
	for (l = 0; l < aress.size(); l++)
		rotted_axyzs.push_back(CBaseFunc::rotateAndTrans(aress[l], this->u, this->t));
	this->tmscore = CBaseFunc::score_fun_once(rotted_axyzs, bress, mts, d02_pro, d02_dna, d02_rna, d02_lig, total_res_num);
	for (l = 0; l < rotted_axyzs.size(); l++)
		delete[] rotted_axyzs[l];
	vector<double*>().swap(rotted_axyzs);
	
	if (NULL != u){
		CBaseFunc::delete2Darr(u, 3);
		delete[] t;
	}
	CBaseFunc::delete2Darr(scomtx, qsize);
}


inline CTMscoreComplex::CTMscoreComplex(const string& qpdb, const string& tpdb, const ALIGN_TYPE& ali_type, const vector<CHAIN_PAIR>* user_chain_pair){
	clock_t start_t, end_t;
	start_t = clock();
	
	this->query = new Complex(qpdb);
	this->templ = new Complex(tpdb);
	this->qsize = query->size();
	this->tsize = templ->size();
	
	this->individual_tmscore_mtx = CBaseFunc::new2Darr(this->qsize, this->tsize);
	
	if (0 == qsize){
		cout << endl;
		cout << "*** Error ***: There is no cared type molecule in the " << qpdb << endl;
		cout << "Please check your inputted file of \"" << qpdb << "\"" << endl;
		cout << endl;
		exit(1);
	}else if (0 == tsize){
		cout << endl;
		cout << "*** Error ***: There is no cared type molecule in the " << tpdb << endl;
		cout << "Please check your inputted file of \"" << tpdb << "\"" << endl;
		exit(1);
	}
	
	this->u = CBaseFunc::new2Darr(3, 3);
	this->t = new double[3];
	this->obj_level_ali = NULL;
	if (NULL != user_chain_pair)
		this->set_obj_level_ali(*user_chain_pair);
	
	this->inv_u = NULL;
	this->inv_t = NULL;
	
	this->qt_match_mtx = NULL;
	this->q_ligAtomMatch_obj_vec = NULL;
	this->t_ligAtomMatch_obj_vec = NULL;
	
	if (1==qsize && 1==tsize){
		prepare_for_monomer();
		align_monomer(false);
	}else{
		prepare_for_multimer(g_use_res_orig_index_or_not, g_use_chain_order_or_not);
		if (g_use_rtmscore_to_search_best){
			// This function is outdated.
			if (g_use_chain_order_or_not){
				if (g_use_res_orig_index_or_not){
					// Without NW-align, Without Greedy-Search
					this->align_multimer_normal_using_not_nwalign_and_not_greadsearch_with_max_rtmscore(false);
				}else{
					// With NW-align, Without Greedy-Search
					this->align_multimer_normal_using_nwalign_and_not_greadsearch_with_max_rtmscore(false);
				}
			}else{
				if (g_use_res_orig_index_or_not){
					// Without NW-align, With Greedy-Search
					this->align_multimer_normal_using_not_nwalign_and_greadsearch_with_max_rtmscore(false);
				}else{
					// With NW-align, With Greedy-Search
					switch (ali_type){
						case SUPERFAST:
						case FAST:
							this->align_multimer_fast_buf_inaccuracy_using_nwalign_and_greadsearch_with_max_rtmscore(ali_type == SUPERFAST);
							break;
						case SLOW:
							this->align_multimer_slow_but_accuracy_using_nwalign_and_greadsearch_with_max_rtmscore(false); 
							break;
						default:
							if (g_fully_cared_ligatom_match)
//								this->align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore_fullycared_ligand(false);
								this->align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore_fullycared_ligand_III(false);
							else this->align_multimer_normal_using_nwalign_and_greadsearch_with_max_rtmscore(false);
							break;
					}
				}
			}
		}else{
			if (g_use_chain_order_or_not){
				if (g_use_res_orig_index_or_not){
					// Without NW-align, Without Greedy-Search
					this->align_multimer_normal_using_not_nwalign_and_not_greadsearch(false);
				}else{
					// With NW-align, Without Greedy-Search
					this->align_multimer_normal_using_nwalign_and_not_greadsearch(false);
				}
			}else{
				if (g_use_res_orig_index_or_not){
					// Without NW-align, With Greedy-Search
					this->align_multimer_normal_using_not_nwalign_and_greadsearch(false);
				}else{
					// With NW-align, With Greedy-Search
					switch (ali_type){
						case SUPERFAST:
						case FAST:
							this->align_multimer_fast_buf_inaccuracy_using_nwalign_and_greadsearch(ali_type == SUPERFAST);
							break;
						case SLOW:
							this->align_multimer_slow_but_accuracy_using_nwalign_and_greadsearch(false);
							break;
						default:
							if (g_fully_cared_ligatom_match)
//								this->align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand(false);
								this->align_multimer_normal_using_nwalign_and_greadsearch_fullycared_ligand_III(false);
							else this->align_multimer_normal_using_nwalign_and_greadsearch(false);
							break;
					}
				}
			}
		}
	}
	
	end_t = clock();
	this->use_seconds = 1.0*(end_t - start_t)/CLOCKS_PER_SEC;
}

inline void CTMscoreComplex::generate_atom_alignment_of_each_ligand_pair_using_index_order(const bool& use_chain_order_or_not){
	if (use_chain_order_or_not){
		int i, j, k;
		qt_match_mtx = new int**[qsize];
		for (i = 0; i < qsize; i++){
			qt_match_mtx[i] = NULL;
			if (i >= tsize) continue;
			
			Molecule* qmol = (*query)[i];
			Molecule* tmol = (*templ)[i];
			const MOLTYPE& qmt = qmol->get_moltype();
			const MOLTYPE& tmt = tmol->get_moltype();
			const string& qseq = qmol->get_seq_str();
			const string& tseq = tmol->get_seq_str();
			const int& isize = qmol->size();
			
			if (LIGAND == qmt){
				if (LIGAND != tmt){
					cout << "ERROR: The types of the " << (i+1) << "-th molecules in the 1st and 2nd complexes are different." << endl;
					cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
					cout << "|----->"<< TOOL_EXE << " pdb1.pdb pdb2.pdb" << endl;
					exit(1);
				}
				
				if (qseq != tseq){
					cout << "ERROR: The atom topology information of the " << (i+1) << " molecules in the 1st and 2nd complexes are different." << endl;
					cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
					cout << "|----->"<< TOOL_EXE << " pdb1.pdb pdb2.pdb" << endl;
					exit(1);
				}
				
				qt_match_mtx[i] = new int*[tsize];
				for (j = 0; j < tsize; j++)
					qt_match_mtx[i][j] = NULL;
					
				qt_match_mtx[i][i] = new int[isize];
				for (k = 0; k < isize; k++)
					qt_match_mtx[i][i][k] = k;
			}
		}
	}else{
		int i, j, k;
		qt_match_mtx = new int**[qsize];
		for (i = 0; i < qsize; i++){
			qt_match_mtx[i] = NULL;
			
			Molecule* imol = (*query)[i];
			const MOLTYPE& imt = imol->get_moltype();
			const string& iseq = imol->get_seq_str();
			const int& isize = imol->size();
			if (LIGAND == imt){
				qt_match_mtx[i] = new int*[tsize];
				for (j = 0; j < tsize; j++){
					qt_match_mtx[i][j] = NULL;
					
					Molecule* jmol = (*templ)[j];
					const MOLTYPE& jmt = jmol->get_moltype();
					const string& jseq = jmol->get_seq_str();
					if (LIGAND == jmt){
						if (iseq == jseq){
							qt_match_mtx[i][j] = new int[isize];
							for (k = 0; k < isize; k++)
								qt_match_mtx[i][j][k] = k;
						}
					}
				}
			}
		}
	}
}

inline void CTMscoreComplex::generate_atom_alignment_of_each_ligand_pair_using_greedysearch(const bool& use_chain_order_or_not){
	if (use_chain_order_or_not){
		int i, j;
		double** u = CBaseFunc::new2Darr(3, 3);
		double* t = new double[3];
		
		q_ligAtomMatch_obj_vec = new LigAtomMatcher*[qsize];
		for (i = 0; i < qsize; i++)
			q_ligAtomMatch_obj_vec[i] = NULL;
		 
		t_ligAtomMatch_obj_vec = new LigAtomMatcher*[tsize];
		for (i = 0; i < tsize; i++)
			t_ligAtomMatch_obj_vec[i] = NULL;
		
		qt_match_mtx = new int**[qsize];
		for (i = 0; i < qsize; i++){
			qt_match_mtx[i] = NULL;
			if (i >= tsize) continue;
			
			Molecule* qmol = (*query)[i];
			Molecule* tmol = (*templ)[i];
			const MOLTYPE& qmt = qmol->get_moltype();
			const MOLTYPE& tmt = tmol->get_moltype();
			const int& isize = qmol->size();
			
			if (LIGAND == qmt){
				if (LIGAND != tmt){
					cout << "ERROR: The types of the " << (i+1) << "th molecules in the 1st and 2nd complexes are different." << endl;
					cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
					cout << "|----->"<< TOOL_EXE << " pdb1.pdb pdb2.pdb" << endl;
					exit(1);
				}
				
				LigAtomMatcher* qamt = new LigAtomMatcher(qmol);
				LigAtomMatcher* tamt = new LigAtomMatcher(tmol);
				
				if (!qamt->match(*tamt)){
					cout << "ERROR: The atom topology information of the " << (i+1) << "-th molecules in the 1st and 2nd complexes are different." << endl;
					cout << "|-- IF you do not know the chain alignment, you can directly use the defualt parameters, for example:" << endl;
					cout << "|----->"<< TOOL_EXE << " pdb1.pdb pdb2.pdb" << endl;
					exit(1);
				}
				
				qt_match_mtx[i] = new int*[tsize];
				for (j = 0; j < tsize; j++)
					qt_match_mtx[i][j] = NULL;
					
				qt_match_mtx[i][i] = new int[isize];
				LigAtomMatcher::quick_identity_atom_align(*qamt, *tamt, qt_match_mtx[i][i], u, t, 4.0);
				q_ligAtomMatch_obj_vec[i] = qamt;
				t_ligAtomMatch_obj_vec[i] = tamt;
			}
		}
		
		CBaseFunc::delete2Darr(u, 3);
		delete[] t;
	}else{
		int i;
		q_ligAtomMatch_obj_vec = new LigAtomMatcher*[qsize];
		for (int i = 0; i < qsize; i++){
			Molecule* imol = (*query)[i];
			const MOLTYPE& imt = imol->get_moltype();
			if (imt == LIGAND)
				q_ligAtomMatch_obj_vec[i] = new LigAtomMatcher(imol);
			else q_ligAtomMatch_obj_vec[i] = NULL;
		}
		
		t_ligAtomMatch_obj_vec = new LigAtomMatcher*[tsize];
		for (int i = 0; i < tsize; i++){
			Molecule* imol = (*templ)[i];
			const MOLTYPE& imt = imol->get_moltype();
			if (imt == LIGAND)
				t_ligAtomMatch_obj_vec[i] = new LigAtomMatcher(imol);
			else t_ligAtomMatch_obj_vec[i] = NULL;
		}
		
		double** u = CBaseFunc::new2Darr(3, 3);
		double* t = new double[3];
		
		qt_match_mtx = new int**[qsize];
		for (int i = 0; i < qsize; i++){
			qt_match_mtx[i] = NULL;
			
			LigAtomMatcher* ith = q_ligAtomMatch_obj_vec[i];
			if (NULL != ith){
				qt_match_mtx[i] = new int*[tsize];
				int isize = ith->size(); 
				for (int j = 0; j < tsize; j++){
					qt_match_mtx[i][j] = NULL;
					
					LigAtomMatcher* jth = t_ligAtomMatch_obj_vec[j];
					if (NULL != jth){
						if (ith->match(*jth)){
							qt_match_mtx[i][j] = new int[isize];
							LigAtomMatcher::quick_identity_atom_align(*ith, *jth, qt_match_mtx[i][j], u, t, 4.0);
						}
					}
				}
			}
		}
		
		CBaseFunc::delete2Darr(u, 3);
		delete[] t;
	}
}

inline int* CTMscoreComplex::get_ij_qt_match_mtx(const int& i, const int& j){
	if (qt_match_mtx[i] == NULL)
		return NULL;
	else return qt_match_mtx[i][j];
}

inline const double& CTMscoreComplex::get_use_seconds(){
	return this->use_seconds;
}

inline void CTMscoreComplex::align_monomer(const bool& fast){
	int i, j, k, l, iali;
	char c1, c2;
	double d0, seqid;
	int identity_num;
	
	Molecule* qmol = (*query)[0];
	Molecule* tmol = (*templ)[0];
	const string& qseq = qmol->get_seq_str();
	const string& tseq = tmol->get_seq_str();
	const MOLTYPE& qmt = qmol->get_moltype();
	const MOLTYPE& tmt = tmol->get_moltype();
	int qmol_len = qmol->size();
	int tmol_len = tmol->size();
	const vector<double*>& qxyz_vec = qmol->get_cared_xyz_vec();
	const vector<double*>& txyz_vec = tmol->get_cared_xyz_vec();
	
	if (qmt != tmt){
		cout << "The molecule typle bewteen query and templ pdb files are different (" << qmol->get_moltype() << "," << tmol->get_moltype() << ")." << endl;
		exit(1);
	}
	
	if (qmt == PROTEIN)
		d0 = CBaseFunc::d0_of_tmscore(qmol_len);
	else if (qmt == DNA)
		d0 = CBaseFunc::d0_of_tmscore_c3prime(qmol_len);
	else if (qmt == RNA)
		d0 = CBaseFunc::d0_of_tmscore_c3prime(qmol_len);
	else d0 = CBaseFunc::d0_of_lsscore(qmol_len);
	
	this->obj_level_ali = new int[1];
	this->obj_level_ali[0] = 0;
	string achain, bchain;
	vector<string> aseq_vec;
	vector<string> bseq_vec;
	if (DETAIL == g_print_result_type){
		achain = this->query->get_chain(0);
		bchain = this->templ->get_chain(this->obj_level_ali[0]);
		
		if (LIGAND == qmt){
			aseq_vec = qmol->get_cared_atomtype_vec_in_lig();
			bseq_vec = tmol->get_cared_atomtype_vec_in_lig();
		}	
	}
	vector<ALIGN_PAIR> __a2b__aa__;
	if (g_use_res_orig_index_or_not){
		identity_num = 0;
		if (LIGAND == qmt){
			
			if (qseq != tseq){
				cout << "The ligand type and the atom order information between query and templ pdb files are not same:" << endl;
				cout << "  Query: " << qseq << endl;
				cout << "  Templ: " << tseq << endl;
				exit(1);
			}
			
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(qxyz_vec, txyz_vec, this->u, this->t, d0, fast);
			rtmscore = tmscore;
			individual_tmscore_mtx[0][0] = tmscore;
			
			if (DETAIL == g_print_result_type){
				for (i = 0; i < qmol_len; i++){
					ALIGN_PAIR ap;
					ap.qchain = achain;
					ap.tchain = bchain;
					ap.qind = i;
					ap.tind = i;
					ap.qoind = qmol->get_ith_orig_index(i);
					ap.toind = tmol->get_ith_orig_index(i);
					ap.qoindsuf = qmol->get_ith_char_following_orig_index_vec(i);
					ap.toindsuf = tmol->get_ith_char_following_orig_index_vec(i);
					ap.qaa = aseq_vec[i];
					ap.taa = bseq_vec[i];
					__a2b__aa__.push_back(ap);
				}
				
				for (i = 0; i < qmol_len; i++){
					__a2b__aa__[i].dis2 = CBaseFunc::distance2(qxyz_vec[i], txyz_vec[i], this->u, this->t);
					this->aa_level_ali.push_back(__a2b__aa__[i]);
				}
			}
		}else {
			int* ali = new int[qmol_len];
			for (i = 0; i < qmol_len; i++){
				k = qmol->get_ith_orig_index(i);
				c1 = qmol->get_ith_char_following_orig_index_vec(i);
				iali = -1;
				for (j = 0; j < tmol_len; j++){
					l = tmol->get_ith_orig_index(j);
					c2 = tmol->get_ith_char_following_orig_index_vec(j);
					if (k == l && c1 == c2){
						iali = j;
						break;
					}
				}
				ali[i] = iali;
				if (-1 != iali) identity_num++;
			}
			
			vector<double*> __qxyz_vec;
			vector<double*> __txyz_vec;
			for (i = 0; i < qmol_len; i++){
				if (-1 != ali[i]){
					__qxyz_vec.push_back(qxyz_vec[i]);
					__txyz_vec.push_back(txyz_vec[ali[i]]);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = i;
						ap.tind = ali[i];
						ap.qoind = qmol->get_ith_orig_index(i);
						ap.toind = tmol->get_ith_orig_index(ali[i]);
						ap.qoindsuf = qmol->get_ith_char_following_orig_index_vec(i);
						ap.toindsuf = tmol->get_ith_char_following_orig_index_vec(ali[i]);
						ap.qaa = qseq[i];
						ap.taa = tseq[ali[i]];
						__a2b__aa__.push_back(ap);
					}
				}
			}
						
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(__qxyz_vec, __txyz_vec, this->u, this->t, d0, fast);
			tmscore = tmscore * identity_num / qmol_len;
			rtmscore = tmscore;
			individual_tmscore_mtx[0][0] = tmscore;
			
			if (DETAIL == g_print_result_type){
				int n = __qxyz_vec.size();
				for (i = 0; i < n; i++){
					__a2b__aa__[i].dis2 = CBaseFunc::distance2(__qxyz_vec[i], __txyz_vec[i], this->u, this->t);
					this->aa_level_ali.push_back(__a2b__aa__[i]);
				}
			}
			delete[] ali;
		}
	}else{ // g_use_res_orig_index_or_not = false
		if (LIGAND == qmt){
			LigAtomMatcher qmatcher(qmol);
			LigAtomMatcher tmatcher(tmol);
			if (!qmatcher.match(tmatcher)){
				cout << "WRONG USAGE: the ligand atoms or their topology information are not same." << endl;
				exit(1);
			}
			
			int* ali = new int[qmol_len];
			LigAtomMatcher::quick_identity_atom_align(qmatcher, tmatcher, ali, this->u, this->t, d0*d0);
			identity_num = qmol_len;
			
			vector<double*> __qxyz_vec;
			vector<double*> __txyz_vec;
			for (i = 0; i < qmol_len; i++){
				if (-1 != ali[i]){
					__qxyz_vec.push_back(qxyz_vec[i]);
					__txyz_vec.push_back(txyz_vec[ali[i]]);
					
					if (DETAIL == g_print_result_type){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = i;
						ap.tind = ali[i];
						ap.qoind = qmol->get_ith_orig_index(i);
						ap.toind = tmol->get_ith_orig_index(ali[i]);
						ap.qoindsuf = qmol->get_ith_char_following_orig_index_vec(i);
						ap.toindsuf = tmol->get_ith_char_following_orig_index_vec(ali[i]);
						ap.qaa = aseq_vec[i];
						ap.taa = bseq_vec[ali[i]];
						__a2b__aa__.push_back(ap);
					}
				}
			}
						
			tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(__qxyz_vec, __txyz_vec, this->u, this->t, d0, fast);
			tmscore = tmscore * identity_num / qmol_len;
			rtmscore = tmscore;
			individual_tmscore_mtx[0][0] = tmscore;
			
			if (DETAIL == g_print_result_type){
				int n = __qxyz_vec.size();
				for (i = 0; i < n; i++){
					__a2b__aa__[i].dis2 = CBaseFunc::distance2(__qxyz_vec[i], __txyz_vec[i], this->u, this->t);
					this->aa_level_ali.push_back(__a2b__aa__[i]);
				}
			}
			
			delete[] ali;
		} else { // LIGAND != qmt
			if (qseq == tseq){
				tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(qxyz_vec, txyz_vec, this->u, this->t, d0, fast);
				rtmscore = tmscore;
				individual_tmscore_mtx[0][0] = tmscore;
				
				if (DETAIL == g_print_result_type){
					for (i = 0; i < qmol_len; i++){
						ALIGN_PAIR ap;
						ap.qchain = achain;
						ap.tchain = bchain;
						ap.qind = i;
						ap.tind = i;
						ap.qoind = qmol->get_ith_orig_index(i);
						ap.toind = tmol->get_ith_orig_index(i);
						ap.qoindsuf = qmol->get_ith_char_following_orig_index_vec(i);
						ap.toindsuf = tmol->get_ith_char_following_orig_index_vec(i);
						
						if (LIGAND == qmt){
							ap.qaa = aseq_vec[i];
							ap.taa = bseq_vec[i];
						}else{
							ap.qaa = qseq[i];
							ap.taa = tseq[i];
						}
						
						__a2b__aa__.push_back(ap);
					}
					
					for (i = 0; i < qmol_len; i++) {
						__a2b__aa__[i].dis2 = CBaseFunc::distance2(qxyz_vec[i], txyz_vec[i], this->u, this->t);
						this->aa_level_ali.push_back(__a2b__aa__[i]);
					}
				}
			}else{ // qseq != tseq
				const vector<char>& qvec_seq = qmol->get_seq_vec();
				const vector<int>&  qoind = qmol->get_orig_index_vec();
				const vector<char>& tvec_seq = tmol->get_seq_vec();
				const vector<int>&  toind = tmol->get_orig_index_vec();
				
				int* __i2j_x_ = new int[qmol_len];
				seqid = CNWalign::nwalign(qvec_seq, qoind, tvec_seq, toind, qmt, __i2j_x_, identity_num);
				if (seqid < g_seqid_cutoff || identity_num < 4) {
					cout << "The molecules of two inputted files are not identical" << endl;
					cout << "  Query: " << qseq << endl;
					cout << "  Templ: " << tseq << endl;
					exit(1);
				}
				
				tmscore = CBaseFunc::cal_rot_tran_from_query_to_templ__(qxyz_vec, txyz_vec, this->u, this->t, d0, __i2j_x_, fast);
				tmscore = tmscore * identity_num / qmol_len;
				rtmscore = tmscore;	
				individual_tmscore_mtx[0][0] = tmscore;
				
				if (DETAIL == g_print_result_type){
					const int* ali = __i2j_x_;
					vector<double*> __qxyz_vec;
					vector<double*> __txyz_vec;
					
					for (i = 0; i < qmol_len; i++){
						if (-1 != ali[i]){
							__qxyz_vec.push_back(qxyz_vec[i]);
							__txyz_vec.push_back(txyz_vec[ali[i]]);
							
							ALIGN_PAIR ap;
							ap.qchain = achain;
							ap.tchain = bchain;
							ap.qind = i;
							ap.tind = ali[i];
							ap.qoind = qmol->get_ith_orig_index(i);
							ap.toind = tmol->get_ith_orig_index(ali[i]);
							ap.qoindsuf = qmol->get_ith_char_following_orig_index_vec(i);
							ap.toindsuf = tmol->get_ith_char_following_orig_index_vec(ali[i]);
							ap.qaa = qseq[i];
							ap.taa = tseq[ali[i]];
							__a2b__aa__.push_back(ap);
						}
					}
					
					int n = __qxyz_vec.size();
					for (i = 0; i < n; i++){
						__a2b__aa__[i].dis2 = CBaseFunc::distance2(__qxyz_vec[i], __txyz_vec[i], this->u, this->t);
						this->aa_level_ali.push_back(__a2b__aa__[i]);
					}
				}
				
				delete[] __i2j_x_;
			}
		}
	}
}

inline CTMscoreComplex::~CTMscoreComplex(){
	int i, j, k, n;
	
	delete query;
	delete templ;
	CBaseFunc::delete2Darr(u, 3);
	delete[] t;
	delete[] obj_level_ali;
	CBaseFunc::delete2Darr(individual_tmscore_mtx, this->qsize);
	
	if (NULL != inv_u){
		CBaseFunc::delete2Darr(inv_u, 3);
		delete[] inv_t;
	}
	
	if (NULL != qt_match_mtx){
		for (i = 0; i < qsize; i++){
			if (NULL != qt_match_mtx[i]){
				for (j = 0; j < tsize; j++){
					if (NULL != qt_match_mtx[i][j])
						delete[] qt_match_mtx[i][j];
				}
				delete[] qt_match_mtx[i];	
			}
		}
		delete[] qt_match_mtx;
	}
	
	if (NULL != q_ligAtomMatch_obj_vec){
		for (i = 0; i < qsize; i++){
			if (NULL != q_ligAtomMatch_obj_vec[i]){
				delete q_ligAtomMatch_obj_vec[i];
			}
		}
		delete[] q_ligAtomMatch_obj_vec;
	}
	
	if (NULL != t_ligAtomMatch_obj_vec){
		for (i = 0; i < tsize; i++){
			if (NULL != t_ligAtomMatch_obj_vec[i]){
				delete t_ligAtomMatch_obj_vec[i];
			}
		}
		delete[] t_ligAtomMatch_obj_vec;
	}
	
	vector<ALIGN_PAIR>().swap(this->aa_level_ali);
}

inline void CTMscoreComplex::set_obj_level_ali(const vector<CHAIN_PAIR>& chainpair){
	if (NULL != this->obj_level_ali)
		delete[] this->obj_level_ali;
	this->obj_level_ali = new int[qsize];
	for (int i = 0; i < qsize; i++)
		this->obj_level_ali[i] = -1;
	
	int* check_templ_ali = new int[tsize];
	for (int i = 0; i < tsize; i++)
		check_templ_ali[i] = -1;	
	
	int n = chainpair.size();
	int chain_ali_num = 0;
	for (int i = 0; i < n; i++){
		CHAIN_PAIR cp = chainpair[i];
		int qind = -1;
		for (int j = 0; j < qsize; j++){
			if (cp.qchain == query->get_chain(j)){
				qind = j;
				break;
			}
		}
		
		if (qind == -1){
			cout << "ERROR : the input chain (" << cp.qchain << ") does not exists in the corresponding input PDB file." << endl;
			exit(1);
		}else if (-1 != this->obj_level_ali[qind]){
			cout << "ERROR : there are two chains aligned the chain of (" << cp.qchain << ")." << endl;
			exit(1);
		}
		
		int tind = -1;
		for (int j = 0; j < tsize; j++){
			if (cp.tchain == templ->get_chain(j)){
				tind = j;
				break;
			}
		}
		
		if (tind == -1){
			cout << "ERROR : the input chain (" << cp.tchain << ") does not exists in the corresponding input PDB file." << endl;
			exit(1);
		}else if (-1 != check_templ_ali[tind]){
			cout << "ERROR : there are two chains aligned the chain of (" << cp.tchain << ")." << endl;
			exit(1);
		}
		
		check_templ_ali[tind] = qind;
		this->obj_level_ali[qind] = tind;
		chain_ali_num++;
	}
	
	if (0 == chain_ali_num){
		cout << "ERROR : You do not give any useful alignment in your input option of '-ia'" << endl;
		exit(1);
	}
	
	delete[] check_templ_ali;
} 

inline CNWalign::CNWalign(const string& aseq, const string& bseq, const MOLTYPE& mol_type) : aseq(aseq), bseq(bseq), mol_type(mol_type){
	this->alen = aseq.size();
	this->blen = bseq.size();
	
	if (this->mol_type == PROTEIN){
		gap_open=-11, gap_extn=-1;
		load_blosum62_for_protein();
	}else{ // dna or rna 
		gap_open=-5; gap_extn =-2;
		load_blosumn_for_nucletide();
	}
	
	identity_ali = NULL;
	identity_num = 0;
	run_needleman_wunsch();
}

inline double CNWalign::nwalign(const vector<char>& aseq, const vector<int>& aoind, const vector<char>& bseq, const vector<int>& boind, const MOLTYPE& mol_type, int* out_ali, int& identity_ali_num){
	int i, j;
	int alen = aseq.size();
	int blen = bseq.size();
	int* aind = new int[alen];
	int* bind = new int[blen];
	
	vector<char> new_aseq;
	vector<char> new_bseq;
	
	int pre_aoind = -1e9;
	for (i= 0; i < alen; i++){
		if (-1e9 != pre_aoind){
			int gapn = aoind[i] - pre_aoind - 1;
			if (0 < gapn){
				if (PROTEIN == mol_type){
					for (j = 0; j < gapn; j++)
						new_aseq.push_back('X');
				} else {
					for (j = 0; j < gapn; j++)
						new_aseq.push_back('x');
				}
			}
		}
		
		pre_aoind = aoind[i];
		aind[i] = new_aseq.size();
		new_aseq.push_back(aseq[i]);
	}
	
	int pre_boind = -1e9;
	for (i= 0; i < blen; i++){
		if (-1e9 != pre_boind) {
			int gapn = boind[i] - pre_boind - 1;
			if (0 < gapn){
				if (PROTEIN == mol_type){
					for (j = 0; j < gapn; j++)
						new_bseq.push_back('X');
				} else {
					for (j = 0; j < gapn; j++)
						new_bseq.push_back('x');
				}
			}
		}
		
		pre_boind = boind[i];
		bind[i] = new_bseq.size();
		new_bseq.push_back(bseq[i]);
	}
	
	int nalen = new_aseq.size();
	int* naind = new int[nalen];
	for (i = 0; i < nalen; i++)
		naind[i] = -1;
	for (i = 0; i < alen; i++)
		naind[aind[i]] = i;
	
	int nblen = new_bseq.size();
	int* nbind = new int[nblen];
	for (i = 0; i < nblen; i++)
		nbind[i] = -1;
	for (i = 0; i < blen; i++)
		nbind[bind[i]] = i;
		
	
	CNWalign nwobj(CBaseFunc::charVec2String(new_aseq), CBaseFunc::charVec2String(new_bseq), mol_type); 
	const int* new_ali = nwobj.get_identity_ali();	
	
	for (i = 0; i < alen; i++)
		out_ali[i] = -1;
	
	for (i = 0; i < nalen; i++){
		if (-1 != new_ali[i] && -1 != naind[i]){
			out_ali[ naind[i] ] = nbind[ new_ali[i] ];
		}
	}
	
	identity_ali_num = 0;
	for (i = 0; i < alen; i++){
		if (-1 != out_ali[i]){
			if (aseq[i] == bseq[out_ali[i]])
				identity_ali_num++;
		}
	}
	
	delete[] naind;
	delete[] nbind;
	delete[] aind;
	delete[] bind;
	
	return 1.*identity_ali_num / (alen<blen?alen:blen);
}

inline CNWalign::~CNWalign(){
	if (NULL != identity_ali){
		delete[] identity_ali;
		identity_ali = NULL;
	}
}

inline void CNWalign::run_needleman_wunsch(){
	int i, j;
	double D, V, H;
	int alen = this->alen+1;
	int blen = this->blen+1;
	
	int* seq1 = new int[alen];
	int* seq2 = new int[blen];
	double** score = CBaseFunc::new2Darr(alen, blen);
	double** val = CBaseFunc::new2Darr(alen+1, blen+1);
	int** idir = CBaseFunc::new2DIntArr(alen+1, blen+1);
	double** preV = CBaseFunc::new2Darr(alen+1, blen+1);
	double** preH = CBaseFunc::new2Darr(alen+1, blen+1);
	int** jpV = CBaseFunc::new2DIntArr(alen+1, blen+1);
	int** jpH = CBaseFunc::new2DIntArr(alen+1, blen+1);
	
	for (i = 1; i < alen; i++)
		for (j = 0; j < imut_seq_len; j++)
			if (aseq[i-1] == imut_seq[j]){
				seq1[i] = j;
				break;
			}
	
	for (i = 1; i < blen; i++)
		for (j = 0; j < imut_seq_len; j++)
			if (bseq[i-1] == imut_seq[j]){
				seq2[i] = j;
				break;
			}
	
	for (i = 1; i < alen; i++)
		for (j = 1; j < blen; j++)
			score[i][j] = imut[seq1[i]][seq2[j]];
	
	val[0][0] = 0;
	val[1][0] = gap_open;
	for (i = 2; i < alen; i++)
		val[i][0] = val[i - 1][0] + gap_extn;
	for (i = 1; i < alen; i++) {
		preV[i][0] = val[i][0]; // not use preV at the beginning
		idir[i][0] = 0; // useless
		jpV[i][0] = 1; // useless
		jpH[i][0] = i; // useless
	}
	val[0][1] = gap_open;
	for (j = 2; j < blen; j++)
		val[0][j] = val[0][j - 1] + gap_extn;
	for (j = 1; j < blen; j++) {
		preH[0][j] = val[0][j];
		idir[0][j] = 0;
		jpV[0][j] = j;
		jpH[0][j] = 1;
	}

	// DP ------------------------------>
	for (j = 1; j < blen; j++) {
		for (i = 1; i < alen; i++) {
			// D=VAL(i-1,j-1)+SCORE(i,j)--------------->
			D = val[i - 1][j - 1] + score[i][j]; // from diagonal, val(i,j) is val(i-1,j-1)
			D = i!=1||j!=1 ? D+g_eps : D;
			
			// H=H+gap_open ------->
			jpH[i][j] = 1;
			double val1 = val[i - 1][j] + gap_open; // gap_open from both D and V
			double val2 = preH[i - 1][j] + gap_extn; // gap_extn from horizontal
			if (val1 > val2) // last step from D or V
			{
				H = val1;
			} else {// last step from H
				H = val2;
				if (i > 1) {
					jpH[i][j] = jpH[i - 1][j] + 1; // record long-gap
				}
			}

			// V=V+gap_open --------->
			jpV[i][j] = 1;
			val1 = val[i][j - 1] + gap_open;
			val2 = preV[i][j - 1] + gap_extn;
			if (val1 > val2) {
				V = val1;
			} else {
				V = val2;
				if (j > 1) {
					jpV[i][j] = jpV[i][j - 1] + 1; // record long-gap
				}
			}

			preH[i][j] = H; // unaccepted H
			preV[i][j] = V; // unaccepted V
			if ((D > H) && (D > V)) {
				idir[i][j] = 1;
				val[i][j] = D;
			} else if (H > V) {
				idir[i][j] = 2;
				val[i][j] = H;
			} else {
				idir[i][j] = 3;
				val[i][j] = V;
			}
		}
	}
	
	this->identity_num = 0;
	this->identity_ali = new int[this->alen];
	for (i = 0; i < this->alen; i++)
		this->identity_ali[i] = -1;
	
	// tracing back the pathway
	i = alen-1;
	j = blen-1;
	while ((i > 0) && (j > 0)) {
		if (idir[i][j] == 1){ // from diagonal
//			if (aseq[i-1] == bseq[j-1]){
//				this->identity_num++;
//				this->identity_ali[i-1] = j - 1;
//			}
			
			// like US-align's option "-TM-score 7", it is not identity
			this->identity_num++;
			this->identity_ali[i-1] = j - 1;
			
			i = i - 1;
			j = j - 1;
		} else if (idir[i][j] == 2){ // from horizonal
			int temp1 = jpH[i][j]; //
			for (int me = 0; me < temp1; me++)
				if (i > 0)
					i = i - 1;
		} else {
			int temp2 = jpV[i][j];
			for (int me = 0; me < temp2; me++)
				if (j > 0)
					j = j - 1;
		}
	}

	delete[] seq1;
	delete[] seq2;
	CBaseFunc::delete2Darr(score, alen);
	CBaseFunc::delete2Darr(val, alen+1);
	CBaseFunc::delete2DIntArr(alen+1, idir);
	CBaseFunc::delete2Darr(preV, alen+1);
	CBaseFunc::delete2Darr(preH, alen+1);
	CBaseFunc::delete2DIntArr(alen+1, jpV);
	CBaseFunc::delete2DIntArr(alen+1, jpH);
	
}

inline const int* CNWalign::get_identity_ali(){
	return identity_ali;	
}

inline const int& CNWalign::operator [](const int& i){
	return identity_ali[i];
}

inline const double CNWalign::identity_div_min_len(){
	int min_len = alen>blen?blen:alen;
	return 1.*identity_num/min_len;
}

inline void CNWalign::load_blosumn_for_nucletide(){
	imut_seq = "acgtux";
	imut_seq_len = imut_seq.size();
	imut[0][0] = 2;
	imut[0][1] = -3;
	imut[0][2] = -3;
	imut[0][3] = -3;
	imut[0][4] = -3;
	imut[0][5] = -3;
	imut[1][0] = -3;
	imut[1][1] = 2;
	imut[1][2] = -3;
	imut[1][3] = -3;
	imut[1][4] = -3;
	imut[1][5] = -3;
	imut[2][0] = -3;
	imut[2][1] = -3;
	imut[2][2] = 2;
	imut[2][3] = -3;
	imut[2][4] = -3;
	imut[2][5] = -3;
	imut[3][0] = -3;
	imut[3][1] = -3;
	imut[3][2] = -3;
	imut[3][3] = 2;
	imut[3][4] = 2;
	imut[3][5] = -3;
	imut[4][0] = -3;
	imut[4][1] = -3;
	imut[4][2] = -3;
	imut[4][3] = 2;
	imut[4][4] = 2;
	imut[4][5] = -3;
	imut[5][0] = -3;
	imut[5][1] = -3;
	imut[5][2] = -3;
	imut[5][3] = -3;
	imut[5][4] = -3;
	imut[5][5] = 0;
}

inline void CNWalign::load_blosum62_for_protein() {
	imut_seq = "ARNDCQEGHILKMFPSTWYVBZX";
	imut_seq_len = imut_seq.size();
	imut[0][0] = 4;
	imut[0][1] = -1;
	imut[0][2] = -2;
	imut[0][3] = -2;
	imut[0][4] = 0;
	imut[0][5] = -1;
	imut[0][6] = -1;
	imut[0][7] = 0;
	imut[0][8] = -2;
	imut[0][9] = -1;
	imut[0][10] = -1;
	imut[0][11] = -1;
	imut[0][12] = -1;
	imut[0][13] = -2;
	imut[0][14] = -1;
	imut[0][15] = 1;
	imut[0][16] = 0;
	imut[0][17] = -3;
	imut[0][18] = -2;
	imut[0][19] = 0;
	imut[0][20] = -2;
	imut[0][21] = -1;
	imut[0][22] = 0;
	imut[1][0] = -1;
	imut[1][1] = 5;
	imut[1][2] = 0;
	imut[1][3] = -2;
	imut[1][4] = -3;
	imut[1][5] = 1;
	imut[1][6] = 0;
	imut[1][7] = -2;
	imut[1][8] = 0;
	imut[1][9] = -3;
	imut[1][10] = -2;
	imut[1][11] = 2;
	imut[1][12] = -1;
	imut[1][13] = -3;
	imut[1][14] = -2;
	imut[1][15] = -1;
	imut[1][16] = -1;
	imut[1][17] = -3;
	imut[1][18] = -2;
	imut[1][19] = -3;
	imut[1][20] = -1;
	imut[1][21] = 0;
	imut[1][22] = -1;
	imut[2][0] = -2;
	imut[2][1] = 0;
	imut[2][2] = 6;
	imut[2][3] = 1;
	imut[2][4] = -3;
	imut[2][5] = 0;
	imut[2][6] = 0;
	imut[2][7] = 0;
	imut[2][8] = 1;
	imut[2][9] = -3;
	imut[2][10] = -3;
	imut[2][11] = 0;
	imut[2][12] = -2;
	imut[2][13] = -3;
	imut[2][14] = -2;
	imut[2][15] = 1;
	imut[2][16] = 0;
	imut[2][17] = -4;
	imut[2][18] = -2;
	imut[2][19] = -3;
	imut[2][20] = 3;
	imut[2][21] = 0;
	imut[2][22] = -1;
	imut[3][0] = -2;
	imut[3][1] = -2;
	imut[3][2] = 1;
	imut[3][3] = 6;
	imut[3][4] = -3;
	imut[3][5] = 0;
	imut[3][6] = 2;
	imut[3][7] = -1;
	imut[3][8] = -1;
	imut[3][9] = -3;
	imut[3][10] = -4;
	imut[3][11] = -1;
	imut[3][12] = -3;
	imut[3][13] = -3;
	imut[3][14] = -1;
	imut[3][15] = 0;
	imut[3][16] = -1;
	imut[3][17] = -4;
	imut[3][18] = -3;
	imut[3][19] = -3;
	imut[3][20] = 4;
	imut[3][21] = 1;
	imut[3][22] = -1;
	imut[4][0] = 0;
	imut[4][1] = -3;
	imut[4][2] = -3;
	imut[4][3] = -3;
	imut[4][4] = 9;
	imut[4][5] = -3;
	imut[4][6] = -4;
	imut[4][7] = -3;
	imut[4][8] = -3;
	imut[4][9] = -1;
	imut[4][10] = -1;
	imut[4][11] = -3;
	imut[4][12] = -1;
	imut[4][13] = -2;
	imut[4][14] = -3;
	imut[4][15] = -1;
	imut[4][16] = -1;
	imut[4][17] = -2;
	imut[4][18] = -2;
	imut[4][19] = -1;
	imut[4][20] = -3;
	imut[4][21] = -3;
	imut[4][22] = -2;
	imut[5][0] = -1;
	imut[5][1] = 1;
	imut[5][2] = 0;
	imut[5][3] = 0;
	imut[5][4] = -3;
	imut[5][5] = 5;
	imut[5][6] = 2;
	imut[5][7] = -2;
	imut[5][8] = 0;
	imut[5][9] = -3;
	imut[5][10] = -2;
	imut[5][11] = 1;
	imut[5][12] = 0;
	imut[5][13] = -3;
	imut[5][14] = -1;
	imut[5][15] = 0;
	imut[5][16] = -1;
	imut[5][17] = -2;
	imut[5][18] = -1;
	imut[5][19] = -2;
	imut[5][20] = 0;
	imut[5][21] = 3;
	imut[5][22] = -1;
	imut[6][0] = -1;
	imut[6][1] = 0;
	imut[6][2] = 0;
	imut[6][3] = 2;
	imut[6][4] = -4;
	imut[6][5] = 2;
	imut[6][6] = 5;
	imut[6][7] = -2;
	imut[6][8] = 0;
	imut[6][9] = -3;
	imut[6][10] = -3;
	imut[6][11] = 1;
	imut[6][12] = -2;
	imut[6][13] = -3;
	imut[6][14] = -1;
	imut[6][15] = 0;
	imut[6][16] = -1;
	imut[6][17] = -3;
	imut[6][18] = -2;
	imut[6][19] = -2;
	imut[6][20] = 1;
	imut[6][21] = 4;
	imut[6][22] = -1;
	imut[7][0] = 0;
	imut[7][1] = -2;
	imut[7][2] = 0;
	imut[7][3] = -1;
	imut[7][4] = -3;
	imut[7][5] = -2;
	imut[7][6] = -2;
	imut[7][7] = 6;
	imut[7][8] = -2;
	imut[7][9] = -4;
	imut[7][10] = -4;
	imut[7][11] = -2;
	imut[7][12] = -3;
	imut[7][13] = -3;
	imut[7][14] = -2;
	imut[7][15] = 0;
	imut[7][16] = -2;
	imut[7][17] = -2;
	imut[7][18] = -3;
	imut[7][19] = -3;
	imut[7][20] = -1;
	imut[7][21] = -2;
	imut[7][22] = -1;
	imut[8][0] = -2;
	imut[8][1] = 0;
	imut[8][2] = 1;
	imut[8][3] = -1;
	imut[8][4] = -3;
	imut[8][5] = 0;
	imut[8][6] = 0;
	imut[8][7] = -2;
	imut[8][8] = 8;
	imut[8][9] = -3;
	imut[8][10] = -3;
	imut[8][11] = -1;
	imut[8][12] = -2;
	imut[8][13] = -1;
	imut[8][14] = -2;
	imut[8][15] = -1;
	imut[8][16] = -2;
	imut[8][17] = -2;
	imut[8][18] = 2;
	imut[8][19] = -3;
	imut[8][20] = 0;
	imut[8][21] = 0;
	imut[8][22] = -1;
	imut[9][0] = -1;
	imut[9][1] = -3;
	imut[9][2] = -3;
	imut[9][3] = -3;
	imut[9][4] = -1;
	imut[9][5] = -3;
	imut[9][6] = -3;
	imut[9][7] = -4;
	imut[9][8] = -3;
	imut[9][9] = 4;
	imut[9][10] = 2;
	imut[9][11] = -3;
	imut[9][12] = 1;
	imut[9][13] = 0;
	imut[9][14] = -3;
	imut[9][15] = -2;
	imut[9][16] = -1;
	imut[9][17] = -3;
	imut[9][18] = -1;
	imut[9][19] = 3;
	imut[9][20] = -3;
	imut[9][21] = -3;
	imut[9][22] = -1;
	imut[10][0] = -1;
	imut[10][1] = -2;
	imut[10][2] = -3;
	imut[10][3] = -4;
	imut[10][4] = -1;
	imut[10][5] = -2;
	imut[10][6] = -3;
	imut[10][7] = -4;
	imut[10][8] = -3;
	imut[10][9] = 2;
	imut[10][10] = 4;
	imut[10][11] = -2;
	imut[10][12] = 2;
	imut[10][13] = 0;
	imut[10][14] = -3;
	imut[10][15] = -2;
	imut[10][16] = -1;
	imut[10][17] = -2;
	imut[10][18] = -1;
	imut[10][19] = 1;
	imut[10][20] = -4;
	imut[10][21] = -3;
	imut[10][22] = -1;
	imut[11][0] = -1;
	imut[11][1] = 2;
	imut[11][2] = 0;
	imut[11][3] = -1;
	imut[11][4] = -3;
	imut[11][5] = 1;
	imut[11][6] = 1;
	imut[11][7] = -2;
	imut[11][8] = -1;
	imut[11][9] = -3;
	imut[11][10] = -2;
	imut[11][11] = 5;
	imut[11][12] = -1;
	imut[11][13] = -3;
	imut[11][14] = -1;
	imut[11][15] = 0;
	imut[11][16] = -1;
	imut[11][17] = -3;
	imut[11][18] = -2;
	imut[11][19] = -2;
	imut[11][20] = 0;
	imut[11][21] = 1;
	imut[11][22] = -1;
	imut[12][0] = -1;
	imut[12][1] = -1;
	imut[12][2] = -2;
	imut[12][3] = -3;
	imut[12][4] = -1;
	imut[12][5] = 0;
	imut[12][6] = -2;
	imut[12][7] = -3;
	imut[12][8] = -2;
	imut[12][9] = 1;
	imut[12][10] = 2;
	imut[12][11] = -1;
	imut[12][12] = 5;
	imut[12][13] = 0;
	imut[12][14] = -2;
	imut[12][15] = -1;
	imut[12][16] = -1;
	imut[12][17] = -1;
	imut[12][18] = -1;
	imut[12][19] = 1;
	imut[12][20] = -3;
	imut[12][21] = -1;
	imut[12][22] = -1;
	imut[13][0] = -2;
	imut[13][1] = -3;
	imut[13][2] = -3;
	imut[13][3] = -3;
	imut[13][4] = -2;
	imut[13][5] = -3;
	imut[13][6] = -3;
	imut[13][7] = -3;
	imut[13][8] = -1;
	imut[13][9] = 0;
	imut[13][10] = 0;
	imut[13][11] = -3;
	imut[13][12] = 0;
	imut[13][13] = 6;
	imut[13][14] = -4;
	imut[13][15] = -2;
	imut[13][16] = -2;
	imut[13][17] = 1;
	imut[13][18] = 3;
	imut[13][19] = -1;
	imut[13][20] = -3;
	imut[13][21] = -3;
	imut[13][22] = -1;
	imut[14][0] = -1;
	imut[14][1] = -2;
	imut[14][2] = -2;
	imut[14][3] = -1;
	imut[14][4] = -3;
	imut[14][5] = -1;
	imut[14][6] = -1;
	imut[14][7] = -2;
	imut[14][8] = -2;
	imut[14][9] = -3;
	imut[14][10] = -3;
	imut[14][11] = -1;
	imut[14][12] = -2;
	imut[14][13] = -4;
	imut[14][14] = 7;
	imut[14][15] = -1;
	imut[14][16] = -1;
	imut[14][17] = -4;
	imut[14][18] = -3;
	imut[14][19] = -2;
	imut[14][20] = -2;
	imut[14][21] = -1;
	imut[14][22] = -2;
	imut[15][0] = 1;
	imut[15][1] = -1;
	imut[15][2] = 1;
	imut[15][3] = 0;
	imut[15][4] = -1;
	imut[15][5] = 0;
	imut[15][6] = 0;
	imut[15][7] = 0;
	imut[15][8] = -1;
	imut[15][9] = -2;
	imut[15][10] = -2;
	imut[15][11] = 0;
	imut[15][12] = -1;
	imut[15][13] = -2;
	imut[15][14] = -1;
	imut[15][15] = 4;
	imut[15][16] = 1;
	imut[15][17] = -3;
	imut[15][18] = -2;
	imut[15][19] = -2;
	imut[15][20] = 0;
	imut[15][21] = 0;
	imut[15][22] = 0;
	imut[16][0] = 0;
	imut[16][1] = -1;
	imut[16][2] = 0;
	imut[16][3] = -1;
	imut[16][4] = -1;
	imut[16][5] = -1;
	imut[16][6] = -1;
	imut[16][7] = -2;
	imut[16][8] = -2;
	imut[16][9] = -1;
	imut[16][10] = -1;
	imut[16][11] = -1;
	imut[16][12] = -1;
	imut[16][13] = -2;
	imut[16][14] = -1;
	imut[16][15] = 1;
	imut[16][16] = 5;
	imut[16][17] = -2;
	imut[16][18] = -2;
	imut[16][19] = 0;
	imut[16][20] = -1;
	imut[16][21] = -1;
	imut[16][22] = 0;
	imut[17][0] = -3;
	imut[17][1] = -3;
	imut[17][2] = -4;
	imut[17][3] = -4;
	imut[17][4] = -2;
	imut[17][5] = -2;
	imut[17][6] = -3;
	imut[17][7] = -2;
	imut[17][8] = -2;
	imut[17][9] = -3;
	imut[17][10] = -2;
	imut[17][11] = -3;
	imut[17][12] = -1;
	imut[17][13] = 1;
	imut[17][14] = -4;
	imut[17][15] = -3;
	imut[17][16] = -2;
	imut[17][17] = 11;
	imut[17][18] = 2;
	imut[17][19] = -3;
	imut[17][20] = -4;
	imut[17][21] = -3;
	imut[17][22] = -2;
	imut[18][0] = -2;
	imut[18][1] = -2;
	imut[18][2] = -2;
	imut[18][3] = -3;
	imut[18][4] = -2;
	imut[18][5] = -1;
	imut[18][6] = -2;
	imut[18][7] = -3;
	imut[18][8] = 2;
	imut[18][9] = -1;
	imut[18][10] = -1;
	imut[18][11] = -2;
	imut[18][12] = -1;
	imut[18][13] = 3;
	imut[18][14] = -3;
	imut[18][15] = -2;
	imut[18][16] = -2;
	imut[18][17] = 2;
	imut[18][18] = 7;
	imut[18][19] = -1;
	imut[18][20] = -3;
	imut[18][21] = -2;
	imut[18][22] = -1;
	imut[19][0] = 0;
	imut[19][1] = -3;
	imut[19][2] = -3;
	imut[19][3] = -3;
	imut[19][4] = -1;
	imut[19][5] = -2;
	imut[19][6] = -2;
	imut[19][7] = -3;
	imut[19][8] = -3;
	imut[19][9] = 3;
	imut[19][10] = 1;
	imut[19][11] = -2;
	imut[19][12] = 1;
	imut[19][13] = -1;
	imut[19][14] = -2;
	imut[19][15] = -2;
	imut[19][16] = 0;
	imut[19][17] = -3;
	imut[19][18] = -1;
	imut[19][19] = 4;
	imut[19][20] = -3;
	imut[19][21] = -2;
	imut[19][22] = -1;
	imut[20][0] = -2;
	imut[20][1] = -1;
	imut[20][2] = 3;
	imut[20][3] = 4;
	imut[20][4] = -3;
	imut[20][5] = 0;
	imut[20][6] = 1;
	imut[20][7] = -1;
	imut[20][8] = 0;
	imut[20][9] = -3;
	imut[20][10] = -4;
	imut[20][11] = 0;
	imut[20][12] = -3;
	imut[20][13] = -3;
	imut[20][14] = -2;
	imut[20][15] = 0;
	imut[20][16] = -1;
	imut[20][17] = -4;
	imut[20][18] = -3;
	imut[20][19] = -3;
	imut[20][20] = 4;
	imut[20][21] = 1;
	imut[20][22] = -1;
	imut[21][0] = -1;
	imut[21][1] = 0;
	imut[21][2] = 0;
	imut[21][3] = 1;
	imut[21][4] = -3;
	imut[21][5] = 3;
	imut[21][6] = 4;
	imut[21][7] = -2;
	imut[21][8] = 0;
	imut[21][9] = -3;
	imut[21][10] = -3;
	imut[21][11] = 1;
	imut[21][12] = -1;
	imut[21][13] = -3;
	imut[21][14] = -1;
	imut[21][15] = 0;
	imut[21][16] = -1;
	imut[21][17] = -3;
	imut[21][18] = -2;
	imut[21][19] = -2;
	imut[21][20] = 1;
	imut[21][21] = 4;
	imut[21][22] = -1;
	imut[22][0] = 0;
	imut[22][1] = -1;
	imut[22][2] = -1;
	imut[22][3] = -1;
	imut[22][4] = -2;
	imut[22][5] = -1;
	imut[22][6] = -1;
	imut[22][7] = -1;
	imut[22][8] = -1;
	imut[22][9] = -1;
	imut[22][10] = -1;
	imut[22][11] = -1;
	imut[22][12] = -1;
	imut[22][13] = -1;
	imut[22][14] = -2;
	imut[22][15] = 0;
	imut[22][16] = 0;
	imut[22][17] = -2;
	imut[22][18] = -1;
	imut[22][19] = -1;
	imut[22][20] = -1;
	imut[22][21] = -1;
	imut[22][22] = -1;
}

inline string CBaseFunc::charVec2String(const vector<char>& vec){
	int size = vec.size();
	char* parr = new char[size+1];
	for (int i = 0; i < size; i++)
		parr[i] = vec[i];
	parr[size] = '\0';
	
	string ans(parr);
	delete[] parr; parr = NULL;
	return ans;
}

inline double CBaseFunc::fabs(const double& det){
	return det>0.?det:-det;
}

inline CBaseFunc::CBaseFunc()
{

}

inline CBaseFunc::~CBaseFunc()
{

}


/*******************************************************
 * @param r : the row number of the 2D matrix
 * @param c : the column number of the 2D matrix
 * @return  : the pointer of the new 2D matrix
 *******************************************************/
inline double** CBaseFunc::new2Darr(const int& r, const int& c){
	double** ans = new double*[r];
	for (int i = 0; i < r; i++){
		ans[i] = new double[c];
		for (int j = 0; j < c; j++){
			ans[i][j] = 0.0;
		}
	}

	return ans;
}

/*******************************************************
 * @param r : the row number of the 2D matrix
 * @function: release the memory of the 2D matrix
 *******************************************************/
inline void CBaseFunc::delete2Darr(double** pMtx, const int& r){
	for (int i = 0; i < r; i++){
		delete[] pMtx[i];
	}
	delete[] pMtx;
	pMtx = NULL;
}

inline bool CBaseFunc::isInit2Darr(double** mtx, const int& r, const int& c){
	int i, j;
	for (i = 0; i < r; i++){
		for (j = 0; j < c; j++){
			if (fabs(mtx[i][j]) > g_eps){
				return false;
			}
		}
	}
	
	return true;
}

inline int** CBaseFunc::new2DIntArr(const int& row, const int& col){
	int **ans=new int*[row];
	for(int i=0;i<row;i++){
		ans[i]=new int[col];
		for(int j=0; j<col; j++)
			ans[i][j] = 0;
	}
	
	return ans;
}

inline void CBaseFunc::delete2DIntArr(const int& n, int ** Arr){
	for(int i = 0; i < n; i++){
		delete [] Arr[i];
	}
	delete[] Arr;
	Arr = NULL;
}

inline bool** CBaseFunc::new2DBoolArr(const int& row, const int& col){
	bool **ans=new bool*[row];
	for(int i=0;i<row;i++){
		ans[i]=new bool[col];
		for(int j=0; j<col; j++)
			ans[i][j] = false;
	}
	
	return ans;
}

inline bool* CBaseFunc::new1Dbool(const int& row){
	bool* ans = new bool[row];
	for (int i = 0; i < row; i++)
		ans[i] = false;
	return ans;
}

inline void CBaseFunc::delete2DBoolArr(const int& n, bool ** Arr){
	for(int i = 0; i < n; i++){
		delete [] Arr[i];
	}
	delete[] Arr;
	Arr = NULL;
}

/*********************************************************************
 * @param acid : the 3 words abbreviation of amino acid
 * @function   : transform the abbreviation of amino acid from 3 words
 *               to 1 word.
 *               e.g., ALA -> A
 *********************************************************************/
inline char CBaseFunc::aminoAcidAbbr3WordsTo1(const string& acid) {
	
	if (0 == acid.compare("ALA"))	return 'A';
	if (0 == acid.compare("CYS"))	return 'C';
	if (0 == acid.compare("ASP"))	return 'D';
	if (0 == acid.compare("GLU"))	return 'E';
	if (0 == acid.compare("PHE"))	return 'F';
	if (0 == acid.compare("GLY"))	return 'G';
	if (0 == acid.compare("HIS"))	return 'H';
	if (0 == acid.compare("ILE"))	return 'I';
	if (0 == acid.compare("LYS"))	return 'K';
	if (0 == acid.compare("LEU"))	return 'L';
	if (0 == acid.compare("MET"))	return 'M';
	if (0 == acid.compare("ASN"))	return 'N';
	if (0 == acid.compare("PRO"))	return 'P';
	if (0 == acid.compare("GLN"))	return 'Q';
	if (0 == acid.compare("ARG"))	return 'R';
	if (0 == acid.compare("SER"))	return 'S';
	if (0 == acid.compare("THR"))	return 'T';
	if (0 == acid.compare("VAL"))	return 'V';
	if (0 == acid.compare("TRP"))	return 'W';
	if (0 == acid.compare("TYR"))	return 'Y';
	
	return 'X';
}

inline char CBaseFunc::nucletideAbbr3WordsTo1(const string& nuc) {
	
	if (0 == nuc.compare("DA"))	return 'a';
	if (0 == nuc.compare("DC"))	return 'c';
	if (0 == nuc.compare("DG"))	return 'g';
	if (0 == nuc.compare("DT"))	return 't';
	if (0 == nuc.compare("DU"))	return 'u';
	if (0 == nuc.compare("A"))	return 'a';
	if (0 == nuc.compare("C"))	return 'c';
	if (0 == nuc.compare("G"))	return 'g';
	if (0 == nuc.compare("T"))	return 't';
	if (0 == nuc.compare("U"))	return 'u';
	if (0 == nuc.compare("RA"))	return 'a';
	if (0 == nuc.compare("RC"))	return 'c';
	if (0 == nuc.compare("RG"))	return 'g';
	if (0 == nuc.compare("RT"))	return 't';
	if (0 == nuc.compare("RU"))	return 'u';

	return 'x';
}

/*********************************************************************
 * @param acid : the 1 word abbreviation of amino acid
 * @function   : transform the abbreviation of amino acid from 1 word
 *               to 3 words.
 *               e.g., A -> ALA
 *********************************************************************/
inline string CBaseFunc::aminoAcidAbbr1WordTo3(const char& acid) {
	
	if ('A' == acid)	return "ALA";
	if ('C' == acid)	return "CYS";
	if ('D' == acid)	return "ASP";
	if ('E' == acid)	return "GLU";
	if ('F' == acid)	return "PHE";
	if ('G' == acid)	return "GLY";
	if ('H' == acid)	return "HIS";
	if ('I' == acid)	return "ILE";
	if ('K' == acid)	return "LYS";
	if ('L' == acid)	return "LEU";
	if ('M' == acid)	return "MET";
	if ('N' == acid)	return "ASN";
	if ('P' == acid)	return "PRO";
	if ('Q' == acid)	return "GLN";
	if ('R' == acid)	return "ARG";
	if ('S' == acid)	return "SER";
	if ('T' == acid)	return "THR";
	if ('V' == acid)	return "VAL";
	if ('W' == acid)	return "TRP";
	if ('Y' == acid)	return "TYR";
	
	return "JUN";
}

/*****************************************************
 * @comments: remove the head and tail spaces
 *****************************************************/
inline string CBaseFunc::stringTrim(const string& str) {
	string ans;
	int start_index = 0;		// inclusive
	int end_index = str.size(); // exclusive
	int i = 0;
	// remove the previous space
	for (start_index = 0; start_index < str.size(); start_index++){
		if (' ' != str[start_index] && '\t' != str[start_index] && '\r' != str[start_index] && '\n' != str[start_index]){
			break;
		}
	}

	for (end_index = str.size(); end_index > start_index; end_index--){
		if (' ' != str[end_index-1] && '\t' != str[end_index-1] && '\r' != str[end_index-1] && '\n' != str[end_index-1]) {
			break;
		}
	}

	for (i = start_index; i < end_index; i++){
		ans += str[i];
	}

	return ans;
}

inline vector<string> CBaseFunc::stringSplit(const string& str, const char& spit){
	int i, size = str.size();
	vector<string> ans;
	vector<char> one;
	for (i = 0; i < size; i++){
		if (spit == str[i]){
			if (one.size() != 0){
				string sub = charVec2String(one);
				ans.push_back(sub);
				vector<char>().swap(one);
			}
		}else{
			one.push_back(str[i]);
		}
	}
	if (one.size() != 0){
		string sub = charVec2String(one);
		ans.push_back(sub);
	}
	
	return ans;
}

inline vector<string> CBaseFunc::stringSplit(const string& str, const char& spit1, const char& spit2){
	int i, size = str.size();
	vector<string> ans;
	vector<char> one;
	for (i = 0; i < size; i++){
		if (spit1 == str[i] || spit2 == str[i]){
			if (one.size() != 0){
				string sub = charVec2String(one);
				ans.push_back(sub);
				vector<char>().swap(one);
			}
		}else{
			one.push_back(str[i]);
		}
	}
	if (one.size() != 0){
		string sub = charVec2String(one);
		ans.push_back(sub);
	}
	
	return ans;
}

inline double CBaseFunc::d0_of_tmscore(const int& Lnorm){
	double d0 = 0.;
		
	if (Lnorm > 15) {
		d0 = 1.24 * pow((Lnorm - 15), (1.0 / 3.0)) - 1.8;
	} else {
		d0 = 0.5;
	}
	if (d0 < 0.5) {
		d0 = 0.5;
	}
	
	return d0;
}

inline double CBaseFunc::d0_of_lsscore(const int& Lnorm){
	double d0 = 0.;
	if (Lnorm < 9){
		d0 = 0.1;
	}else{
		d0 = 0.55 * pow(Lnorm - 9, 1.0/3.0) + 0.15;
	}
	
	if (d0 < 0.1)
		d0 = 0.1;
	if (d0 > 4.0)
		d0 = 4.0;
		
	return d0;
} 

inline double CBaseFunc::d0_of_tmscore_c3prime(const int& Lnorm){
	double d0 = 0.;
	if(Lnorm<=11) d0=0.3;
    else if(Lnorm>11&&Lnorm<=15) d0=0.4;
    else if(Lnorm>15&&Lnorm<=19) d0=0.5;
    else if(Lnorm>19&&Lnorm<=23) d0=0.6;
    else if(Lnorm>23&&Lnorm<30)  d0=0.7;
    else d0=(0.6*pow((Lnorm*1.0-0.5), 1.0/2)-2.5);
    
    return d0;
}

inline double CBaseFunc::__2merExchangedSearch(double** scoMtx, const int& rownum, const int& colnum, int* alivec, int* transpose_alivec, const double& prev_score, const int& iter_max){
	int row_ind, col_ind, iter, exchange_row_ind, exchange_col_ind, new_col_ind, rela_another_row_ind;
	double maxScoreChange, valChange, prev_sco;
	
	double score = prev_score;
	for (iter = 0; iter < iter_max; iter++){
		prev_sco = score;
		if (iter % 2 != 0){
			// direction: 0 -> rownum
			for (row_ind = 0; row_ind < rownum; row_ind++){
				
				col_ind = alivec[row_ind];
				maxScoreChange = 0.;
				exchange_row_ind = -1;
				exchange_col_ind = -1;
				for (new_col_ind = 0; new_col_ind < colnum; new_col_ind++){
					if (new_col_ind == col_ind) continue;
					rela_another_row_ind = transpose_alivec[new_col_ind];
					
					// calculate the changed value 
					valChange = scoMtx[row_ind][new_col_ind];
					if (-1 != col_ind)
						valChange = scoMtx[row_ind][new_col_ind] - scoMtx[row_ind][col_ind];
					
					if (-1 != rela_another_row_ind){
						if (-1 != col_ind)
							valChange += scoMtx[rela_another_row_ind][col_ind] - scoMtx[rela_another_row_ind][new_col_ind];
						else valChange -= scoMtx[rela_another_row_ind][new_col_ind];
					}
					
					if (valChange > maxScoreChange){
						maxScoreChange = valChange;
						exchange_row_ind = rela_another_row_ind;
						exchange_col_ind = new_col_ind;
					}
				}
			
				if (maxScoreChange > 1e-9){
					score += maxScoreChange;
					alivec[row_ind] = exchange_col_ind;
					if (-1 != exchange_row_ind)
						alivec[exchange_row_ind] = col_ind;
						
					if (-1 != exchange_col_ind)
						transpose_alivec[exchange_col_ind] = row_ind;
					if (-1 != col_ind)
						transpose_alivec[col_ind] = exchange_row_ind;
				}
			}//for (row_ind
			
		}else{
			// direction: rownum -> 0 
			for (row_ind = rownum-1; row_ind >= 0; row_ind--){
				col_ind = alivec[row_ind];
				
				maxScoreChange = 0.;
				exchange_row_ind = -1;
				exchange_col_ind = -1;
				for (new_col_ind = 0; new_col_ind < colnum; new_col_ind++){
					if (new_col_ind == col_ind) continue;
					rela_another_row_ind = transpose_alivec[new_col_ind];
					
					// calculate the changed value 
					valChange = scoMtx[row_ind][new_col_ind];
					if (-1 != col_ind)
						valChange = scoMtx[row_ind][new_col_ind] - scoMtx[row_ind][col_ind];
					
					if (-1 != rela_another_row_ind){
						if (-1 != col_ind)
							valChange += scoMtx[rela_another_row_ind][col_ind] - scoMtx[rela_another_row_ind][new_col_ind];
						else valChange -= scoMtx[rela_another_row_ind][new_col_ind];
					}
					
					if (valChange > maxScoreChange){
						maxScoreChange = valChange;
						exchange_row_ind = rela_another_row_ind;
						exchange_col_ind = new_col_ind;
					}
				}
			
				if (maxScoreChange > 1e-9){
					score += maxScoreChange;
					alivec[row_ind] = exchange_col_ind;
					if (-1 != exchange_row_ind)
						alivec[exchange_row_ind] = col_ind;
						
					if (-1 != exchange_col_ind)
						transpose_alivec[exchange_col_ind] = row_ind;
					if (-1 != col_ind)
						transpose_alivec[col_ind] = exchange_row_ind;
				}
			}//for (row_ind
		}
		
		if (fabs(score - prev_sco) < 1e-9) break; 
	}
	
	return score;
}

inline double CBaseFunc::__3merExchangedSearch(double** scoMtx, const int& rownum, const int& colnum, int* alivec, int* transpose_alivec, const double& prev_score, const int& iter_max){
	int iter, row_ind, col_ind, jth_exchange_row_ind, jth_exchange_col_ind, kth_exchange_row_ind, kth_exchange_col_ind,
			 jth_new_col_ind, jth_rela_another_row_ind, kth_new_col_ind, kth_rela_another_row_ind;
	double maxScoreChange, valChange, prev_sco;
	
	double score = prev_score;
	for (iter = 0; iter < iter_max; iter++){
		prev_sco = score;
		for (row_ind = 0; row_ind < rownum; row_ind++){
			col_ind = alivec[row_ind];
			
			maxScoreChange = 0.;
			jth_exchange_row_ind = -1;
			jth_exchange_col_ind = -1;
			kth_exchange_row_ind = -1;
			kth_exchange_col_ind = -1;
			for (jth_new_col_ind = 0; jth_new_col_ind < colnum; jth_new_col_ind++){ // j does not mean the column index
				if (jth_new_col_ind == col_ind) continue;
				jth_rela_another_row_ind = transpose_alivec[jth_new_col_ind];
				if (-1 == jth_rela_another_row_ind) continue;
					
				for (kth_new_col_ind = 0; kth_new_col_ind < colnum; kth_new_col_ind++){// k does not mean the column index
					if (kth_new_col_ind == col_ind || kth_new_col_ind == jth_new_col_ind) continue;
					kth_rela_another_row_ind = transpose_alivec[kth_new_col_ind];
							
					// calculate the changed value
					valChange = scoMtx[row_ind][jth_new_col_ind];
					if (-1 != col_ind)
						valChange = scoMtx[row_ind][jth_new_col_ind] - scoMtx[row_ind][col_ind];
								
					valChange += scoMtx[jth_rela_another_row_ind][kth_new_col_ind] - scoMtx[jth_rela_another_row_ind][jth_new_col_ind];
							
					if (-1 != kth_rela_another_row_ind){
						if (-1 != col_ind)
								valChange += scoMtx[kth_rela_another_row_ind][col_ind] - scoMtx[kth_rela_another_row_ind][kth_new_col_ind];
						else valChange -= scoMtx[kth_rela_another_row_ind][kth_new_col_ind];
					}
							
					if (valChange > maxScoreChange){
						maxScoreChange = valChange;
						jth_exchange_row_ind = jth_rela_another_row_ind;
						jth_exchange_col_ind = jth_new_col_ind;
						kth_exchange_row_ind = kth_rela_another_row_ind;
						kth_exchange_col_ind = kth_new_col_ind;
					}
				}// for_k
			}//for_j
			
			if (maxScoreChange > 1e-9){
				score += maxScoreChange;
				
				if (scoMtx[row_ind][jth_exchange_col_ind] > 1e-9) {
					alivec[row_ind] = jth_exchange_col_ind;
					transpose_alivec[jth_exchange_col_ind] = row_ind;
				}else {
					alivec[row_ind] = -1;
					transpose_alivec[jth_exchange_col_ind] = -1;
				}
				
				if (scoMtx[jth_exchange_row_ind][kth_exchange_col_ind] > 1e-9) {
					alivec[jth_exchange_row_ind] = kth_exchange_col_ind;
					transpose_alivec[kth_exchange_col_ind] = jth_exchange_row_ind;
				}else {
					alivec[jth_exchange_row_ind] = -1;
					transpose_alivec[kth_exchange_col_ind] = -1;
				}
				
				if (-1 != kth_exchange_row_ind)
					if (-1 == col_ind || scoMtx[kth_exchange_row_ind][col_ind] <= 1e-9)
						alivec[kth_exchange_row_ind] = -1;
					else alivec[kth_exchange_row_ind] = col_ind;
						
				if (-1 != col_ind)
					if (-1 == kth_exchange_row_ind || scoMtx[kth_exchange_row_ind][col_ind] <= 1e-9)
						transpose_alivec[col_ind] = -1;
					else transpose_alivec[col_ind] = kth_exchange_row_ind;
				
//					System.out.println("=============================================================");
//					double sum = 0;
//					for (int i = 0; i < alivec.length; i++) {
//						if (-1 != alivec[i]) {
//							System.out.println(i + " : " + alivec[i] + " : " + scoMtx[i][alivec[i]]);
//							sum += scoMtx[i][alivec[i]];
//						}
//					}
//					System.out.println("sum = " + sum + ", score = " + score + ", " + maxScoreChange);
			}
		}//for (row_ind
		
		if (fabs(score - prev_sco) < 1e-9) break; 
	}
	
	return score;
}

// Please make sure that rownum > colnum, it will save several seconds O(n^3)
inline double CBaseFunc::greedySearch(double** scoMtx, const int& rownum, const int& colnum, int* alivec, int* transpose_alivec){
	int i, corrRowInd, corrColInd, colInd, rr, cc;
	double maxValue, value;
	
	for (i = 0; i < rownum; i++) alivec[i] = -1;
	for (i = 0; i < colnum; i++) transpose_alivec[i] = -1;
	
	int minnum = rownum > colnum ? colnum : rownum;
	double score = 0.0;
	for (i = 0; i < minnum; i++){	// i means the ith alignment, not row index
		maxValue = 0.0;
		corrRowInd = -1;
		corrColInd = -1;
		for (rr = 0; rr < rownum; rr++){	// rr means the rr-th row in rowSortedInds
			if (-1 == alivec[rr]){
				colInd = -1;
				value = 0.;
				for (cc = 0; cc < colnum; cc++){
					if (-1 == transpose_alivec[cc] && value < scoMtx[rr][cc]){
						colInd = cc;
						value = scoMtx[rr][cc];
					}
				}
					
				if (maxValue < value){
					maxValue = value;
					corrRowInd = rr;
					corrColInd = colInd;
				}
			}
		}
		
		score += maxValue;
		if (-1 != corrRowInd)
			alivec[corrRowInd] = corrColInd;
		if (-1 != corrColInd)
			transpose_alivec[corrColInd] = corrRowInd;
	}
	
	return score;
}

inline double CBaseFunc::cal_rot_tran_from_query_to_templ__(
	const vector<double*>& query,
	const vector<double*>& templ,
	double** out_u,
	double* out_t,
	const double& user_d0,
	const bool& fast)
{
	int i, j, k;
	int nmax, nseq, n_ali;
	n_ali = nseq = query.size();
	nmax = nseq + 1;
	
	double d, d2, d0, d02;
	
	int n_cut; // ![1,n_ali],align residues for the score
	double score;
	double score_max;
	double d0_search;
	
	int* L_ini = new int[nmax];
	int* i_ali = new int[nmax];
	double** u = new2Darr(4, 4);
	double* t = new double[4];
	double* iq = new double[nmax];
	double* xa = new double[nmax];
	double* ya = new double[nmax];
	double* za = new double[nmax];
	double* xt = new double[nmax];
	double* yt = new double[nmax];
	double* zt = new double[nmax];
	double* xb = new double[nmax];
	double* yb = new double[nmax];
	double* zb = new double[nmax];
	double** r_1 = new2Darr(4, nmax);
	double** r_2 = new2Darr(4, nmax);
	int* k_ali = new int[nmax];
	int* k_ali0 = new int[nmax];
	
//		seq1A = "*" + p1.getSeq();
	for (int kk = 0; kk < nseq; kk++) {
		xa[kk+1] = query[kk][0];
		ya[kk+1] = query[kk][1];
		za[kk+1] = query[kk][2];
	}
	
//		seq1B = "*" + p2.getSeq();
	for (int kk = 0; kk < nseq; kk++) {
		xb[kk+1] = templ[kk][0];
		yb[kk+1] = templ[kk][1];
		zb[kk+1] = templ[kk][2];
	}
	
	int ka0 = 0;
	d0 = user_d0;
	d02 = d0*d0;
	
	// *** d0_search ----->
	d0_search = d0;
	if (d0_search > 8)
		d0_search = 8;
	if (d0_search < 4.5)
		d0_search = 4.5;
		
	// *** iterative parameters ----->
	//int n_it = 20; // !maximum number of iterations
	int n_it = g_maxmum_number_of_iterations;
	
	//int n_init_max = 6; // !maximum number of L_init
	int n_init_max = g_maxmum_number_of_L_init;
	
	int n_init = 0;
	int L_ini_min = 4;
	if (n_ali < 4)
		L_ini_min = n_ali;
	
	bool flag1 = false;
	for (i = 1; i <= (n_init_max - 1); i++) {
		n_init = n_init + 1;
		L_ini[n_init] = n_ali / (int) pow(2, (n_init - 1));
		if (L_ini[n_init] <= L_ini_min) {
			L_ini[n_init] = L_ini_min;
			flag1 = true;
			break;
		}
	}
	if (flag1 == false) {
		n_init = n_init + 1;
		L_ini[n_init] = L_ini_min;
	}

	score_max = 0; // !TM-score
	int LL, ka, i_init, L_init, iL_max, iL;
	for (i_init = 1; i_init <= n_init; i_init++) { // 333
		L_init = L_ini[i_init];
		iL_max = n_ali - L_init + 1;
		for (iL = 1; iL <= iL_max; fast? iL+=2 : iL++) { // 300 !on aligned residues,
			LL = 0;
			ka = 0;
			for (i = 1; i <= L_init; i++) {
				k = iL + i - 1; // ![1,n_ali] common
				r_1[1][i] = xa[k];
				r_1[2][i] = ya[k];
				r_1[3][i] = za[k];
				r_2[1][i] = xb[k];
				r_2[2][i] = yb[k];
				r_2[3][i] = zb[k];
				ka = ka + 1;
				k_ali[ka] = k;
				LL = LL + 1;
			}
			u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
			for (j = 1; j <= nseq; j++) {
				xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
				yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
				zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
			}
			d = d0_search - 1;
			score_fun(xt, yt, zt, xb, yb, zb, d*d, n_cut, n_ali, i_ali, d02, nseq, score); 
			// iteration
			if (score_max < score) {
				score_max = score;
				ka0 = ka;
				for (i = 1; i <= ka0; i++) {
					k_ali0[i] = k_ali[i];
				}
				
				for (i = 1; i < 4; i++) {
					for (j = 1; j < 4; j++) {
						out_u[i-1][j-1] = u[i][j];
					}
					out_t[i-1] = t[i];
				}
			}
			
			// *** iteration for extending
			// ---------------------------------->
			d = d0_search + 1;
			for (int it = 1; it <= n_it; it++) {
				LL = 0;
				ka = 0;
				for (i = 1; i <= n_cut; i++) {
					int m = i_ali[i]; // ![1,n_ali]
					r_1[1][i] = xa[m];
					r_1[2][i] = ya[m];
					r_1[3][i] = za[m];
					r_2[1][i] = xb[m];
					r_2[2][i] = yb[m];
					r_2[3][i] = zb[m];
					ka = ka + 1;
					k_ali[ka] = m;
					LL = LL + 1;
				}
				u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
				for (j = 1; j <= nseq; j++) {
					xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
					yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
					zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
				}
				score_fun(xt, yt, zt, xb, yb, zb, d*d, n_cut, n_ali, i_ali, d02, nseq, score);
				if (score_max < score) {
					score_max = score;
					ka0 = ka;
					for (i = 1; i <= ka; i++) {
						k_ali0[i] = k_ali[i];
					}
					
					for (i = 1; i < 4; i++) {
						for (j = 1; j < 4; j++) {
							out_u[i-1][j-1] = u[i][j];
						}
						out_t[i-1] = t[i];
					}
				}
				if (it == n_it) {
					break;
				}
				if (n_cut == ka) { // then
					int neq = 0;
					for (i = 1; i <= n_cut; i++) {
						if (i_ali[i] == k_ali[i]) {
							neq = neq + 1;
						}
					}
					if (n_cut == neq) {
						break;
					}
				}
			} 
		}
	}
	
	delete[] L_ini;
	delete[] i_ali;
	delete2Darr(u, 4);
	delete[] t;
	delete[] iq;
	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] xt;
	delete[] yt;
	delete[] zt;
	delete[] xb;
	delete[] yb;
	delete[] zb;
	delete2Darr(r_1, 4);
	delete2Darr(r_2, 4);
	delete[] k_ali;
	delete[] k_ali0;
	
	return score_max;
}

inline double CBaseFunc::cal_rot_tran_from_query_to_templ__II(
	const vector<double*>& query,
	const vector<double*>& templ,
	double** out_u,
	double* out_t,
	const double& user_d0,
	const bool& fast)
{
	int i, j, k;
	int nmax, nseq, n_ali;
	n_ali = nseq = query.size();
	nmax = nseq + 1;
	
	double d, d2, d0, d02;
	
	int n_cut; // ![1,n_ali],align residues for the score
	double score;
	double score_max;
	double d0_search;
	
	int* L_ini = new int[nmax];
	int* i_ali = new int[nmax];
	double** u = new2Darr(4, 4);
	double* t = new double[4];
	double* iq = new double[nmax];
	double* xa = new double[nmax];
	double* ya = new double[nmax];
	double* za = new double[nmax];
	double* xt = new double[nmax];
	double* yt = new double[nmax];
	double* zt = new double[nmax];
	double* xb = new double[nmax];
	double* yb = new double[nmax];
	double* zb = new double[nmax];
	double** r_1 = new2Darr(4, nmax);
	double** r_2 = new2Darr(4, nmax);
	int* k_ali = new int[nmax];
	int* k_ali0 = new int[nmax];
	
//		seq1A = "*" + p1.getSeq();
	for (int kk = 0; kk < nseq; kk++) {
		xa[kk+1] = query[kk][0];
		ya[kk+1] = query[kk][1];
		za[kk+1] = query[kk][2];
	}
	
//		seq1B = "*" + p2.getSeq();
	for (int kk = 0; kk < nseq; kk++) {
		xb[kk+1] = templ[kk][0];
		yb[kk+1] = templ[kk][1];
		zb[kk+1] = templ[kk][2];
	}
	
	int ka0 = 0;
	d0 = user_d0;
	d02 = d0*d0;
	
	// *** d0_search ----->
	d0_search = d0;
	if (d0_search > 8)
		d0_search = 8;
	if (d0_search < 4.5)
		d0_search = 4.5;
		
	// *** iterative parameters ----->
	//int n_it = 20; // !maximum number of iterations
	int n_it = g_maxmum_number_of_iterations;
	
	//int n_init_max = 6; // !maximum number of L_init
	int n_init_max = g_maxmum_number_of_L_init;
	
	int n_init = 0;
	int L_ini_min = 4;
	if (n_ali < 4)
		L_ini_min = n_ali;
	
	bool flag1 = false;
	for (i = 1; i <= (n_init_max - 1); i++) {
		n_init = n_init + 1;
		L_ini[n_init] = n_ali / (int) pow(2, (n_init - 1));
		if (L_ini[n_init] <= L_ini_min) {
			L_ini[n_init] = L_ini_min;
			flag1 = true;
			break;
		}
	}
	if (flag1 == false) {
		n_init = n_init + 1;
		L_ini[n_init] = L_ini_min;
	}

	score_max = 0; // !TM-score
	int LL, ka, i_init, L_init, iL_max, iL;
	for (i_init = 1; i_init <= n_init; i_init++) { // 333
		L_init = L_ini[i_init];
		iL_max = n_ali - L_init + 1;
		for (iL = 1; iL <= iL_max; fast? iL+=40 : iL++) { // 300 !on aligned residues,
			LL = 0;
			ka = 0;
			for (i = 1; i <= L_init; i++) {
				k = iL + i - 1; // ![1,n_ali] common
				r_1[1][i] = xa[k];
				r_1[2][i] = ya[k];
				r_1[3][i] = za[k];
				r_2[1][i] = xb[k];
				r_2[2][i] = yb[k];
				r_2[3][i] = zb[k];
				ka = ka + 1;
				k_ali[ka] = k;
				LL = LL + 1;
			}
			u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
			for (j = 1; j <= nseq; j++) {
				xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
				yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
				zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
			}
			d = d0_search - 1;
			score_fun(xt, yt, zt, xb, yb, zb, d*d, n_cut, n_ali, i_ali, d02, nseq, score); 
			// iteration
			if (score_max < score) {
				score_max = score;
				ka0 = ka;
				for (i = 1; i <= ka0; i++) {
					k_ali0[i] = k_ali[i];
				}
				
				for (i = 1; i < 4; i++) {
					for (j = 1; j < 4; j++) {
						out_u[i-1][j-1] = u[i][j];
					}
					out_t[i-1] = t[i];
				}
			}
			
			// *** iteration for extending
			// ---------------------------------->
			d = d0_search + 1;
			for (int it = 1; it <= n_it; it++) {
				LL = 0;
				ka = 0;
				for (i = 1; i <= n_cut; i++) {
					int m = i_ali[i]; // ![1,n_ali]
					r_1[1][i] = xa[m];
					r_1[2][i] = ya[m];
					r_1[3][i] = za[m];
					r_2[1][i] = xb[m];
					r_2[2][i] = yb[m];
					r_2[3][i] = zb[m];
					ka = ka + 1;
					k_ali[ka] = m;
					LL = LL + 1;
				}
				u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
				for (j = 1; j <= nseq; j++) {
					xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
					yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
					zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
				}
				score_fun(xt, yt, zt, xb, yb, zb, d*d, n_cut, n_ali, i_ali, d02, nseq, score);
				if (score_max < score) {
					score_max = score;
					ka0 = ka;
					for (i = 1; i <= ka; i++) {
						k_ali0[i] = k_ali[i];
					}
					
					for (i = 1; i < 4; i++) {
						for (j = 1; j < 4; j++) {
							out_u[i-1][j-1] = u[i][j];
						}
						out_t[i-1] = t[i];
					}
				}
				if (it == n_it) {
					break;
				}
				if (n_cut == ka) { // then
					int neq = 0;
					for (i = 1; i <= n_cut; i++) {
						if (i_ali[i] == k_ali[i]) {
							neq = neq + 1;
						}
					}
					if (n_cut == neq) {
						break;
					}
				}
			} 
		}
	}
	
	delete[] L_ini;
	delete[] i_ali;
	delete2Darr(u, 4);
	delete[] t;
	delete[] iq;
	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] xt;
	delete[] yt;
	delete[] zt;
	delete[] xb;
	delete[] yb;
	delete[] zb;
	delete2Darr(r_1, 4);
	delete2Darr(r_2, 4);
	delete[] k_ali;
	delete[] k_ali0;
	
	return score_max;
}

inline double CBaseFunc::cal_rot_tran_from_query_to_templ__(const vector<double*>& query, const vector<double*>& templ, const vector<MOLTYPE>& moltypes,
			 double** out_u, double* out_t, const double& user_d0_pro, const double& user_d0_dna, const double& user_d0_rna, const double& user_d0_lig, const bool& fast){
	
	int i, j, k;
	int nmax, nseq, n_ali;
	n_ali = nseq = query.size();
	nmax = nseq + 1;
	
	double d_pro, d2_pro, d0_pro, d02_pro, d0_search_pro;
	double d_dna, d2_dna, d0_dna, d02_dna, d0_search_dna;
	double d_rna, d2_rna, d0_rna, d02_rna, d0_search_rna;
	double d_lig, d2_lig, d0_lig, d02_lig, d0_search_lig;
	
	int n_cut; // ![1,n_ali],align residues for the score
	double score;
	double score_max;
	
	int* L_ini = new int[nmax];
	int* i_ali = new int[nmax];
	double** u = new2Darr(4, 4);
	double* t = new double[4];
	double* iq = new double[nmax];
	MOLTYPE* mt = new MOLTYPE[nmax];
	double* xa = new double[nmax];
	double* ya = new double[nmax];
	double* za = new double[nmax];
	double* xt = new double[nmax];
	double* yt = new double[nmax];
	double* zt = new double[nmax];
	double* xb = new double[nmax];
	double* yb = new double[nmax];
	double* zb = new double[nmax];
	double** r_1 = new2Darr(4, nmax);
	double** r_2 = new2Darr(4, nmax);
	int* k_ali = new int[nmax];
	int* k_ali0 = new int[nmax];
	
//		seq1A = "*" + p1.getSeq();
	for (int kk = 0; kk < nseq; kk++) {
		xa[kk+1] = query[kk][0];
		ya[kk+1] = query[kk][1];
		za[kk+1] = query[kk][2];
		xb[kk+1] = templ[kk][0];
		yb[kk+1] = templ[kk][1];
		zb[kk+1] = templ[kk][2];
		mt[kk+1] = moltypes[kk];
	}
	
	int ka0 = 0;
	d0_pro = user_d0_pro;
	d0_dna = user_d0_dna;
	d0_rna = user_d0_rna;
	d0_lig = user_d0_lig;
	
	d02_pro = d0_pro*d0_pro;
	d02_dna = d0_dna*d0_dna;
	d02_rna = d0_rna*d0_rna;
	d02_lig = d0_lig*d0_lig;
	
	// *** d0_search ----->
	d0_search_pro = d0_pro;
	if (d0_search_pro > 8)
		d0_search_pro = 8;
	if (d0_search_pro < 4.5)
		d0_search_pro = 4.5;
	
	d0_search_dna = d0_dna;
	if (d0_search_dna > 8)
		d0_search_dna = 8;
	if (d0_search_dna < 4.5)
		d0_search_dna = 4.5;
	
	d0_search_rna = d0_rna;
	if (d0_search_rna > 8)
		d0_search_rna = 8;
	if (d0_search_rna < 4.5)
		d0_search_rna = 4.5;
		
	d0_search_lig = d0_lig;
	if (d0_search_lig > 4)
		d0_search_lig = 4;
	if (d0_search_lig < 4.5)
		d0_search_lig = 4.5;
		
	// *** iterative parameters ----->
	//int n_it = 20; // !maximum number of iterations
	int n_it = g_maxmum_number_of_iterations;
	
	//int n_init_max = 6; // !maximum number of L_init
	int n_init_max = g_maxmum_number_of_L_init;
	
	int n_init = 0;
	int L_ini_min = 4;
	if (n_ali < 4)
		L_ini_min = n_ali;
	
	bool flag1 = false;
	for (i = 1; i <= (n_init_max - 1); i++) {
		n_init = n_init + 1;
		L_ini[n_init] = n_ali / (int) pow(2, (n_init - 1));
		if (L_ini[n_init] <= L_ini_min) {
			L_ini[n_init] = L_ini_min;
			flag1 = true;
			break;
		}
	}
	if (flag1 == false) {
		n_init = n_init + 1;
		L_ini[n_init] = L_ini_min;
	}

	score_max = 0; // !TM-score
	int LL, ka, i_init, L_init, iL_max, iL;
	for (i_init = 1; i_init <= n_init; i_init++) { // 333
		L_init = L_ini[i_init];
		iL_max = n_ali - L_init + 1;
		for (iL = 1; iL <= iL_max; fast? iL+=2 : iL++) { // 300 !on aligned residues,
			LL = 0;
			ka = 0;
			for (i = 1; i <= L_init; i++) {
				k = iL + i - 1; // ![1,n_ali] common
				r_1[1][i] = xa[k];
				r_1[2][i] = ya[k];
				r_1[3][i] = za[k];
				r_2[1][i] = xb[k];
				r_2[2][i] = yb[k];
				r_2[3][i] = zb[k];
				ka = ka + 1;
				k_ali[ka] = k;
				LL = LL + 1;
			}
			u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
			for (j = 1; j <= nseq; j++) {
				xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
				yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
				zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
			}
			d_pro = d0_search_pro - 1;
			d_dna = d0_search_dna - 1;
			d_rna = d0_search_rna - 1;
			d_lig = d0_search_lig - 1;
			score_fun(xt, yt, zt, xb, yb, zb, mt, d_pro*d_pro, d_dna*d_dna, d_rna*d_rna, d_lig*d_lig, n_cut, n_ali, i_ali, d02_pro, d02_dna, d02_rna, d02_lig, nseq, score); 
			// iteration
			if (score_max < score) {
				score_max = score;
				ka0 = ka;
				for (i = 1; i <= ka0; i++) {
					k_ali0[i] = k_ali[i];
				}
				
				for (i = 1; i < 4; i++) {
					for (j = 1; j < 4; j++) {
						out_u[i-1][j-1] = u[i][j];
					}
					out_t[i-1] = t[i];
				}
			}
			
			// *** iteration for extending
			// ---------------------------------->
			d_pro = d0_search_pro + 1;
			d_dna = d0_search_dna + 1;
			d_rna = d0_search_rna + 1;
			d_lig = d0_search_lig + 1;
			for (int it = 1; it <= n_it; it++) {
				LL = 0;
				ka = 0;
				for (i = 1; i <= n_cut; i++) {
					int m = i_ali[i]; // ![1,n_ali]
					r_1[1][i] = xa[m];
					r_1[2][i] = ya[m];
					r_1[3][i] = za[m];
					r_2[1][i] = xb[m];
					r_2[2][i] = yb[m];
					r_2[3][i] = zb[m];
					ka = ka + 1;
					k_ali[ka] = m;
					LL = LL + 1;
				}
				u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
				for (j = 1; j <= nseq; j++) {
					xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
					yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
					zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
				}
				score_fun(xt, yt, zt, xb, yb, zb, mt, d_pro*d_pro, d_dna*d_dna, d_rna*d_rna, d_lig*d_lig, n_cut, n_ali, i_ali, d02_pro, d02_dna, d02_rna, d02_lig, nseq, score);
				if (score_max < score) {
					score_max = score;
					ka0 = ka;
					for (i = 1; i <= ka; i++) {
						k_ali0[i] = k_ali[i];
					}
					
					for (i = 1; i < 4; i++) {
						for (j = 1; j < 4; j++) {
							out_u[i-1][j-1] = u[i][j];
						}
						out_t[i-1] = t[i];
					}
				}
				if (it == n_it) {
					break;
				}
				if (n_cut == ka) { // then
					int neq = 0;
					for (i = 1; i <= n_cut; i++) {
						if (i_ali[i] == k_ali[i]) {
							neq = neq + 1;
						}
					}
					if (n_cut == neq) {
						break;
					}
				}
			} 
		}
	}
	
	delete[] L_ini;
	delete[] i_ali;
	delete2Darr(u, 4);
	delete[] t;
	delete[] iq;
	delete[] mt;
	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] xt;
	delete[] yt;
	delete[] zt;
	delete[] xb;
	delete[] yb;
	delete[] zb;
	delete2Darr(r_1, 4);
	delete2Darr(r_2, 4);
	delete[] k_ali;
	delete[] k_ali0;
		
	return score_max;
}


inline double CBaseFunc::cal_rot_tran_from_query_to_templ__II(const vector<double*>& query, const vector<double*>& templ, const vector<MOLTYPE>& moltypes,
			 double** out_u, double* out_t, const double& user_d0_pro, const double& user_d0_dna, const double& user_d0_rna, const double& user_d0_lig, const bool& fast){
	
	int i, j, k;
	int nmax, nseq, n_ali;
	n_ali = nseq = query.size();
	nmax = nseq + 1;
	
	double d_pro, d2_pro, d0_pro, d02_pro, d0_search_pro;
	double d_dna, d2_dna, d0_dna, d02_dna, d0_search_dna;
	double d_rna, d2_rna, d0_rna, d02_rna, d0_search_rna;
	double d_lig, d2_lig, d0_lig, d02_lig, d0_search_lig;
	
	int n_cut; // ![1,n_ali],align residues for the score
	double score;
	double score_max;
	
	int* L_ini = new int[nmax];
	int* i_ali = new int[nmax];
	double** u = new2Darr(4, 4);
	double* t = new double[4];
	double* iq = new double[nmax];
	MOLTYPE* mt = new MOLTYPE[nmax];
	double* xa = new double[nmax];
	double* ya = new double[nmax];
	double* za = new double[nmax];
	double* xt = new double[nmax];
	double* yt = new double[nmax];
	double* zt = new double[nmax];
	double* xb = new double[nmax];
	double* yb = new double[nmax];
	double* zb = new double[nmax];
	double** r_1 = new2Darr(4, nmax);
	double** r_2 = new2Darr(4, nmax);
	int* k_ali = new int[nmax];
	int* k_ali0 = new int[nmax];
	
//		seq1A = "*" + p1.getSeq();
	for (int kk = 0; kk < nseq; kk++) {
		xa[kk+1] = query[kk][0];
		ya[kk+1] = query[kk][1];
		za[kk+1] = query[kk][2];
		xb[kk+1] = templ[kk][0];
		yb[kk+1] = templ[kk][1];
		zb[kk+1] = templ[kk][2];
		mt[kk+1] = moltypes[kk];
	}
	
	int ka0 = 0;
	d0_pro = user_d0_pro;
	d0_dna = user_d0_dna;
	d0_rna = user_d0_rna;
	d0_lig = user_d0_lig;
	
	d02_pro = d0_pro*d0_pro;
	d02_dna = d0_dna*d0_dna;
	d02_rna = d0_rna*d0_rna;
	d02_lig = d0_lig*d0_lig;
	
	// *** d0_search ----->
	d0_search_pro = d0_pro;
	if (d0_search_pro > 8)
		d0_search_pro = 8;
	if (d0_search_pro < 4.5)
		d0_search_pro = 4.5;
	
	d0_search_dna = d0_dna;
	if (d0_search_dna > 8)
		d0_search_dna = 8;
	if (d0_search_dna < 4.5)
		d0_search_dna = 4.5;
	
	d0_search_rna = d0_rna;
	if (d0_search_rna > 8)
		d0_search_rna = 8;
	if (d0_search_rna < 4.5)
		d0_search_rna = 4.5;
		
	d0_search_lig = d0_lig;
	if (d0_search_lig > 4)
		d0_search_lig = 4;
	if (d0_search_lig < 4.5)
		d0_search_lig = 4.5;
		
	// *** iterative parameters ----->
	//int n_it = 20; // !maximum number of iterations
	int n_it = g_maxmum_number_of_iterations;
	
	//int n_init_max = 6; // !maximum number of L_init
	int n_init_max = g_maxmum_number_of_L_init;
	
	int n_init = 0;
	int L_ini_min = 4;
	if (n_ali < 4)
		L_ini_min = n_ali;
	
	bool flag1 = false;
	for (i = 1; i <= (n_init_max - 1); i++) {
		n_init = n_init + 1;
		L_ini[n_init] = n_ali / (int) pow(2, (n_init - 1));
		if (L_ini[n_init] <= L_ini_min) {
			L_ini[n_init] = L_ini_min;
			flag1 = true;
			break;
		}
	}
	if (flag1 == false) {
		n_init = n_init + 1;
		L_ini[n_init] = L_ini_min;
	}

	score_max = 0; // !TM-score
	int LL, ka, i_init, L_init, iL_max, iL;
	for (i_init = 1; i_init <= n_init; i_init++) { // 333
		L_init = L_ini[i_init];
		iL_max = n_ali - L_init + 1;
		for (iL = 1; iL <= iL_max; fast? iL+=40 : iL++) { // 300 !on aligned residues,
			LL = 0;
			ka = 0;
			for (i = 1; i <= L_init; i++) {
				k = iL + i - 1; // ![1,n_ali] common
				r_1[1][i] = xa[k];
				r_1[2][i] = ya[k];
				r_1[3][i] = za[k];
				r_2[1][i] = xb[k];
				r_2[2][i] = yb[k];
				r_2[3][i] = zb[k];
				ka = ka + 1;
				k_ali[ka] = k;
				LL = LL + 1;
			}
			u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
			for (j = 1; j <= nseq; j++) {
				xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
				yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
				zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
			}
			d_pro = d0_search_pro - 1;
			d_dna = d0_search_dna - 1;
			d_rna = d0_search_rna - 1;
			d_lig = d0_search_lig - 1;
			score_fun(xt, yt, zt, xb, yb, zb, mt, d_pro*d_pro, d_dna*d_dna, d_rna*d_rna, d_lig*d_lig, n_cut, n_ali, i_ali, d02_pro, d02_dna, d02_rna, d02_lig, nseq, score); 
			// iteration
			if (score_max < score) {
				score_max = score;
				ka0 = ka;
				for (i = 1; i <= ka0; i++) {
					k_ali0[i] = k_ali[i];
				}
				
				for (i = 1; i < 4; i++) {
					for (j = 1; j < 4; j++) {
						out_u[i-1][j-1] = u[i][j];
					}
					out_t[i-1] = t[i];
				}
			}
			
			// *** iteration for extending
			// ---------------------------------->
			d_pro = d0_search_pro + 1;
			d_dna = d0_search_dna + 1;
			d_rna = d0_search_rna + 1;
			d_lig = d0_search_lig + 1;
			for (int it = 1; it <= n_it; it++) {
				LL = 0;
				ka = 0;
				for (i = 1; i <= n_cut; i++) {
					int m = i_ali[i]; // ![1,n_ali]
					r_1[1][i] = xa[m];
					r_1[2][i] = ya[m];
					r_1[3][i] = za[m];
					r_2[1][i] = xb[m];
					r_2[2][i] = yb[m];
					r_2[3][i] = zb[m];
					ka = ka + 1;
					k_ali[ka] = m;
					LL = LL + 1;
				}
				u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
				for (j = 1; j <= nseq; j++) {
					xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
					yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
					zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
				}
				score_fun(xt, yt, zt, xb, yb, zb, mt, d_pro*d_pro, d_dna*d_dna, d_rna*d_rna, d_lig*d_lig, n_cut, n_ali, i_ali, d02_pro, d02_dna, d02_rna, d02_lig, nseq, score);
				if (score_max < score) {
					score_max = score;
					ka0 = ka;
					for (i = 1; i <= ka; i++) {
						k_ali0[i] = k_ali[i];
					}
					
					for (i = 1; i < 4; i++) {
						for (j = 1; j < 4; j++) {
							out_u[i-1][j-1] = u[i][j];
						}
						out_t[i-1] = t[i];
					}
				}
				if (it == n_it) {
					break;
				}
				if (n_cut == ka) { // then
					int neq = 0;
					for (i = 1; i <= n_cut; i++) {
						if (i_ali[i] == k_ali[i]) {
							neq = neq + 1;
						}
					}
					if (n_cut == neq) {
						break;
					}
				}
			} 
		}
	}
	
	delete[] L_ini;
	delete[] i_ali;
	delete2Darr(u, 4);
	delete[] t;
	delete[] iq;
	delete[] mt;
	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] xt;
	delete[] yt;
	delete[] zt;
	delete[] xb;
	delete[] yb;
	delete[] zb;
	delete2Darr(r_1, 4);
	delete2Darr(r_2, 4);
	delete[] k_ali;
	delete[] k_ali0;
		
	return score_max;
}


inline double CBaseFunc::cal_rot_tran_from_query_to_templ__for_rTMscore(
			const int& aligned_chain_num, 
			const int& chain_num, 
            const vector<double*>& query, const vector<double*>& templ, 
			map<int, int>& chain_index_corr_to_query__aa_num,
			map<int, MOLTYPE>& chain_index_corr_to_query__moltype,
            map<int, double>& chain_index_corr_to_query__d0,
            map<int, double>& chain_index_corr_to_query__d02,
			const vector<int>& chain_index_corr_to_query, 
			double** out_u, double* out_t, const bool& fast){
	
	int i, j, k;
	int nmax, nseq, n_ali;
	n_ali = nseq = query.size();
	nmax = nseq + 1;
	
	map<int, double>::iterator it;
	map<int, double> chain_index_corr_to_query__d2;
	
	int n_cut; // ![1,n_ali],align residues for the score
	double score;
	double score_max;
	
	int* L_ini = new int[nmax];
	int* i_ali = new int[nmax];
	double** u = new2Darr(4, 4);
	double* t = new double[4];
	double* iq = new double[nmax];
	double* xa = new double[nmax];
	double* ya = new double[nmax];
	double* za = new double[nmax];
	int* __chain_index = new int[nmax];
	double* xt = new double[nmax];
	double* yt = new double[nmax];
	double* zt = new double[nmax];
	double* xb = new double[nmax];
	double* yb = new double[nmax];
	double* zb = new double[nmax];
	double** r_1 = new2Darr(4, nmax);
	double** r_2 = new2Darr(4, nmax);
	int* k_ali = new int[nmax];
	int* k_ali0 = new int[nmax];
	
//		seq1A = "*" + p1.getSeq();
	for (int kk = 0; kk < nseq; kk++) {
		xa[kk+1] = query[kk][0];
		ya[kk+1] = query[kk][1];
		za[kk+1] = query[kk][2];
		xb[kk+1] = templ[kk][0];
		yb[kk+1] = templ[kk][1];
		zb[kk+1] = templ[kk][2];
		__chain_index[kk+1] = chain_index_corr_to_query[kk];
	}
	
	int ka0 = 0;
	
	// *** d0_search ----->
	map<int, double> d0_search;
	for (it = chain_index_corr_to_query__d0.begin(); it != chain_index_corr_to_query__d0.end(); it++){
		int cind = it->first;
		MOLTYPE& mt = chain_index_corr_to_query__moltype[it->first];
		switch (mt){
			case PROTEIN:
				d0_search[cind] = chain_index_corr_to_query__d0[cind];
				if (d0_search[cind] > 8)
					d0_search[cind] = 8;
				if (d0_search[cind] < 4.5)
					d0_search[cind] = 4.5;
				break;
			case DNA:
			case RNA:
				d0_search[cind] = chain_index_corr_to_query__d0[cind];
				if (d0_search[cind] > 8)
					d0_search[cind] = 8;
				if (d0_search[cind] < 4.5)
					d0_search[cind] = 4.5;
				break;
			case LIGAND:
				d0_search[cind] = chain_index_corr_to_query__d0[cind];
				if (d0_search[cind] > 8)
					d0_search[cind] = 8;
				if (d0_search[cind] < 4.5)
					d0_search[cind] = 4.5;
				break;
					
		}
	}
		
	// *** iterative parameters ----->
	//int n_it = 20; // !maximum number of iterations
	int n_it = g_maxmum_number_of_iterations;
	
	//int n_init_max = 6; // !maximum number of L_init
	int n_init_max = g_maxmum_number_of_L_init;
	
	int n_init = 0;
	int L_ini_min = 4;
	if (n_ali < 4)
		L_ini_min = n_ali;
	
	bool flag1 = false;
	for (i = 1; i <= (n_init_max - 1); i++) {
		n_init = n_init + 1;
		L_ini[n_init] = n_ali / (int) pow(2, (n_init - 1));
		if (L_ini[n_init] <= L_ini_min) {
			L_ini[n_init] = L_ini_min;
			flag1 = true;
			break;
		}
	}
	if (flag1 == false) {
		n_init = n_init + 1;
		L_ini[n_init] = L_ini_min;
	}

	score_max = 0; // !TM-score
	int LL, ka, i_init, L_init, iL_max, iL;
	for (i_init = 1; i_init <= n_init; i_init++) { // 333
		L_init = L_ini[i_init];
		iL_max = n_ali - L_init + 1;
		for (iL = 1; iL <= iL_max; fast ? iL+=2 : iL++) { // 300 !on aligned residues,
			LL = 0;
			ka = 0;
			for (i = 1; i <= L_init; i++) {
				k = iL + i - 1; // ![1,n_ali] common
				r_1[1][i] = xa[k];
				r_1[2][i] = ya[k];
				r_1[3][i] = za[k];
				r_2[1][i] = xb[k];
				r_2[2][i] = yb[k];
				r_2[3][i] = zb[k];
				ka = ka + 1;
				k_ali[ka] = k;
				LL = LL + 1;
			}
			u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
			for (j = 1; j <= nseq; j++) {
				xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
				yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
				zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
			}
			
			map<int, double>().swap(chain_index_corr_to_query__d2);
			for (it = d0_search.begin(); it != d0_search.end(); it++){
				int cind = it->first;
				MOLTYPE& mt = chain_index_corr_to_query__moltype[it->first];
				switch (mt){
					case PROTEIN:
						chain_index_corr_to_query__d2[cind] = d0_search[cind] - 1;
						break;
					case DNA:
					case RNA:
						chain_index_corr_to_query__d2[cind] = d0_search[cind] - 1;
						break;
					case LIGAND:
						chain_index_corr_to_query__d2[cind] = d0_search[cind] - 1;
						break;
				}
				chain_index_corr_to_query__d2[cind] = chain_index_corr_to_query__d2[cind]*chain_index_corr_to_query__d2[cind];
			} 
			
			score_fun_rtmsco(aligned_chain_num, chain_num, xt, yt, zt, chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d2, chain_index_corr_to_query__d02, __chain_index, xb, yb, zb, n_cut, n_ali, i_ali, score); 
			// iteration
			if (score_max < score) {
				score_max = score;
				ka0 = ka;
				for (i = 1; i <= ka0; i++) {
					k_ali0[i] = k_ali[i];
				}
				
				for (i = 1; i < 4; i++) {
					for (j = 1; j < 4; j++) {
						out_u[i-1][j-1] = u[i][j];
					}
					out_t[i-1] = t[i];
				}
			}
			
			// *** iteration for extending
			// ---------------------------------->
			map<int, double>().swap(chain_index_corr_to_query__d2);
			for (it = d0_search.begin(); it != d0_search.end(); it++){
				int cind = it->first;
				MOLTYPE& mt = chain_index_corr_to_query__moltype[it->first];
				switch (mt){
					case PROTEIN:
						chain_index_corr_to_query__d2[cind] = d0_search[cind] + 1;
						break;
					case DNA:
					case RNA:
						chain_index_corr_to_query__d2[cind] = d0_search[cind] + 1;
						break;
					case LIGAND:
						chain_index_corr_to_query__d2[cind] = d0_search[cind] + 1;
						break;
				}
				chain_index_corr_to_query__d2[cind] = chain_index_corr_to_query__d2[cind]*chain_index_corr_to_query__d2[cind];
			} 
			
			for (int it = 1; it <= n_it; it++) {
				LL = 0;
				ka = 0;
				for (i = 1; i <= n_cut; i++) {
					int m = i_ali[i]; // ![1,n_ali]
					r_1[1][i] = xa[m];
					r_1[2][i] = ya[m];
					r_1[3][i] = za[m];
					r_2[1][i] = xb[m];
					r_2[2][i] = yb[m];
					r_2[3][i] = zb[m];
					ka = ka + 1;
					k_ali[ka] = m;
					LL = LL + 1;
				}
				u3b(r_1, r_2, LL, 1, u, t); // !u rotate r_1 to r_2
				for (j = 1; j <= nseq; j++) {
					xt[j] = t[1] + u[1][1] * xa[j] + u[1][2] * ya[j] + u[1][3] * za[j];
					yt[j] = t[2] + u[2][1] * xa[j] + u[2][2] * ya[j] + u[2][3] * za[j];
					zt[j] = t[3] + u[3][1] * xa[j] + u[3][2] * ya[j] + u[3][3] * za[j];
				}
				score_fun_rtmsco(aligned_chain_num, chain_num, xt, yt, zt,  chain_index_corr_to_query__aa_num, chain_index_corr_to_query__d2, chain_index_corr_to_query__d02, __chain_index, xb, yb, zb, n_cut, n_ali, i_ali, score);
				if (score_max < score) {
					score_max = score;
					ka0 = ka;
					for (i = 1; i <= ka; i++) {
						k_ali0[i] = k_ali[i];
					}
					
					for (i = 1; i < 4; i++) {
						for (j = 1; j < 4; j++) {
							out_u[i-1][j-1] = u[i][j];
						}
						out_t[i-1] = t[i];
					}
				}
				if (it == n_it) {
					break;
				}
				if (n_cut == ka) { // then
					int neq = 0;
					for (i = 1; i <= n_cut; i++) {
						if (i_ali[i] == k_ali[i]) {
							neq = neq + 1;
						}
					}
					if (n_cut == neq) {
						break;
					}
				}
			} 
		}
	}
	
	delete[] L_ini;
	delete[] i_ali;
	delete2Darr(u, 4);
	delete[] t;
	delete[] iq;
	delete[] xa;
	delete[] ya;
	delete[] za;
	delete[] __chain_index;
	delete[] xt;
	delete[] yt;
	delete[] zt;
	delete[] xb;
	delete[] yb;
	delete[] zb;
	delete2Darr(r_1, 4);
	delete2Darr(r_2, 4);
	delete[] k_ali;
	delete[] k_ali0;
		
	return score_max;
}

inline void CBaseFunc::u3b(double** x, double** y, const int& n, const int& mode, double** u, double* t) {
	// cccccccccccccccc Calculate sum of (r_d-r_m)^2
	// cccccccccccccccccccccccccc
	// c w - w(m) is weight for atom pair c m (given)
	// c x - x(i,m) are coordinates of atom c m in set x (given)
	// c y - y(i,m) are coordinates of atom c m in set y (given)
	// c n - n is number of atom pairs (given)
	// c mode - 0:calculate rms only (given)
	// c 1:calculate rms,u,t (takes longer)
	// c rms - sum of w*(ux+t-y)**2 over all atom pairs (result)
	// c u - u(i,j) is rotation matrix for best superposition (result)
	// c t - t(i) is translation vector for best superposition (result)
	// c ier - 0: a unique optimal superposition has been determined(result)
	// c -1: superposition is not unique but optimal
	// c -2: no result obtained because of negative weights w
	// c or all weights equal to zero.
	// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	int ier = -1;
	if (n < 1) {
		return ;
	}
	ier = -2;
	
	int i, j, m, m1, l, k;
	double rms, e0, d, h, g;
	double cth, sth, sqrth, p, det, sigma;
	double* xc = new double[3 + 1];
	double* yc = new double[3 + 1];
	double** a = new2Darr(3 + 1, 3 + 1);
	double** b = new2Darr(3 + 1, 3 + 1);
	double** r = new2Darr(3 + 1, 3 + 1);
	double* e = new double[3 + 1];
	double* rr = new double[6 + 1];
	double* ss = new double[6 + 1];
	
	double sqrt3 = 1.73205080756888, tol = 0.01;
	int ip[] = { -100, 1, 2, 4, 2, 3, 5, 4, 5, 6 };
	int ip2312[] = { -100, 2, 3, 1, 2 };
	int a_failed = 0, b_failed = 0;
	double epsilon = 0.000000001;

	// initializtation
	rms = 0;
	e0 = 0;
	for (i = 1; i <= 3; i++) {
		xc[i] = 0.0;
		yc[i] = 0.0;
		t[i] = 0.0;
		for (j = 1; j <= 3; j++) {
			u[i][j] = 0.0;
			r[i][j] = 0.0;
			a[i][j] = 0.0;
			if (i == j) {
				u[i][j] = 1.0;
				a[i][j] = 1.0;
			}
		}
	}
	
	// compute centers for vector sets x, y
	for (m = 1; m <= n; m++) {
		for (i = 1; i <= 3; i++) {
			xc[i] = xc[i] + x[i][m];
			yc[i] = yc[i] + y[i][m];
		}
	}
	for (i = 1; i <= 3; i++) {
		xc[i] = xc[i] / n;
		yc[i] = yc[i] / n;
	}
	// compute e0 and matrix r
	for (m = 1; m <= n; m++) {
		for (i = 1; i <= 3; i++) {
			e0 = e0 + (x[i][m] - xc[i]) * (x[i][m] - xc[i]) + (y[i][m] - yc[i]) * (y[i][m] - yc[i]);
			d = y[i][m] - yc[i];
			for (j = 1; j <= 3; j++) {
				r[i][j] = r[i][j] + d * (x[j][m] - xc[j]);
			}
		}
	}
	// compute determinat of matrix r
	det = r[1][1] * ((r[2][2] * r[3][3]) - (r[2][3] * r[3][2]))
			- r[1][2] * ((r[2][1] * r[3][3]) - (r[2][3] * r[3][1]))
			+ r[1][3] * ((r[2][1] * r[3][2]) - (r[2][2] * r[3][1]));

	sigma = det;
	// compute tras(r)*r
	m = 0;
	for (j = 1; j <= 3; j++) {
		for (i = 1; i <= j; i++) {
			m = m + 1;
			rr[m] = r[1][i] * r[1][j] + r[2][i] * r[2][j] + r[3][i] * r[3][j];
		}
	}

	double spur = (rr[1] + rr[3] + rr[6]) / 3.0;
	double cof = (((((rr[3] * rr[6] - rr[5] * rr[5]) + rr[1] * rr[6]) - rr[4] * rr[4]) + rr[1] * rr[3])
			- rr[2] * rr[2]) / 3.0;
	det = det * det;

	for (i = 1; i <= 3; i++) {
		e[i] = spur;
	}
	
	if (spur > 0) {
		d = spur * spur;
		h = d - cof;
		g = (spur * cof - det) / 2.0 - spur * h;
		if (h > 0) {
			sqrth = sqrt(h);
			d = h * h * h - g * g;
			if (d < 0.0) {
				d = 0.0;
			}
			d = atan2(sqrt(d), -g) / 3.0;
			cth = sqrth * cos(d);
			sth = sqrth * sqrt3 * sin(d);
			e[1] = (spur + cth) + cth;
			e[2] = (spur - cth) + sth;
			e[3] = (spur - cth) - sth;
			if (mode != 0) {// compute a
				for (l = 1; l <= 3; l = l + 2) {
					d = e[l];
					ss[1] = (d - rr[3]) * (d - rr[6]) - rr[5] * rr[5];
					ss[2] = (d - rr[6]) * rr[2] + rr[4] * rr[5];
					ss[3] = (d - rr[1]) * (d - rr[6]) - rr[4] * rr[4];
					ss[4] = (d - rr[3]) * rr[4] + rr[2] * rr[5];
					ss[5] = (d - rr[1]) * rr[5] + rr[2] * rr[4];
					ss[6] = (d - rr[1]) * (d - rr[3]) - rr[2] * rr[2];

					if (fabs(ss[0]) <= epsilon)
						ss[0] = 0.0;
					if (fabs(ss[1]) <= epsilon)
						ss[1] = 0.0;
					if (fabs(ss[2]) <= epsilon)
						ss[2] = 0.0;
					if (fabs(ss[3]) <= epsilon)
						ss[3] = 0.0;
					if (fabs(ss[4]) <= epsilon)
						ss[4] = 0.0;
					if (fabs(ss[5]) <= epsilon)
						ss[5] = 0.0;

					if (fabs(ss[1]) >= fabs(ss[3])) {
						j = 1;
						if (fabs(ss[1]) < fabs(ss[6])) {
							j = 3;
						}
					} else if (fabs(ss[3]) >= fabs(ss[6])) {
						j = 2;
					} else {
						j = 3;
					}

					d = 0.0;
					j = 3 * (j - 1);
					for (i = 1; i <= 3; i++) {
						k = ip[i + j];
						a[i][l] = ss[k];
						d = d + ss[k] * ss[k];
					}
					// if( d > 0.0 ) d = 1.0 / sqrt(d);
					if (d > 0)
						d = 1.0 / sqrt(d);
					else
						d = 0.0;
					for (i = 1; i <= 3; i++) {
						a[i][l] = a[i][l] * d;
					}
				} // for l

				d = a[1][1] * a[1][3] + a[2][1] * a[2][3] + a[3][1] * a[3][3];

				if ((e[1] - e[2]) > (e[2] - e[3])) {
					m1 = 3;
					m = 1;
				} else {
					m1 = 1;
					m = 3;
				}
				p = 0;
				for (i = 1; i <= 3; i++) {
					a[i][m1] = a[i][m1] - d * a[i][m];
					p = p + a[i][m1] * a[i][m1];
				}

				if (p <= tol) {
					p = 1.0;
					for (i = 1; i <= 3; i++) {
						if (p < fabs(a[i][m])) {
							continue;
						}
						p = fabs(a[i][m]);
						j = i;
					}
					k = ip2312[j];
					l = ip2312[j + 1];
					p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
					if (p > tol) {
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m] / p;
						a[l][m1] = a[k][m] / p;
					} else {// goto 40
						a_failed = 1;
					}
				} // if p<=tol
				else {
					p = 1.0 / sqrt(p);
					for (i = 1; i <= 3; i++) {
						a[i][m1] = a[i][m1] * p;
					}
				} // else p<=tol
				if (a_failed != 1) {
					a[1][2] = a[2][3] * a[3][1] - a[2][1] * a[3][3];
					a[2][2] = a[3][3] * a[1][1] - a[3][1] * a[1][3];
					a[3][2] = a[1][3] * a[2][1] - a[1][1] * a[2][3];
				}
			} // if(mode!=0)
		} // h>0

		// compute b anyway
		if (mode != 0 && a_failed != 1) {// a is computed correctly
			// compute b
			for (l = 1; l <= 2; l++) {
				d = 0.0;
				for (i = 1; i <= 3; i++) {
					b[i][l] = r[i][1] * a[1][l] + r[i][2] * a[2][l] + r[i][3] * a[3][l];
					d = d + b[i][l] * b[i][l];
				}
				// if( d > 0 ) d = 1.0 / sqrt(d);
				if (d > 0)
					d = 1.0 / sqrt(d);
				else
					d = 0.0;
				for (i = 1; i <= 3; i++) {
					b[i][l] = b[i][l] * d;
				}
			}
			d = b[1][1] * b[1][2] + b[2][1] * b[2][2] + b[3][1] * b[3][2];
			p = 0.0;

			for (i = 1; i <= 3; i++) {
				b[i][2] = b[i][2] - d * b[i][1];
				p = p + b[i][2] * b[i][2];
			}

			if (p <= tol) {
				p = 1.0;
				for (i = 1; i <= 3; i++) {
					if (p < fabs(b[i][1])) {
						continue;
					}
					p = fabs(b[i][1]);
					j = i;
				}
				k = ip2312[j];
				l = ip2312[j + 1];
				p = sqrt(b[k][1] * b[k][1] + b[l][1] * b[l][1]);
				if (p > tol) {
					b[j][2] = 0.0;
					b[k][2] = -b[l][1] / p;
					b[l][2] = b[k][1] / p;
				} else {
					b_failed = 1;
				}
			} // if( p <= tol )
			else {
				p = 1.0 / sqrt(p);
				for (i = 1; i <= 3; i++) {
					b[i][2] = b[i][2] * p;
				}
			}
			if (b_failed != 1) {
				b[1][3] = b[2][1] * b[3][2] - b[2][2] * b[3][1];
				b[2][3] = b[3][1] * b[1][2] - b[3][2] * b[1][1];
				b[3][3] = b[1][1] * b[2][2] - b[1][2] * b[2][1];
				// compute u
				for (i = 1; i <= 3; i++) {
					for (j = 1; j <= 3; j++) {
						u[i][j] = b[i][1] * a[j][1] + b[i][2] * a[j][2] + b[i][3] * a[j][3];
					}
				}
			}

			// compute t
			for (i = 1; i <= 3; i++) {
				t[i] = ((yc[i] - u[i][1] * xc[1]) - u[i][2] * xc[2]) - u[i][3] * xc[3];
			}
		} // if(mode!=0 && a_failed!=1)
	} // spur>0
	else // just compute t and errors
	{
		// compute t
		for (i = 1; i <= 3; i++) {
			t[i] = ((yc[i] - u[i][1] * xc[1]) - u[i][2] * xc[2]) - u[i][3] * xc[3];
		}
	} // else spur>0

	// compute rms
	for (i = 1; i <= 3; i++) {
		if (e[i] < 0)
			e[i] = 0;
		e[i] = sqrt(e[i]);
	}
	ier = 0;
	if (e[2] <= (e[1] * 1.0e-05))
		ier = -1;
	d = e[3];
	if (sigma < 0.0) {
		d = -d;
		if ((e[2] - e[3]) <= (e[1] * 1.0e-05))
			ier = -1;
	}
	d = (d + e[2]) + e[1];
	rms = (e0 - d) - d;
	if (rms < 0.00000000001)
		rms = 0.0;

	delete[] xc;
	delete[] yc;
	CBaseFunc::delete2Darr(a, 3+1);
	CBaseFunc::delete2Darr(b, 3+1);
	CBaseFunc::delete2Darr(r, 3+1);
	delete[] e;
	delete[] rr;
	delete[] ss;
}

inline void CBaseFunc::score_fun(double* xt, double* yt, double* zt, double* xb, double* yb, double* zb, const double& d2, int& n_cut, int& n_ali, int* i_ali, const double& d02, const int& nseq, double& score) {
	double d_tmp2 = d2;
	double score_sum = 0; // !TMscore
	double dis2 = 0;
	int k; 
	double xbias, ybias, zbias;
	double d_tmp; 
	
	while (true) {
		n_cut = 0; // !number of residue-pairs dis<d, for iteration
		score_sum = 0; // !TMscore
		dis2 = 0;
		for (k = 1; k <= n_ali; k++) {
			xbias = xt[k] - xb[k];
			ybias = yt[k] - yb[k];
			zbias = zt[k] - zb[k];
			dis2 = xbias*xbias + ybias*ybias + zbias*zbias;
			
			// *** for iteration:
			if (dis2 < d_tmp2) {
				n_cut++;
				i_ali[n_cut] = k; // ![1,n_ali], mark the residue-pairs in dis<d
			}
			
			// *** for TM-score:
			score_sum += 1. / (1 + dis2/d02);
		}
		
		if (n_cut < 3 && n_ali > 3){
			d_tmp = sqrt(d_tmp2) + 0.5;
			d_tmp2 = d_tmp*d_tmp;
		}
			
		else break;
	}
	
	score = score_sum / nseq; // !TM-score
}

inline void CBaseFunc::score_fun_rtmsco(
			const int& aligned_chain_num,
			const int& chain_num,
			double* xt, double* yt, double* zt, 
			map<int, int>& chain_index_corr_to_query__aa_num,
            map<int, double>& chain_index_corr_to_query__d2,
            map<int, double>& chain_index_corr_to_query__d02,
			int* chain_index_corr_to_query, 
			double* xb, double* yb, double* zb,
			int& n_cut, int& n_ali, int* i_ali, double& score) {
	map<int, double> chain_scores;
	map<int, double>::iterator it;
	
	map<int, double> d_tmp2;
	for (it = chain_index_corr_to_query__d2.begin(); it != chain_index_corr_to_query__d2.end(); it++)
		d_tmp2[it->first] = it->second;
	
	int chain_index;
	double sco = 0;
	double dis2 = 0;
	int k; 
	double xbias, ybias, zbias;
	double d_tmp;
	
	while (true) {
		n_cut = 0; // !number of residue-pairs dis<d, for iteration
		map<int, double>().swap(chain_scores); // !TMscore
		dis2 = 0;
		for (k = 1; k <= n_ali; k++) {
			chain_index = chain_index_corr_to_query[k];
			
			xbias = xt[k] - xb[k];
			ybias = yt[k] - yb[k];
			zbias = zt[k] - zb[k];
			dis2 = xbias*xbias + ybias*ybias + zbias*zbias;
			
			// *** for iteration:
			if (dis2 < d_tmp2[chain_index]) {
				n_cut++;
				i_ali[n_cut] = k; // ![1,n_ali], mark the residue-pairs in dis<d
			}
			
			// *** for TM-score:
			sco = 1. / (1 + dis2/chain_index_corr_to_query__d02[chain_index]);
			
			it = chain_scores.find(chain_index);
			if (chain_scores.end() == it){
				chain_scores[chain_index] = sco;
			}else{
				it->second += sco;
			}
		}
		
		if (n_cut < 3 && n_ali > 3){
			d_tmp = sqrt(d_tmp2[chain_index]) + 0.5;
			d_tmp2[chain_index] = d_tmp * d_tmp;
		}
		else break;
	}
	
	score = 0.;
	for (it = chain_scores.begin(); it != chain_scores.end(); it++){
		chain_index = it->first;
		sco = it->second / chain_index_corr_to_query__aa_num[chain_index];
		score += 1./ sco;
	}
	
	score = aligned_chain_num * aligned_chain_num / (chain_num * score) ; // rTM-score
}

inline double CBaseFunc::score_fun_rtmsco_once(
		const int& aligned_chain_num,
		const int& chain_num,
		const vector<double*>& rotted_qxyz, 
		const vector<double*>& txyz, 
		const vector<int>& chain_index_corr_to_query, 
		map<int, int>& chain_index_corr_to_query__aa_num,
		map<int, double>& chain_index_corr_to_query__d02) {
	map<int, double> chain_scores;
	map<int, double>::iterator it;
	
	int ali_len = rotted_qxyz.size();
	int chain_index;
	double sco = 0;
	double dis2 = 0;
	int k; 
	double xbias, ybias, zbias;
	double d_tmp;
	
	dis2 = 0;
	for (k = 0; k < ali_len; k++) {
		chain_index = chain_index_corr_to_query[k];
	
		xbias = rotted_qxyz[k][0] - txyz[k][0];
		ybias = rotted_qxyz[k][1] - txyz[k][1];
		zbias = rotted_qxyz[k][2] - txyz[k][2];
		dis2 = xbias*xbias + ybias*ybias + zbias*zbias;
		
		// *** for TM-score:
		sco = 1. / (1 + dis2/chain_index_corr_to_query__d02[chain_index]);
		
		it = chain_scores.find(chain_index);
		if (chain_scores.end() == it){
			chain_scores[chain_index] = sco;
		}else{
			it->second += sco;
		}			
	}
	
	double score = 0.;
	for (it = chain_scores.begin(); it != chain_scores.end(); it++){
		chain_index = it->first;
		sco = it->second / chain_index_corr_to_query__aa_num[chain_index];
		score += 1./ sco;
	}
	
	return aligned_chain_num * aligned_chain_num / (chain_num * score) ; // rTM-score
}

inline double CBaseFunc::score_fun_once(
		const vector<double*>& rotted_qxyz, 
		const vector<double*>& txyz, 
		const vector<MOLTYPE>& mts,
		const double& d0_pro2, const double& d0_dna2, const double& d0_rna2, const double& d0_lig2,
		const int& total_res_num) {
	int ali_len = rotted_qxyz.size();
	double score_sum = 0; // !TMscore
	double dis2 = 0;
	int k; 
	double xbias, ybias, zbias;
	
	for (k = 0; k < ali_len; k++) {
		xbias = rotted_qxyz[k][0] - txyz[k][0];
		ybias = rotted_qxyz[k][1] - txyz[k][1];
		zbias = rotted_qxyz[k][2] - txyz[k][2];
		dis2 = xbias*xbias + ybias*ybias + zbias*zbias;
		
		if (PROTEIN == mts[k]){
			score_sum += 1. / (1 + dis2/d0_pro2);
		}else if (DNA == mts[k]){
			score_sum += 1. / (1 + dis2/d0_dna2);
		}else if (RNA == mts[k]){
			score_sum += 1. / (1 + dis2/d0_rna2);
		}else {
			score_sum += 1. / (1. + dis2/d0_lig2);
		}
	}
		
	return score_sum / total_res_num; // !TM-score
}

inline void CBaseFunc::score_fun(double* xt, double* yt, double* zt, double* xb, double* yb, double* zb, MOLTYPE* mt, 
								const double& d_pro2, const double& d_dna2, const double& d_rna2, const double& d_lig2, 
								int& n_cut, int& n_ali, int* i_ali, 
								const double& d0_pro2, const double& d0_dna2, const double& d0_rna2, const double& d0_lig2, 
								const int& nseq, double& score) {
	double d_tmp_pro2 = d_pro2;
	double d_tmp_dna2 = d_dna2;
	double d_tmp_rna2 = d_rna2;
	double d_tmp_lig2 = d_lig2;
	
	double score_sum = 0; // !TMscore
	double dis2 = 0;
	int k; 
	double xbias, ybias, zbias;
	double d_tmp;
	
	while (true) {
		n_cut = 0; // !number of residue-pairs dis<d, for iteration
		score_sum = 0; // !TMscore
		dis2 = 0;
		for (k = 1; k <= n_ali; k++) {
			xbias = xt[k] - xb[k];
			ybias = yt[k] - yb[k];
			zbias = zt[k] - zb[k];
			dis2 = xbias*xbias + ybias*ybias + zbias*zbias;
			
			if (PROTEIN == mt[k]){
				// *** for iteration:
				if (dis2 < d_tmp_pro2) {
					n_cut++;
					i_ali[n_cut] = k; // ![1,n_ali], mark the residue-pairs in dis<d
				}
				
				// *** for TM-score:
				score_sum += 1. / (1 + dis2/d0_pro2);
			}else if (DNA == mt[k]){
				// *** for iteration:
				if (dis2 < d_tmp_dna2) {
					n_cut++;
					i_ali[n_cut] = k; // ![1,n_ali], mark the residue-pairs in dis<d
				}
				
				// *** for TM-score:
				score_sum += 1. / (1 + dis2/d0_dna2);
			}else if (RNA == mt[k]){
				// *** for iteration:
				if (dis2 < d_tmp_rna2) {
					n_cut++;
					i_ali[n_cut] = k; // ![1,n_ali], mark the residue-pairs in dis<d
				}
				
				// *** for TM-score:
				score_sum += 1. / (1 + dis2/d0_rna2);
			}else {
				// *** for iteration:
				if (dis2 < d_tmp_lig2) {
					n_cut++;
					i_ali[n_cut] = k; // ![1,n_ali], mark the residue-pairs in dis<d
				}
				
				// *** for TM-score:
				score_sum += 1. / (1. + dis2/d0_lig2);
			}
		}
		
		if (n_cut < 3 && n_ali > 3){
			d_tmp = sqrt(d_tmp_pro2) + 0.5;
			d_tmp_pro2 = d_tmp * d_tmp;
			
			d_tmp = sqrt(d_tmp_dna2) + 0.5;
			d_tmp_dna2 = d_tmp * d_tmp;
			
			d_tmp = sqrt(d_tmp_rna2) + 0.5;
			d_tmp_rna2 = d_tmp * d_tmp;
			
			d_tmp = sqrt(d_tmp_lig2) + 0.5;
			d_tmp_lig2 = d_tmp * d_tmp;
		}
		else break;
	}
	
	score = score_sum / nseq; // !TM-score
}

inline Molecule::Molecule(const MOLTYPE& moltype, const vector<string>* contents) {
	this->m_moltype = moltype;
	if (moltype == PROTEIN)
		load_one_protein(contents);
	else if (moltype == RNA || moltype == DNA)
		load_one_nucletide(contents);
	else if (moltype == LIGAND)
		load_one_ligand(contents);
	else {
		cout << "Oops, you are not the molecule types of protein, rna, dna, ligand..." << endl;
		exit(1);
	}
	
	if (this->m_moltype != LIGAND)
		m_seq_str = CBaseFunc::charVec2String(m_seq_vec);
	else {
//		int i, n = m_cared_atomtype_vec.size();
//		for (i = 0; i < n-1; i++){
//			m_seq_str += m_cared_atomtype_vec[i] + " ";
//		}
//		if (n > 0) m_seq_str += m_cared_atomtype_vec[n-1];

		int i, n = m_cared_atomsimpletype_vec.size();
		for (i = 0; i < n-1; i++){
			m_seq_str += m_cared_atomsimpletype_vec[i] + " ";
		}
		if (n > 0) m_seq_str += m_cared_atomsimpletype_vec[n-1];
	}
}

inline void Molecule::load_one_protein(const vector<string>* contents){
	int size = contents->size();
	int line = 0;
	int llen; 
	string lline;
	double x, y, z;
	double* xyz = NULL;
	char cstr[50];
	int res_orig_index;
	char char_following_res_orig_index;
	
	while (size > line){
		lline = (*contents)[line];
		CBaseFunc::toUpperString(lline);
		
		llen = lline.size();
		
		if (llen < 54) {
			line++;
			continue;
		}
		
		strcpy(cstr, (lline.substr(30, 8)).c_str());
		sscanf(cstr, "%lf", &x);
		strcpy(cstr, (lline.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", &y);
		strcpy(cstr, (lline.substr(46, 8)).c_str());
		sscanf(cstr, "%lf", &z);
		
		xyz = new double[3];
		xyz[0] = x;
		xyz[1] = y;
		xyz[2] = z;
		
		if (0 == lline.compare(12, 4, g_atomtype_res)){
			strcpy(cstr, (lline.substr(22, 4)).c_str());
			sscanf(cstr, "%d", &res_orig_index);
			char_following_res_orig_index = lline[26];
			
			if (0 == m_orig_index_vec.size() || res_orig_index != m_orig_index_vec[m_orig_index_vec.size()-1] 
					|| char_following_res_orig_index != m_char_following_orig_index_vec[m_char_following_orig_index_vec.size()-1]){
				m_cared_xyz_vec.push_back(xyz);
				m_orig_index_vec.push_back(res_orig_index);
				m_char_following_orig_index_vec.push_back(char_following_res_orig_index);
				m_seq_vec.push_back(CBaseFunc::aminoAcidAbbr3WordsTo1(lline.substr(17, 3)));
			}else{
				cout << "Warning! Duplicated residue " << res_orig_index << "." << endl; 
				
				bool is_not_clash_CA = true;
				if (m_cared_xyz_vec.size() != 0){
					double dis2 = CBaseFunc::distance2(xyz, m_cared_xyz_vec[m_cared_xyz_vec.size()-1]);
					if (dis2 < 14/*~=3.75*3.75*/)
						is_not_clash_CA = false;
				}
				
				if (is_not_clash_CA){
					m_cared_xyz_vec.push_back(xyz);
					m_orig_index_vec.push_back(res_orig_index);
					m_char_following_orig_index_vec.push_back(char_following_res_orig_index);
					m_seq_vec.push_back(CBaseFunc::aminoAcidAbbr3WordsTo1(lline.substr(17, 3)));
				}
			} 
		}
		
		m_all_info_vec.push_back(lline);
		m_all_xyz_vec.push_back(xyz);
		
		line++;
	}
}

inline void Molecule::load_one_nucletide(const vector<string>* contents){
	int size = contents->size();
	int line = 0;
	int llen; 
	string lline;
	double x, y, z;
	double* xyz = NULL;
	char cstr[50];
	int res_orig_index;
	char char_following_res_orig_index;
	
	while (size > line){
		lline = (*contents)[line];
		CBaseFunc::toUpperString(lline);
		
		llen = lline.size();
		
		if (llen < 54) {
			line++;
			continue;
		}
				
		strcpy(cstr, (lline.substr(30, 8)).c_str());
		sscanf(cstr, "%lf", &x);
		strcpy(cstr, (lline.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", &y);
		strcpy(cstr, (lline.substr(46, 8)).c_str());
		sscanf(cstr, "%lf", &z);
		
		xyz = new double[3];
		xyz[0] = x;
		xyz[1] = y;
		xyz[2] = z;
		
		if (0 == lline.compare(12, 4, g_atomtype_nuc)){
			strcpy(cstr, (lline.substr(22, 4)).c_str());
			sscanf(cstr, "%d", &res_orig_index);
			char_following_res_orig_index = lline[26];
			
			if (0 == m_orig_index_vec.size() || res_orig_index != m_orig_index_vec[m_orig_index_vec.size()-1] 
					|| char_following_res_orig_index != m_char_following_orig_index_vec[m_char_following_orig_index_vec.size()-1]){
				m_cared_xyz_vec.push_back(xyz);
				m_orig_index_vec.push_back(res_orig_index);
				m_char_following_orig_index_vec.push_back(char_following_res_orig_index);
				m_seq_vec.push_back(CBaseFunc::nucletideAbbr3WordsTo1(CBaseFunc::stringTrim(lline.substr(17, 3))));
			}else{
				cout << "Warning! Duplicated residue " << res_orig_index << "." << endl; 
				
				bool is_not_clash_C3 = true;
				if (m_cared_xyz_vec.size() != 0){
					double dis2 = CBaseFunc::distance2(xyz, m_cared_xyz_vec[m_cared_xyz_vec.size()-1]);
					if (dis2 < 14/*~=3.75*3.75*/)
						is_not_clash_C3 = false;
				}
				
				if (is_not_clash_C3){
					m_cared_xyz_vec.push_back(xyz);
					m_orig_index_vec.push_back(res_orig_index);
					m_char_following_orig_index_vec.push_back(char_following_res_orig_index);
					m_seq_vec.push_back(CBaseFunc::nucletideAbbr3WordsTo1(CBaseFunc::stringTrim(lline.substr(17, 3))));
				}
			} 
		}
		
		m_all_info_vec.push_back(lline);
		m_all_xyz_vec.push_back(xyz);
		
		line++;
	}
}

inline void Molecule::load_one_ligand(const vector<string>* contents){
	int size = contents->size();
	int line = 0;
	int llen; 
	string lline;
	double x, y, z;
	double* xyz = NULL;
	char cstr[50];
	int res_orig_index;
	string atomtype;
	string atomsimpletype;
	
	while (size > line){
		lline = (*contents)[line];
		llen = lline.size();
		
		if (llen < 54) {
			line++;
			continue;
		}
				
		strcpy(cstr, (lline.substr(30, 8)).c_str());
		sscanf(cstr, "%lf", &x);
		strcpy(cstr, (lline.substr(38, 8)).c_str());
		sscanf(cstr, "%lf", &y);
		strcpy(cstr, (lline.substr(46, 8)).c_str());
		sscanf(cstr, "%lf", &z);
		
		xyz = new double[3];
		xyz[0] = x;
		xyz[1] = y;
		xyz[2] = z;
		
		strcpy(cstr, (lline.substr(22, 4)).c_str());
		sscanf(cstr, "%d", &res_orig_index);
		
		atomtype = CBaseFunc::stringTrim(lline.substr(12, 4));
		atomsimpletype = CBaseFunc::stringTrim(lline.substr(76, 2));
		if (atomsimpletype != "H"){
			m_orig_index_vec.push_back(res_orig_index);
			m_char_following_orig_index_vec.push_back(' ');
			
			m_cared_xyz_vec.push_back(xyz);
			CBaseFunc::toUpperString(atomtype);
			m_cared_atomtype_vec.push_back(atomtype);
			CBaseFunc::toUpperString(atomsimpletype);
			m_cared_atomsimpletype_vec.push_back(atomsimpletype);
		}
		
		m_all_info_vec.push_back(lline);
		m_all_xyz_vec.push_back(xyz);
		
		line++;
	}
}

inline const MOLTYPE& Molecule::get_moltype() {
	return m_moltype;
}

inline const bool Molecule::is_same_molytpe(const Molecule& mol) {
	return mol.m_moltype == this->m_moltype;
}

inline const double* Molecule::operator [] (const int& i) {
	return m_cared_xyz_vec[i];
}

inline const int Molecule::size() {
	return m_cared_xyz_vec.size();
}

inline Molecule::~Molecule() {
	int i;
	int n = this->m_all_xyz_vec.size();
	for (i = 0; i < n; i++) {
		if (NULL != m_all_xyz_vec[i]) {
			delete[] m_all_xyz_vec[i];
			m_all_xyz_vec[i] = NULL;
		}
	}
	
	vector<string>().swap(m_all_info_vec);
	vector<double*>().swap(m_all_xyz_vec);
	vector<double*>().swap(m_cared_xyz_vec);
	vector<int>().swap(m_orig_index_vec);
	vector<char>().swap(m_seq_vec);
	vector<string>().swap(this->m_cared_atomtype_vec);
}

inline Complex::Complex(const string& pdb) {
	load(pdb);
}

inline Complex::~Complex(){
	int i = 0;
	int n = m_mols.size();
	for (i = 0; i < n; i++) {
		if (NULL != m_mols[i]) {
			delete m_mols[i];
			m_mols[i] = NULL;
		}
	}
	vector<Molecule*>().swap(m_mols);
}

inline void Complex::load(const string& path){
	if (PDB == g_file_format)
		load_pdb(path);
	else if (CIF == g_file_format) 
		load_cif(path);
	else {
		if (is_cif_func(path)){
			load_cif(path);
		}else{
			load_pdb(path);
		}
	}
}

inline void Complex::load_pdb(const string& pdb){
	int i, n;
	
	map<char, vector<string>* > big_molecule_chain_contents_map;
	map<char, MOLTYPE> big_molecule_chain_moltype;
	vector<char> chain_orders;
	
	map<string, vector<string>* > small_molecule_chain_contents_map;
	vector<string> ind_orders;
	
	map<char, vector<string>* >::iterator it;
	map<string, vector<string>* >::iterator smit;
	vector<string>* chain_contents;
	string line;
	MOLTYPE moltype;
	char chain;
	string hetatm_chain_oindex; 
	string hetatm;
	
	// load all atom lines in pdb
	ifstream fin(pdb.c_str());
	while (getline(fin, line)) {
		//just load the first model
		if (0 == line.compare(0, 3, "END"))
			break;
		
		//load all atom and hetatm lines
		if (0 == line.compare(0, 4, "ATOM")) {
			chain = line[21];
			it = big_molecule_chain_contents_map.find(chain);
			if (it != big_molecule_chain_contents_map.end()){
				(it->second)->push_back(line);
			}else{
				MOLTYPE mt = parse_atom_mol_type(line);
				chain_contents = new vector<string>();
				chain_contents->push_back(line);
				big_molecule_chain_contents_map[chain] = chain_contents;
				big_molecule_chain_moltype[chain] = mt;
				chain_orders.push_back(chain);
			}
		}else if (0 == line.compare(0, 6, "HETATM")) {
			hetatm = line.substr(17, 3);
			if (g_is_load_H2O || hetatm != "HOH"){
				hetatm_chain_oindex = CBaseFunc::stringTrim(hetatm) + "[" + CBaseFunc::stringTrim(line.substr(21, 5)) + "]";
			
				smit = small_molecule_chain_contents_map.find(hetatm_chain_oindex);
				if (smit != small_molecule_chain_contents_map.end()){
					(smit->second)->push_back(line);
				}else{
					chain_contents = new vector<string>();
					chain_contents->push_back(line);
					small_molecule_chain_contents_map[hetatm_chain_oindex] = chain_contents;
					
					ind_orders.push_back(hetatm_chain_oindex);
				}
			}
		}
	}
	fin.close();
	
	n = chain_orders.size();
	for (i = 0; i < n; i++) {
		chain = chain_orders[i];
		chain_contents = big_molecule_chain_contents_map[chain];
		moltype = big_molecule_chain_moltype[chain];
		
		Molecule* mol = new Molecule(moltype, chain_contents);
		m_mols.push_back(mol);
		
		string chstr(1, chain);
		this->m_chains.push_back(chstr);
		
		if (mol->size() >= 4 && ((PROTEIN==moltype && g_is_sup_protein) || (DNA==moltype && g_is_sup_dna) || (RNA==moltype && g_is_sup_rna)))
			this->m_avaid_inds.push_back(i);
	}
	
	int nbigmol = m_mols.size();
	n = ind_orders.size();
	for (i = 0; i < n; i++) {
		hetatm_chain_oindex = ind_orders[i];
		chain_contents = small_molecule_chain_contents_map[hetatm_chain_oindex];
		
		Molecule* mol = new Molecule(LIGAND, chain_contents);
		m_mols.push_back(mol);
		
		this->m_chains.push_back(hetatm_chain_oindex);
		
		if (g_is_sup_ligand && mol->size() >= 4)
			this->m_avaid_inds.push_back(i+nbigmol);
	}
	
	// delete the new objects
	for (map<char, vector<string>* >::iterator it = big_molecule_chain_contents_map.begin(); it != big_molecule_chain_contents_map.end(); it++){
		vector<string>().swap(*(it->second));
		delete (it->second);
	}
	map<char, vector<string>* >().swap(big_molecule_chain_contents_map);
	
	for (map<string, vector<string>* >::iterator it = small_molecule_chain_contents_map.begin(); it != small_molecule_chain_contents_map.end(); it++){
		vector<string>().swap(*(it->second));
		delete (it->second);
	}
	map<string, vector<string>* >().swap(small_molecule_chain_contents_map);
}

inline string Complex::transfer_ATOM_line_from_cif_to_pdb(const string &cifLine){
	//ATOM   1    N  N     . ASP A 1 1   ? 5.873   -9.253  -23.702 1.00 47.90 ? ? ? ? ? ? 467  ASP A N     1 
	//ATOM   1    N N     . GLY A 1 12  ? -26.680 -5.524  12.927  1.00 84.97  ? ? ? ? ? ? 12  GLY A N     1 
	size_t start = cifLine.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    
    std::string line = cifLine.substr(start);
    
    std::vector<std::string> fields;
    std::string field;
    bool inField = false;
    
    for (size_t i = 0; i < line.length(); i++) {
        if (line[i] == ' ' || line[i] == '\t') {
            if (inField) {
                fields.push_back(field);
                field.clear();
                inField = false;
            }
        } else {
            field += line[i];
            inField = true;
        }
    }
    
    if (!field.empty()) {
        fields.push_back(field);
    }
    
    // 检查字段数量
    if (fields.size() < 15) {
        return "ERROR: Line has too few fields: " + cifLine;
    }
    
    // 提取字段 - 根据常见的mmCIF格式
    // 索引可能因具体文件而异，这里使用常见的顺序
    std::string atomNum, atomType, atomName, altLoc, resName, chainId, resSeq, iCode;
    std::string x, y, z, occupancy, tempFactor;
    
    // 尝试不同的字段索引方案
    if (fields.size() >= 26) {
        // 可能是格式1，有更多字段
        atomNum = fields[1];
        atomType = fields[2];
        atomName = fields[3];
        altLoc = fields[4];
        resName = fields[5];
        chainId = fields[6];
        // fields[7] 是 entity_id
        resSeq = fields[8];
        iCode = fields[9];
        x = fields[10];
        y = fields[11];
        z = fields[12];
        occupancy = fields[13];
        tempFactor = fields[14];
    } else if (fields.size() >= 15) {
        // 可能是简化格式
        atomNum = fields[1];
        atomType = fields[2];
        
        // 检查fields[3]是否是合理的原子名称
        if (fields[3].length() <= 4 && 
            (fields[3][0] == 'N' || fields[3][0] == 'C' || fields[3][0] == 'O' || 
             fields[3][0] == 'S' || fields[3][0] == 'H')) {
            atomName = fields[3];
            altLoc = (fields.size() > 4) ? fields[4] : " ";
            resName = (fields.size() > 5) ? fields[5] : "UNK";
            chainId = (fields.size() > 6) ? fields[6] : "A";
            resSeq = (fields.size() > 8) ? fields[8] : "1";
            iCode = (fields.size() > 9) ? fields[9] : " ";
            x = (fields.size() > 10) ? fields[10] : "0.000";
            y = (fields.size() > 11) ? fields[11] : "0.000";
            z = (fields.size() > 12) ? fields[12] : "0.000";
            occupancy = (fields.size() > 13) ? fields[13] : "1.00";
            tempFactor = (fields.size() > 14) ? fields[14] : "0.00";
        }
    }
    
    if (altLoc == "." || altLoc == "?") altLoc = " ";
    if (iCode == "." || iCode == "?") iCode = " ";
    
    std::string formattedAtomName;
    if (atomName.length() == 1) {
        formattedAtomName = " " + atomName + "  ";
    } else if (atomName.length() == 2) {
        // 检查是否元素符号在开头
        if (atomName[0] >= 'A' && atomName[0] <= 'Z' && 
            atomName[1] >= 'a' && atomName[1] <= 'z') {
            formattedAtomName = " " + atomName + " ";
        } else {
            formattedAtomName = " " + atomName + " ";
        }
    } else if (atomName.length() == 3) {
        formattedAtomName = " " + atomName;
    } else if (atomName.length() == 4) {
        formattedAtomName = atomName;
    } else {
        formattedAtomName = " " + atomName.substr(0, 3);
    }
    
    // 格式化残基名称
    std::string formattedResName = resName;
    if (formattedResName.length() == 1) {
        formattedResName = formattedResName + "  ";
    } else if (formattedResName.length() == 2) {
        formattedResName = formattedResName + " ";
    }
    
    // 构建PDB格式行
    std::ostringstream pdbStream;
    pdbStream << "ATOM  ";
    
    // 原子序号 (1-5列)
    pdbStream << std::setw(5) << std::right << atomNum;
    pdbStream << " ";
    
    // 原子名称 (13-16列)
    pdbStream << formattedAtomName;
    
    // 交替位置标识符 (17列)
    pdbStream << altLoc[0];
    
    // 残基名称 (18-20列)
    pdbStream << formattedResName;
    
    // 空格 (21列)
    pdbStream << " ";
    
    // 链标识符 (22列)
    pdbStream << chainId[0];
    
    // 残基序列号 (23-26列)
    pdbStream << std::setw(4) << std::right << resSeq;
    
    // 插入码 (27列)
    pdbStream << iCode[0];
    
    // 3个空格 (28-30列)
    pdbStream << "   ";
    
    // 坐标 (31-54列)
    try {
        double xVal = std::stod(x);
        double yVal = std::stod(y);
        double zVal = std::stod(z);
        
        pdbStream << std::fixed << std::setprecision(3);
        pdbStream << std::setw(8) << std::right << xVal;
        pdbStream << std::setw(8) << std::right << yVal;
        pdbStream << std::setw(8) << std::right << zVal;
    } catch (...) {
        pdbStream << std::setw(8) << "0.000" 
                  << std::setw(8) << "0.000" 
                  << std::setw(8) << "0.000";
    }
    
    // 占有率和温度因子 (55-66列)
    try {
        double occVal = std::stod(occupancy);
        double tempVal = std::stod(tempFactor);
        
        pdbStream << std::fixed << std::setprecision(2);
        pdbStream << std::setw(6) << std::right << occVal;
        pdbStream << std::setw(6) << std::right << tempVal;
    } catch (...) {
        pdbStream << std::setw(6) << "1.00" 
                  << std::setw(6) << "0.00";
    }
    
    // 10个空格 (67-76列)
    pdbStream << "          ";
    
    // 元素符号 (77-78列)
    if (atomType.length() >= 2) {
        pdbStream << atomType.substr(0, 2);
    } else {
        pdbStream << std::setw(2) << std::right << atomType;
    }
    
    return pdbStream.str();
}

inline string Complex::transfer_HETATM_line_from_cif_to_pdb(const string &cifLine){
	size_t start = cifLine.find_first_not_of(" \t");
    if (start == std::string::npos) return "";
    
    std::string line = cifLine.substr(start);
    
    std::vector<std::string> fields;
    std::string field;
    bool inField = false;
    
    for (size_t i = 0; i < line.length(); i++) {
        if (line[i] == ' ' || line[i] == '\t') {
            if (inField) {
                fields.push_back(field);
                field.clear();
                inField = false;
            }
        } else {
            field += line[i];
            inField = true;
        }
    }
    
    if (!field.empty()) {
        fields.push_back(field);
    }
    
    // 检查字段数量
    if (fields.size() < 15) {
        return "ERROR: Line has too few fields: " + cifLine;
    }
    
    // 提取字段 - 根据常见的mmCIF格式
    // 索引可能因具体文件而异，这里使用常见的顺序
    std::string atomNum, atomType, atomName, altLoc, resName, chainId, resSeq, iCode;
    std::string x, y, z, occupancy, tempFactor;
    
    // 尝试不同的字段索引方案
    if (fields.size() >= 26) {
        // 可能是格式1，有更多字段
        atomNum = fields[1];
        atomType = fields[2];
        atomName = fields[3];
        altLoc = fields[4];
        resName = fields[5];
        chainId = fields[6];
        // fields[7] 是 entity_id
        resSeq = fields[8];
        iCode = fields[9];
        x = fields[10];
        y = fields[11];
        z = fields[12];
        occupancy = fields[13];
        tempFactor = fields[14];
    } else if (fields.size() >= 15) {
        // 可能是简化格式
        atomNum = fields[1];
        atomType = fields[2];
        
        // 检查fields[3]是否是合理的原子名称
        if (fields[3].length() <= 4 && 
            (fields[3][0] == 'N' || fields[3][0] == 'C' || fields[3][0] == 'O' || 
             fields[3][0] == 'S' || fields[3][0] == 'H')) {
            atomName = fields[3];
            altLoc = (fields.size() > 4) ? fields[4] : " ";
            resName = (fields.size() > 5) ? fields[5] : "UNK";
            chainId = (fields.size() > 6) ? fields[6] : "A";
            resSeq = (fields.size() > 8) ? fields[8] : "1";
            iCode = (fields.size() > 9) ? fields[9] : " ";
            x = (fields.size() > 10) ? fields[10] : "0.000";
            y = (fields.size() > 11) ? fields[11] : "0.000";
            z = (fields.size() > 12) ? fields[12] : "0.000";
            occupancy = (fields.size() > 13) ? fields[13] : "1.00";
            tempFactor = (fields.size() > 14) ? fields[14] : "0.00";
        }
    }
    
    if (altLoc == "." || altLoc == "?") altLoc = " ";
    if (iCode == "." || iCode == "?") iCode = " ";
    
    std::string formattedAtomName;
    if (atomName.length() == 1) {
        formattedAtomName = " " + atomName + "  ";
    } else if (atomName.length() == 2) {
        // 检查是否元素符号在开头
        if (atomName[0] >= 'A' && atomName[0] <= 'Z' && 
            atomName[1] >= 'a' && atomName[1] <= 'z') {
            formattedAtomName = " " + atomName + " ";
        } else {
            formattedAtomName = " " + atomName + " ";
        }
    } else if (atomName.length() == 3) {
        formattedAtomName = " " + atomName;
    } else if (atomName.length() == 4) {
        formattedAtomName = atomName;
    } else {
        formattedAtomName = " " + atomName.substr(0, 3);
    }
    
    // 格式化残基名称
    std::string formattedResName = resName;
    if (formattedResName.length() == 1) {
        formattedResName = formattedResName + "  ";
    } else if (formattedResName.length() == 2) {
        formattedResName = formattedResName + " ";
    }
    
    // 构建PDB格式行
    std::ostringstream pdbStream;
    pdbStream << "HETATM";
    
    // 原子序号 (1-5列)
    pdbStream << std::setw(5) << std::right << atomNum;
    pdbStream << " ";
    
    // 原子名称 (13-16列)
    pdbStream << formattedAtomName;
    
    // 交替位置标识符 (17列)
    pdbStream << altLoc[0];
    
    // 残基名称 (18-20列)
    pdbStream << formattedResName;
    
    // 空格 (21列)
    pdbStream << " ";
    
    // 链标识符 (22列)
    pdbStream << chainId[0];
    
    // 残基序列号 (23-26列)
    pdbStream << std::setw(4) << std::right << resSeq;
    
    // 插入码 (27列)
    pdbStream << iCode[0];
    
    // 3个空格 (28-30列)
    pdbStream << "   ";
    
    // 坐标 (31-54列)
    try {
        double xVal = std::stod(x);
        double yVal = std::stod(y);
        double zVal = std::stod(z);
        
        pdbStream << std::fixed << std::setprecision(3);
        pdbStream << std::setw(8) << std::right << xVal;
        pdbStream << std::setw(8) << std::right << yVal;
        pdbStream << std::setw(8) << std::right << zVal;
    } catch (...) {
        pdbStream << std::setw(8) << "0.000" 
                  << std::setw(8) << "0.000" 
                  << std::setw(8) << "0.000";
    }
    
    // 占有率和温度因子 (55-66列)
    try {
        double occVal = std::stod(occupancy);
        double tempVal = std::stod(tempFactor);
        
        pdbStream << std::fixed << std::setprecision(2);
        pdbStream << std::setw(6) << std::right << occVal;
        pdbStream << std::setw(6) << std::right << tempVal;
    } catch (...) {
        pdbStream << std::setw(6) << "1.00" 
                  << std::setw(6) << "0.00";
    }
    
    // 10个空格 (67-76列)
    pdbStream << "          ";
    
    // 元素符号 (77-78列)
    if (atomType.length() >= 2) {
        pdbStream << atomType.substr(0, 2);
    } else {
        pdbStream << std::setw(2) << std::right << atomType;
    }
    
    return pdbStream.str();
}

inline bool Complex::is_cif_func(const string& path){
	bool ans = false;
	string line; 
	ifstream fin(path.c_str());
	while (getline(fin, line)) {
		//just load the first model
		if (0 == line.compare(0, 14, "_audit_conform") 
			|| 0 == line.compare(0, 16, "_struct.entry_id")
			|| 0 == line.compare(0, 10, "_atom_site")){
			ans = true;
			break;
		}
	}
	fin.close();
	return ans; 
}

inline void Complex::load_cif(const string& cif){
	int i, n;
	
	map<char, vector<string>* > big_molecule_chain_contents_map;
	map<char, MOLTYPE> big_molecule_chain_moltype;
	vector<char> chain_orders;
	
	map<string, vector<string>* > small_molecule_chain_contents_map;
	vector<string> ind_orders;
	
	map<char, vector<string>* >::iterator it;
	map<string, vector<string>* >::iterator smit;
	vector<string>* chain_contents;
	string line;
	MOLTYPE moltype;
	char chain;
	string hetatm_chain_oindex; 
	string hetatm;
	
	// load all atom lines in cif
	ifstream fin(cif.c_str());
	bool is_atom_start = false;
	while (getline(fin, line)) {
		//just load the first model
		if (is_atom_start && 0 == line.compare(0, 1, "#"))
			break;
		
		//load all atom and hetatm lines
		if (0 == line.compare(0, 4, "ATOM")) {
			is_atom_start = true;
			line = transfer_ATOM_line_from_cif_to_pdb(line);
			chain = line[21];
			it = big_molecule_chain_contents_map.find(chain);
			if (it != big_molecule_chain_contents_map.end()){
				(it->second)->push_back(line);
			}else{
				MOLTYPE mt = parse_atom_mol_type(line);
				chain_contents = new vector<string>();
				chain_contents->push_back(line);
				big_molecule_chain_contents_map[chain] = chain_contents;
				big_molecule_chain_moltype[chain] = mt;
				chain_orders.push_back(chain);
			}
		}else if (0 == line.compare(0, 6, "HETATM")) {
			is_atom_start = true;
			line = transfer_HETATM_line_from_cif_to_pdb(line); 
			
			hetatm = line.substr(17, 3);
			if (g_is_load_H2O || hetatm != "HOH"){
				hetatm_chain_oindex = CBaseFunc::stringTrim(hetatm) + "[" + CBaseFunc::stringTrim(line.substr(21, 5)) + "]";
			
				smit = small_molecule_chain_contents_map.find(hetatm_chain_oindex);
				if (smit != small_molecule_chain_contents_map.end()){
					(smit->second)->push_back(line);
				}else{
					chain_contents = new vector<string>();
					chain_contents->push_back(line);
					small_molecule_chain_contents_map[hetatm_chain_oindex] = chain_contents;
					
					ind_orders.push_back(hetatm_chain_oindex);
				}
			}
		}
	}
	fin.close();
	
	n = chain_orders.size();
	for (i = 0; i < n; i++) {
		chain = chain_orders[i];
		chain_contents = big_molecule_chain_contents_map[chain];
		moltype = big_molecule_chain_moltype[chain];
		
		Molecule* mol = new Molecule(moltype, chain_contents);
		m_mols.push_back(mol);
		
		string chstr(1, chain);
		this->m_chains.push_back(chstr);
		
		if (mol->size() >= 4 && ((PROTEIN==moltype && g_is_sup_protein) || (DNA==moltype && g_is_sup_dna) || (RNA==moltype && g_is_sup_rna)))
			this->m_avaid_inds.push_back(i);
	}
	
	int nbigmol = m_mols.size();
	n = ind_orders.size();
	for (i = 0; i < n; i++) {
		hetatm_chain_oindex = ind_orders[i];
		chain_contents = small_molecule_chain_contents_map[hetatm_chain_oindex];
		
		Molecule* mol = new Molecule(LIGAND, chain_contents);
		m_mols.push_back(mol);
		
		this->m_chains.push_back(hetatm_chain_oindex);
		
		if (g_is_sup_ligand && mol->size() >= 4)
			this->m_avaid_inds.push_back(i+nbigmol);
	}
	
	// delete the new objects
	for (map<char, vector<string>* >::iterator it = big_molecule_chain_contents_map.begin(); it != big_molecule_chain_contents_map.end(); it++){
		vector<string>().swap(*(it->second));
		delete (it->second);
	}
	map<char, vector<string>* >().swap(big_molecule_chain_contents_map);
	
	for (map<string, vector<string>* >::iterator it = small_molecule_chain_contents_map.begin(); it != small_molecule_chain_contents_map.end(); it++){
		vector<string>().swap(*(it->second));
		delete (it->second);
	}
	map<string, vector<string>* >().swap(small_molecule_chain_contents_map);
}

inline const MOLTYPE Complex::parse_atom_mol_type(const string& line){
	string threeWordAA = CBaseFunc::stringTrim(line.substr(17, 3));
	int len = threeWordAA.size();
	if (3 == len){
		return PROTEIN;
	}else{
		if ('D' == threeWordAA[0])
			return DNA;
		else return RNA;
	}
}

inline const string Molecule::get_seq_str(){
	return this->m_seq_str;
}

inline const vector<char>& Molecule::get_seq_vec(){
	return m_seq_vec;
}

inline double* CBaseFunc::rotateAndTrans(double* xyz, double** u, double* t){
	double x = xyz[0];
	double y = xyz[1];
	double z = xyz[2];
	double xt = t[0] + u[0][0]*x + u[0][1]*y + u[0][2]*z;
	double yt = t[1] + u[1][0]*x + u[1][1]*y + u[1][2]*z;
	double zt = t[2] + u[2][0]*x + u[2][1]*y + u[2][2]*z;
	
	double* ans = new double[3];
	ans[0] = xt;
	ans[1] = yt;
	ans[2] = zt;
	
	return ans;
}

inline Molecule* Complex::operator [](const int& i){
	return m_mols[m_avaid_inds[i]];
}

inline const int Complex::size(){
	return m_avaid_inds.size();
}

inline const int Complex::all_size(){
	return m_mols.size();
}

inline Molecule* Complex::get_ith_mol_in_all(const int& i){
	return m_mols[i];
}

inline const vector<string>& Molecule::to_str(){
	return m_all_info_vec;
}

inline const vector<string> Molecule::to_str(double** u, double* t){
	int i = 0;
	double* oxyz = NULL;
	double* nxyz = NULL;
	string oinfo;
	char xstr[9];
	char ystr[9];
	char zstr[9];
	
	vector<string> ans; 
	
	const int N = m_all_info_vec.size();
	for (i = 0; i < N; i++){
		oinfo = m_all_info_vec[i];
		oxyz = m_all_xyz_vec[i];
		nxyz = CBaseFunc::rotateAndTrans(oxyz, u, t);

		sprintf(xstr, "%8.3f", nxyz[0]);
		sprintf(ystr, "%8.3f", nxyz[1]);
		sprintf(zstr, "%8.3f", nxyz[2]);
		
		ans.push_back(oinfo.substr(0, 30) + xstr + ystr + zstr + oinfo.substr(54));
		delete[] nxyz;
	}
	
	return ans;
}

inline const vector<string> Complex::to_str(){
	vector<string> ans;
	int i, n = m_mols.size();
	for (i = 0; i < n; i++){
		vector<string> ith = m_mols[i]->to_str();
		ans.insert(ans.end(), ith.begin(), ith.end());
		ans.push_back("TER");
	}
	
	return ans;
}

inline const vector<string> Complex::to_str(double** u, double* t){
	vector<string> ans;
	int i, n = m_mols.size();
	for (i = 0; i < n; i++){
		vector<string> ith = m_mols[i]->to_str(u, t);
		ans.insert(ans.end(), ith.begin(), ith.end());
		ans.push_back("TER");
	}
	
	return ans;
}

inline void Complex::save(const string& path){
	ofstream fout(path.c_str());
	vector<string> cc = to_str();
	int i, n = cc.size();
	for (i = 0; i < n; i++){
		fout << cc[i] << endl;
	}
	fout.close(); 
}

inline void Complex::save(const string& path, double** u, double* t){
	ofstream fout(path.c_str());
	vector<string> cc = to_str(u, t);
	int i, n = cc.size();
	for (i = 0; i < n; i++){
		fout << cc[i] << endl;
	}
	fout.close(); 
}

inline bool CBaseFunc::is_same(const vector<char>& a, const vector<char>& b){
	int i, n = a.size();
	bool ans = true;
	if (n != b.size())
		ans = false;
	else{
		for (i = 0; i < n; i++){
			if (a[i] != b[i]){
				ans = false;
				break;				
			}
		}	
	}
	
	return ans;
}

inline const vector<double*>& Molecule::get_cared_xyz_vec(){
	return this->m_cared_xyz_vec;
}

inline const int& CNWalign::get_identity_num(){
	return identity_num;
}

inline double CBaseFunc::cal_rot_tran_from_query_to_templ__(const vector<double*>& query, const vector<double*>& templ, double** out_u, double* out_t, const double& user_d0, const int* q2t, const bool& fast){
	int i, j, L = query.size();
	
	vector<double*> pse_qxyz_vec;
	vector<double*> pse_txyz_vec;
	for (i = 0; i < L; i++){
		j = q2t[i];
		if (-1 != j){
			pse_qxyz_vec.push_back(query[i]);
			pse_txyz_vec.push_back(templ[j]);
		}
	}
	
	return CBaseFunc::cal_rot_tran_from_query_to_templ__(pse_qxyz_vec, pse_txyz_vec, out_u, out_t, user_d0, fast); 
}

inline const int& Molecule::get_ith_orig_index(const int& i){
	return this->m_orig_index_vec[i];
}

inline const char& Molecule::get_ith_char_following_orig_index_vec(const int& i){
	return this->m_char_following_orig_index_vec[i];
}

inline const vector<int>& Molecule::get_orig_index_vec(){
	return this->m_orig_index_vec;
}

inline const double& CTMscoreComplex::get_tmscore(){
	return this->tmscore;
}

inline double** CTMscoreComplex::get_u(){
	return this->u;
}

inline const double* CTMscoreComplex::get_t(){
	return this->t;
}

inline double CBaseFunc::distance2(const double* a, const double* b){
	double xbias = a[0] - b[0];
	double ybias = a[1] - b[1];
	double zbias = a[2] - b[2];
	return xbias*xbias + ybias*ybias + zbias*zbias;
}

inline double CBaseFunc::distance2(const double* q, const double* t, double** qu, double* qt){
	double x = q[0];
	double y = q[1];
	double z = q[2];
	double xt = qt[0] + qu[0][0]*x + qu[0][1]*y + qu[0][2]*z;
	double yt = qt[1] + qu[1][0]*x + qu[1][1]*y + qu[1][2]*z;
	double zt = qt[2] + qu[2][0]*x + qu[2][1]*y + qu[2][2]*z;
	
	double xbias = xt - t[0];
	double ybias = yt - t[1];
	double zbias = zt - t[2];
	return xbias*xbias + ybias*ybias + zbias*zbias;
}

inline double CBaseFunc::rough_score(const MOLTYPE& mt, const vector<double*>& q, const vector<double*>& t, const int* q2t){
	int i;
	int qsize = q.size();
	int tsize = t.size();
	int size = qsize>tsize ? qsize : tsize;
	double d0, dis2;
	
	if (mt == PROTEIN){
		d0 = CBaseFunc::d0_of_tmscore(size);
	}else if (mt != LIGAND){
		d0 = CBaseFunc::d0_of_lsscore(size);
	}else{
		d0 = CBaseFunc::d0_of_tmscore_c3prime(size);
	}
	double d02 = d0*d0;
	
	double tmsco = 0.;
	for (i = 0; i < qsize; i++){
		if (-1 != q2t[i]){
			dis2 = CBaseFunc::distance2(q[i], t[q2t[i]]);
			tmsco += 1./(1. + dis2/d02);
		}
	}
	
	return tmsco / size;
}

inline void CBaseFunc::u3b(const vector<double*>& query, const vector<double*>& templ, double** out_u, double* out_t){
	int n = query.size();
	double** x = new2Darr(3+1, n+1);
	double** y = new2Darr(3+1, n+1);
	double** u = new2Darr(3+1, 3+1);
	double* t = new double[3+1];
	int mode = 1;
	int ier = 0;
	
	for (int i = 0; i < n; i++) {
		x[1][i+1] = query[i][0];
		x[2][i+1] = query[i][1];
		x[3][i+1] = query[i][2];
		
		y[1][i+1] = templ[i][0];
		y[2][i+1] = templ[i][1];
		y[3][i+1] = templ[i][2];
	}
	
	u3b(x, y, n, mode, u, t);
	
	for (int i = 1; i < 4; i++){
		for (int j = 1; j < 4; j++){
			out_u[i-1][j-1] = u[i][j];
		}
		out_t[i-1] = t[i];
	}
	
	delete2Darr(x, 3+1);
	delete2Darr(y, 3+1);
	delete2Darr(u, 3+1);
	delete[] t;
}

inline double CBaseFunc::u3b_func(const vector<double*>& query, const vector<double*>& templ, const vector<MOLTYPE>& mts, const double& d02_pro, const double& d02_dna, const double& d02_rna, const double& d02_lig){
	int n = query.size();
	double** x = new2Darr(3+1, n+1);
	double** y = new2Darr(3+1, n+1);
	double** u = new2Darr(3+1, 3+1);
	double* t = new double[3+1];
	double** u33 = new2Darr(3, 3);
	double* t3 = new double[3];
	int mode = 1;
	int ier = 0;
	
	for (int i = 0; i < n; i++) {
		x[1][i+1] = query[i][0];
		x[2][i+1] = query[i][1];
		x[3][i+1] = query[i][2];
		
		y[1][i+1] = templ[i][0];
		y[2][i+1] = templ[i][1];
		y[3][i+1] = templ[i][2];
	}
	
	u3b(x, y, n, mode, u, t);
	
	for (int i = 1; i < 4; i++){
		for (int j = 1; j < 4; j++){
			u33[i-1][j-1] = u[i][j];
		}
		t3[i-1] = t[i];
	} 
	
	double tmscore = 0.;
	for (int i = 0; i < n; i++) {
		double* newx = rotateAndTrans(query[i], u33, t3);
		double dis2 = distance2(newx, templ[i]);
		const MOLTYPE& mt = mts[i];	
		
		double d02 = 0.;
		if (mt == PROTEIN){
			d02 = d02_pro;
		}else if (mt == DNA){
			d02 = d02_dna;
		}else if (mt == RNA){
			d02 = d02_rna;
		}else{
			d02 = d02_lig;
		}
		
		tmscore += 1./(1. + dis2/d02);
		
		delete[] newx;
	}
	
	delete2Darr(u33, 3);
	delete2Darr(x, 3+1);
	delete2Darr(y, 3+1);
	delete2Darr(u, 3+1);
	delete[] t;
	delete[] t3;
	
	return tmscore;
}

inline const string& Complex::get_chain(const int& i){
	return this->m_chains[m_avaid_inds[i]];
}

inline Molecule* Complex::get_mol_base_on_chain(const string& chain){
	int size = m_chains.size();
	for (int i = 0; i < size; i++){
		if (chain == m_chains[i])
			return m_mols[i];
	}
	return NULL; 
}

inline void CBaseFunc::inverseRotAndTransMtx(double** input_rot_mtx, double* input_trans, double** output_rot_mtx, double* output_trans) {
	inverse33Mtx(input_rot_mtx, output_rot_mtx);
	
	output_trans[0] = -(input_trans[0]*output_rot_mtx[0][0] + input_trans[1]*output_rot_mtx[0][1] + input_trans[2]*output_rot_mtx[0][2]);
	output_trans[1] = -(input_trans[0]*output_rot_mtx[1][0] + input_trans[1]*output_rot_mtx[1][1] + input_trans[2]*output_rot_mtx[1][2]);
	output_trans[2] = -(input_trans[0]*output_rot_mtx[2][0] + input_trans[1]*output_rot_mtx[2][1] + input_trans[2]*output_rot_mtx[2][2]);
}

inline double CBaseFunc::yuzishiForInverse33Mtx(double** A, int i, int j){
	double M[2][2];
	for (int r = 0, rr = 0; r < 3; r++){
		if (r == i) continue;
		for (int c = 0, cc = 0; c < 3; c++){
			if (c == j)continue;
			M[rr][cc] = A[r][c];
			cc++;
		}
		rr++;
	}

	double ans = M[0][0]*M[1][1] - M[0][1]*M[1][0];
	
	if (1 == (i+j)%2){
		return -ans;
	}else{
		return ans;
	}
}

inline bool CBaseFunc::inverse33Mtx(double** A, double** invA){	
	double det = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][1] + A[0][2]*A[1][0]*A[2][1] 
		         - A[2][0]*A[1][1]*A[0][2] - A[1][0]*A[0][1]*A[2][2] - A[0][0]*A[2][1]*A[1][2];

	if (0.00000001 > fabs(det))
		return false;

	invA[0][0] = yuzishiForInverse33Mtx(A, 0, 0)/det; invA[0][1] = yuzishiForInverse33Mtx(A, 1, 0)/det; invA[0][2] = yuzishiForInverse33Mtx(A, 2, 0)/det; 
	invA[1][0] = yuzishiForInverse33Mtx(A, 0, 1)/det; invA[1][1] = yuzishiForInverse33Mtx(A, 1, 1)/det; invA[1][2] = yuzishiForInverse33Mtx(A, 2, 1)/det; 
	invA[2][0] = yuzishiForInverse33Mtx(A, 0, 2)/det; invA[2][1] = yuzishiForInverse33Mtx(A, 1, 2)/det; invA[2][2] = yuzishiForInverse33Mtx(A, 2, 2)/det;
	
	double norm = 0.0;
	for (int i = 0; i < 3; i++){
		norm += A[0][i]*invA[i][0] + A[1][i]*invA[i][1] + A[2][i]*invA[i][2];
	}
	norm /= 3.0;

	invA[0][0] /= norm; invA[0][1] /= norm; invA[0][2] /= norm; 
	invA[1][0] /= norm; invA[1][1] /= norm; invA[1][2] /= norm; 
	invA[2][0] /= norm; invA[2][1] /= norm; invA[2][2] /= norm;

	return true;
}

inline void CTMscoreComplex::save_invUt(const string& savepath){
	ofstream fout(savepath.c_str());
	vector<string> inv_ut_vec = formatInvRotMtx();
	int i, n = inv_ut_vec.size();
	for (i = 0; i < n; i++)
		fout << inv_ut_vec[i] << endl;
	fout.close();
}

inline void CTMscoreComplex::save_superposition_ditances(const string& savepath){
	ofstream fout(savepath.c_str());
	fout << "Detail superposition information: " << endl;
	fout << "----------------------------------------------------------------------------" << endl;
	fout << "               Structure 1                                Structure 2       " << endl;
	fout << "        --------------------------                --------------------------" << endl;
	fout << "ind.    mol.          oind.   ant*   distance     mol.          oind.    ant" << endl;  
	fout << "----------------------------------------------------------------------------" << endl;
	char buf[200];
	int i, j, k=0, n = this->aa_level_ali.size();
	for (j = 0; j < tsize; j++){
		string tchain = templ->get_chain(j);
		for (i = 0; i < n; i++){
			ALIGN_PAIR ap = this->aa_level_ali[i];
			if (0 != strcmp(tchain.c_str(), ap.tchain.c_str()))
				continue;
			
			double dis = sqrt(ap.dis2);
			
			if (dis < 10.){
				sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.6f    %10s    %5d%c   %3s", 
						k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
			}else if (dis < 100){
				sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.5f    %10s    %5d%c   %3s", 
						k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
			}else if (dis < 1000){
				sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.4f    %10s    %5d%c   %3s", 
						k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
			}else if (dis < 10000){
				sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.3f    %10s    %5d%c   %3s", 
						k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
			}else if (dis < 100000){
				sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.2f    %10s    %5d%c   %3s", 
						k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
			}else if (dis < 1000000){
				sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.1f    %10s    %5d%c   %3s", 
						k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
			}
			
			fout << buf << endl;
		}
	}	
	fout << "------------------------------------------------------------------" << endl;
	fout << "* 'ant' means amino acid type, nucletide type, or atom type, resp-" << endl;
	fout << "  ectively."<<endl; 
	fout << "------------------------------------------------------------------" << endl;
	fout.close();
}

inline void CTMscoreComplex::print_result(){
	if (CBaseFunc::isInit2Darr(this->individual_tmscore_mtx, this->qsize, this->tsize)){
		cout << "There is an issue for \"individual_tmscore_mtx\", please send the two complex structures you inputted to hj@ism.cams.cn." << endl;
		cout << "We will be very grateful to you!!!" << endl;
		exit(1);
	} 
	
	int i, j, k, n;
	
	CBaseFunc::print_logo();
	
	int* inv_obj_level_ali = new int[tsize];
	for (i = 0; i < tsize; i++)
		inv_obj_level_ali[i] = -1;
	for (i = 0; i < qsize; i++){
		if (-1 != obj_level_ali[i]){
			inv_obj_level_ali[obj_level_ali[i]] = i;
		}
	}
	
	bool is_all_qchain_matched = true;
	for (i = 0; i < qsize; i++){
		if (-1 == obj_level_ali[i]){
			is_all_qchain_matched = false;
			break;	
		}
	}
	bool is_all_tchain_matched = qsize==tsize ? is_all_qchain_matched : false;
	
	if (!is_all_qchain_matched || !is_all_tchain_matched) {
		int q_unmatched_num = 0;
		for (i = 0; i < qsize; i++){
			if (-1 == obj_level_ali[i]){
				q_unmatched_num++;
			}
		}
		int t_unmatched_num = 0;
		for (i = 0; i < tsize; i++){
			if (-1 == inv_obj_level_ali[i]){
				t_unmatched_num++;
			}
		}
		cout << endl;
		cout << "*****************************************************************************************" << endl;
		cout << "WARNING:" << endl; 
		cout << "-->> THE CHAINS/MOLECULES OF TWO INPUTTED STRUCTURES CAN NOT MATCH EACH OTHER PERFECTLY." << endl;
		if (0 != t_unmatched_num) {
			if (1 == t_unmatched_num){
				cout << "-->> UNMATCHED CHAIN/MOLECULE IN Structure 1: ";
			}else {
				cout << "-->> UNMATCHED CHAINS/MOLECULES IN Structure 1: ";
			}
			for (i = 0; i < tsize; i++){
				if (-1 == inv_obj_level_ali[i]){
					cout << this->templ->get_chain(i) << ' ';
				}
			}
			cout << endl;	
		}
		if (0 != q_unmatched_num) {
			if (1 == q_unmatched_num){
				cout << "-->> UNMATCHED CHAIN/MOLECULE IN Structure 2: ";
			}else {
				cout << "-->> UNMATCHED CHAINS/MOLECULES IN Structure 2: ";
			}
			for (i = 0; i < qsize; i++){
				if (-1 == obj_level_ali[i]){
					cout << this->query->get_chain(i) << ' ';
				}
			}
			cout << endl;	
		}
		cout << "*****************************************************************************************" << endl;
		cout << endl;
		cout << endl;
	} 
		
	int qpro_num = 0;
	int qdna_num = 0;
	int qrna_num = 0;
	int qlig_num = 0;
	int tpro_num = 0;
	int tdna_num = 0;
	int trna_num = 0;
	int tlig_num = 0;
	vector<string> qpro_chains;
	vector<string> qdna_chains;
	vector<string> qrna_chains;
	vector<string> qlig_types;
	vector<string> tpro_chains;
	vector<string> tdna_chains;
	vector<string> trna_chains;
	vector<string> tlig_types;
	vector<int> qpro_aanums;
	vector<int> qdna_nucnums;
	vector<int> qrna_nucnums;
	vector<int> qlig_atnums;
	vector<int> tpro_aanums;
	vector<int> tdna_nucnums;
	vector<int> trna_nucnums;
	vector<int> tlig_atnums;
	int qprot_aanum_in_tot = 0;
	int qdna_ntnum_in_tot = 0;
	int qrna_ntnum_in_tot = 0;
	int qlig_atnum_in_tot = 0;
	int tprot_aanum_in_tot = 0;
	int tdna_ntnum_in_tot = 0;
	int trna_ntnum_in_tot = 0;
	int tlig_atnum_in_tot = 0;
	
	n = this->query->size();
	for (i = 0; i < n; i++){
		Molecule* mol = (*(this->query))[i];
		MOLTYPE mt = mol->get_moltype();
		switch (mt){
			case PROTEIN:
				qpro_num++;
				qpro_chains.push_back(this->query->get_chain(i));
				qpro_aanums.push_back(mol->size());
				qprot_aanum_in_tot += mol->size();
				break;
			case DNA:
				qdna_num++;
				qdna_chains.push_back(this->query->get_chain(i));
				qdna_nucnums.push_back(mol->size());
				qdna_ntnum_in_tot += mol->size();
				break;
			case RNA:
				qrna_num++;
				qrna_chains.push_back(this->query->get_chain(i));
				qrna_nucnums.push_back(mol->size());
				qrna_ntnum_in_tot += mol->size();
				break;
			case LIGAND:
				qlig_num++;
				qlig_types.push_back(this->query->get_chain(i));
				qlig_atnums.push_back(mol->size());
				qlig_atnum_in_tot += mol->size();
				break; 
		}
	}
	
	n = this->templ->size();
	for (i = 0; i < n; i++){
		Molecule* mol = (*(this->templ))[i];
		MOLTYPE mt = mol->get_moltype();
		switch (mt){
			case PROTEIN:
				tpro_num++;
				tpro_chains.push_back(this->templ->get_chain(i));
				tpro_aanums.push_back(mol->size());
				tprot_aanum_in_tot += mol->size();
				break;
			case DNA:
				tdna_num++;
				tdna_chains.push_back(this->templ->get_chain(i));
				tdna_nucnums.push_back(mol->size());
				tdna_ntnum_in_tot += mol->size();
				break;
			case RNA:
				trna_num++;
				trna_chains.push_back(this->templ->get_chain(i));
				trna_nucnums.push_back(mol->size());
				trna_ntnum_in_tot += mol->size();
				break;
			case LIGAND:
				tlig_num++;
				tlig_types.push_back(this->templ->get_chain(i));
				tlig_atnums.push_back(mol->size());
				tlig_atnum_in_tot += mol->size(); 
				break; 
		}
	}
	
	char buf[5000];
	cout << "Information of Structure 1: " << endl;
	if (0 != tpro_num){
		sprintf(buf, "  >Protein number: %d, amino acid number: %d, protein chain/molecule(s): ", tpro_num, tprot_aanum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << tpro_chains[0] << " (" << tpro_aanums[0] << " a.a.)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = tpro_chains.size();
		for (i = 1; i < n; i++){
			cout << buf << tpro_chains[i] << " (" << tpro_aanums[i] << " a.a.)" << endl;
		}
	}
	if (0 != tdna_num){
		sprintf(buf, "  >DNA number: %d, base pair number: %d, DNA chain/molecule(s): ", tdna_num, tdna_ntnum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << tdna_chains[0] << " (" << tdna_nucnums[0] << " bp)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = tdna_chains.size();
		for (i = 1; i < n; i++){
			cout << buf << tdna_chains[i] << " (" << tdna_nucnums[i] << " bp)" << endl;
		}
	} 
	if (0 != trna_num){
		sprintf(buf, "  >RNA number: %d, nucletide number: %d, RNA chain/molecule(s): ", trna_num, trna_ntnum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << trna_chains[0] << " (" << trna_nucnums[0] << " nt)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = trna_chains.size();
		for (i = 1; i < n; i++){
			cout << buf << trna_chains[i] << " (" << trna_nucnums[i] << " nt)" << endl;
		}
	}
	if (0 != tlig_num){
		sprintf(buf, "  >Ligand number: %d, heavy atom number: %d, Ligand molecules: ", tlig_num, tlig_atnum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << tlig_types[0] << " (" << tlig_atnums[0] << " atoms)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = tlig_types.size();
		for (i = 1; i < n; i++){
			cout << buf << tlig_types[i] << " (" << tlig_atnums[i] << " atoms)" << endl;
		}
	}
	cout << endl;
	cout << "Information of Structure 2: " << endl;
	if (0 != qpro_num){
		sprintf(buf, "  >Protein number: %d, amino acid number: %d, protein chain/molecule(s): ", qpro_num, qprot_aanum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << qpro_chains[0] << " (" << qpro_aanums[0] << " a.a.)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = qpro_chains.size();
		for (i = 1; i < n; i++){
			cout << buf << qpro_chains[i] << " (" << qpro_aanums[i] << " a.a.)" << endl;
		}
	} 
	if (0 != qdna_num){
		sprintf(buf, "  >DNA number: %d, base pair number: %d, DNA chain/molecule(s): ", qdna_num, qdna_ntnum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << qdna_chains[0] << " (" << qdna_nucnums[0] << " bp)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = qdna_chains.size();
		for (i = 1; i < n; i++){
			cout << buf << qdna_chains[i] << " (" << qdna_nucnums[i] << " bp)" << endl;
		}
	} 
	if (0 != qrna_num){
		sprintf(buf, "  >RNA number: %d, nucletide number: %d, RNA chain/molecule(s): ", qrna_num, qrna_ntnum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << qrna_chains[0] << " (" << qrna_nucnums[0] << " nt)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = qrna_chains.size();
		for (i = 1; i < n; i++){
			cout << buf << qrna_chains[i] << " (" << qrna_nucnums[i] << " nt)" << endl;
		}
	}
	if (0 != qlig_num){
		sprintf(buf, "  >Ligand number: %d, heavy atom number: %d, Ligand molecule(s): ", qlig_num, qlig_atnum_in_tot);
		int spacenum = strlen(buf);
		cout << buf << qlig_types[0] << " (" << qlig_atnums[0] << " atoms)" << endl;
		for (i = 0; i < spacenum; i++)
			buf[i] = ' ';
		buf[i] = '\0';
		
		n = qlig_types.size();
		for (i = 1; i < n; i++){
			cout << buf << qlig_types[i] << " (" << qlig_atnums[i] << " atoms)" << endl;
		}
	}
	cout << endl << endl;
	
	cout << "Molecule mapping information: " << endl;
	cout << "  >Molecule in Structure 1: ";
	for (i = 0; i < tsize; i++){
		if (-1 != inv_obj_level_ali[i]){
			cout << this->templ->get_chain(i) << ' ';
		}
	}
	cout << endl;
	cout << "  >Molecule in Structure 2: ";
	for (i = 0; i < tsize; i++){
		if (-1 != inv_obj_level_ali[i]){
			cout << this->query->get_chain(inv_obj_level_ali[i]) << ' ';
		}
	}
	cout << endl;
	
	map<string, double> chain_tmscore;
	map<string, double> chain_d02;
	
	for (i = 0; i < qsize; i++){
		if (-1 != obj_level_ali[i]){
			double d02 = this->chain_index_corr_to_query__d02[i];
			chain_d02[query->get_chain(i)] = d02;
			chain_tmscore[query->get_chain(i)] = 0.;
		}
	}
	n = this->aa_level_ali.size();
	for (i = 0; i < n; i++){
		ALIGN_PAIR ap = this->aa_level_ali[i];
		double dis2 = ap.dis2;
		double d02 = chain_d02[ap.qchain];
		chain_tmscore[ap.qchain] += 1./(1. + dis2 / d02);
	}
	cout << endl << endl;
	
	cout << "TM-score of each pair of mapped molecules calculated using individual rotation matrix" << endl;
	for (i = 0; i < tsize; i++){
		if (-1 != inv_obj_level_ali[i]){
			double individual_tmsco = this->individual_tmscore_mtx[inv_obj_level_ali[i]][i];
			sprintf(buf, "  >TM-score between %s(in Structure 1) and %s(in Structure 2): %8.6f", templ->get_chain(i).c_str(), query->get_chain(inv_obj_level_ali[i]).c_str(), individual_tmsco);
			cout << buf << endl; 
		}
	}
	cout << endl << endl;
	
	cout << "TM-score of each pair of mapped molecules calculated using complex rotation matrix" << endl;
	for (i = 0; i < tsize; i++){
		if (-1 != inv_obj_level_ali[i]){
			double tmsco = chain_tmscore[query->get_chain(inv_obj_level_ali[i])];
			int num = this->chain_index_corr_to_query__aa_num[inv_obj_level_ali[i]];
			sprintf(buf, "  >TM-score between %s(in Structure 1) and %s(in Structure 2): %8.6f", templ->get_chain(i).c_str(), query->get_chain(inv_obj_level_ali[i]).c_str(), tmsco/num);
			cout << buf << endl;
			
			chain_tmscore[query->get_chain(inv_obj_level_ali[i])] = tmsco/num;
		}
	}
	cout << endl << endl;
	
	cout << "rTM-score(n): calculated using the highest n TM-scores in the previous region" << endl;
	n = chain_tmscore.size();
	vector<string> chains;
	double* tmscore_arr = new double[n];
	map<string, double>::iterator it;
	for  (i = 0, it = chain_tmscore.begin(); it != chain_tmscore.end(); i++, it++){
		chains.push_back(it->first);
		tmscore_arr[i] = it->second;
	}
	int* sinds = CSort::quickDescendSortIndex(n, tmscore_arr, 0, n-1);
	int max_rTM_level = n;
	double corr_rTM_of_max_rTM_level = 0.;
	for (i = 1; i <= n; i++){
		double rTM = 0.;
		for (j = 0; j < i; j++){
			rTM += 1./ tmscore_arr[sinds[j]];
		}
		rTM = i/rTM;
		
		corr_rTM_of_max_rTM_level = rTM;
		sprintf(buf, "  >rTM-score(%d): %8.6f", i, rTM);
		cout << buf << endl; 
	}
//	n = templ->size(); // n = query->size();
//	for (; i <= n; i++){
//		sprintf(buf, "  >rTM-score(%d): %8.6f (donot matched )", i, 0.);
//		cout << buf << endl; 
//	}
	delete[] sinds;
	delete[] tmscore_arr;
	cout << endl << endl;
	
	cout << "Final metrics for complex structures of Structure 1 and Structure 2:" << endl;
	double itmscore = calcluate_itmscore(g_interact_dis_cut_pow2);
	if (0 > itmscore){
		if (is_all_qchain_matched && is_all_tchain_matched)
			sprintf(buf, "  >TM-score: %8.6f, rTM-score: %8.6f, iTM-score: n/a", tmscore, rtmscore);
		else sprintf(buf, "  >TM-score: %8.6f, rTM-score: n/a, iTM-score: n/a", tmscore);
		cout << buf << endl;
	}else {
		if (is_all_qchain_matched && is_all_tchain_matched){
			sprintf(buf, "  >TM-score: %8.6f, rTM-score: %8.6f, iTM-score: %8.6f", tmscore, rtmscore, itmscore);
			cout << buf << endl;
		}else{
			sprintf(buf, "  >TM-score: %8.6f, rTM-score: n/a, iTM-score: %8.6f", tmscore, itmscore);
			cout << buf << endl;
		} 
	}
	cout << endl << endl;
	
	cout << "Complex Rotation Matrix (rotate Structure 1 to Structure 2):" << endl;
	vector<string> inv_ut_vec = formatInvRotMtx();
	n = inv_ut_vec.size();
	for (i = 0; i < n; i++)
		cout << inv_ut_vec[i] << endl;
	
	if (g_is_output_in_detail){
		cout << endl << endl;
		cout << "Detail superposition information: " << endl;
	
		if (0 == qlig_num || 0 == tlig_num){
			cout << "------------------------------------------------------------------" << endl;
			cout << "             Structure 1                          Structure 2     " << endl;
			cout << "        ---------------------                ---------------------" << endl;
			cout << "ind.    mol.     oind.   ant*   distance     mol.     oind.    ant" << endl;  
			cout << "------------------------------------------------------------------" << endl;
		}else{
			cout << "----------------------------------------------------------------------------" << endl;
			cout << "              Structure 1                                Structure 2        " << endl;
			cout << "        --------------------------                --------------------------" << endl;
			cout << "ind.    mol.          oind.   ant*   distance     mol.          oind.    ant" << endl;  
			cout << "----------------------------------------------------------------------------" << endl;
		}
		
		n = this->aa_level_ali.size();
		k = 0;
		for (j = 0; j < tsize; j++){
			string tchain = templ->get_chain(j);
			for (i = 0; i < n; i++){
				ALIGN_PAIR ap = this->aa_level_ali[i];
				if (0 != strcmp(tchain.c_str(), ap.tchain.c_str()))
					continue;
				
				double dis = sqrt(ap.dis2);
				if (0 == qlig_num || 0 == tlig_num){
					if (dis < 10.){
						sprintf(buf, "%4d    %5s    %5d%c   %3s    %8.6f    %5s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
					}else if (dis < 100){
						sprintf(buf, "%4d    %5s    %5d%c   %3s    %8.5f    %5s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());	
					}else if (dis < 1000){
						sprintf(buf, "%4d    %5s    %5d%c   %3s    %8.4f    %5s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());
					}else if (dis < 10000){
						sprintf(buf, "%4d    %5s    %5d%c   %3s    %8.3f    %5s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());
					}else if (dis < 100000){
						sprintf(buf, "%4d    %5s    %5d%c   %3s    %8.2f    %5s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());
					}else if (dis < 1000000){
						sprintf(buf, "%4d    %5s    %5d%c   %3s    %8.1f    %5s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());
					}
				}else{
					if (dis < 10.){
						sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.6f    %10s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
					}else if (dis < 100){
						sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.5f    %10s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
					}else if (dis < 1000){
						sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.4f    %10s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
					}else if (dis < 10000){
						sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.3f    %10s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
					}else if (dis < 100000){
						sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.2f    %10s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
					}else if (dis < 1000000){
						sprintf(buf, "%4d    %10s    %5d%c   %3s    %8.1f    %10s    %5d%c   %3s", 
								k+1, ap.tchain.c_str(), ap.toind, ap.toindsuf, ap.taa.c_str(), dis, ap.qchain.c_str(), ap.qoind, ap.qoindsuf, ap.qaa.c_str());			
					}
				}
				
				cout << buf << endl;
				k++;
			}	
		}
		
		if (0 == qlig_num || 0 == tlig_num){
			cout << "------------------------------------------------------------------" << endl;
			cout << "* 'ant' means amino acid type, nucletide type, or atom type, resp-" << endl;
			cout << "  ectively."<<endl; 
			cout << "------------------------------------------------------------------" << endl;
		}else if  (0 != qlig_num && 0 != tlig_num){
			cout << "----------------------------------------------------------------------------" << endl;
			cout << "* 'ant' means amino acid type, nucletide type, or atom type, respectively." << endl;
			cout << "----------------------------------------------------------------------------" << endl;
		}
	}
	cout << endl << endl;
	cout << "Taking " << setprecision(8) << use_seconds << " seconds in total." << endl;
	
	
	delete[] inv_obj_level_ali;
}


inline vector<string> CTMscoreComplex::formatRotMtx(){
	vector<string> ans;
	
	char buf[1000];
	sprintf(buf, "i%18s %15s %15s %15s", "t[i]", "u[i][0]", "u[i][1]", "u[i][2]");
	ans.push_back(string(buf));
	for (int k = 0; k < 3; k++){
		sprintf(buf, "%d%18.10f %15.10f %15.10f %15.10f", k, t[k], u[k][0], u[k][1], u[k][2]);
		ans.push_back(string(buf));
	}
//	ans.push_back(string("\nCode for rotating Structure 2 from (x,y,z) to (X,Y,Z):"));
//	ans.push_back(string("for(k=0; k<L; k++)"));
//	ans.push_back(string("{"));
//	ans.push_back(string("   X[k] = t[0] + u[0][0]*x[k] + u[0][1]*y[k] + u[0][2]*z[k]"));
//	ans.push_back(string("   Y[k] = t[1] + u[1][0]*x[k] + u[1][1]*y[k] + u[1][2]*z[k]"));
//	ans.push_back(string("   Z[k] = t[2] + u[2][0]*x[k] + u[2][1]*y[k] + u[2][2]*z[k]"));
//	ans.push_back(string("}"));

	return ans;
}

inline vector<string> CTMscoreComplex::formatInvRotMtx(){
	invUt();
	
	vector<string> ans;
	
	char buf[1000];
	sprintf(buf, "i%18s %15s %15s %15s", "t[i]", "u[i][0]", "u[i][1]", "u[i][2]");
	ans.push_back(string(buf));
	for (int k = 0; k < 3; k++){
		sprintf(buf, "%d%18.10f %15.10f %15.10f %15.10f", k, inv_t[k], inv_u[k][0], inv_u[k][1], inv_u[k][2]);
		ans.push_back(string(buf));
	}
//	ans.push_back(string("\nCode for rotating Structure 1 from (x,y,z) to (X,Y,Z):"));
//	ans.push_back(string("for(k=0; k<L; k++)"));
//	ans.push_back(string("{"));
//	ans.push_back(string("   X[k] = t[0] + u[0][0]*x[k] + u[0][1]*y[k] + u[0][2]*z[k]"));
//	ans.push_back(string("   Y[k] = t[1] + u[1][0]*x[k] + u[1][1]*y[k] + u[1][2]*z[k]"));
//	ans.push_back(string("   Z[k] = t[2] + u[2][0]*x[k] + u[2][1]*y[k] + u[2][2]*z[k]"));
//	ans.push_back(string("}"));

	return ans;
}

inline void CTMscoreComplex::save_roted_query(const string& savepath){
	query->save(savepath, u, t);
}

inline void CTMscoreComplex::invUt(){
	if (NULL == inv_u){
		inv_u = CBaseFunc::new2Darr(3, 3);
		inv_t = new double[3];
		CBaseFunc::inverseRotAndTransMtx(u, t, inv_u, inv_t);
	}
}

inline void CTMscoreComplex::save_roted_templ(const string& savepath){
	invUt();
	templ->save(savepath, inv_u, inv_t);
}


inline void CTMscoreComplex::save_sup_pdb_for_web(const string& savepath){
	int i, n;
	
	vector<string> q  = query->to_str();
	invUt();
	vector<string> rt = templ->to_str(inv_u, inv_t);
	
	char buf[200];
	
	ofstream fout(savepath.c_str());
	n = q.size();
	for (i = 0; i < n && i < 50000; i++){
		const string& c = q[i];
		if (0 == c.compare(0, 4, "ATOM")
			|| 0 == c.compare(0, 6, "HETATM")){
			if (c.length() > 22){
				sprintf(buf, "%5d", i);
				fout << c.substr(0, 6) << buf << c.substr(11, 10) << 'N' << c.substr(22) << endl;
			}else fout << c << endl; 			
		}
	}
	
	n = rt.size();
	for (i = 0; i < n && i < 50000; i++){
		const string& c = rt[i];
		if (0 == c.compare(0, 4, "ATOM")
			|| 0 == c.compare(0, 6, "HETATM")){
			if (c.length() > 22){
				sprintf(buf, "%5d", i+50000);
				fout << c.substr(0, 6) << buf << c.substr(11, 10) << 'M' << c.substr(22) << endl;
			}else fout << c << endl;
		}
	}
	
	fout.close();
}


inline const vector<CHAIN_PAIR>* CBaseFunc::parse_matched_chain_pairs(const string& matched_chain_pair_path){
	string line; 
	// load all atom lines in pdb
	ifstream fin(matched_chain_pair_path.c_str());
	vector<CHAIN_PAIR>* ans = new vector<CHAIN_PAIR>();
	while (getline(fin, line)) {
		// ignore the comment lines
		if (0 == line.compare(0, 1, "#"))
			continue;
		line = CBaseFunc::stringTrim(line);
		if (2 > line.size())
			continue;
		
		vector<string> items = CBaseFunc::stringSplit(line, ' ', '\t');
		if (2 != items.size()){
			cout << "ERROR : The format of line \""<< line << "\" in \"" << matched_chain_pair_path << "\" is not correct" << endl;
			exit(1);
		}
		CHAIN_PAIR p;
		p.qchain = items[0];
		p.tchain = items[1];
		
		ans->push_back(p);
	}
	fin.close();
	
	return ans;
}

inline vector<string> CBaseFunc::string2stringVec(const string& str){
	vector<string> ans;
	int n = str.size();
	for (int i = 0; i < n; i++){
		ans.push_back(str.substr(i, i+1));
	}
	
	return ans;
}

inline void CBaseFunc::toUpperString(string &str){  
   transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper);  
}

inline string CBaseFunc::eraseAll(const string &str, char ch){
    string s(str.c_str());
    int index = 0;
    if( !s.empty())
    {
        while( (index = s.find(ch, index)) != string::npos)
        {
            s.erase(index,1);
        }
    }
		
		return s;
}

inline string CBaseFunc::eraseAll(const string &str, const char* arr, int len){
    string s(str.c_str());
    for (int i = 0; i < len; i++){
		s = eraseAll(s, arr[i]);
	}
	return s;
}

inline int*** CBaseFunc::new3DIntArr(int row, int col, int thd){
	int ***ans=new int**[row];
	for(int i=0;i<row;i++){
		ans[i]=new int*[col];
		for (int j=0; j<col;j++){
			ans[i][j]=new int[thd];
			for (int k=0; k<thd; k++){
				ans[i][j][k] = 0;
			}
		}
	}
	
	return ans;
}

inline const string& Molecule::get_cared_atomsimpletype_in_lig(const int& i){
	return m_cared_atomsimpletype_vec[i];
}

inline const vector<string>& Molecule::get_cared_atomtype_vec_in_lig(){
	return m_cared_atomtype_vec;
}

inline void CBaseFunc::delete3DIntArr(int row, int col, int*** Arr){
	for(int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			delete[] Arr[i][j];
		}
		delete[] Arr[i];
	}
	delete[] Arr;
	Arr = NULL;
}

inline const bool LigAtomMatcher::is_same_important(const int& i, const LigAtomMatcher& other, const int& other_i){
	vector<vector<string>* >& ith = *(attypes_in_roads[i]);
	vector<vector<string>* >& other_ith = *(other.attypes_in_roads[other_i]);
	
//	cout << "DEBUG START" << endl;
//	int imax = 0;
//	for (int i = 0; i < ith.size(); i++){
//		if (imax < (*ith[i]).size())
//			imax = (*ith[i]).size(); 
//	}
//	int omax = 0;
//	for (int i = 0; i < other_ith.size(); i++){
//		if (omax < (*other_ith[i]).size())
//			omax = (*other_ith[i]).size(); 
//	}
//	if (imax == omax){
//		for (int i = 0; i < ith.size(); i++){
//			for (int j = 0; j < (*ith[i]).size(); j++){
//				cout << (*ith[i])[j] << ' ';
//			}
//			cout << endl;
//		}
//		cout << "======" << endl;
//		for (int i = 0; i < other_ith.size(); i++){
//			for (int j = 0; j < (*other_ith[i]).size(); j++){
//				cout << (*other_ith[i])[j] << ' ';
//			}
//			cout << endl;
//		}
//	}
//	cout << "DEBUG END" << endl;
	
	return CVector::is_same_contents(ith, other_ith, g_lig_match_in_order_or_content);
}

inline void LigAtomMatcher::extract_atgs(){
	int i, k, n, m;
	bool is_new_group;
	for (i = 0; i < atnum; i++){
		is_new_group = true;
		n = atgs.size();
		for (k = 0; k < n; k++){
			vector<int>& kth_atg = *atgs[k];
			if (is_same_important(kth_atg[0], i)){
				kth_atg.push_back(i);
				is_new_group = false;
				break;
			}
		}
		
		if (is_new_group){
			vector<int>* p_atg = new vector<int>();
			p_atg->push_back(i);
			atgs.push_back(p_atg);	
		}
	}
}

LigAtomMatcher::LigAtomMatcher(Molecule* p_lig):lig(*p_lig){
	atnum = lig.size();
	is_same_imp_mtx = CBaseFunc::new2DBoolArr(atnum, atnum);
	
	string* attypes = new string[atnum];
	double* vdwr = new double[atnum];
	double** adjmtx = CBaseFunc::new2Darr(atnum, atnum);
	
	AtomVarDerWaalRadius vdwrobj;
	for (int i = 0; i < atnum; i++){
		attypes[i] = p_lig->get_cared_atomsimpletype_in_lig(i);
		vdwr[i] = vdwrobj[attypes[i]];
	}
	
	for (int i = 0; i < atnum; i++){
		for (int j = i+1; j < atnum; j++){
			double vdwij = (vdwr[i]>vdwr[j] ? vdwr[i] : vdwr[j])+0.06;
			double dis2 = CBaseFunc::distance2(lig[i], lig[j]);
			if (dis2 < vdwij*vdwij)
				adjmtx[i][j] = adjmtx[j][i] = 0.5*(vdwr[i]+vdwr[j]);
		}
	}
	
	ShortestRoad sr(adjmtx, atnum);
	int*** all_roads = sr.getAllRoad8Floyd();
	for (int i = 0; i < atnum; i++)
		attypes_in_roads.push_back(extract_ith_roads(all_roads[i], attypes));
	
	for (int i = 0; i < atnum; i++){
		vector<vector<string>* >& iroads = *attypes_in_roads[i];
		for (int j = i+1; j < atnum; j++){
			vector<vector<string>* >& jroads = *attypes_in_roads[j];
			is_same_imp_mtx[i][j] = is_same_imp_mtx[j][i] = is_same_roads(iroads, jroads);
		}
	}
	
	extract_atgs();
	
	delete[] attypes;
	delete[] vdwr;
	CBaseFunc::delete2Darr(adjmtx, atnum);
}

LigAtomMatcher:: ~LigAtomMatcher(){
	int n = atgs.size();
	for (int i = 0; i < n; i++){
		if (NULL != atgs[i])
			delete atgs[i];
	}
	
	n = attypes_in_roads.size();
	for (int i = 0; i < n; i++){
		vector<vector<string>* >& ith = *attypes_in_roads[i];
		if (NULL != attypes_in_roads[i]){
			int m = ith.size();
			for (int j = 0; j < m; j++){
				if (NULL != ith[j]){
					delete ith[j];
				}
			}
			delete attypes_in_roads[i];
		}
		
	}
	
	CBaseFunc::delete2DBoolArr(atnum, is_same_imp_mtx);
}

inline bool LigAtomMatcher::is_same_roads(const vector<vector<string>* >& arvec, const vector<vector<string>* >& brvec){
	return CVector::is_same_contents(arvec, brvec, g_lig_match_in_order_or_content);
}

inline vector<vector<string>* >* LigAtomMatcher::extract_ith_roads(int** ar, const string* attypes){
	vector<vector<string>* >* arvec = new vector<vector<string>* >();
	for (int i = 0; i < atnum; i++){
		vector<string>* road = new vector<string>();
		for (int k = 0; k < atnum; k++){
			if (-1 == ar[i][k]) break;
			road->push_back(attypes[ar[i][k]]);
		}
		
		if (0 < road->size()){
			arvec->push_back(road);	
		} else delete road;
	}
	
	if (0 >= arvec->size()){
		delete[] arvec;
		return NULL;
	}else{
		return arvec;
	}
}

inline const bool& LigAtomMatcher::is_same_important(const int& i, const int& j){
	return is_same_imp_mtx[i][j];
}

inline double CBaseFunc::gdt_ts(const vector<double>& dis2_vec){
	const double cutoffs_pow2[] = {1.*1., 2.*2., 4.*4., 8.*8.};
	const int cutoffs_len = 4;
	int qlen = dis2_vec.size();
	double GDT_Px[cutoffs_len];
	for (int k = 0; k < cutoffs_len; k++) {
		GDT_Px[k] = 0.;
		for (int i = 0; i < qlen; i++) {
			if (dis2_vec[i] <= cutoffs_pow2[k])
				GDT_Px[k]++;
		}
	}
	
	double ans = 0.;
	for (int k = 0; k < cutoffs_len; k++) {
		ans += GDT_Px[k]/qlen*100;
	}
	
	return ans/cutoffs_len;
}

inline double CBaseFunc::gdt_ha(const vector<double>& dis2_vec){
	const double cutoffs_pow2[] = {0.5*0.5, 1.*1., 2.*2., 4.*4.};
	const int cutoffs_len = 4;
	int qlen = dis2_vec.size();
	double GDT_Px[cutoffs_len];
	for (int k = 0; k < cutoffs_len; k++) {
		GDT_Px[k] = 0.;
		for (int i = 0; i < qlen; i++) {
			if (dis2_vec[i] <= cutoffs_pow2[k])
				GDT_Px[k]++;
		}
	}
	
	double ans = 0.;
	for (int k = 0; k < cutoffs_len; k++) {
		ans += GDT_Px[k]/qlen*100;
	}
	
	return ans/cutoffs_len;
}
inline double CBaseFunc::rmsd(const vector<double>& dis2_vec){
	double ans = 0.;
	int qlen = dis2_vec.size();
	for (int i = 0; i < qlen; i++) {
		ans += dis2_vec[i];
	}
	ans /= qlen;
	return sqrt(ans);
}

inline void CBaseFunc::print_help(const char* arg){
	CBaseFunc::print_logo();
	
	cout << "TM-score-Comp (TM-scoreC) is a quick and accurate algorithm for measuring quality of " << endl
	     << "complex structure predictions of proteins, nucleic acids, and small molecule ligands." << endl
	     << endl;
	cout << " Usage 1: " << arg << " structure_1.pdb structure_2.pdb [Options]" << endl
	     << " Usage 2: " << arg << " structure_1.cif structure_2.cif [Options]" << endl << endl
		 << " Options:" << endl
		 << "   -mol      Select molecule types in the inputted pdb files to superimpose. "<< endl 
		 << "               all (default): superimpose all molecules in inputs." << endl
		 << "               prt          : superimpose proteins in inputs only." << endl
		 << "               dna          : superimpose DNAs in inputs only." << endl
		 << "               rna          : superimpose RNAs in inputs only." << endl
		 << "               lig          : superimpose small molecule ligands in inputs only." << endl
		 << "               p+d          : superimpose proteins and DNAs in inputs." << endl
		 << "               p+r          : superimpose proteins and RNAs in inputs." << endl
		 << "               p+l          : superimpose proteins and ligands in inputs." << endl
		 << "               d+r          : superimpose DNAs and RNAs in inputs." << endl
		 << "               d+l          : superimpose DNAs and ligands in inputs." << endl
		 << "               r+l          : superimpose RNAs and ligands in inputs." << endl
		 << "               pdr          : superimpose proteins, DNAs and RNA in inputs." << endl
		 << "               pdl          : superimpose proteins, DNAs and ligands in inputs." << endl
		 << "               prl          : superimpose proteins, RNAs and ligands in inputs." << endl
		 << "               drl          : superimpose DNAs, RNAs and ligands in inputs." << endl
		 << "   -s        Select TM-score or rTM-score to search the best superposition." << endl
		 << "               t (default and strongly suggested): use TM-score " << endl
		 << "               r                                 : use rTM-score " << endl
		 << "   -d0       The TM-score scaled by an assigned d0, e.g., '-d0 3.5' reports MaxSub" << endl
		 << "             score, where d0 is 3.5 Angstrom." << endl
		 << "   -da       Is molecule order information in two inputted pdb files matched? Note " << endl
		 << "             that, only one of the options of '-ia', and '-da y' can be applied at " << endl 
		 << "             the same time." << endl
		 << "               y          : using molecule order in inputted files to generate molecule mapping" << endl
		 << "               n (default): genrating molecule mapping automatically" << endl
		 << "   -ffm      Setting the file formats of two inputted structures:" << endl
		 << "               pdb           : the formats of two inputted files are both 'PDB' format (see " << endl
		 << "                               https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)." << endl
		 << "               mmcif         : the formats of two inputted files are both 'mmCIF' format (see " << endl
		 << "                               https://mmcif.wwpdb.org/docs/tutorials/content/atomic-description.html)." << endl
		 << "               auto (default): the formats of two inputted files are different, this program will check " << endl
		 << "                               them to be one of 'PDB' and 'mmCIF' automatically." << endl
		 << "   -ia       Input molecule mapping txt file path (see detail in README). Note that" << endl
		 << "             only one of the options of '-ia', and '-da y' can be applied at once" << endl 
		 << "   -ri       Are residue/nucletide/atom indexes of the homologous protein/DNA(RNA)/ligand" << endl
		 << "             molecules in the inputted pdb files matched?" << endl
		 << "               y          : using the matched residue/nucleotide/atom indexes to generate" << endl
		 << "                            the residue/nucletide/atom assignment directly for all homolous " << endl
		 << "                            protein/DNA(RNA)/ligand molecules." << endl
 		 << "               n (default): using global sequence aligment tool, i.e., NW-align, to generate" << endl
 		 << "                            residue/nucleotide alignment in homologous protein/DNA(RNA) molecules" << endl 
 		 << "                            once and using the atom mapping algorithm (see Manuscript) to generated" << endl
		 << "                            the atom mapping for ligands." << endl 
		 << "   -sid      When option of \"-ri\" is set to 'n', this option can be selected to " << endl
		 << "             set the global sequence identity cutoff (**default: 0.7**, range from 0 to" << endl
		 << "             1), which is used to ensure wheter two molecule sequence is homologous " << endl
		 << "             or not. When NWalign-outputted global sequence identity value is larger" << endl
		 << "             than the inputted cutoff, the corresponding sequences will be looked as " << endl 
		 << "             homologous with each other. In TM-score-Comp, only homologous molecules "<< endl
		 << "             can matched each other." << endl
		 << "   -wt       Whether load water molecules in inputted pdb files or not?" << endl
		 << "               y          : load water molecules in inputted pdb files." << endl
 		 << "               n (default): do not load water molecules in inputted pdb files." << endl
 		 << "   -odis     Whether  directly output the distance information in detail between each "<<endl
		 << "             pair of the mapped elements of residues for protein, nucletides for DNAs "<<endl
		 <<"              and RNAs and atoms for ligands in the command line window." << endl
 		 << "               y          : output the distances of mapped residues/nucletides/atoms in detial" << endl
 		 << "               n (default): donot output the distances of mapped residues/nucletides/atoms." << endl
		 << "   -mode     Modes of TM-score-Comp, which can imfact the speed and accuracy of the final superimposed " << endl
		 << "             results." << endl 
		 << "               fast            : in average, the speed and accuracy are faster and slightly lower than " << endl
		 << "                                 those of 'normal', respectively" << endl 
		 << "               normal (default): the normal model of TM-score-Comp" << endl
		 << "   -atom-nuc 4-character atom name used to represent a nucletide. Default is \" C3'\" for " << endl
		 << "             DNA/RNA (note the space before C3')." << endl
		 << "             PLEASE MAKE SURE THAT THE ATOM TYPE EXISTS IN EVERY NUCLETIDE." << endl
		 << "   -atom-res 4-character atom name used to represent a residue. Default is \" CA \" for " << endl
		 << "             protein (note the spaces before and after CA)." << endl
		 << "             PLEASE MAKE SURE THAT THE ATOM TYPE EXISTS IN EVERY AMINO ACID RESIDUE." << endl
		 << "   -nit      Set maximum iteration number in the original TM-score program," << endl
		 << "             range from 20 to positive infinity, default is 20." << endl
		 << "             ** NOTE THAT, THIS OPTION MAY INCREACE THE RUNING TIME, PLEASE MAKE SURE" << endl
		 << "             THAT YOU WANT MORE ITERATION TIMES, WHICH MAY GIVE YOU A HIGH SCORE VALUE." << endl
		 << "   -nLinit   Set maximum number of L_init in the original TM-score program," << endl 
		 << "             range from 6 to positive infinity, default is 6." << endl
		 << "             ** NOTE THAT, THIS OPTION MAY INCREACE THE RUNING TIME, PLEASE MAKE SURE" << endl
		 << "             THAT YOU WANT MORE L_init NUMBER, WHICH MAY GIVE YOU A HIGH SCORE VALUE." << endl
		 << "   -clig     Whether re-mapping ligand atom pairs in every scoring time or just at begining." << endl 
		 << "             Note that, this option is just worked on \"-mode normal\", and cannot worked with " << endl
		 << "             options of \"-ia\", \"-da y\", and \"-ri y\"" << endl
		 << "               y (defualt): re-mapping ligand atom pair in every scoring time." << endl
		 << "               n          : re-mapping ligand atom pair just at begining." << endl
		 << "   -o        Save the Structure 1 with PDB format after superposition to file." << endl
		 << "   -srm      Save the complex rotation matrix to file." << endl
		 << "   -ssp      Save the detail superposition information to file." << endl
		 << "   -h        print this help" << endl << endl
		 << " Example usages:" << endl
		 << "    "<< arg <<" predicted.pdb native.pdb" << endl
		 << "    "<< arg <<" predicted.cif native.cif" << endl
		 << "    "<< arg <<" predicted.pdb native.pdb -mol all" << endl
		 << "    "<< arg <<" predicted.pdb native.pdb -da y" << endl
		 << "    "<< arg <<" predicted.pdb native.pdb -ia chain_align.txt" << endl
		 << "    "<< arg <<" predicted.pdb native.pdb -ri y" << endl
		 << "    "<< arg <<" predicted.pdb native.pdb -sid 0.9" << endl
		 << "    "<< arg <<" predicted.pdb native.pdb -o rotatted_predicted.pdb" << endl
		 << "    "<< arg <<" -h"<< endl
		 << "    "<< arg <<" -v" << endl << endl;
	exit(1);
}

//inline void CBaseFunc::print_logo(){
//	cout 
//	<<"================================================================================" << endl
//	<<" ______,__ __                            ___                       _                 " << endl
//	<<"(_) | /|  |  |                          / (_)                     | |                " << endl
//	<<"    |  |  |  |   ,   __   __   ,_    _ |      __   _  _  _     _  | |  _             " << endl
//	<<"  _ |  |  |  |  / \\_/    /  \\_/  |  |/ |     /  \\_/ |/ |/ |  |/ \\_|/  |/  /\\/   " << endl
//	<<" (_/   |  |  |_/ \\/ \\___/\\__/    |_/|__/\\___/\\__/   |  |  |_/|__/ |__/|__/ /\\_/" << endl
//	<<"                                                            /|                       " << endl
//	<<"                                                            \\|                      " << endl
//	<<"Version of TM-score-Complex (TMscoreC): " << VERSION << endl
//	<<"Please email comments and suggestions to Jun Hu (hj@ism.cams.cn)" << endl
//	<<"================================================================================" << endl << endl;
//}

inline void CBaseFunc::print_logo(){
	cout 
	<<"========================================================================================" << endl
	<<"             ______,__ __                            ___                      " << endl
	<<"            (_) | /|  |  |                          / (_)                     " << endl
	<<"                |  |  |  |   ,   __   __   ,_    _ |      __   _  _  _     _  " << endl
	<<"              _ |  |  |  |  / \\_/    /  \\_/  |  |/ |     /  \\_/ |/ |/ |  |/ \\" << endl
	<<"             (_/   |  |  |_/ \\/ \\___/\\__/    |_/|__/\\___/\\__/   |  |  |_/|__/ " << endl
	<<"                                                                        /|                       " << endl
	<<"                                                                        \\|                      " << endl
	<<"Version of TM-score-Comp (TMscoreC): " << VERSION << endl
	<<"Please email comments and suggestions to Jun Hu (hj@ism.cams.cn)" << endl
	<<"========================================================================================" << endl << endl;
}
