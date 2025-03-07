/*
 * globals.hpp
 *
 *  Created on: June 20, 2012
 *      Author: hwchoi
 */

#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <cmath>

#include <boost/foreach.hpp>
#include <boost/array.hpp>
#include "PreyClass.hpp"
#include "BaitClass.hpp"
#include "InterClass.hpp"
#include "UIClass.hpp"
typedef unsigned Count_t;
// typedef double Quant_t;
class Quant_t;
namespace saint {
	inline bool isnan(const Quant_t);
}
class Quant_t {
	friend Quant_t exp(Quant_t);
	friend bool saint::isnan(const Quant_t x);
 	double quant;
public:
 	Quant_t(double x) : quant(x){ }
 	Quant_t() : quant(0){ }
 	inline Quant_t& operator= (double a) {
		quant = a;
		return *this;
	}
 	inline operator double& () {
 		// if(std::isnan(quant)) throw runtime_error("invalid");
 		return quant;
 	}
 	inline operator const double& () const {
 		// if(std::isnan(quant)) throw runtime_error("invalid const");
 		return quant;
 	}
	inline bool operator<(Quant_t r) const {
		// if(std::isnan(quant) || std::isnan(r.quant)) cout<<"comparing with NAN"<<flush;
		return quant < r.quant;
	}
};
#include "Fastmat.hpp"
#include "Stats_fns.hpp"


using namespace std;

const unsigned iter_Z_d = 10;

namespace saint {
	const double nan = std::numeric_limits<double>::quiet_NaN();
	// inline bool isnan(const double x) { return std::isnan(x); }
	inline bool isnan(const Quant_t x) { return std::isnan(x.quant); }
}
inline Quant_t exp(Quant_t a){
	if(saint::isnan(a)) return saint::nan;
	return std::exp(a);
}
vector<string> splitString(const string & input_txt);
void getFileDimensions(string inputFile, int &nr, int &nc);

deque<InterClass> parseInterFile(const string & inputFile, size_t &ninter);

void mapRowCol( deque<InterClass> &IDATA, const deque<PreyClass> &PDATA, const deque<BaitClass> &BDATA, const map<string, BaitClass>&);

std::set<string> parsePreyFile( deque<PreyClass> &PDATA, const string & inputFile, size_t &nprey);

map<string, BaitClass> parseBaitFile( deque<BaitClass> &BDATA, const string & inputFile, size_t &nip);
void sortBaitData( deque<BaitClass> &BDATA );

void createList( deque<UIClass> &UIDATA, const deque<InterClass> &IDATA, const deque<BaitClass> &BDATA, size_t &nuinter );

int get_nexpr( deque<BaitClass> &BDATA );
int get_nctrl( const deque<BaitClass> &BDATA );

inline std::string detect_line_ending(const string& file_name) {
	std::ifstream infile(file_name);
	char tmp;
	for(int i = 0;i < 1e3; i++){
		infile.get(tmp);
		if(tmp == '\r') {		// Mac or Windows
			infile.get(tmp);
			if(tmp == '\n' )	// Windows
				return "\r\n";
			return "\r";		// Mac
		}
		if(tmp == '\n')			// Unix
			return "\n";
	}
	return "";
}

template<typename T>
size_t container_size(const T& container) {
	return container.size();
}

inline istream& skip_line(istream& is, const string& eol, const unsigned times = 1){
	for(unsigned i=0; i<times; i++)
		is.ignore(std::numeric_limits<std::streamsize>::max(), eol[eol.size()-1]);
	return is;
}

inline double quant_log_pdf(Quant_t x, const double mu=0, const double sigma=1) {
	if(saint::isnan(x)) {
		throw "isnan";
		/*
		  const double tmp = log_pnorm(t,mu,sigma);
		  if (tmp == -std::numeric_limits<double>::infinity()) // prevent underflow
		  return log(std::numeric_limits<double>::min());
		  return tmp;
		*/
		// return log_pnorm(x,mu,sigma);
	} else {
		if(abs(x-mu)/sigma>5)
			x=(signbit(static_cast<double>(x))?-1:1)*5*sigma+mu;
		return normal_log_pdf(x,mu,sigma);
	}
}

struct Model_data {
	// data
	const size_t nprey;
	const size_t nbait;

	const Quantmat test_mat_DATA;
	const Fastmat<valarray<Quant_t> > test_mat;
	const Fastmat<valarray<Quant_t> > test_mat1;
	const vector<unsigned char> n_rep_vec;
	const vector<vector<size_t> > p2p_mapping;
	// const double t;
	const vector<Quant_t> t;
	const size_t n_ctrl_ip;

	// parameters
	Fastmat<valarray1<bool> > Z;

	double beta0, beta1, gamma;

	const vector<double> eta;
	vector<double> d;
	const valarray<double> sd_true;
	const valarray<double> sd_false;


	void print_data_summary(const deque<PreyClass> &PDATA,const vector<string> &ubait, const deque<UIClass> &UIDATA) const;

	double llikelihood() const; // full log likelihood
	double loglikelihood_Z(const size_t i, const size_t j, const size_t rep, const Fastmat<vector<boost::array<double, 2> > >& pre_calc_loglik) const;
	// double llik_theta(const std::vector<double> &x, std::vector<double> & /*grad*/, const size_t i);
	boost::array<Fastmat<double>, 3> calculateScore() const;
	vector<double> get_MRF_parameters() const {
		return {beta0,beta1,gamma};
	}
	void set_MRF_parameters(double beta0_, double beta1_, double gamma_){
		beta0=beta0_, beta1=beta1_, gamma=gamma_;
	}
	void set_MRF_parameters(const vector<double>& a){
		if(a.size()!=3)
			throw runtime_error("3 MRF parameters needed");
		beta0=a[0], beta1=a[1], gamma=a[2];
	}
	void print_MRF_parameters(ostream& out = std::cout) const {
		out << "";
	}
	void wrt_MRF();
	void wrt_MRF_gamma_0();
	void wrt_d();
	void icm_Z();
	void list_output(const string& filename, const vector<string>& prey_list, const deque<UIClass>& UIDATA, const Quantmat& ctrl_mat_DATA, vector<size_t>& ip_idx_to_bait_no, const Fastmat<double>& average_score, const Fastmat<double>& maximum_score, const Fastmat<double>& topo_average_score, const Fastmat<double>& topo_maximum_score, const Fastmat<double>& topo_odds_score) const;
private:
	Fastmat<vector<boost::array<double, 2> > > precalculate_densities() const;
	double llik_MRF_gamma_0( const std::vector<double> &x, std::vector<double>& /*grad*/) const;
	double llik_MRF( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum) const;
	double llik_gamma( const std::vector<double>& x, std::vector<double>& /*grad*/, const Fastmat<double>& gsum, const Fastmat<vector<boost::array<double, 2> > >& log_densities) const;
};

struct Options {
	double f;
	unsigned R;
	int L;
	std::vector<string> input_files;
};

#include "../va/va.hpp"

#endif /* GLOBALS_HPP_ */
