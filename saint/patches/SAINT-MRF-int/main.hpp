#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/array.hpp>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "globals.hpp"


using namespace std;

namespace std1 {
	template<typename _CharT, typename _Traits, typename _Allocator>
	std::basic_istream<_CharT,_Traits>&
	getline(std::basic_istream<_CharT,_Traits>& __is,
			std::basic_string<_CharT,_Traits,_Allocator>& __str, const std::basic_string<_CharT,_Traits,_Allocator>& string_delim) {
		switch(string_delim.size()){
		case 1:
			return std::getline(__is, __str, string_delim[0]);
		case 2:
		{
			std::basic_istream<_CharT,_Traits>& tmp = std::getline(__is, __str, string_delim[1]);
			if(__str.size() > 0){
				
				if(__str[__str.size() - 1] == string_delim[0])
					__str.erase(__str.size() - 1, 1);
			}
			return tmp;
		}
		default:
			throw runtime_error("invalid delimiter");
		}
	}
}

double str2dbl(string ch) {
	double ret;
	std::istringstream iss(ch);
	iss >> ret;
	return ret;
}

// Function returns a vector of strings for the given string object.
// The string is split on white space characters.
std::vector<string> splitString(const string& input_txt) {
	std::vector<string> strVec;
	boost::algorithm::split(strVec, input_txt, boost::is_any_of("\t"));
	return strVec;
}


// Function gets the dimensions of the input file for creation of the matrix
void getFileDimensions(string inputFile, int& nr, int& nc) {

	// open the file and count number lines
	ifstream inF;
	int lineCtr = 0;
	int numCols = 0;
	int numRows = 0;
	string line;
	std::vector<string> v;

	inF.open(inputFile, ios::in);

	if(!inF.is_open()) {
		throw runtime_error("\nERROR! Unable to open " + inputFile + ". File not found\n\n");
		cerr << "\nERROR! Unable to open " << inputFile << ". File not found\n\n";
		exit(0);
	}
	while( !inF.eof() ) {
		line = "";
		getline(inF, line);
		if(line.length() > 10) lineCtr++;
	}
	numRows = lineCtr - 1; // the first row is the header line, we ignore it
	inF.clear(); // clear eof bit
	inF.seekg(0, ios::beg); // set file pointer to beginning of file again

	// count number of columns in file
	line = "";
	getline(inF, line);
	v = splitString(line);
	numCols = v.size() - 1; // first element is a row header

	inF.close(); // close input file stream

	nr = numRows;
	nc = numCols;
}

/**************** Interaction *****************/
deque<InterClass> parseInterFile(const string& inputFile, size_t& ninter) {
	cout << "Parsing interaction file " << inputFile << " ..." << std::flush;
	deque<InterClass> IDATA;
	ifstream inF(inputFile);
	const string line_ending = detect_line_ending(inputFile);
	ninter = 0;
	while( !inF.eof() ) {
		string line = "";
		std1::getline(inF, line, line_ending);

		vector<string> curLineVec = splitString(line);
		if( curLineVec.size() < 4) continue;
		IDATA.push_back( InterClass() );
		InterClass& tmp = IDATA.back();
		// InterClass tmp;
		tmp.set_ipId( curLineVec.at(0) );
		tmp.set_baitId( curLineVec.at(1) );
		tmp.set_preyId( curLineVec.at(2) );
		tmp.set_quant( log(atof(curLineVec.at(3).c_str())) ); // This line is changed to take log of the quant. atof converts a string to a double. log is base e.
		// IDATA.push_back( tmp );
		ninter++;
	}
	inF.close();
	cout << "done." << endl;
	return IDATA;
}

void mapRowCol(deque<InterClass>& IDATA, const deque<PreyClass>& PDATA, const deque<BaitClass>& BDATA, const map<string, BaitClass>& bait_Id_map) {
	for(deque<InterClass>::iterator m = IDATA.begin(); m != IDATA.end(); m++) {
		// lookup prey and get rowId
		// deque<PreyClass>::const_iterator mp =
		// 	find_if(PDATA.begin(), PDATA.end(),
		// 			[&m](const PreyClass& a) {return a.get_preyId() == m->get_preyId();});
		deque<PreyClass>::const_iterator mp = PDATA.begin();
		for(; mp != PDATA.end(); mp++)
			if(mp->get_preyId() == m->get_preyId())
				break;

		if(mp == PDATA.end())
			throw runtime_error("prey " + m->get_preyId() + " not found");

		m->set_rowId( mp->get_rowId() );
		// lookup bait and get colId
		// deque<BaitClass>::const_iterator mb =
		// 	find_if(BDATA.begin(), BDATA.end(),
		// 			[&m](const BaitClass& a) {return a.get_ipId() == m->get_ipId();});
		deque<BaitClass>::const_iterator mb = BDATA.begin();
		for(; mb != BDATA.end(); mb++)
			if(mb->get_ipId() == m->get_ipId())
				break;
		if(mb == BDATA.end())
			throw runtime_error("bait not found");
		m->set_colId( mb->get_colId() );

		m->is_ctrl = bait_Id_map.at(m->get_ipId()).get_isCtrl();
	}
}



/****************     Prey     *****************/

set<string> parsePreyFile(deque<PreyClass>& PDATA, const string& inputFile, size_t& nprey) {
	// This function parses the prey file and returns the number of preys

	cout << "Parsing prey file" << inputFile << " ..." << std::flush; // flush is used to print the output to the screen

	ifstream inF(inputFile); // open the file
	const string line_ending = detect_line_ending(inputFile); // detect line ending

	nprey = 0; // initialize number of preys to 0

	set<string> prey_Id_set; // Initialize a set to store prey IDs
	string line; // initialize a string to store each line of the file
	while(std1::getline(inF, line, line_ending)) { // read each line of the file
		vector<string> curLineVec = splitString(line); // split the line on white space characters
		string prey_id, gene_name; // initialize strings to store prey ID and gene name
		double prey_length = 0.0; // initialize a double to store prey length
		double imputed_mu = 0.0; // initialize a double to store imputed mu
		prey_id = curLineVec.at(0); // first column is prey ID

		
		using boost::conversion::try_lexical_convert; // try_lexical_convert is used to convert a string to a number
		using boost::lexical_cast; // lexical_cast is used to convert a string to a number
		if(curLineVec.size()==1){ // if there is only one column, assume it is a gene name
			printf("Warning: prey file has only one column. Assuming it is a gene name.\n");
			gene_name = curLineVec.at(0); // gene name is the first column
		}else if(curLineVec.size()==2){ // if there are two columns
			const bool is_prey_length = try_lexical_convert(curLineVec.at(1) ,prey_length); // try to convert the second column to a number
			// if second column is not a number, assume it is a gene name
			gene_name = is_prey_length? curLineVec.at(0) : curLineVec.at(1); // if second column is a number, first column is prey ID and second column is gene name
		}else if(curLineVec.size()>2){
			// if there are three columns, first column is prey ID, second column is gene name, and third name is imputed mu
			gene_name = curLineVec.at(1); // second column is gene name
			
			imputed_mu = lexical_cast<double>(curLineVec.at(2)); // third column is imputed mu
			
			// prey_length = lexical_cast<double>(curLineVec.at(1)); // second column is prey length
			// gene_name = curLineVec.at(2); // third column is gene name
		}

		// prey_length = curLineVec.size()>=2 ? atof(curLineVec.at(1).c_str()):0;
		// gene_name = curLineVec.size()<3 ? curLineVec.at(0) : curLineVec.at(2);

		PreyClass tmp; // create a PreyClass object
		tmp.set_rowId( nprey ); // set rowId to be the current prey number
		tmp.set_preyId( prey_id ); prey_Id_set.insert(prey_id); // set preyId to be the current prey ID
		// tmp.set_preyLength( prey_length ); // set preyLength to be the current prey length
		tmp.set_preyGeneId( gene_name ); // set preyGeneId to be the current prey gene name
		tmp.set_imputedMu( imputed_mu ); // set imputedMu to be the current imputed mu

		PDATA.push_back( tmp ) ; // add the PreyClass object to the PDATA deque

		nprey++; // increment the number of preys

	}
	if(nprey != prey_Id_set.size()) // check for duplicate preys
		throw runtime_error("duplicate preys in prey file");
	inF.close();

	// return nprey;
	cout << "done." << endl;
	return prey_Id_set; // return the set of prey IDs
}

/****************     Bait     *****************/
map<string, BaitClass> parseBaitFile(deque<BaitClass>& BDATA, const string& inputFile, size_t& nip, size_t& n_test_ip) {
	cout << "Parsing prey file " <<  inputFile << " ..." << std::flush;

	ifstream inF;
	const string line_ending = detect_line_ending(inputFile);

	nip = 0;

	inF.open(inputFile.c_str(), ios::in);
	map<string, BaitClass> bait_Id_map;
	n_test_ip = 0;
	size_t n_ctrl_ip = 0;
	while( !inF.eof() ) {
		string line;
		std1::getline(inF, line, line_ending);

		if(line == "") continue;
		vector<string> curLineVec = splitString(line);

		if(curLineVec.size() != 3) {
			cout << nip << endl;
			cout << curLineVec.size() << endl;
			cout << curLineVec[0] << endl;
			throw runtime_error("Bait file must have 3 columns");
		}

		string node1, node2, node3; // two nodes to an interaction
		node1 = curLineVec.at(0);
		node2 = curLineVec.at(1);
		node3 = curLineVec.at(2);

		BaitClass tmp;
		if(node3 != "T" && node3 != "C")
			throw runtime_error("3rd column of bait file must be 'T' or 'C'");
		tmp.set_colId( node3 == "T" ? n_test_ip++ : n_ctrl_ip++ );
		tmp.set_ipId( node1 );
		tmp.set_baitId( node2 );
		tmp.set_isCtrl( node3 != "T" );


		BDATA.push_back( tmp );
		bait_Id_map[node1] = tmp;

		nip++;
	}

	inF.close();
	cout << "done." << endl;
	return bait_Id_map;
}



void sortBaitData( deque<BaitClass>& BDATA ) {

	deque<BaitClass> newBDATA;

	deque<BaitClass>::iterator m;

	for(m = BDATA.begin(); m != BDATA.end(); m++)
		if(m->get_isCtrl() == true)
			newBDATA.push_back( *m );

	for(m = BDATA.begin(); m != BDATA.end(); m++)
		if(m->get_isCtrl() == false)
			newBDATA.push_back( *m );

	if(BDATA.size() != newBDATA.size())
		cerr << "Bait sorting error" << endl;

	BDATA = deque<BaitClass>(newBDATA);
	// BDATA = deque<BaitClass>();
	// m = newBDATA.begin();
	// while( m != newBDATA.end() ) {
	//	BDATA.push_back( *m );
	//	m++;
	// }
}

/********************* unique interaction *********************/
void createList( deque<UIClass>& UIDATA, const deque<InterClass>& IDATA, const deque<BaitClass>& BDATA, const deque<PreyClass>& PDATA,const size_t nprey, const size_t nbait, const vector<string>& ubait, vector<size_t>& ip_idx_to_bait_no, size_t& nuinter ) {
	//ip_idx_to_bait_no = {};
	ip_idx_to_bait_no.clear();
	{
		map<string, size_t> ubait_map;
		for(size_t i = 0; i < ubait.size(); ++i)
			ubait_map[ubait[i]] = i;
		BOOST_FOREACH(const auto& bait , BDATA){
			if( !bait.get_isCtrl() )
				ip_idx_to_bait_no.push_back( ubait_map.at(bait.get_baitId()) );
		}
	}
	UIDATA.clear();
	// UIDATA must be a deque. deque::push_back() does not invalidate references.
	Fastmat<UIClass*> UImat(nprey, nbait);
	BOOST_FOREACH(const auto& inter , IDATA)
		if(!inter.is_ctrl){
			const size_t j = inter.get_colId();
			const size_t b = ip_idx_to_bait_no[j];
			const size_t i = inter.get_rowId();
			if(!UImat(i,b)){
				UIDATA.push_back({});
				UIClass& tmp = UIDATA.back();
				tmp.set_baitId( inter.get_baitId() );
				tmp.set_preyId( inter.get_preyId() );
				tmp.set_preyGeneId( PDATA[i].get_preyGeneId() );
				tmp.set_rowId( i );
				// tmp.add_colId( j );
				UImat(i,b) = &tmp;
			}
			UImat(i,b)->add_colId( j );
		}
	nuinter = UIDATA.size();
	return;


	// original code below does the same thing as above but slower
	UIDATA.clear();
	nuinter = 0;
	map<string, const BaitClass*> bait_map;
	BOOST_FOREACH(const auto& bait , BDATA)
		bait_map[bait.get_ipId()] = &bait;

	for(auto m = IDATA.begin(); m != IDATA.end(); m++) {
		auto bm = bait_map.find(m->get_ipId());
		if(bm == bait_map.end())
			throw runtime_error("ipId not found");
		bool tmpIsCtrl = bm->second->get_isCtrl();
		if(tmpIsCtrl == false) {
			deque<UIClass>::iterator um;
			for(um = UIDATA.begin(); um != UIDATA.end(); um++)
				if( (um->get_preyId() == m->get_preyId())
					&&
					(um->get_baitId() == m->get_baitId()) )
					break;
			// if identical, the just add colId
			// if not, then create new entry for ui
			if(um == UIDATA.end()) { // nothing with the same bait and prey IDs found
				// UIClass& tmp = UIDATA.back();
				UIClass tmp;
				tmp.set_baitId( m->get_baitId() );
				tmp.set_preyId( m->get_preyId() );
				tmp.set_rowId( m->get_rowId() );
				tmp.add_colId( m->get_colId() );
				UIDATA.push_back( tmp );
				nuinter++;
			}
			else {
				um->add_colId( m->get_colId() );
			}
		}
	}
	// return nuinter;
}

int get_nexpr( deque<BaitClass>& BDATA ) {

	deque<BaitClass>::iterator m;
	int ret = 0;

	for(m = BDATA.begin(); m != BDATA.end(); m++) {
		if( m->get_isCtrl() == false ) ret++;
	}

	return ret;

}

int get_nctrl(const deque<BaitClass>& BDATA ) {
	int ret = 0;
	for(auto m = BDATA.begin(); m != BDATA.end(); m++)
		if( m->get_isCtrl() == true ) ret++;
	return ret;
}

void createMatrixData( Quantmat& test_mat_DATA, Quantmat& ctrl_mat_DATA, const deque<InterClass>& IDATA) {
	for(deque<InterClass>::const_iterator m = IDATA.begin(); m != IDATA.end(); m++) {
		int r = m->get_rowId(), c = m->get_colId();
		Quant_t q = m->get_quant();
		m->is_ctrl ? ctrl_mat_DATA(r, c) = q : test_mat_DATA(r, c) = q;
	}
}
void zero_self_interactions(Quantmat& test_mat_DATA, const deque<UIClass>& UIDATA) {
	BOOST_FOREACH(const auto& UI , UIDATA) {
		if(UI.get_preyGeneId() == UI.get_baitId()) {
			BOOST_FOREACH(const int col , UI.get_colId())
				test_mat_DATA(UI.get_rowId(), col) = saint::nan;
		}
	}
}

void Model_data::print_data_summary(const deque<PreyClass>& PDATA,const vector<string>& ubait, const deque<UIClass>& UIDATA) const {
	set<string> uprey;
	for(size_t i = 0; i<nprey; ++i)
		uprey.insert( PDATA.at(i).get_preyGeneId() );
	if(0)
		if(PDATA.size() != uprey.size())
			throw std::runtime_error("prey list has duplicates!");
	cout << "number of:\npreys: " << uprey.size() << '\n';
	cout << "unique baits: " << ubait.size() << '\n';
	cout << "ctrl IPs: " << n_ctrl_ip << '\n';
	cout << "non zero unique interactions: "<< (UIDATA.size()) << '\n';
	vector<string> tmp;
	BOOST_FOREACH(const string& ub , ubait)
		if(uprey.count( ub ) == 0)
			tmp.push_back(ub);
	if(tmp.size() > 0){
		cout << "The following baits are not preys:" << '\n';
		BOOST_FOREACH(const string& t , tmp) cout<<t<<'\n';
		cout << "-----" << '\n';
	}
}

// class inherited from ofstream to convert nan to “.”
struct outputstream : public std::ofstream {
	using std::ofstream::operator<<;
	inline outputstream& operator<<(const Quant_t value){
		if(saint::isnan(value))
			return static_cast<outputstream&>(static_cast<std::ofstream&>(*this) << ".");
		return static_cast<outputstream&>(static_cast<std::ofstream&>(*this) << value);
	}
	inline outputstream& operator<<(double value){
		if(std::isnan(value))
			return static_cast<outputstream&>(static_cast<std::ofstream&>(*this) << ".");
		return static_cast<outputstream&>(static_cast<std::ofstream&>(*this) << value);
	}
};
template<typename T>
inline outputstream& operator<<(outputstream& out, T value){
 	return static_cast<outputstream&>(static_cast<std::ofstream&>(out) << value);
}

void Model_data::list_output(const string& filename, const vector<string>& prey_list, const deque<UIClass>& UIDATA, const Quantmat& ctrl_mat_DATA, vector<size_t>& ip_idx_to_bait_no, const Fastmat<double>& average_score, const Fastmat<double>& maximum_score, const Fastmat<double>& topo_average_score, const Fastmat<double>& topo_maximum_score, const Fastmat<double>& topo_odds_score) const {
	boost::array<string, 18> column_names = {{"Bait", "Prey", "PreyGene", "Intensity", "IntensitySum", "AvgIntensity", "NumReplicates", "ctrlIntensity", "AvgP", "MaxP", "TopoAvgP", "TopoMaxP","SaintScore","OddsScore","FoldChange", "BFDR", "ImputedValue", "boosted_by"}};
	const string delim = "\t";
	outputstream f;
	// std::ofstream f;
	f.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	f.open(filename);
	f << std::setprecision(3) << std::fixed;
	// write column names
	for(size_t i = 0; i < column_names.size() - 1; ++i)
		f << column_names[i] << delim;
	f << column_names.back() << '\n';
	// write data
	BOOST_FOREACH(const auto& ui , UIDATA) { // for each unique interaction
		size_t row = ui.get_rowId(),
			baitcol = ip_idx_to_bait_no[ui.get_colId().front()];
		f << ui.get_baitId() << delim // Bait column
		  << ui.get_preyId() << delim // Prey column
		  << ui.get_preyGeneId() << delim; // PreyGene column
		const valarray<Quant_t> quant_orig = exp(test_mat(row, baitcol)); // There is no definition for test_mat anywhere in this file. I assume it is a member of Model_data.
		valarray<Quant_t> quant = quant_orig; // This is the quant data that will be modified
		// const valarray<double>& quant = reinterpret_cast<const valarray<double>&>(quant_orig);
		for(size_t i=0;i<quant.size();++i) // replace nan with exp(t)
			if(saint::isnan(quant[i])) // if quant[i] is nan
				quant[i] = exp(t[i]); // replace it with exp(t)
		double AvgQuant = quant.sum()/quant.size(); // average of quant data
		f << quant_orig[0]; // Intensity column
		for(size_t col = 1; col < quant_orig.size(); ++col) // for each column
			f << '|' << quant_orig[col]; // write the intensity
		f << delim;
		f << quant.sum() << delim // IntensitySum column
		  << AvgQuant << delim // AvgIntensity column
		  << quant.size() << delim; // NumReplicates column
		const valarray<Quant_t> ctrlCounts = exp(ctrl_mat_DATA(row,__));
		f << ctrlCounts[0]; // ctrlIntensity column
		for(size_t j = 1;j < ctrlCounts.size(); ++j)
			f << "|" <<ctrlCounts[j];
		f << delim;
		double avg_p = average_score(row, baitcol); // This is the average score for the prey
		unsigned denom = 0;
		double numer = 0;
		for(unsigned i=0; i < average_score.size(); ++i) { // for each prey
			bool tmp = average_score[i] > avg_p; // if prey i has a higher average score than the current prey
			numer += average_score[i]*tmp; // add the average score of prey i to the numerator if prey i has a higher average score than the current prey
			denom += tmp; // increment the denominator if prey i has a higher average score than the current prey
		}
		double BFDR = (denom==0 ? 0 : 1-numer/denom); // BFDR is the number of preys with a higher average score than the current prey divided by the total number of preys
		double max_p = maximum_score(row, baitcol); // This is the maximum score for the prey
		double topo_avg_p  = topo_average_score(row, baitcol); // This is the topological average score for the prey
		// topo_avg_p = max(topo_avg_p, avg_p);
		double topo_max_p = topo_maximum_score(row, baitcol); // This is the topological maximum score for the prey
		// topo_max_p = max(topo_max_p, max_p);
		f
			//<< exp(eta[row]) << delim << exp(eta[row] + d[row]) << delim
			// << sd_true[row] << delim
			// << sd_false[row] << delim
			<< avg_p << delim // AvgP column
			<< max_p << delim // MaxP column
			<< topo_avg_p << delim // TopoAvgP column
			<< topo_max_p << delim // TopoMaxP column
			<< max(topo_avg_p, avg_p) << delim // SaintScore column
			<< topo_odds_score(row, baitcol) << delim
			// << AvgQuant/eta[row] << delim
			<< AvgQuant/exp(eta[row]) << delim //FoldChange
			<< BFDR << delim
			<< exp(t[row]) << delim; // ImputedValue column
		BOOST_FOREACH(const unsigned l , p2p_mapping[row]){
			if(Z(l, baitcol).mean() > 0)
				f << prey_list[l] << '|'; // boosted_by column
		}
		f << '\n';
	}
}