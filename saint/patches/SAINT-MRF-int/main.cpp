/*
  Copyright (C) <2011>  <Hyungwon Choi, Damian Fermin>
  For troubleshooting, contact hyung_won_choi@nuhs.edu.sg.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You can obtain a copy of the GNU General Public License from
  <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <iomanip>

#include "BaitClass.hpp"
#include "PreyClass.hpp"
#include "InterClass.hpp"
#include "UIClass.hpp"
#include "main.hpp"
#include "Mapping.hpp"
#include "Stats.hpp"

using namespace std;


// ICM
Model_data icms(Model_data& dp, const bool with_gamma) {
	double oldllik, newllik = dp.llikelihood(); // This line defines two doubles and sets them to the log likelihood of the statistical model.
	if(!with_gamma) dp.gamma = 0;else dp.gamma=0; // This line sets the gamma parameter of the statistical model to zero if the gamma parameter is not used. Otherwise, it sets the gamma parameter of the statistical model to the value of the gamma parameter.
	for (unsigned iter = 1;iter <= 15; ++iter) { // This loop iterates through the number of iterations.
		dp.print_MRF_parameters(); // This line prints the parameters of the statistical model.
		oldllik = newllik; // This line sets the old log likelihood to the new log likelihood.
		// But aren't the old and new log likelihoods the same?  - Yes, they are the same.
		dp.icm_Z(); // This line performs the ICM algorithm on the statistical model.
		// dp.wrt_d();//TODO
		if(with_gamma)
			dp.wrt_MRF();
		else
			dp.wrt_MRF_gamma_0();
		newllik = dp.llikelihood();
		if( (newllik >= oldllik) &&
			(exp(newllik - oldllik) -1 < 1e-3) ) break;
	}
	dp.print_MRF_parameters();
	return dp;
}

boost::array<Fastmat<double>, 3> computation(const bool with_gamma, Model_data& dp, ostream& parameters_output_file) {
	// This function computes the scores for the statistical model. It also returns the scores.
	const auto dp0 = icms(dp, with_gamma); // This line defines a Model_data object and sets it to the statistical model. The icms function is called to compute the scores for the statistical model. The scores are stored in the Model_data object.
	const auto all_scores = dp0.calculateScore();
	// cout << "MRF parameters estimates: ";
	dp0.print_MRF_parameters();
	// write estimated MRF parameters to a file
	dp0.print_MRF_parameters(parameters_output_file);
	return all_scores;
}

#include <boost/program_options.hpp>
Options program_arguments(int argc, char* argv[]) {
	namespace po = boost::program_options;
	double f;
	unsigned R;
	int L;
	string config_file;

	// Declare a group of options that will be
	// allowed only on command line
	po::options_description generic("Generic options");
	generic.add_options()
		("version,v", "print version string")
		("help,h", "produce help message")
		("config,c", po::value<string>(&config_file),//->default_value("saintMRF.cfg"),
		 "name of a file of a configuration.")
		;

	// Declare a group of options that will be
	// allowed both on command line and in
	// config file
	po::options_description config("Configuration");
	config.add_options()
		("frequency,f", po::value<double>(&f)->default_value(0.5),
		 "frequency")
		("ncounts,R", po::value<unsigned>(&R)->default_value(100),
		 "number of maximum counts to use for probability calculation")
		("ncontrols,L", po::value<int>(&L)->default_value(100),
		 "number of controls for each prey that is used in the calculations")
		;

	// Hidden options, will be allowed both on command line and
	// in config file, but will not be shown to the user.
	po::options_description hidden("Hidden options");
	hidden.add_options()
		("input-file", po::value< vector<string> >(), "input file")
		;


	po::options_description cmdline_options;
	cmdline_options.add(generic).add(config).add(hidden);

	po::options_description config_file_options;
	config_file_options.add(config).add(hidden);

	po::options_description visible("Allowed options");
	visible.add(generic).add(config);

	po::positional_options_description p;
	p.add("input-file", -1);

	po::variables_map vm;
	store(po::command_line_parser(argc, argv).
		  options(cmdline_options).positional(p).run(), vm);
	notify(vm);
	if(vm.count("help") || vm.count("version")) {
		if (vm.count("version")) {
			cout << "SAINT express version 3.6.3\n";
		}
		if (vm.count("help")) {
			cout << "\nUSAGE: SAINTexpress-int [OPTIONS] [<interaction data> <prey data> <bait data>] [known interaction data]\n";
			cout << "OUTPUT: a text file 'list.txt'. See manual for details on the output.\n\n"
				"Example uasge 1: SAINTexpress-int inter.dat prey.dat bait.dat\n"
				 << "Example uasge 2: SAINTexpress-int -L4 inter.dat prey.dat bait.dat\n"
				"\tOnly the highest 4 control counts will be used.\n"
				 << "Example uasge 3: SAINTexpress-int -L4 inter.dat prey.dat bait.dat GO.txt\n"
				"\tA file containing known interactions 'GO.txt' is utilized in the computation of the scores.\n"
				 << "Example uasge 4: SAINTexpress-int -R2 inter.dat prey.dat bait.dat\n"
				"\tOnly the 2 replicates in each interaction with the highest counts is involved in the computation of the scores.\n"
				 << "Get version: SAINTexpress-int -v\n";
			cout << visible << "\n";
		}
		exit(0);
	}
	if(config_file=="") {
		config_file = "saintMRF.cfg";
		// cout << "Config file argument not supplied. Using " << config_file << " as config file\n";
	}
	ifstream ifs(config_file.c_str());
	if (!ifs)
		// cout << "Cannot open config file: " << config_file << endl
		;
	else
	{
		cout << "Using config file: " << config_file << endl;
		store(parse_config_file(ifs, config_file_options), vm);
		notify(vm);
	}
	vector<string> input_files_vec;
	// boost::array<string, 4> input_files;
	if (vm.count("input-file"))
	{
		input_files_vec=vm["input-file"].as< vector<string> >();
		if(input_files_vec.size()==4)
			;
		// input_files_vec = vm["input-file"].as< boost::array<string, 4> >();
		else if(input_files_vec.size() == 3)
			;
		else {
			cerr << "four input files needed.\n";
			exit(1);
		}
	}
	else {
		input_files_vec = {"inter.dat", "prey.dat", "bait.dat"};
		cout << "Input files not supplied, using defaults.\n";
	}
	cout << "Input files are: "
		 << input_files_vec[0] << ", "
		 << input_files_vec[1] << ", "
		 << input_files_vec[2];
	if(input_files_vec.size()==4)
		cout << ", " << input_files_vec[3];
	cout << endl;
	if(f<=0 || 1<=f) {
		cerr << "Frequency (f) must be between 0 to 1." << endl;
		exit(0);
	}
	if(R < 1) {
		cerr << "number of maximum counts to use for probability calculation (n) must be a positive integer." << endl;
		exit(0);
	}
	if(L < 1) {
		cerr << "number of controls (L) must be a positive integer." << endl;
		exit(0);
	}
	// cout << "f is " << f << "\n";
	// cout << "n is " << R << "\n";
	// cout << "L is " << L << "\n";
	return {f, R, L, input_files_vec};
}


Options opts;
int main(int argc, char* argv[]) {
	

	// define input files
	/*const Options*/ opts = program_arguments(argc, argv);
	const string inputInter = opts.input_files[0]; // interaction file
	const string inputPrey = opts.input_files[1]; // prey file
	const string inputBait = opts.input_files[2];
	const string inputGO = opts.input_files.size()==4 ? // GO file - optional
		opts.input_files[3]:"";
	
	// define variables
	size_t ninter; // number of interactions
	size_t nuinter; // number of unique interactions
	size_t nbait; // number of baits
	size_t nprey; // number of preys
	size_t nip; // number of IPs

	/********************/
	/* Reading the data */
	/********************/

	cerr << "Interaction file: \"" << inputInter << "\"" <<endl;
	cerr << "Prey file: \"" << inputPrey << "\"" << endl;
	cerr << "Bait file: \"" << inputBait << "\"" << endl;
	cerr << "GO file: \"" << inputGO << "\"" << endl;

	/* Reading in the data */
	// get prey data
	deque<PreyClass> PDATA; // This line defines a deque of PreyClass objects. A Deque is a double-ended queue. A queue is a FIFO (first in first out) data structure. A deque is a queue that can be accessed from both ends.
	parsePreyFile(PDATA, inputPrey, nprey); // This function parses the prey file and stores the data in the deque PDATA. It also returns the number of preys.
	


	// get bait data
	deque<BaitClass> BDATA; // This line defines a deque of BaitClass objects.
	size_t n_test_ip; // number of test IPs
	const map<string, BaitClass> bait_Id_map = parseBaitFile(BDATA, inputBait, nip, n_test_ip); // This function parses the bait file and stores the data in the deque BDATA. It also returns a map that maps bait IDs to BaitClass objects.
	// get interaction data
	deque<InterClass> IDATA = parseInterFile(inputInter, ninter); // This function parses the interaction file and stores the data in the deque IDATA. It also returns the number of interactions.
	cout << "Setting matrix indices for each interaction..." << flush;
	mapRowCol( IDATA, PDATA, BDATA, bait_Id_map);  // mapping interaction data to matrix data. This function sets the prey and bait indices for each interaction.
	cout << "done." << endl;

	/***************************/
	/* Creating data structure */
	/***************************/

	// Matrix data format
	// MatrixClass QDATA(nprey, nip);
	const size_t n_ctrl_ip = get_nctrl( BDATA ); // number of control IPs
	// Quantmat test_mat_DATA(nprey, nip - n_ctrl_ip);
	// Quantmat ctrl_mat_DATA(nprey, n_ctrl_ip);
	Quantmat test_mat_DATA(nprey, nip - n_ctrl_ip, saint::nan); // This line defines a Quantmat object. A Quantmat object is a matrix of Quant_t objects. A Quant_t object is a double. The Quantmat object is a matrix of doubles.
	Quantmat ctrl_mat_DATA(nprey, n_ctrl_ip, saint::nan); // This line defines a Quantmat object. A Quantmat object is a matrix of Quant_t objects. A Quant_t object is a double. The Quantmat object is a matrix of doubles.
	// createMatrixData( QDATA, IDATA );
	cout << "Creating matrix..." << flush;
	createMatrixData( test_mat_DATA, ctrl_mat_DATA, IDATA); // This function creates the matrix data from the interaction data. It also returns the number of unique interactions.
	cout << "done." << endl;
	vector<double> tmp; // This line defines a vector of doubles. It is used to store the values of the matrix data for each interaction.
	// BOOST_FOREACH(const Quant_t v, test_mat_DATA.mat){
	for(auto i=begin(test_mat_DATA.mat); i!=end(test_mat_DATA.mat);++i){ // This loop iterates through the matrix data.
		const Quant_t v = *i; // This line defines a Quant_t object and sets it to the value of the matrix data for the current interaction.
		if(!saint::isnan(v) && -1e5<v) // filter out log(0)
			tmp.push_back(v); // This line adds the value of the matrix data for the current interaction to the vector tmp.
	}
	// BOOST_FOREACH(const Quant_t v, ctrl_mat_DATA.mat){
	for(auto i=begin(ctrl_mat_DATA.mat); i!=end(ctrl_mat_DATA.mat); ++i){ // This loop iterates through the matrix data.
		const Quant_t v = *i; // This line defines a Quant_t object and sets it to the value of the matrix data for the current interaction.
		if(!saint::isnan(v) && -1e5<v) // filter out log(0)
			tmp.push_back(v); // This line adds the value of the matrix data for the current interaction to the vector tmp.
	}
	test_mat_DATA.mat -= mean(tmp); // This line subtracts the mean of the values of the matrix data from the matrix data.
	test_mat_DATA.mat /= sqrt(var1(tmp)); // This line divides the values of the matrix data by the square root of the variance of the values of the matrix data.
	ctrl_mat_DATA.mat -= mean(tmp); // This line subtracts the mean of the values of the matrix data from the matrix data.
	ctrl_mat_DATA.mat /= sqrt(var1(tmp)); // This line divides the values of the matrix data by the square root of the variance of the values of the matrix data.	
	// Unique interaction

	// Add a print statement to see what the imputed mu value is for the first prey
	cout << "The imputed mu of the first prey is " << PDATA.at(0).get_imputedMu() << endl;

	// I need to iterate through the prey data and normalize the imputed mu values in the same way as the matrix data.
	for(auto i=begin(PDATA); i!=end(PDATA); ++i){ // This loop iterates through the prey data.
		PreyClass& p = *i; // This line defines a reference to a PreyClass object and sets it to the current prey.
		if(!saint::isnan(p.get_imputedMu()) && -1e5<p.get_imputedMu()) // This line checks if the imputed mu value for the current prey is not NaN and is greater than -1e5.
			p.set_imputedMu( (log(p.get_imputedMu()) - mean(tmp)) / sqrt(var1(tmp)) ); // This line normalizes the imputed mu value for the current prey.
			// if statement for if i is 0
			if (i == begin(PDATA)) {
				cout << p.get_imputedMu() << endl;
				cout << mean(tmp) << endl;
				cout << sqrt(var1(tmp)) << endl;
				cout << (p.get_imputedMu() - mean(tmp)) / sqrt(var1(tmp)) << endl;
			}
	}

	// Add print statement for after normalization
	cout << "The imputed mu of the first prey is " << PDATA.at(0).get_imputedMu() << endl;

	// unique baits (test baits only);
	const vector<string> ubait = uniqueBait( BDATA, nbait ); // This line defines a vector of strings and sets it to the unique baits.
	// ip index to bait index
	vector<size_t> ip_idx_to_bait_no; // This line defines a vector of size_t objects. It is used to map IP indices to bait indices.
	cout << "Creating a list of unique interactions..." << flush;
	deque<UIClass> UIDATA; // This line defines a deque of UIClass objects. A UIClass object is a unique interaction object.
	createList( UIDATA, IDATA, BDATA, PDATA, nprey, nbait, ubait, ip_idx_to_bait_no, nuinter ); // This function creates a list of unique interactions. It also returns the number of unique interactions.
	cout << "done." << endl;
	// zero out reflexive interactions
	zero_self_interactions(test_mat_DATA, UIDATA); // This function sets the values of the matrix data for reflexive interactions to zero.

	/*************************/
	/*   Mapping proteins    */
	/* Statistical Analysis  */
	/*************************/
	vector<string> prey_list; // This line defines a vector of strings. It is used to store the prey IDs.
	for(size_t i = 0; i<nprey; ++i) // This loop iterates through the preys.
		prey_list.push_back( PDATA.at(i).get_preyId() ); // This line adds the prey ID for the current prey to the vector prey_list.

	// p2p_mapping is a vector of vectors of size_t objects. It is used to map prey indices to prey indices.	
	const vector<vector<size_t> > p2p_mapping = (inputGO==""? vector<vector<size_t> >(prey_list.size(),vector<size_t>(0)) : createP2Pmap(prey_list, inputGO)); // This line defines a vector of vectors of size_t objects and sets it to the protein to protein mapping. If the GO file is not empty, then the protein to protein mapping is created from the GO file. Otherwise, the protein to protein mapping is set to an empty vector.
	
	//Model_data dp = statModel(p2p_mapping, ubait, test_mat_DATA, ctrl_mat_DATA, ip_idx_to_bait_no, nprey, nbait,PDATA); // This function creates a statistical model. It also returns a Model_data object.
	Model_data dp = statModel(p2p_mapping, ubait, test_mat_DATA, ctrl_mat_DATA, ip_idx_to_bait_no, nprey, nbait,PDATA);
	
	// dp.print_data_summary(PDATA, ubait, UIDATA); // print data summary
	// ofstream parameters_output_file("parameters.txt");
	ostream& parameters_output_file(std::cout); // This line defines a reference to an ostream object and sets it to the standard output stream.
	Model_data dp_gamma0 = dp; // This line defines a Model_data object and sets it to the statistical model.

	// This is the line I need to look at in more detail.
	const boost::array<Fastmat<double>, 3> scores = computation(false/*gamma set to zero*/, dp_gamma0, parameters_output_file);	// This line computes the scores for the statistical model. It also returns the scores.
	Model_data dp_MRF = dp; // This line defines a Model_data object and sets it to the statistical model.
	const auto ctor=Fastmat<double>(nprey,nbait); // This line defines a Fastmat<double> object and sets it to a matrix of doubles with dimensions nprey by nbait.
	boost::array<Fastmat<double>, 3> topo_scores{{ctor,ctor,ctor}}; // This line defines a boost::array of Fastmat<double> objects and sets it to a boost::array of matrices of doubles with dimensions nprey by nbait.
	if(inputGO != "")
		topo_scores = computation(true/*gamma > 0*/, dp_MRF, parameters_output_file); 		// Topo, gamma > 0
	else
		topo_scores = scores;
	// parameters_output_file.close();
	// list_output will print the topo_scores as maximum of scores and topo_scores
	dp_MRF.list_output("list.txt", prey_list, UIDATA, ctrl_mat_DATA, ip_idx_to_bait_no, scores[0], scores[1], topo_scores[0], topo_scores[1], topo_scores[2]); // This function prints the list of unique interactions to a file.
	// matrix_output("matrix.txt", PDATA, /*UIDATA,*/ test_mat_DATA, ctrl_mat_DATA, ip_idx_to_bait_no, topo_scores[0], scores[0], ubait, dp);
	
	return 0;
}
