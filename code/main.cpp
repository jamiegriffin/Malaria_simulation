
#include <sys/stat.h>
#include <omp.h>

#include "library.h"

#include "parameter_map.cpp"
#include "classes.h"
#include "mosq.h"
#include "people.h"
#include "simulation.h"

#include "equilibrium/malaria_parms_demog.cpp"
#include "equilibrium/malaria_parms.cpp"
#include "equilibrium/malaria_state.cpp"

#include "events/event_manager.cpp"
#include "events/event.cpp"

#include "mosq.cpp"
#include "people.cpp"
#include "simulation.cpp"

int main(int argc, char *argv[]){

	const int min_arg=5;
	if(argc-1 < min_arg)
		error_crit("There should be at least " + to_string(min_arg) + " arguments on the command line, there were " + to_string(argc-1));

	const string arg5 = argv[5];
	const string ext5=arg5.substr(arg5.size()-4, 4);
	const bool demog=(ext5==".txt");
	
	if(demog && argc-1 < min_arg+2)
		error_crit("There should be at least " + to_string(min_arg+2) + " arguments on the command line, there were " + to_string(argc-1));

	int arg=1;
	const string root=string(argv[arg++]) + "/";
	const string default_parms=root + string(argv[arg++]);
	const string site_file=root + string(argv[arg++]);
	const string demog_file=demog ? (root + string(argv[arg++])) : string();
	const string pop_file=demog ? (root + string(argv[arg++])) : string();
	const string output_root=root + string(argv[arg++]) + "/";
	const string name=argv[arg++];

	const string output_file=output_root + name + ".txt";
	ofstream out0;
	const string out_cout=output_root  + name + ".cout";

	out0.open(out_cout.c_str(), ios::out);
	streambuf *old_cout = cout.rdbuf(out0.rdbuf());
	streambuf *old_cerr = cerr.rdbuf(out0.rdbuf());

	ifstream input(default_parms.c_str());
	if(!input)
		error_crit("Can't open file "+default_parms+" for input");
	map<string, string> file_map;
	while(!input.eof()){
		string s;
		string t;
		input >> s >> t;
		if(!input.eof()){
			if(!input.good())
				error_crit("Error in importing data from "+default_parms);
			file_map.insert(pair<string, string>(s, t));
		}
	}
	input.close();

	int arg1=arg;
	while(arg1<argc){
		const string s=argv[arg1++];
		if(s=="add")
			continue;
		if(arg1==argc)
			error_crit(
				"Not counting the word add, there should be an even number of additional command line arguments, name1 value1 name2 value2 ...,  after the first " 
					+ to_string(min_arg));
		const string f=argv[arg1++];
		const int i=s.find("file");
		if(i!=string::npos){
			if(file_map.find(s)==file_map.end())
				error_crit(s + " is not one of the input list of files");
			file_map[s]=f;
		}
	}
	
	for (map<string, string>::iterator it = file_map.begin(); it != file_map.end(); it++) {
		if(it->first != "output_file" && it->first != "saved_file")
			import_to_map(root + it->second);
	}
	import_to_map(site_file, true);

	cout << '\n' << datestr() << '\t' << timestr() << '\n';
	cout << "executable\t" <<  argv[0] << '\n';
	cout << "output file\t" << output_file << '\n';
	cout << "root\t" << root << '\n';
	cout << "default_parms\t" << default_parms << '\n';	
	cout << "model_parms\t" << root + file_map["model_file"] << '\n';
	cout << "site_file\t" << site_file << '\n';
		
	cout << "\nParameters in input files\n";
	for(Parameter_map::const_iterator p=parameter_map.begin(); p != parameter_map.end(); p++)
		cout << p->first << '\t' << p->second.first << '\n';
	if(arg<argc){
		cout << '\n';
		cout << "\nParameters modified on the command line\n";
	}


	bool to_add=false;
	while(arg<argc){
		const string s=argv[arg++];
		if(s=="add"){
			to_add=true;
			continue;
		}
		if(arg==argc)
			error_crit(
				"Not counting the word add, there should be an even number of additional command line arguments, name1 value1 name2 value2 ...,  after the first " 
					+ to_string(min_arg));
		const string f=argv[arg++];
		const int i=s.find("file");
		if(i==string::npos){
			const double x=atof(f.c_str());
			to_add ? add_to_map(s, x) : modify_map(s, x);			
			cout << s << '\t' << x << endl;
		}
	}

	//*******************************************************************************
	const bool overwrite=from_map_bool("overwrite");
	if(!overwrite){
		struct stat file_info;
		if(stat(output_file.c_str(),&file_info)==0)
			error_crit("output file " + output_file + " already exists and overwrite = 0, so exiting program");
	}

	const bool parallel_omp=from_map_bool("parallel", true);
	omp_set_num_threads(parallel_omp ? omp_get_max_threads() : 1);

	//*******************************************************************************
	//*******************************************************************************
	const double begin_omp = omp_get_wtime();

	add_to_map("prop_gamb_ss", 1-(from_map("prop_arab", 0, 1, 0) + from_map("prop_fun", 0, 1, 0)));
	modify_itn_cov();

	const int itn_years=find_map_years("itn_cov");
	bool mult_baseline = itn_years>0;
	if(!mult_baseline && from_map_bool("det_baseline", false))
		modify_map("num_runs", 1);
	
	const int recalculate = from_map_int("recalculate", 0, 6, 0);
	if (recalculate>0) {
		const bool baseline_int = from_map_bool("baseline_int", true);
		if (recalculate>=4 && recalculate<=6)
			modify_map("baseline_int", 0);
		const double total_M_new = recalculate_M(mult_baseline, root, demog_file, pop_file);
		modify_map("baseline_int", baseline_int);
		modify_map("total_M", total_M_new);

		if (recalculate==1 || recalculate==4) {
			const double time_taken = omp_get_wtime()-begin_omp;
			ofstream out(output_file.c_str(), ios::out);
			if (recalculate==1 || recalculate==4)
				out << from_map("total_M") << '\t' << total_M_new << '\t' << time_taken << '\n';
		}
		if (recalculate==3 || recalculate==6)
			change_site_file(site_file);
	}
	
	if(recalculate==0 || recalculate==2 || recalculate==5){
		const int num_runs=from_map_int("num_runs");
		const int num_per_yr=from_map_int("output_per_yr");
		const int output_type=from_map_int("output_type", 0, 1);
		const int init_output=from_map_int("init_output", 0, 20);
		const int first_output=output_type==0 ? -init_output : min(from_map_int("output_year0", -20), -init_output); 

		int final_output_year=-1;
		while(in_map("output_year" + to_string(++final_output_year)))
			;
		final_output_year--;
		const int final_run=output_type==0 ? from_map_int("final_run", 0) : 
			ceil(from_map("output_year" + to_string(final_output_year), 0));

		modify_map("init_output", -first_output);
		modify_map("final_run", final_run);

		const string output_vars_file = root + file_map["output_file"];
		vector<Output_var> output_vars;
		vector<string> quantity;
		vector<string> var_names;
		vector<double> age0, age1;
		import_file(output_vars_file, quantity, age0, age1, var_names);
		for(int i=0; i<quantity.size(); i++)
			output_vars.push_back(Output_var(quantity[i], age0[i], age1[i]));

		const int num_vars=quantity.size();
		const int num_num=7;
		vector<vector<vector<double> > > to_output(num_runs, vector<vector<double> >(num_vars + num_num + 1, 
											vector<double>((final_run - first_output + 1)*num_per_yr + 1)));

		for(int r=0; r<num_runs; r++)
			for(int t=0; t<to_output[r][0].size(); t++)
				to_output[r][0][t]=-999*dy;
		string saved_file;
		if(file_map.find("saved_file")==file_map.end()){
			saved_file=output_file;
			saved_file.replace(saved_file.find(".txt"), 4, "_saved.bin");
		}
		else
			saved_file=output_root + file_map["saved_file"];
		
		if(from_map_int("output_header", 0, 2)<2)
			run_final(to_output, output_vars, root, demog_file, pop_file, saved_file);

		output_sim(to_output, output_file, var_names, name, final_output_year);
	}
	//*******************************************************************************

	cout << omp_get_wtime()-begin_omp << " seconds\n";
	
	cout.rdbuf(old_cout);
	cerr.rdbuf(old_cerr);
}
