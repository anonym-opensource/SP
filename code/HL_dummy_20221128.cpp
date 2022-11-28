#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <string>
#include <list>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <typeinfo>
#include <unordered_set>
#include <unordered_map>
#include <climits>
#include <math.h>
#include <thread>
#include <chrono>
#include <shared_mutex>
using namespace std;


/*header files in the Boost library: https://www.boost.org/ */
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

/*some other header files*/
#include <ThreadPool.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_save_read.h>
#include <graph_hash_of_mixed_weighted/HL/two_hop_v2/graph_hash_of_mixed_weighted_CT_v2.h>
#include <text mining/binary_save_read_vector.h>
#include <text mining/binary_save_read_vector_of_vectors.h>
#include <graph_v_of_v_idealID/graph_v_of_v_idealID_shortest_paths.h>

/*querying*/

pair<double, double> querying_element(graph_hash_of_mixed_weighted_CT_v2_case_info* case_info, int s, int t) {

	try {
		auto begin = std::chrono::high_resolution_clock::now();
		CT_extract_distance(*case_info, s, t);
		auto end = std::chrono::high_resolution_clock::now();
		double query_dis_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s


		begin = std::chrono::high_resolution_clock::now();
		vector<pair<int, int>> path;
		CT_extract_path(*case_info, s, t, path);
		end = std::chrono::high_resolution_clock::now();
		double query_path_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

		return { query_dis_time , query_path_time };
	}
	catch (...) {
		cout << "querying_element throw error!" << endl;
		return { 0,0 };
	}
}

pair<double, double> querying(graph_hash_of_mixed_weighted_CT_v2_case_info& case_info, vector<pair<int, int>>& query_list) {

	double query_dis_avg_time = 0, query_path_avg_time = 0;
	int query_times = query_list.size();

	ThreadPool pool(80); // thread_num for querying
	std::vector< std::future<pair<double, double>> > results; // return typename: xxx
	graph_hash_of_mixed_weighted_CT_v2_case_info* case_info_p = &case_info;
	for (int i = 0; i < query_times; i++) {
		int s = query_list[i].first, t = query_list[i].second;
		results.emplace_back(
			pool.enqueue([s, t, case_info_p] { // pass const type value j to thread; [] can be empty
				return querying_element(case_info_p, s, t);
				})
		);
	}
	for (auto&& result : results) {
		pair<double, double> r = result.get();
		query_dis_avg_time += r.first;
		query_path_avg_time += r.second;
	}
	query_dis_avg_time = query_dis_avg_time / (double)query_times;
	query_path_avg_time = query_path_avg_time / (double)query_times;

	return { query_dis_avg_time , query_path_avg_time };
}


/*GST & MDC*/
class GTS_MDC_solve_info {
public:
	double GST_time = 0, MDC_time = 0, GST_HL_time = 0, MDC_HL_time = 0;
};

GTS_MDC_solve_info solve_GST_MDC_element(graph_hash_of_mixed_weighted_CT_v2_case_info* case_info, graph_v_of_v_idealID* ideal_g_p, vector<int>* MDC_query_groups_p) {

	double time_SP, time_GST, time_GST_HL, time_MDC, time_MDC_HL;

	graph_v_of_v_idealID& ideal_g = *ideal_g_p;
	vector<int>& MDC_query_groups = *MDC_query_groups_p;

	/*find SPs*/
	auto begin = std::chrono::high_resolution_clock::now();
	int T_num = MDC_query_groups.size();
	int T_min_ID = 0, T_min_size = INT_MAX;
	vector<vector<float>> distances(T_num);
	vector<vector<int>> predecessors(T_num);
	for (int i = 0; i < T_num; i++) {
		int degree = ideal_g[MDC_query_groups[i]].size();
		if (degree < T_min_size) {
			T_min_size = degree;
			T_min_ID = i;
		}
		graph_v_of_v_idealID_shortest_paths(ideal_g, MDC_query_groups[i], distances[i], predecessors[i]);
	}
	auto end = std::chrono::high_resolution_clock::now();
	time_SP = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	if (T_min_size == 0) {
		GTS_MDC_solve_info x;
		return x;
	}

	/*GST*/
	if (1) {
		begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted solu_min;
		double solu_g_cost_min = 1e5;
		for (int i = 0; i < T_min_size; i++) {
			int v = ideal_g[MDC_query_groups[T_min_ID]][i].first;
			graph_hash_of_mixed_weighted solu_g;
			double solu_g_cost = 0;
			for (int j = 0; j < T_num; j++) {
				if (j != T_min_ID) {
					int s = v;
					int t = predecessors[j][s];
					while (t != s && t != MDC_query_groups[j]) {
						double ec = graph_v_of_v_idealID_edge_weight(ideal_g, s, t);
						graph_hash_of_mixed_weighted_add_edge(solu_g, s, t, ec);
						solu_g_cost += ec;
						s = t;
						t = predecessors[j][s];
					}
				}
			}
			if (solu_g_cost_min > solu_g_cost) {
				solu_g_cost_min = solu_g_cost;
				solu_min = solu_g;
			}
		}
		end = std::chrono::high_resolution_clock::now();
		time_GST = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	/*GST_HL*/
	if (1) {
		begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted solu_min;
		double solu_g_cost_min = 1e5;
		for (int i = 0; i < T_min_size; i++) {
			int v = ideal_g[MDC_query_groups[T_min_ID]][i].first;
			graph_hash_of_mixed_weighted solu_g;
			double solu_g_cost = 0;
			for (int j = 0; j < T_num; j++) {
				if (j != T_min_ID) {
					vector<pair<int, int>> path;
					CT_extract_path(*case_info, v, MDC_query_groups[j], path);
					for (int k = path.size() - 1; k >= 0; k--) {
						int s = path[k].first, t = path[k].second;
						if (s != MDC_query_groups[j] && t != MDC_query_groups[j] && s < INT_MAX && t < INT_MAX) {
							double ec = graph_v_of_v_idealID_edge_weight(ideal_g, s, t);
							graph_hash_of_mixed_weighted_add_edge(solu_g, s, t, ec);
						}
					}
				}
			}
			if (solu_g_cost_min > solu_g_cost) {
				solu_g_cost_min = solu_g_cost;
				solu_min = solu_g;
			}
		}
		end = std::chrono::high_resolution_clock::now();
		time_GST_HL = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	/*MDC*/
	if (1) {
		begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted solu_min;
		double d_min = std::numeric_limits<double>::max();
		int d_min_v = 0;
		for (int i = 0; i < T_min_size; i++) {
			int v = ideal_g[MDC_query_groups[T_min_ID]][i].first;
			double d = 0;
			for (int j = 0; j < T_num; j++) {
				if (j != T_min_ID) {
					double dis = distances[j][v];
					if (d < dis) {
						d = dis;
					}
				}
			}
			if (d_min > d) {
				d_min = d;
				d_min_v = v;
			}
		}
		graph_hash_of_mixed_weighted solu_g;
		for (int j = 0; j < T_num; j++) {
			if (j != T_min_ID) {
				int s = d_min_v;
				int t = predecessors[j][s];
				while (t != s && t != MDC_query_groups[j]) {
					double ec = graph_v_of_v_idealID_edge_weight(ideal_g, s, t);
					graph_hash_of_mixed_weighted_add_edge(solu_g, s, t, ec);
					s = t;
					t = predecessors[j][s];
				}
			}
		}
		end = std::chrono::high_resolution_clock::now();
		time_MDC = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	/*MDC_HL*/
	if (1) {
		begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted solu_min;
		double d_min = std::numeric_limits<double>::max();
		int d_min_v = 0;
		for (int i = 0; i < T_min_size; i++) {
			int v = ideal_g[MDC_query_groups[T_min_ID]][i].first;
			double d = 0;
			for (int j = 0; j < T_num; j++) {
				if (j != T_min_ID) {
					double dis = CT_extract_distance(*case_info, MDC_query_groups[j], v);
					if (d < dis) {
						d = dis;
					}
				}
			}
			if (d_min > d) {
				d_min = d;
				d_min_v = v;
			}
		}
		graph_hash_of_mixed_weighted solu_g;
		for (int j = 0; j < T_num; j++) {
			if (j != T_min_ID) {
				vector<pair<int, int>> path;
				CT_extract_path(*case_info, d_min_v, MDC_query_groups[j], path);
				for (int k = path.size() - 1; k >= 0; k--) {
					int s = path[k].first, t = path[k].second;
					if (s != MDC_query_groups[j] && t != MDC_query_groups[j] && s < INT_MAX && t < INT_MAX) {
						double ec = graph_v_of_v_idealID_edge_weight(ideal_g, s, t);
						graph_hash_of_mixed_weighted_add_edge(solu_g, s, t, ec);
					}
				}
			}
		}
		end = std::chrono::high_resolution_clock::now();
		time_MDC_HL = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	GTS_MDC_solve_info x;
	x.GST_time = time_SP + time_GST;
	x.GST_HL_time = time_GST_HL;
	x.MDC_time = time_SP + time_MDC;
	x.MDC_HL_time = time_MDC_HL;
	return x;
}

GTS_MDC_solve_info solve_GST_MDC(graph_hash_of_mixed_weighted_CT_v2_case_info& case_info, graph_hash_of_mixed_weighted& input_graph, int max_N_ID, vector<vector<int>>& MDC_query_list) {

	graph_v_of_v_idealID ideal_g = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID_2(input_graph, max_N_ID);

	GTS_MDC_solve_info avg;
	int N = MDC_query_list.size();

	ThreadPool pool(50); // thread_num for querying
	std::vector< std::future<GTS_MDC_solve_info> > results; // return typename: xxx
	graph_hash_of_mixed_weighted_CT_v2_case_info* case_info_p = &case_info;
	graph_v_of_v_idealID* ideal_g_p = &ideal_g;
	for (int i = 0; i < N; i++) {
		vector<int>* MDC_query_groups_p = &MDC_query_list[i];
		results.emplace_back(
			pool.enqueue([case_info_p, ideal_g_p, MDC_query_groups_p] { // pass const type value j to thread; [] can be empty
				return solve_GST_MDC_element(case_info_p, ideal_g_p, MDC_query_groups_p);
				})
		);
	}
	for (auto&& result : results) {
		GTS_MDC_solve_info r = result.get();
		avg.GST_time += r.GST_time;
		avg.MDC_time += r.MDC_time;
		avg.GST_HL_time += r.GST_HL_time;
		avg.MDC_HL_time += r.MDC_HL_time;
	}

	avg.GST_time = avg.GST_time / (double)N;
	avg.MDC_time = avg.MDC_time / (double)N;
	avg.GST_HL_time = avg.GST_HL_time / (double)N;
	avg.MDC_HL_time = avg.MDC_HL_time / (double)N;

	return avg;
}



/*exp*/

void exp_element(string data_name, int ec_type, int thread_num, int d, long long int max_bit_size, double max_run_time_seconds) {

	cout << "start indexing " << data_name << " " << ec_type << endl;

	string path = "2022_HL_dummy_data";

	graph_hash_of_mixed_weighted input_graph;
	vector<pair<int, int>> query_list;
	vector<vector<int>> MDC_query_list;
	int max_non_dummy_ID = 0; // the number of non-dummy vertices - 1
	int max_N_ID = 0; // actually should be the number of both dummy and non-dummy vertices
	if (data_name == "amazon") {
		max_non_dummy_ID = 548551;
		max_N_ID = 574511;
		binary_read_vector(path + "//amazon//final//amazon_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//amazon//final//amazon_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//amazon//final//amazon_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//amazon//final//amazon_graph_random_ec.binary");
		}
	}
	else if (data_name == "gplus") {
		max_non_dummy_ID = 107595;
		max_N_ID = 111719;
		binary_read_vector(path + "//gplus//final//gplus_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//gplus//final//gplus_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//gplus//final//gplus_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//gplus//final//gplus_graph_random_ec.binary");
		}
	}
	else if (data_name == "movielens") {
		max_non_dummy_ID = 62422;
		max_N_ID = 62444;
		binary_read_vector(path + "//movielens//final//movielens_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//movielens//final//movielens_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//movielens//final//movielens_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//movielens//final//movielens_graph_random_ec.binary");
		}
	}
	else if (data_name == "dblp") {
		max_non_dummy_ID = 1248890;
		max_N_ID = 1372828;
		binary_read_vector(path + "//dblp//final//dblp_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//dblp//final//dblp_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//dblp//final//dblp_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//dblp//final//dblp_graph_random_ec.binary");
		}
	}
	else if (data_name == "reddit") {
		max_non_dummy_ID = 4262833;
		max_N_ID = 5409492;
		binary_read_vector(path + "//reddit//final//reddit_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//reddit//final//reddit_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//reddit//final//reddit_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//reddit//final//reddit_graph_random_ec.binary");
		}
	}
	else if (data_name == "musae") {
		max_non_dummy_ID = 19108;
		max_N_ID = 32293;
		binary_read_vector(path + "//musae//final//musae_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//musae//final//musae_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//musae//final//musae_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//musae//final//musae_graph_random_ec.binary");
		}
	}
	else if (data_name == "deezer") {
		max_non_dummy_ID = 143883;
		max_N_ID = 143973;
		binary_read_vector(path + "//deezer//final//deezer_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//deezer//final//deezer_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//deezer//final//deezer_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//deezer//final//deezer_graph_random_ec.binary");
		}
	}
	else if (data_name == "Pokec") {
		max_non_dummy_ID = 1139631;
		max_N_ID = 1170890;
		binary_read_vector(path + "//Pokec//final//Pokec_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//Pokec//final//Pokec_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//Pokec//final//Pokec_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//Pokec//final//Pokec_graph_random_ec.binary");
		}
	}
	else if (data_name == "twitch") {
		max_non_dummy_ID = 34117;
		max_N_ID = 37282;
		binary_read_vector(path + "//twitch//final//twitch_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//twitch//final//twitch_MDC_query_list.binary", MDC_query_list);
		if (ec_type == 0) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//twitch//final//twitch_graph_Jacard_ec.binary");
		}
		else if (ec_type == 1) {
			input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//twitch//final//twitch_graph_random_ec.binary");
		}
	}
	else if (data_name == "github") {
	max_non_dummy_ID = 37699;
	max_N_ID = 41706;
	binary_read_vector(path + "//github//final//github_query_list.binary", query_list);
	binary_read_vector_of_vectors(path + "//github//final//github_MDC_query_list.binary", MDC_query_list);
	if (ec_type == 0) {
		input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//github//final//github_graph_Jacard_ec.binary");
	}
	else if (ec_type == 1) {
		input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//github//final//github_graph_random_ec.binary");
	}
	}
	else if (data_name == "example") {
		max_non_dummy_ID = 8e3;
		max_N_ID = 1e4 + 1;
		binary_read_vector(path + "//example_query_list.binary", query_list);
		binary_read_vector_of_vectors(path + "//example_MDC_query_list.binary", MDC_query_list);
		input_graph = graph_hash_of_mixed_weighted_binary_read(path + "//example_graph.binary");
	}

	/*output*/
	ofstream outputFile;
	outputFile.precision(8);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	string save_name = "exp_" + data_name + "_ectype_" + to_string(ec_type) + "_T_" + to_string(thread_num) + "_d_" + to_string(d);
	outputFile.open(save_name + ".csv");
	outputFile << "data,ec_type,thread_num,d,"
		<< "CT-PSL*_MB,CT-PSL*_time,CT-PSL*_query_dis_time,CT-PSL*_query_path_time,"
		<< "CT-Cano-PSL*_MB,CT-Cano-PSL*_time,CT-Cano-PSL*_query_dis_time,CT-Cano-PSL*_query_path_time,CT-Cano-PSL*_cano_time,CT-Cano-PSL*_cano_ratio,"
		<< "CT-Cano-PLL_MB,CT-Cano-PLL_time,CT-Cano-PLL_query_dis_time,CT-Cano-PLL_query_path_time,CT-Cano-PLL_cano_time,CT-Cano-PLL_cano_ratio,"
		<< "CT-Cano-partial-PLL_MB,CT-Cano-partial-PLL_time,CT-Cano-partial-PLL_query_dis_time,CT-Cano-partial-PLL_query_path_time,CT-Cano-partial-PLL_cano_time,CT-Cano-partial-PLL_cano_ratio,"
		<< "GST_time,GST_HL_time,MDC_time,MDC_HL_time" << endl;

	outputFile << data_name << "," << ec_type << "," << thread_num << "," << d << "," << std::flush;

	bool PSL_reach_limits = false;
	string reach_limits_type = "reach_limits_type unknown";
	/* CT-PSL* */
	if (1) {
		cout << "start CT-PSL*" << endl;
		graph_hash_of_mixed_weighted_CT_v2_case_info case_info;
		case_info.two_hop_case_info.use_2019R1 = 1;
		case_info.two_hop_case_info.use_2019R2 = 1;
		case_info.two_hop_case_info.use_canonical_repair = 0;
		case_info.d = d;
		case_info.use_PLL = 0;
		case_info.thread_num = thread_num;
		case_info.max_bit_size = max_bit_size;
		case_info.max_run_time_seconds = max_run_time_seconds;

		bool catch_error = false;
		try {
			CT_v2(input_graph, max_N_ID, case_info);
		}
		catch (string s) {
			PSL_reach_limits = true;
			cout << s << endl;
			reach_limits_type = s;
			clear_gloval_values_CT();
			catch_error = true;
		}
		cout << "finish CT-PSL*" << endl;
		if (catch_error) {
			if (reach_limits_type == reach_limit_error_string_MB) {
				outputFile << "0,0,0,0," << std::flush;
			}
			else if (reach_limits_type == reach_limit_error_string_time) {
				outputFile << "-1,-1,-1,-1," << std::flush;
			}
			else {
				outputFile << "1e5,1e5,1e5,1e5," << std::flush;
			}
		}
		else {
			pair<double, double> query_results = querying(case_info, query_list);

			outputFile << (double)case_info.compute_label_bit_size() / 1024 / 1024 << "," << case_info.time_total << "," <<
				query_results.first << "," << query_results.second << "," << std::flush;
		}
		case_info.record_all_details(save_name + "_CT-PSL*");
		case_info.clear_labels();
	}

	/* CT-Cano-PSL* */
	if (1) {
		cout << "start CT-Cano-PSL*" << endl;
		graph_hash_of_mixed_weighted_CT_v2_case_info case_info;
		case_info.two_hop_case_info.use_2019R1 = 1;
		case_info.two_hop_case_info.use_2019R2 = 1;
		case_info.two_hop_case_info.use_canonical_repair = 1;
		case_info.d = d;
		case_info.use_PLL = 0;
		case_info.thread_num = thread_num;
		case_info.max_bit_size = max_bit_size;
		case_info.max_run_time_seconds = max_run_time_seconds;

		bool catch_error = false;
		try {
			if (PSL_reach_limits) {
				throw(reach_limits_type);
			}
			CT_v2(input_graph, max_N_ID, case_info);
		}
		catch (string s) {
			cout << s << endl;
			reach_limits_type = s;
			clear_gloval_values_CT();
			catch_error = true;
		}
		cout << "finish CT-Cano-PSL*" << endl;
		if (catch_error) {
			if (reach_limits_type == reach_limit_error_string_MB) {
				outputFile << "0,0,0,0,0,0," << std::flush;
			}
			else if (reach_limits_type == reach_limit_error_string_time) {
				outputFile << "-1,-1,-1,-1,-1,-1," << std::flush;
			}
			else {
				outputFile << "1e5,1e5,1e5,1e5,1e5,1e5," << std::flush;
			}
		}
		else {
			pair<double, double> query_results = querying(case_info, query_list);

			outputFile << (double)case_info.compute_label_bit_size() / 1024 / 1024 << "," << case_info.time_total << "," <<
				query_results.first << "," << query_results.second << "," << case_info.two_hop_case_info.time_canonical_repair1 + case_info.two_hop_case_info.time_canonical_repair2
				<< "," << case_info.two_hop_case_info.canonical_repair_remove_label_ratio << "," << std::flush;
		}
		case_info.record_all_details(save_name + "_CT-Cano-PSL*");
		case_info.clear_labels();
	}

	/* CT-Cano-PLL */
	if (1) {
		cout << "start CT-Cano-PLL" << endl;
		graph_hash_of_mixed_weighted_CT_v2_case_info case_info;
		case_info.two_hop_case_info.use_2019R1 = 1;
		case_info.two_hop_case_info.use_2019R2 = 1;
		case_info.two_hop_case_info.use_canonical_repair = 1;
		case_info.d = d;
		case_info.use_PLL = 1;
		case_info.thread_num = thread_num;
		case_info.max_bit_size = max_bit_size;
		case_info.max_run_time_seconds = max_run_time_seconds;

		bool catch_error = false;
		try {
			CT_v2(input_graph, max_N_ID, case_info);
		}
		catch (string s) {
			cout << s << endl;
			reach_limits_type = s;
			clear_gloval_values_CT();
			catch_error = true;
		}
		cout << "finish CT-Cano-PLL" << endl;
		if (catch_error) {
			if (reach_limits_type == reach_limit_error_string_MB) {
				outputFile << "0,0,0,0,0,0," << std::flush;
			}
			else if (reach_limits_type == reach_limit_error_string_time) {
				outputFile << "-1,-1,-1,-1,-1,-1," << std::flush;
			}
			else {
				outputFile << "1e5,1e5,1e5,1e5,1e5,1e5," << std::flush;
			}
		}
		else {
			pair<double, double> query_results = querying(case_info, query_list);

			outputFile << (double)case_info.compute_label_bit_size() / 1024 / 1024 << "," << case_info.time_total << "," <<
				query_results.first << "," << query_results.second << "," << case_info.two_hop_case_info.time_canonical_repair1 + case_info.two_hop_case_info.time_canonical_repair2
				<< "," << case_info.two_hop_case_info.canonical_repair_remove_label_ratio << "," << std::flush;
		}
		case_info.record_all_details(save_name + "_CT-Cano-PLL");
		case_info.clear_labels();
	}

	/* CT-Cano-partial-PLL */
	if (1) {
		cout << "start CT-Cano-partial-PLL" << endl;
		graph_hash_of_mixed_weighted_CT_v2_case_info case_info;
		case_info.two_hop_case_info.use_2019R1 = 1;
		case_info.two_hop_case_info.use_2019R2 = 1;
		case_info.two_hop_case_info.use_canonical_repair = 1;
		case_info.d = d;
		case_info.use_PLL = 1;
		case_info.two_hop_case_info.use_dummy_dij_search_in_PLL = 1;
		case_info.max_non_dummy_ID = max_non_dummy_ID;
		case_info.thread_num = thread_num;
		case_info.max_bit_size = max_bit_size;
		case_info.max_run_time_seconds = max_run_time_seconds;

		bool catch_error = false;
		try {
			CT_v2(input_graph, max_N_ID, case_info);
		}
		catch (string s) {
			cout << s << endl;
			reach_limits_type = s;
			clear_gloval_values_CT();
			catch_error = true;
		}
		cout << "finish CT-Cano-partial-PLL" << endl;
		if (catch_error) {
			if (reach_limits_type == reach_limit_error_string_MB) {
				outputFile << "0,0,0,0,0,0,0,0,0,0," << std::flush;
			}
			else if (reach_limits_type == reach_limit_error_string_time) {
				outputFile << "-1,-1,-1,-1,-1,-1,-1,-1,-1,-1," << std::flush;
			}
			else {
				outputFile << "1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5,1e5," << std::flush;
			}
		}
		else {
			pair<double, double> query_results = querying(case_info, query_list);

			GTS_MDC_solve_info x = solve_GST_MDC(case_info, input_graph, max_N_ID, MDC_query_list);

			outputFile << (double)case_info.compute_label_bit_size() / 1024 / 1024 << "," << case_info.time_total << "," <<
				query_results.first << "," << query_results.second << "," << case_info.two_hop_case_info.time_canonical_repair1 + case_info.two_hop_case_info.time_canonical_repair2
				<< "," << case_info.two_hop_case_info.canonical_repair_remove_label_ratio
				<< "," << x.GST_time << "," << x.GST_HL_time << "," << x.MDC_time << "," << x.MDC_HL_time << "," << std::endl;
		}
		case_info.record_all_details(save_name + "_CT-Cano-partial-PLL");
		case_info.clear_labels();
	}

}

void HL_main_exp() {

	vector<string> used_datas = { "musae", "twitch", "github", "amazon", "reddit", "dblp" };

	long long int max_bit_size = pow(1024, 3) * 300; 
	double max_run_time_seconds = 3600 * 24; 
	int thread_num = 80;
	int d = 50; // querying is too slow when d is large
	 
	/*Jacard & random*/
	if (1) {
		for (int i = 0; i < used_datas.size(); i++) {
			exp_element(used_datas[i], 0, thread_num, d, max_bit_size, max_run_time_seconds);
			exp_element(used_datas[i], 1, thread_num, d, max_bit_size, max_run_time_seconds);
		}
	}
}







int main()
{
	cout << "Start running..." << endl;
	auto begin = std::chrono::high_resolution_clock::now();
	/*the two values below are for #include <graph_hash_of_mixed_weighted.h>*/
	graph_hash_of_mixed_weighted_turn_on_value = 1e3;
	graph_hash_of_mixed_weighted_turn_off_value = 1e1;

	HL_main_exp();


	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	cout << "END    runningtime: " << runningtime << "s" << endl;
}

