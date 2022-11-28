# Indexing Shortest Paths Between Vertices and Vertex Groups in Weighted Knowledge Graphs

The introduction of these files are as follows. 


# Datasets

There are six datasets: Musae, Twitch, Github, Amazon, DBLP, Reddit (zipped together; should be unzipped first). 

There are sie files for each dataset. For example, the Amazon dataset contains the following files: 
1. amazon_edges.txt: the edges, as well as random and Jacard edge weights.
2. amazon_vertex_names.txt: the vertices as well as vertex names.
3. amazon_graph_Jacard_ec.binary: the binary file of an adjacency list of the Amazon graph with Jacard weights. This file will be read by the following codes.
4. amazon_graph_random_ec.binary: the binary file of an adjacency list of the Amazon graph with random weights. This file will be read by the following codes.
5. amazon_MDC_query_list.binary: the binary file of the randomly generated query list for the indexing algorithms in the main experiments (each query is a pair of dummy and non-dummy vertices). This file will be read by the following codes.
6. amazon_query_list.binary: the binary file of the randomly generated query list for the GST and MDC algorithms in the main experiments (each query is a set of dummy vertices, respresenting the set of vertex groups in the GST and MDC problems). This file will be read by the following codes.





# C++ codes 

The C++ source codes are in <b>HL_dummy_202*****.cpp</b>. 

It is recommended to fold all the regions of codes for easy reading (by pressing Ctrl+M+O in Visual Studio). 

Running these codes requires some header files (e.g. #include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>; see cppheader_202*****.zip) and the Boost library: https://www.boost.org/ (e.g. #include <boost/random.hpp>) . 

We use the above datasets to generate the above binary files for the experiments. The reason for generating binary files is that it is much faster to read binary files than to read raw data files. 

After making the header files and binary files ready, <b>all the experiments in our paper can be conducted by running the function "HL_main_exp()"</b>, Make sure there is enough memory (500 GB memory).

To read these C++ codes in detail, it is recommended to start from "HL_main_exp()", and then go to "exp_element()". More detailed codes in other regions can then be traced, such as the port of proposed algorithms in "graph_hash_of_mixed_weighted/HL/two_hop_v2/graph_hash_of_mixed_weighted_CT_v2.h"


