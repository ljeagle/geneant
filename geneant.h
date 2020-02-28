#include <string>
#include <vector>

using std::vector;
using std::string;

int toIdx(char c);
vector<string> pairwise(string s1, string s2);
vector< vector<string> > pairwise_batch(string query_seq, vector<string> target_seqs);
