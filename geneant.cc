#include <stdlib.h>
#include <sys/time.h>
#include <fstream>
#include <iostream>

#include "geneant.h"

using namespace std;

int MIN_INT = -1000000000;

int mx(int a, int b){
	return a>b? a: b;
}

int fsb[5][5] = {{5,-1,-2,-1,-3}, {-1,5,-3,-2,-4}, {-2,-3,5,-2,-2}, {-1,-2,-2,5,-1}, {-3,-4,-2,-1,MIN_INT}};

int toIdx(char c){
	switch(c){
	case 'a': return 0;
	case 'c': return 1;
	case 'g': return 2;
	case 't': return 3;
	default: return 4;
	}
}

vector<string> pairwise(string s1, string s2) {
	bool exchanged = false;
	if (s2.length() < s1.length())
	{
		string tmp = s1;
		s1 = s2;
		s2 = tmp;
		exchanged = true;
	}

	int len1 = s1.length(), len2 = s2.length();
	int WINDOW_SIZE = 150 * 2, MAX_DISTANCE = WINDOW_SIZE / 2;

	int **fen = new int *[len1 + 1];  
	for(int i=0; i<len1 + 1; ++i)  
	{
		fen[i] = new int[len2 + 1];  
	}  

	string alignment_rule[WINDOW_SIZE + 1];
	for (int i = 0; i <= WINDOW_SIZE; i++)  alignment_rule[i].reserve(len1 + len2);

	alignment_rule[MAX_DISTANCE] = "";
	for (int i = MAX_DISTANCE + 1; i <= WINDOW_SIZE; i++)
	{
		alignment_rule[i] = alignment_rule[i - 1] + "+";
	} 

	string aux_alignment_rule[WINDOW_SIZE + 1];

	for (int i = 0; i <= WINDOW_SIZE; i++)  aux_alignment_rule[i].reserve(len1 + len2);
	aux_alignment_rule[MAX_DISTANCE] = "";
	for (int i = MAX_DISTANCE + 1; i <= WINDOW_SIZE; i++)
	{
		aux_alignment_rule[i] = aux_alignment_rule[i - 1] + "-";
	} 

	string reg_str = ""; 

	for (int i = 0; i <= len1; i++)  {
		for (int j = 0; j <= len2; j++)  {
			if (i == 0 && j == 0){
				fen[i][j] = 0;
				continue;
			}
			if (abs(i - j) > MAX_DISTANCE)   
			{ 
				if (i == len1 && j > i)  
				{ 
					reg_str = alignment_rule[WINDOW_SIZE] + string(len2 - (i + MAX_DISTANCE), '+');  
					break;
				}  
				continue; 
			}
			int temp = MIN_INT, gap1 = MIN_INT, gap2 = MIN_INT, gap3 = MIN_INT;
			bool ret = false;
			if (i > 0) {
				gap1 = fen[i-1][j] + fsb[toIdx(s1[i-1])][4];
				temp = mx(temp, gap1);
			}
			if (j > 0) 
			{
				gap2 = fen[i][j-1] + fsb[4][toIdx(s2[j-1])];
				temp = mx(temp, gap2);
			}
			if (i > 0 && j > 0) 
			{
				int value = fsb[toIdx(s1[i-1])][toIdx(s2[j-1])];
				if (value == 5)  { ret = true;}
				gap3 = fen[i-1][j-1] + value;
				temp = mx(temp, gap3);
			}

			fen[i][j] = temp;


			if (i == 0)  continue;

			int pos = j - i + MAX_DISTANCE;
			if (j == 0)   { alignment_rule[pos] = aux_alignment_rule[i + MAX_DISTANCE]; continue; }
			if (temp == gap1)
			{
				if (pos < WINDOW_SIZE)  alignment_rule[pos] = alignment_rule[pos + 1] + "-"; 
			}  else if (temp == gap2)
			{
				if (pos > 0)   alignment_rule[pos] = alignment_rule[pos - 1] + "+"; 
			}  else if (temp == gap3)
			{
				if (ret)  alignment_rule[pos] = alignment_rule[pos] + "=";
				else   alignment_rule[pos] = alignment_rule[pos] + "~";
			}
			if (i == len1 && j == len2)  
			{
				reg_str = alignment_rule[pos];
			}

		}
	}
	vector<string> result(2);
	for (int i = 0; i < 2; i++)
	{
		string tmp = "";
		for (int j = 0, k = 0; j < reg_str.length(); j++)
		{
			if (reg_str[j] == '+')
			{
				if (i == 0) tmp = tmp + "*";
				else   tmp = tmp + s2[k++]; 
			}  else if (reg_str[j] == '-')
			{
				if (i == 0)  tmp = tmp + s1[k++];
				else  tmp = tmp + "*";
			}  else
			{
				if (i == 0)  tmp = tmp + s1[k++];
				else  tmp = tmp + s2[k++];
			}
		}
		if(exchanged) {
			result[1-i] = tmp;
		}
		else {
			result[i] = tmp;
		}
	}

	int count = 0;
	for (int i = 0; i < reg_str.length(); i++)
	{
		if (reg_str[i] == '=') count++;
	}
	double similarity = (double)count/(double)reg_str.length();

	for (int i = 0; i < len1 + 1; ++i)
	{
		delete fen[i];
	}
	return result;
}


vector< vector<string> > pairwise_batch(string query_seq, vector<string> target_seqs) {
	vector< vector<string> > aligned_seqs(target_seqs.size());
#pragma omp parallel for
	for(size_t i = 0; i < target_seqs.size(); i++) {
		aligned_seqs[i] = pairwise(query_seq, target_seqs[i]);
	} 
	return aligned_seqs;
}
