#ifndef _DEGREE_H_
#define _DEGREE_H_

#include<vector>
#include"gurobi_c++.h"
using namespace std;

class Degree {
	int bi[31][32];
public:
	GRBEnv env;
	float max_time = 0;
	float min_time = 100000;

	Degree();

	int upperBoundUnivariate(int r, int n, int d, vector<int>& k);
	void computeT(vector<bool>& flag, vector<int>& k, int r, int n, bool& isFull, int curLayer, int left,
		int& totalLayer, vector<int>& b);
	int buildModelUnivariate(vector<vector<bool>>& A, vector<vector<bool>>& B, int r, int n, int d);
	void initializeVector(GRBModel& model, vector<GRBVar>& v, int cols);
	void initializeVector(GRBModel& model, vector<vector<GRBVar>>& v, int rows, int cols);
	void modAddition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb,
		vector<GRBVar>& c, vector<GRBVar>& c1, vector<GRBVar>& q, vector<GRBVar>& d, int n);
	void addition(GRBModel& model, vector<GRBVar>& a, vector<GRBVar>& b, vector<GRBVar>& c, vector<GRBVar>& s,
		int la, int sumLen);
	int solveGeneralOP(vector<int>& Z, int pi, int r, int n);
	int MaxNum(vector<int>& k, int r, int n, int d);
	bool isCondition(vector<int>& k, int r, int n, int d);
	bool checkEquivalence(vector<int>& k, int n);
	void iterativeEnu(vector<int>& k, int curLayer, int totalLayer, int r1, int r2, int n, int d,
		int& validNum);
	void outputTime();
	void searchParameters(int r1, int r2, int n, int d, int w);


};



#endif