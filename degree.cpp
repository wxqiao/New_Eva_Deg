#include"degree.h"
#include<vector>
#include<set>
#include<algorithm>
#include<ctime>
using namespace std;

//compute T_{r,w}^R and construct flag (flag[j]=1<==>j in T_{r,w}^R) (r>=2)
// let flag=[0,...,0],flag.size=n; k=[h_1,...,h_w]; isFull=false; curLayer=0;left=r-2;totalLayer=w
void Degree::computeT(vector<bool>& flag, vector<int>& k, int r, int n, bool& isFull, int curLayer, int left,
	int& totalLayer, vector<int>& b) {
	if (curLayer == totalLayer - 1) {
		b.push_back(left);
		int sum = 0;
		int pos = 0;
		for (int i = 0; i < totalLayer; i++) {
			sum += b[i];
			pos = (pos + k[i] * b[i]) % n;
		}
		if (sum != r - 2) {
			cout << "program errors" << endl;
		}
		flag[pos] = 1;
		b.pop_back();
		isFull = true;
		for (int i = 0; i < n; i++) {
			if (flag[i] == false) {
				isFull = false;
				break;
			}
		}
	}
	else {
		for (int i = 0; i <= left; i++) {
			b.push_back(i);
			computeT(flag, k, r, n, isFull, curLayer + 1, left - i, totalLayer, b);
			b.pop_back();
			if (isFull) {
				break;
			}
		}
	}
}

Degree::Degree() {
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 6);

	int row = 31, col = 32;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			bi[i][j] = 0;
		}
	}
	bi[0][0] = 1;
	for (int i = 1; i < row; i++) {
		bi[i][0] = 1;
		for (int j = 1; j < i + 1; j++) {
			bi[i][j] = bi[i - 1][j - 1] + bi[i - 1][j];
		}
	}
}

//Model for the maximal number of different elements
int Degree::buildModelUnivariate(vector<vector<bool>>& A, vector<vector<bool>>& B, int r, int n, int d) {
	//clock_t t;
	//t = clock();

	vector<vector<GRBVar>> Aup, Adw, Bup, Bdw;
	vector<GRBVar> X;
	GRBModel model = GRBModel(env);
	initializeVector(model, X, n);
	initializeVector(model, Aup, r, n);
	initializeVector(model, Adw, r, n);
	initializeVector(model, Bup, r, n);
	initializeVector(model, Bdw, r, n);

	for (int i = 1; i < r; i++) {
		for (int j = 0; j < n; j++) {
			model.addConstr(Aup[i][j] == Adw[i - 1][(j + d) % n]);
			model.addConstr(Bup[i][j] == Bdw[i - 1][(j + d) % n]);
		}
	}

	for (int i = 0; i < n; i++) {
		model.addConstr(Aup[0][i] == 0);
		model.addConstr(Adw[r - 1][i] == 0);
		model.addConstr(Bup[0][i] == 0);
		model.addConstr(Bdw[r - 1][i] == 0);
	}

	for (int i = 0; i < r; i++) {
		for (int j = 0; j < n; j++) {
			if (A[i][j] == 0) {
				model.addConstr(Aup[i][j] == 0);
				model.addConstr(Adw[i][j] == 0);
			}
			if (B[i][j] == 0) {
				model.addConstr(Bup[i][j] == 0);
				model.addConstr(Bdw[i][j] == 0);
			}
		}
	}

	for (int i = 0; i < r-1; i++) {
		GRBLinExpr sumA = 0;
		GRBLinExpr sumB = 0;
		for (int j = 0; j < n; j++) {
			sumA = sumA + Adw[i][j];
			sumB = sumB + Bdw[i][j];
		}
		model.addConstr(sumA <= bi[r - 2][i]);
		model.addConstr(sumB <= bi[r - 2][i]);
	}

	for (int i = 0; i<n; i++) {
		GRBLinExpr sum = 0;
		for (int j = 0; j < r; j++) {
			model.addConstr(X[i] >= Aup[j][i]);
			model.addConstr(X[i] >= Adw[j][i]);
			model.addConstr(X[i] >= Bup[j][i]);
			model.addConstr(X[i] >= Bdw[j][i]);
			sum = sum + Aup[j][i] + Adw[j][i] + Bup[j][i] + Bdw[j][i];
		}
		model.addConstr(X[i] <= sum);
	}

	GRBLinExpr obj = 0;
	for (int i = 0; i < n; i++) {
		obj = obj + X[i];
	}

	model.setObjective(obj, GRB_MAXIMIZE);
	model.optimize();

	//t = clock() - t;

	if (model.get(GRB_IntAttr_Status) == 3) {
		return -2;
	}

	/*float ft = ((float)t) / CLOCKS_PER_SEC;
	if (ft > max_time)
		max_time = ft;
	if (ft < min_time)
		min_time = ft;*/

	return model.getObjective().getValue();
}

void Degree::initializeVector(GRBModel& model, vector<GRBVar>& v, int cols) {
	v.clear();
	v.resize(cols);
	for (int i = 0; i < cols; i++)
		v[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
}

void Degree::initializeVector(GRBModel& model, vector<vector<GRBVar>>& v, int rows, int cols) {
	v.clear();
	v.resize(rows);
	for (int i = 0; i < rows; i++) {
		v[i].resize(cols);
		for (int j = 0; j < cols; j++) {
			v[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
		}
	}
}


void Degree::modAddition(GRBModel& model, vector<GRBVar>& a, int sa, vector<GRBVar>& b, int sb, 
	vector<GRBVar>& c, vector<GRBVar>& c1, vector<GRBVar>& q, vector<GRBVar>& d, int n) {
	model.addConstr(c[0] == 0);
	for (int i = 0; i < n; i++) {
		model.addConstr(a[(i - sa + n) % n] + b[(i - sb + n) % n] + c[i] == q[i] + 2 * c[i + 1]);
	}
	model.addConstr(c1[0] == c[n]);
	for (int i = 0; i < n; i++) {
		model.addConstr(2 * c1[i + 1] + d[i] == q[i] + c1[i]);
	}
}

void Degree::addition(GRBModel& model, vector<GRBVar>& a, vector<GRBVar>& b, vector<GRBVar>& c, vector<GRBVar>& s,
	int la, int sumLen) {
	model.addConstr(c[0] == 0);
	/*if (n < sumLen) {
		cout << "errors! Please modify addition()!" << endl;
	}*/
	for (int i = 0; i < la; i++) {
		model.addConstr(a[i] + b[i] + c[i] == s[i] + 2 * c[i + 1]);
	}
	for (int i = la; i < sumLen; i++) {
		model.addConstr(b[i] + c[i] == s[i] + 2 * c[i + 1]);
	}
}

//solve the model for upperbound
int Degree::solveGeneralOP(vector<int>& Z, int pi, int r, int n) {
	//clock_t  t;
	//t = clock();

	GRBModel model = GRBModel(env);
	int k = Z.size();

	int sumLen = ceil(log2f(n));
	sumLen = sumLen + r + 1;//we need to sum up at most n variables smaller than 2^r.

	vector<vector<GRBVar> > a;
	vector<vector<GRBVar> > ms;//modular sum
	vector<vector<GRBVar> > cs;//common sum
	vector<vector<GRBVar> > mcarry;//carry for the modular addition
	vector<vector<GRBVar> > ccarry;//carry for the common addition
	vector<GRBVar> u;//use to denote the choice of pi positions from Z

	initializeVector(model, a, k, n);
	initializeVector(model, ms, k * 2 + 1, n);
	initializeVector(model, cs, k + 1, sumLen);
	initializeVector(model, mcarry, k * 2, n + 1);
	initializeVector(model, ccarry, k, sumLen + 1);
	initializeVector(model, u, k);

	for (int i = n - 1; i >= r + 1; i--) {
		for (int j = 0; j < k; j++) {
			model.addConstr(a[j][i] == 0);
		}
	}

	for (int i = 0; i < n; i++) {
		model.addConstr(ms[0][i] == 0);
	}

	for (int i = 0; i < sumLen; i++) {
		model.addConstr(cs[0][i] == 0);
	}

	for (int i = 0; i < k; i++) {
		modAddition(model, ms[i * 2], 0, a[i], Z[i], mcarry[i * 2], mcarry[i * 2 + 1], ms[i * 2 + 1], ms[i * 2 + 2], n);
		//addition(model, cs[i], a[i], ccarry[i], cs[i + 1], r+1, sumLen);
		addition(model, a[i], cs[i], ccarry[i], cs[i + 1], r + 1, sumLen);
	}

	//extra condition on the common addition sum
	for (int i = sumLen - 1; i >= r + 1; i--) {
		model.addConstr(cs[k][i] == 0);
	}
	for (int i = 0; i < r; i++) {
		model.addConstr(cs[k][i] + cs[k][r] <= 1);
	}

	//extra condition on the choice
	GRBLinExpr uSum = 0;
	for (int i = 0; i < k; i++) {
		uSum = uSum + u[i];
		GRBLinExpr sum = 0;
		for (int j = 0; j < n; j++) {
			sum = sum + a[i][j];
			model.addConstr(u[i] >= a[i][j]);
		}
		model.addConstr(u[i] <= sum);
	}
	model.addConstr(uSum <= pi);

	GRBLinExpr obj = 0;
	for (int i = 0; i < n; i++) {
		obj = obj + ms[k * 2][i];
	}

	model.setObjective(obj, GRB_MAXIMIZE);
	model.optimize();

	//t = clock() - t;

	/*float ft = ((float)t) / CLOCKS_PER_SEC;
	if (ft > max_time)
		max_time = ft;
	if (ft < min_time)
		min_time = ft;*/

	int res = model.getObjective().getValue();

	return res;
}




//upperbound for the output of the r-th round (r>=2)
int Degree::upperBoundUnivariate(int r, int n, int d, vector<int>& k) {
	clock_t st;
	st = clock();

	int w = k.size();

	int max = 0;

	vector<vector<bool>> A(r);
	for (int s = 0; s < A.size(); s++)
		A[s].resize(n);
	vector<vector<bool>> B(r);
	for (int t = 0; t < B.size(); t++)
		B[t].resize(n);

	vector<bool> flag;
	flag.clear();
	for (int i = 0; i < n; i++)
		flag.push_back(0);
	bool isFull = false;
	vector<int> b;
	b.clear();
	computeT(flag, k, r, n, isFull, 0, r - 2, w, b);

	for (int i = 0; i < w; i++) {	
		//construct A(r) for classA of k[i]
		for (int s = 0; s < r; s++) {
			for (int t = 0; t < n; t++) {
				A[s][(t + (r - s) * d + k[i]) % n] = flag[t];
			}
		}

		for (int j = 0; j < w; j++) {
			//construct B(r) for classB of k[j]
			for (int s = 0; s < r; s++) {
				for (int t = 0; t < n; t++) {
					B[s][(t + (r -1- s) * d + k[j]) % n] = flag[t];
				}
			}
			//construct Z for A and B
			vector<int> Z;
			Z.clear();
			for (int s = 0; s < n; s++) {
				bool isFind = false;
				for (int t = 0; t < r; t++) {
					if (A[t][s] || B[t][s]) {
						isFind = true;
						break;
					}

				}
				if (isFind) {
					Z.push_back(s);
				}
			}

			int pi = buildModelUnivariate(A, B, r, n, d);

			int res = solveGeneralOP(Z, pi, r, n);

			if (res > max) {
				max = res;
			}
			
		}
	}

	st = clock() - st;
	float ft = ((float)st) / CLOCKS_PER_SEC;
	cout<<" - solving time:" << ft<<"s";

	return max;
}



//Compute pi (the maximal number of different elements in Class)
int Degree::MaxNum(vector<int>& k, int r, int n, int d) {
	clock_t ct;
	ct = clock();

	int w = k.size();

	int max = 0;

	vector<vector<bool>> A(r);
	for (int s = 0; s < A.size(); s++)
		A[s].resize(n);
	vector<vector<bool>> B(r);
	for (int t = 0; t < B.size(); t++)
		B[t].resize(n);

	vector<bool> flag;
	flag.clear();
	for (int i = 0; i < n; i++)
		flag.push_back(0);
	bool isFull = false;
	vector<int> b;
	b.clear();
	computeT(flag, k, r, n, isFull, 0, r - 2, w, b);

	for (int i = 0; i < w; i++) {
		//construct A(r) for classA of k[i]
		for (int s = 0; s < r; s++) {
			for (int t = 0; t < n; t++) {
				A[s][(t + (r - s) * d + k[i]) % n] = flag[t];
			}
		}

		for (int j = 0; j < w; j++) {
			//construct B(r) for classB of k[j]
			for (int s = 0; s < r; s++) {
				for (int t = 0; t < n; t++) {
					B[s][(t + (r - 1 - s) * d + k[j]) % n] = flag[t];
				}
			}

			int pi = buildModelUnivariate(A, B, r, n, d);

			if (pi > max) {
				max = pi;
			}

		}
	}

	ct = clock() - ct;
	float ft = ((float)ct) / CLOCKS_PER_SEC;

	if (ft > max_time)
		max_time = ft;
	if (ft < min_time)
		min_time = ft;

	return max;
	
}


//test if the condition holds at r-th round
bool Degree::isCondition(vector<int>& k, int r, int n, int d) {
	int e = log2(n);
	int max = n;
	if (e >= r)
		max = 1 << r;

	int pi = MaxNum(k, r, n, d);

	if (pi == max) {
		return 1;
	}
	return 0;

}


// true:ÒÑ´æÔÚ
bool Degree::checkEquivalence(vector<int>& k, int n) {
	vector<int> eq(k.size());
	
	for (int i = 1; i < k.size(); i++) {
		for (int j = 0; j < k.size(); j++) {
			eq[j] = (k[j] - k[i] + n) % n;
		}
		sort(eq.begin(), eq.end());
		bool isFind = false;
		for (int i = 1; i < k.size(); i++) {
			if (eq[i] == k[i]) {
				continue;
			}
			else if (eq[i] > k[i]) {
				isFind = false;
				break;
			}
			else {
				isFind = true;
				break;
			}
		}
		if (isFind)
			return true;
	}
	return false;
}

void Degree::iterativeEnu(vector<int>& k, int curLayer, int totalLayer, int r1, int r2, int n, int d, 
	int& validNum) {
	if (curLayer == totalLayer) {
		
		bool isEqui = checkEquivalence(k, n);
		if (!isEqui) {
			

			bool isValid_r1 = isCondition(k, r1, n, d);

			if (isValid_r1) {
				bool isValid_r2 = isCondition(k, r2, n, d);

				if (isValid_r2) {
					validNum++;

					for (int j = 0; j < k.size(); j++) {
						cout << k[j] << " ";

					}

					cout << ": ";
					cout << "The conditions hold!    ";
					outputTime();

					if (validNum > 61) {
						system("pause");
					}

				}

			}






			
		}
	}
	else {
		int start = k[k.size() - 1] + 1;
		for (int i = start; i < n; i++) {
			k.push_back(i);
			iterativeEnu(k, curLayer + 1, totalLayer, r1, r2, n, d, validNum);
			k.pop_back();
		}
	}
}

void Degree::outputTime() {
	cout << "min_time:" << min_time << "seconds;       ";
	cout << "max_time:" << max_time << "seconds" << endl;
}

void Degree::searchParameters(int r1, int r2, int n, int d, int w) {
	vector<int> k;
	k.clear();
	int validNum = 0;
	k.push_back(0);
	iterativeEnu(k, 1, w, r1, r2, n, d, validNum);
	cout << "total valid candidates:" << validNum << endl;
}


