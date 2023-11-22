#include<iostream>
#include<vector>
#include"degree.h"
using namespace std;

void result1(Degree& degree) {
	cout << "find optimal (h1,h2,h3,h4) for (n,d)=(129,1):" << endl;
	degree.searchParameters(7, 8, 129, 1, 4);
	degree.outputTime();
	degree.min_time = 100000;
	degree.max_time = 0;

}

void result2(Degree& degree) {
	cout << "find optimal (h1,h2,h3) for (n,d)=(63,32):" << endl;
	degree.searchParameters(5, 6, 63, 32, 3);
	degree.outputTime();
	degree.min_time = 100000;
	degree.max_time = 0;

}



void result3(Degree& degree) {
	cout << "test upper bounds for different (h1,h2) where (n,d)=(129,1):" << endl;
	int arr[5][2] = { {0,6},{0,9}, {0,36}, {0,42}, {0,63} };
	vector<int> h(2);
	for (int i = 0; i < 5; i++) {
		cout << "h1, h2:";
		for (int j = 0; j < 2; j++) {
			h[j] = arr[i][j];
		}
		cout << h[0] << " " << h[1] << endl;
		int n = 129, d = 1;
		for (int r = 2; r < 30; r++) {
			int res = degree.upperBoundUnivariate(r, n, d, h);
			cout << "-->round " << r << ":" << res << endl;
			if (res == n)
				break;
		}
		cout << endl;
	}

}

void result4(Degree& degree) {
	cout << "test upper bounds for different (h1,h2) where (n,d)=(63,32):" << endl;
	int arr[5][2] = { {0,3},{0,6},{0,36} };
	vector<int> h(2);
	for (int i = 0; i < 3; i++) {
		cout << "h1, h2:";
		for (int j = 0; j < 2; j++) {
			h[j] = arr[i][j];
		}
		cout << h[0] << " " << h[1] << endl;
		int n = 63, d = 32;
		for (int r = 2; r < 30; r++) {
			int res = degree.upperBoundUnivariate(r, n, d, h);
			cout << "--round " << r << ":" << res << endl;
			//degree.min_time = 100000;
			//degree.max_time = 0;
			if (res == n)
				break;
		}
		cout << endl;
	}
}


int main() {
	Degree degree;
	int cmd = 1;
	cout << "please input your command:" << endl;
	cout << "0----->exit" << endl;
	cout << "1----->find optimal (h1,h2,h3,h4) for (n,d)=(129,1)" << endl;
	cout << "2----->find optimal (h1,h2,h3) for (n,d)=(63,32)" << endl;
	cout << "3----->test upper bounds for different (h1,h2) where (n,d)=(129,1)" << endl;
	cout << "4----->test upper bounds for different (h1,h2) where (n,d)=(63,32)" << endl;
	cout << "5----->test upper bound for your choice of(h1, ..., h_w) and (n, d, r)" << endl;
	cout << "6----->compute the maximal number pi for your choice of(h1, ..., h_w) and (n, d, r)" << endl;
	cout << "7----->test upper bounds for your choice of(h1, ..., h_w) and (n, d)" << endl;


	while (cin >> cmd) {
		if (cmd == 1)
			result1(degree);

		else if (cmd == 2)
			result2(degree);

		else if (cmd == 3)
			result3(degree);

		else if (cmd == 4)
			result4(degree);

		else if (cmd == 5) {
			
			int n = 0, d = 0, r = 0, w = 0;
			cout << "input n,d,r,w (seperated with space):";
			cin >> n >> d >> r >> w;
			vector<int> h(w);
			cout << "please input h1, ..., h" << w << " (seperated with space):";


			for (int i = 0; i < w; i++) {
				cin >> h[i];
			}

			int res = degree.upperBoundUnivariate(r, n, d, h);
			cout << "-->round " << r << ":" << res << endl;
			cout << endl;
			
		}

		else if (cmd == 6) {

			int n = 0, d = 0, r = 0, w = 0;
			cout << "input n,d,r,w (seperated with space):";
			cin >> n >> d >> r >> w;
			vector<int> h(w);
			cout << "please input h1, ..., h" << w << " (seperated with space):";


			for (int i = 0; i < w; i++) {
				cin >> h[i];
			}

			int pi=degree.MaxNum(h, r, n, d);
			cout << "the maximal number:" << pi << "   ";
			cout << "-solving time:" << degree.max_time << "s"<<endl;
			cout << endl;

			degree.min_time = 100000;
			degree.max_time = 0;

		}

		else if (cmd == 7) {

			int n = 0, d = 0, r = 0, w = 0;
			cout << "input n,d,w (seperated with space):";
			cin >> n >> d >> w;
			vector<int> h(w);
			cout << "please input h1, ..., h" << w << " (seperated with space):";


			for (int i = 0; i < w; i++) {
				cin >> h[i];
			}

			for (int r = 2; r < 30; r++) {
				int res = degree.upperBoundUnivariate(r, n, d, h);
				cout << "-->round " << r << ":" << res << endl;
				if (res == n)
					break;
			}
			cout << endl;

		}



		
		else {
			break;
		}


		cout << "please input your command:" << endl;
		cout << "0----->exit" << endl;
		cout << "1----->find optimal (h1,h2,h3,h4) for (n,d)=(129,1)" << endl;
		cout << "2----->find optimal (h1,h2,h3) for (n,d)=(63,32)" << endl;
		cout << "3----->test upper bounds for different (h1,h2) where (n,d)=(129,1)" << endl;
		cout << "4----->test upper bounds for different (h1,h2) where (n,d)=(63,32)" << endl;
		cout << "5----->test upper bound for your choice of(h1, ..., h_w) and (n, d, r)" << endl;
		cout << "6----->compute the maximal number pi for your choice of(h1, ..., h_w) and (n, d, r)" << endl;
		cout << "7----->test upper bounds for your choice of(h1, ..., h_w) and (n, d)" << endl;

	}

}