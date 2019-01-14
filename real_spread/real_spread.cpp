#include <iostream>
#include <random>
#include <vector>
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <io.h>

using namespace std;

class Network
{
private:
	int nsize;
	int realsize;
	vector<vector<int> > n; //�ڽӱ�
	vector<double> beta;
	vector<double> si_beta; //si���ʱ�
	vector<double> th_beta; //th���ʱ�
	vector<bool> infected; //�ڵ��Ƿ񱻸�Ⱦ
	vector<int> ninfected; //�ڵ��Ⱦ�ھ���
	int infectednum; //��Ⱦ����
	vector<double> ratiolist;
	double ratio; //��Ⱦ��
	vector<double> timelist; //ʱ���
	vector<double> threshold; //��ֵ��
	vector<bool> susceptible; //�ڵ��Ƿ�ɱ���Ⱦ������Ⱦ���У�
	vector<int> suslist; //�ɱ���Ⱦ�ڵ��
	vector<string> files;
	vector<vector<int> > kk;
	vector<int> nk; //�ڵ�Ķ�
	default_random_engine e;
	double si_totr; //SI������
	double th_totr; //��ֵ������
	int node; //��ǰ����Ⱦ�ڵ�
	double addr; //�ۼ�����
	double rt; //�����(0,totr)
	int neighbor; //�ھӽڵ�
	double t; //��ʱ��
	double dt; //��Ⱦһ���ڵ����ӵ�ʱ��
	double bt; //�´����ʱ��
	int ti; //��ti��ʱ����
	clock_t totstart, totend; //��ʱ��
	clock_t start, end; //����ʵ���ʱ��
	int begin;//��ſ�ʼ�ڵ�
	vector<int> spread;//��Ⱦ����

public:
	void createER(string filename) {  //����er����
		cout << "N:";
		cin >> nsize;
		double k;//ƽ����
		cout << "<k>:";
		cin >> k;
		int M = nsize*k / 2;
		e.seed(time(0));
		int n1, n2;
		n.clear();
		n.resize(nsize + 1);
		uniform_int_distribution<int> u(1, nsize);
		ofstream out(filename);
		if (out.is_open()) {
			int edge = 0;
			while (edge<M)
			{
				n1 = (u(e));
				n2 = (u(e));
				if ((n1 != n2) && (disconnected(n1, n2))) {
					n[n1].push_back(n2);
					n[n2].push_back(n1);
					out << n1 << " " << n2 << endl;
					edge++;
				}
			}
		}
		out.close();
	}

	void getnk() {
		nk.clear();
		nk.resize(n.size());
		for (int i = 0; i < nk.size(); i++) {
			nk[i] = n[i].size();
		}
	}

	void getAllFile(const char* path) {
		_finddata_t file;
		intptr_t lf;
		//�����ļ���·��
		if ((lf = _findfirst(path, &file)) == -1)  //����*.*��Ϊ*.txt�����������ļ��������е�txt�ļ���
			cout << "Not Found!" << endl;
		else {
			//����ļ���
			files.push_back(file.name);
			while (_findnext(lf, &file) == 0) {
				files.push_back(file.name);
			}
		}
		_findclose(lf);
	}

	void countdegree(const char* path, string outpath) {
		getnk();
		vector<int> degree = { 10,20,30,40,50,60,70,80 };
		kk.clear();
		kk.resize(degree.size());
		char path2[30];
		strcpy(path2, path);
		strcat(path2, "*.*");
		kk.clear();
		getAllFile(path2);
		string fpath = path;
		string filename;
		string line;
		int sum;
		int m;
		for (int i = 0; i < files.size(); i++) {//��ȡ�����ļ�
			filename = fpath + files[i];
			ifstream in(filename);
			while (getline(in, line)) 
			{
				istringstream is(line);
				is >> node >> m;
				switch (nk[node])
				{
				case 10:
					kk[0].push_back(m);
					break;
				case 20:
					kk[1].push_back(m);
					break;
				case 30:
					kk[2].push_back(m);
					break;
				case 40:
					kk[3].push_back(m);
					break;
				case 50:
					kk[4].push_back(m);
					break;
				case 60:
					kk[5].push_back(m);
					break;
				case 70:
					kk[6].push_back(m);
					break;
				case 80:
					kk[7].push_back(m);
					break;
				default:
					break;
				}
			}
		}
		for (int i = 0; i < kk.size(); i++) {
			ofstream out(outpath + "_" + to_string(degree[i]) + ".txt");
			sort(kk[i].begin(), kk[i].end());
			m = kk[i][0];
			sum = 1;
			for (int j = 1; j < kk[i].size(); j++) {
				if (m == kk[i][j]) {
					sum++;
				}
				else {
					out << (double)m / degree[i] << " " << ((double)sum / kk[i].size())*degree[i] << endl;
					m = kk[i][j];
					sum = 1;
				}
			}
			out << (double)m / degree[i] << " " << ((double)sum / kk[i].size())*degree[i] << endl;
			out.close();
		}

	}

	bool disconnected(int n1, int n2) {
		for (int i = 0; i < n[n1].size(); i++) {
			if (n[n1][i] == n2) {
				return false;
			}
		}
		return true;
	}

	void readdata(string filename) {
		ifstream in(filename);
		string line;
		if (in) // �и��ļ�
		{
			int n1, n2;
			n.clear();
			nsize = 0;
			while (getline(in, line)) // line�в�����ÿ�еĻ��з�
			{
				istringstream is(line);
				is >> n1 >> n2;
				if ((n1 >= nsize) || (n2 >= nsize)) {
					nsize = ((n1 > n2) ? n1 : n2) + 1;
					n.resize(nsize);
				}
				n[n1].push_back(n2);
				n[n2].push_back(n1);
			}
			realsize = 0;
			for (int i = 0; i < nsize; i++) {
				if (n[i].size() != 0) {
					realsize++;
				}
				else
				{
					cout << i << endl;
				}
			}
			e.seed(time(0));
			srand((unsigned)time(NULL));
			if (n[0].size() != 0) {
				begin = 0;
			}
			else
			{
				begin = 1;
			}
		}
		else // û�и��ļ�
		{
			cout << "no such file" << endl;
			return;
		}
		cout << "data loaded" << endl;
	}
	
	void readSpread(string filename) {
		ifstream in(filename);
		string line;
		spread.clear();
		if (in) // �и��ļ�
		{
			while (getline(in, line)) // line�в�����ÿ�еĻ��з�
			{
				istringstream is(line);
				is >> node;
				spread.push_back(node);
			}
		}
		else // û�и��ļ�
		{
			cout << "no such file" << endl;
			return;
		}
		cout << "spread loaded" << endl;
	}

	void initSpread(int x0) {
		susceptible.clear();
		susceptible.resize(nsize, false);
		suslist.clear();
		infected.clear();
		infected.resize(nsize, false);
		ninfected.clear();
		ninfected.resize(nsize, 0);
		infectednum = 0;
		si_totr = 0;
		th_totr = 0;
		set_seeds_random(x0);
		t = 0;
		bt = 0;
		dt = 0;
		ti = 0;
	}

	void set_seeds_random(int x0) {
		int seedi;
		uniform_int_distribution<int> u(0, nsize - 1);
		for (int i = 0; i < x0; ) {
			seedi = u(e);
			if ((!infected[seedi]) && (n[seedi].size()>0)) { //δ����Ⱦ
				infected[seedi] = true;
				infectednum++;
				i++;
			}
		}  
//		infected[10] = true;
//		infectednum++; 
	}


	void set_beta() { //����ÿ���ڵ�ĸ�Ⱦ����
		uniform_real_distribution<double> v(0, 1);
		for (int i = 0; i < nsize; i++) {
			if (v(e) <= 0.5) {
				beta[i] = 0.65;
			}
			else
			{
				beta[i] = 0.65;
			}
		}
             /*for (int i = 0; i < nsize; i++) {
			beta[i] = 2;
		}*/
	}

	void set_si_beta(double sibeta) { //����si�ڵ�ĸ�Ⱦ����
		si_beta.clear();
		si_beta.resize(nsize, sibeta);
	}

	void set_th_beta(double thbeta) { //����th�ڵ�ĸ�Ⱦ����
		th_beta.clear();
		th_beta.resize(nsize, thbeta);
	} 

	void set_threshold(double th) { //����ÿ���ڵ����ֵ
		threshold.clear();
		threshold.resize(nsize, 0);
		for (int i = 0; i < nsize; i++) {
			threshold[i] = n[i].size()*th;
		}
	}

	void init(double T, double endt) {
		ratiolist.clear();
		ratiolist.resize(endt / T + 1, 0);
		beta.clear();
		beta.resize(nsize, 0);
		timelist.clear();
		timelist.resize(endt / T + 1, 0);
		set_beta();
		for (int i = 0; i < timelist.size(); i++) {
			timelist[i] = i*T;
		}
	}

	void showtime(clock_t start, clock_t end) {
		cout << (end - start)/CLOCKS_PER_SEC << "s" << endl;
	}
	
	void real_spreading(int rp, string spreadfile, string real_neighbor){
		readSpread(spreadfile);
    	         ninfected.clear();
		ninfected.resize(nsize, 0);
    	         ofstream out3(real_neighbor + "_" + to_string(rp)); //output result
    	         for (int i = 0; i < spread.size(); i++) {
	             node =  spread[i];
	             out3 << node << " " << ninfected[node] << endl;
        	             for (int j = 0; j < n[node].size(); j++) {
    			neighbor = n[node][j];
    			ninfected[neighbor]++;
			}
		}	
	}


	void SI(int repeat, int seeds, double T, double endt, string SI_ratio, string SI_infected_neighbor) {
		totstart = clock();
		init(T, endt);
		ofstream out1(SI_ratio); //��Ⱦ�����
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << ":" << endl;
			ofstream out2(SI_infected_neighbor + "_" + to_string(rp));//�Ѹ�Ⱦ�ھ��������
			initSpread(seeds);
			for (int i = 0; i < nsize; i++) {//�������ӽڵ㣬�����Ⱦ�ڵ㣬��������
				if (infected[i]) { //
					for (int k = 0; k < n[i].size(); k++) {
						if (infected[n[i][k]]) {
							ninfected[i]++;
						}
					}
					out2 << i << " " << ninfected[i] << endl;
					for (int j = 0; j < n[i].size(); j++) {
						neighbor = n[i][j];
						if (!infected[neighbor]) {
							ninfected[neighbor]++;
							si_totr += beta[neighbor];
							if (!susceptible[neighbor]) {
								suslist.push_back(neighbor);
								susceptible[neighbor] = true;
							}
						}
					}
				}
			}
			ratiolist[ti] += (double)infectednum / (nsize - begin);
			bt += T;
			while (true)
			{
				addr = 0;
				uniform_real_distribution<double> u(0, si_totr);
				rt = u(e);
				if (suslist.size() > 0) {
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						addr += ninfected[node] * beta[node];
						if (addr > rt) {
							suslist.erase(suslist.begin() + i);
							break;
						}
					}
					//cout << node << " " << ninfected[node] << endl;
					out2 << node << " " << ninfected[node] << endl;
					infectednum++;
					susceptible[node] = false;
					infected[node] = true;
					dt = -log(rt / si_totr) / si_totr; //����ʱ��
					si_totr -= ninfected[node] * beta[node]; //�Ƴ��ýڵ�����
					t += dt;
					ninfected[node] = 0;
					for (int i = 0; i < n[node].size(); i++) {//����ñ���Ⱦ�ڵ���ھ�
						neighbor = n[node][i];
						if (!infected[neighbor]) { //��Ӹýڵ�δ����Ⱦ�ھ�Ϊ��Ⱦ�ڵ㣬��������
							ninfected[neighbor]++;
							si_totr += beta[neighbor];
							if (!susceptible[neighbor]) {
								susceptible[neighbor] = true;
								suslist.push_back(neighbor);
							}
						}
					}
				}
				else
				{
					break;
				}
				if (t > endt) {
					break;
				}
				while (bt <= t) {
					ti++;
					ratiolist[ti] += (double)(infectednum - 1) / (nsize - begin);
					bt += T;
				}

			}
			out2.close();
			while (ti < ratiolist.size()-1) {
				ti++;
				ratiolist[ti] += (double)infectednum / (nsize - begin);
			}
			end = clock();
			showtime(start, end);
		}
		for (int i = 0; i < ratiolist.size(); i++) {
			out1 << timelist[i] << " " << ratiolist[i] / repeat << endl;
		}
		out1.close();
		totend = clock();
		cout << "total time: " ;
		showtime(totstart, totend);
	}

	void TH(int repeat, int seeds, double th, double T, double endt, string TH_ratio, string TH_infected_neighbor) {
		totstart = clock();
		init(T, endt);
		set_threshold(th);
		ofstream out1(TH_ratio); //��Ⱦ�����
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << " ";
			ofstream out2(TH_infected_neighbor + "_" + to_string(rp)); //�Ѹ�Ⱦ�ھ��������
			initSpread(seeds);
			for (int i = 0; i < nsize; i++) {//�������ӽڵ�i�������Ⱦ�ڵ㣬��������
				if (infected[i]) { //
					for (int k = 0; k < n[i].size(); k++) {
						neighbor = n[i][k];
						if (infected[neighbor]) {//�ھ��Ѹ�Ⱦ
							ninfected[i]++;
						}
						else{//�ھ�δ��Ⱦ
							ninfected[neighbor]++;
						}
					}
					out2 << i << " " << ninfected[i] << endl;
				}
			}
			for (int i = 0; i < nsize; i++) {
				if ((!infected[i])&&(n[i].size()>0)&&(ninfected[i] >= threshold[i])) {
					suslist.push_back(i);
					susceptible[i] = true;
					th_totr += beta[i];
				}
			}
			ratiolist[ti] += (double)infectednum / (nsize - begin);
			bt += T;
			while (true)
			{
				addr = 0;
				uniform_real_distribution<double> u(0, th_totr);
				rt = u(e);
				if (suslist.size() > 0) {
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						addr += beta[node];
						if (addr > rt) {
							suslist.erase(suslist.begin() + i);
							break;
						}
					}
					if (ninfected[node] > n[node].size()) {
						cout << node << "-" << ninfected[node] << ">" << n[node].size() << endl;
					}
					out2 << node << " " << ninfected[node] << endl;
					infectednum++;
					susceptible[node] = false;
					infected[node] = true;
					dt = -log(rt / th_totr) / th_totr; //����ʱ��
					t += dt;;
					th_totr -= beta[node]; //�Ƴ��ýڵ�����
					for (int i = 0; i < n[node].size(); i++) {
						neighbor = n[node][i];
						if (!infected[neighbor]) { //���¸ýڵ�δ����Ⱦ�ھӵĸ�Ⱦ�ھ���
							ninfected[neighbor]++; //�ھӵĸ�Ⱦ�ھ�����1
							if (!susceptible[neighbor]) {//��������Ⱦ�б���
								if (ninfected[neighbor] >= threshold[neighbor]) { //��Ⱦ�ھ�����������ֵ���������Ⱦ�б���������
									th_totr += beta[neighbor];
									susceptible[neighbor] = true;
									suslist.push_back(neighbor);
								}
							}
						}
					}
				}
				else
				{
					break;
				}

				if (t > endt) {
					break;
				}
				while (bt <= t) {
					ti++;
					ratiolist[ti] += (double)(infectednum - 1) / (nsize - begin);
					bt += T;
				}
			}
			out2.close();
			while (ti<ratiolist.size()-1 ) {
				ti++;
				ratiolist[ti] += (double)infectednum / (nsize - begin);
			}
			end = clock();
			showtime(start, end);
		}
		for (int i = 0; i < ratiolist.size(); i++) {
			out1 << timelist[i] << " " << ratiolist[i] / repeat << endl;
		}
		out1.close();
		totend = clock();
		cout << "total time: ";
		showtime(totstart, totend);
	}
};

int main()
{
	Network n;
//	n.createER("ER.txt"); // create ER network
	n.readdata("user_relation.txt"); // read network
	/*for (int i = 0; i < 26994; i++) {
		string datafile = "weiboid"+ to_string(i)+".txt";
		n.real_spreading(i, "./real_spread/" + datafile, "./SI/SI_infected_neighbor");
	}*/
	n.countdegree("./real_spread/", "./real_spread/");
//	n.SI(200, 50, 0.1, 20, "./SI/SI_ratio.txt", "./SI/SI_infected_neighbor"); //(�ظ��������������������ʱ��������ʱ�䣬��Ⱦ���ļ�����Ⱦ�ھ������ļ�)
//	n.TH(200, 1, 0.001, 0.1, 20, "./TH/TH_ratio.txt", "./TH/TH_infected_neighbor"); //(�ظ�������������������ֵ�����ʱ��������ʱ�䣬��Ⱦ���ļ�����Ⱦ�ھ������ļ�)
	system("pause");
}