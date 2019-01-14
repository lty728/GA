#include <iostream>
#include <random>
#include <vector>
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>

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
	vector<int> infecttimes;//����Ⱦ����
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
				n1 = u(e);
				n2 = u(e);
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

	void initSpread(int x0) {
		susceptible.clear();
		susceptible.resize(nsize, false);
		suslist.clear();
		infected.clear();
		infected.resize(nsize, false);
		ninfected.clear();
		ninfected.resize(nsize, 0);
		infecttimes.clear();
		infecttimes.resize(nsize, 0);
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
	}


	void set_beta(double sbeta) { //����ÿ���ڵ�ĸ�Ⱦ����
		beta.clear();
		beta.resize(nsize, sbeta);
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
		set_beta(0.65);
		for (int i = 0; i < timelist.size(); i++) {
			timelist[i] = i*T;
		}
	}

	void showtime(clock_t start, clock_t end) {
		cout << end - start << "ms" << endl;
	}

	void SI(int repeat, int seeds, int it,double T, double endt, string SI_ratio, string SI_infected_neighbor) {
		totstart = clock();
		init(T, endt);
		ofstream out1(SI_ratio); //��Ⱦ�����
		int pos;
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << ":" << " ";
			ofstream out2(SI_infected_neighbor + "_" + to_string(rp) + ".txt");//�Ѹ�Ⱦ�ھ��������
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
				if (suslist.size() > 0) {
					addr = 0;
					uniform_real_distribution<double> u(0, si_totr);
					rt = u(e);
					dt = -log(rt / si_totr) / si_totr; //����ʱ��
					t += dt;
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						addr += ninfected[node] * beta[node];
						if (addr > rt) {
							infecttimes[node]++;
							pos = i;//��¼λ��
							break;
						}
					}
					if (infecttimes[node] >= it) {
						suslist.erase(suslist.begin() + pos);
						out2 << node << " " << ninfected[node] << endl;
						infectednum++;
						susceptible[node] = false;
						infected[node] = true;
						si_totr -= ninfected[node] * beta[node]; //�Ƴ��ýڵ�����
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
					ratiolist[ti] += (double)infectednum / (nsize - begin);
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

	void TH(int repeat, int seeds, int it, double th, double T, double endt, string TH_ratio, string TH_infected_neighbor) {
		totstart = clock();
		init(T, endt);
		set_threshold(th);
		ofstream out1(TH_ratio); //��Ⱦ�����
		int pos;
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << ":" << " ";
			ofstream out2(TH_infected_neighbor + "_" + to_string(rp)  + ".txt"); //�Ѹ�Ⱦ�ھ��������
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
				if (suslist.size() > 0) {
					addr = 0;
					uniform_real_distribution<double> u(0, th_totr);
					rt = u(e);
					dt = -log(rt / th_totr) / th_totr; //����ʱ��
					t += dt;;
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						addr += beta[node];
						if (addr > rt) {
							pos = i;
							infecttimes[node]++;
							break;
						}
					}
					if (infecttimes[node] >= it) {
						suslist.erase(suslist.begin() + pos);
						out2 << node << " " << ninfected[node] << endl;
						infectednum++;
						susceptible[node] = false;
						infected[node] = true;
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
					ratiolist[ti] += (double)infectednum / (nsize - begin);
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

	void hybrid(int repeat, int seeds, double sibeta, double thbeta, double th, double T, double endt, string hy_ratio, string hy_infected_neighbor) {
		totstart = clock();
		init(T, endt);
		set_si_beta(sibeta);
		set_th_beta(thbeta);
		int th_num; //��ֵ��Ⱦ�ڵ���
		ofstream out1(hy_ratio); //��Ⱦ�����
		set_threshold(th);
		uniform_real_distribution<double> u(0, 1);
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << " ";
			ofstream out2(hy_infected_neighbor + "_" + to_string(rp) + ".txt"); //�Ѹ�Ⱦ�ھ��������
			initSpread(seeds);
			th_num = 0;
			for (int i = 0; i < nsize; i++) {//�������ӽڵ㣬�����Ⱦ�ڵ㣬��������
				if (infected[i]) { //
					for (int k = 0; k < n[i].size(); k++) {
						neighbor = n[i][k];
						if (infected[neighbor]) {
							ninfected[i]++;
						}
						else {
							ninfected[neighbor]++;
							si_totr += si_beta[neighbor];
							if (!susceptible[neighbor]) { //�ھ�δ������Ⱦ
								suslist.push_back(neighbor);
								susceptible[neighbor] = true;
							}
						}
					}
					out2 << i << " " << ninfected[i] << endl;
				}
			}
			for (int i = 0; i < suslist.size(); i++) {
				node = suslist[i];
				if (ninfected[node] >= threshold[node]) { //�ھӴﵽ��ֵ
					th_num++;
					th_totr += th_beta[node];
					si_totr -= ninfected[node] * si_beta[node];
				}
			}
			ratiolist[ti] += (double)infectednum / (nsize - begin);
			bt += T;
			while (suslist.size()>0)
			{
				addr = 0;
				if (th_num > 0) { //�дﵽ��ֵ�Ľڵ�
					if (th_totr <= 0) {
						cout << endl;
					}				
					rt = u(e);
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						if (ninfected[node] >= threshold[node]) {
							addr += th_beta[node];
							if (addr >= rt*th_totr) {
								suslist.erase(suslist.begin() + i);
								break;
							}
						}
					}
					th_num--;
					dt = -log(rt) / (th_totr + si_totr); //����ʱ��
					th_totr -= th_beta[node]; //�Ƴ��ýڵ����� 
					out2 << node << " " << ninfected[node] << endl;
				}
				else { //û�дﵽ��ֵ�Ľڵ�
					rt = u(e);
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						addr += si_beta[node] * ninfected[node];
						if (addr >= rt*si_totr) {
							suslist.erase(suslist.begin() + i);
							break;
						}
					}
					dt = -log(rt) / si_totr; //����ʱ��
					si_totr -= si_beta[node] * ninfected[node]; //�Ƴ��ýڵ�����
					out2 << node << " " << ninfected[node] << endl;
				}
				for (int i = 0; i < n[node].size(); i++) { //���±���Ⱦ�ڵ���ھӵ���Ϣ
					neighbor = n[node][i];
					if (!infected[neighbor]) { //���¸ýڵ�δ����Ⱦ�ھӵĸ�Ⱦ�ھ���
						ninfected[neighbor]++; //�ھӵĸ�Ⱦ�ھ�����1
						if (!susceptible[neighbor]) {//��������Ⱦ�б���
							suslist.push_back(neighbor);
							susceptible[neighbor] = true;
							if (ninfected[neighbor] >= threshold[neighbor]) { //��Ⱦ�ھ�����������ֵ���������ֵ����
								th_totr += th_beta[neighbor];
								th_num++;
							}
							else {
								si_totr += si_beta[neighbor];
							}
						}
						else { //����Ⱦ�б���
							if (((ninfected[neighbor] - 1) < threshold[neighbor]) && (ninfected[neighbor] >= threshold[neighbor])) { //ԭ��δ����ֵ���ڱ��θ�Ⱦ��ﵽ��ֵ
								th_totr += th_beta[neighbor];
								th_num++;
								si_totr -= (ninfected[neighbor] - 1)*si_beta[neighbor];
							}
							else if (ninfected[neighbor]<threshold[neighbor]) {//���θ�Ⱦ����δ����ֵ
								si_totr += si_beta[neighbor];
							}
							else {}//֮ǰ���Ѵ���ֵ
						}
					}
				}
				infectednum++;
				susceptible[node] = false;
				infected[node] = true;
				t += dt;
				if (t > endt) {
					cout << endl;
					break;
				}
				while (bt <= t) {
					ti++;
					ratiolist[ti] += (double)infectednum / (nsize - begin);
					bt += T;
				}
			}
			out2.close();
			while (ti<ratiolist.size() - 1) {
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
	//n.createER("F:/data/ER_1000_6.txt"); //����һ��ER����
	n.readdata("F:/data/ER_1000_6.txt"); //��ȡ��������
	//n.SI(100, 1, 3, 0.1, 20, "F:/data/result/SI/SI_ratio.txt", "F:/data/result/SI/SI_infected_neighbor"); //(�ظ�������������������Ⱦ���������ʱ��������ʱ�䣬��Ⱦ���ļ�����Ⱦ�ھ������ļ�)
	n.TH(100, 1, 3, 0.1, 0.1, 40, "F:/data/result/TH/TH_ratio.txt", "F:/data/result/TH/TH_infected_neighbor"); //(�ظ�������������������ֵ�����ʱ��������ʱ�䣬��Ⱦ���ļ�����Ⱦ�ھ������ļ�)
	//n.hybrid(10, 1, 1, 100, 0.5, 0.05, 5, "F:/data/result/hy2/hy2_ratio.txt", "F:/data/result/hy2/hy2_infected_neighbor");
	//int repeat, int seeds, double sibeta, double thbeta, double th, double T, double endt
	//system("pause");
}

