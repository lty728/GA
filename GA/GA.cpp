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
	vector<vector<int> > n; //�ڽӱ�
	vector<double> beta; //���ʱ�
	vector<double> nr; //ÿ���ڵ�ĵ�ǰ����
	vector<bool> infected; //�ڵ��Ƿ񱻸�Ⱦ
	int infectednum; //��Ⱦ����
	vector<double> ratiolist;
	double ratio; //��Ⱦ��
	vector<double> timelist; //ʱ���
	vector<double> threshold; //��ֵ��
	vector<bool> susceptible; //�ڵ��Ƿ�ɱ���Ⱦ
	vector<int> suslist; //�ɱ���Ⱦ�ڵ�� 
	default_random_engine e;
	double totr; //������

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
		nr.clear();
		nr.resize(nsize, 0);
		infectednum = 0;
		totr = 0; 
		set_seeds_random(x0);
	}

	void set_seeds_random(int x0) {
		e.seed(time(0));
		int seedi;
		int neighbor;
		uniform_int_distribution<int> u(0, nsize - 1);
		for (int i = 0; i < x0; ) {
			seedi = u(e);
			if ((!infected[seedi])&&(n[seedi].size()>0)) { //δ����Ⱦ
				infected[seedi] = true;
				infectednum++;
				if (susceptible[seedi]) { //��ǰ�ɱ���Ⱦ
					susceptible[seedi] = false;
				}
				i++;
			}
		}
		ratio = (double)infectednum / nsize;
		for (int i = 0; i < nsize; i++) {
			if (susceptible[i]) {
				suslist.push_back(i);
				for (int j = 0; j < n[i].size(); j++) {
					neighbor = n[i][j];
					if (infected[neighbor]) {
						nr[i] += beta[i];
						totr += beta[i];
					}
				}
			}
		}
	}


	void set_beta() { //����ÿ���ڵ�ĸ�Ⱦ����
		e.seed(time(0));
		uniform_real_distribution<double> v(0, 1);
		for (int i = 0; i < nsize; i++) {
			if (v(e) <= 0.5) {
				beta[i] = 0.5;
			}
			else
			{
				beta[i] = 1;
			}
		}
	}

	void set_threshold(double th) { //����ÿ���ڵ����ֵ
		threshold.clear();
		threshold.resize(nsize, 0);
		for (int i = 0; i < nsize; i++) {
			threshold[i] = n[i].size()*th;
		}
	}

	void SI(int repeat,int seeds, double T, double endt,string SI_ratio,string SI_infected_neighbor) {;
		clock_t totstart, totend;
		totstart = clock();
		ratiolist.clear();
		ratiolist.resize(endt / T + 1, 0);
		beta.clear();
		beta.resize(nsize, 0);
		timelist.clear();
		timelist.resize(endt / T + 1, 0);
		set_beta();
		ofstream out1(SI_ratio); //��Ⱦ�����
		int node; //��ǰ����Ⱦ�ڵ�
		double addr; //�ۼ�����
		double rt; //�����(0,totr)
		int neighbor; //�ھӽڵ�
		int infectedneighbor; //�Ѹ�Ⱦ�ھ�����
		double t; //��ʱ��
		double dt;
		int ti; //��ti��ʱ����
		double tempr;
		uniform_real_distribution<double> v(0, 1);
		clock_t start, end;
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp <<" ";
			ofstream out2(SI_infected_neighbor + "_" + to_string(rp));//�Ѹ�Ⱦ�ھ��������
			initSpread(seeds);
			e.seed(time(0));
			t = 0;
			dt = 0;
			ti = 0;
			for (int i = 0; i < nsize; i++) {//�������ӽڵ㣬�����Ⱦ�ڵ㣬��������
				if (infected[i]) { //
					for (int j = 0; j < n[i].size(); j++) {
						neighbor = n[i][j];
						if (!infected[neighbor]) {
							nr[neighbor] += beta[neighbor];
							totr += beta[neighbor];
							if (!susceptible[neighbor]) {
								suslist.push_back(neighbor);
								susceptible[neighbor] = true;
							}
						}
					}
				}
			}
			ratiolist[ti] += (double)infectednum / nsize;
			while (true)
			{
				if (t > endt) {
					break;
				}
				addr = 0;
				uniform_real_distribution<double> u(0, totr);
				rt = u(e);
				if (suslist.size() > 0) {
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						addr += nr[node];
						if (addr > rt) {
							suslist.erase(suslist.begin() + i);
							break;
						}
					}
					infectednum++;
					susceptible[node] = false;
					infected[node] = true;
					tempr = totr;
					totr -= nr[node]; //�Ƴ��ýڵ�����
					nr[node] = 0;
					infectedneighbor = 0;
					for (int i = 0; i < n[node].size(); i++) {//����ñ���Ⱦ�ڵ���ھ�
						neighbor = n[node][i];
						if (!infected[neighbor]) { //��Ӹýڵ�δ����Ⱦ�ھ�Ϊ��Ⱦ�ڵ㣬��������
							nr[neighbor] += beta[neighbor];
							totr += beta[neighbor];
							if (!susceptible[neighbor]) {
								susceptible[neighbor] = true;
								suslist.push_back(neighbor);
							}
						}
						else {
							infectedneighbor++;
						}
					}
					out2 << node << " " << infectedneighbor << endl;
				}
				else
				{
					break;
				}
				dt += -log(v(e)) / tempr; //����ʱ��
				if (dt >= T) {
					t += dt;
					ti++;
					dt = 0;
					ratio = (double)infectednum / nsize;
					ratiolist[ti] += ratio;
					timelist[ti] += t;
				}
				
			}
			out2.close();
			end = clock();
			cout << end - start << "ms" << endl;
		}
		for (int i = 0; i < ratiolist.size(); i++) {
			out1 << timelist[i] / repeat << " " << ratiolist[i] / repeat << endl;
		}
		out1.close();
		totend = clock();
		cout << "total time:" << totend - totstart << "ms" << endl;
	}

	void TH(int repeat, int seeds, double th, double T, double endt, string TH_ratio, string TH_infected_neighbor) {
		clock_t totstart, totend;
		totstart = clock();
		vector<int> ninfected(nsize,0); //�ڵ��Ⱦ�ھ���
		ofstream out1(TH_ratio); //��Ⱦ�����
		set_threshold(th);		
		ratiolist.clear();
		ratiolist.resize(endt / T + 1, 0);
		timelist.clear();
		timelist.resize(endt / T + 1, 0);
		beta.clear();
		beta.resize(nsize, 0);
		set_beta();
		double t; //��ʱ��
		int node; //��ǰ����Ⱦ�ڵ�
		double addr; //�ۼ�����
		double rt; //�����(0,1)
		int neighbor; //����Ⱦ�ڵ��ھ�
		int infectedneighbor; //�Ѹ�Ⱦ�ھ�����
		double dt;
		int ti;
		clock_t start, end;
		double tempr;
		uniform_real_distribution<double> v(0, 1);
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << " ";
			ofstream out2(TH_infected_neighbor + "_" + to_string(rp)); //�Ѹ�Ⱦ�ھ��������
			initSpread(seeds);
			e.seed(time(0));
			t = 0;
			dt = 0;
			ti = 0;
			for (int i = 0; i < nsize; i++) {//�������ӽڵ㣬�����Ⱦ�ڵ㣬��������
				if (infected[i]) { //
					for (int j = 0; j < n[i].size(); j++) {
						neighbor = n[i][j];
						if (!infected[neighbor]) {
							nr[neighbor] += beta[neighbor];
							ninfected[neighbor]++;
							if (ninfected[neighbor] >= threshold[neighbor]) {
								suslist.push_back(neighbor);
								susceptible[neighbor] = true;
								totr += beta[neighbor];
							}
						}
					}
				}
			}
			ratiolist[ti] += (double)infectednum / nsize;
			while (true)
			{
				if (t > endt) {
					break;
				}
				addr = 0;
				uniform_real_distribution<double> u(0, totr);
				rt = u(e);
				if (suslist.size() > 0) {
					for (int i = 0; i < suslist.size(); i++) {//ѡ��һ���ڵ��Ⱦ
						node = suslist[i];
						addr += nr[node];
						if (addr > rt) {
							suslist.erase(suslist.begin() + i);
							break;
						}
					}
					infectednum++;					
					susceptible[node] = false;
					infected[node] = true;
					tempr = totr;
					totr -= nr[node]; //�Ƴ��ýڵ�����
					nr[node] = 0;
					infectedneighbor = 0;
					for (int i = 0; i < n[node].size(); i++) {
						neighbor = n[node][i];
						if (!infected[neighbor]) { //���¸ýڵ�δ����Ⱦ�ھ�����
							nr[neighbor] += beta[neighbor];
							ninfected[neighbor]++;
							if (ninfected[neighbor] >= threshold[neighbor]) { //��Ⱦ�ھ�����������ֵ���������Ⱦ�б���������
								totr += beta[neighbor];
								susceptible[neighbor] = true;
								suslist.push_back(neighbor);
							}
						}
						else {
							infectedneighbor++;
						}
					}
					out2 << node << " " << infectedneighbor << endl;
				}
				else
				{
					break;
				}
				double k = v(e);
				dt += -log(k) / tempr; //����ʱ��
				if (dt > 100) {
					cout << k << " " << tempr << endl;
				}
				if (dt >= T) {
					t += dt;
					ti++;
					dt = 0;
					ratio = (double)infectednum / nsize;
					ratiolist[ti] += ratio;
					timelist[ti] += t;
				}

			}	
			out2.close();
			end = clock();
			cout << end - start << "ms" << endl;
		}	
		for (int i = 0; i < ratiolist.size(); i++) {
			out1 << timelist[i] / repeat << " " << ratiolist[i] / repeat << endl;
		}
		out1.close();
		totend = clock();
		cout << "total time:" << totend - totstart << "ms" << endl;
	}



};

int main()
{
	Network n;
	//n.createER("F:/data/ER_100k.txt"); //����һ��ER����
	n.readdata("F:/data/ER_100k.txt"); //��ȡ��������
	//n.SI(5,1000, 0.03, 3, "F:/data/SI/SI_ratio", "F:/data/SI/SI_infected_neighbor"); //(�ظ��������������������ʱ��������ʱ�䣬��Ⱦ���ļ�����Ⱦ�ھ������ļ�)
	n.TH(5,1000, 0.05, 0.03, 3, "F:/data/TH/TH_ratio", "F:/data/TH/TH_infected_neighbor"); //(�ظ�������������������ֵ�����ʱ��������ʱ�䣬��Ⱦ���ļ�����Ⱦ�ھ������ļ�)
	system("pause");
}

