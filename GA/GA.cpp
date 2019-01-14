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
	vector<vector<int> > n; //邻接表
	vector<double> beta; //速率表
	vector<double> nr; //每个节点的当前速率
	vector<bool> infected; //节点是否被感染
	int infectednum; //感染数量
	vector<double> ratiolist;
	double ratio; //感染率
	vector<double> timelist; //时间表
	vector<double> threshold; //阈值表
	vector<bool> susceptible; //节点是否可被感染
	vector<int> suslist; //可被感染节点表 
	default_random_engine e;
	double totr; //总速率

public:
	void createER(string filename) {  //构造er网络
		cout << "N:";
		cin >> nsize;
		double k;//平均度
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
		if (in) // 有该文件  
		{
			int n1, n2;
			n.clear();
			nsize = 0;
			while (getline(in, line)) // line中不包括每行的换行符  
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
		else // 没有该文件  
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
			if ((!infected[seedi])&&(n[seedi].size()>0)) { //未被感染
				infected[seedi] = true;
				infectednum++;
				if (susceptible[seedi]) { //当前可被感染
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


	void set_beta() { //设置每个节点的感染速率
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

	void set_threshold(double th) { //设置每个节点的阈值
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
		ofstream out1(SI_ratio); //感染率输出
		int node; //当前待感染节点
		double addr; //累加速率
		double rt; //随机数(0,totr)
		int neighbor; //邻居节点
		int infectedneighbor; //已感染邻居数量
		double t; //总时间
		double dt;
		int ti; //第ti次时间间隔
		double tempr;
		uniform_real_distribution<double> v(0, 1);
		clock_t start, end;
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp <<" ";
			ofstream out2(SI_infected_neighbor + "_" + to_string(rp));//已感染邻居数量输出
			initSpread(seeds);
			e.seed(time(0));
			t = 0;
			dt = 0;
			ti = 0;
			for (int i = 0; i < nsize; i++) {//处理种子节点，添加易染节点，更新速率
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
					for (int i = 0; i < suslist.size(); i++) {//选出一个节点感染
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
					totr -= nr[node]; //移除该节点速率
					nr[node] = 0;
					infectedneighbor = 0;
					for (int i = 0; i < n[node].size(); i++) {//处理该被感染节点的邻居
						neighbor = n[node][i];
						if (!infected[neighbor]) { //添加该节点未被感染邻居为易染节点，更新速率
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
				dt += -log(v(e)) / tempr; //计算时间
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
		vector<int> ninfected(nsize,0); //节点感染邻居数
		ofstream out1(TH_ratio); //感染率输出
		set_threshold(th);		
		ratiolist.clear();
		ratiolist.resize(endt / T + 1, 0);
		timelist.clear();
		timelist.resize(endt / T + 1, 0);
		beta.clear();
		beta.resize(nsize, 0);
		set_beta();
		double t; //总时间
		int node; //当前待感染节点
		double addr; //累加速率
		double rt; //随机数(0,1)
		int neighbor; //待感染节点邻居
		int infectedneighbor; //已感染邻居数量
		double dt;
		int ti;
		clock_t start, end;
		double tempr;
		uniform_real_distribution<double> v(0, 1);
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << " ";
			ofstream out2(TH_infected_neighbor + "_" + to_string(rp)); //已感染邻居数量输出
			initSpread(seeds);
			e.seed(time(0));
			t = 0;
			dt = 0;
			ti = 0;
			for (int i = 0; i < nsize; i++) {//处理种子节点，添加易染节点，更新速率
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
					for (int i = 0; i < suslist.size(); i++) {//选出一个节点感染
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
					totr -= nr[node]; //移除该节点速率
					nr[node] = 0;
					infectedneighbor = 0;
					for (int i = 0; i < n[node].size(); i++) {
						neighbor = n[node][i];
						if (!infected[neighbor]) { //更新该节点未被感染邻居速率
							nr[neighbor] += beta[neighbor];
							ninfected[neighbor]++;
							if (ninfected[neighbor] >= threshold[neighbor]) { //感染邻居数量大于阈值，则加入易染列表，更新速率
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
				dt += -log(k) / tempr; //计算时间
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
	//n.createER("F:/data/ER_100k.txt"); //生成一个ER网络
	n.readdata("F:/data/ER_100k.txt"); //读取网络数据
	//n.SI(5,1000, 0.03, 3, "F:/data/SI/SI_ratio", "F:/data/SI/SI_infected_neighbor"); //(重复次数，种子数量，输出时间间隔，总时间，感染率文件，感染邻居数量文件)
	n.TH(5,1000, 0.05, 0.03, 3, "F:/data/TH/TH_ratio", "F:/data/TH/TH_infected_neighbor"); //(重复次数，种子数量，阈值，输出时间间隔，总时间，感染率文件，感染邻居数量文件)
	system("pause");
}

