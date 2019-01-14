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
	vector<vector<int> > n; //邻接表
	vector<double> beta;
	vector<double> si_beta; //si速率表
	vector<double> th_beta; //th速率表
	vector<bool> infected; //节点是否被感染
	vector<int> ninfected; //节点感染邻居数
	int infectednum; //感染数量
	vector<double> ratiolist;
	double ratio; //感染率
	vector<double> timelist; //时间表
	vector<double> threshold; //阈值表
	vector<bool> susceptible; //节点是否可被感染（在易染表中）
	vector<int> suslist; //可被感染节点表 
	vector<int> infecttimes;//被感染次数
	default_random_engine e;
	double si_totr; //SI总速率
	double th_totr; //阈值总速率
	int node; //当前待感染节点
	double addr; //累加速率
	double rt; //随机数(0,totr)
	int neighbor; //邻居节点
	double t; //总时间
	double dt; //感染一个节点增加的时间
	double bt; //下次输出时间
	int ti; //第ti次时间间隔
	clock_t totstart, totend; //总时间
	clock_t start, end; //单次实验的时间
	int begin;//编号开始节点

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
			if ((!infected[seedi]) && (n[seedi].size()>0)) { //未被感染
				infected[seedi] = true;
				infectednum++;
				i++;
			}
		}
	}


	void set_beta(double sbeta) { //设置每个节点的感染速率
		beta.clear();
		beta.resize(nsize, sbeta);
	}

	void set_si_beta(double sibeta) { //设置si节点的感染速率
		si_beta.clear();
		si_beta.resize(nsize, sibeta);
	}

	void set_th_beta(double thbeta) { //设置th节点的感染速率
		th_beta.clear();
		th_beta.resize(nsize, thbeta);
	}

	void set_threshold(double th) { //设置每个节点的阈值
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
		ofstream out1(SI_ratio); //感染率输出
		int pos;
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << ":" << " ";
			ofstream out2(SI_infected_neighbor + "_" + to_string(rp) + ".txt");//已感染邻居数量输出
			initSpread(seeds);
			for (int i = 0; i < nsize; i++) {//处理种子节点，添加易染节点，更新速率
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
					dt = -log(rt / si_totr) / si_totr; //计算时间
					t += dt;
					for (int i = 0; i < suslist.size(); i++) {//选出一个节点感染
						node = suslist[i];
						addr += ninfected[node] * beta[node];
						if (addr > rt) {
							infecttimes[node]++;
							pos = i;//记录位置
							break;
						}
					}
					if (infecttimes[node] >= it) {
						suslist.erase(suslist.begin() + pos);
						out2 << node << " " << ninfected[node] << endl;
						infectednum++;
						susceptible[node] = false;
						infected[node] = true;
						si_totr -= ninfected[node] * beta[node]; //移除该节点速率
						for (int i = 0; i < n[node].size(); i++) {//处理该被感染节点的邻居
							neighbor = n[node][i];
							if (!infected[neighbor]) { //添加该节点未被感染邻居为易染节点，更新速率
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
		ofstream out1(TH_ratio); //感染率输出
		int pos;
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << ":" << " ";
			ofstream out2(TH_infected_neighbor + "_" + to_string(rp)  + ".txt"); //已感染邻居数量输出
			initSpread(seeds);
			for (int i = 0; i < nsize; i++) {//处理种子节点i，添加易染节点，更新速率
				if (infected[i]) { //
					for (int k = 0; k < n[i].size(); k++) {
						neighbor = n[i][k];
						if (infected[neighbor]) {//邻居已感染
							ninfected[i]++;
						}
						else{//邻居未感染
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
					dt = -log(rt / th_totr) / th_totr; //计算时间
					t += dt;;
					for (int i = 0; i < suslist.size(); i++) {//选出一个节点感染
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
						th_totr -= beta[node]; //移除该节点速率
						for (int i = 0; i < n[node].size(); i++) {
							neighbor = n[node][i];
							if (!infected[neighbor]) { //更新该节点未被感染邻居的感染邻居数
								ninfected[neighbor]++; //邻居的感染邻居数加1
								if (!susceptible[neighbor]) {//若不在易染列表中
									if (ninfected[neighbor] >= threshold[neighbor]) { //感染邻居数量大于阈值，则加入易染列表，更新速率
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
		int th_num; //阈值易染节点数
		ofstream out1(hy_ratio); //感染率输出
		set_threshold(th);
		uniform_real_distribution<double> u(0, 1);
		for (int rp = 1; rp <= repeat; rp++) {
			start = clock();
			cout << rp << " ";
			ofstream out2(hy_infected_neighbor + "_" + to_string(rp) + ".txt"); //已感染邻居数量输出
			initSpread(seeds);
			th_num = 0;
			for (int i = 0; i < nsize; i++) {//处理种子节点，添加易染节点，更新速率
				if (infected[i]) { //
					for (int k = 0; k < n[i].size(); k++) {
						neighbor = n[i][k];
						if (infected[neighbor]) {
							ninfected[i]++;
						}
						else {
							ninfected[neighbor]++;
							si_totr += si_beta[neighbor];
							if (!susceptible[neighbor]) { //邻居未加入易染
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
				if (ninfected[node] >= threshold[node]) { //邻居达到阈值
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
				if (th_num > 0) { //有达到阈值的节点
					if (th_totr <= 0) {
						cout << endl;
					}				
					rt = u(e);
					for (int i = 0; i < suslist.size(); i++) {//选出一个节点感染
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
					dt = -log(rt) / (th_totr + si_totr); //计算时间
					th_totr -= th_beta[node]; //移除该节点速率 
					out2 << node << " " << ninfected[node] << endl;
				}
				else { //没有达到阈值的节点
					rt = u(e);
					for (int i = 0; i < suslist.size(); i++) {//选出一个节点感染
						node = suslist[i];
						addr += si_beta[node] * ninfected[node];
						if (addr >= rt*si_totr) {
							suslist.erase(suslist.begin() + i);
							break;
						}
					}
					dt = -log(rt) / si_totr; //计算时间
					si_totr -= si_beta[node] * ninfected[node]; //移除该节点速率
					out2 << node << " " << ninfected[node] << endl;
				}
				for (int i = 0; i < n[node].size(); i++) { //更新被感染节点的邻居的信息
					neighbor = n[node][i];
					if (!infected[neighbor]) { //更新该节点未被感染邻居的感染邻居数
						ninfected[neighbor]++; //邻居的感染邻居数加1
						if (!susceptible[neighbor]) {//若不在易染列表中
							suslist.push_back(neighbor);
							susceptible[neighbor] = true;
							if (ninfected[neighbor] >= threshold[neighbor]) { //感染邻居数量大于阈值，则更新阈值速率
								th_totr += th_beta[neighbor];
								th_num++;
							}
							else {
								si_totr += si_beta[neighbor];
							}
						}
						else { //在易染列表中
							if (((ninfected[neighbor] - 1) < threshold[neighbor]) && (ninfected[neighbor] >= threshold[neighbor])) { //原本未达阈值，在本次感染后达到阈值
								th_totr += th_beta[neighbor];
								th_num++;
								si_totr -= (ninfected[neighbor] - 1)*si_beta[neighbor];
							}
							else if (ninfected[neighbor]<threshold[neighbor]) {//本次感染后仍未达阈值
								si_totr += si_beta[neighbor];
							}
							else {}//之前就已达阈值
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
	//n.createER("F:/data/ER_1000_6.txt"); //生成一个ER网络
	n.readdata("F:/data/ER_1000_6.txt"); //读取网络数据
	//n.SI(100, 1, 3, 0.1, 20, "F:/data/result/SI/SI_ratio.txt", "F:/data/result/SI/SI_infected_neighbor"); //(重复次数，种子数量，感染次数，输出时间间隔，总时间，感染率文件，感染邻居数量文件)
	n.TH(100, 1, 3, 0.1, 0.1, 40, "F:/data/result/TH/TH_ratio.txt", "F:/data/result/TH/TH_infected_neighbor"); //(重复次数，种子数量，阈值，输出时间间隔，总时间，感染率文件，感染邻居数量文件)
	//n.hybrid(10, 1, 1, 100, 0.5, 0.05, 5, "F:/data/result/hy2/hy2_ratio.txt", "F:/data/result/hy2/hy2_infected_neighbor");
	//int repeat, int seeds, double sibeta, double thbeta, double th, double T, double endt
	//system("pause");
}

