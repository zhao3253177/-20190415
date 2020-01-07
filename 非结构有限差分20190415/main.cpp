#include<iostream>
#include<fstream>
#include<cmath>
#include"const.h"
#include"functions.h"
#include"init.h"
#include"output.h"
#include"patition.h"
#include<stdlib.h>
#include<ctime>
#include<iomanip>
#include <vector>
#include<ctime>
using namespace std;
vector <mesh> A0;
vector<vector <mesh>> A;
vector<vector <int>> ad;
vector <vector<mesh>> Ar;//记录上一时刻的物理量；

double dt;
double t_sim = 0;
int step = 0;
double xL, xR, yU, yD;
vector <polygon_mesh> M0;
vector<vector <polygon_mesh>> M;
clock_t start;
clock_t finish;
double res;
int main()
{
	srand((unsigned)time(NULL));
	using namespace MeshPara;
	using ConstPara::t_end;
	//init_mesh1();
	init_mesh();
	init_polygon_mesh();
	getType();
	reorderMesh();
	partition_Point();

	partition_Mesh();
	findNeiborSec();
	getNeibor();

	findConnectPoint();
	change_mesh();
	//remesh_bound();//将边界x或y坐标移动到和内部点相同，改善边界条件

	//out_neighbor();
	reorder_neighbor();
	coordinate_trans();
	//init_Ar();
	//init_flow_shockwaveCross();
	initFlow();
	out_M("mesh/" + methodType + "/step = " + to_string(step));
	out_res();
	start = clock();
	while (t_sim < t_end)
	{
		//judge();
		Ar = A;

		if (methodType == "F")//装配法 动网格
		{
			movemesh();
			get_dt();
			if (t_sim + dt > t_end)
				dt = t_end - t_sim;
			coordinate_trans();
			clear_Vm();
			update_IN();
			for (int j = 0; j < A[0].size(); j++)
			{
				if (A[0][j].type == "U")
					A[0][j].ρ = A[0][j - 80].ρ, A[0][j].p = A[0][j - 80].p, A[0][j].u = A[0][j - 80].u, A[0][j].v = A[0][j - 80].v;
				else if (A[0][j].type == "D")
					A[0][j].ρ = A[0][j + 80].ρ, A[0][j].p = A[0][j + 80].p, A[0][j].u = A[0][j + 80].u, A[0][j].v = A[0][j + 80].v;
				else if (A[0][j].type == "SHOCK")
				{
					if (abs(A[0][j].y - yD) < 1e-10)
						A[0][j].ρ = A[0][j + 1].ρ, A[0][j].p = A[0][j + 1].p, A[0][j].u = A[0][j + 1].u, A[0][j].v = A[0][j + 1].v;
					else if (abs(A[0][j].y - yU) < 1e-10)
						A[0][j].ρ = A[0][j - 1].ρ, A[0][j].p = A[0][j - 1].p, A[0][j].u = A[0][j - 1].u, A[0][j].v = A[0][j - 1].v;
				}
			}
			update_bound_shockwave_fitting();//更新激波边界条件，此处为装配法核心
			update_Vm();//更新网格点运动速度
			update_bound();//更新边界条件
			for (int j = 0; j < A[0].size(); j++)
			{
				if (A[0][j].type == "U")
					A[0][j].ρ = A[0][j - 80].ρ, A[0][j].p = A[0][j - 80].p, A[0][j].u = A[0][j - 80].u, A[0][j].v = A[0][j - 80].v;
				else if (A[0][j].type == "D")
					A[0][j].ρ = A[0][j + 80].ρ, A[0][j].p = A[0][j + 80].p, A[0][j].u = A[0][j + 80].u, A[0][j].v = A[0][j + 80].v;
				else if (A[0][j].type == "SHOCK")
				{
					if (abs(A[0][j].y - yD) < 1e-10)
						A[0][j].ρ = A[0][j + 1].ρ, A[0][j].p = A[0][j + 1].p, A[0][j].u = A[0][j + 1].u, A[0][j].v = A[0][j + 1].v;
					else if (abs(A[0][j].y - yU) < 1e-10)
						A[0][j].ρ = A[0][j - 1].ρ, A[0][j].p = A[0][j - 1].p, A[0][j].u = A[0][j - 1].u, A[0][j].v = A[0][j - 1].v;
				}
			}
		}
		else//捕捉 静止网格
		{
			get_dt();

			if (t_sim + dt > t_end)
				dt = t_end - t_sim;
			update_IN();
			update_bound_shockwave();//捕捉法的边界条件处理
			if (t_sim == 0)
			{
				using namespace Normal;
				int i, j;
				for (i = 0; i < A.size(); i++)
				{
					for (j = 0; j < A[i].size(); j++)
					{
						ρ2 = 1-0.2;
						if (A[i][j].x < 1)
						{
							A[i][j].ρ = ρ1;
							A[i][j].u = u1;
							A[i][j].v = v1;
							A[i][j].p = p1;
						}
						else
						{
							A[i][j].ρ = ρ2;
							A[i][j].u = u2;
							A[i][j].v = v2;
							A[i][j].p = p2;
						}
					}
				}

			}
		}


		step++;
		t_sim = t_sim + dt;
		if (step % 10 == 0)//输出临时结果
		{
			//if (abs(res - compute_res()) < 1e-10)
			//	break;
			//else
			res = compute_res();
			out_res();

			if (step % 10 == 0)
			{
				out_M("mesh/" + methodType + "/step = " + to_string(step));

				cout << "step = " << step << "   dt = " << dt << "   t = " << t_sim << "   CPU_t = " << (double)(double(clock()) - start) / CLOCKS_PER_SEC << "   res = " << res << endl;

				if (step % 10 == 0)
				{
					if (step % 10 == 0)
					{
						ofstream fout;
						fout.precision(15);
						int i, j;
						if (FlowType == "normal" || FlowType == "oblique")
						{
							fout.open("data/" + FlowType + methodType + "y=0.081step" + to_string(step) + ".dat");
							fout << "variables = x,rho,u,v,p" << endl;
							//fout << "solutiontime =" << t_sim << endl;
							for (i = 0; i < A.size(); i++)
							{
								for (j = 0; j < A[i].size(); j++)
								{
									if (abs(A[i][j].y - 0.081/*0.159089*/) < dx / 2)
										fout << A[i][j].x << "   " << A[i][j].ρ << "   " << A[i][j].u << "   " << A[i][j].v << "   " << A[i][j].p << endl;
								}
							}
						}
						else if (FlowType == "intersection")
						{

							fout.open("data/" + FlowType + methodType + "_x=0.01685 step" + to_string(step) + ".dat");
							fout << "variables = y,rho" << endl;
							for (i = 0; i < A.size(); i++)
							{
								for (j = 0; j < A[i].size(); j++)
								{
									if (abs(A[i][j].x - 0.01685) < dx / 2)
										fout << A[i][j].y << "   " << A[i][j].ρ << endl;
								}
							}
							fout.close();
							fout.clear();
							fout.open("data/" + FlowType + methodType + "_x=0.0016 step" + to_string(step) + ".dat");
							fout << "variables = y,rho" << endl;
							for (i = 0; i < A.size(); i++)
							{
								for (j = 0; j < A[i].size(); j++)
								{
									if (abs(A[i][j].x - 0.0016) < dx / 2)
										fout << A[i][j].y << "   " << A[i][j].ρ << endl;
								}
							}

						}
					}
				}

			}
		}


		//if (t_sim > 0.02)
		//	break;
		if (step > 100000)
			break;
		//if (step > 1000 && res < 1e-10)
		//	break;
	}
	out_M("mesh/" + methodType + "/step = " + to_string(step));
	system("PAUSE");
}