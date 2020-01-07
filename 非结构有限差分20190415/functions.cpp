#include"const.h"
#include"functions.h"
#include<omp.h>
#include"shockwave.h"
#include<iostream>
#include"Prandtl-Meyer.h"
using std::vector;
using namespace ConstPara;
using namespace MeshPara;

double distance(mesh& a, mesh& b)
{
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	return sqrt(dx * dx + dy * dy);
}
void get_dt()//��t
{
	extern vector<vector <mesh>> A;
	double max�� = 0, max�� = 0;
	extern double dt;
	double t;
	int i, j, k;
	double max1, max2;
	double S��, S��, c, u��, u��;

	dt = t_end;
	max1 = max2 = 0;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].type != "IN")
				continue;
			max�� = max�� = 0;
			if (A[i][j].neibor.size() > 3)
			{
				S�� = sqrt(A[i][j].��x[0] * A[i][j].��x[0] + A[i][j].��y[0] * A[i][j].��y[0]);
				S�� = sqrt(A[i][j].��x[0] * A[i][j].��x[0] + A[i][j].��y[0] * A[i][j].��y[0]);
				c = sqrt(�� * A[i][j].p / A[i][j].��);
				u�� = A[i][j].u * A[i][j].��x[0] + A[i][j].v * A[i][j].��y[0];
				u�� = A[i][j].u * A[i][j].��x[0] + A[i][j].v * A[i][j].��y[0];
				//max�� = max(S��, u��);
				//max�� = max(S��, u��);
				max�� = abs(u��) + c * S��;
				max�� = abs(u��) + c * S��;
				max1 = max(max1, max��);
				max2 = max(max2, max��);

			}
			else
			{
				for (k = 0; k < 3; k++)
				{
					S�� = sqrt(A[i][j].��x[k] * A[i][j].��x[k] + A[i][j].��y[k] * A[i][j].��y[k]);
					S�� = sqrt(A[i][j].��x[k] * A[i][j].��x[k] + A[i][j].��y[k] * A[i][j].��y[k]);
					c = sqrt(�� * A[i][j].p / A[i][j].��);
					u�� = A[i][j].u * A[i][j].��x[k] + A[i][j].v * A[i][j].��y[k];
					u�� = A[i][j].u * A[i][j].��x[k] + A[i][j].v * A[i][j].��y[k];
					//max�� = max(S��, u��);
					//max�� = max(S��, u��);
					max�� += abs(u��) + c * S��;
					max�� += abs(u��) + c * S��;
				}
				max�� = max�� / 3;
				max�� = max�� / 3;
				max1 = max(max1, max��);
				max2 = max(max2, max��);
			}
		}
	}
	t = CFL / (max1 + max2);
	//t = CFL / (max�� + max��);
	dt = min(dt, t);

}
//void update_AfromU()
//{
//	extern vector<vector <mesh>> A;
//	extern vector<vector<vector <double>>> U;
//	int i, j;
//	extern int step;
//
//	for (i = 0; i < A.size(); i++)
//	{
//#pragma omp parallel
//
//		for (j = 0; j < A[i].size(); j++)
//		{
//			A[i][j].�� = U[i][j][0];
//			A[i][j].u = U[i][j][1] / U[i][j][0];
//			A[i][j].v = U[i][j][2] / U[i][j][0];
//			A[i][j].p = (�� - 1) * (U[i][j][3] - 0.5 * A[i][j].�� * (A[i][j].u * A[i][j].u + A[i][j].v * A[i][j].v));
//			A[i][j].step = step;
//		}
//	}
//
//}
void update_bound_uniform()
{
	extern vector<vector <mesh>> A;
	int i, j;
	using namespace Init;

#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].type == "IN" || A[i][j].type == "N")
				continue;
			A[i][j].�� = ��0;
			A[i][j].u = u0;
			A[i][j].v = v0;
			A[i][j].p = p0;
			//if (A[i].type == "L")
			//{
			//	A[i].�� = ��1;
			//	A[i].u = u1;
			//	A[i].v = v1;
			//	A[i].p = p1;

			//}
			if (A[i][j].type != "L" || A[i][j].type != "R")
				A[i][j].v = 0;
		}
	}
}
void update_bound_shockwave()
{
	extern vector <mesh> A0;
	extern vector<vector <mesh>> A;
	int i, j;
	if (FlowType == "normal")
	{
		using namespace Normal;
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				//��1 = 3.8571, p1 = 10.3333, u1 = 2.6293, v1 = 0;
				//��2 = 1 + 0.2 * sin(5 * A[i][j].x), p2 = 1, u2 = 0, v2 = 0;
				//��2 = 3 + 0.5 * sin(5 * A[i][j].x);

				if (A[i][j].type == "IN")
					continue;
				else if (A[i][j].type == "U")
				{
					A[i][j].�� = A[i][j - Xnum].��;
					A[i][j].u = A[i][j - Xnum].u;
					A[i][j].v = A[i][j - Xnum].v;
					A[i][j].p = A[i][j - Xnum].p;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].�� = (A[i][A[i][j].neibor[0]].�� + A[i][A[i][j].neibor[1]].��) / 2;
						A[i][j].u = (A[i][A[i][j].neibor[0]].u + A[i][A[i][j].neibor[1]].u) / 2;
						A[i][j].v = (A[i][A[i][j].neibor[0]].v + A[i][A[i][j].neibor[1]].v) / 2;
						A[i][j].p = (A[i][A[i][j].neibor[0]].p + A[i][A[i][j].neibor[1]].p) / 2;
					}
				}
				else if (A[i][j].type == "D")
				{
					A[i][j].�� = A[i][j + Xnum].��;
					A[i][j].u = A[i][j + Xnum].u;
					A[i][j].v = A[i][j + Xnum].v;
					A[i][j].p = A[i][j + Xnum].p;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].�� = (A[i][A[i][j].neibor[0]].�� + A[i][A[i][j].neibor[1]].��) / 2;
						A[i][j].u = (A[i][A[i][j].neibor[0]].u + A[i][A[i][j].neibor[1]].u) / 2;
						A[i][j].v = (A[i][A[i][j].neibor[0]].v + A[i][A[i][j].neibor[1]].v) / 2;
						A[i][j].p = (A[i][A[i][j].neibor[0]].p + A[i][A[i][j].neibor[1]].p) / 2;
					}

				}
				else if (A[i][j].type == "L")
				{
					A[i][j].�� = ��1;
					A[i][j].u = u1;
					A[i][j].v = v1;
					A[i][j].p = p1;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].�� = (A[i][A[i][j].neibor[0]].�� + A[i][A[i][j].neibor[1]].��) / 2;
						A[i][j].u = (A[i][A[i][j].neibor[0]].u + A[i][A[i][j].neibor[1]].u) / 2;
						A[i][j].v = (A[i][A[i][j].neibor[0]].v + A[i][A[i][j].neibor[1]].v) / 2;
						A[i][j].p = (A[i][A[i][j].neibor[0]].p + A[i][A[i][j].neibor[1]].p) / 2;
					}

				}
				else if (A[i][j].type == "R")
				{
					A[i][j].�� = A[i][j - 1].��;
					A[i][j].u = A[i][j - 1].u;
					A[i][j].v = A[i][j - 1].v;
					A[i][j].p = A[i][j - 1].p;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].�� = (A[i][A[i][j].neibor[0]].�� + A[i][A[i][j].neibor[1]].��) / 2;
						A[i][j].u = (A[i][A[i][j].neibor[0]].u + A[i][A[i][j].neibor[1]].u) / 2;
						A[i][j].v = (A[i][A[i][j].neibor[0]].v + A[i][A[i][j].neibor[1]].v) / 2;
						A[i][j].p = (A[i][A[i][j].neibor[0]].p + A[i][A[i][j].neibor[1]].p) / 2;
					}

				}

			}
		}
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "U" || A[i][j].type == "D")
				{
					if (A[i][j].neibor.size() == 2)
					{
						if (abs(A[i][j].x - A[i][A[i][j].neibor[0]].x) < 1e-10 && abs(A[i][j].y - A[i][A[i][j].neibor[0]].y) < 1e-10)
						{
							A[i][j].�� = A[i][A[i][j].neibor[0]].��;
							A[i][j].u = A[i][A[i][j].neibor[0]].u;
							A[i][j].v = A[i][A[i][j].neibor[0]].v;
							A[i][j].p = A[i][A[i][j].neibor[0]].p;
						}
						else
						{
							A[i][j].�� = A[i][A[i][j].neibor[1]].��;
							A[i][j].u = A[i][A[i][j].neibor[1]].u;
							A[i][j].v = A[i][A[i][j].neibor[1]].v;
							A[i][j].p = A[i][A[i][j].neibor[1]].p;

						}
					}
				}

			}
		}




	}
	else if (FlowType == "intersection")
	{
		using namespace ShockwaveCross;
		double ��12 = -30 / 180.0 * ConstPara::pi;
		double ��13 = 30 / 180.0 * ConstPara::pi;
		mesh A12, A13, C;
		Line L12, L13;
		A12.x = -dx / 2, A12.y = A0[Pnum - 1].y - 2 * dy * 5 * Ynum / 45;
		A13.x = -dx / 2, A13.y = 2 * dy * 5 * Ynum / 45;
		L12 = getLine(��12, A12);
		L13 = getLine(��13, A13);
		C = getCrossPoint(L12, L13);
		L12 = getLine(��2, C);
		L13 = getLine(��3, C);
		mesh A1, A2, A3;
		A1.�� = ��1, A1.p = p1, A1.u = u1, A1.v = v1;
		get_down(A1, A2, ��2);
		get_down(A1, A3, ��3);

		��2 = A2.��, p2 = A2.p, u2 = A2.u, v2 = A2.v;
		��3 = A3.��, p3 = A3.p, u3 = A3.u, v3 = A3.v;
#pragma omp parallel
		for (i = 0; i < A.size(); i++)
		{


			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "IN")
					continue;
				else if (A[i][j].type == "L")
				{
					if (A[i][j].x * L12.A + A[i][j].y * L12.B + L12.C < 0)
					{
						A[i][j].�� = ��2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;
					}
					else if (A[i][j].x * L13.A + A[i][j].y * L13.B + L13.C > 0)
					{
						A[i][j].�� = ��3;
						A[i][j].u = u3;
						A[i][j].v = v3;
						A[i][j].p = p3;
					}
					else
					{
						A[i][j].�� = ��1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
				}
				else if (A[i][j].type == "R")
				{
					A[i][j].�� = A[i][j - 1].��;
					A[i][j].u = A[i][j - 1].u;
					A[i][j].v = A[i][j - 1].v;
					A[i][j].p = A[i][j - 1].p;

				}
				else if (A[i][j].type == "U")
				{
					A[i][j].�� = A[i][j - Xnum].��;
					A[i][j].u = A[i][j - Xnum].u;
					A[i][j].v = A[i][j - Xnum].v;
					A[i][j].p = A[i][j - Xnum].p;

				}
				else if (A[i][j].type == "D")
				{
					A[i][j].�� = A[i][j + Xnum].��;
					A[i][j].u = A[i][j + Xnum].u;
					A[i][j].v = A[i][j + Xnum].v;
					A[i][j].p = A[i][j + Xnum].p;
				}
			}
		}
	}
	else if (FlowType == "oblique")
	{
		using namespace Oblique;
		mesh M1;
		Line L;
		vector<mesh> t;
		M1.x = 0.0014, M1.y = 0;

		L = getLine(��, M1);

		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "IN")
					continue;
				else if (A[i][j].type == "U")
				{

					A[i][j].�� = A[i][j - Xnum].��;
					A[i][j].u = A[i][j - Xnum].u;
					A[i][j].v = A[i][j - Xnum].v;
					A[i][j].p = A[i][j - Xnum].p;

				}
				else if (A[i][j].type == "D")
				{
					if (A[i][j].x > 0.002)
					{
						A[i][j].�� = ��2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;

					}
					else
					{
						A[i][j].�� = A[i][j + Xnum].��;
						A[i][j].u = A[i][j + Xnum].u;
						A[i][j].v = A[i][j + Xnum].v;
						A[i][j].p = A[i][j + Xnum].p;

					}

				}
				else if (A[i][j].type == "L")
				{
					A[i][j].�� = ��1;
					A[i][j].u = u1;
					A[i][j].v = v1;
					A[i][j].p = p1;
				}
				else if (A[i][j].type == "R")
				{
					A[i][j].�� = A[i][j - 1].��;
					A[i][j].u = A[i][j - 1].u;
					A[i][j].v = A[i][j - 1].v;
					A[i][j].p = A[i][j - 1].p;
				}

			}
		}
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "U" || A[i][j].type == "D")
				{
					if (A[i][j].neibor.size() == 2)
					{
						if (abs(A[i][j].x - A[i][A[i][j].neibor[0]].x) < 1e-10 && abs(A[i][j].y - A[i][A[i][j].neibor[0]].y) < 1e-10)
						{
							A[i][j].�� = A[i][A[i][j].neibor[0]].��;
							A[i][j].u = A[i][A[i][j].neibor[0]].u;
							A[i][j].v = A[i][A[i][j].neibor[0]].v;
							A[i][j].p = A[i][A[i][j].neibor[0]].p;
						}
						else
						{
							A[i][j].�� = A[i][A[i][j].neibor[1]].��;
							A[i][j].u = A[i][A[i][j].neibor[1]].u;
							A[i][j].v = A[i][A[i][j].neibor[1]].v;
							A[i][j].p = A[i][A[i][j].neibor[1]].p;

						}
					}
				}

			}
		}

	}

}
void update_bound_shockwave_fitting()//�����߽磬����װ�䷨
{
	extern vector<vector <mesh>> A;
	extern vector<vector <mesh>> Ar;
	extern double xL, xR, yU, yD;
	extern int step;
	int i, j, k, m, n;
	double Mas, Ma1 = 0, Ma2;
	double Vs = 0, V1, V1n, V1t, V2n, V2t;
	double c1, c2, p2, ��2;
	double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
	double �� = 0;
	double ��1, ��2, ��3;
	double R_x, R_y;//�ο���ֵ��������Vs
	int m1, n1;
	mesh up, down;
	mesh upr, downr;
	double RR, RL, unb, ab;
	double ��20, p20, u20, v20, c20;
	double ��21, p21, u21, v21, c21;
	if (FlowType == "oblique" || FlowType == "normal")
		R_x = Ar[0][4000].x, R_y = Ar[0][4000].y;
	else if (FlowType == "intersection")
		R_x = Ar[0][0].x, R_y = Ar[0][0].y;
	for (i = 0; i < Ar.size(); i++)
	{

		for (j = 0; j < Ar[i].size(); j++)
		{
			if (A[i][j].step == step + 1 || (Ar[i][j].type != "SHOCK" && Ar[i][j].type != "DISCON"))
				continue;
			�� = get_��(R_x, R_y, Ar[i][j].x, Ar[i][j].y);//�����Ƕ�

			if ((FlowType == "oblique" || FlowType == "normal") && j == 4000)
			{
				for (int k = 0; k < A[i][j].neibor.size(); k++)
				{
					m = A[i][j].neibor[k];
					if (A[i][m].type == "SHOCK")
						x1 = A[i][m].x, y1 = A[i][m].y;
				}
				�� = get_��(R_x, R_y, x1, y1);//�����Ƕ�
				//�� = pi / 2;
			}
			�� = pi / 2;
			m = A[i][j].neiborsec;
			n = A[i][j].neiborsec_ad;
			if (A[i][j].p < A[m][n].p)
			{
				up = A[i][j];
				down = A[m][n];
				upr = Ar[i][j];
				downr = Ar[m][n];
			}
			else
			{
				down = A[i][j];
				up = A[m][n];
				downr = Ar[i][j];
				upr = Ar[m][n];
			}
			if (up.type == "SHOCK")
			{

				m1 = n1 = 0;
				for (k = 0; k < down.neibor.size(); k++)
				{
					m1 = down.neibor[k];
					if (A[i][m1].type == "SHOCK" || A[i][m1].type == "CENTER")
						continue;
					else
					{
						n1 = m1;
						break;
					}
				}

				Ma1 = get_Mu(upr, Ar[i][n1], ��);
				//std::cout << Ma1 << std::endl;
				c1 = get_c(upr.��, upr.p);
				c2 = get_c(Ar[i][n1].��, Ar[i][n1].p);
				V1n = get_un(upr, ��);
				V1t = get_ut(upr, ��);
				V2n = get_un(Ar[i][n1], ��);
				p2 = upr.p * get_p2p1(Ma1);
				��2 = upr.�� * get_��2��1(Ma1);

				down.p = p2;
				down.�� = ��2;
				RR = (V2n + 2 * c2 / (�� - 1));

				if (�� >= 0)
				{
					Vs = (upr.u * sin(��) - upr.v * cos(��)) - Ma1 * c1;
					Vs = -Vs;
					V2n = get_udn(upr, downr, Ma1, Vs, ��);
					V2t = V1t;
				}
				else if (�� < 0)
				{
					Vs = (upr.u * sin(-��) + upr.v * cos(-��)) - Ma1 * c1;
					Vs = -Vs;
					V2n = get_udn(upr, downr, Ma1, Vs, ��);
					V2t = V1t;
				}
				c21 = get_c(��2, p2);
				RL = (V2n - 2 * c21 / (�� - 1));
				unb = (RR + RL) / 2;
				ab = (�� - 1) * (RR - RL) / 4;
			}
			else if (up.type == "DISCON")
			{
				c1 = get_c(upr.��, upr.p);
				c2 = get_c(downr.��, downr.p);
				V1n = get_un(upr, ��);
				V2n = get_un(downr, ��);
				V1t = get_ut(upr, ��);
				V2t = get_ut(downr, ��);
				double J1 = V1n + 2 * c1 / (�� - 1);
				double J2 = V2n - 2 * c2 / (�� - 1);
				double S1 = upr.p / pow(upr.��, ��);
				double S2 = downr.p / pow(downr.��, ��);
				double pb, pb1, pb0;
				pb1 = upr.p;
				pb0 = downr.p;
				pb = (upr.p + downr.p) / 2;
				while (abs(pb - pb1) > 1e-10)
				{
					pb0 = pb1;
					pb1 = pb;
					double fxn_1 = 2 * sqrt(�� * pow(pb0, (�� - 1) / ��) * pow(S1, 1 / ��)) / (�� - 1) + 2 * sqrt(�� * pow(pb0, (�� - 1) / ��) * pow(S2, 1 / ��)) / (�� - 1) - (J2 - J1);
					double fxn = 2 * sqrt(�� * pow(pb1, (�� - 1) / ��) * pow(S1, 1 / ��)) / (�� - 1) + 2 * sqrt(�� * pow(pb1, (�� - 1) / ��) * pow(S2, 1 / ��)) / (�� - 1) - (J2 - J1);
					pb = pb1 - fxn / (fxn - fxn_1) * (pb1 - pb0);
				}
				//���β�������
				up.p = pb;
				up.�� = pow((up.p / S1), 1 / ��);
				c1 = get_c(up.��, up.p);
				V1n = 2 * c1 / (�� - 1) - J1;
				//���β�������
				down.p = pb;
				down.�� = pow((down.p / S2), 1 / ��);
				c2 = get_c(down.��, down.p);
				V2n = 2 * c2 / (�� - 1) + J2;
				Vs = (V1n + V2n) / 2;
			}
			if (�� >= 0)
			{
				up.um = down.um = Vs * sin(��);
				up.vm = down.vm = -Vs * cos(��);
				//up.u = V1t * cos(��) + V1n * sin(��);
				//up.v = V1t * sin(��) - V1n * cos(��);

				//down.u = V2t * cos(��) + V2n * sin(��);
				//down.v = V2t * sin(��) - V2n * cos(��);
				//A[i][n1].u = unb;
				//A[i][n1].�� = ��2 * pow((ab * ab / (c21 * c21)), 1 / (�� - 1));
				//A[i][n1].p = ab * ab * down.�� / ��;
				double p20 = down.p;
				double ��20 = down.��;
				down.p = p2;
				down.�� = ��2 * pow(down.p / p20 ,1 / ��);
				double c = V2n * V2n / 2 + �� * down.p / ((�� - 1) * down.��);
				down.u = sqrt(2 * (c - �� * down.p / ((�� - 1) * down.��)));

				//down.p = up.p;
				//down.�� = up.��;
				//down.u = up.u;
				//down.v = up.v;
			}
			else
			{
				up.um = down.um = Vs * sin(-��);
				up.vm = down.vm = Vs * cos(-��);
				//up.u = V1t * cos(-��) + V1n * sin(-��);
				//up.v = -V1t * sin(-��) + V1n * cos(-��);
				down.u = V2t * cos(-��) + V2n * sin(-��);
				down.v = -V2t * sin(-��) + V2n * cos(-��);

			}
			if (FlowType == "intersection" && i == 0)
			{
				up.um = down.um = -down.um;
				up.vm = down.vm = -down.vm;
			}
			if (FlowType == "intersection")
				up.um = down.um = 0;
			else if (FlowType == "oblique" || FlowType == "normal")
				up.vm = down.vm = 0;

			up.step = down.step = step + 1;
			if (Ar[i][j].p < Ar[m][n].p)
			{
				A[i][j] = up;
				A[m][n] = down;
			}
			else
			{
				A[i][j] = down;
				A[m][n] = up;
			}
		}
	}
}
double get_��(double x1, double y1, double x2, double y2)//��ֱ����x��ļн�
{
	double ��;
	if (x1 == x2)
	{
		//if (y2 < y1)
		//	�� = -pi / 2;
		//else
		//	�� = pi / 2;
		�� = pi / 2;
	}
	else if (y1 == y2)
		�� = 0;
	else
	{
		�� = atan(abs((y2 - y1) / (x2 - x1)));
		if ((x2 > x1 && y2 < y1) || (x2 < x1 && y2 > y1))
			�� = -��;
	}
	return ��;
}
void update_Vm()
{
	extern double xL, xR, yU, yD;
	extern vector<vector<mesh>> A;
	extern vector<mesh>A0;
	int i, j, k, n;
	double d1, d2;
	extern int step;
	int s = 1, m;
	double um, vm;
	if (FlowType == "oblique" || FlowType == "normal")
	{
		double x;
		if (FlowType == "oblique")
		{
			A[0][231].um = A[0][231].vm = 0;
			A[1][1419].um = A[1][1419].vm = 0;
		}

		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "SHOCK")
					continue;
				um = 0, vm = 0, x = 0;
				for (k = 0; k < A[i][j].moveConnct.size(); k++)
				{
					m = A[i][j].moveConnct[k];
					um += A[i][m].um;
					vm += A[i][m].vm;
					x += A[i][m].x;
				}
				um = um / A[i][j].moveConnct.size();
				vm = vm / A[i][j].moveConnct.size();
				x = x / A[i][j].moveConnct.size();
				if (i == 0)
				{
					d1 = abs(xL - A[i][j].x);
					d2 = abs(x - xL);
					A[i][j].um = d1 / d2 * um;
					A[i][j].vm = d1 / d2 * vm;
					//A[i][j].um = A[i][j].x / x * um;
					//A[i][j].vm = A[i][j].x / x * vm;
				}
				if (i == 1)
				{
					d1 = abs(x - A[i][j].x);
					d2 = abs(x - xR);
					A[i][j].um = (d2 - d1) / d2 * um;
					A[i][j].vm = (d2 - d1) / d2 * vm;

					//A[i][j].um = (A0[Xnum - 1].x - A[i][j].x) / (A0[Xnum - 1].x - x) * um;
					//A[i][j].vm = (A0[Xnum - 1].x - A[i][j].x) / (A0[Xnum - 1].x - x) * vm;
				}
			}
		}
		//for (i = 0; i < A.size(); i++)
		//{
		//	for (j = 0; j < A[i].size(); j++)
		//	{
		//		if (A[i][j].type == "SHOCK")
		//			continue;
		//		um = 0, vm = 0, x = 0;
		//		for (k = 0; k < A[i][j].moveConnct.size(); k++)
		//		{
		//			m = A[i][j].moveConnct[k];
		//			um += A[i][m].um;
		//			vm += A[i][m].vm;
		//			x += A[i][m].x;
		//		}
		//		um = um / A[i][j].moveConnct.size();
		//		vm = vm / A[i][j].moveConnct.size();
		//		x = x / A[i][j].moveConnct.size();
		//		if (i == 0)
		//		{
		//			A[i][j].um = A[i][j].x / x * um;
		//			A[i][j].vm = A[i][j].x / x * vm;
		//		}
		//		if (i == 1)
		//		{
		//			A[i][j].um = (A0[Xnum - 1].x - A[i][j].x) / (A0[Xnum - 1].x - x) * um;
		//			A[i][j].vm = (A0[Xnum - 1].x - A[i][j].x) / (A0[Xnum - 1].x - x) * vm;
		//		}
		//	}
		//}
	}
	else if (FlowType == "intersection")
	{
		double x, y, x1, y1, x2, y2;
		double um1, vm1, um2, vm2;
		double d1, d2;
		int n1, n2;
		for (i = 0; i < A.size(); i++)
		{
			A[i][0].um = A[i][0].vm = 0;
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "SHOCK")
				{
					continue;
				}
				else if (A[i][j].type == "DISCON")
				{
					continue;
				}
				else if (A[i][j].type == "CENTER")
				{
					A[i][j].um = A[i][j].vm = 0;
				}
				else if (A[i][j].y == yU || A[i][j].y == yD)
					A[i][j].um = A[i][j].vm = 0;
				else
				{

					if (i == 0 || i == 3 || i == 4)
					{
						A[i][j].vm = 0;
						for (k = 0; k < A[i][j].moveConnct.size(); k++)
						{
							m = A[i][j].moveConnct[k];
							if (A[i][m].y > A[i][j].y)
								y1 = A[i][m].y, vm1 = A[i][m].vm;
							else
								y2 = A[i][m].y, vm2 = A[i][m].vm;
						}
					}

					else if (i == 1)
					{
						y1 = yU;
						y2 = A[i][A[i][j].moveConnct[0]].y;
						vm1 = 0;
						vm2 = A[i][A[i][j].moveConnct[0]].vm;
					}
					else
					{
						y1 = A[i][A[i][j].moveConnct[0]].y;
						y2 = yD;
						vm1 = A[i][A[i][j].moveConnct[0]].vm;
						vm2 = 0;
					}
					y = A[i][j].y;
					double con = (y1 - y) / (y - y2);
					A[i][j].vm = (y1 + vm1 + con * (y2 + vm2) - (y + con * y)) / (con + 1);
					A[i][j].um = 0;
				}



			}

		}


	}

}
void clear_Vm()//���������˶��ٶ�
{
	extern vector<vector<mesh>> A;
	int i, j;
#pragma omp parallel
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			A[i][j].um = 0;
			A[i][j].vm = 0;
		}
	}
}
void update_bound()
{
	extern vector<vector<mesh>> A;

	int i, j, k, m, n;
	if (methodType == "C")
	{
		if (FlowType == "normal")
		{
			using namespace Normal;
#pragma omp parallel

			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j].type == "IN" || A[i][j].type == "SHOCK")
						continue;
					if (i == 0)
					{
						A[i][j].�� = ��1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
					if (i == 1)
					{
						A[i][j].�� = ��2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;
					}

				}
			}
		}
		if (FlowType == "oblique")
		{
			using namespace Oblique;
#pragma omp parallel
			for (i = 0; i < A.size(); i++)
			{


				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j].type == "IN" || A[i][j].type == "SHOCK")
						continue;
					if (i == 0)
					{
						A[i][j].�� = ��1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
					if (i == 1)
					{
						A[i][j].�� = ��2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;
					}
					if (A[i][j].type == "R")
					{
						A[i][j].�� = A[i][j - 1].��;
						A[i][j].u = A[i][j - 1].u;
						A[i][j].v = A[i][j - 1].v;
						A[i][j].p = A[i][j - 1].p;
					}

				}
			}
		}
	}
	else if (methodType == "F")
	{
		if (FlowType == "normal")
		{
			using namespace Normal;
			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					//��1 = 3.8571, p1 = 10.3333, u1 = 2.6293, v1 = 0;
					//��2 = 0.5 + 0.1 * sin(5 * A[i][j].x), p2 = 1, u2 = 0, v2 = 0;
					��2 = 3 +1*sin(5 * A[i][j].x);
					//��2 = 1 + 0.1 * sin(5 * A[i][j].x);
					if (i == 1)
					{
						A[i][j].�� = ��2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;
						continue;
					}
					if (A[i][j].type == "IN" || A[i][j].type == "SHOCK" || A[i][j].type == "DISCON")
						continue;

					else if (A[i][j].type == "U" || A[i][j].type == "D")
					{
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							m = A[i][j].neibor[k];
							if (A[i][m].type == "SHOCK")
								continue;
							else if (abs(A[i][j].y - A[i][m].y) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j].�� = A[i][n].��;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
					else if (A[i][j].type == "L")
					{
						A[i][j].�� = ��1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
					else if (A[i][j].type == "R")
					{
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							m = A[i][j].neibor[k];
							if (A[i][m].type == "SHOCK")
								continue;
							else if (abs(A[i][j].x - A[i][m].x) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j].�� = A[i][n].��;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
				}
			}

		}
		else if (FlowType == "oblique")
		{
			using namespace Oblique;

			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j].type == "IN" || A[i][j].type == "SHOCK" || A[i][j].type == "DISCON")
						continue;
					else if (A[i][j].type == "U")
					{
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							m = A[i][j].neibor[k];
							if (A[i][m].type == "SHOCK")
								continue;
							else if (abs(A[i][j].y - A[i][m].y) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j].�� = A[i][n].��;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
					else if (A[i][j].type == "D")
					{
						if (i == 0)
						{
							A[i][j].�� = ��1;
							A[i][j].u = u1;
							A[i][j].v = v1;
							A[i][j].p = p1;
						}
						else
						{
							A[i][j].�� = ��2;
							A[i][j].u = u2;
							A[i][j].v = v2;
							A[i][j].p = p2;

						}

					}
					else if (A[i][j].type == "L")
					{
						A[i][j].�� = ��1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
					else if (A[i][j].type == "R")
					{
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							m = A[i][j].neibor[k];
							if (A[i][m].type == "SHOCK")
								continue;
							else if (abs(A[i][j].x - A[i][m].x) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j].�� = A[i][n].��;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
				}
			}
		}
		else if (FlowType == "intersection")
		{
			using namespace ShockwaveCross;
			for (i = 0; i < A.size(); i++)
			{

				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j].type == "IN")
						continue;
					else if (A[i][j].type == "U" || A[i][j].type == "D")
					{
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							m = A[i][j].neibor[k];
							if (A[i][m].type == "SHOCK")
								continue;
							else if (abs(A[i][j].y - A[i][m].y) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j].�� = A[i][n].��;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
					else if (A[i][j].type == "L")
					{
						if (i == 0)
						{
							A[i][j].�� = ��1;
							A[i][j].u = u1;
							A[i][j].v = v1;
							A[i][j].p = p1;
						}
						else if (i == 1)
						{
							A[i][j].�� = ��2;
							A[i][j].u = u2;
							A[i][j].v = v2;
							A[i][j].p = p2;
						}
						else if (i == 2)
						{
							A[i][j].�� = ��3;
							A[i][j].u = u3;
							A[i][j].v = v3;
							A[i][j].p = p3;

						}
					}
					else if (A[i][j].type == "R")
					{
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							m = A[i][j].neibor[k];
							if (A[i][m].type == "SHOCK")
								continue;
							else if (abs(A[i][j].x - A[i][m].x) > 1e-10)
							{
								n = m;
								break;
							}
							else
								n = m;
						}
						A[i][j].�� = A[i][n].��;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
				}
			}

		}
	}
}
void update_bound_shockwaveCross()
{
	using namespace ShockwaveCross;
	extern vector<vector <mesh>> A;
	int i, j;
#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].type == "IN")
				continue;
			else if (A[i][j].type == "U")
			{
				A[i][j].�� = ��2;
				A[i][j].u = u2;
				A[i][j].v = v2;
				A[i][j].p = p2;
			}
			else if (A[i][j].type == "L")
			{
				A[i][j].�� = ��1;
				A[i][j].u = u1;
				A[i][j].v = v1;
				A[i][j].p = p1;

			}
			else if (A[i][j].type == "D")
			{
				A[i][j].�� = ��3;
				A[i][j].u = u3;
				A[i][j].v = v3;
				A[i][j].p = p3;
			}
			else if (A[i][j].type == "R")
			{
				A[i][j].�� = A[i][j - 1].��;
				A[i][j].u = A[i][j - 1].u;
				A[i][j].v = A[i][j - 1].v;
				A[i][j].p = A[i][j - 1].p;
			}

		}
	}
}
void update_bound_Prandtl_Meyer()
{
	using namespace Prandtl_Meyer;
	extern vector<vector <mesh>> A;
	int i, j;
#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].type == "IN")
				continue;
			else if (A[i][j].type == "L")
			{
				A[i][j].�� = ��1;
				A[i][j].u = u1;
				A[i][j].v = v1;
				A[i][j].p = p1;
			}
			else if (A[i][j].type == "DL")
			{
				A[i][j].�� = ��1;
				A[i][j].u = u1;
				A[i][j].v = v1;
				A[i][j].p = p1;
			}
			else if (A[i][j].type == "DR")
			{
				A[i][j].�� = ��2;
				A[i][j].u = u2;
				A[i][j].v = v2;
				A[i][j].p = p2;
			}
			else if (A[i][j].type == "R")
			{
				A[i][j].�� = A[i][j - 1].��;
				A[i][j].u = A[i][j - 1].u;
				A[i][j].v = A[i][j - 1].v;
				A[i][j].p = A[i][j - 1].p;
			}
		}
	}
}



//void get_F()
//{
//	extern vector <mesh> A;
//	extern double F[Pnum][4];
//	for (int i = 0; i < Pnum; i++)
//	{
//		F[i][0] = A[i].��*A[i].u;
//		F[i][1] = A[i].��*A[i].u*A[i].u + A[i].p;
//		F[i][2] = A[i].��*A[i].u*A[i].v;
//		F[i][3] = A[i].u*(A[i].p / (�� - 1) + 0.5*A[i].��*(A[i].u*A[i].u + A[i].v*A[i].v) + A[i].p);
//	}
//}
//void get_G()
//{
//	extern vector <mesh> A;
//	extern double G[Pnum][4];
//	for (int i = 0; i < Pnum; i++)
//	{
//		G[i][0] = A[i].��*A[i].v; extern vector<vector<vector <double>>> Utr;
//		G[i][1] = A[i].��*A[i].u *A[i].v;
//		G[i][2] = A[i].��*A[i].u*A[i].v + A[i].p;
//		G[i][3] = A[i].v*(A[i].p / (�� - 1) + 0.5*A[i].��*(A[i].u*A[i].u + A[i].v*A[i].v) + A[i].p);
//	}
//}

//void update_INtrans()
//{
//	extern vector <mesh> A0;
//	extern vector<vector<vector <double>>> U;
//	extern vector<vector<vector <double>>> Utr;
//	extern int i;
//
//#pragma omp parallel 
//
//	for (int i = 0; i < A0.size(); i++)
//	{
//		if (A0[i].type != "IN")
//			continue;
//		if (A0[i].neibor.size() != 3)
//		{
//			Utr[i][0][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][0][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][0][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][0][3] = U[0][i][3] * A0[i].J[0];
//		}
//		else
//		{
//			Utr[i][0][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][0][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][0][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][0][3] = U[0][i][3] * A0[i].J[0];
//
//			Utr[i][1][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][1][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][1][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][1][3] = U[0][i][3] * A0[i].J[0];
//
//			Utr[i][2][0] = U[0][i][0] * A0[i].J[0];
//			Utr[i][2][1] = U[0][i][1] * A0[i].J[0];
//			Utr[i][2][2] = U[0][i][2] * A0[i].J[0];
//			Utr[i][2][3] = U[0][i][3] * A0[i].J[0];
//		}
//	}
//}
double compute_res()//����в�
{
	extern vector<vector <mesh>> A;
	extern  vector < vector <mesh>> Ar;
	extern double dt;
	int i, j;
	double res = 0;
	int n = A.size() * A[0].size();
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			res += abs(A[i][j].�� - Ar[i][j].��) / Ar[i][j].��;
		}
	}
	return res / n;
}
void record()
{
	extern vector<vector <mesh>> A;
	extern  vector < vector <mesh>> Ar;
	extern double t_sim;
	int i, j;
	if (t_sim == 0)
	{
		for (i = 0; i < A.size(); i++)
		{
			Ar.push_back(A[i]);
		}
	}
	else
	{
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				Ar[i][j] = A[i][j];
			}
		}

	}
}
void judge()//�ж���Χ�������Ƿ�ı䣬���û�ı�Ļ���������������٣���
{
	extern vector <vector<mesh>> A;
	extern vector<vector <mesh>> Ar;
	int judge;
	extern int step;
	if (step == 0)
		std::cout << "���㿪ʼ��" << std::endl;
	else
		for (int i = 0; i < A.size(); i++)
		{
			for (int j = 0; j < A[i].size(); j++)
			{
				judge = 0;
				for (int k = 0; k < A[i][j].neibor.size(); k++)
				{
					int n = A[i][j].neibor[k];
					if (A[i][n].p != Ar[i][n].p && A[i][n].u != Ar[i][n].u && A[i][n].v != Ar[i][n].v && A[i][n].�� != Ar[i][n].��)
						judge = 1;
				}
				if (judge == 0)
					A[i][j].change = "N";
				else
					A[i][j].change = "Y";
			}
		}
}

void update_IN()//�����ڲ���
{
	extern vector<vector <mesh>> A;
	extern vector<vector <mesh>> Ar;
	double U[4];
	double temp[4] = { 0 };
	extern double dt;
	int i, j, k;
	Flux F1, F2, G1, G2;
	Flux F1_, F2_, G1_, G2_;
	int n1, n2, n3, n4;

	double d1, d2;
	//#pragma omp parallel for private(F1,F2,G1,G2,temp)
	for (i = 0; i < Ar.size(); i++)
	{
		for (j = 0; j < Ar[i].size(); j++)
		{
			if (Ar[i][j].type == "L" || Ar[i][j].type == "R" || Ar[i][j].type == "U" || Ar[i][j].type == "D" /*|| A[i][j].change == "N"*/ || A[i][j].type == "SHOCK")
				continue;
			U[0] = Ar[i][j].��;
			U[1] = Ar[i][j].�� * Ar[i][j].u;
			U[2] = Ar[i][j].�� * Ar[i][j].v;
			U[3] = 0.5 * Ar[i][j].�� * (Ar[i][j].u * Ar[i][j].u + Ar[i][j].v * Ar[i][j].v) + Ar[i][j].p / (�� - 1);

			if (Ar[i][j].neibor.size() == 3)
			{
				F1 = F2 = G1 = G2 = { 0 };

				for (k = 0; k < 12; k++)
				{
					n1 = method[k][0];
					n2 = method[k][1];
					n3 = method[k][2];
					n4 = method[k][3];

					F1_ = HLLC_��2(Ar[i][Ar[i][j].neibor[n3]], Ar[i][j], Ar[i][j], k);
					F2_ = HLLC_��2(Ar[i][j], Ar[i][Ar[i][j].neibor[n1]], Ar[i][j], k);
					G1_ = HLLC_��2(Ar[i][Ar[i][j].neibor[n4]], Ar[i][j], Ar[i][j], k);
					G2_ = HLLC_��2(Ar[i][j], Ar[i][Ar[i][j].neibor[n2]], Ar[i][j], k);
					F1.f1 += F1_.f1, F1.f2 += F1_.f2, F1.f3 += F1_.f3, F1.f4 += F1_.f4;
					F2.f1 += F2_.f1, F2.f2 += F2_.f2, F2.f3 += F2_.f3, F2.f4 += F2_.f4;
					G1.f1 += G1_.f1, G1.f2 += G1_.f2, G1.f3 += G1_.f3, G1.f4 += G1_.f4;
					G2.f1 += G2_.f1, G2.f2 += G2_.f2, G2.f3 += G2_.f3, G2.f4 += G2_.f4;
				}
				F1.f1 /= k, F1.f2 /= k, F1.f3 /= k, F1.f4 /= k;
				F2.f1 /= k, F2.f2 /= k, F2.f3 /= k, F2.f4 /= k;
				G1.f1 /= k, G1.f2 /= k, G1.f3 /= k, G1.f4 /= k;
				G2.f1 /= k, G2.f2 /= k, G2.f3 /= k, G2.f4 /= k;

				U[0] = U[0] - dt * (F2.f1 - F1.f1) - dt * (G2.f1 - G1.f1);
				U[1] = U[1] - dt * (F2.f2 - F1.f2) - dt * (G2.f2 - G1.f2);
				U[2] = U[2] - dt * (F2.f3 - F1.f3) - dt * (G2.f3 - G1.f3);
				U[3] = U[3] - dt * (F2.f4 - F1.f4) - dt * (G2.f4 - G1.f4);

				//if (k == 0)
				//{
				//	F1 = HLLC_��(Ar[i][Ar[i][j].neibor[2]], Ar[i][j], Ar[i][j], k);
				//	F2 = HLLC_��(Ar[i][j], Ar[i][Ar[i][j].neibor[0]], Ar[i][j], k);
				//	G1 = HLLC_��(Ar[i][Ar[i][j].neibor[2]], Ar[i][j], Ar[i][j], k);
				//	G2 = HLLC_��(Ar[i][j], Ar[i][Ar[i][j].neibor[1]], Ar[i][j], k);
				//	d1 = distance(Ar[i][Ar[i][j].neibor[2]], Ar[i][Ar[i][j].neibor[0]]) / 2;
				//	d2 = distance(Ar[i][Ar[i][j].neibor[2]], Ar[i][Ar[i][j].neibor[1]]) / 2;
				//}
				//else if (k == 1)
				//{
				//	F1 = HLLC_��(Ar[i][Ar[i][j].neibor[0]], Ar[i][j], Ar[i][j], k);
				//	F2 = HLLC_��(Ar[i][j], Ar[i][Ar[i][j].neibor[1]], Ar[i][j], k);
				//	G1 = HLLC_��(Ar[i][Ar[i][j].neibor[0]], Ar[i][j], Ar[i][j], k);
				//	G2 = HLLC_��(Ar[i][j], Ar[i][Ar[i][j].neibor[2]], Ar[i][j], k);
				//	d1 = distance(Ar[i][Ar[i][j].neibor[0]], Ar[i][Ar[i][j].neibor[1]]) / 2;
				//	d2 = distance(Ar[i][Ar[i][j].neibor[0]], Ar[i][Ar[i][j].neibor[2]]) / 2;
				//}
				//else
				//{
				//	F1 = HLLC_��(Ar[i][Ar[i][j].neibor[1]], Ar[i][j], Ar[i][j], k);
				//	F2 = HLLC_��(Ar[i][j], Ar[i][Ar[i][j].neibor[2]], Ar[i][j], k);
				//	G1 = HLLC_��(Ar[i][Ar[i][j].neibor[1]], Ar[i][j], Ar[i][j], k);
				//	G2 = HLLC_��(Ar[i][j], Ar[i][Ar[i][j].neibor[0]], Ar[i][j], k);
				//	d1 = distance(Ar[i][Ar[i][j].neibor[1]], Ar[i][Ar[i][j].neibor[2]]) / 2;
				//	d2 = distance(Ar[i][Ar[i][j].neibor[1]], Ar[i][Ar[i][j].neibor[0]]) / 2;

				//}
				//temp[0] += U[0] - dt * (F2.f1 - F1.f1) / d1 - dt * (G2.f1 - G1.f1) / d2;
				//temp[1] += U[1] - dt * (F2.f2 - F1.f2) / d1 - dt * (G2.f2 - G1.f2) / d2;
				//temp[2] += U[2] - dt * (F2.f3 - F1.f3) / d1 - dt * (G2.f3 - G1.f3) / d2;
				//temp[3] += U[3] - dt * (F2.f4 - F1.f4) / d1 - dt * (G2.f4 - G1.f4) / d2;
			//}
			//U[0] = temp[0] / 3;
			//U[1] = temp[1] / 3;
			//U[2] = temp[2] / 3;
			//U[3] = temp[3] / 3;
			}

			else if (A[i][j].neibor.size() == 4)
			{
				F1 = HLLC_��2(Ar[i][Ar[i][j].neibor[2]], Ar[i][j], Ar[i][j], 0);
				F2 = HLLC_��2(Ar[i][j], Ar[i][Ar[i][j].neibor[0]], Ar[i][j], 0);
				G1 = HLLC_��2(Ar[i][Ar[i][j].neibor[3]], Ar[i][j], Ar[i][j], 0);
				G2 = HLLC_��2(Ar[i][j], Ar[i][Ar[i][j].neibor[1]], Ar[i][j], 0);

				U[0] = U[0] - dt * (F2.f1 - F1.f1) - dt * (G2.f1 - G1.f1);
				U[1] = U[1] - dt * (F2.f2 - F1.f2) - dt * (G2.f2 - G1.f2);
				U[2] = U[2] - dt * (F2.f3 - F1.f3) - dt * (G2.f3 - G1.f3);
				U[3] = U[3] - dt * (F2.f4 - F1.f4) - dt * (G2.f4 - G1.f4);
			}
			else if (A[i][j].neibor.size() > 4)
			{
				F1 = HLLC_��(A[i][A[i][j].neibor1[2]], A[i][j], A[i][j], 0);
				F2 = HLLC_��(A[i][j], A[i][A[i][j].neibor1[0]], A[i][j], 0);
				G1 = HLLC_��(A[i][A[i][j].neibor1[3]], A[i][j], A[i][j], 0);
				G2 = HLLC_��(A[i][j], A[i][A[i][j].neibor1[1]], A[i][j], 0);
				U[0] = U[0] - dt * (F2.f1 - F1.f1 + G2.f1 - G1.f1);
				U[1] = U[1] - dt * (F2.f2 - F1.f2 + G2.f2 - G1.f2);
				U[2] = U[2] - dt * (F2.f3 - F1.f3 + G2.f3 - G1.f3);
				U[3] = U[3] - dt * (F2.f4 - F1.f4 + G2.f4 - G1.f4);
			}
			else if (A[i][j].neibor.size() == 2)
			{
				if (A[i][j].type == "SHOCK")
				{
					double m, n;
					m = A[i][j].neibor[0];
					n = A[i][j].neibor[1];

					A[i][j].�� = (A[i][m].�� + A[i][n].��) / 2;
					A[i][j].u = (A[i][m].u + A[i][n].u) / 2;
					A[i][j].v = (A[i][m].v + A[i][n].v) / 2;
					A[i][j].p = (A[i][m].p + A[i][n].p) / 2;
				}
				else
					continue;
			}
			else
				std::cout << i << "   " << j << "   " << "somthing wrong in updat_u" << std::endl;
			A[i][j].�� = U[0];
			A[i][j].u = U[1] / U[0];
			A[i][j].v = U[2] / U[0];
			A[i][j].p = (�� - 1) * (U[3] - 0.5 * A[i][j].�� * (A[i][j].u * A[i][j].u + A[i][j].v * A[i][j].v));
		}
	}

}

//void choose_U(int i)
//{
//	extern vector <mesh> A;
//	extern vector <mesh> A;
//	extern vector<vector <double>> U;
//	extern double U0[Pnum][4];
//	extern double U1[Pnum][4];
//	extern double U2[Pnum][4];
//	double ��00, u00, v00, p00;
//	double ��, u, v, p;
//	double ��0, u0, v0, p0;
//	double ��1, u1, v1, p1;
//	double ��2, u2, v2, p2;
//	��00 = A[i].��;
//	u00 = A[i].u;
//	v00 = A[i].v;
//	p00 = A[i].p;
//
//	��0 = U0[i][0] / A[i].J[0];
//	u0 = U0[i][1] / U0[i][0];
//	v0 = U0[i][2] / U0[i][0];
//	p0 = (�� - 1)*(U0[i][3] / A[i].J[0] - 0.5*A[i].��*(A[i].u*A[i].u + A[i].v*A[i].v));
//
//	��1 = U1[i][0] / A[i].J[1];
//	u1 = U1[i][1] / U1[i][0];
//	v1 = U1[i][2] / U1[i][0];
//	p1 = (�� - 1)*(U1[i][3] / A[i].J[1] - 0.5*A[i].��*(A[i].u*A[i].u + A[i].v*A[i].v));
//	//��0 = ��1;
//	//u0 = u1;
//	//v0 = v1;
//	//p0 = p1;
//
//	��2 = U2[i][0] / A[i].J[2];
//	u2 = U2[i][1] / U2[i][0];
//	v2 = U2[i][2] / U2[i][0];
//	p2 = (�� - 1)*(U2[i][3] / A[i].J[1] - 0.5*A[i].��*(A[i].u*A[i].u + A[i].v*A[i].v));
//	�� = absmax(absmax(��0, ��1), ��2);
//	u = absmin(absmin(u0, u1), u2);
//	v = absmax(absmax(v0, v1), v2);
//	p = absmax(absmax(p0, p1), p2);
//	U[i][0] = ��;
//	U[i][1] = �� * u;
//	U[i][2] = �� * v;
//	U[i][3] = 0.5*��*(u*u + v * v) + p / (�� - 1);
//}
//void choose_U(int i)
//{
//	extern vector <mesh> A;
//
//	extern double U[Pnum][4];
//	extern double U0[Pnum][4];
//	extern double U1[Pnum][4];
//	extern double U2[Pnum][4];
//
//	double u0 = abs(U0[i][1] / U0[i][0]);
//	double u1 = abs(U1[i][1] / U1[i][0]);
//	double u2 = abs(U2[i][1] / U2[i][0]);
//	if (u0 < u1&&u0 < u2)
//	{
//		U[i][0] = U0[i][0] / A[i].J[0];
//		U[i][1] = U0[i][1] / A[i].J[0];
//		U[i][2] = U0[i][2] / A[i].J[0];
//		U[i][3] = U0[i][3] / A[i].J[0];
//	}
//	else if (u1 < u0&&u1 < u2)
//	{
//		U[i][0] = U1[i][0] / A[i].J[1];
//		U[i][1] = U1[i][1] / A[i].J[1];
//		U[i][2] = U1[i][2] / A[i].J[1];
//		U[i][3] = U1[i][3] / A[i].J[1];
//	}
//	else
//	{
//		U[i][0] = U2[i][0] / A[i].J[2];
//		U[i][1] = U2[i][1] / A[i].J[2];
//		U[i][2] = U2[i][2] / A[i].J[2];
//		U[i][3] = U2[i][3] / A[i].J[2];
//	}
//
//}
double absmax(double a, double b)
{
	if (abs(a) >= abs(b))
		return a;
	else
		return b;
}
double absmin(double a, double b)
{
	if (abs(a) <= abs(b))
		return a;
	else
		return b;

}
double max(double a, double b)
{
	if (a >= b)
		return a;
	else
		return b;
}
double min(double a, double b)
{
	if (a <= b)
		return a;
	else
		return b;

}
double get_��(mesh A, mesh B)//��������������x��ļн�
{
	double dy = abs(A.y - B.y);
	double dx = abs(A.x - B.x);
	double �� = atan(dy / dx);
	return ��;
}
void reorder_neighbor()
{
	extern vector<vector <mesh>> A;
	double maxy, miny, maxx, minx;
	int t;
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].type != "IN" || A[i][j].neibor.size() < 3)
				continue;
			if (A[i][j].neibor.size() == 3)
			{
				int n1 = A[i][j].neibor[0];
				int n2 = A[i][j].neibor[1];
				int n3 = A[i][j].neibor[2];
				if (j - 1 == n2)
					t = n2, n2 = n3, n3 = t;
				A[i][j].neibor[0] = n1;
				A[i][j].neibor[1] = n2;
				A[i][j].neibor[2] = n3;
			}
			if (A[i][j].neibor.size() == 4)
			{
				int n1 = A[i][j].neibor[0];
				int n2 = A[i][j].neibor[1];
				int n3 = A[i][j].neibor[2];
				int n4 = A[i][j].neibor[3];
				double maxy = max(max(max(A[i][n1].y, A[i][n2].y), A[i][n3].y), A[i][n4].y);
				if (A[i][n1].y == maxy)
					t = n1, n1 = n2, n2 = t;
				else if (A[i][n3].y == maxy)
					t = n3, n3 = n2, n2 = t;
				else if (A[i][n4].y == maxy)
					t = n4, n4 = n2, n2 = t;
				double miny = min(min(A[i][n1].y, A[i][n3].y), A[i][n4].y);
				if (A[i][n1].y == miny)
					t = n1, n1 = n4, n4 = t;
				else if (A[i][n3].y == miny)
					t = n3, n3 = n4, n4 = t;
				maxx = max(A[i][n1].x, A[i][n3].x);
				if (A[i][n3].x == maxx)
					t = n1, n1 = n3, n3 = t;
				A[i][j].neibor[0] = n1;
				A[i][j].neibor[1] = n2;
				A[i][j].neibor[2] = n3;
				A[i][j].neibor[3] = n4;
			}
			if (A[i][j].neibor.size() > 4)
			{
				int n1 = A[i][j].neibor1[0];
				int n2 = A[i][j].neibor1[1];
				int n3 = A[i][j].neibor1[2];
				int n4 = A[i][j].neibor1[3];
				double maxy = max(max(max(A[i][n1].y, A[i][n2].y), A[i][n3].y), A[i][n4].y);
				if (A[i][n1].y == maxy)
					t = n1, n1 = n2, n2 = t;
				else if (A[i][n3].y == maxy)
					t = n3, n3 = n2, n2 = t;
				else if (A[i][n4].y == maxy)
					t = n4, n4 = n2, n2 = t;
				double miny = min(min(A[i][n1].y, A[i][n3].y), A[i][n4].y);
				if (A[i][n1].y == miny)
					t = n1, n1 = n4, n4 = t;
				else if (A[i][n3].y == miny)
					t = n3, n3 = n4, n4 = t;
				maxx = max(A[i][n1].x, A[i][n3].x);
				if (A[i][n3].x == maxx)
					t = n1, n1 = n3, n3 = t;
				A[i][j].neibor1[0] = n1;
				A[i][j].neibor1[1] = n2;
				A[i][j].neibor1[2] = n3;
				A[i][j].neibor1[3] = n4;
			}

		}

	}
}
double area(mesh A, mesh B, mesh C, mesh D)//�������ĵ㹹���ı������ 
{
	return 0.5* abs(A.x * B.y + B.x * C.y + C.x * D.y + D.x * A.y - B.x * A.y - C.x * B.y - D.x * C.y - A.x * D.y);
}

void movemesh()//�����ƶ�
{
	extern double dt;
	extern vector<vector<mesh>> A;
	int i, j;
#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			A[i][j].x = A[i][j].x + A[i][j].um * dt;
			A[i][j].y = A[i][j].y + A[i][j].vm * dt;
		}
	}
	//extern double t_sim;
	//extern int step;
	//extern vector <mesh> A0;
	//double L0 = A0[Xnum - 1].x0 - A0[0].x0;
	//double L1 = L0 / 2;
	////double L = L1 + 0.5*L1*sin(20 * pi*t_sim);
	//double L = L1 + 0.5 * L1 * sin(pi * step / 2000);
	//for (int i = 0; i < A0.size(); i++)
	//{
	//	if (A0[i].x0 < L1)
	//		A0[i].x = A0[i].x0 * (L / L1);
	//	else
	//		A0[i].x = A0[Xnum - 1].x0 - (A0[Xnum - 1].x0 - A0[i].x0) * ((L0 - L) / L1);
	//}
}
mesh getCrossPoint(Line L1, Line L2)
{
	mesh L;
	L.x = (L2.B * L1.C - L1.B * L2.C) / (L2.A * L1.B - L1.A * L2.B);
	L.y = (L2.A * L1.C - L1.A * L2.C) / (L1.A * L2.B - L2.A * L1.B);
	return L;
}

void findNeiborSec()
{
	extern vector<vector<mesh>> A;
	int i, j, k, m;
//#pragma omp parallel

	for (i = 0; i < A.size(); i++)
	{
		if (A.size() == 1)
			continue;
		for (j = 0; j < A[i].size(); j++)
		{
			for (k = i + 1; k < A.size(); k++)
			{
				for (m = 0; m < A[k].size(); m++)
				{
					if (A[i][j].x == A[k][m].x && A[i][j].y == A[k][m].y)
					{
						A[i][j].neiborsec = k;
						A[i][j].neiborsec_ad = m;
						A[k][m].neiborsec = i;
						A[k][m].neiborsec_ad = j;
					}
				}
			}

		}
	}
}
Line getLine(mesh A, mesh B)
{
	Line L;
	if (A.x == B.x)
	{
		if (A.y == B.y)
		{
			std::cout << "something wrong in getLine !" << std::endl;
			L.A = 0, L.B = 0, L.C = 0;
		}
		else
			L.A = 1, L.B = 0, L.C = -A.x;
	}
	else if (A.y == B.y)
	{
		L.A = 0, L.B = 1, L.C = -A.y;
	}
	else
		L.A = (B.y - A.y) / (B.x - A.x), L.B = -1, L.C = (B.x * A.y - A.x * B.y) / (B.x - A.x);
	return L;
}
Line getLine(double ��, mesh A)
{
	double k = tan(��);
	double b = A.y - k * A.x;
	Line L;
	L.A = k, L.B = -1, L.C = b;
	return L;

}