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
void get_dt()//求Δt
{
	extern vector<vector <mesh>> A;
	double maxξ = 0, maxη = 0;
	extern double dt;
	double t;
	int i, j, k;
	double max1, max2;
	double Sξ, Sη, c, uξ, uη;

	dt = t_end;
	max1 = max2 = 0;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].type != "IN")
				continue;
			maxξ = maxη = 0;
			if (A[i][j].neibor.size() > 3)
			{
				Sξ = sqrt(A[i][j].ξx[0] * A[i][j].ξx[0] + A[i][j].ξy[0] * A[i][j].ξy[0]);
				Sη = sqrt(A[i][j].ηx[0] * A[i][j].ηx[0] + A[i][j].ηy[0] * A[i][j].ηy[0]);
				c = sqrt(γ * A[i][j].p / A[i][j].ρ);
				uξ = A[i][j].u * A[i][j].ξx[0] + A[i][j].v * A[i][j].ξy[0];
				uη = A[i][j].u * A[i][j].ηx[0] + A[i][j].v * A[i][j].ηy[0];
				//maxξ = max(Sξ, uξ);
				//maxη = max(Sη, uη);
				maxξ = abs(uξ) + c * Sξ;
				maxη = abs(uξ) + c * Sη;
				max1 = max(max1, maxξ);
				max2 = max(max2, maxη);

			}
			else
			{
				for (k = 0; k < 3; k++)
				{
					Sξ = sqrt(A[i][j].ξx[k] * A[i][j].ξx[k] + A[i][j].ξy[k] * A[i][j].ξy[k]);
					Sη = sqrt(A[i][j].ηx[k] * A[i][j].ηx[k] + A[i][j].ηy[k] * A[i][j].ηy[k]);
					c = sqrt(γ * A[i][j].p / A[i][j].ρ);
					uξ = A[i][j].u * A[i][j].ξx[k] + A[i][j].v * A[i][j].ξy[k];
					uη = A[i][j].u * A[i][j].ηx[k] + A[i][j].v * A[i][j].ηy[k];
					//maxξ = max(Sξ, uξ);
					//maxη = max(Sη, uη);
					maxξ += abs(uξ) + c * Sξ;
					maxη += abs(uξ) + c * Sη;
				}
				maxξ = maxξ / 3;
				maxη = maxη / 3;
				max1 = max(max1, maxξ);
				max2 = max(max2, maxη);
			}
		}
	}
	t = CFL / (max1 + max2);
	//t = CFL / (maxξ + maxη);
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
//			A[i][j].ρ = U[i][j][0];
//			A[i][j].u = U[i][j][1] / U[i][j][0];
//			A[i][j].v = U[i][j][2] / U[i][j][0];
//			A[i][j].p = (γ - 1) * (U[i][j][3] - 0.5 * A[i][j].ρ * (A[i][j].u * A[i][j].u + A[i][j].v * A[i][j].v));
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
			A[i][j].ρ = ρ0;
			A[i][j].u = u0;
			A[i][j].v = v0;
			A[i][j].p = p0;
			//if (A[i].type == "L")
			//{
			//	A[i].ρ = ρ1;
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
				//ρ1 = 3.8571, p1 = 10.3333, u1 = 2.6293, v1 = 0;
				//ρ2 = 1 + 0.2 * sin(5 * A[i][j].x), p2 = 1, u2 = 0, v2 = 0;
				//ρ2 = 3 + 0.5 * sin(5 * A[i][j].x);

				if (A[i][j].type == "IN")
					continue;
				else if (A[i][j].type == "U")
				{
					A[i][j].ρ = A[i][j - Xnum].ρ;
					A[i][j].u = A[i][j - Xnum].u;
					A[i][j].v = A[i][j - Xnum].v;
					A[i][j].p = A[i][j - Xnum].p;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].ρ = (A[i][A[i][j].neibor[0]].ρ + A[i][A[i][j].neibor[1]].ρ) / 2;
						A[i][j].u = (A[i][A[i][j].neibor[0]].u + A[i][A[i][j].neibor[1]].u) / 2;
						A[i][j].v = (A[i][A[i][j].neibor[0]].v + A[i][A[i][j].neibor[1]].v) / 2;
						A[i][j].p = (A[i][A[i][j].neibor[0]].p + A[i][A[i][j].neibor[1]].p) / 2;
					}
				}
				else if (A[i][j].type == "D")
				{
					A[i][j].ρ = A[i][j + Xnum].ρ;
					A[i][j].u = A[i][j + Xnum].u;
					A[i][j].v = A[i][j + Xnum].v;
					A[i][j].p = A[i][j + Xnum].p;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].ρ = (A[i][A[i][j].neibor[0]].ρ + A[i][A[i][j].neibor[1]].ρ) / 2;
						A[i][j].u = (A[i][A[i][j].neibor[0]].u + A[i][A[i][j].neibor[1]].u) / 2;
						A[i][j].v = (A[i][A[i][j].neibor[0]].v + A[i][A[i][j].neibor[1]].v) / 2;
						A[i][j].p = (A[i][A[i][j].neibor[0]].p + A[i][A[i][j].neibor[1]].p) / 2;
					}

				}
				else if (A[i][j].type == "L")
				{
					A[i][j].ρ = ρ1;
					A[i][j].u = u1;
					A[i][j].v = v1;
					A[i][j].p = p1;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].ρ = (A[i][A[i][j].neibor[0]].ρ + A[i][A[i][j].neibor[1]].ρ) / 2;
						A[i][j].u = (A[i][A[i][j].neibor[0]].u + A[i][A[i][j].neibor[1]].u) / 2;
						A[i][j].v = (A[i][A[i][j].neibor[0]].v + A[i][A[i][j].neibor[1]].v) / 2;
						A[i][j].p = (A[i][A[i][j].neibor[0]].p + A[i][A[i][j].neibor[1]].p) / 2;
					}

				}
				else if (A[i][j].type == "R")
				{
					A[i][j].ρ = A[i][j - 1].ρ;
					A[i][j].u = A[i][j - 1].u;
					A[i][j].v = A[i][j - 1].v;
					A[i][j].p = A[i][j - 1].p;
					if (A[i][j].neibor.size() == 2)
					{
						A[i][j].ρ = (A[i][A[i][j].neibor[0]].ρ + A[i][A[i][j].neibor[1]].ρ) / 2;
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
							A[i][j].ρ = A[i][A[i][j].neibor[0]].ρ;
							A[i][j].u = A[i][A[i][j].neibor[0]].u;
							A[i][j].v = A[i][A[i][j].neibor[0]].v;
							A[i][j].p = A[i][A[i][j].neibor[0]].p;
						}
						else
						{
							A[i][j].ρ = A[i][A[i][j].neibor[1]].ρ;
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
		double θ12 = -30 / 180.0 * ConstPara::pi;
		double θ13 = 30 / 180.0 * ConstPara::pi;
		mesh A12, A13, C;
		Line L12, L13;
		A12.x = -dx / 2, A12.y = A0[Pnum - 1].y - 2 * dy * 5 * Ynum / 45;
		A13.x = -dx / 2, A13.y = 2 * dy * 5 * Ynum / 45;
		L12 = getLine(θ12, A12);
		L13 = getLine(θ13, A13);
		C = getCrossPoint(L12, L13);
		L12 = getLine(β2, C);
		L13 = getLine(β3, C);
		mesh A1, A2, A3;
		A1.ρ = ρ1, A1.p = p1, A1.u = u1, A1.v = v1;
		get_down(A1, A2, β2);
		get_down(A1, A3, β3);

		ρ2 = A2.ρ, p2 = A2.p, u2 = A2.u, v2 = A2.v;
		ρ3 = A3.ρ, p3 = A3.p, u3 = A3.u, v3 = A3.v;
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
						A[i][j].ρ = ρ2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;
					}
					else if (A[i][j].x * L13.A + A[i][j].y * L13.B + L13.C > 0)
					{
						A[i][j].ρ = ρ3;
						A[i][j].u = u3;
						A[i][j].v = v3;
						A[i][j].p = p3;
					}
					else
					{
						A[i][j].ρ = ρ1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
				}
				else if (A[i][j].type == "R")
				{
					A[i][j].ρ = A[i][j - 1].ρ;
					A[i][j].u = A[i][j - 1].u;
					A[i][j].v = A[i][j - 1].v;
					A[i][j].p = A[i][j - 1].p;

				}
				else if (A[i][j].type == "U")
				{
					A[i][j].ρ = A[i][j - Xnum].ρ;
					A[i][j].u = A[i][j - Xnum].u;
					A[i][j].v = A[i][j - Xnum].v;
					A[i][j].p = A[i][j - Xnum].p;

				}
				else if (A[i][j].type == "D")
				{
					A[i][j].ρ = A[i][j + Xnum].ρ;
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

		L = getLine(β, M1);

		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "IN")
					continue;
				else if (A[i][j].type == "U")
				{

					A[i][j].ρ = A[i][j - Xnum].ρ;
					A[i][j].u = A[i][j - Xnum].u;
					A[i][j].v = A[i][j - Xnum].v;
					A[i][j].p = A[i][j - Xnum].p;

				}
				else if (A[i][j].type == "D")
				{
					if (A[i][j].x > 0.002)
					{
						A[i][j].ρ = ρ2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;

					}
					else
					{
						A[i][j].ρ = A[i][j + Xnum].ρ;
						A[i][j].u = A[i][j + Xnum].u;
						A[i][j].v = A[i][j + Xnum].v;
						A[i][j].p = A[i][j + Xnum].p;

					}

				}
				else if (A[i][j].type == "L")
				{
					A[i][j].ρ = ρ1;
					A[i][j].u = u1;
					A[i][j].v = v1;
					A[i][j].p = p1;
				}
				else if (A[i][j].type == "R")
				{
					A[i][j].ρ = A[i][j - 1].ρ;
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
							A[i][j].ρ = A[i][A[i][j].neibor[0]].ρ;
							A[i][j].u = A[i][A[i][j].neibor[0]].u;
							A[i][j].v = A[i][A[i][j].neibor[0]].v;
							A[i][j].p = A[i][A[i][j].neibor[0]].p;
						}
						else
						{
							A[i][j].ρ = A[i][A[i][j].neibor[1]].ρ;
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
void update_bound_shockwave_fitting()//激波边界，用于装配法
{
	extern vector<vector <mesh>> A;
	extern vector<vector <mesh>> Ar;
	extern double xL, xR, yU, yD;
	extern int step;
	int i, j, k, m, n;
	double Mas, Ma1 = 0, Ma2;
	double Vs = 0, V1, V1n, V1t, V2n, V2t;
	double c1, c2, p2, ρ2;
	double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
	double θ = 0;
	double θ1, θ2, θ3;
	double R_x, R_y;//参考点值，用于求Vs
	int m1, n1;
	mesh up, down;
	mesh upr, downr;
	double RR, RL, unb, ab;
	double ρ20, p20, u20, v20, c20;
	double ρ21, p21, u21, v21, c21;
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
			θ = get_θ(R_x, R_y, Ar[i][j].x, Ar[i][j].y);//激波角度

			if ((FlowType == "oblique" || FlowType == "normal") && j == 4000)
			{
				for (int k = 0; k < A[i][j].neibor.size(); k++)
				{
					m = A[i][j].neibor[k];
					if (A[i][m].type == "SHOCK")
						x1 = A[i][m].x, y1 = A[i][m].y;
				}
				θ = get_θ(R_x, R_y, x1, y1);//激波角度
				//θ = pi / 2;
			}
			θ = pi / 2;
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

				Ma1 = get_Mu(upr, Ar[i][n1], θ);
				//std::cout << Ma1 << std::endl;
				c1 = get_c(upr.ρ, upr.p);
				c2 = get_c(Ar[i][n1].ρ, Ar[i][n1].p);
				V1n = get_un(upr, θ);
				V1t = get_ut(upr, θ);
				V2n = get_un(Ar[i][n1], θ);
				p2 = upr.p * get_p2p1(Ma1);
				ρ2 = upr.ρ * get_ρ2ρ1(Ma1);

				down.p = p2;
				down.ρ = ρ2;
				RR = (V2n + 2 * c2 / (γ - 1));

				if (θ >= 0)
				{
					Vs = (upr.u * sin(θ) - upr.v * cos(θ)) - Ma1 * c1;
					Vs = -Vs;
					V2n = get_udn(upr, downr, Ma1, Vs, θ);
					V2t = V1t;
				}
				else if (θ < 0)
				{
					Vs = (upr.u * sin(-θ) + upr.v * cos(-θ)) - Ma1 * c1;
					Vs = -Vs;
					V2n = get_udn(upr, downr, Ma1, Vs, θ);
					V2t = V1t;
				}
				c21 = get_c(ρ2, p2);
				RL = (V2n - 2 * c21 / (γ - 1));
				unb = (RR + RL) / 2;
				ab = (γ - 1) * (RR - RL) / 4;
			}
			else if (up.type == "DISCON")
			{
				c1 = get_c(upr.ρ, upr.p);
				c2 = get_c(downr.ρ, downr.p);
				V1n = get_un(upr, θ);
				V2n = get_un(downr, θ);
				V1t = get_ut(upr, θ);
				V2t = get_ut(downr, θ);
				double J1 = V1n + 2 * c1 / (γ - 1);
				double J2 = V2n - 2 * c2 / (γ - 1);
				double S1 = upr.p / pow(upr.ρ, γ);
				double S2 = downr.p / pow(downr.ρ, γ);
				double pb, pb1, pb0;
				pb1 = upr.p;
				pb0 = downr.p;
				pb = (upr.p + downr.p) / 2;
				while (abs(pb - pb1) > 1e-10)
				{
					pb0 = pb1;
					pb1 = pb;
					double fxn_1 = 2 * sqrt(γ * pow(pb0, (γ - 1) / γ) * pow(S1, 1 / γ)) / (γ - 1) + 2 * sqrt(γ * pow(pb0, (γ - 1) / γ) * pow(S2, 1 / γ)) / (γ - 1) - (J2 - J1);
					double fxn = 2 * sqrt(γ * pow(pb1, (γ - 1) / γ) * pow(S1, 1 / γ)) / (γ - 1) + 2 * sqrt(γ * pow(pb1, (γ - 1) / γ) * pow(S2, 1 / γ)) / (γ - 1) - (J2 - J1);
					pb = pb1 - fxn / (fxn - fxn_1) * (pb1 - pb0);
				}
				//上游参数更新
				up.p = pb;
				up.ρ = pow((up.p / S1), 1 / γ);
				c1 = get_c(up.ρ, up.p);
				V1n = 2 * c1 / (γ - 1) - J1;
				//下游参数更新
				down.p = pb;
				down.ρ = pow((down.p / S2), 1 / γ);
				c2 = get_c(down.ρ, down.p);
				V2n = 2 * c2 / (γ - 1) + J2;
				Vs = (V1n + V2n) / 2;
			}
			if (θ >= 0)
			{
				up.um = down.um = Vs * sin(θ);
				up.vm = down.vm = -Vs * cos(θ);
				//up.u = V1t * cos(θ) + V1n * sin(θ);
				//up.v = V1t * sin(θ) - V1n * cos(θ);

				//down.u = V2t * cos(θ) + V2n * sin(θ);
				//down.v = V2t * sin(θ) - V2n * cos(θ);
				//A[i][n1].u = unb;
				//A[i][n1].ρ = ρ2 * pow((ab * ab / (c21 * c21)), 1 / (γ - 1));
				//A[i][n1].p = ab * ab * down.ρ / γ;
				double p20 = down.p;
				double ρ20 = down.ρ;
				down.p = p2;
				down.ρ = ρ2 * pow(down.p / p20 ,1 / γ);
				double c = V2n * V2n / 2 + γ * down.p / ((γ - 1) * down.ρ);
				down.u = sqrt(2 * (c - γ * down.p / ((γ - 1) * down.ρ)));

				//down.p = up.p;
				//down.ρ = up.ρ;
				//down.u = up.u;
				//down.v = up.v;
			}
			else
			{
				up.um = down.um = Vs * sin(-θ);
				up.vm = down.vm = Vs * cos(-θ);
				//up.u = V1t * cos(-θ) + V1n * sin(-θ);
				//up.v = -V1t * sin(-θ) + V1n * cos(-θ);
				down.u = V2t * cos(-θ) + V2n * sin(-θ);
				down.v = -V2t * sin(-θ) + V2n * cos(-θ);

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
double get_θ(double x1, double y1, double x2, double y2)//求直线与x轴的夹角
{
	double θ;
	if (x1 == x2)
	{
		//if (y2 < y1)
		//	θ = -pi / 2;
		//else
		//	θ = pi / 2;
		θ = pi / 2;
	}
	else if (y1 == y2)
		θ = 0;
	else
	{
		θ = atan(abs((y2 - y1) / (x2 - x1)));
		if ((x2 > x1 && y2 < y1) || (x2 < x1 && y2 > y1))
			θ = -θ;
	}
	return θ;
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
void clear_Vm()//清除各点的运动速度
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
						A[i][j].ρ = ρ1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
					if (i == 1)
					{
						A[i][j].ρ = ρ2;
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
						A[i][j].ρ = ρ1;
						A[i][j].u = u1;
						A[i][j].v = v1;
						A[i][j].p = p1;
					}
					if (i == 1)
					{
						A[i][j].ρ = ρ2;
						A[i][j].u = u2;
						A[i][j].v = v2;
						A[i][j].p = p2;
					}
					if (A[i][j].type == "R")
					{
						A[i][j].ρ = A[i][j - 1].ρ;
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
					//ρ1 = 3.8571, p1 = 10.3333, u1 = 2.6293, v1 = 0;
					//ρ2 = 0.5 + 0.1 * sin(5 * A[i][j].x), p2 = 1, u2 = 0, v2 = 0;
					ρ2 = 3 +1*sin(5 * A[i][j].x);
					//ρ2 = 1 + 0.1 * sin(5 * A[i][j].x);
					if (i == 1)
					{
						A[i][j].ρ = ρ2;
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
						A[i][j].ρ = A[i][n].ρ;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
					else if (A[i][j].type == "L")
					{
						A[i][j].ρ = ρ1;
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
						A[i][j].ρ = A[i][n].ρ;
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
						A[i][j].ρ = A[i][n].ρ;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
					else if (A[i][j].type == "D")
					{
						if (i == 0)
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
					else if (A[i][j].type == "L")
					{
						A[i][j].ρ = ρ1;
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
						A[i][j].ρ = A[i][n].ρ;
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
						A[i][j].ρ = A[i][n].ρ;
						A[i][j].u = A[i][n].u;
						A[i][j].v = A[i][n].v;
						A[i][j].p = A[i][n].p;
					}
					else if (A[i][j].type == "L")
					{
						if (i == 0)
						{
							A[i][j].ρ = ρ1;
							A[i][j].u = u1;
							A[i][j].v = v1;
							A[i][j].p = p1;
						}
						else if (i == 1)
						{
							A[i][j].ρ = ρ2;
							A[i][j].u = u2;
							A[i][j].v = v2;
							A[i][j].p = p2;
						}
						else if (i == 2)
						{
							A[i][j].ρ = ρ3;
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
						A[i][j].ρ = A[i][n].ρ;
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
				A[i][j].ρ = ρ2;
				A[i][j].u = u2;
				A[i][j].v = v2;
				A[i][j].p = p2;
			}
			else if (A[i][j].type == "L")
			{
				A[i][j].ρ = ρ1;
				A[i][j].u = u1;
				A[i][j].v = v1;
				A[i][j].p = p1;

			}
			else if (A[i][j].type == "D")
			{
				A[i][j].ρ = ρ3;
				A[i][j].u = u3;
				A[i][j].v = v3;
				A[i][j].p = p3;
			}
			else if (A[i][j].type == "R")
			{
				A[i][j].ρ = A[i][j - 1].ρ;
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
				A[i][j].ρ = ρ1;
				A[i][j].u = u1;
				A[i][j].v = v1;
				A[i][j].p = p1;
			}
			else if (A[i][j].type == "DL")
			{
				A[i][j].ρ = ρ1;
				A[i][j].u = u1;
				A[i][j].v = v1;
				A[i][j].p = p1;
			}
			else if (A[i][j].type == "DR")
			{
				A[i][j].ρ = ρ2;
				A[i][j].u = u2;
				A[i][j].v = v2;
				A[i][j].p = p2;
			}
			else if (A[i][j].type == "R")
			{
				A[i][j].ρ = A[i][j - 1].ρ;
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
//		F[i][0] = A[i].ρ*A[i].u;
//		F[i][1] = A[i].ρ*A[i].u*A[i].u + A[i].p;
//		F[i][2] = A[i].ρ*A[i].u*A[i].v;
//		F[i][3] = A[i].u*(A[i].p / (γ - 1) + 0.5*A[i].ρ*(A[i].u*A[i].u + A[i].v*A[i].v) + A[i].p);
//	}
//}
//void get_G()
//{
//	extern vector <mesh> A;
//	extern double G[Pnum][4];
//	for (int i = 0; i < Pnum; i++)
//	{
//		G[i][0] = A[i].ρ*A[i].v; extern vector<vector<vector <double>>> Utr;
//		G[i][1] = A[i].ρ*A[i].u *A[i].v;
//		G[i][2] = A[i].ρ*A[i].u*A[i].v + A[i].p;
//		G[i][3] = A[i].v*(A[i].p / (γ - 1) + 0.5*A[i].ρ*(A[i].u*A[i].u + A[i].v*A[i].v) + A[i].p);
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
double compute_res()//计算残差
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
			res += abs(A[i][j].ρ - Ar[i][j].ρ) / Ar[i][j].ρ;
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
void judge()//判断周围物理量是否改变，如果没改变的话就跳过，计算加速？？
{
	extern vector <vector<mesh>> A;
	extern vector<vector <mesh>> Ar;
	int judge;
	extern int step;
	if (step == 0)
		std::cout << "计算开始！" << std::endl;
	else
		for (int i = 0; i < A.size(); i++)
		{
			for (int j = 0; j < A[i].size(); j++)
			{
				judge = 0;
				for (int k = 0; k < A[i][j].neibor.size(); k++)
				{
					int n = A[i][j].neibor[k];
					if (A[i][n].p != Ar[i][n].p && A[i][n].u != Ar[i][n].u && A[i][n].v != Ar[i][n].v && A[i][n].ρ != Ar[i][n].ρ)
						judge = 1;
				}
				if (judge == 0)
					A[i][j].change = "N";
				else
					A[i][j].change = "Y";
			}
		}
}

void update_IN()//更新内部点
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
			U[0] = Ar[i][j].ρ;
			U[1] = Ar[i][j].ρ * Ar[i][j].u;
			U[2] = Ar[i][j].ρ * Ar[i][j].v;
			U[3] = 0.5 * Ar[i][j].ρ * (Ar[i][j].u * Ar[i][j].u + Ar[i][j].v * Ar[i][j].v) + Ar[i][j].p / (γ - 1);

			if (Ar[i][j].neibor.size() == 3)
			{
				F1 = F2 = G1 = G2 = { 0 };

				for (k = 0; k < 12; k++)
				{
					n1 = method[k][0];
					n2 = method[k][1];
					n3 = method[k][2];
					n4 = method[k][3];

					F1_ = HLLC_Χ2(Ar[i][Ar[i][j].neibor[n3]], Ar[i][j], Ar[i][j], k);
					F2_ = HLLC_Χ2(Ar[i][j], Ar[i][Ar[i][j].neibor[n1]], Ar[i][j], k);
					G1_ = HLLC_Υ2(Ar[i][Ar[i][j].neibor[n4]], Ar[i][j], Ar[i][j], k);
					G2_ = HLLC_Υ2(Ar[i][j], Ar[i][Ar[i][j].neibor[n2]], Ar[i][j], k);
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
				//	F1 = HLLC_Χ(Ar[i][Ar[i][j].neibor[2]], Ar[i][j], Ar[i][j], k);
				//	F2 = HLLC_Χ(Ar[i][j], Ar[i][Ar[i][j].neibor[0]], Ar[i][j], k);
				//	G1 = HLLC_Υ(Ar[i][Ar[i][j].neibor[2]], Ar[i][j], Ar[i][j], k);
				//	G2 = HLLC_Υ(Ar[i][j], Ar[i][Ar[i][j].neibor[1]], Ar[i][j], k);
				//	d1 = distance(Ar[i][Ar[i][j].neibor[2]], Ar[i][Ar[i][j].neibor[0]]) / 2;
				//	d2 = distance(Ar[i][Ar[i][j].neibor[2]], Ar[i][Ar[i][j].neibor[1]]) / 2;
				//}
				//else if (k == 1)
				//{
				//	F1 = HLLC_Χ(Ar[i][Ar[i][j].neibor[0]], Ar[i][j], Ar[i][j], k);
				//	F2 = HLLC_Χ(Ar[i][j], Ar[i][Ar[i][j].neibor[1]], Ar[i][j], k);
				//	G1 = HLLC_Υ(Ar[i][Ar[i][j].neibor[0]], Ar[i][j], Ar[i][j], k);
				//	G2 = HLLC_Υ(Ar[i][j], Ar[i][Ar[i][j].neibor[2]], Ar[i][j], k);
				//	d1 = distance(Ar[i][Ar[i][j].neibor[0]], Ar[i][Ar[i][j].neibor[1]]) / 2;
				//	d2 = distance(Ar[i][Ar[i][j].neibor[0]], Ar[i][Ar[i][j].neibor[2]]) / 2;
				//}
				//else
				//{
				//	F1 = HLLC_Χ(Ar[i][Ar[i][j].neibor[1]], Ar[i][j], Ar[i][j], k);
				//	F2 = HLLC_Χ(Ar[i][j], Ar[i][Ar[i][j].neibor[2]], Ar[i][j], k);
				//	G1 = HLLC_Υ(Ar[i][Ar[i][j].neibor[1]], Ar[i][j], Ar[i][j], k);
				//	G2 = HLLC_Υ(Ar[i][j], Ar[i][Ar[i][j].neibor[0]], Ar[i][j], k);
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
				F1 = HLLC_Χ2(Ar[i][Ar[i][j].neibor[2]], Ar[i][j], Ar[i][j], 0);
				F2 = HLLC_Χ2(Ar[i][j], Ar[i][Ar[i][j].neibor[0]], Ar[i][j], 0);
				G1 = HLLC_Υ2(Ar[i][Ar[i][j].neibor[3]], Ar[i][j], Ar[i][j], 0);
				G2 = HLLC_Υ2(Ar[i][j], Ar[i][Ar[i][j].neibor[1]], Ar[i][j], 0);

				U[0] = U[0] - dt * (F2.f1 - F1.f1) - dt * (G2.f1 - G1.f1);
				U[1] = U[1] - dt * (F2.f2 - F1.f2) - dt * (G2.f2 - G1.f2);
				U[2] = U[2] - dt * (F2.f3 - F1.f3) - dt * (G2.f3 - G1.f3);
				U[3] = U[3] - dt * (F2.f4 - F1.f4) - dt * (G2.f4 - G1.f4);
			}
			else if (A[i][j].neibor.size() > 4)
			{
				F1 = HLLC_Χ(A[i][A[i][j].neibor1[2]], A[i][j], A[i][j], 0);
				F2 = HLLC_Χ(A[i][j], A[i][A[i][j].neibor1[0]], A[i][j], 0);
				G1 = HLLC_Υ(A[i][A[i][j].neibor1[3]], A[i][j], A[i][j], 0);
				G2 = HLLC_Υ(A[i][j], A[i][A[i][j].neibor1[1]], A[i][j], 0);
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

					A[i][j].ρ = (A[i][m].ρ + A[i][n].ρ) / 2;
					A[i][j].u = (A[i][m].u + A[i][n].u) / 2;
					A[i][j].v = (A[i][m].v + A[i][n].v) / 2;
					A[i][j].p = (A[i][m].p + A[i][n].p) / 2;
				}
				else
					continue;
			}
			else
				std::cout << i << "   " << j << "   " << "somthing wrong in updat_u" << std::endl;
			A[i][j].ρ = U[0];
			A[i][j].u = U[1] / U[0];
			A[i][j].v = U[2] / U[0];
			A[i][j].p = (γ - 1) * (U[3] - 0.5 * A[i][j].ρ * (A[i][j].u * A[i][j].u + A[i][j].v * A[i][j].v));
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
//	double ρ00, u00, v00, p00;
//	double ρ, u, v, p;
//	double ρ0, u0, v0, p0;
//	double ρ1, u1, v1, p1;
//	double ρ2, u2, v2, p2;
//	ρ00 = A[i].ρ;
//	u00 = A[i].u;
//	v00 = A[i].v;
//	p00 = A[i].p;
//
//	ρ0 = U0[i][0] / A[i].J[0];
//	u0 = U0[i][1] / U0[i][0];
//	v0 = U0[i][2] / U0[i][0];
//	p0 = (γ - 1)*(U0[i][3] / A[i].J[0] - 0.5*A[i].ρ*(A[i].u*A[i].u + A[i].v*A[i].v));
//
//	ρ1 = U1[i][0] / A[i].J[1];
//	u1 = U1[i][1] / U1[i][0];
//	v1 = U1[i][2] / U1[i][0];
//	p1 = (γ - 1)*(U1[i][3] / A[i].J[1] - 0.5*A[i].ρ*(A[i].u*A[i].u + A[i].v*A[i].v));
//	//ρ0 = ρ1;
//	//u0 = u1;
//	//v0 = v1;
//	//p0 = p1;
//
//	ρ2 = U2[i][0] / A[i].J[2];
//	u2 = U2[i][1] / U2[i][0];
//	v2 = U2[i][2] / U2[i][0];
//	p2 = (γ - 1)*(U2[i][3] / A[i].J[1] - 0.5*A[i].ρ*(A[i].u*A[i].u + A[i].v*A[i].v));
//	ρ = absmax(absmax(ρ0, ρ1), ρ2);
//	u = absmin(absmin(u0, u1), u2);
//	v = absmax(absmax(v0, v1), v2);
//	p = absmax(absmax(p0, p1), p2);
//	U[i][0] = ρ;
//	U[i][1] = ρ * u;
//	U[i][2] = ρ * v;
//	U[i][3] = 0.5*ρ*(u*u + v * v) + p / (γ - 1);
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
double get_β(mesh A, mesh B)//求出两个网格点与x轴的夹角
{
	double dy = abs(A.y - B.y);
	double dx = abs(A.x - B.x);
	double β = atan(dy / dx);
	return β;
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
double area(mesh A, mesh B, mesh C, mesh D)//求任意四点构成四边形面积 
{
	return 0.5* abs(A.x * B.y + B.x * C.y + C.x * D.y + D.x * A.y - B.x * A.y - C.x * B.y - D.x * C.y - A.x * D.y);
}

void movemesh()//网格移动
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
Line getLine(double θ, mesh A)
{
	double k = tan(θ);
	double b = A.y - k * A.x;
	Line L;
	L.A = k, L.B = -1, L.C = b;
	return L;

}