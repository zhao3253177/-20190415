//分区文件，用于激波装配法，将原网格分成若干区域
#include<iostream>
#include"const.h"
#include"functions.h"
#include<vector>
#include<string>
#include<fstream>
using std::vector;
extern vector <mesh> A0;
extern vector <polygon_mesh> M0;
using namespace ShockwaveCross;
using MeshPara::Pnum;

extern vector<vector <mesh>> A;
extern vector<vector <polygon_mesh>> M;

Line L12;
Line L13;
mesh C;
int findAd(mesh A, vector <mesh> An)
{
	int s = -1;
	for (int i = 0; i < An.size(); i++)
	{
		if (A.x == An[i].x && A.y == An[i].y)
			s = i;
		else
			continue;
	}
	return s;
}

void partition_Point()//对已有的网格点进行分区
{
	using namespace MeshPara;

	if (methodType == "C")
	{
		vector<mesh> t;
		int i;
		A.push_back(t);
		for (i = 0; i < A0.size(); i++)
			A[0].push_back(A0[i]);
	}
	else if (methodType == "F")
	{
		if (FlowType == "normal" || FlowType == "oblique")
		{
			vector<mesh> t;
			int i;
			A.push_back(t);
			A.push_back(t);
			for (i = 0; i < A0.size(); i++)
			{
				if (A0[i].x < 1.003)
					A0[i].section = 1, A[0].push_back(A0[i]);
				else
					A0[i].section = 2, A[1].push_back(A0[i]);
			}
		}
		else if (FlowType == "intersection")
		{
			double θ12 = -30 / 180.0 * ConstPara::pi;
			double θ13 = 30 / 180.0 * ConstPara::pi;
			mesh A12;
			mesh A13;
			vector<mesh> t;
			A12.x = -dx / 2, A12.y = A0[Pnum - 1].y - 2 * dy * 5 * Ynum / 43;
			A13.x = -dx / 2, A13.y = 2 * dy * 5 * Ynum / 43;
			L12 = getLine(θ12, A12);
			L13 = getLine(θ13, A13);
			C = getCrossPoint(L12, L13);
			C.type = "CENTER";
			A.push_back(t);
			A.push_back(t);
			A.push_back(t);
			A.push_back(t);
			A.push_back(t);
			A[0].push_back(C);
			A[1].push_back(C);
			A[2].push_back(C);
			A[3].push_back(C);
			A[4].push_back(C);
			int i;

			for (i = 0; i < A0.size(); i++)
			{
				if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C > 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C < 0)
					A0[i].section = 1, A[0].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C < 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C < 0)
					A0[i].section = 2, A[1].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C > 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C > 0)
					A0[i].section = 3, A[2].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C < 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C > 0 && A0[i].y > C.y && abs(A0[i].y - C.y) > 1e-10)
					A0[i].section = 4, A[3].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C < 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C > 0 && A0[i].y < C.y && abs(A0[i].y - C.y) > 1e-10)
					A0[i].section = 5, A[4].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C == 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C < 0)
					A0[i].section = 12, A[0].push_back(A0[i]), A[1].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C == 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C > 0)
					A0[i].section = 35, A[2].push_back(A0[i]), A[4].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C > 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C == 0)
					A0[i].section = 13, A[0].push_back(A0[i]), A[2].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C < 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C == 0)
					A0[i].section = 24, A[1].push_back(A0[i]), A[3].push_back(A0[i]);
				else if (A0[i].x * L12.A + A0[i].y * L12.B + L12.C < 0 && A0[i].x * L13.A + A0[i].y * L13.B + L13.C > 0 && abs(A0[i].y - C.y) <= 1e-10)
					A0[i].section = 45, A0[i].type = "DISCON", A[3].push_back(A0[i]), A[4].push_back(A0[i]);
				else
					std::cout << "something wrong in partition !" << std::endl;
			}
		}

	}


}
void partition_Mesh()//对多边形网格进行划分以及增加点
{

	int i, j, s;
	int  n0, n1, n2, n3, n4, n5;
	polygon_mesh m;
	vector <polygon_mesh> t;
	int  ad0, ad1, ad2, ad3, ad4, ad5;
	using namespace MeshPara;
	if (methodType == "C")
	{
		M.push_back(t);
		for (i = 0; i < M0.size(); i++)
		{
			M[0].push_back(M0[i]);
		}
	}
	else if (methodType == "F")
	{
		if (FlowType == "normal" || FlowType == "oblique")
		{
			M.push_back(t);
			M.push_back(t);
			for (i = 0; i < M0.size(); i++)
			{
				n0 = M0[i].node[0];
				n1 = M0[i].node[1];
				n2 = M0[i].node[2];
				n3 = M0[i].node[3];
				n4 = M0[i].node[4];
				n5 = M0[i].node[5];
				if (A0[n0].section == 1 && A0[n1].section == 1 && A0[n2].section == 1 && A0[n3].section == 1 && A0[n4].section == 1 && A0[n5].section == 1)
				{
					ad0 = findAd(A0[n0], A[0]);
					ad1 = findAd(A0[n1], A[0]);
					ad2 = findAd(A0[n2], A[0]);
					ad3 = findAd(A0[n3], A[0]);
					ad4 = findAd(A0[n4], A[0]);
					ad5 = findAd(A0[n5], A[0]);
					M[0].push_back(m);
					s = int(M[0].size()) - 1;
					M[0][s].node.push_back(ad0);
					M[0][s].node.push_back(ad1);
					M[0][s].node.push_back(ad2);
					M[0][s].node.push_back(ad3);
					M[0][s].node.push_back(ad4);
					M[0][s].node.push_back(ad5);

					M[0][s].face_start.push_back(ad0);
					M[0][s].face_start.push_back(ad1);
					M[0][s].face_start.push_back(ad2);
					M[0][s].face_start.push_back(ad3);
					M[0][s].face_start.push_back(ad4);
					M[0][s].face_start.push_back(ad5);

					M[0][s].face_end.push_back(ad1);
					M[0][s].face_end.push_back(ad2);
					M[0][s].face_end.push_back(ad3);
					M[0][s].face_end.push_back(ad4);
					M[0][s].face_end.push_back(ad5);
					M[0][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 2 && A0[n1].section == 2 && A0[n2].section == 2 && A0[n3].section == 2 && A0[n4].section == 2 && A0[n5].section == 2)
				{
					ad0 = findAd(A0[n0], A[1]);
					ad1 = findAd(A0[n1], A[1]);
					ad2 = findAd(A0[n2], A[1]);
					ad3 = findAd(A0[n3], A[1]);
					ad4 = findAd(A0[n4], A[1]);
					ad5 = findAd(A0[n5], A[1]);
					M[1].push_back(m);
					s = int(M[1].size()) - 1;
					M[1][s].node.push_back(ad0);
					M[1][s].node.push_back(ad1);
					M[1][s].node.push_back(ad2);
					M[1][s].node.push_back(ad3);
					M[1][s].node.push_back(ad4);
					M[1][s].node.push_back(ad5);

					M[1][s].face_start.push_back(ad0);
					M[1][s].face_start.push_back(ad1);
					M[1][s].face_start.push_back(ad2);
					M[1][s].face_start.push_back(ad3);
					M[1][s].face_start.push_back(ad4);
					M[1][s].face_start.push_back(ad5);

					M[1][s].face_end.push_back(ad1);
					M[1][s].face_end.push_back(ad2);
					M[1][s].face_end.push_back(ad3);
					M[1][s].face_end.push_back(ad4);
					M[1][s].face_end.push_back(ad5);
					M[1][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 1 && A0[n1].section == 2 && A0[n2].section == 2 && A0[n3].section == 2 && A0[n4].section == 1 && A0[n5].section == 1)
				{
					mesh m1, m2;
					m1.x = (A0[n0].x + A0[n1].x) / 2, m1.y = (A0[n0].y + A0[n1].y) / 2;
					m2.x = (A0[n3].x + A0[n4].x) / 2, m2.y = (A0[n3].y + A0[n4].y) / 2;
					//m1.x = m2.x = 0.0071;
					ad0 = findAd(A0[n0], A[0]);
					ad1 = findAd(m1, A[0]);
					if (findAd(m1, A[0]) == -1)
					{
						A[0].push_back(m1);
						ad1 = A[0].size() - 1;
						A[0][ad1].type = "SHOCK";
					}
					ad2 = findAd(m2, A[0]);
					if (findAd(m2, A[0]) == -1)
					{
						A[0].push_back(m2);
						ad2 = A[0].size() - 1;
						A[0][ad2].type = "SHOCK";

					}
					ad3 = findAd(A0[n4], A[0]);
					ad4 = findAd(A0[n5], A[0]);
					M[0].push_back(m);
					s = int(M[0].size()) - 1;
					M[0][s].node.push_back(ad0);
					M[0][s].node.push_back(ad1);
					M[0][s].node.push_back(ad2);
					M[0][s].node.push_back(ad3);
					M[0][s].node.push_back(ad4);

					M[0][s].face_start.push_back(ad0);
					M[0][s].face_start.push_back(ad1);
					M[0][s].face_start.push_back(ad2);
					M[0][s].face_start.push_back(ad3);
					M[0][s].face_start.push_back(ad4);

					M[0][s].face_end.push_back(ad1);
					M[0][s].face_end.push_back(ad2);
					M[0][s].face_end.push_back(ad3);
					M[0][s].face_end.push_back(ad4);
					M[0][s].face_end.push_back(ad0);
					ad0 = findAd(m1, A[1]);
					if (findAd(m1, A[1]) == -1)
					{
						A[1].push_back(m1);
						ad0 = A[1].size() - 1;
						A[1][ad0].type = "SHOCK";

					}

					ad1 = findAd(A0[n1], A[1]);
					ad2 = findAd(A0[n2], A[1]);
					ad3 = findAd(A0[n3], A[1]);
					ad4 = findAd(m2, A[1]);
					if (findAd(m2, A[1]) == -1)
					{
						A[1].push_back(m2);
						ad4 = A[1].size() - 1;
						A[1][ad4].type = "SHOCK";

					}
					M[1].push_back(m);
					s = int(M[1].size()) - 1;
					M[1][s].node.push_back(ad0);
					M[1][s].node.push_back(ad1);
					M[1][s].node.push_back(ad2);
					M[1][s].node.push_back(ad3);
					M[1][s].node.push_back(ad4);

					M[1][s].face_start.push_back(ad0);
					M[1][s].face_start.push_back(ad1);
					M[1][s].face_start.push_back(ad2);
					M[1][s].face_start.push_back(ad3);
					M[1][s].face_start.push_back(ad4);

					M[1][s].face_end.push_back(ad1);
					M[1][s].face_end.push_back(ad2);
					M[1][s].face_end.push_back(ad3);
					M[1][s].face_end.push_back(ad4);
					M[1][s].face_end.push_back(ad0);

				}

			}

		}
		else if (FlowType == "intersection")
		{
			Line L1, L2;
			M.push_back(t);
			M.push_back(t);
			M.push_back(t);
			M.push_back(t);
			M.push_back(t);
			for (i = 0; i < M0.size(); i++)
			{
				n0 = M0[i].node[0];
				n1 = M0[i].node[1];
				n2 = M0[i].node[2];
				n3 = M0[i].node[3];
				n4 = M0[i].node[4];
				n5 = M0[i].node[5];

				if (A0[n0].section == 1 && A0[n1].section == 1 && A0[n2].section == 1 && A0[n3].section == 1 && A0[n4].section == 1 && A0[n5].section == 1)
				{
					ad0 = findAd(A0[n0], A[0]);
					ad1 = findAd(A0[n1], A[0]);
					ad2 = findAd(A0[n2], A[0]);
					ad3 = findAd(A0[n3], A[0]);
					ad4 = findAd(A0[n4], A[0]);
					ad5 = findAd(A0[n5], A[0]);
					M[0].push_back(m);
					s = int(M[0].size()) - 1;
					M[0][s].node.push_back(ad0);
					M[0][s].node.push_back(ad1);
					M[0][s].node.push_back(ad2);
					M[0][s].node.push_back(ad3);
					M[0][s].node.push_back(ad4);
					M[0][s].node.push_back(ad5);

					M[0][s].face_start.push_back(ad0);
					M[0][s].face_start.push_back(ad1);
					M[0][s].face_start.push_back(ad2);
					M[0][s].face_start.push_back(ad3);
					M[0][s].face_start.push_back(ad4);
					M[0][s].face_start.push_back(ad5);

					M[0][s].face_end.push_back(ad1);
					M[0][s].face_end.push_back(ad2);
					M[0][s].face_end.push_back(ad3);
					M[0][s].face_end.push_back(ad4);
					M[0][s].face_end.push_back(ad5);
					M[0][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 2 && A0[n1].section == 2 && A0[n2].section == 2 && A0[n3].section == 2 && A0[n4].section == 2 && A0[n5].section == 2)
				{
					ad0 = findAd(A0[n0], A[1]);
					ad1 = findAd(A0[n1], A[1]);
					ad2 = findAd(A0[n2], A[1]);
					ad3 = findAd(A0[n3], A[1]);
					ad4 = findAd(A0[n4], A[1]);
					ad5 = findAd(A0[n5], A[1]);
					M[1].push_back(m);
					s = int(M[1].size()) - 1;
					M[1][s].node.push_back(ad0);
					M[1][s].node.push_back(ad1);
					M[1][s].node.push_back(ad2);
					M[1][s].node.push_back(ad3);
					M[1][s].node.push_back(ad4);
					M[1][s].node.push_back(ad5);

					M[1][s].face_start.push_back(ad0);
					M[1][s].face_start.push_back(ad1);
					M[1][s].face_start.push_back(ad2);
					M[1][s].face_start.push_back(ad3);
					M[1][s].face_start.push_back(ad4);
					M[1][s].face_start.push_back(ad5);

					M[1][s].face_end.push_back(ad1);
					M[1][s].face_end.push_back(ad2);
					M[1][s].face_end.push_back(ad3);
					M[1][s].face_end.push_back(ad4);
					M[1][s].face_end.push_back(ad5);
					M[1][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 3 && A0[n1].section == 3 && A0[n2].section == 3 && A0[n3].section == 3 && A0[n4].section == 3 && A0[n5].section == 3)
				{
					ad0 = findAd(A0[n0], A[2]);
					ad1 = findAd(A0[n1], A[2]);
					ad2 = findAd(A0[n2], A[2]);
					ad3 = findAd(A0[n3], A[2]);
					ad4 = findAd(A0[n4], A[2]);
					ad5 = findAd(A0[n5], A[2]);
					M[2].push_back(m);
					s = int(M[2].size()) - 1;
					M[2][s].node.push_back(ad0);
					M[2][s].node.push_back(ad1);
					M[2][s].node.push_back(ad2);
					M[2][s].node.push_back(ad3);
					M[2][s].node.push_back(ad4);
					M[2][s].node.push_back(ad5);

					M[2][s].face_start.push_back(ad0);
					M[2][s].face_start.push_back(ad1);
					M[2][s].face_start.push_back(ad2);
					M[2][s].face_start.push_back(ad3);
					M[2][s].face_start.push_back(ad4);
					M[2][s].face_start.push_back(ad5);

					M[2][s].face_end.push_back(ad1);
					M[2][s].face_end.push_back(ad2);
					M[2][s].face_end.push_back(ad3);
					M[2][s].face_end.push_back(ad4);
					M[2][s].face_end.push_back(ad5);
					M[2][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 4 && A0[n1].section == 4 && A0[n2].section == 4 && A0[n3].section == 4 && A0[n4].section == 4 && A0[n5].section == 4)
				{
					ad0 = findAd(A0[n0], A[3]);
					ad1 = findAd(A0[n1], A[3]);
					ad2 = findAd(A0[n2], A[3]);
					ad3 = findAd(A0[n3], A[3]);
					ad4 = findAd(A0[n4], A[3]);
					ad5 = findAd(A0[n5], A[3]);
					M[3].push_back(m);
					s = int(M[3].size()) - 1;
					M[3][s].node.push_back(ad0);
					M[3][s].node.push_back(ad1);
					M[3][s].node.push_back(ad2);
					M[3][s].node.push_back(ad3);
					M[3][s].node.push_back(ad4);
					M[3][s].node.push_back(ad5);

					M[3][s].face_start.push_back(ad0);
					M[3][s].face_start.push_back(ad1);
					M[3][s].face_start.push_back(ad2);
					M[3][s].face_start.push_back(ad3);
					M[3][s].face_start.push_back(ad4);
					M[3][s].face_start.push_back(ad5);

					M[3][s].face_end.push_back(ad1);
					M[3][s].face_end.push_back(ad2);
					M[3][s].face_end.push_back(ad3);
					M[3][s].face_end.push_back(ad4);
					M[3][s].face_end.push_back(ad5);
					M[3][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 5 && A0[n1].section == 5 && A0[n2].section == 5 && A0[n3].section == 5 && A0[n4].section == 5 && A0[n5].section == 5)
				{
					ad0 = findAd(A0[n0], A[4]);
					ad1 = findAd(A0[n1], A[4]);
					ad2 = findAd(A0[n2], A[4]);
					ad3 = findAd(A0[n3], A[4]);
					ad4 = findAd(A0[n4], A[4]);
					ad5 = findAd(A0[n5], A[4]);
					M[4].push_back(m);
					s = int(M[4].size()) - 1;
					M[4][s].node.push_back(ad0);
					M[4][s].node.push_back(ad1);
					M[4][s].node.push_back(ad2);
					M[4][s].node.push_back(ad3);
					M[4][s].node.push_back(ad4);
					M[4][s].node.push_back(ad5);

					M[4][s].face_start.push_back(ad0);
					M[4][s].face_start.push_back(ad1);
					M[4][s].face_start.push_back(ad2);
					M[4][s].face_start.push_back(ad3);
					M[4][s].face_start.push_back(ad4);
					M[4][s].face_start.push_back(ad5);

					M[4][s].face_end.push_back(ad1);
					M[4][s].face_end.push_back(ad2);
					M[4][s].face_end.push_back(ad3);
					M[4][s].face_end.push_back(ad4);
					M[4][s].face_end.push_back(ad5);
					M[4][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 1 && A0[n1].section == 1 && A0[n2].section == 2 && A0[n3].section == 2 && A0[n4].section == 2 && A0[n5].section == 1)
				{
					ad0 = findAd(A0[n0], A[0]);
					ad1 = findAd(A0[n1], A[0]);

					L1 = getLine(A0[n1], A0[n2]);
					L2 = getLine(A0[n4], A0[n5]);
					if (findAd(getCrossPoint(L1, L12), A[0]) == -1)
					{
						A[0].push_back(getCrossPoint(L1, L12));
						ad2 = A[0].size() - 1;
						A[0][ad2].type = "SHOCK";

					}
					else
						ad2 = findAd(getCrossPoint(L1, L12), A[0]);

					if (findAd(getCrossPoint(L2, L12), A[0]) == -1)
					{
						A[0].push_back(getCrossPoint(L2, L12));
						ad3 = A[0].size() - 1;
						A[0][ad3].type = "SHOCK";

					}
					else
						ad3 = findAd(getCrossPoint(L2, L12), A[0]);
					ad4 = findAd(A0[n5], A[0]);
					M[0].push_back(m);
					s = int(M[0].size()) - 1;
					M[0][s].node.push_back(ad0);
					M[0][s].node.push_back(ad1);
					M[0][s].node.push_back(ad2);
					M[0][s].node.push_back(ad3);
					M[0][s].node.push_back(ad4);

					M[0][s].face_start.push_back(ad0);
					M[0][s].face_start.push_back(ad1);
					M[0][s].face_start.push_back(ad2);
					M[0][s].face_start.push_back(ad3);
					M[0][s].face_start.push_back(ad4);

					M[0][s].face_end.push_back(ad1);
					M[0][s].face_end.push_back(ad2);
					M[0][s].face_end.push_back(ad3);
					M[0][s].face_end.push_back(ad4);
					M[0][s].face_end.push_back(ad0);


					if (findAd(getCrossPoint(L1, L12), A[1]) == -1)
					{
						A[1].push_back(getCrossPoint(L1, L12));
						ad0 = A[1].size() - 1;
						A[1][ad0].type = "SHOCK";

					}
					else
						ad0 = findAd(getCrossPoint(L1, L12), A[1]);

					if (findAd(getCrossPoint(L2, L12), A[1]) == -1)
					{
						A[1].push_back(getCrossPoint(L2, L12));
						ad4 = A[1].size() - 1;
						A[1][ad4].type = "SHOCK";

					}
					else
						ad4 = findAd(getCrossPoint(L2, L12), A[1]);
					ad1 = findAd(A0[n2], A[1]);
					ad2 = findAd(A0[n3], A[1]);
					ad3 = findAd(A0[n4], A[1]);
					M[1].push_back(m);
					s = int(M[1].size()) - 1;
					M[1][s].node.push_back(ad0);
					M[1][s].node.push_back(ad1);
					M[1][s].node.push_back(ad2);
					M[1][s].node.push_back(ad3);
					M[1][s].node.push_back(ad4);

					M[1][s].face_start.push_back(ad0);
					M[1][s].face_start.push_back(ad1);
					M[1][s].face_start.push_back(ad2);
					M[1][s].face_start.push_back(ad3);
					M[1][s].face_start.push_back(ad4);

					M[1][s].face_end.push_back(ad1);
					M[1][s].face_end.push_back(ad2);
					M[1][s].face_end.push_back(ad3);
					M[1][s].face_end.push_back(ad4);
					M[1][s].face_end.push_back(ad0);

				}
				else if (A0[n0].section == 3 && A0[n1].section == 3 && A0[n2].section == 5 && (A0[n3].section == 5 || A0[n3].section == 45) && (A0[n4].section == 5 || A0[n4].section == 45) && A0[n5].section == 3)
				{
					ad0 = findAd(A0[n0], A[2]);
					ad1 = findAd(A0[n1], A[2]);

					L1 = getLine(A0[n1], A0[n2]);
					L2 = getLine(A0[n4], A0[n5]);
					if (findAd(getCrossPoint(L1, L12), A[2]) == -1)
					{
						A[2].push_back(getCrossPoint(L1, L12));
						ad2 = A[2].size() - 1;
						A[2][ad2].type = "SHOCK";

					}
					else
						ad2 = findAd(getCrossPoint(L1, L12), A[2]);

					if (findAd(getCrossPoint(L2, L12), A[2]) == -1)
					{
						A[2].push_back(getCrossPoint(L2, L12));
						ad3 = A[2].size() - 1;
						A[2][ad3].type = "SHOCK";

					}
					else
						ad3 = findAd(getCrossPoint(L2, L12), A[2]);
					ad4 = findAd(A0[n5], A[2]);
					M[2].push_back(m);
					s = int(M[2].size()) - 1;
					M[2][s].node.push_back(ad0);
					M[2][s].node.push_back(ad1);
					M[2][s].node.push_back(ad2);
					M[2][s].node.push_back(ad3);
					M[2][s].node.push_back(ad4);

					M[2][s].face_start.push_back(ad0);
					M[2][s].face_start.push_back(ad1);
					M[2][s].face_start.push_back(ad2);
					M[2][s].face_start.push_back(ad3);
					M[2][s].face_start.push_back(ad4);

					M[2][s].face_end.push_back(ad1);
					M[2][s].face_end.push_back(ad2);
					M[2][s].face_end.push_back(ad3);
					M[2][s].face_end.push_back(ad4);
					M[2][s].face_end.push_back(ad0);


					if (findAd(getCrossPoint(L1, L12), A[4]) == -1)
					{
						A[4].push_back(getCrossPoint(L1, L12));
						ad0 = A[4].size() - 1;
						A[4][ad0].type = "SHOCK";

					}
					else
						ad0 = findAd(getCrossPoint(L1, L12), A[4]);

					if (findAd(getCrossPoint(L2, L12), A[4]) == -1)
					{
						A[4].push_back(getCrossPoint(L2, L12));
						ad4 = A[4].size() - 1;
						A[4][ad4].type = "SHOCK";

					}
					else
						ad4 = findAd(getCrossPoint(L2, L12), A[4]);
					ad1 = findAd(A0[n2], A[4]);
					ad2 = findAd(A0[n3], A[4]);
					ad3 = findAd(A0[n4], A[4]);
					M[4].push_back(m);
					s = int(M[4].size()) - 1;
					M[4][s].node.push_back(ad0);
					M[4][s].node.push_back(ad1);
					M[4][s].node.push_back(ad2);
					M[4][s].node.push_back(ad3);
					M[4][s].node.push_back(ad4);

					M[4][s].face_start.push_back(ad0);
					M[4][s].face_start.push_back(ad1);
					M[4][s].face_start.push_back(ad2);
					M[4][s].face_start.push_back(ad3);
					M[4][s].face_start.push_back(ad4);

					M[4][s].face_end.push_back(ad1);
					M[4][s].face_end.push_back(ad2);
					M[4][s].face_end.push_back(ad3);
					M[4][s].face_end.push_back(ad4);
					M[4][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 3 && A0[n1].section == 3 && A0[n2].section == 3 && A0[n3].section == 1 && A0[n4].section == 1 && A0[n5].section == 1)
				{
					L1 = getLine(A0[n0], A0[n5]);
					L2 = getLine(A0[n2], A0[n3]);

					if (findAd(getCrossPoint(L1, L13), A[0]) == -1)
					{
						A[0].push_back(getCrossPoint(L1, L13));
						ad0 = A[0].size() - 1;
						A[0][ad0].type = "SHOCK";

					}
					else
						ad0 = findAd(getCrossPoint(L1, L13), A[0]);

					if (findAd(getCrossPoint(L2, L13), A[0]) == -1)
					{
						A[0].push_back(getCrossPoint(L2, L13));
						ad1 = A[0].size() - 1;
						A[0][ad1].type = "SHOCK";

					}
					else
						ad1 = findAd(getCrossPoint(L2, L13), A[0]);

					ad2 = findAd(A0[n3], A[0]);
					ad3 = findAd(A0[n4], A[0]);
					ad4 = findAd(A0[n5], A[0]);
					M[0].push_back(m);
					s = int(M[0].size()) - 1;
					M[0][s].node.push_back(ad0);
					M[0][s].node.push_back(ad1);
					M[0][s].node.push_back(ad2);
					M[0][s].node.push_back(ad3);
					M[0][s].node.push_back(ad4);

					M[0][s].face_start.push_back(ad0);
					M[0][s].face_start.push_back(ad1);
					M[0][s].face_start.push_back(ad2);
					M[0][s].face_start.push_back(ad3);
					M[0][s].face_start.push_back(ad4);

					M[0][s].face_end.push_back(ad1);
					M[0][s].face_end.push_back(ad2);
					M[0][s].face_end.push_back(ad3);
					M[0][s].face_end.push_back(ad4);
					M[0][s].face_end.push_back(ad0);

					if (findAd(getCrossPoint(L2, L13), A[2]) == -1)
					{
						A[2].push_back(getCrossPoint(L2, L13));
						ad3 = A[2].size() - 1;
						A[2][ad3].type = "SHOCK";

					}
					else
						ad3 = findAd(getCrossPoint(L2, L13), A[2]);

					if (findAd(getCrossPoint(L1, L13), A[2]) == -1)
					{
						A[2].push_back(getCrossPoint(L1, L13));
						ad4 = A[2].size() - 1;
						A[2][ad4].type = "SHOCK";

					}
					else
						ad4 = findAd(getCrossPoint(L1, L13), A[2]);

					ad0 = findAd(A0[n0], A[2]);
					ad1 = findAd(A0[n1], A[2]);
					ad2 = findAd(A0[n2], A[2]);
					M[2].push_back(m);
					s = int(M[2].size()) - 1;
					M[2][s].node.push_back(ad0);
					M[2][s].node.push_back(ad1);
					M[2][s].node.push_back(ad2);
					M[2][s].node.push_back(ad3);
					M[2][s].node.push_back(ad4);

					M[2][s].face_start.push_back(ad0);
					M[2][s].face_start.push_back(ad1);
					M[2][s].face_start.push_back(ad2);
					M[2][s].face_start.push_back(ad3);
					M[2][s].face_start.push_back(ad4);

					M[2][s].face_end.push_back(ad1);
					M[2][s].face_end.push_back(ad2);
					M[2][s].face_end.push_back(ad3);
					M[2][s].face_end.push_back(ad4);
					M[2][s].face_end.push_back(ad0);
				}
				else if ((A0[n0].section == 4 || A0[n0].section == 45) && (A0[n1].section == 4 || A0[n1].section == 45) && A0[n2].section == 4 && A0[n3].section == 2 && A0[n4].section == 2 && A0[n5].section == 2)
				{
					L1 = getLine(A0[n0], A0[n5]);
					L2 = getLine(A0[n2], A0[n3]);

					if (findAd(getCrossPoint(L1, L13), A[1]) == -1)
					{
						A[1].push_back(getCrossPoint(L1, L13));
						ad0 = A[1].size() - 1;
						A[1][ad0].type = "SHOCK";

					}
					else
						ad0 = findAd(getCrossPoint(L1, L13), A[1]);

					if (findAd(getCrossPoint(L2, L13), A[1]) == -1)
					{
						A[1].push_back(getCrossPoint(L2, L13));
						ad1 = A[1].size() - 1;
						A[1][ad1].type = "SHOCK";

					}
					else
						ad1 = findAd(getCrossPoint(L2, L13), A[1]);

					ad2 = findAd(A0[n3], A[1]);
					ad3 = findAd(A0[n4], A[1]);
					ad4 = findAd(A0[n5], A[1]);
					M[1].push_back(m);
					s = int(M[1].size()) - 1;
					M[1][s].node.push_back(ad0);
					M[1][s].node.push_back(ad1);
					M[1][s].node.push_back(ad2);
					M[1][s].node.push_back(ad3);
					M[1][s].node.push_back(ad4);

					M[1][s].face_start.push_back(ad0);
					M[1][s].face_start.push_back(ad1);
					M[1][s].face_start.push_back(ad2);
					M[1][s].face_start.push_back(ad3);
					M[1][s].face_start.push_back(ad4);

					M[1][s].face_end.push_back(ad1);
					M[1][s].face_end.push_back(ad2);
					M[1][s].face_end.push_back(ad3);
					M[1][s].face_end.push_back(ad4);
					M[1][s].face_end.push_back(ad0);

					if (findAd(getCrossPoint(L2, L13), A[3]) == -1)
					{
						A[3].push_back(getCrossPoint(L2, L13));
						ad3 = A[3].size() - 1;
						A[3][ad3].type = "SHOCK";

					}
					else
						ad3 = findAd(getCrossPoint(L2, L13), A[3]);

					if (findAd(getCrossPoint(L1, L13), A[3]) == -1)
					{
						A[3].push_back(getCrossPoint(L1, L13));
						ad4 = A[3].size() - 1;
						A[3][ad4].type = "SHOCK";

					}
					else
						ad4 = findAd(getCrossPoint(L1, L13), A[3]);

					ad0 = findAd(A0[n0], A[3]);
					ad1 = findAd(A0[n1], A[3]);
					ad2 = findAd(A0[n2], A[3]);
					M[3].push_back(m);
					s = int(M[3].size()) - 1;
					M[3][s].node.push_back(ad0);
					M[3][s].node.push_back(ad1);
					M[3][s].node.push_back(ad2);
					M[3][s].node.push_back(ad3);
					M[3][s].node.push_back(ad4);

					M[3][s].face_start.push_back(ad0);
					M[3][s].face_start.push_back(ad1);
					M[3][s].face_start.push_back(ad2);
					M[3][s].face_start.push_back(ad3);
					M[3][s].face_start.push_back(ad4);

					M[3][s].face_end.push_back(ad1);
					M[3][s].face_end.push_back(ad2);
					M[3][s].face_end.push_back(ad3);
					M[3][s].face_end.push_back(ad4);
					M[3][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 5 && A0[n1].section == 5 && A0[n2].section == 45 && A0[n3].section == 4 && A0[n4].section == 4 && A0[n5].section == 45)
				{

					ad0 = findAd(A0[n2], A[3]);
					ad1 = findAd(A0[n3], A[3]);
					ad2 = findAd(A0[n4], A[3]);
					ad3 = findAd(A0[n5], A[3]);
					M[3].push_back(m);
					s = int(M[3].size()) - 1;
					M[3][s].node.push_back(ad0);
					M[3][s].node.push_back(ad1);
					M[3][s].node.push_back(ad2);
					M[3][s].node.push_back(ad3);

					M[3][s].face_start.push_back(ad0);
					M[3][s].face_start.push_back(ad1);
					M[3][s].face_start.push_back(ad2);
					M[3][s].face_start.push_back(ad3);

					M[3][s].face_end.push_back(ad1);
					M[3][s].face_end.push_back(ad2);
					M[3][s].face_end.push_back(ad3);
					M[3][s].face_end.push_back(ad0);


					ad0 = findAd(A0[n0], A[4]);
					ad1 = findAd(A0[n1], A[4]);
					ad2 = findAd(A0[n2], A[4]);
					ad3 = findAd(A0[n5], A[4]);
					M[4].push_back(m);
					s = int(M[4].size()) - 1;
					M[4][s].node.push_back(ad0);
					M[4][s].node.push_back(ad1);
					M[4][s].node.push_back(ad2);
					M[4][s].node.push_back(ad3);

					M[4][s].face_start.push_back(ad0);
					M[4][s].face_start.push_back(ad1);
					M[4][s].face_start.push_back(ad2);
					M[4][s].face_start.push_back(ad3);

					M[4][s].face_end.push_back(ad1);
					M[4][s].face_end.push_back(ad2);
					M[4][s].face_end.push_back(ad3);
					M[4][s].face_end.push_back(ad0);
				}
				else if (A0[n0].section == 5 && A0[n1].section == 5 && A0[n2].section == 5 && A0[n3].section == 45 && A0[n4].section == 45 && A0[n5].section == 5)
				{
					ad0 = findAd(A0[n0], A[4]);
					ad1 = findAd(A0[n1], A[4]);
					ad2 = findAd(A0[n2], A[4]);
					ad3 = findAd(A0[n3], A[4]);
					ad4 = findAd(A0[n4], A[4]);
					ad5 = findAd(A0[n5], A[4]);
					M[4].push_back(m);
					s = int(M[4].size()) - 1;
					M[4][s].node.push_back(ad0);
					M[4][s].node.push_back(ad1);
					M[4][s].node.push_back(ad2);
					M[4][s].node.push_back(ad3);
					M[4][s].node.push_back(ad4);
					M[4][s].node.push_back(ad5);

					M[4][s].face_start.push_back(ad0);
					M[4][s].face_start.push_back(ad1);
					M[4][s].face_start.push_back(ad2);
					M[4][s].face_start.push_back(ad3);
					M[4][s].face_start.push_back(ad4);
					M[4][s].face_start.push_back(ad5);

					M[4][s].face_end.push_back(ad1);
					M[4][s].face_end.push_back(ad2);
					M[4][s].face_end.push_back(ad3);
					M[4][s].face_end.push_back(ad4);
					M[4][s].face_end.push_back(ad5);
					M[4][s].face_end.push_back(ad0);

				}
				else if (A0[n0].section == 45 && A0[n1].section == 45 && A0[n2].section == 4 && A0[n3].section == 4 && A0[n4].section == 4 && A0[n5].section == 4)
				{
					ad0 = findAd(A0[n0], A[3]);
					ad1 = findAd(A0[n1], A[3]);
					ad2 = findAd(A0[n2], A[3]);
					ad3 = findAd(A0[n3], A[3]);
					ad4 = findAd(A0[n4], A[3]);
					ad5 = findAd(A0[n5], A[3]);
					M[3].push_back(m);
					s = int(M[3].size()) - 1;
					M[3][s].node.push_back(ad0);
					M[3][s].node.push_back(ad1);
					M[3][s].node.push_back(ad2);
					M[3][s].node.push_back(ad3);
					M[3][s].node.push_back(ad4);
					M[3][s].node.push_back(ad5);

					M[3][s].face_start.push_back(ad0);
					M[3][s].face_start.push_back(ad1);
					M[3][s].face_start.push_back(ad2);
					M[3][s].face_start.push_back(ad3);
					M[3][s].face_start.push_back(ad4);
					M[3][s].face_start.push_back(ad5);

					M[3][s].face_end.push_back(ad1);
					M[3][s].face_end.push_back(ad2);
					M[3][s].face_end.push_back(ad3);
					M[3][s].face_end.push_back(ad4);
					M[3][s].face_end.push_back(ad5);
					M[3][s].face_end.push_back(ad0);

				}
				else if (A0[n0].section == 3 && A0[n1].section == 3 && A0[n2].section == 45 && A0[n3].section == 2 && A0[n4].section == 2 && A0[n5].section == 1)
				{
					L1 = getLine(A0[n0], A0[n5]);
					L2 = getLine(A0[n4], A0[n5]);

					if (findAd(getCrossPoint(L1, L13), A[0]) == -1)
					{
						A[0].push_back(getCrossPoint(L1, L13));
						ad0 = A[0].size() - 1;
						A[0][ad0].type = "SHOCK";

					}
					else
						ad0 = findAd(getCrossPoint(L1, L13), A[0]);

					ad1 = 0;

					if (findAd(getCrossPoint(L2, L12), A[0]) == -1)
					{
						A[0].push_back(getCrossPoint(L2, L12));
						ad2 = A[0].size() - 1;
						A[0][ad2].type = "SHOCK";

					}
					else
						ad2 = findAd(getCrossPoint(L1, L12), A[0]);
					ad3 = findAd(A0[n5], A[0]);
					M[0].push_back(m);
					s = int(M[0].size()) - 1;
					M[0][s].node.push_back(ad0);
					M[0][s].node.push_back(ad1);
					M[0][s].node.push_back(ad2);
					M[0][s].node.push_back(ad3);
					M[0][s].face_start.push_back(ad0);
					M[0][s].face_start.push_back(ad1);
					M[0][s].face_start.push_back(ad2);
					M[0][s].face_start.push_back(ad3);
					M[0][s].face_end.push_back(ad1);
					M[0][s].face_end.push_back(ad2);
					M[0][s].face_end.push_back(ad3);
					M[0][s].face_end.push_back(ad0);

					//M[1]
					ad0 = 0;
					ad2 = findAd(A0[n3], A[1]);
					ad3 = findAd(A0[n4], A[1]);
					L1 = getLine(A0[n2], A0[n3]);
					L2 = getLine(A0[n4], A0[n5]);
					if (findAd(getCrossPoint(L1, L13), A[1]) == -1)
					{
						A[1].push_back(getCrossPoint(L1, L13));
						ad1 = A[1].size() - 1;
						A[1][ad1].type = "SHOCK";

					}
					else
						ad1 = findAd(getCrossPoint(L1, L13), A[1]);

					if (findAd(getCrossPoint(L2, L12), A[1]) == -1)
					{
						A[1].push_back(getCrossPoint(L2, L12));
						ad4 = A[1].size() - 1;
						A[1][ad4].type = "SHOCK";

					}
					else
						ad4 = findAd(getCrossPoint(L1, L12), A[1]);
					M[1].push_back(m);
					s = int(M[1].size()) - 1;
					M[1][s].node.push_back(ad0);
					M[1][s].node.push_back(ad1);
					M[1][s].node.push_back(ad2);
					M[1][s].node.push_back(ad3);
					M[1][s].node.push_back(ad4);

					M[1][s].face_start.push_back(ad0);
					M[1][s].face_start.push_back(ad1);
					M[1][s].face_start.push_back(ad2);
					M[1][s].face_start.push_back(ad3);
					M[1][s].face_start.push_back(ad4);

					M[1][s].face_end.push_back(ad1);
					M[1][s].face_end.push_back(ad2);
					M[1][s].face_end.push_back(ad3);
					M[1][s].face_end.push_back(ad4);
					M[1][s].face_end.push_back(ad0);

					//M[2]
					ad0 = findAd(A0[n0], A[2]);
					ad1 = findAd(A0[n1], A[2]);
					L1 = getLine(A0[n1], A0[n2]);
					L2 = getLine(A0[n0], A0[n5]);
					if (findAd(getCrossPoint(L1, L12), A[2]) == -1)
					{
						A[2].push_back(getCrossPoint(L1, L12));
						ad2 = A[2].size() - 1;
						A[2][ad2].type = "SHOCK";

					}
					else
						ad2 = findAd(getCrossPoint(L1, L12), A[2]);

					ad3 = 0;

					if (findAd(getCrossPoint(L2, L13), A[2]) == -1)
					{
						A[2].push_back(getCrossPoint(L2, L13));
						ad4 = A[2].size() - 1;
						A[2][ad4].type = "SHOCK";

					}
					else
						ad4 = findAd(getCrossPoint(L2, L13), A[2]);
					M[2].push_back(m);
					s = int(M[2].size()) - 1;
					M[2][s].node.push_back(ad0);
					M[2][s].node.push_back(ad1);
					M[2][s].node.push_back(ad2);
					M[2][s].node.push_back(ad3);
					M[2][s].node.push_back(ad4);

					M[2][s].face_start.push_back(ad0);
					M[2][s].face_start.push_back(ad1);
					M[2][s].face_start.push_back(ad2);
					M[2][s].face_start.push_back(ad3);
					M[2][s].face_start.push_back(ad4);

					M[2][s].face_end.push_back(ad1);
					M[2][s].face_end.push_back(ad2);
					M[2][s].face_end.push_back(ad3);
					M[2][s].face_end.push_back(ad4);
					M[2][s].face_end.push_back(ad0);
					//M[3]
					ad0 = 0;
					ad1 = findAd(A0[n2], A[3]);
					L1 = getLine(A0[n2], A0[n3]);
					if (findAd(getCrossPoint(L1, L13), A[3]) == -1)
					{
						A[3].push_back(getCrossPoint(L1, L13));
						ad2 = A[3].size() - 1;
						A[3][ad2].type = "SHOCK";

					}
					else
						ad2 = findAd(getCrossPoint(L1, L13), A[3]);

					M[3].push_back(m);
					s = int(M[3].size()) - 1;
					M[3][s].node.push_back(ad0);
					M[3][s].node.push_back(ad1);
					M[3][s].node.push_back(ad2);
					M[3][s].face_start.push_back(ad0);
					M[3][s].face_start.push_back(ad1);
					M[3][s].face_start.push_back(ad2);
					M[3][s].face_end.push_back(ad1);
					M[3][s].face_end.push_back(ad2);
					M[3][s].face_end.push_back(ad0);
					//M[4]
					ad0 = 0;
					ad1 = findAd(A0[n2], A[4]);
					L1 = getLine(A0[n1], A0[n2]);
					if (findAd(getCrossPoint(L1, L12), A[4]) == -1)
					{
						A[4].push_back(getCrossPoint(L1, L12));
						ad2 = A[4].size() - 1;
						A[4][ad2].type = "SHOCK";

					}
					else
						ad2 = findAd(getCrossPoint(L1, L12), A[4]);

					M[4].push_back(m);
					s = int(M[4].size()) - 1;
					M[4][s].node.push_back(ad0);
					M[4][s].node.push_back(ad1);
					M[4][s].node.push_back(ad2);
					M[4][s].face_start.push_back(ad0);
					M[4][s].face_start.push_back(ad1);
					M[4][s].face_start.push_back(ad2);
					M[4][s].face_end.push_back(ad1);
					M[4][s].face_end.push_back(ad2);
					M[4][s].face_end.push_back(ad0);


				}

				else
					std::cout << "something wrong in partitionMesh !" << std::endl;
			}
		}

	}


}

void out_M1(std::string name)
{
	int face = 0;
	for (int i = 0; i < M[0].size(); i++)
	{
		face += M[0][i].node.size();
	}
	std::ofstream fout;
	using std::endl;
	fout.open(name + ".dat");
	int i;
	//fout << "FILETYPE = GRID" << endl;
	//fout << "VARIABLES = \"X\", \"Y\"" << endl;
	fout << "VARIABLES =  \"X\", \"Y\"" << std::endl;

	fout << "ZONE T=\"Test\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "Nodes = " << A[0].size() << endl;
	fout << "Elements = " << M[0].size() << endl;
	fout << "Faces = " << face << endl;
	fout << "NumConnectedBoundaryFaces=0 " << endl;
	fout << "TotalNumBoundaryConnections=0 " << endl;

	fout.scientific;
	for (i = 0; i < A[0].size(); i++)
	{
		fout << A[0][i].x << endl;
	}
	fout << endl;

	for (i = 0; i < A[0].size(); i++)
	{
		fout << A[0][i].y << endl;
	}
	fout << endl;

	for (i = 0; i < M[0].size(); i++)
	{
		for (int j = 0; j < M[0][i].face_start.size(); j++)
		{
			fout << M[0][i].face_start[j] + 1 << "   " << M[0][i].face_end[j] + 1 << endl;
		}
	}
	fout << endl;
	for (i = 0; i < M[0].size(); i++)
	{
		for (int j = 0; j < M[0][i].face_start.size(); j++)
			fout << i + 1 << "  ";
		fout << endl;
	}
	fout << endl;
	for (i = 0; i < M[0].size(); i++)
	{
		for (int j = 0; j < M[0][i].face_start.size(); j++)
			fout << 0 << "  ";
		fout << endl;
	}
}
void out_M2(std::string name)
{
	int face = 0;
	for (int i = 0; i < M[1].size(); i++)
	{
		face += M[1][i].node.size();
	}
	std::ofstream fout;
	using std::endl;
	fout.open(name + ".dat");
	int i;
	//fout << "FILETYPE = GRID" << endl;
	//fout << "VARIABLES = \"X\", \"Y\"" << endl;
	fout << "VARIABLES =  \"X\", \"Y\"" << std::endl;

	fout << "ZONE T=\"Test\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "Nodes = " << A[1].size() << endl;
	fout << "Elements = " << M[1].size() << endl;
	fout << "Faces = " << face << endl;
	fout << "NumConnectedBoundaryFaces=0 " << endl;
	fout << "TotalNumBoundaryConnections=0 " << endl;

	fout.scientific;
	for (i = 0; i < A[1].size(); i++)
	{
		fout << A[1][i].x << endl;
	}
	fout << endl;

	for (i = 0; i < A[1].size(); i++)
	{
		fout << A[1][i].y << endl;
	}
	fout << endl;

	for (i = 0; i < M[1].size(); i++)
	{
		for (int j = 0; j < M[1][i].face_start.size(); j++)
		{
			fout << M[1][i].face_start[j] + 1 << "   " << M[1][i].face_end[j] + 1 << endl;
		}
	}
	fout << endl;
	for (i = 0; i < M[1].size(); i++)
	{
		for (int j = 0; j < M[1][i].face_start.size(); j++)
			fout << i + 1 << "  ";
		fout << endl;
	}
	fout << endl;
	for (i = 0; i < M[1].size(); i++)
	{
		for (int j = 0; j < M[1][i].face_start.size(); j++)
			fout << 0 << "  ";
		fout << endl;
	}

}

