#include"const.h"
#include"functions.h"
#include"shockwave.h"
#include<iostream>
#include<ctime>
#include<stdlib.h>
#include<vector>
#include"Prandtl-Meyer.h"
#include"init.h"




using std::vector;
using namespace MeshPara;
void init_mesh()
{
	extern vector <mesh> A0;
	extern vector<vector<int>> ad;
	mesh t;
	vector<int> a;
	int i, j;
	i = 0;
	extern double xL, xR, yU, yD;

	if (FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal")
	{
		if (meshType == 10 || meshType == 11 || meshType == 12)
		{
			while (i < Pnum)
			{
				for (j = 0; j < Xnum; j++)
				{
					A0.push_back(t);

					if (j == 0)
					{
						if (i == 0)
							A0[i].x = dx / 2, A0[i].y = 0;
						else if (A0[i - Xnum].x == 0)
							A0[i].x = dx / 2, A0[i].y = A0[i - Xnum].y + dy;
						else
							A0[i].x = 0, A0[i].y = A0[i - Xnum].y + dy;
					}
					else if (A0[i - j].x == 0)
					{
						if (j % 2 == 0)
							A0[i].x = A0[i - 1].x + dx, A0[i].y = A0[i - 1].y;
						else
							A0[i].x = A0[i - 1].x + 2 * dx, A0[i].y = A0[i - 1].y;
					}
					else
					{
						if (j % 2 == 0)
							A0[i].x = A0[i - 1].x + 2 * dx, A0[i].y = A0[i - 1].y;
						else
							A0[i].x = A0[i - 1].x + dx, A0[i].y = A0[i - 1].y;
					}
					i++;
				}
			}
			xL = min(A0[0].x, A0[Xnum].x);
			xR = max(A0[Xnum - 1].x, A0[Xnum - 1 + Xnum].x);
			yU = A0[Pnum - 1].y;
			yD = A0[0].y;
		}
		else
		{
			xL = 0;
			xR = 2 * dx + 3 * dx * (Xnum - 2) / 2;
			yU = dy * (Ynum - 1);
			yD = 0;
			double dx1 = (xR - xL) / (Xnum - 1);
			double dy1 = (yU - yD) / (Ynum - 1);
			while (i < Pnum)
			{
				for (j = 0; j < Xnum; j++)
				{
					A0.push_back(t);

					if (j == 0)
					{
						if (i == 0)
							A0[i].x = 0, A0[i].y = 0;
						else
							A0[i].x = 0, A0[i].y = A0[i - Xnum].y + dy1;
					}
					else
						A0[i].x = A0[i - 1].x + dx1, A0[i].y = A0[i - 1].y;
					i++;
				}
			}
		}
		//for (i = 0; i < Pnum; i++)
		//{
		//	if (A0[i].y > yD)
		//		A0[i].y -= dy;
		//	if (abs(A0[i].y - yU) < 1e-10)
		//		A0[i].y -= 2 * dy;
		//}
		//xL = min(A0[0].x, A0[Xnum].x);
		//xR = max(A0[Xnum - 1].x, A0[Xnum - 1 + Xnum].x);
		//yU = A0[Pnum - 1].y;
		//yD = A0[0].y;

	}
	if (FlowType == "Prandtl-Meyer")
	{
		using namespace Prandtl_Meyer;
		using ConstPara::��;
		double LineLength = 0.003;
		//��߽ṹ��������
		double dx;
		double dl = LineLength / Ynum;
		double ex = cos(��1);
		double ey = sin(��1);
		double xend, yend;
		while (i < Pnum)
		{
			if (i == 0)
			{
				xend = LineLength;
				yend = 0;
			}
			else
			{
				xend = dl * ex + A0[i - 1].x;
				yend = dl * ey + A0[i - 1].y;
			}
			dx = xend / (Xnum - 1);
			for (j = 0; j < Xnum; j++)
			{
				A0.push_back(t);
				if (j == 0)
				{
					if (i == 0)
						A0[i].x = 0, A0[i].y = 0;
					else
						A0[i].x = 0, A0[i].y = yend;
				}
				else
				{
					A0[i].x = A0[i - 1].x + dx, A0[i].y = A0[i - 1].y;
				}
				i++;
			}
		}
		//�����������
		int k = 0;
		int s = 0;
		while (k < (Xnum - 1) * (Ynum - 1))
		{
			for (j = 0; j < Xnum - 1; j++)
			{
				ad.push_back(a);
				ad[k].push_back(k + s);
				ad[k].push_back(k + 1 + s);
				ad[k].push_back(k + 1 + Xnum + s);
				ad[k].push_back(k + Xnum + s);
				if (j == Xnum - 2)
					s++;
				k++;
			}
		}
		//�ұ߽ṹ��������
		ex = cos(��2 - ��2);
		ey = sin(��2 - ��2);
		double xstart, ystart;
		double dl2 = 1.5 * LineLength * sin(��2) / Ynum;
		s = 0;
		while (i < 2 * Pnum - 2)
		{
			if (i == Pnum)
			{
				xstart = LineLength;
				ystart = 0;
				xend = LineLength + 2 * LineLength * cos(��2);
				yend = -2 * LineLength * sin(��2);
			}
			else
			{
				if (s == 1)
				{
					xstart = 1.3 * dl * ex + A0[Xnum - 1].x;
					ystart = 1.3 * dl * ey + A0[Xnum - 1].y;
				}
				else
				{
					xstart = 1.3 * dl * ex + A0[i - Xnum].x;
					ystart = 1.3 * dl * ey + A0[i - Xnum].y;
				}
				xend = dl2 * sin(��2) + A0[i - 1].x;
				yend = dl2 * cos(��2) + A0[i - 1].y;
			}
			double d = sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart)) / (Xnum - 1);
			double ex1 = (xend - xstart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
			double ey1 = (yend - ystart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
			for (j = 0; j < Xnum; j++)
			{
				if (j == 0 && i == Pnum)
					continue;
				A0.push_back(t);
				if (j == 0)
					A0[i].x = xstart, A0[i].y = ystart;
				else
				{
					if (i == Pnum)
						A0[i].x = A0[Xnum - 1].x + d * ex1, A0[i].y = A0[Xnum - 1].y + d * ey1;
					else
						A0[i].x = A0[i - 1].x + d * ex1, A0[i].y = A0[i - 1].y + d * ey1;
				}
				i++;
				if (i >= 2 * Pnum - 1)
					break;
			}
			s++;
		}
		//�ұ���������
		k = ad.size();
		s = Xnum + Ynum - 1;
		while (k < 2 * (Xnum - 1) * (Ynum - 1))
		{
			for (j = 0; j < Xnum - 1; j++)
			{
				ad.push_back(a);
				if (k == (Xnum - 1) * (Ynum - 1))
					ad[k].push_back(Xnum - 1);
				else
					ad[k].push_back(k + s - 1);
				ad[k].push_back(k + s);
				ad[k].push_back(k + Xnum + s);
				ad[k].push_back(k + Xnum + s - 1);
				if (j == Xnum - 2)
					s++;
				k++;
			}
		}
		//�м�нǷǽṹ����
		for (j = 0; j < Ynum; j++)
		{
			if (j == 0)
				continue;
			if (j == 1)
			{
				ad.push_back(a);
				ad[ad.size() - 1].push_back(2 * Xnum - 1);
				ad[ad.size() - 1].push_back(Xnum - 1);
				ad[ad.size() - 1].push_back(Xnum - 1 + Pnum);
			}
			else
			{
				//�������ʼ��ͽ�����
				xstart = A0[(j + 1) * Xnum - 1].x;
				ystart = A0[(j + 1) * Xnum - 1].y;
				xend = A0[j * Xnum - 1 + Pnum].x;
				yend = A0[j * Xnum - 1 + Pnum].y;
				//��ʼ�㵽������ĵȷֳ��ȣ��Լ���λ����������xyֵ
				double d = sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart)) / j;
				double ex1 = (xend - xstart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
				double ey1 = (yend - ystart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
				int size = A0.size();
				for (k = 0; k < j - 1; k++)
				{
					A0.push_back(t);
					A0[A0.size() - 1].x = xstart + (k + 1) * d * ex1;
					A0[A0.size() - 1].y = ystart + (k + 1) * d * ey1;
					if (j == 2)
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back((j + 1) * Xnum - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);
						ad[ad.size() - 1].push_back(j * Xnum - 1 + Pnum);
					}
					else if (k == 0)
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back((j + 1) * Xnum - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back(size - j + 2);
					}
					else if (k == j - 2)
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);
						ad[ad.size() - 1].push_back(j * Xnum - 1 + Pnum);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back(A0.size() - 2);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);

					}
					else
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back(A0.size() - 2);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(A0.size() - 1);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad[ad.size() - 1].push_back(size - j + 2 + k);
					}
				}
			}

		}

	}
}

void reorderMesh()//���߽��Ϊ����
{
	extern vector <mesh> A0;
	extern vector<vector<int>> ad;
	extern double xL, xR, yU, yD;
	xL = min(A0[0].x, A0[Xnum].x);
	xR = max(A0[Xnum - 1].x, A0[Xnum - 1 + Xnum].x);
	yU = A0[Pnum - 1].y;
	yD = A0[0].y;
	using ConstPara::pi;
	int i = 0;
	if ((FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal") && meshType != 20)
	{
		for (i = 0; i < A0.size(); i++)
		{
			if (abs(A0[i].x - xL) <= 1e-10 || abs(A0[i].x - xR) <= 1e-10)
				continue;
			if (abs(A0[i].x - (xL + dx * cos(pi / 3))) <= 1e-10)
				A0[i].x = xL, A0[i].type = "L";
			else if (abs(A0[i].x - (xR - dx * cos(pi / 3))) <= 1e-10)
				A0[i].x = xR, A0[i].type = "R";
			else if (abs(A0[i].y - (yD + dx * sin(pi / 3))) <= 1e-10)
				A0[i].y = yD, A0[i].type = "D";
			else if (abs(A0[i].y - (yU - dx * sin(pi / 3))) <= 1e-10)
				A0[i].y = yU, A0[i].type = "U";
		}
	}
}
void remesh_bound()//���߽�x��y�����ƶ������ڲ�����ͬ�����Ʊ߽�����
{
	extern vector<vector <mesh>> A;
	int i, j, k, m;
	if (FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal")
	{
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "D" || A[i][j].type == "U")
				{
					for (k = 0; k < A[i][j].neibor.size(); k++)
					{
						m = A[i][j].neibor[k];
						if (abs(A[i][j].y - A[i][m].y) < 1e-10)
							continue;
						else
							A[i][j].x = A[i][m].x;
					}
				}
			}
		}
	}
}
void init_mesh1()
{
	extern vector <mesh> A0;
	int i, j;
	i = 0;
	mesh t;
	while (i < Pnum)
	{
		for (j = 0; j < Xnum; j++)
		{
			A0.push_back(t);

			if (j == 0)
			{
				if (i == 0)
					A0[i].x = 0, A0[i].y = 0;
				else if (A0[i - Xnum].x == 0)
					A0[i].x = 0, A0[i].y = A0[i - Xnum].y + dy;
				else
					A0[i].x = 0, A0[i].y = A0[i - Xnum].y + dy;
			}
			else if (A0[i - j].x == 0)
			{
				if (j % 2 == 0)
					A0[i].x = A0[i - 1].x + dx, A0[i].y = A0[i - 1].y;
				else
					A0[i].x = A0[i - 1].x + dx, A0[i].y = A0[i - 1].y;
			}
			else
			{
				if (j % 2 == 0)
					A0[i].x = A0[i - 1].x + dx, A0[i].y = A0[i - 1].y;
				else
					A0[i].x = A0[i - 1].x + dx, A0[i].y = A0[i - 1].y;
			}
			A0[i].x0 = A0[i].x;
			A0[i].y0 = A0[i].y;
			i++;
		}
	}
}
void brick_mesh()
{
	extern vector<vector <mesh>> A;
	int i, j, k, m, n;
	for (i = 0; i < A.size(); i++)
	{
		n = 0;
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].type != "IN")
				continue;
			n = 0;
			for (k = 0; k < A[i][j].neibor.size(); k++)
			{
				m = A[i][j].neibor[k];
				if (m == j + 1)
				{
					n++;
					break;
				}
			}
			if (n == 1)
				A[i][j].x = A[i][j].x - dx / 4;
			else
				A[i][j].x = A[i][j].x + dx / 4;
		}
	}
}
void change_mesh()
{
	extern vector<vector <mesh>> A;
	extern double xL, xR;
	using ConstPara::random;
	//srand((unsigned)time(NULL));//�����Ŷ�
	int i, j;
	double r;
	double sumx, sumy;
	int m, n;
	int k;

	if (meshType == 10)
	{
		if (FlowType == "intersection")
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					if (abs(A[i][j].x - 0.01685) <= dx / 2)
					{
						A[i][j].x = 0.01685;
					}
				}
			}

	}
	else if (meshType == 11)
	{
		if (FlowType == "normal" || FlowType == "oblique")
		{
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j].x == xL || A[i][j].x == xR || A[i][j].type == "L" || A[i][j].type == "R" || A[i][j].type == "D" || A[i][j].type == "U" || A[i][j].type == "CENTER")
						continue;
					r = 0.6 * dx;

					if (abs(A[i][j].y - 0.00415692) < 1e-7)
					{
						if (A[i][j].neibor[1] == i - 1)
							A[i][j].x += rand() / double(RAND_MAX) * random * 2 * r - rand() / double(RAND_MAX) * random * r;
						else
							A[i][j].x += rand() / double(RAND_MAX) * random * r - rand() / double(RAND_MAX) * random * 2 * r;
					}
					else if (A[i][j].type == "SHOCK" || A[i][j].type == "DISCON")
					{
						sumx = sumy = 0;
						r = 0.2 * dx;

						m = A[i][j].neiborsec;
						n = A[i][j].neiborsec_ad;
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							sumx += A[i][A[i][j].neibor[k]].x;
							sumy += A[i][A[i][j].neibor[k]].y;
						}
						for (k = 0; k < A[m][n].neibor.size(); k++)
						{
							sumx += A[m][A[m][n].neibor[k]].x;
							sumy += A[m][A[m][n].neibor[k]].y;
						}
						A[m][n].x = A[i][j].x = sumx / (A[i][j].neibor.size() + A[m][n].neibor.size())
							+ rand() / double(RAND_MAX) * random * 6 * r - rand() / double(RAND_MAX) * 6 * random * r;
						A[m][n].y = A[i][j].y = sumy / (A[i][j].neibor.size() + A[m][n].neibor.size())
							/*+ rand() / double(RAND_MAX) * random * r - rand() / double(RAND_MAX) * random * r*/;
					}
					else
					{
						if (A[i][j].neibor[1] == i - 1)
							A[i][j].x += rand() / double(RAND_MAX) * random * 2 * r - rand() / double(RAND_MAX) * random * r;
						else
							A[i][j].x += rand() / double(RAND_MAX) * random * r - rand() / double(RAND_MAX) * random * 2 * r;
						A[i][j].y += rand() / double(RAND_MAX) * random * 2 * r - rand() / double(RAND_MAX) * random * 2 * r;
					}
				}
			}
		}
		else if (FlowType == "intersection")
		{
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					if (A[i][j].x == xL || A[i][j].x == xR || A[i][j].type == "L" || A[i][j].type == "R" || A[i][j].type == "D" || A[i][j].type == "U" || A[i][j].type == "CENTER")
						continue;
					r = 0.5 * dx;
					//if (abs(A[i][j].x - 0.01685) <= dx / 2)
					//{
					//	A[i][j].x = 0.01685;
					//	//A[i][j].y += rand() / double(RAND_MAX) * random * 2 * r - rand() / double(RAND_MAX) * random * 2 * r;
					//}
					//if (A[i][j].y >= A[i][0].y)
					//{
					//	continue;
					//}
					//else if (abs(A[i][j].x - 0.01685) <= dx / 2)
					//{
					//	A[i][j].x = 0.01685;
					//	//A[i][j].y += rand() / double(RAND_MAX) * random * 1 * r - rand() / double(RAND_MAX) * random * 1 * r;
					//}

					/*else*/ if (A[i][j].type == "SHOCK" || A[i][j].type == "DISCON")
					{
						sumx = sumy = 0;
						r = 0.2 * dx;

						m = A[i][j].neiborsec;
						n = A[i][j].neiborsec_ad;
						for (k = 0; k < A[i][j].neibor.size(); k++)
						{
							sumx += A[i][A[i][j].neibor[k]].x;
							sumy += A[i][A[i][j].neibor[k]].y;
						}
						for (k = 0; k < A[m][n].neibor.size(); k++)
						{
							sumx += A[m][A[m][n].neibor[k]].x;
							sumy += A[m][A[m][n].neibor[k]].y;
						}
						A[m][n].x = A[i][j].x = sumx / (A[i][j].neibor.size() + A[m][n].neibor.size())
							+ rand() / double(RAND_MAX) * random * 3 * r - rand() / double(RAND_MAX) * 3 * random * r;
						A[m][n].y = A[i][j].y = sumy / (A[i][j].neibor.size() + A[m][n].neibor.size())
							+ rand() / double(RAND_MAX) * random * 1.5 * r - rand() / double(RAND_MAX) * 1.5 * random * r;
					}
					else
					{
						if (A[i][j].neibor[1] == i - 1)
							A[i][j].x += rand() / double(RAND_MAX) * random * 2 * r - rand() / double(RAND_MAX) * random * r;
						else
							A[i][j].x += rand() / double(RAND_MAX) * random * r - rand() / double(RAND_MAX) * random * 2 * r;
						A[i][j].y += rand() / double(RAND_MAX) * random * 2 * r - rand() / double(RAND_MAX) * random * 2 * r;
					}
				}
			}
		}

	}
	else if (meshType == 12)
		brick_mesh();



}
void init_polygon_mesh()
{

	extern vector <mesh> A0;
	extern vector <polygon_mesh> M0;
	polygon_mesh t;
	int i, j;
	i = 0;
	int k = 0;
	if (meshType == 10 || meshType == 11 || meshType == 12)
	{
		while (i < Pnum)
		{

			if (i + 1 + 2 * Xnum > Pnum - 1)
				break;

			k++;
			for (int j = 0; j < Xnum; j++)
			{
				if (abs(abs(A0[i].x - A0[i + 1].x) - dx) > 1e-10)
				{
					i++;
					continue;
				}
				if ((i + k - 1) % 2 != 0)
				{
					i++;
					continue;
				}
				if (i + 1 + 2 * Xnum > Pnum - 1)
					continue;
				M0.push_back(t);
				int s = int(M0.size()) - 1;
				M0[s].node.push_back(i);
				M0[s].node.push_back(i + 1);
				M0[s].node.push_back(i + 1 + Xnum);
				M0[s].node.push_back(i + 1 + 2 * Xnum);
				M0[s].node.push_back(i + 2 * Xnum);
				M0[s].node.push_back(i + Xnum);

				M0[s].face_start.push_back(i);
				M0[s].face_start.push_back(i + 1);
				M0[s].face_start.push_back(i + 1 + Xnum);
				M0[s].face_start.push_back(i + 1 + 2 * Xnum);
				M0[s].face_start.push_back(i + 2 * Xnum);
				M0[s].face_start.push_back(i + Xnum);

				M0[s].face_end.push_back(i + 1);
				M0[s].face_end.push_back(i + 1 + Xnum);
				M0[s].face_end.push_back(i + 1 + 2 * Xnum);
				M0[s].face_end.push_back(i + 2 * Xnum);
				M0[s].face_end.push_back(i + Xnum);
				M0[s].face_end.push_back(i);
				i++;
			}
		}
	}
	else
	{
		extern double xL, xR, yU, yD;
		int s;
		while (i < Pnum)
		{
			if (abs(A0[i].x - xR) < 1e-10 || abs(A0[i].y - yU) < 1e-10)
			{
				i++;
				continue;
			}
			M0.push_back(t);
			s = int(M0.size()) - 1;
			M0[s].node.push_back(i);
			M0[s].node.push_back(i + 1);
			M0[s].node.push_back(i + 1 + Xnum);
			M0[s].node.push_back(i + Xnum);

			M0[s].face_start.push_back(i);
			M0[s].face_start.push_back(i + 1);
			M0[s].face_start.push_back(i + 1 + Xnum);
			M0[s].face_start.push_back(i + Xnum);

			M0[s].face_end.push_back(i + 1);
			M0[s].face_end.push_back(i + 1 + Xnum);
			M0[s].face_end.push_back(i + Xnum);
			M0[s].face_end.push_back(i);
			i++;
		}

	}
}
void getNeibor()
{
	extern vector<vector<int>> ad;
	extern vector<vector <mesh>> A;
	extern vector<vector <polygon_mesh>> M;

	int i, j, k, m, s;
	//if (methodType == "Fitting")
	//{
		//���������ڲ�ͬ�ķ�������Ը÷�����ĳһ�㣬����˷����Ķ��������ڵ��¼�˸õ꣬��˶��������������õ�Ϊ���ڽڵ�
	for (i = 0; i < M.size(); i++)
	{
		for (j = 0; j < M[i].size(); j++)
		{
			for (k = 0; k < M[i][j].face_start.size(); k++)
			{
				int n1, n2, s;
				n1 = M[i][j].face_start[k];
				n2 = M[i][j].face_end[k];
				s = 0;
				for (m = 0; m < A[i][n1].neibor.size(); m++)
				{
					if (n2 == A[i][n1].neibor[m])
						s++;
				}
				if (s == 0)
					A[i][n1].neibor.push_back(n2);
				s = 0;
				for (m = 0; m < A[i][n2].neibor.size(); m++)
				{
					if (n1 == A[i][n2].neibor[m])
						s++;
				}
				if (s == 0)
					A[i][n2].neibor.push_back(n1);
			}
		}
	}



	//for (i = 0; i < A.size(); i++)//��ͬ����
	//{
	//	for (j = 0; j < A[i].size(); j++)//ÿ�������ڲ�ͬ�ڵ�
	//	{
	//		for (k = 0; k < M[i].size(); k++)//��ͬ�����Ķ��������
	//		{
	//			for (m = 0; m < M[i][k].node.size(); m++)//ÿ�����������ڵ���Ϣ
	//			{
	//				if (j == M[i][k].node[m])
	//				{
	//					int n1, n2, n11 = -1, n22 = -1;
	//					if (m == 0)
	//					{
	//						n1 = M[i][k].node[M[i][k].node.size() - 1];
	//						n2 = M[i][k].node[1];
	//					}
	//					else if (m == M[i][k].node.size() - 1)
	//					{
	//						n1 = M[i][k].node[m - 1];
	//						n2 = M[i][k].node[0];
	//					}
	//					else
	//					{
	//						n1 = M[i][k].node[m - 1];
	//						n2 = M[i][k].node[m + 1];
	//					}
	//					for (s = 0; s < A[i][j].neibor.size(); s++)
	//					{
	//						if (n1 == A[i][j].neibor[s])
	//							n11 = 1;
	//						if (n2 == A[i][j].neibor[s])
	//							n22 = 1;
	//					}
	//					if (n11 == -1)
	//						A[i][j].neibor.push_back(n1);
	//					if (n22 == -1)
	//						A[i][j].neibor.push_back(n2);
	//				}
	//			}
	//		}

	//	}
	//}



	//for (i = 0; i < A.size(); i++)//��ͬ����
	//{
	//	for (j = 0; j < A[i].size(); j++)//ÿ�������ڲ�ͬ�ڵ�
	//	{
	//		if (A[i][j].type != "IN")
	//			continue;
	//		s = 0;
	//		for (k = 0; k < A[i][j].neibor.size(); k++)
	//		{
	//			m = A[i][j].neibor[k];
	//			if (abs(A[i][m].x - A[i][j + 1].x) < 1e-6 && abs(A[i][m].y - A[i][j + 1].y) < 1e-6)
	//				s++;
	//		}
	//		if (s == 0)
	//			A[i][j].neibor.push_back(j + 1);
	//		else
	//			A[i][j].neibor.push_back(j - 1);
	//	}
	//}

	//}
	//else if (methodType == "Capturing")
	//{
	//	if (FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal")
	//	{
	//		int temp[4];
	//		i = 0;
	//		while (i < A[0].size())
	//		{
	//			for (j = 0; j < Xnum; j++)
	//			{
	//				if (i < Xnum)
	//				{
	//					double d1 = distance(A[0][i + 1], A[0][i]);
	//					if (abs(d1 - dx) < 1e-6)
	//						A[0][i].neibor.push_back(i + 1);
	//					else
	//						A[0][i].neibor.push_back(i - 1);
	//					A[0][i].neibor.push_back(i + Xnum);
	//				}
	//				else if (i >= A[0].size() - Xnum)
	//				{
	//					A[0][i].neibor.push_back(i - Xnum);
	//					double d1 = distance(A[0][i - 1], A[0][i]);
	//					if (abs(d1 - dx) < 1e-6)
	//						A[0][i].neibor.push_back(i - 1);
	//					else
	//						A[0][i].neibor.push_back(i + 1);
	//				}
	//				else
	//				{
	//					A[0][i].neibor.push_back(i - Xnum);
	//					double d1 = distance(A[0][i - 1], A[0][i]);
	//					if (abs(d1 - dx) < 1e-6)
	//						A[0][i].neibor.push_back(i - 1);
	//					else
	//						A[0][i].neibor.push_back(i + 1);
	//					A[0][i].neibor.push_back(i + Xnum);
	//				}

	//				i++;
	//			}
	//		}


	//	}
	//	int m, n;
	//	if (FlowType == "Prandtl-Meyer")
	//	{
	//		for (i = 0; i < A[0].size(); i++)
	//		{
	//			for (j = 0; j < ad.size(); j++)
	//			{
	//				for (k = 0; k < ad[j].size(); k++)
	//				{
	//					if (i != ad[j][k])
	//						continue;
	//					if (i == 109)
	//						i = 109;
	//					if (k == 0)
	//					{
	//						n = 0;
	//						for (int m = 0; m < A[0][i].neibor.size(); m++)
	//						{
	//							if (ad[j][k + 1] == A[0][i].neibor[m])
	//								n++;
	//						}
	//						if (n == 0)
	//							A[0][i].neibor.push_back(ad[j][k + 1]);
	//						n = 0;
	//						for (int m = 0; m < A[0][i].neibor.size(); m++)
	//						{
	//							if (ad[j][ad[j].size() - 1] == A[0][i].neibor[m])
	//								n++;
	//						}
	//						if (n == 0)

	//							A[0][i].neibor.push_back(ad[j][ad[j].size() - 1]);
	//					}
	//					else if (k == ad[j].size() - 1)
	//					{
	//						n = 0;
	//						for (int m = 0; m < A[0][i].neibor.size(); m++)
	//						{
	//							if (ad[j][0] == A[0][i].neibor[m])
	//								n++;
	//						}
	//						if (n == 0)

	//							A[0][i].neibor.push_back(ad[j][0]);
	//						n = 0;
	//						for (int m = 0; m < A[0][i].neibor.size(); m++)
	//						{
	//							if (ad[j][k - 1] == A[0][i].neibor[m])
	//								n++;
	//						}
	//						if (n == 0)

	//							A[0][i].neibor.push_back(ad[j][k - 1]);
	//					}
	//					else
	//					{
	//						n = 0;
	//						for (int m = 0; m < A[0][i].neibor.size(); m++)
	//						{
	//							if (ad[j][k - 1] == A[0][i].neibor[m])
	//								n++;
	//						}
	//						if (n == 0)

	//							A[0][i].neibor.push_back(ad[j][k - 1]);
	//						n = 0;
	//						for (int m = 0; m < A[0][i].neibor.size(); m++)
	//						{
	//							if (ad[j][k + 1] == A[0][i].neibor[m])
	//								n++;
	//						}
	//						if (n == 0)

	//							A[0][i].neibor.push_back(ad[j][k + 1]);
	//					}

	//				}
	//			}
	//		}
	//		for (i = 0; i < A[0].size(); i++)
	//		{
	//			if (A[0][i].neibor.size() > 4)
	//			{
	//				double S = 0, Smax = 0;
	//				for (int n = 0; n < A[0][i].neibor.size(); n++)
	//				{
	//					for (int j = n + 1; j < A[0][i].neibor.size(); j++)
	//					{
	//						for (int k = j + 1; k < A[0][i].neibor.size(); k++)
	//						{
	//							for (int m = k + 1; m < A[0][i].neibor.size(); m++)
	//							{
	//								S = max(Smax, area(A[0][A[0][i].neibor[n]], A[0][A[0][i].neibor[j]], A[0][A[0][i].neibor[k]], A[0][A[0][i].neibor[m]]));
	//								if (S != Smax)
	//								{
	//									A[0][i].neibor1[0] = A[0][i].neibor[n];
	//									A[0][i].neibor1[1] = A[0][i].neibor[j];
	//									A[0][i].neibor1[2] = A[0][i].neibor[k];
	//									A[0][i].neibor1[3] = A[0][i].neibor[m];
	//									Smax = S;
	//								}
	//							}

	//						}

	//					}

	//				}
	//			}

	//		}

	//	}
	//}

}
void getType()
{
	extern vector <mesh> A0;
	extern double xL, xR, yU, yD;


	if (FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal")
	{
		int i;
		for (i = 0; i < Pnum; i++)
		{
			if (abs(A0[i].x - xL) < 1e-10)
			{
				//if (A[i].y == yD)
				//	A[i].type = "LD";
				//else if (A[i].y == yU)
				//	A[i].type = "LU";
				//else
				A0[i].type = "L";
			}
			else if (abs(A0[i].x - xR) < 1e-10)
			{
				//if (A[i].y == yD)
				//	A[i].type = "RD";
				//else if (A[i].y == yU)
				//	A[i].type = "RU";
				//else
				A0[i].type = "R";
			}
			else if (abs(A0[i].y - yD) < 1e-10)
				A0[i].type = "D";
			else if (abs(A0[i].y - yU) < 1e-10)
				A0[i].type = "U";
			else
				A0[i].type = "IN";
		}

	}
	if (FlowType == "Prandtl-Meyer")
	{
		int i;
		for (i = 0; i < Pnum; i++)
		{
			if (A0[i].x == 0)
				A0[i].type = "L";
			else if (A0[i].y == 0)
				A0[i].type = "DL";//�±߽������
			else
				A0[i].type = "IN";

		}
		for (i = Pnum; i < 2 * Pnum - 1; i++)
		{
			if (i < Pnum + Xnum - 1)
				A0[i].type = "DR";//�±߽���Ұ��
			else if ((i - (Pnum + Xnum - 2)) % Xnum == 0)
				A0[i].type = "R";
			else
				A0[i].type = "IN";
		}
		for (i = 2 * Pnum - 1; i < A0.size(); i++)
		{
			A0[i].type = "IN";
		}
		A0[Xnum - 1].type = "IN";
	}
}

void initFlow()
{
	if (methodType == "F")
	{
		int i, j;
		extern vector<vector <mesh>> A;
		if (FlowType == "normal")
		{
			using namespace Normal;
			mesh A1, A2;
			A1.�� = ��1, A1.p = p1, A1.u =u1, A1.v = v1;
			get_down(A1, A2, ConstPara::pi / 2);
			��2 = A2.��, p2 = A2.p, u2 = A2.u, v2 = A2.v;
			//u2 += 1;
			//u1 += 1;
			��1 = A2.��, p1 = A2.p, u1 = A2.u, v1 = A2.v;
			��2 = A1.��, p2 = A1.p, u2 = A1.u, v2 = A1.v;
			u1 = u2 - u1;
			u2 = 0;
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					//��1 = 3.8571, p1 = 10.3333, u1 = 2.6293, v1 = 0;
					//��2 = 1 + 0.2 * sin(5 * A[i][j].x), p2 = 1, u2 = 0, v2 = 0;
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

				}

			}
		}
		if (FlowType == "oblique")
		{
			using namespace Oblique;
			mesh A1, A2;
			A1.�� = ��1, A1.p = p1, A1.u = u1, A1.v = v1;
			get_down(A1, A2, ��);
			��2 = A2.��, p2 = A2.p, u2 = A2.u, v2 = A2.v;

			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
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

				}

			}
		}
		else if (FlowType == "intersection")
		{
			using namespace ShockwaveCross;
			using namespace ConstPara;
			double ��40 = 45 * pi / 180, ��41 = 50 * pi / 180, ��42 = 60 * pi / 180;
			double p41, p40;
			double p51, p50;
			double ��41, ��40;
			double ��51, ��50;
			double ��51, ��50;
			double f��41, f��40;
			double un20, un21, Mu20, Mu21;
			double un30, un31, Mu30, Mu31;
			mesh A1, A2, A3;
			A1.�� = ��1, A1.p = p1, A1.u = u1, A1.v = v1;
			get_down(A1, A2, ��2/*-0.05*pi/180*/ );
			get_down(A1, A3, ��3);

			��2 = A2.��, p2 = A2.p, u2 = A2.u, v2 = A2.v;
			��3 = A3.��, p3 = A3.p, u3 = A3.u, v3 = A3.v;

			//A2.�� = ��2, A2.p = p2, A2.u = u2, A2.v = v2;
			//A3.�� = ��3, A3.p = p3, A3.u = u3, A3.v = v3;
			mesh A40, A41, A50, A51;
			while (abs(��42 - ��41) > 1e-15)
			{
				while (��42 > -��2 || ��42 < ��2)
				{
					if (abs(��42) > 1e5)
					{
						��42 = ��42 / 1e5;
						continue;
					}
					if (��42 > -��2)
						��42 = ��42 + ��2;
					else if (��42 < ��2)
						��42 = ��42 - ��2;
					//if (��42 > pi / 2)
					//	��42 = ��42 - pi / 2;
					//	else if (��42 < -pi / 2)
					//	��42 = ��42 + pi / 2;
				}
				if (abs(��42) < 1e-7)
					��42 = 30 * pi / 180;

				��40 = (��41 + ��42) / 2;
				��41 = ��42;

				get_down(A2, A40, ��40);
				A50.p = A40.p;
				��50 = get_��(A3, A50.p, -1);
				get_down(A3, A50, ��50);

				get_down(A2, A41, ��41);
				A51.p = A41.p;
				��51 = get_��(A3, A51.p, -1);
				get_down(A3, A51, ��51);

				��40 = get_��(A40.u, A40.v);
				��41 = get_��(A41.u, A41.v);
				��50 = get_��(A50.u, A50.v);
				��51 = get_��(A51.u, A51.v);
				f��41 = ��41 - ��51;
				f��40 = ��40 - ��50;
				if (abs(f��41 - f��40) <= 1e-20)
				{
					��42 = ��41;
					break;
				}
				��42 = ��41 - f��41 / (f��41 - f��40) * (��41 - ��40);
			}
			��4 = ��42;
			get_down(A2, A41, ��4);
			A51.p = A41.p;
			��5 = get_��(A3, A51.p, -1);
			get_down(A3, A51, ��5);
			��4 = get_��(A41.u, A41.v);
			��5 = get_��(A51.u, A51.v);

			double un2 = u2 * sin(��4) - v2 * cos(��4);
			double ut2 = u2 * cos(��4) + v2 * sin(��4);
			double Mu2 = get_Ma(un2, 0, ��2, p2);
			��4 = ��2 * get_��2��1(Mu2);
			p4 = p2 * get_p2p1(Mu2);
			double Md4 = sqrt((Mu2 * Mu2 + 2 / (�� - 1)) / (2 * �� * Mu2 * Mu2 / (�� - 1) - 1));
			double c4 = sqrt(�� * p4 / ��4);
			double un4 = Md4 * c4;
			double ut4 = ut2;
			u4 = un4 * sin(��4) + ut4 * cos(��4);
			v4 = -un4 * cos(��4) + ut4 * sin(��4);


			double un3 = -u3 * sin(��5) + v3 * cos(��5);
			double ut3 = u3 * cos(��5) + v3 * sin(��5);
			double Mu3 = get_Ma(un3, 0, ��3, p3);
			��5 = ��3 * get_��2��1(Mu3);
			p5 = p3 * get_p2p1(Mu3);
			double Md5 = sqrt((Mu3 * Mu3 + 2 / (�� - 1)) / (2 * �� * Mu3 * Mu3 / (�� - 1) - 1));
			double c5 = sqrt(�� * p5 / ��5);
			double un5 = Md5 * c5;
			double ut5 = ut3;
			u5 = un5 * sin(-��5) + ut5 * cos(-��5);
			v5 = un5 * cos(-��5) - ut5 * sin(-��5);
			std::cout.precision(20);
			std::cout << "p2= " << p2 << "   p3= " << p3 << std::endl;
			std::cout << "u2= " << u2 << "   u3= " << u3 << std::endl;
			std::cout << "v2= " << v2 << "   v3= " << v3 << std::endl;
			std::cout << "rho2= " << ��2 << "   rho3= " << ��3 << std::endl;
			std::cout << "��4= " << ��4 * 180 / pi << "   ��5= " << ��5 * 180 / pi << std::endl;
			std::cout << "��4= " << ��4 * 180 / pi << "   ��5= " << ��5 * 180 / pi << std::endl;
			std::cout << "p4= " << p4 << "   p5= " << p5 << std::endl;
			std::cout << "u4= " << u4 << "   u5= " << u5 << std::endl;
			std::cout << "v4= " << v4 << "   v5= " << v5 << std::endl;
			std::cout << "rho4= " << ��4 << "   rho5= " << ��5 << std::endl;
			std::cout << std::endl;
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					A[i][j].um = A[i][j].vm = 0;
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
						//if (A[i][j].type == "L")
						//{
						A[i][j].�� = ��3;
						A[i][j].u = u3;
						A[i][j].v = v3;
						A[i][j].p = p3;

						//}
						//else
						//{
						//	A[i][j].�� = ��3;
						//	A[i][j].u = u3 - 5;
						//	A[i][j].v = v3 -0.5;
						//	A[i][j].p = p3;
						//}
					}
					else if (i == 3)
					{
						A[i][j].�� = ��4;
						A[i][j].u = u4;
						A[i][j].v = v4;
						A[i][j].p = p4;
					}
					else if (i == 4)
					{
						A[i][j].�� = ��5;
						A[i][j].u = u5;
						A[i][j].v = v5;
						A[i][j].p = p5;
					}
					else
					{
						std::cout << "something wrong in intersection!" << std::endl;
					}

				}

			}

		}
	}
	else if (methodType == "C")
	{
		if (FlowType == "normal")
			init_flow_normal();
		else if (FlowType == "oblique")
			init_flow_shockwave();
		else if (FlowType == "intersection")
			init_flow_shockwaveCross();
	}
}
void init_flow_normal()
{
	extern vector<vector <mesh>> A;
	using namespace Normal;
	int i, j;
	mesh A1, A2;
	A1.�� = ��1, A1.p = p1, A1.u = u1, A1.v = v1;
	get_down(A1, A2, ConstPara::pi / 2);
	��2 = A2.��, p2 = A2.p, u2 = A2.u, v2 = A2.v;
	std::cout.precision(20);
	std::cout << "p1= " << p1 << "   p2= " << p2 << std::endl;
	std::cout << "u1= " << u1 << "   u2= " << u2 << std::endl;
	std::cout << "v1= " << v1 << "   v2= " << v2 << std::endl;
	std::cout << "rho1= " << ��1 << "   rho2= " << ��2 << std::endl;
	//u2 += 1;
	//u1 += 1;
	��1 = A2.��, p1 = A2.p, u1 = A2.u, v1 = A2.v;
	��2 = A1.��, p2 = A1.p, u2 = A1.u, v2 = A1.v;
	u1 = u2 - u1;
	u2 = 0;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			��1 = 3.8571, p1 = 10.3333, u1 = 2.6293, v1 = 0;
			��2 = 1 + 0.2 /** sin(5 * A[i][j].x)*/, p2 = 1, u2 = 0, v2 = 0;
			//��2 = 3 + 0.5 * sin(5 * A[i][j].x);

			if (A[i][j].x < 1)
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
	}

}
void init_flow_uniform()//������
{
	extern vector<vector <mesh>> A;
	int i;
	using namespace Init;
	for (i = 0; i < A[0].size(); i++)
	{
		A[0][i].�� = ��0;
		A[0][i].u = u0;
		A[0][i].v = v0;
		A[0][i].p = p0;
	}
}
void init_flow_shockwave()//б����
{
	extern vector<vector <mesh>> A;
	using namespace Oblique;
	mesh A1, A2;
	A1.�� = ��1, A1.p = p1, A1.u = u1, A1.v = v1;
	get_down(A1, A2, ��);
	��2 = A2.��, p2 = A2.p, u2 = A2.u, v2 = A2.v;
	using std::cout;
	using std::endl;
	cout.precision(10);
	cout << "��1=" << ��1 << ", u1= " << u1 << ", v1=" << v1 << ",  p1=" << p1 << endl;
	cout << "��2=" << ��2 << ", u2= " << u2 << ", v2=" << v2 << ",  p2=" << p2 << endl;

	int i, j;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j].x < 0.002)
			{
				A[i][j].�� = ��2;
				A[i][j].u = u2;
				A[i][j].v = v2;
				A[i][j].p = p2;
			}
			else
			{
				A[i][j].�� = ��2;
				A[i][j].u = u2;
				A[i][j].v = v2;
				A[i][j].p = p2;
			}
		}
		//if (A[0][i].type == "IN")
		//{
		//	A[0][i].�� = ��1;
		//	A[0][i].u = u1;
		//	A[0][i].v = v1;
		//	A[0][i].p = p1;
		//}
		//else if (A[0][i].type == "L")
		//{
		//	A[0][i].�� = ��1;
		//	A[0][i].u = u1;
		//	A[0][i].v = v1;
		//	A[0][i].p = p1;
		//}
		//else if (A[0][i].type == "U")
		//{
		//	if (A[0][i].x < A[0][startpoint].x)
		//	{
		//		A[0][i].�� = ��1;
		//		A[0][i].u = u1;
		//		A[0][i].v = v1;
		//		A[0][i].p = p1;
		//	}
		//	else
		//	{
		//		A[0][i].�� = ��2;
		//		A[0][i].u = u2;
		//		A[0][i].v = v2;
		//		A[0][i].p = p2;
		//	}
		//}
		//else if (A[0][i].type == "D" || A[0][i].type == "R")
		//{
		//	if (A[0][i].x < A[0][startpoint].x + dy * Ynum / tan(��))
		//	{
		//		A[0][i].�� = ��1;
		//		A[0][i].u = u1;
		//		A[0][i].v = v1;
		//		A[0][i].p = p1;
		//	}
		//	else
		//	{
		//		A[0][i].�� = ��2;
		//		A[0][i].u = u2;
		//		A[0][i].v = v2;
		//		A[0][i].p = p2;
		//	}
		//}
		//else
		//	std::cout << "something wrong!" << std::endl;
	}
}
void init_flow_shockwaveCross()//ͬ�༤���ཻ
{
	extern vector<vector <mesh>> A;
	extern vector <mesh> A0;
	double Ma1;
	using namespace ShockwaveCross;
	using namespace Init;
	int i;
	using namespace ConstPara;
	Line L12;
	Line L13;
	mesh C;

	mesh A1, A2, A3;
	A1.�� = ��1, A1.p = p1, A1.u = u1, A1.v = v1;
	get_down(A1, A2, ��2);
	get_down(A1, A3, ��3);
	��2 = A2.��, p2 = A2.p, u2 = A2.u, v2 = A2.v;
	��3 = A3.��, p3 = A3.p, u3 = A3.u, v3 = A3.v;
	double ��12 = -30 / 180.0 * ConstPara::pi;
	double ��13 = 30 / 180.0 * ConstPara::pi;
	mesh A12;
	mesh A13;
	vector<mesh> t;
	A12.x = -dx / 2, A12.y = A0[Pnum - 1].y - 2 * dy * 5 * Ynum / 43;
	A13.x = -dx / 2, A13.y = 2 * dy * 5 * Ynum / 43;
	L12 = getLine(��12, A12);
	L13 = getLine(��13, A13);
	C = getCrossPoint(L12, L13);
	L12 = getLine(��2, C);
	L13 = getLine(��3, C);
	for (i = 0; i < Pnum; i++)
	{
		if (A[0][i].type == "L")
		{
			if (A[0][i].x * L12.A + A[0][i].y * L12.B + L12.C < 0)
			{
				A[0][i].�� = ��2;
				A[0][i].u = u2;
				A[0][i].v = v2;
				A[0][i].p = p2;
			}
			else if (A[0][i].x * L13.A + A[0][i].y * L13.B + L13.C > 0)
			{
				A[0][i].�� = ��3;
				A[0][i].u = u3;
				A[0][i].v = v3;
				A[0][i].p = p3;
			}
			else
			{
				A[0][i].�� = ��1;
				A[0][i].u = u1;
				A[0][i].v = v1;
				A[0][i].p = p1;
			}
		}
		else
		{
			A[0][i].�� = ��0;
			A[0][i].u = u0;
			A[0][i].v = v0;
			A[0][i].p = p0;
		}
	}
}
//void init_U()
//{
//	extern vector<vector <mesh>> A;
//	extern vector<vector<vector <double>>> U;
//	vector<vector <double>> u0;
//	vector <double> u1;
//	int i, j;
//	for (i = 0; i < A.size(); i++)
//	{
//
//		U.push_back(u0);
//		for (j = 0; j < A[i].size(); j++)
//		{
//			U[i].push_back(u1);
//			U[i][j].push_back(A[i][j].��);
//			U[i][j].push_back(A[i][j].�� * A[i][j].u);
//			U[i][j].push_back(A[i][j].�� * A[i][j].v);
//			U[i][j].push_back(0.5 * A[i][j].�� * (A[i][j].u * A[i][j].u + A[i][j].v * A[i][j].v) + A[i][j].p / (ConstPara::�� - 1));
//		}
//	}
//}
void findConnectPoint()
{
	extern vector<vector <mesh>> A;
	int i, j, k;
	if (FlowType == "oblique" || FlowType == "normal")
	{
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "SHOCK")
					continue;
				for (k = 0; k < A[i].size(); k++)
				{
					if (A[i][k].type != "SHOCK")
						continue;
					if (abs(A[i][j].y - A[i][k].y) <= dy + 1e-10)
						A[i][j].moveConnct.push_back(k);
				}

			}
		}
	}
	else if (FlowType == "intersection")
	{

		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j].type == "SHOCK" || A[i][j].type == "DISCON" || A[i][j].type == "CENTER")
					continue;
				for (k = 0; k < A[i].size(); k++)
				{
					if (A[i][k].type == "SHOCK" || A[i][k].type == "DISCON")
					{
						if (abs(A[i][k].x - A[i][j].x) <= dx / 2 + 1e-10)
							A[i][j].moveConnct.push_back(k);
					}

				}

			}
		}



	}
}
//void init_Ar()
//{
//	extern vector <mesh> A0;
//	extern vector <mesh> Ar;
//	for (int i = 0; i < Pnum; i++)
//	{
//		Ar[i].x = A0[i].x;
//		Ar[i].y = A0[i].y;
//		Ar[i].�� = A0[i].��;
//		Ar[i].u = A0[i].u;
//		Ar[i].v = A0[i].v;
//		Ar[i].p = A0[i].p;
//		Ar[i].neibor = A0[i].neibor;
//	}
//}