#include"const.h"
#include"functions.h"
#include<vector>

void coordinate_trans()
{
	extern std::vector<std::vector<mesh>> A;
	int i, j,k;
	int n1, n2, n3, n4;
	extern double t_sim;
	extern double dt;
	extern std::vector < std::vector <mesh>> Ar;

	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{

			if (A[i][j].neibor.size() < 3)//С��3�ĵ��޷�������任
				continue;
			if (t_sim == 0)
			{
				if (A[i][j].neibor.size() == 3)
				{
					using MeshPara::method;
					for (k = 0; k < 12; k++)
					{
						n1 = method[k][0];
						n2 = method[k][1];
						n3 = method[k][2];
						n4 = method[k][3];
						A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[n1]].x - A[i][A[i][j].neibor[n3]].x));
						A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[n1]].y - A[i][A[i][j].neibor[n3]].y));
						A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[n2]].x - A[i][A[i][j].neibor[n4]].x));
						A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[n2]].y - A[i][A[i][j].neibor[n4]].y));
						A[i][j].J.push_back(1 / (A[i][j].x��[k] * A[i][j].y��[k] - A[i][j].x��[k] * A[i][j].y��[k]));
						A[i][j].��x.push_back(A[i][j].y��[k] * A[i][j].J[k]);
						A[i][j].��y.push_back(-A[i][j].x��[k] * A[i][j].J[k]);
						A[i][j].��x.push_back(-A[i][j].y��[k] * A[i][j].J[k]);
						A[i][j].��y.push_back(A[i][j].x��[k] * A[i][j].J[k]);
						A[i][j].x��.push_back(0);
						A[i][j].y��.push_back(0);
						A[i][j].��t.push_back(0);
						A[i][j].��t.push_back(0);
					}
					//A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[0]].x - A[i][A[i][j].neibor[2]].x));
					//A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[0]].y - A[i][A[i][j].neibor[2]].y));
					//A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[1]].x - A[i][A[i][j].neibor[2]].x));
					//A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[1]].y - A[i][A[i][j].neibor[2]].y));
					//A[i][j].J.push_back(1 / (A[i][j].x��[0] * A[i][j].y��[0] - A[i][j].x��[0] * A[i][j].y��[0]));
					//A[i][j].��x.push_back(A[i][j].y��[0] * A[i][j].J[0]);
					//A[i][j].��y.push_back(-A[i][j].x��[0] * A[i][j].J[0]);
					//A[i][j].��x.push_back(-A[i][j].y��[0] * A[i][j].J[0]);
					//A[i][j].��y.push_back(A[i][j].x��[0] * A[i][j].J[0]);
					//A[i][j].x��.push_back(0);
					//A[i][j].y��.push_back(0);
					//A[i][j].��t.push_back(0);
					//A[i][j].��t.push_back(0);

					//A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[1]].x - A[i][A[i][j].neibor[0]].x));
					//A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[1]].y - A[i][A[i][j].neibor[0]].y));
					//A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[2]].x - A[i][A[i][j].neibor[0]].x));
					//A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[2]].y - A[i][A[i][j].neibor[0]].y));
					//A[i][j].J.push_back(1 / (A[i][j].x��[1] * A[i][j].y��[1] - A[i][j].x��[1] * A[i][j].y��[1]));
					//A[i][j].��x.push_back(A[i][j].y��[1] * A[i][j].J[1]);
					//A[i][j].��y.push_back(-A[i][j].x��[1] * A[i][j].J[1]);
					//A[i][j].��x.push_back(-A[i][j].y��[1] * A[i][j].J[1]);
					//A[i][j].��y.push_back(A[i][j].x��[1] * A[i][j].J[1]);
					//A[i][j].x��.push_back(0);
					//A[i][j].y��.push_back(0);
					//A[i][j].��t.push_back(0);
					//A[i][j].��t.push_back(0);

					//A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[2]].x - A[i][A[i][j].neibor[1]].x));
					//A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[2]].y - A[i][A[i][j].neibor[1]].y));
					//A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[0]].x - A[i][A[i][j].neibor[1]].x));
					//A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[0]].y - A[i][A[i][j].neibor[1]].y));
					//A[i][j].J.push_back(1 / (A[i][j].x��[2] * A[i][j].y��[2] - A[i][j].x��[2] * A[i][j].y��[2]));
					//A[i][j].��x.push_back(A[i][j].y��[2] * A[i][j].J[2]);
					//A[i][j].��y.push_back(-A[i][j].x��[2] * A[i][j].J[2]);
					//A[i][j].��x.push_back(-A[i][j].y��[2] * A[i][j].J[2]);
					//A[i][j].��y.push_back(A[i][j].x��[2] * A[i][j].J[2]);
					//A[i][j].x��.push_back(0);
					//A[i][j].y��.push_back(0);
					//A[i][j].��t.push_back(0);
					//A[i][j].��t.push_back(0);
				}
				else if (A[i][j].neibor.size() == 4)
				{
					A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[0]].x - A[i][A[i][j].neibor[2]].x));
					A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[0]].y - A[i][A[i][j].neibor[2]].y));
					A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor[1]].x - A[i][A[i][j].neibor[3]].x));
					A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor[1]].y - A[i][A[i][j].neibor[3]].y));
					A[i][j].J.push_back(1 / (A[i][j].x��[0] * A[i][j].y��[0] - A[i][j].x��[0] * A[i][j].y��[0]));
					A[i][j].��x.push_back(A[i][j].y��[0] * A[i][j].J[0]);
					A[i][j].��y.push_back(-A[i][j].x��[0] * A[i][j].J[0]);
					A[i][j].��x.push_back(-A[i][j].y��[0] * A[i][j].J[0]);
					A[i][j].��y.push_back(A[i][j].x��[0] * A[i][j].J[0]);
					A[i][j].x��.push_back(0);
					A[i][j].y��.push_back(0);
					A[i][j].��t.push_back(0);
					A[i][j].��t.push_back(0);
				}
				else if (A[i][j].neibor.size() > 4)
				{
					A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor1[0]].x - A[i][A[i][j].neibor1[2]].x));
					A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor1[0]].y - A[i][A[i][j].neibor1[2]].y));
					A[i][j].x��.push_back(0.5 * (A[i][A[i][j].neibor1[1]].x - A[i][A[i][j].neibor1[3]].x));
					A[i][j].y��.push_back(0.5 * (A[i][A[i][j].neibor1[1]].y - A[i][A[i][j].neibor1[3]].y));
					A[i][j].J.push_back(1 / (A[i][j].x��[0] * A[i][j].y��[0] - A[i][j].x��[0] * A[i][j].y��[0]));
					A[i][j].��x.push_back(A[i][j].y��[0] * A[i][j].J[0]);
					A[i][j].��y.push_back(-A[i][j].x��[0] * A[i][j].J[0]);
					A[i][j].��x.push_back(-A[i][j].y��[0] * A[i][j].J[0]);
					A[i][j].��y.push_back(A[i][j].x��[0] * A[i][j].J[0]);
					A[i][j].x��.push_back(0);
					A[i][j].y��.push_back(0);
					A[i][j].��t.push_back(0);
					A[i][j].��t.push_back(0);

				}

			}
			else//t_sim��������ʱ����Ҫ����任��ӦΪ������
			{
				if (A[i][j].neibor.size() == 4)
				{
					A[i][j].x��[0] = 0.5 * (A[i][A[i][j].neibor[0]].x - A[i][A[i][j].neibor[2]].x);
					A[i][j].y��[0] = 0.5 * (A[i][A[i][j].neibor[0]].y - A[i][A[i][j].neibor[2]].y);
					A[i][j].x��[0] = 0.5 * (A[i][A[i][j].neibor[1]].x - A[i][A[i][j].neibor[3]].x);
					A[i][j].y��[0] = 0.5 * (A[i][A[i][j].neibor[1]].y - A[i][A[i][j].neibor[3]].y);
					A[i][j].J[0] = 1 / (A[i][j].x��[0] * A[i][j].y��[0] - A[i][j].x��[0] * A[i][j].y��[0]);
					A[i][j].��x[0] = A[i][j].y��[0] * A[i][j].J[0];
					A[i][j].��y[0] = -A[i][j].x��[0] * A[i][j].J[0];
					A[i][j].��x[0] = -A[i][j].y��[0] * A[i][j].J[0];
					A[i][j].��y[0] = A[i][j].x��[0] * A[i][j].J[0];
					A[i][j].x��[0] = (A[i][j].x - Ar[i][j].x) / dt;
					A[i][j].y��[0] = (A[i][j].y - Ar[i][j].y) / dt;
					A[i][j].��t[0] = -A[i][j].x��[0] * A[i][j].��x[0] - A[i][j].y��[0] * A[i][j].��y[0];
					A[i][j].��t[0] = -A[i][j].x��[0] * A[i][j].��x[0] - A[i][j].y��[0] * A[i][j].��y[0];
				}
				if (A[i][j].neibor.size() > 4)
				{
					A[i][j].x��[0] = 0.5 * (A[i][A[i][j].neibor1[0]].x - A[i][A[i][j].neibor1[2]].x);
					A[i][j].y��[0] = 0.5 * (A[i][A[i][j].neibor1[0]].y - A[i][A[i][j].neibor1[2]].y);
					A[i][j].x��[0] = 0.5 * (A[i][A[i][j].neibor1[1]].x - A[i][A[i][j].neibor1[3]].x);
					A[i][j].y��[0] = 0.5 * (A[i][A[i][j].neibor1[1]].y - A[i][A[i][j].neibor1[3]].y);
					A[i][j].J[0] = 1 / (A[i][j].x��[0] * A[i][j].y��[0] - A[i][j].x��[0] * A[i][j].y��[0]);
					A[i][j].��x[0] = A[i][j].y��[0] * A[i][j].J[0];
					A[i][j].��y[0] = -A[i][j].x��[0] * A[i][j].J[0];
					A[i][j].��x[0] = -A[i][j].y��[0] * A[i][j].J[0];
					A[i][j].��y[0] = A[i][j].x��[0] * A[i][j].J[0];
					A[i][j].x��[0] = (A[i][j].x - Ar[i][j].x) / dt;
					A[i][j].y��[0] = (A[i][j].y - Ar[i][j].y) / dt;
					A[i][j].��t[0] = -A[i][j].x��[0] * A[i][j].��x[0] - A[i][j].y��[0] * A[i][j].��y[0];
					A[i][j].��t[0] = -A[i][j].x��[0] * A[i][j].��x[0] - A[i][j].y��[0] * A[i][j].��y[0];

				}
				if (A[i][j].neibor.size() == 3)
				{
					A[i][j].x��[0] = 0.5 * (A[i][A[i][j].neibor[0]].x - A[i][A[i][j].neibor[2]].x);
					A[i][j].y��[0] = 0.5 * (A[i][A[i][j].neibor[0]].y - A[i][A[i][j].neibor[2]].y);
					A[i][j].x��[0] = 0.5 * (A[i][A[i][j].neibor[1]].x - A[i][A[i][j].neibor[2]].x);
					A[i][j].y��[0] = 0.5 * (A[i][A[i][j].neibor[1]].y - A[i][A[i][j].neibor[2]].y);
					A[i][j].J[0] = 1 / (A[i][j].x��[0] * A[i][j].y��[0] - A[i][j].x��[0] * A[i][j].y��[0]);
					A[i][j].��x[0] = A[i][j].y��[0] * A[i][j].J[0];
					A[i][j].��y[0] = -A[i][j].x��[0] * A[i][j].J[0];
					A[i][j].��x[0] = -A[i][j].y��[0] * A[i][j].J[0];
					A[i][j].��y[0] = A[i][j].x��[0] * A[i][j].J[0];
					A[i][j].x��[0] = (A[i][j].x - Ar[i][j].x) / dt;
					A[i][j].y��[0] = (A[i][j].y - Ar[i][j].y) / dt;
					A[i][j].��t[0] = -A[i][j].x��[0] * A[i][j].��x[0] - A[i][j].y��[0] * A[i][j].��y[0];
					A[i][j].��t[0] = -A[i][j].x��[0] * A[i][j].��x[0] - A[i][j].y��[0] * A[i][j].��y[0];

					A[i][j].x��[1] = 0.5 * (A[i][A[i][j].neibor[1]].x - A[i][A[i][j].neibor[0]].x);
					A[i][j].y��[1] = 0.5 * (A[i][A[i][j].neibor[1]].y - A[i][A[i][j].neibor[0]].y);
					A[i][j].x��[1] = 0.5 * (A[i][A[i][j].neibor[2]].x - A[i][A[i][j].neibor[0]].x);
					A[i][j].y��[1] = 0.5 * (A[i][A[i][j].neibor[2]].y - A[i][A[i][j].neibor[0]].y);
					A[i][j].J[1] = 1 / (A[i][j].x��[1] * A[i][j].y��[1] - A[i][j].x��[1] * A[i][j].y��[1]);
					A[i][j].��x[1] = A[i][j].y��[1] * A[i][j].J[1];
					A[i][j].��y[1] = -A[i][j].x��[1] * A[i][j].J[1];
					A[i][j].��x[1] = -A[i][j].y��[1] * A[i][j].J[1];
					A[i][j].��y[1] = A[i][j].x��[1] * A[i][j].J[1];
					A[i][j].x��[1] = (A[i][j].x - Ar[i][j].x) / dt;
					A[i][j].y��[1] = (A[i][j].y - Ar[i][j].y) / dt;
					A[i][j].��t[1] = -A[i][j].x��[1] * A[i][j].��x[1] - A[i][j].y��[1] * A[i][j].��y[1];
					A[i][j].��t[1] = -A[i][j].x��[1] * A[i][j].��x[1] - A[i][j].y��[1] * A[i][j].��y[1];

					A[i][j].x��[2] = 0.5 * (A[i][A[i][j].neibor[2]].x - A[i][A[i][j].neibor[1]].x);
					A[i][j].y��[2] = 0.5 * (A[i][A[i][j].neibor[2]].y - A[i][A[i][j].neibor[1]].y);
					A[i][j].x��[2] = 0.5 * (A[i][A[i][j].neibor[0]].x - A[i][A[i][j].neibor[1]].x);
					A[i][j].y��[2] = 0.5 * (A[i][A[i][j].neibor[0]].y - A[i][A[i][j].neibor[1]].y);
					A[i][j].J[2] = 1 / (A[i][j].x��[2] * A[i][j].y��[2] - A[i][j].x��[2] * A[i][j].y��[2]);
					A[i][j].��x[2] = A[i][j].y��[2] * A[i][j].J[2];
					A[i][j].��y[2] = -A[i][j].x��[2] * A[i][j].J[2];
					A[i][j].��x[2] = -A[i][j].y��[2] * A[i][j].J[2];
					A[i][j].��y[2] = A[i][j].x��[2] * A[i][j].J[2];
					A[i][j].x��[2] = (A[i][j].x - Ar[i][j].x) / dt;
					A[i][j].y��[2] = (A[i][j].y - Ar[i][j].y) / dt;
					A[i][j].��t[2] = -A[i][j].x��[2] * A[i][j].��x[2] - A[i][j].y��[2] * A[i][j].��y[2];
					A[i][j].��t[2] = -A[i][j].x��[2] * A[i][j].��x[2] - A[i][j].y��[2] * A[i][j].��y[2];

				}

			}
		}
	}
}
