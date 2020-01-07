#pragma once
#include<cmath>
#include<string>
#include<vector>
using std::vector;
using std::string;
namespace ConstPara//����
{
	const double CFL = 0.5;
	const double t_end = 100;//����ʱ��
	const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628;
	const double �� = 1.4;//���ȱ�
	const double random = 0.495;
}
namespace Init//��ʼ����
{
	extern double ��0;
	extern double v0;
	extern double p0;
	extern double u0;
}
namespace Normal//��ʼ����
{
	extern double ��1;
	extern double v1;
	extern double p1;
	extern double u1;
	extern double ��2;
	extern double v2;
	extern double p2;
	extern double u2;
}

namespace MeshPara
{
	const string FlowType = "normal";
	//Ŀǰ�У����Ͳ�Prandtl-Meyer;������normal��б����oblique �������ཻintersection,����other
	const int meshType = 20;
	//10������������11�Ŷ�������������12ש������
	//20�����νṹ����
	const string methodType = "C";
	//���㷽����������׽C-Capturing������װ��F-Fitting
	const int method[12][4] = {
{1,2,0,0},//����3
{0,1,2,2},
{2,0,1,1},

{0,1,1,2},//����2
{1,2,2,0},
{2,0,0,1},


{0,0,1,2},//����1
{1,1,2,0},
{2,2,0,1},

{0,1,2,0},//����4
{2,0,1,2},
{1,2,0,1},

//{2,1,1,0},//����22
//{0,2,2,1},
//{1,0,0,2},
//
//{2,1,0,0},//����32
//{1,0,2,2},
//{0,2,1,1},
//
//{0,0,2,1},//����12
//{1,1,0,2},
//{2,2,1,0},
//
//{0,2,1,0},//����42
//{2,1,0,2},
//{1,0,2,1},

	};

	const int Xnum = 1000;//x�����ĸ���20,40
	const int Ynum = 25;//y�����ĸ���43,87
	const int Pnum = Xnum * Ynum;//һ����ĸ���
	const double dx = 0.0167/2.5;
	const double dy = sqrt(3.0) / 2 * dx;
}

struct mesh
{
	double x;
	double y;
	double ��;//�ܶ�
	double u;//ˮƽ�ٶ�
	double v;//��ֱ�ٶ�
	double p;//ѹǿ
	double x0;//������ĳ�ʼ����
	double y0;
	double um;//����˶��ٶ�
	double vm;
	int section;//����
	int step = 0;//�����õ��ںβ�����µ�
	int neiborsec = -1;//���ڷ����������ֵΪ�������ʾ�õ㲻�Ǽ����㣬��֮Ϊ���ڷ������õļ�����
	int neiborsec_ad = -1;//���ڷ����ĵ�ַ
	vector <double> ��x;
	vector <double> ��y;
	vector <double> ��t;
	vector <double> ��x;
	vector <double> ��y;
	vector <double> ��t;
	vector <double> x��;
	vector <double> x��;
	vector <double> x��;
	vector <double> y��;
	vector <double> y��;
	vector <double> y��;
	vector <double> J;//�ſɱ�����ʽ
	vector <int>neibor;//��¼�ø������ڸ��λ����Ϣ
	vector <int>moveConnct;//�˶������㣬�������˶����moveConnect���й�
	int neibor1[4] = { 0 };//�����ڸ�㳬���ĸ���ɸѡ�ĸ���˴�������任�ͼ������õ������ڸ�㶼�ô˴�
	string type;//������ͣ���Ϊ�������ұ߽��Լ��ڲ���U��D,L,R,IN,������SHOCK���Ӵ���ϵ�DISCON��contact discontinuity��,�����ཻ��CENTER
	string change = "Y";//��Χ�����������Ƿ�ı䣬û�С�N�����ı䡰Y�������ڼ������
};
struct Flux
{
	double f1;
	double f2;
	double f3;
	double f4;
};
struct Coor//��ʾ����
{
	double x;
	double y;
};
struct Line//��ʾֱ��Ax+Bx+C=0
{
	double A;
	double B;
	double C;
};
struct polygon_mesh
{
	//��ֻ��¼�����ţ�����������������
	vector <int> node;
	vector <int>face_start;
	vector <int>face_end;
	int section;
};
namespace Oblique//б����
{
	extern double ��;
	extern double ��;
	extern int startpoint;
	extern double  ��1, p1, Ma1, u1, v1;
	extern double  ��2, p2, Ma2, u2, v2;
}
namespace Prandtl_Meyer//��������Ү������
{
	extern double  ��0, p0, ��0, ��0, ��0, ��0, Ma0, u0, v0;
	extern double  ��1, p1, ��1, ��1, ��1, ��1, Ma1, u1, v1;
	extern double  ��2, p2, ��2, ��2, ��2, ��2, Ma2, u2, v2;
}
namespace ShockwaveCross//�����ཻ
{

	extern double  ��1, p1, ��1, ��1, Ma1, u1, v1;
	extern double  ��2, p2, ��2, ��2, Ma2, u2, v2;
	extern double  ��3, p3, ��3, ��3, Ma3, u3, v3;
	extern double  ��4, p4, ��4, ��4, Ma4, u4, v4;
	extern double  ��5, p5, ��5, ��5, Ma5, u5, v5;
}

