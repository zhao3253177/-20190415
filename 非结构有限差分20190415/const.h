#pragma once
#include<cmath>
#include<string>
#include<vector>
using std::vector;
using std::string;
namespace ConstPara//常数
{
	const double CFL = 0.5;
	const double t_end = 100;//结束时刻
	const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628;
	const double γ = 1.4;//比热比
	const double random = 0.495;
}
namespace Init//初始网格
{
	extern double ρ0;
	extern double v0;
	extern double p0;
	extern double u0;
}
namespace Normal//初始网格
{
	extern double ρ1;
	extern double v1;
	extern double p1;
	extern double u1;
	extern double ρ2;
	extern double v2;
	extern double p2;
	extern double u2;
}

namespace MeshPara
{
	const string FlowType = "normal";
	//目前有，膨胀波Prandtl-Meyer;正激波normal，斜激波oblique ，激波相交intersection,其他other
	const int meshType = 20;
	//10正六边形网格，11扰动正六边形网格，12砖块网格，
	//20长方形结构网格
	const string methodType = "C";
	//计算方法，激波捕捉C-Capturing，激波装配F-Fitting
	const int method[12][4] = {
{1,2,0,0},//方法3
{0,1,2,2},
{2,0,1,1},

{0,1,1,2},//方法2
{1,2,2,0},
{2,0,0,1},


{0,0,1,2},//方法1
{1,1,2,0},
{2,2,0,1},

{0,1,2,0},//方法4
{2,0,1,2},
{1,2,0,1},

//{2,1,1,0},//方法22
//{0,2,2,1},
//{1,0,0,2},
//
//{2,1,0,0},//方法32
//{1,0,2,2},
//{0,2,1,1},
//
//{0,0,2,1},//方法12
//{1,1,0,2},
//{2,2,1,0},
//
//{0,2,1,0},//方法42
//{2,1,0,2},
//{1,0,2,1},

	};

	const int Xnum = 1000;//x方向点的个数20,40
	const int Ynum = 25;//y方向点的个数43,87
	const int Pnum = Xnum * Ynum;//一共点的个数
	const double dx = 0.0167/2.5;
	const double dy = sqrt(3.0) / 2 * dx;
}

struct mesh
{
	double x;
	double y;
	double ρ;//密度
	double u;//水平速度
	double v;//竖直速度
	double p;//压强
	double x0;//动网格的初始网格
	double y0;
	double um;//格点运动速度
	double vm;
	int section;//分区
	int step = 0;//表明该点在何步骤更新的
	int neiborsec = -1;//相邻分区，如果该值为负，则表示该点不是激波点，反之为相邻分区共用的激波点
	int neiborsec_ad = -1;//相邻分区的地址
	vector <double> ξx;
	vector <double> ξy;
	vector <double> ξt;
	vector <double> ηx;
	vector <double> ηy;
	vector <double> ηt;
	vector <double> xξ;
	vector <double> xη;
	vector <double> xτ;
	vector <double> yξ;
	vector <double> yη;
	vector <double> yτ;
	vector <double> J;//雅可比行列式
	vector <int>neibor;//记录该格点的相邻格点位置信息
	vector <int>moveConnct;//运动关联点，即本点运动与该moveConnect点有关
	int neibor1[4] = { 0 };//当相邻格点超过四个，筛选四个入此处，坐标变换和计算所用到的相邻格点都用此处
	string type;//格点类型，分为上下左右边界以及内部，U，D,L,R,IN,激波点SHOCK，接触间断点DISCON（contact discontinuity）,激波相交点CENTER
	string change = "Y";//周围网格物理量是否改变，没有“N”，改变“Y”；用于计算加速
};
struct Flux
{
	double f1;
	double f2;
	double f3;
	double f4;
};
struct Coor//表示坐标
{
	double x;
	double y;
};
struct Line//表示直线Ax+Bx+C=0
{
	double A;
	double B;
	double C;
};
struct polygon_mesh
{
	//都只记录点的序号，用于输出多边形网格
	vector <int> node;
	vector <int>face_start;
	vector <int>face_end;
	int section;
};
namespace Oblique//斜激波
{
	extern double β;
	extern double δ;
	extern int startpoint;
	extern double  ρ1, p1, Ma1, u1, v1;
	extern double  ρ2, p2, Ma2, u2, v2;
}
namespace Prandtl_Meyer//普朗特麦耶尔流动
{
	extern double  ρ0, p0, θ0, λ0, δ0, μ0, Ma0, u0, v0;
	extern double  ρ1, p1, θ1, λ1, δ1, μ1, Ma1, u1, v1;
	extern double  ρ2, p2, θ2, λ2, δ2, μ2, Ma2, u2, v2;
}
namespace ShockwaveCross//激波相交
{

	extern double  ρ1, p1, β1, δ1, Ma1, u1, v1;
	extern double  ρ2, p2, β2, δ2, Ma2, u2, v2;
	extern double  ρ3, p3, β3, δ3, Ma3, u3, v3;
	extern double  ρ4, p4, β4, δ4, Ma4, u4, v4;
	extern double  ρ5, p5, β5, δ5, Ma5, u5, v5;
}

