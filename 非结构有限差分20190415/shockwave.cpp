#include"const.h"
#include<cmath>
#include"shockwave.h"
//���㼤����ز���Ǯ��𢡶��������ѧ��p228-241
using namespace ConstPara;
double get_c(double ��, double p)//���ٹ�ʽ
{
	return sqrt(�� * p / ��);
}
double get_Ma(double u, double v, double ��, double p)//�������
{
	double a = get_c(��, p);
	double velocity = sqrt(u * u + v * v);
	return velocity / a;
}
//��������ϵʽ
//7-99
double get_Mas(double p1, double p2)
{
	double Mas = 1 + (�� + 1) * (p2 / p1 - 1) / (2 * ��);
	return sqrt(Mas);
}
//7-105
double get_p2p1(double Ma1)
{
	return  2 * ��* Ma1* Ma1 / (�� + 1) - (�� - 1) / (�� + 1);
}
//7-106
double get_��2��1(double Ma1)
{
	return  (�� + 1) * Ma1 * Ma1 / ((�� - 1) * Ma1 * Ma1 + 2);
}
double get_Ma2(double Ma1)
{
	return sqrt((Ma1 * Ma1 + 2 / (�� - 1)) / (2 * �� * Ma1 * Ma1 / (�� - 1) - 1));
}
//�����������  �޶���,��ʿ����p57 3.5
double get_Mu(mesh & U, mesh & D, double ��)
{
	using namespace ConstPara;
	double �� = 1e-10;
	double Mu = 2, Mu1 = 20;
	double udn, uun, udt, uut;
	uun = get_un(U, ��);
	udn = get_un(D, ��);
	uut = get_ut(U, ��);
	udt = get_ut(D, ��);
	double cu = sqrt(�� * U.p / U.��);
	double J, Fx, fx;
	double a, b;
	while (abs(Mu - Mu1) >= ��)
	{
		Mu = Mu1;
		double J = 2 * sqrt(�� * D.p / D.��) / (�� - 1) + udn;
		a = 2 * �� * Mu1 * Mu1 - (�� - 1);
		b = (�� - 1) * Mu1 * Mu1 + 2;
		Fx = sqrt(a * b) / (Mu1 * (�� - 1)) + (Mu1 * Mu1 - 1) / Mu1 - (�� + 1) * (J - uun) / (2 * cu);
		fx = (4 * �� * Mu1 * b + 2 * (�� - 1) * Mu1 * a) / (2 * Mu1 * (�� - 1) * sqrt(a * b)) - sqrt(a * b) / (Mu1 * Mu1 * (�� - 1)) + 1 + 1 / (Mu1 * Mu1);
		Mu1 = Mu1 - Fx / fx;
	}
	return Mu;
}
//3.4
double get_udn(mesh & U, mesh & D, double Ma1, double Vs, double ��)
{
	double ��1 = get��fromMa(Ma1);
	double ��2 = 1 / ��1;
	double Ma2 = getMafrom��(��2);
	double c2 = sqrt(�� * D.p / D.��);
	double V2 = Ma2 * c2;
	return Vs-V2;
	//double a = ((�� - 1) * Ma1 + 2) / (2 * �� * Ma1 - (�� - 1));
	//double b = (a * �� * D.p) / D.��;
	//double cu = sqrt(�� * U.p / U.��);
	//double uun = U.u * sin(��) - U.v * cos(��);
	//double udn = sqrt(b) - uun +cu * Ma1 ;
	//return udn;
}

double get��fromMa(double Ma)
{
	double a = 2 + (�� - 1) * Ma * Ma;
	double b = 0;
	double c = -(�� + 1) * Ma * Ma;
	return (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
}
double getMafrom��(double ��)
{
	double ��2 = �� * ��;
	return sqrt((2 * ��2 / (�� + 1)) / (1 - (�� - 1) / (�� + 1) * ��2));
}

//����Ϊб����,����ֻ֤�ʺϲ�ǰ�ٶ�ˮƽ��v=0,������ͨ����
double get_Ma2(double Ma1, double ��)//б���������������Ǯ��𢡶��������ѧ��p241,7-131
{
	double Ma2 = (Ma1 * Ma1 + 2 / (�� - 1)) / (2 * �� / (�� - 1) * Ma1 * Ma1 * sin(��) * sin(��) - 1) + 2 / (�� - 1) * Ma1 * Ma1 * cos(��) * cos(��) / (Ma1 * Ma1 * sin(��) * sin(��) + 2 / (�� - 1));
	Ma2 = sqrt(Ma2);
	return Ma2;
}
//double get_��(double Ma1, double ��)//���������
//{
//
//	double tan�� = (Ma1 * Ma1 * sin(��) * sin(��) - 1) / ((Ma1 * Ma1 * ((�� + 1) / 2 - sin(��) * sin(��)) + 1) * tan(��));
//	double �� = atan(tan��);
//	//if (�� < 0)
//	//	�� = �� + pi;
//	return ��;
//}
double get_��from��(double Ma1, double ��)
{
	//�ҽط� �������ѧ�뼼��p132
	double �� = 40 * pi / 180, ��1 = 50 * pi / 180, ��0 = 60 * pi / 180;
	double f��1, f��0;
	if (�� == 0.92846675631531839)
		�� = 0.92846675631531839;
	while (abs(�� - ��1) > 1e-10)
	{
		��0 = ��1;
		��1 = ��;
		while (��1 > pi / 2 || ��1 < -pi / 2)
		{
			if (��1 > pi / 2)
				��1 = ��1 - pi / 2;
			else if (��1 < -pi / 2)
				��1 = ��1 + pi / 2;
		}
		while (��0 > pi / 2 || ��0 < -pi / 2)
		{
			if (��0 > pi / 2)
				��0 = ��0 - pi / 2;
			else if (��0 < -pi / 2)
				��0 = ��0 + pi / 2;
		}

		f��1 = tan(��) - (Ma1 * Ma1 * sin(��1) * sin(��1) - 1) / ((Ma1 * Ma1 * ((�� + 1) / 2 - sin(��1) * sin(��1)) + 1) * tan(��1));
		f��0 = tan(��) - (Ma1 * Ma1 * sin(��0) * sin(��0) - 1) / ((Ma1 * Ma1 * ((�� + 1) / 2 - sin(��0) * sin(��0)) + 1) * tan(��0));
		�� = ��1 - f��1 / (f��1 - f��0) * (��1 - ��0);
	}
	return ��;
}
double get_ufromMa2(double Ma2, double ��2, double p2, double ��)
{
	double a2 = sqrt(�� * p2 / ��2);
	double velocity = Ma2 * a2;
	double u = velocity * cos(��);
	return u;
}
double get_vfromMa2(double Ma2, double ��2, double p2, double ��)
{
	double a2 = sqrt(�� * p2 / ��2);
	double velocity = Ma2 * a2;
	double v = velocity * sin(��);
	return v;
}
double get_p2(double Ma1, double ��, double p1)
{
	double p2 = (2 * �� / (�� + 1) * Ma1 * Ma1 * sin(��) * sin(��) - (�� - 1) / (�� + 1)) * p1;
	return p2;
}
double get_��2(double Ma1, double ��, double ��1)
{
	double ��2 = ((�� + 1) * Ma1 * Ma1 * sin(��) * sin(��) / ((�� - 1) * Ma1 * Ma1 * sin(��) * sin(��) + 2)) * ��1;
	return ��2;
}


//����Ϊͨ����б�����㷨
double get_un(mesh & A, double ��)
{
	if (�� >= 0)
		return A.u * sin(��) - A.v * cos(��);
	else
		return A.u * sin(-��) + A.v * cos(-��);
}
double get_ut(mesh & A, double ��)
{
	if (�� >= 0)
		return A.u * cos(��) + A.v * sin(��);
	else
		return A.u * cos(-��) - A.v * sin(-��);
}

void get_down(mesh & U, mesh & D, double ��)//���β���
{
	double udn, uun, udt, uut;
	uun = get_un(U, ��);
	uut = get_ut(U, ��);
	double au = get_c(U.��, U.p);
	double Mu = uun / au;
	D.p = U.p * get_p2p1(Mu);
	D.�� = U.�� * get_��2��1(Mu);
	double Md = get_Ma2(Mu);
	double ad = get_c(D.��, D.p);
	udn = Md * ad;
	udt = uut;
	if (�� >= 0)
	{
		D.u = udt * cos(��) + udn * sin(��);
		D.v = udt * sin(��) - udn * cos(��);
	}
	else
	{
		D.u = udt * cos(-��) + udn * sin(-��);
		D.v = -udt * sin(-��) + udn * cos(-��);
	}
}
double get_��(double u, double v)
{
	double �� = atan(abs(v) / abs(u));
	if (u > 0)
	{
		if (v > 0)
			return ��;
		else if (v < 0)
			return -��;
		else
			return 0;
	}
	else if (u < 0)
	{
		if (v > 0)
			return -��;
		else if (v < 0)
			return ��;
		else
			return 0;
	}
	else
	{
		if (v > 0)
			return pi / 2;
		else if (v < 0)
			return -pi / 2;
		else
			return 0;
	}
}
double get_��(mesh & U, double p2, int type)//type=-1��ʾ����Ħ½�Ϊ��������Ϊ��
{
	double �� = 40 * pi / 180, ��1 = 50 * pi / 180, ��0 = 60 * pi / 180;
	double uun0, uun1;
	double Mu0, Mu1;
	double p20, p21;
	double f��1, f��0;
	while (abs(�� - ��1) > 1e-15)
	{
		while (�� > pi / 2 || �� < -pi / 2)
		{
			if (�� > pi / 2)
				�� = �� - pi / 2;
			else if (�� < -pi / 2)
				�� = �� + pi / 2;
		}
		��0 = ��1;
		��1 = ��;
		if (type == -1)
		{
			if (��0 > 0)
				��0 = -��0;
			if (��1 > 0)
				��1 = -��1;
		}
		else if (type == 1)
		{
			if (��0 < 0)
				��0 = -��0;
			if (��1 < 0)
				��1 = -��1;
		}
		uun0 = get_un(U, ��0);
		uun1 = get_un(U, ��1);
		Mu0 = get_Ma(uun0, 0, U.��, U.p);
		Mu1 = get_Ma(uun1, 0, U.��, U.p);
		p20 = U.p * get_p2p1(Mu0);
		p21 = U.p * get_p2p1(Mu1);
		f��0 = p20 - p2;
 		f��1 = p21 - p2;
		if (f��0 == f��1)
			break;
		else
			�� = ��1 - f��1 / (f��1 - f��0) * (��1 - ��0);
	}
	return ��;

}