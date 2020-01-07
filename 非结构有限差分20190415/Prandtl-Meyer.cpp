#include<cmath>
#include"const.h"
#include"shockwave.h"
using ConstPara::��;
double get��from��(double ��)
{
	double k = (�� - 1) / (�� + 1);
	double ��2 = �� * ��;
	return sqrt(1 / k)*atan(sqrt((k*(��2 - 1)) / (1 - k * ��2))) - atan(sqrt((��2 - 1) / (1 - k * ��2)));
}
double get��fromMa(double Ma)
{
	return asin(1 / Ma);
}
double get��fromMa(double Ma)
{
	double �� = get��fromMa(Ma);
	return sqrt((�� + 1) / (�� - 1))*asin(sqrt((�� - 1)*(��*�� - 1) / 2));
}
double get��from��(double ��)
{
	double �� = 1e-6;
	double ��1 = 1, ��2 = sqrt((�� + 1) / (�� - 1));
	double ��M = (��1 + ��2) / 2;
	double ��1, ��M=0, ��2;
	while (abs(��M - ��) > 1e-6)
	{
		��M = (��1 + ��2) / 2;
		��1 = get��from��(��1);
		��M = get��from��(��M);
		��2 = get��from��(��2);
		if ((��1 - ��)*(��M - ��) < 0)
			��1 = ��1, ��2 = ��M;
		else
			��1 = ��M, ��2 = ��2;
	}
	return (��1 + ��2) / 2;
}
//����Ϊ��������ʽ�������������Ͳ�����
double getTfrom��(double ��,double T0)
{
	double temp = 1 - (�� - 1)*��*�� / (�� + 1);
	return temp*T0;
}
double getT0from��andT(double ��, double T)
{
	double temp = 1 - (�� - 1)*��*�� / (�� + 1);
	return T / temp;
}
double getpfrom��(double ��, double p0)
{
	double temp = 1 - (�� - 1)*��*�� / (�� + 1);
	temp = pow(temp, �� / (�� - 1));
	return temp * p0;
}
double getp0from��andp(double ��, double p)
{
	double temp = 1 - (�� - 1)*��*�� / (�� + 1);
	temp = pow(temp, �� / (�� - 1));
	return p / temp;
}
double get��from��(double ��, double ��0)
{
	double temp = 1 - (�� - 1)*��*�� / (�� + 1);
	temp = pow(temp, 1 / (�� - 1));
	return temp * ��0;
}
double get��0from��and��(double ��, double ��)
{
	double temp = 1 - (�� - 1)*��*�� / (�� + 1);
	temp = pow(temp, 1 / (�� - 1));
	return �� / temp;
}
