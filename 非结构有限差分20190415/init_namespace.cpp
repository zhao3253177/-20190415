#include"const.h"
#include"shockwave.h"
#include"Prandtl-Meyer.h"
#include"patition.h"
using namespace MeshPara;
using ConstPara::t_end;


namespace Init//��ʼ����
{
	double ��0 = 1;
	double v0 = 0;
	double p0 = 1;
	double u0 = 1 * sqrt(ConstPara::�� * p0 / ��0);
}
namespace Normal
{
	//double ��1 = 1;
	//double v1 = 0;
	//double p1 = 1;
	//double u1 = 1.1 * sqrt(ConstPara::�� * p1 / ��1);
	//double ��2;
	//double v2;
	//double p2;
	//double u2;
	double ��1 = 5;
	double v1 = 0;
	double p1 = 4;
	double u1 = 3 * sqrt(ConstPara::�� * p1 / ��1);
	double ��2 = 3.0919952786742946;
	double v2 = 0;
	double p2 = 5.5743027276648576;
	double u2 = 2.2088635344975316;

}
namespace Oblique//б����
{
	double �� = 40* ConstPara::pi / 180;
	double ��1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 10* sqrt(ConstPara::�� * p1 / ��1);
	int startpoint = MeshPara::Xnum * 3 / 10;
	double Ma1 = get_Ma(u1, v1, ��1, p1);
	double Ma2 = get_Ma2(Ma1, ��);
	double p2 = get_p2(Ma1, ��, p1);
	double ��2 = get_��2(Ma1, ��, ��1);
	double �� = get_��(Ma1, ��);
	double u2 = get_ufromMa2(Ma2, ��2, p2, ��);
	double v2 = get_vfromMa2(Ma2, ��2, p2, ��);
}
namespace Prandtl_Meyer//��������Ү������
{
	double ��1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 1.3 * sqrt(ConstPara::�� * p1 / ��1);
	double Ma1 = get_Ma(u1, v1, ��1, p1);
	double ��1 = get��fromMa(Ma1);
	double ��1 = get��from��(��1);
	double ��1 = get��fromMa(Ma1);
	double ��1 = get��fromMa(Ma1);
	double ��2 = 10.0 / 180 * ConstPara::pi;
	double ��2 = get��from��(��2);
	double Ma2 = getMafrom��(��2);
	double ��2 = get��fromMa(Ma2);
	double ��2 = get��fromMa(Ma2);
	double p0 = getp0from��andp(��1, p1);
	double ��0 = get��0from��and��(��1, ��1);
	double p2 = getpfrom��(��2, p0);
	double ��2 = get��from��(��2, ��0);
	double u2 = Ma2 * sqrt(p2 * ConstPara::�� / ��2) * cos(��2);
	double v2 = -Ma2 * sqrt(p2 * ConstPara::�� / ��2) * sin(��2);

}
namespace ShockwaveCross//�����ཻ
{
	double ��1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 5* sqrt(ConstPara::�� * p1 / ��1);
	double ��2 = -41 / 180.0 * ConstPara::pi;
	double ��3 = 30/ 180.0 * ConstPara::pi;
	double Ma1 = get_Ma(u1, v1, ��1, p1);

	double Ma2 = get_Ma2(Ma1, ��2);
	double p2 = get_p2(Ma1, ��2, p1);
	double ��2 = get_��2(Ma1, ��2, ��1);
	double ��2 = get_��(Ma1, ��2);
	double u2 = get_ufromMa2(Ma2, ��2, p2, ��2);
	double v2 = get_vfromMa2(Ma2, ��2, p2, ��2);

	double Ma3 = get_Ma2(Ma1, ��3);
	double p3 = get_p2(Ma1, ��3, p1);
	double ��3 = get_��2(Ma1, ��3, ��1);
	double ��3 = get_��(Ma1, ��3);
	double u3 = get_ufromMa2(Ma3, ��3, p3, ��3);
	double v3 = get_vfromMa2(Ma3, ��3, p3, ��3);

	double p4 = 1;
	double ��4 = 1;
	double ��4 = 1;
	double u4 = 1;
	double v4 = 1;
	double ��4 = 1;

	double p5 = 1;
	double ��5 = 1;
	double ��5 = 1;
	double u5 = 1;
	double v5 = 1;
	double ��5 = 1;
}
