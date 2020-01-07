#include"const.h"
#include"shockwave.h"
#include"Prandtl-Meyer.h"
#include"patition.h"
using namespace MeshPara;
using ConstPara::t_end;


namespace Init//³õÊ¼Íø¸ñ
{
	double ¦Ñ0 = 1;
	double v0 = 0;
	double p0 = 1;
	double u0 = 1 * sqrt(ConstPara::¦Ã * p0 / ¦Ñ0);
}
namespace Normal
{
	//double ¦Ñ1 = 1;
	//double v1 = 0;
	//double p1 = 1;
	//double u1 = 1.1 * sqrt(ConstPara::¦Ã * p1 / ¦Ñ1);
	//double ¦Ñ2;
	//double v2;
	//double p2;
	//double u2;
	double ¦Ñ1 = 5;
	double v1 = 0;
	double p1 = 4;
	double u1 = 3 * sqrt(ConstPara::¦Ã * p1 / ¦Ñ1);
	double ¦Ñ2 = 3.0919952786742946;
	double v2 = 0;
	double p2 = 5.5743027276648576;
	double u2 = 2.2088635344975316;

}
namespace Oblique//Ð±¼¤²¨
{
	double ¦Â = 40* ConstPara::pi / 180;
	double ¦Ñ1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 10* sqrt(ConstPara::¦Ã * p1 / ¦Ñ1);
	int startpoint = MeshPara::Xnum * 3 / 10;
	double Ma1 = get_Ma(u1, v1, ¦Ñ1, p1);
	double Ma2 = get_Ma2(Ma1, ¦Â);
	double p2 = get_p2(Ma1, ¦Â, p1);
	double ¦Ñ2 = get_¦Ñ2(Ma1, ¦Â, ¦Ñ1);
	double ¦Ä = get_¦Ä(Ma1, ¦Â);
	double u2 = get_ufromMa2(Ma2, ¦Ñ2, p2, ¦Ä);
	double v2 = get_vfromMa2(Ma2, ¦Ñ2, p2, ¦Ä);
}
namespace Prandtl_Meyer//ÆÕÀÊÌØÂóÒ®¶ûÁ÷¶¯
{
	double ¦Ñ1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 1.3 * sqrt(ConstPara::¦Ã * p1 / ¦Ñ1);
	double Ma1 = get_Ma(u1, v1, ¦Ñ1, p1);
	double ¦Ë1 = get¦ËfromMa(Ma1);
	double ¦Ä1 = get¦Äfrom¦Ë(¦Ë1);
	double ¦Ì1 = get¦ÌfromMa(Ma1);
	double ¦È1 = get¦ÈfromMa(Ma1);
	double ¦Ä2 = 10.0 / 180 * ConstPara::pi;
	double ¦Ë2 = get¦Ëfrom¦Ä(¦Ä2);
	double Ma2 = getMafrom¦Ë(¦Ë2);
	double ¦Ì2 = get¦ÌfromMa(Ma2);
	double ¦È2 = get¦ÈfromMa(Ma2);
	double p0 = getp0from¦Ëandp(¦Ë1, p1);
	double ¦Ñ0 = get¦Ñ0from¦Ëand¦Ñ(¦Ë1, ¦Ñ1);
	double p2 = getpfrom¦Ë(¦Ë2, p0);
	double ¦Ñ2 = get¦Ñfrom¦Ë(¦Ë2, ¦Ñ0);
	double u2 = Ma2 * sqrt(p2 * ConstPara::¦Ã / ¦Ñ2) * cos(¦Ä2);
	double v2 = -Ma2 * sqrt(p2 * ConstPara::¦Ã / ¦Ñ2) * sin(¦Ä2);

}
namespace ShockwaveCross//¼¤²¨Ïà½»
{
	double ¦Ñ1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 5* sqrt(ConstPara::¦Ã * p1 / ¦Ñ1);
	double ¦Â2 = -41 / 180.0 * ConstPara::pi;
	double ¦Â3 = 30/ 180.0 * ConstPara::pi;
	double Ma1 = get_Ma(u1, v1, ¦Ñ1, p1);

	double Ma2 = get_Ma2(Ma1, ¦Â2);
	double p2 = get_p2(Ma1, ¦Â2, p1);
	double ¦Ñ2 = get_¦Ñ2(Ma1, ¦Â2, ¦Ñ1);
	double ¦Ä2 = get_¦Ä(Ma1, ¦Â2);
	double u2 = get_ufromMa2(Ma2, ¦Ñ2, p2, ¦Ä2);
	double v2 = get_vfromMa2(Ma2, ¦Ñ2, p2, ¦Ä2);

	double Ma3 = get_Ma2(Ma1, ¦Â3);
	double p3 = get_p2(Ma1, ¦Â3, p1);
	double ¦Ñ3 = get_¦Ñ2(Ma1, ¦Â3, ¦Ñ1);
	double ¦Ä3 = get_¦Ä(Ma1, ¦Â3);
	double u3 = get_ufromMa2(Ma3, ¦Ñ3, p3, ¦Ä3);
	double v3 = get_vfromMa2(Ma3, ¦Ñ3, p3, ¦Ä3);

	double p4 = 1;
	double ¦Ñ4 = 1;
	double ¦Ä4 = 1;
	double u4 = 1;
	double v4 = 1;
	double ¦Â4 = 1;

	double p5 = 1;
	double ¦Ñ5 = 1;
	double ¦Ä5 = 1;
	double u5 = 1;
	double v5 = 1;
	double ¦Â5 = 1;
}
