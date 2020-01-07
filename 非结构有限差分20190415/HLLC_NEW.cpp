#include<cmath>
#include"const.h"
#include<iostream>
#include<omp.h>
#include"functions.h"
using namespace std;
using namespace ConstPara;
//2άŷ����������ɢ���ŵ���������������ѧ�̡̳�p433 13.3.7
Flux HLLC_��(mesh& CL, mesh& CR, mesh& C, int method)//�������Ҳ��㣬����任�ο���
{
	double ��x = C.��x[method];
	double ��y = C.��y[method];
	double ��t = C.��t[method];
	double J = C.J[method];
	double D�� = sqrt(��x * ��x + ��y * ��y);
	double ��1 = ��x / D��;
	double ��2 = ��y / D��;
	double ��3 = ��t / D��;
	double ��cL = CL.u * ��1 + CL.v * ��2 + ��3;
	double ��cR = CR.u * ��1 + CR.v * ��2 + ��3;
	double aL = sqrt(�� * CL.p / CL.��);
	double aR = sqrt(�� * CR.p / CR.��);
	double a = sqrt(CL.��);
	double b = sqrt(CR.��);
	double SL = ��cL - aL;
	double SR = ��cR + aR;

	double EL = CL.p / (�� - 1) + 0.5 * CL.�� * CL.u * CL.u + 0.5 * CL.�� * CL.v * CL.v;
	double ER = CR.p / (�� - 1) + 0.5 * CR.�� * CR.u * CR.u + 0.5 * CR.�� * CR.v * CR.v;
	Flux FL, DU;
	FL.f1 = (CL.�� * ��cL);
	FL.f2 = (CL.�� * CL.u * ��cL + ��1 * CL.p);
	FL.f3 = (CL.�� * CL.v * ��cL + ��2 * CL.p);
	FL.f4 = (��cL * (EL + CL.p) - ��3 * CL.p);
	DU.f1 = (CR.�� * ��cR);
	DU.f2 = (CR.�� * CR.u * ��cR + ��1 * CR.p);
	DU.f3 = (CR.�� * CR.v * ��cR + ��2 * CR.p);
	DU.f4 = (��cR * (ER + CR.p) - ��3 * CR.p);
	double SM = (CL.p - CR.p + CL.�� * ��cL * (��cL - SL) + CR.�� * ��cR * (SR - ��cR)) / (CR.�� * (SR - ��cR) + CL.�� * (��cL - SL));
	double pM = 0.5 * (CL.�� * (��cL - SL) * (��cL - SM) + CR.�� * (��cR - SR) * (��cR - SM) + CR.p + CL.p);
	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
	double ��L = 1 / (SL - SM);
	double ��R = 1 / (SR - SM);

	double ��LM = ��L * CL.�� * (SL - ��cL);
	double ��uLS = ��L * ((SL - ��cL) * (CL.�� * CL.u) + (pM - CL.p) * ��1);
	double ��vLS = ��L * ((SL - ��cL) * (CL.�� * CL.v) + (pM - CL.p) * ��2);
	double eLS = ��L * ((SL - ��cL) * EL - CL.p * ��cL + pM * SM - (pM - CL.p) * ��3);

	double ��RM = ��R * CR.�� * (SR - ��cR);
	double ��uRS = ��R * ((SR - ��cR) * (CR.�� * CR.u) + (pM - CR.p) * ��1);
	double ��vRS = ��R * ((SR - ��cR) * (CR.�� * CR.v) + (pM - CR.p) * ��2);
	double eRS = ��R * ((SR - ��cR) * ER - CR.p * ��cR + pM * SM - (pM - CR.p) * ��3);

	Flux FML, FMR;
	FML.f1 = SM * ��LM;
	FML.f2 = ��uLS * SM + pM * ��1;
	FML.f3 = ��vLS * SM + pM * ��2;
	FML.f4 = (eLS + pM) * SM - pM * ��3;
	FMR.f1 = SM * ��RM;
	FMR.f2 = ��uRS * SM + pM * ��1;
	FMR.f3 = ��vRS * SM + pM * ��2;
	FMR.f4 = (eRS + pM) * SM - pM * ��3;

	Flux F_HLLC;
	if (SL > 0)
		F_HLLC = FL;
	else if (SL <= 0 && SM > 0)
		F_HLLC = FML;
	else if (SM <= 0 && SR >= 0)
		F_HLLC = FMR;
	else
		F_HLLC = DU;
	F_HLLC.f1 = F_HLLC.f1 * D��;
	F_HLLC.f2 = F_HLLC.f2 * D��;
	F_HLLC.f3 = F_HLLC.f3 * D��;
	F_HLLC.f4 = F_HLLC.f4 * D��;

	return F_HLLC;
}
Flux HLLC_��2(mesh & CL, mesh & CR, mesh & C, int method)//�������Ҳ��㣬����任�ο���
{
	double ��x = C.��x[method];
	double ��y = C.��y[method];
	double ��t = C.��t[method];
	double J = C.J[method];
	double D�� = sqrt(��x * ��x + ��y * ��y);
	double ��1 = ��x / D��;
	double ��2 = ��y / D��;
	double ��3 = ��t / D��;
	double aL = sqrt(�� * CL.p / CL.��);
	double aR = sqrt(�� * CR.p / CR.��);
	double a = sqrt(CL.��);
	double b = sqrt(CR.��);
	double ubl = CL.u * ��1 + CL.v * ��2 + ��3;
	double ubr = CR.u * ��1 + CR.v * ��2 + ��3;
	double ubp = (a * ubl + b * ubr) / (a + b);
	double aM = (a * aL * aL + b * aR * aR) / (a + b) + 0.5 * (a * b / ((a + b) * (a + b))) * ((CL.u - CR.u) * (CL.u - CR.u) + (CL.v - CR.v) * (CL.v - CR.v));


	double conl = min(ubl - aL, ubp - aM);
	double conr = max(ubr + aR, ubp + aM);
	double s1 = CR.�� * (conr - ubr);
	double s2 = CL.�� * (ubl - conl);
	double ss = (s1 * ubr + s2 * ubl + CL.p - CR.p) / (s1 + s2);
	double ps = s2 * (ubl - ss) + CL.p;
	Flux FL, FR;
	FL.f1 = CL.��;
	FL.f2 = CL.�� * CL.u;
	FL.f3 = CL.�� * CL.v;
	FL.f4 = CL.p * �� / (�� - 1) + 0.5 * CL.�� * CL.u * CL.u + 0.5 * CL.�� * CL.v * CL.v; ;
	FR.f1 = CR.��;
	FR.f2 = CR.�� * CR.u;
	FR.f3 = CR.�� * CR.v;
	FR.f4 = CR.p * �� / (�� - 1) + 0.5 * CR.�� * CR.u * CR.u + 0.5 * CR.�� * CR.v * CR.v;
	double par1 = (ubl - conl) / (ss - conl);
	double par2 = (ps - CL.p) / (ss - conl);

	Flux FML, FMR;
	FML.f1 = par1 * FL.f1;
	FML.f2 = par1 * FL.f2 - par2 * ��1;
	FML.f3 = par1 * FL.f3 - par2 * ��2;
	FML.f4 = par1 * FL.f4 + par2 * ��3;
	double par3 = (ubr - conr) / (ss - conr);
	double par4 = (ps - CR.p) / (ss - conr);

	FMR.f1 = par3 * FR.f1;
	FMR.f2 = par3 * FR.f2 - par4 * ��1;
	FMR.f3 = par3 * FR.f3 - par4 * ��2;
	FMR.f4 = par3 * FR.f4 + par4 * ��3;

	Flux F_HLLC;
	if (conl > 0)
	{
		F_HLLC.f1 = FL.f1 * ubl;
		F_HLLC.f2 = FL.f2 * ubl + ��1 * CL.p;
		F_HLLC.f3 = FL.f3 * ubl + ��2 * CL.p;
		F_HLLC.f4 = FL.f4 * ubl - ��3 * CL.p;
	}
	else if (conl <= 0 && ss > 0)
	{
		F_HLLC.f1 = FML.f1 * ss;
		F_HLLC.f2 = FML.f2 * ss + ��1 * ps;
		F_HLLC.f3 = FML.f3 * ss + ��2 * ps;
		F_HLLC.f4 = FML.f4 * ss - ��3 * ps;
	}
	else if (ss <= 0 && conr >= 0)
	{
		F_HLLC.f1 = FMR.f1 * ss;
		F_HLLC.f2 = FMR.f2 * ss + ��1 * ps;
		F_HLLC.f3 = FMR.f3 * ss + ��2 * ps;
		F_HLLC.f4 = FMR.f4 * ss - ��3 * ps;
	}
	else
	{
		F_HLLC.f1 = FR.f1 * ubr;
		F_HLLC.f2 = FR.f2 * ubr + ��1 * CR.p;
		F_HLLC.f3 = FR.f3 * ubr + ��2 * CR.p;
		F_HLLC.f4 = FR.f4 * ubr - ��3 * CR.p;
	}
	F_HLLC.f1 = F_HLLC.f1 * D��;
	F_HLLC.f2 = F_HLLC.f2 * D��;
	F_HLLC.f3 = F_HLLC.f3 * D��;
	F_HLLC.f4 = F_HLLC.f4 * D��;

	return F_HLLC;
}

Flux HLLC_��(mesh& CD, mesh& CU, mesh& C, int method)
{
	double ��x = C.��x[method];
	double ��y = C.��y[method];
	double ��t = C.��t[method];
	double J = C.J[method];

	double D�� = sqrt(��x * ��x + ��y * ��y);
	double ��1 = ��x / D��;
	double ��2 = ��y / D��;
	double ��3 = ��t / D��;
	double ��cU = CU.u * ��1 + CU.v * ��2;
	double ��cD = CD.u * ��1 + CD.v * ��2;
	double aU = sqrt(�� * CU.p / CU.��);
	double aD = sqrt(�� * CD.p / CD.��);
	double EU = CU.p / (�� - 1) + 0.5 * CU.�� * CU.u * CU.u + 0.5 * CU.�� * CU.v * CU.v;
	double ED = CD.p / (�� - 1) + 0.5 * CD.�� * CD.u * CD.u + 0.5 * CD.�� * CD.v * CD.v;
	Flux GD, GU;
	GD.f1 = (CD.�� * ��cD);
	GD.f2 = (CD.�� * CD.u * ��cD + ��1 * CD.p);
	GD.f3 = (CD.�� * CD.v * ��cD + ��2 * CD.p);
	GD.f4 = (��cD * (ED + CD.p) - ��3 * CD.p);

	GU.f1 = (CU.�� * ��cU);
	GU.f2 = (CU.�� * CU.u * ��cU + ��1 * CD.p);
	GU.f3 = (CU.�� * CU.v * ��cU + ��2 * CD.p);
	GU.f4 = (��cU * (EU + CU.p) - ��3 * CD.p);

	double SD = ��cD - aD;
	double SU = ��cU + aU;
	double SM = (CD.p - CU.p - CD.�� * ��cD * (SD - ��cD) + CU.�� * ��cU * (SU - ��cU)) / (CU.�� * (SU - ��cU) - CD.�� * (SD - ��cD));
	double pM = 0.5 * (CD.�� * (��cD - SD) * (��cD - SM) + CU.�� * (��cU - SU) * (��cU - SM) + CU.p + CD.p);
	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
	double ��D = 1 / (SD - SM);
	double ��U = 1 / (SU - SM);

	double ��DM = ��D * CD.�� * (SD - ��cD);
	double ��uDS = ��D * ((SD - ��cD) * (CD.�� * CD.u) + (pM - CD.p) * ��1);
	double ��vDS = ��D * ((SD - ��cD) * (CD.�� * CD.v) + (pM - CD.p) * ��2);
	double eDS = ��D * ((SD - ��cD) * ED - CD.p * ��cD + pM * SM - (pM - CD.p) * ��3);

	double ��UM = ��U * CU.�� * (SU - ��cU);
	double ��uUS = ��U * ((SU - ��cU) * (CU.�� * CU.u) + (pM - CU.p) * ��1);
	double ��vUS = ��U * ((SU - ��cU) * (CU.�� * CU.v) + (pM - CU.p) * ��2);
	double eUS = ��U * ((SU - ��cU) * EU - CU.p * ��cU + pM * SM - (pM - CU.p) * ��3);

	Flux  FMD, FMU;

	FMD.f1 = SM * ��DM;
	FMD.f2 = ��uDS * SM + pM * ��1;
	FMD.f3 = ��vDS * SM + pM * ��2;
	FMD.f4 = (eDS + pM) * SM - pM * ��3;
	FMU.f1 = SM * ��UM;
	FMU.f2 = ��uUS * SM + pM * ��1;
	FMU.f3 = ��vUS * SM + pM * ��2;
	FMU.f4 = (eUS + pM) * SM - pM * ��3;

	Flux G_HLLC;
	if (SD >= 0)
		G_HLLC = GD;
	else if (SD <= 0 && SM >= 0)
		G_HLLC = FMD;
	else if (SM <= 0 && SU >= 0)
		G_HLLC = FMU;
	else
		G_HLLC = GU;
	G_HLLC.f1 = G_HLLC.f1 * D��;
	G_HLLC.f2 = G_HLLC.f2 * D��;
	G_HLLC.f3 = G_HLLC.f3 * D��;
	G_HLLC.f4 = G_HLLC.f4 * D��;

	return G_HLLC;
}
Flux HLLC_��2(mesh& CL, mesh& CR, mesh& C, int method)
{
	double ��x = C.��x[method];
	double ��y = C.��y[method];
	double ��t = C.��t[method];
	double J = C.J[method];
	double D�� = sqrt(��x * ��x + ��y * ��y);
	double ��1 = ��x / D��;
	double ��2 = ��y / D��;
	double ��3 = ��t / D��;
	double aL = sqrt(�� * CL.p / CL.��);
	double aR = sqrt(�� * CR.p / CR.��);
	double a = sqrt(CL.��);
	double b = sqrt(CR.��);
	double ubl = CL.u * ��1 + CL.v * ��2 + ��3;
	double ubr = CR.u * ��1 + CR.v * ��2 + ��3;
	double ubp = (a * ubl + b * ubr) / (a + b);
	double aM = (a * aL * aL + b * aR * aR) / (a + b) + 0.5 * (a * b / ((a + b) * (a + b))) * ((CL.u - CR.u) * (CL.u - CR.u) + (CL.v - CR.v) * (CL.v - CR.v));


	double conl = min(ubl - aL, ubp - aM);
	double conr = max(ubr + aR, ubp + aM);
	double s1 = CR.�� * (conr - ubr);
	double s2 = CL.�� * (ubl - conl);
	double ss = (s1 * ubr + s2 * ubl + CL.p - CR.p) / (s1 + s2);
	double ps = s2 * (ubl - ss) + CL.p;
	Flux FL, FR;
	FL.f1 = CL.��;
	FL.f2 = CL.�� * CL.u;
	FL.f3 = CL.�� * CL.v;
	FL.f4 = CL.p * �� / (�� - 1) + 0.5 * CL.�� * CL.u * CL.u + 0.5 * CL.�� * CL.v * CL.v; ;
	FR.f1 = CR.��;
	FR.f2 = CR.�� * CR.u;
	FR.f3 = CR.�� * CR.v;
	FR.f4 = CR.p * �� / (�� - 1) + 0.5 * CR.�� * CR.u * CR.u + 0.5 * CR.�� * CR.v * CR.v;
	double par1 = (ubl - conl) / (ss - conl);
	double par2 = (ps - CL.p) / (ss - conl);

	Flux FML, FMR;
	FML.f1 = par1 * FL.f1;
	FML.f2 = par1 * FL.f2 - par2 * ��1;
	FML.f3 = par1 * FL.f3 - par2 * ��2;
	FML.f4 = par1 * FL.f4 + par2 * ��3;
	double par3 = (ubr - conr) / (ss - conr);
	double par4 = (ps - CR.p) / (ss - conr);

	FMR.f1 = par3 * FR.f1;
	FMR.f2 = par3 * FR.f2 - par4 * ��1;
	FMR.f3 = par3 * FR.f3 - par4 * ��2;
	FMR.f4 = par3 * FR.f4 + par4 * ��3;

	Flux F_HLLC;
	if (conl > 0)
	{
		F_HLLC.f1 = FL.f1 * ubl;
		F_HLLC.f2 = FL.f2 * ubl + ��1 * CL.p;
		F_HLLC.f3 = FL.f3 * ubl + ��2 * CL.p;
		F_HLLC.f4 = FL.f4 * ubl - ��3 * CL.p;
	}
	else if (conl <= 0 && ss > 0)
	{
		F_HLLC.f1 = FML.f1 * ss;
		F_HLLC.f2 = FML.f2 * ss + ��1 * ps;
		F_HLLC.f3 = FML.f3 * ss + ��2 * ps;
		F_HLLC.f4 = FML.f4 * ss - ��3 * ps;
	}
	else if (ss <= 0 && conr >= 0)
	{
		F_HLLC.f1 = FMR.f1 * ss;
		F_HLLC.f2 = FMR.f2 * ss + ��1 * ps;
		F_HLLC.f3 = FMR.f3 * ss + ��2 * ps;
		F_HLLC.f4 = FMR.f4 * ss - ��3 * ps;
	}
	else
	{
		F_HLLC.f1 = FR.f1 * ubr;
		F_HLLC.f2 = FR.f2 * ubr + ��1 * CR.p;
		F_HLLC.f3 = FR.f3 * ubr + ��2 * CR.p;
		F_HLLC.f4 = FR.f4 * ubr - ��3 * CR.p;
	}
	F_HLLC.f1 = F_HLLC.f1 * D�� ;
	F_HLLC.f2 = F_HLLC.f2 * D�� ;
	F_HLLC.f3 = F_HLLC.f3 * D�� ;
	F_HLLC.f4 = F_HLLC.f4 * D�� ;

	return F_HLLC;
}
