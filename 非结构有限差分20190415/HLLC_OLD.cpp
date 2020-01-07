//#include<cmath>
//#include"const.h"
//#include<iostream>
//#include<omp.h>
//
//using namespace std;
////2άŷ����������ɢ���ŵ���������������ѧ�̡̳�p433 13.3.7
//double * HLLC_��(mesh CL, mesh CR, mesh C, int method)//�������Ҳ��㣬����任�ο���
//{
//	double ��x = C.��x[method];//�˴����Ǧ�x,�Ǧ�x*J-1=y��
//	double ��y = C.��y[method];
//	double J = C.J[method];
//	double D�� = sqrt(��x*��x + ��y * ��y);
//	double ��cL = (��x * CL.u + ��y * CL.v) / D��;
//	double ��cR = (��x * CR.u + ��y * CR.v) / D��;
//	double aL = sqrt(��*CL.p / CL.��);
//	double aR = sqrt(��*CR.p / CR.��);
//	double EL = CL.p / (�� - 1) + 0.5*CL.�� * CL.u * CL.u + 0.5*CL.�� * CL.v * CL.v;
//	double ER = CR.p / (�� - 1) + 0.5*CR.�� * CR.u * CR.u + 0.5*CR.�� * CR.v * CR.v;
//	double FL[4], FR[4];
//
//	FL[0] = (CL.��*��cL);
//	FL[1] = (CL.��*CL.u*��cL + ��x * CL.p);
//	FL[2] = (CL.��*CL.v*��cL + ��y * CL.p);
//	FL[3] = (��cL * (EL + CL.p));
//	FR[0] = (CR.��*��cR);
//	FR[1] = (CR.��*CR.u*��cR + ��x * CR.p);
//	FR[2] = (CR.��*CR.v*��cR + ��y * CR.p);
//	FR[3] = (��cR * (ER + CR.p));
//	double SL = ��cL - aL;
//	double SR = ��cR + aR;
//	double SM = (CL.p - CR.p + CL.�� * ��cL * (��cL - SL) + CR.�� * ��cR* (SR - ��cR)) / (CR.�� * (SR - ��cR) - CL.�� * (SL - ��cL));
//	double pM = 0.5*(CL.�� * (��cL - SL)*(��cL - SM) + CR.�� * (��cR - SR)*(��cR - SM) + CR.p + CL.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ��L = 1 / (SL - SM);
//	double ��R = 1 / (SR - SM);
//
//	double ��LM = ��L * CL.��*(SL - ��cL);
//	double ��uLS = ��L * ((SL - ��cL)*(CL.��*CL.u) + (pM - CL.p)*��x);
//	double ��vLS = ��L * ((SL - ��cL)*(CL.��*CL.v) + (pM - CL.p)*��y);
//	double eLS = ��L * ((SL - ��cL)*EL - CL.p*��cL + pM * SM);
//
//	double ��RM = ��R * CR.��*(SR - ��cR);
//	double ��uRS = ��R * ((SR - ��cR)*(CR.��*CR.u) + (pM - CR.p)*��x);
//	double ��vRS = ��R * ((SR - ��cR)*(CR.��*CR.v) + (pM - CR.p)*��y);
//	double eRS = ��R * ((SR - ��cR)*ER - CR.p*��cR + pM * SM);
//
//	double  FML[4], FMR[4];
//	FML[0] = SM * ��LM;
//	FML[1] = ��uLS * SM + pM * ��x;
//	FML[2] = ��vLS * SM + pM * ��y;
//	FML[3] = (eLS + pM)*SM;
//	FMR[0] = SM * ��RM;
//	FMR[1] = ��uRS * SM + pM * ��x;
//	FMR[2] = ��vRS * SM + pM * ��y;
//	FMR[3] = (eRS + pM)*SM;
//
//	double *F_HLLC = new double[4];
//	if (SL > 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FL[i];
//	else if (SL <= 0 && SM > 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FML[i];
//	else if (SM <= 0 && SR >= 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FMR[i];
//	else
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FR[i];
//	for (int i = 0; i < 4; i++)
//		F_HLLC[i] = F_HLLC[i] *D��;
//	return F_HLLC;
//}
//
//double * HLLC_��(mesh CD, mesh CU, mesh C, int method)
//{
//	double ��x = C.��x[method];
//	double ��y = C.��y[method];
//	double J = C.J[method];
//	double D�� = sqrt(��x*��x + ��y * ��y);
//	double ��cU = (C.��x[method] * CU.u + C.��y[method] * CD.v) / D��;
//	double ��cD = (C.��x[method] * CD.u + C.��y[method] * CU.v) / D��;
//	double aU = sqrt(��*CU.p / CU.��);
//	double aD = sqrt(��*CD.p / CD.��);
//	double EU = CU.p / (�� - 1) + 0.5*CU.�� * CU.u * CU.u + 0.5*CU.�� * CU.v * CU.v;
//	double ED = CD.p / (�� - 1) + 0.5*CD.�� * CD.u * CD.u + 0.5*CD.�� * CD.v * CD.v;
//	double GD[4], GU[4];
//
//	GD[0] = (CD.��*��cD);
//	GD[1] = (CD.��*CD.u*��cD + ��x * CD.p);
//	GD[2] = (CD.��*CD.v*��cD + ��y * CD.p);
//	GD[3] = (��cD * (ED + CD.p));
//
//	GU[0] = (CU.��*��cU);
//	GU[1] = (CU.��*CU.u*��cU + ��x * CD.p);
//	GU[2] = (CU.��*CU.v*��cU + ��y * CD.p);
//	GU[3] = (��cU * (EU + CU.p));
//
//	double SD = ��cD - aD;
//	double SU = ��cU + aU;
//	double SM = (CD.p - CU.p - CD.�� * ��x * (SD - ��cD) + CU.�� * ��x* (SU - ��cU)) / (CU.�� * (SU - ��cU) - CD.�� * (SD - ��cD));
//	double pM = 0.5*(CD.�� * (��cD - SD)*(��cD - SM) + CU.�� * (��cU - SU)*(��cU - SM) + CU.p + CD.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ��D = 1 / (SD - SM);
//	double ��U = 1 / (SU - SM);
//
//	double ��DM = ��D * CD.��*(SD - ��cD);
//	double ��uDS = ��D * ((SD - ��cD)*(CD.��*CD.u) + (pM - CD.p)*��x);
//	double ��vDS = ��D * ((SD - ��cD)*(CD.��*CD.v) + (pM - CD.p)*��y);
//	double eDS = ��D * ((SD - ��cD)*ED - CD.p*��cD + pM * SM);
//
//	double ��UM = ��U * CU.��*(SU - ��cU);
//	double ��uUS = ��U * ((SU - ��cU)*(CU.��*CU.u) + (pM - CU.p)*��x);
//	double ��vUS = ��U * ((SU - ��cU)*(CU.��*CU.v) + (pM - CU.p)*��y);
//	double eUS = ��U * ((SU - ��cU)*EU - CU.p*��cU + pM * SM);
//
//	double  FMD[4], FMU[4];
//
//	FMD[0] = SM * ��DM;
//	FMD[1] = ��uDS * SM + pM * ��x;
//	FMD[2] = ��vDS * SM + pM * ��y;
//	FMD[3] = (eDS + pM)*SM;
//	FMU[0] = SM * ��UM;
//	FMU[1] = ��uUS * SM + pM * ��x;
//	FMU[2] = ��vUS * SM + pM * ��y;
//	FMU[3] = (eUS + pM)*SM;
//
//	double *G_HLLC = new double[4];
//	if (SD >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = GD[i];
//	else if (SD <= 0 && SM >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = FMD[i];
//	else if (SM <= 0 && SU >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = FMU[i];
//	else
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = GU[i];
//	for (int i = 0; i < 4; i++)
//		G_HLLC[i] = G_HLLC[i] *D��;
//
//	return G_HLLC;
//}
