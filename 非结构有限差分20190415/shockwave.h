#pragma once
double get_c(double ��, double p);//���ٹ�ʽ

double get_Mas(double p1, double p2);
double get_p2p1(double Ma1);
double get_��2��1(double Ma1);
double get_Ma2(double Ma1);
double get_un(mesh& A, double ��);
double get_ut(mesh& A, double ��);
void get_down(mesh& U, mesh& D, double ��);//���β���

double get_Ma(double u, double v, double ��, double p);//�������
double get_Ma2(double Ma1, double ��);//б���������������Ǯ��𢡶��������ѧ��p241,7-131
//double get_��(double Ma1, double ��);//���������
double get_��(double u, double v);
double get_��(mesh& U, double p2, int type);

double get_��from��(double Ma1, double ��);
double get_ufromMa2(double Ma2, double ��2, double p2, double ��);
double get_vfromMa2(double Ma2, double ��2, double p2, double ��);
double get_p2(double Ma1, double ��,double p1);
double get_��2(double Ma1, double ��,double ��1);
double get_Mu(mesh& U, mesh& D, double ��);
double get_udn(mesh& U, mesh& D, double Ma1, double Vs, double ��);
double get��fromMa(double Ma);
double getMafrom��(double ��);
