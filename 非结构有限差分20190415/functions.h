#pragma once
#include"const.h"
double distance(mesh& a, mesh& b);
void get_dt();
//void update_AfromU();
//void choose_U(int i);
//void get_F();
//void get_G();
void coordinate_trans();
Flux HLLC_��(mesh& CL, mesh& CR, mesh& C, int method);
Flux HLLC_��2(mesh& CL, mesh& CR, mesh& C, int method);//�������Ҳ��㣬����任�ο���

Flux HLLC_��(mesh& CD, mesh& CU, mesh& C, int method);
Flux HLLC_��2(mesh& CL, mesh& CR, mesh& C, int method);
void record();
void judge();
void update_IN();
void update_bound_uniform();
void update_bound_shockwave();
void update_bound_shockwaveCross();
void update_bound_Prandtl_Meyer();
double get_��(mesh A, mesh B);//��������������x��ļн�
void reorder_neighbor();
double area(mesh A, mesh B, mesh C, mesh D);//�������ĵ㹹���ı������ 
void movemesh();
void findNeiborSec();
void update_bound_shockwave_fitting();//�����߽磬����װ�䷨
void update_Vm();
void clear_Vm();
void update_bound();
double get_��(double x1, double y1, double x2, double y2);//��ֱ����x��ļн�

double max(double a, double b);
double min(double a, double b);
double absmax(double a, double b);
double absmin(double a, double b);

mesh getCrossPoint(Line L1, Line L2);
Line getLine(mesh A, mesh B);
Line getLine(double ��, mesh A);
double compute_res();//����в�
