#ifndef MPSPARTICLE
#define MPSPARTICLE

#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。

class mpsparticle
{
public:
	int id;
	int firstID;
	double r[DIMENTION];
	double r0[DIMENTION];
	double u[DIMENTION];
	double n[DIMENTION];	//単位法線ベクトル
	double old_A[DIMENTION];
	double F[DIMENTION]; //粒子に働く電磁力
	double potential[DIMENTION];//表面張力項//m/s2
	


	double P;//圧力
	double h;//エンタルピー
	double PND;//reを用いた粒子数密度
	double PND2;//re2を用いた粒子数密度 
	double val;//その都度適当な値の格納に利用
	double T; //温度[K]
	double L; //基準長さ 可変解像度粒子法関係で用いる。一定解像度の場合はleと同値
	double vis;//粘性
	double ep;//相当ひずみ率
	double sigma;//相当流動応力
	
	int type;	//FLUID INWALL OUTWALL
	int materialID;	//材質番号 1ならdensity,2ならdensity2などを使用
	int surface;//0:内部 1:表面
	int index;//格納されている格子の番号
	int N;    //re内に存在する周辺粒子数
	int N2;   //re2内に存在する周辺粒子数
	int N3;   //re3内に存在する周辺粒子数
	int NEI[300];//re内に存在する周辺粒子番号
	int NEI2[300];
	int NEI3[450]; 

	int toBEM;	//FEMに転送するかどうか。ON or OFF
	int color;	//色　1なら赤2なら青　FSWでの恒久的な色付け用

	double heat_gene_before1;
	double heat_generation;

	double dir_Pst;
	double dir_Pem;
	int trace;

	int division_flag;
};


#endif