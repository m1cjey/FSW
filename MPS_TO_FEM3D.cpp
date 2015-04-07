#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納
#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
#include"PART.h"		//class PART定義
#include"BEMclass.h"	//BEM2D関係のclass 定義
#include"FEM3Dclass.h"
#include<omp.h>
#include<vector>
#include"function.h"

//直方体解析領域作成関数
void make_cube_region(mpsconfig *CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2]);

//液滴
void MPSTOFEM3D_droplet(mpsconfig *CON,int *node_num,vector<point3D> &NODE,vector<mpsparticle> &PART, int fluid_number, int particle_number);

//MPS_TO_FEM3Dmain関数
void MPS_TO_FEM3Dmain(mpsconfig *CON,int *node_num,vector<point3D> &NODE,vector<mpsparticle> &PART, int fluid_number, int particle_number)
{
	int model=CON->get_model_number();
	if(model==3) MPSTOFEM3D_droplet(CON,node_num,NODE,PART,  fluid_number,  particle_number);//液滴

	ofstream fp("node_c.dat");
	for(int i=1;i<=*node_num;i++) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	fp.close();
	cout<<"MPSTOFEM3D完了 節点数="<<*node_num<<endl;
}

//液滴
void MPSTOFEM3D_droplet(mpsconfig *CON,int *node_num,vector<point3D> &NODE,vector<mpsparticle> &PART, int fluid_number, int particle_number)
{
	int num=0;						//節点数
	double le=CON->get_distancebp();
    double err=1e-10;
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);			//NODEは節点番号1からスタートするから、ここでひとつ確保しておく

	////流体粒子出力
    for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON || i%4==0)
		//if(PART[i].toFEM==ON)//解析中、FEMに渡す粒子番号を統一したい
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=FLUID;
			NODE[num].particleID=i;				//節点iに対応する粒子番号はi
			NODE[num].remesh=ON;			//リメッシュON
		}
    }////////////*/

	if(CON->get_region_shape()==0)
	{
		int divN[3];
		divN[A_X]=10;
		divN[A_Y]=10;
		divN[A_Z]=10;
		double regionX[2]={CON->get_XL(),CON->get_XR()};
		double regionY[2]={CON->get_YD(),CON->get_YU()};
		double regionZ[2]={CON->get_ZD(),CON->get_ZU()};
		int count=num;//現時点での節点数
		
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//解析領域は直方体

		//次に、いま作成した解析領域の少し内側に節点を追加する。これはFINE3Dで境界条件面によけいな節点が追加されることを防ぐためである。
		double dX=(regionX[1]-regionX[0])/divN[A_X];				//先ほど作成した解析領域のメッシュ長さdX
		double dY=(regionY[1]-regionY[0])/divN[A_Y];				//先ほど作成した解析領域のメッシュ長さdY
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//先ほど作成した解析領域のメッシュ長さdZ
		//以下のように定義した寸法の箱領域を作成すれば、解析領域と新しい箱との隙間には良質な正四面体に近いメッシュが作成される
		regionX[0]=CON->get_XL()+dX*sqrt(3.0)*0.5; regionX[1]=CON->get_XR()-dX*sqrt(3.0)*0.5;
		regionY[0]=CON->get_YD()+dY*sqrt(3.0)*0.5; regionY[1]=CON->get_YU()-dY*sqrt(3.0)*0.5;
		regionZ[0]=CON->get_ZD()+dZ*sqrt(3.0)*0.5; regionZ[1]=CON->get_ZU()-dZ*sqrt(3.0)*0.5;

		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//内側の箱作成

		for(int i=count+1;i<=num;i++)				//境界条件
		{
			double Z=NODE[i].r[A_Z];
			if(Z>CON->get_ZU()-err) NODE[i].boundary_condition=2;
			else if(Z<CON->get_ZD()+err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh領域作成
		regionX[0]=CON->get_XL()*0.5; regionX[1]=CON->get_XR()*0.5;
		regionY[0]=CON->get_YD()*0.5; regionY[1]=CON->get_YU()*0.5;
		regionZ[0]=CON->get_ZD()*0.7; regionZ[1]=CON->get_ZU()*0.7;
		count=num;//現時点での節点数
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//解析領域は直方体

		for(int i=count+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//境界条件
			NODE[i].remesh=ON;
		}
		
	}

	

	*node_num=num;

}

//直方体解析領域作成関数
void make_cube_region(mpsconfig *CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2])
{
	int node_num=*node;
	//divN[3];								//各辺の分割数
	
	double Xmin=regionX[0];				//解析領域
	double Xmax=regionX[1];
	double Ymin=regionY[0];
	double Ymax=regionY[1];
	double Zmin=regionZ[0];
	double Zmax=regionZ[1];

	double divL[3];								//分割幅
	divL[A_X]=(Xmax-Xmin)/divN[A_X];
	divL[A_Y]=(Ymax-Ymin)/divN[A_Y];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];

	point3D NODE01;

	//底面
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin;					//解析領域の底面
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	//上面
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmax;					//解析領域の上面
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}

	
	//側面Y
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];		
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymax;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}

	//側面
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmax;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	*node=node_num;
}
