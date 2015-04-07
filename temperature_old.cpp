#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納
#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
#include"PART.h"		//class PART定義
#include"BEMclass.h"	//BEM2D関係のclass 定義
#include"FEM3Dclass.h"	//FEM3D関係のclass 定義
#include<omp.h>
#include<vector>
#include"function.h"

///温度場計算関数
void calc_Temperature(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t)
{
	cout<<"温度場計算------";

	//ここではｴﾝﾀﾙﾋﾟｰを主体にして計算する。すなわち、
	//Dh/Dt=kΔT+Q k:熱伝導率
	//ちなみに式変形すれば
	//DT/Dt=k/(ρC)ΔT+Q/(ρC)となる

	unsigned int timeA=GetTickCount();					//計算開始時刻
	int dim=CON->get_dimention();						//解析次元
	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;							//ラプラシアン影響半径

	double V=get_volume(CON);				//粒子の体積
	double *density=new double [particle_number];//各粒子の密度格納
	double *Cp=new double [particle_number];		//比熱
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1)
		{
			density[i]=CON->get_density();
			Cp[i]=CON->get_Cp();
		}
		else if(PART[i].materialID==2) 
		{
			density[i]=CON->get_density2();
			Cp[i]=CON->get_Cp2();
		}
	}
	for(int i=fluid_number;i<particle_number;i++)
	{
		density[i]=CON->get_wall_density();
		Cp[i]=CON->get_wall_Cp();
	}
	double *mass=new double [fluid_number];	//各粒子の質量格納
	for(int i=0;i<fluid_number;i++) mass[i]=density[i]*V;
		
	double *MP=new double [fluid_number];		//融点
	double *latent_H=new double [fluid_number];		//潜熱
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1)
		{	
			MP[i]=CON->get_MP();
			latent_H[i]=CON->get_latent_H();
		}
		else if(PART[i].materialID==2)
		{
			MP[i]=CON->get_MP2();
			latent_H[i]=CON->get_latent_H2();
		}
	}

	double roomT=CON->get_roomT();						//室温[K]
	double wall_mass=V*CON->get_wall_density();			//壁粒子の質量

    double *T_laplacian = new double [particle_number];
    double *T=new double [particle_number];//温度
    double *k=new double [particle_number];//各粒子の熱伝導率
   
    ///エンタルピーから各粒子の温度T[i]を求める
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//融解開始点のエンタルピー
		double hs1=hs0+latent_H[i]*mass[i];			//融解終了点のエンタルピー

        if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];			//固体
		else if(hs0<=PART[i].h && PART[i].h<=hs1) T[i]=MP[i];	//融点
		else if(hs1<PART[i].h)											//液体
		{       
			T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];
			if( PART[i].type==SOLID ||PART[i].type==BOSOLID ||PART[i].type>CFD)//相変態
			{
				PART[i].type=FLUID;
				if(PART[i].PND>n0*CON->get_beta()) PART[i].surface=OFF;
				else PART[i].surface=ON;
			}  
		}
    }
	for(int i=fluid_number;i<particle_number;i++)//壁粒子
    {
        T[i]=PART[i].h/wall_mass/CON->get_wall_Cp();//固体
    }
    ///////////////
	
    ///////各粒子の熱伝導率k[i]を求める
	double maxk=0;//最大熱伝導率
	double densityCp=CON->get_density()*CON->get_Cp();
    for(int i=0;i<fluid_number;i++)
    {
		if(PART[i].materialID==1) k[i]=CON->get_k();
		else if(PART[i].materialID==2) k[i]=CON->get_k2();
		if(k[i]>maxk) maxk=k[i];
    }
	for(int i=fluid_number;i<particle_number;i++)
    {
		k[i]=CON->get_wall_k();//壁の熱伝導率
		if(k[i]>maxk) maxk=k[i];
    }
	maxk/=densityCp;//熱拡散率
	if(maxk*dt/(le*le)>0.2) cout<<"熱伝導率に関してdtが大きすぎる dt<"<<0.2*le*le/maxk<<"にしてください"<<endl;

    ////温度の拡散計算
    if(CON->get_insulate()==0)//断熱
    {
        for(int i=0;i<particle_number;i++)
        {
            T_laplacian[i]=0;//初期化
			if(PART[i].type!=INWALL && PART[i].type!=OUTWALL)
			{
				double lam=0;//正確なλ
				double W=0;//粒子数密度
				for(int k1=0;k1<PART[i].N2;k1++)
				{       
					int j=PART[i].NEI2[k1]; 
					if(PART[j].type!=INWALL && PART[j].type!=OUTWALL)
					{ 
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);//粒子間の距離
						double w=kernel(R,dis);//重み関数
						W+=w;
						lam+=dis*dis*w;
						double dT=T[j]-T[i];//温度差
						if(CON->get_T_laplacian()==0 || CON->get_T_laplacian()==1)
						{
							T_laplacian[i]+=dT*w*k[i]*k[j]/(k[i]+k[j])*2;
						}
						else if(CON->get_T_laplacian()==2)
						{ 
							T_laplacian[i]+=dT*w/(dis*dis)*k[i]*k[j]/(k[i]+k[j])*2;
						}  
					}
				} 
	    
				if(W!=0) lam/=W;
				else if(W==0) lam=lamda;
				int d=CON->get_dimention();
				if(CON->get_T_laplacian()==0)  T_laplacian[i]=T_laplacian[i]*2*d/(n0*lamda);
				else if(CON->get_T_laplacian()==1 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/(W*lam);
				else if(CON->get_T_laplacian()==2 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/W;
			}
		}
    }
    else if(CON->get_insulate()==1)//非断熱
    {
        for(int i=0;i<particle_number;i++)
        {
            T_laplacian[i]=0;//初期化
			//if(PART[i].type!=INWALL && PART[i].type!=OUTWALL)
			{
				double lam=0;//正確なλ
				double W=0;//粒子数密度
				for(int k1=0;k1<PART[i].N2;k1++)
				{       
					int j=PART[i].NEI2[k1]; 
		
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);//粒子間の距離
					double w=kernel(R,dis);//重み関数
					W+=w;
					lam+=dis*dis*w;
					double dT=T[j]-T[i];//温度差

					if(CON->get_T_laplacian()==0 || CON->get_T_laplacian()==1)
					{
						T_laplacian[i]+=dT*w*k[i]*k[j]/(k[i]+k[j])*2;///(density[i]*Cp[i]);
					}
					else if(CON->get_T_laplacian()==2)
					{ 
						T_laplacian[i]+=dT*w/(dis*dis)*k[i]*k[j]/(k[i]+k[j])*2;
					}
				}
				if(W!=0) lam/=W;
				else if(W==0) lam=lamda;
				int d=CON->get_dimention();
				if(CON->get_T_laplacian()==0)  T_laplacian[i]=T_laplacian[i]*2*d/(n0*lamda);
				else if(CON->get_T_laplacian()==1 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/(W*lam);
				else if(CON->get_T_laplacian()==2 &&W!=0) T_laplacian[i]=T_laplacian[i]*2*d/W;
				
			}
		}
    }
    //////////////*/
  
    for(int i=0;i<particle_number;i++) PART[i].h+=T_laplacian[i]*dt*V;//エンタルピー更新

	if(CON->get_air_cool()==ON)//空気との熱伝達を考慮する場合
	{
		double h=50;			//空気の熱伝達率　単位は[W/m^2/K]. 値はきちんと調べること。ヌッセル数とかと関係あり？
		double S;//粒子の断面積
		if(dim==2) S=le;
		else
		{
			if(CON->get_model_set_way()==0) S=le*le;//正方格子の場合
			else if(CON->get_model_set_way()==1) S=sqrt(3.0)/2*le*le;//細密格子の場合
		}

		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].surface==ON) PART[i].h+=h*(roomT-T[i])*S*dt;
		}
	}//////////*/
   
    
    ////温度プロット
	double height=CON->get_TplotZ();
    plot_T(CON ,PART,particle_number,T,height);
	if(CON->get_T_AVS()>0)
	{
		if(t%CON->get_T_AVS()==0 || t==1)output_temperature_avs(CON,PART, t, particle_number, fluid_number,T, height);
	}


	if(CON->get_dimention()==3)//XZ平面のTも出力
	{
		ofstream fh("T_XZ.dat");
		for(int i=0;i<particle_number;i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fh<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<T[i]<<endl;
		fh.close();
	}
    
    delete [] T;
    delete [] k;
    delete [] T_laplacian;

	delete [] density;
	delete [] mass;
	delete [] Cp;
	delete [] MP;
	delete [] latent_H;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//温度プロット関数
void plot_T(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,double *T,double height)
{
	ofstream t("T_XY.dat");
	double le=CON->get_distancebp();
	double L=le*0.5;

	if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{
			double x=PART[i].r[A_X];
			double y=PART[i].r[A_Y];
			double Z=T[i];
			t<<x<<" "<<y<<" "<<Z<<endl;
		}
	}
	else if(CON->get_dimention()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			double z=PART[i].r[A_Z];
			if(z>height-L && z<height+L)
			{
				double x=PART[i].r[A_X];
				double y=PART[i].r[A_Y];
				double Z=T[i];
				t<<x<<" "<<y<<" "<<Z<<endl;
			}
		}
	}

	t.close();
}

//温度AVSファイル出力関数
void output_temperature_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number,double *T,double height)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	t=1;//いまはわざと毎ステップ上書き

	//sprintf(filename,"pressure/pressure%d",t);//フォルダを作成して管理する場合はこちら
	sprintf(filename,"T_XZ%d",t);//他のファイルと同じ階層に生成するならこちら
	ofstream fout(filename);
	if(!fout)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n++;
			}
		}
	}
	else if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{
			double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
			double y=PART[i].r[A_Y]*1.0E+05;
			double z=PART[i].r[A_Z]*1.0E+05;
			double P=T[i];
			fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
			n++;
		}
	}
	fout.close();
	//sprintf(filename,"pressure/pressure%d.fld",t);//フォルダを作成して管理する場合はこちら
	sprintf(filename,"T_XZ%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	ofstream fout2(filename);
	if(!fout2)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	fout2 << "dim1=" << n <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=1" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=temperature" << endl << endl;
	//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout2 << "variable 1 file=T_XZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    1 file=T_XZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    2 file=T_XZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    3 file=T_XZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2.close();
}
