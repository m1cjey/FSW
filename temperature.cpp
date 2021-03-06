#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
//#include"define.h"	//#define 格納
//#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
//#include"PART.h"		//class PART定義
#include"BEMclass.h"	//BEM2D関係のclass 定義
//#include"FEM3Dclass.h"	//FEM3D関係のclass 定義
//#include<omp.h>
//#include<vector>
#include"function.h"

///温度場計算関数
void calc_Temperature(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t)
{
	int counth=0;
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].heat_generation>0) counth++;
	}
	//cout<<"heatが0でない粒子数="<<counth<<endl;
	//for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) PART[i].heat_generation=10;
	cout<<"温度場計算------";

	//ここではｴﾝﾀﾙﾋﾟｰを主体にして計算する。すなわち、
	//Dh/Dt=kΔT+Q k:熱伝導率
	//ちなみに式変形すれば
	//DT/Dt=k/(ρC)ΔT+Q/(ρC)となる

	unsigned int timeA=GetTickCount();					//計算開始時刻
	int dim=CON->get_dimention();						//解析次元
	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;							//ラプラシアン影響半径

	double V=CON->get_particle_volume();			//粒子の体積
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


	
	//発熱量を計算する
	for(int i=0;i<particle_number;i++)
	{
		PART[i].h+=0.5*dt*(PART[i].heat_generation+PART[i].heat_gene_before1)*V;		//1ステップ前の発熱と合わせて台形則で近似
		//PART[i].h+=0.5*dt*(PART[i].heat_generation+PART[i].heat_gene_before1)*1;//*V;		//1ステップ前の発熱と合わせて台形則で近似
	}

    double *T_laplacian = new double [particle_number];
    double *T=new double [particle_number];//温度
	int *insulate=new int [particle_number];//壁粒子が断熱されているかどうか　局所的に断熱状態にしたいときなどに用いる
    double *k=new double [particle_number];//各粒子の熱伝導率

	for(int i=0;i<particle_number;i++) insulate[i]=OFF;
	
   
    ///エンタルピーから各粒子の温度T[i]を求める
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//融解開始点のエンタルピー
		double hs1=hs0+latent_H[i]*mass[i];			//融解終了点のエンタルピー

		if(PART[i].h<hs0)
		{
			//cout<<"固体"<<endl;
			T[i]=PART[i].h/mass[i]/Cp[i];			//固体
		}
		else if(hs0<=PART[i].h && PART[i].h<hs1)
		{
			//cout<<"融解中?"<<endl;
			T[i]=MP[i];	//融点
			
		}
		else if(hs1<=PART[i].h)											//液体
		{
			//cout<<"液体"<<endl;
			T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];
			if(CON->get_material()==H2O && T[i]>=393) T[i]=393;//現在沸騰は考えないので、もし温度が100度越えすれば100度に戻す(気化熱の実装の方が先決か?）
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
		if(CON->get_model_number()==26)//釜の場合、壁の場所に応じて物性地をかえる
		{
			double air_mass=V*1.205;//空気の重さ
			double air_Cp=1006;//空気の比熱
			if(CON->get_dimention()==2)
			{
				if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08)//上壁
				{
					T[i]=PART[i].h/air_mass/air_Cp;//固体
				}
				else T[i]=PART[i].h/wall_mass/CON->get_wall_Cp();//固体
			}
			else if(CON->get_dimention()==3)
			{
			}
		}
        else T[i]=PART[i].h/wall_mass/CON->get_wall_Cp();//固体
    }
    ///////////////

	////fswにおいて、底部以外の壁面の温度条件を変える
	if(CON->get_model_number()==19)
	{
		double width0=18*1e-3;//25*1e-3;//18*1e-3;//流体が占める幅
		double Width=width0;		//横18mm+壁粒子を左右に4粒子分
		double Height=6*1e-3;	//高さ6mm+壁粒子を下に4粒子分
		double Depth=36e-3;//18e-3;//30*1e-3;//27*1e-3+6*le*A;	//奥行き27mm+壁粒子を左右に4粒子分
		double depth0=9e-3;				//ツール中心と、手前の壁との距離
	
		//double Z_mod=le*B+(1-sqrt(2.0)/2)*le;//ツールが流体へめり込むのを防ぐための調整量　この値だけツール以外の物体が下がる
		double Z_mod=le+(1-sqrt(2.0)/2)*le;
			
		///////////////////box粒子の材質決定

		for(int i=fluid_number;i<particle_number;i++) 

		{
			if(PART[i].r[A_X]>(Width*0.5)+0.2*le)//右の壁
			{
				T[i]=CON->get_roomT();
				insulate[i]=OFF;
			}
			else if(PART[i].r[A_X]<-(Width*0.5)-0.2*le)//左の壁
			{
				T[i]=CON->get_roomT();
				insulate[i]=OFF;
			}
				
			else if(PART[i].r[A_Y]<-(depth0)-0.2*le)//手前の壁
			{
				T[i]=CON->get_roomT();
				insulate[i]=OFF;
			}
			else if(PART[i].r[A_Y]>(Depth-depth0)+0.2*le)//奥の壁
			{
					T[i]=CON->get_roomT();
					insulate[i]=OFF;
			}
			if(PART[i].r[A_X]<=(Width*0.5)+0.2*le && PART[i].r[A_X]>=-(Width*0.5)-0.2*le)
			{
				if(PART[i].r[A_Y]>=-(depth0)-0.2*le && PART[i].r[A_Y]<=(Depth-depth0)+0.2*le)//下の壁
				{
					if(PART[i].r[A_Z]<-0.2*le-Z_mod)
					{
						insulate[i]=ON;	
					}
				}
			}
		}
	}

	//境界の温度そのものが実験などで既にもとまっている場合、それを対応する境界に位置する粒子に与える
	//set_temperature_boundary(CON, PART, fluid_number, particle_number, T,dt,t);
	
	//////////////////////*/
	
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
		if(CON->get_model_number()==26)//釜の場合、壁の場所に応じて物性地をかえる
		{
			if(CON->get_dimention()==2)
			{
				if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08)//上壁
				{
					k[i]=0.0257;//空気の熱伝導率
				}
				else k[i]=CON->get_wall_k();//壁の熱伝導率
			}
			else if(CON->get_dimention()==3)
			{
			}
		}
    }

	maxk/=densityCp;//熱拡散率
	if(maxk*dt/(le*le)>0.2) cout<<"熱伝導率に関してdtが大きすぎる dt<"<<0.2*le*le/maxk<<"にしてください"<<endl;	//理論的限界は0.5


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
	else if(CON->get_insulate()==2)//局所的に断熱
    {
        for(int i=0;i<particle_number;i++)
        {
            T_laplacian[i]=0;//初期化
			if(insulate[i]==OFF)
			{
				double lam=0;//正確なλ
				double W=0;//粒子数密度
				for(int k1=0;k1<PART[i].N2;k1++)
				{       
					int j=PART[i].NEI2[k1]; 
					if(insulate[j]==OFF)
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
							T_laplacian[i]+=dT*w*k[i]*k[j]/(k[i]+k[j])*2;///(density[i]*Cp[i]);
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
    //////////////*/
  
	//エンタルピー更新
    for(int i=0;i<particle_number;i++) PART[i].h+=T_laplacian[i]*dt*V;

	if(CON->get_air_cool()==ON)//空気との熱伝達を考慮する場合
	{
		double h=5;//50;			//空気の熱伝達率　単位は[W/m^2/K]. 値はきちんと調べること。ヌッセル数とかと関係あり？
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
			else if(CON->get_model_number()==26 && PART[i].heat_generation>0) PART[i].h+=h*(roomT-T[i])*S*dt;//釜において、発熱量が存在する(=釜の最外部)も罰熱があるとする
		}
	}//////////*/
    
	double Tmax=0;
	for(int i=0;i<particle_number;i++) if(T[i]>Tmax) Tmax=T[i];
	double Tminw=100000;
	for(int i=0;i<particle_number;i++) if(T[i]<Tminw) if(PART[i].type!=FLUID) Tminw=T[i];
	//cout<<"最大温度="<<Tmax<<"[K]"<<endl;
	double Tmaxf=0;
	for(int i=0;i<fluid_number;i++) if(T[i]>Tmaxf) Tmaxf=T[i];
	cout<<"最大温度="<<Tmax<<"[K]"<<" 流体最大温度"<<Tmaxf<<"[K]"<<endl;

	//発熱量を初期化
	for(int i=0;i<particle_number;i++)
	{
		PART[i].heat_gene_before1=PART[i].heat_generation;	//1step前の情報へと格納しなおし
		PART[i].heat_generation=0;							//初期化　発熱量を増やすときは+=の記述にすると、いろんな関数内で発熱が起こったときにも対処できる//現在は渦電流損のみで、舞ステップ求めてい
	}
	
	//粒子へ温度を記憶
	for(int i=0;i<particle_number;i++)
	{
		//if(PART[i].type==FLUID) PART[i].T=T[i];
		 PART[i].T=T[i];
		//if(PART[i].T<CON->get_initialT()) cout<<"初期温度より低い?"<<endl; 
	}

	////温度プロット
	double height=CON->get_TplotZ();
    plot_T(CON ,PART,particle_number,T,height);
	if(CON->get_T_AVS()>0)
	{
		//if(t%CON->get_T_AVS()==0 || t==1)output_temperature_avs(CON,PART, t, particle_number, fluid_number,T, height);
		//if(t%CON->get_T_AVS()==0 || t==1)output_temperature_avs2(CON,PART, t, particle_number, fluid_number,T, height);//06/03
	}

	///microavs用_全部出力
	/*
	if(t==1 || t%10==0)
	{
		int n=0;
		char filename[30];
		sprintf_s(filename,"PART.T%d",t);//他のファイルと同じ階層に生成するならこちら
		ofstream fout5(filename);
		if(!fout5)
		{
			cout << "cannot open" << filename << endl;
			exit(EXIT_FAILURE);
		}
		if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type!=FLUID)
				{
					double x=PART[i].r[A_X];	//rは非常に小さい値なので10^5倍しておく
					double y=PART[i].r[A_Y];//*1.0E+05;
					double z=PART[i].r[A_Z];//*1.0E+05;
					double T=PART[i].T;
					fout5 << T << "\t" << x << "\t" << y << "\t" << z << endl;
					//fout5 << insulate[i] << "\t" << x << "\t" << y << "\t" << z << endl; こっちだと、断熱、非断熱の条件設定を可視化できる
					n++;
				}
			}
		}
		fout5.close();
		
		sprintf_s(filename,"PART.T%d.fld",t);//他のファイルと同じ階層に生成するならこちら
		ofstream fout6(filename);
		if(!fout6)
		{
			cout << "cannot open" << filename << endl;
			exit(EXIT_FAILURE);
		}

		fout6 << "# AVS field file" << endl;
		fout6 << "ndim=1" << endl;
		fout6 << "dim1=" << n <<endl;
		fout6 << "nspace=3" << endl;
		fout6 << "veclen=1" << endl;
		fout6 << "data=float" << endl;
		fout6 << "field=irregular" << endl;
		fout6 << "label=T" << endl << endl;
		fout6 << "variable 1 file=PART.T" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout6 << "coord    1 file=PART.T" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout6 << "coord    2 file=PART.T" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout6 << "coord    3 file=PART.T" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout6.close();
	}*/	//05/07
	
	if(CON->get_dimention()==3)//XZ平面のTも出力
	{
		
		ofstream fh("T_XZ.dat");
		for(int i=0;i<particle_number;i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fh<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<T[i]<<endl;
		fh.close();
	}
    
    delete [] T;
	delete [] insulate;
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
	int flag_out_f=0;
	int flag_out_b=0;
	char filename[30];
	char filename_n[30];
	char filename_f[30];
	char filename_b[30];
	int n=0,nn=0,nf=0,nb=0;

	double le=CON->get_distancebp();
	double cross_section=CON->get_speed_face_p();
	int output_face=CON->get_speed_face();
	int output_face_n=0;

	double shold_R=4.25*1e-3;
	if(CON->get_tool_angle()>0)	
	{
		double	angle=CON->get_tool_angle();
		shold_R*=cos(angle);
	}

	//出力面の移動
	if(CON->get_process_type()==1||CON->get_process_type()==2)
	{
		double dt=CON->get_dt();
		int dwell_step=CON->get_dwelling_time()/dt;
		int change_step=CON->get_change_step();
		double speed2=CON->get_move_speed2();

		if(t>dwell_step+change_step)	output_face+=speed2*dt*(t-dwell_step-change_step);
	}


	//FSW方向転換
	if(CON->get_process_type()==2||CON->get_process_type()==0)
	{
		double probe_H=4*1e-3;
		double t_base=probe_H/CON->get_move_speed();		
		double TIME=CON->get_dt()*(t-1);
		double t_dw=CON->get_dwelling_time();

		if(TIME>=t_base)
		{
			if(CON->get_output_forward()==ON)	flag_out_f=ON;
			if(CON->get_output_backward()==ON)	flag_out_b=ON;
		}
	}
	else
	{
			if(CON->get_output_forward()==ON)	flag_out_f=ON;
			if(CON->get_output_backward()==ON)	flag_out_b=ON;
	}			
	//t=1;//いまはわざと毎ステップ上書き
	
	//sprintf_s(filename,"pressure/pressure%d",t);//フォルダを作成して管理する場合はこちら
	if(output_face==0)
	{
		sprintf_s(filename,"T_YZ%d",t);//他のファイルと同じ階層に生成するならこちら
		if(CON->get_output_another_face()==ON)
		{
			sprintf_s(filename_n,"T_XZ%d",t);
			output_face_n=1;
		}
		if(flag_out_f==ON)	sprintf_s(filename_f,"T_YZ_forward%d",t);
		if(flag_out_b==ON)	sprintf_s(filename_b,"T_YZ_backward%d",t);
	}
	else if(output_face==1)
	{
	//	sprintf_s(filename,"T_XZ%d",t);
		if(CON->get_output_another_face()==ON)
		{
			sprintf_s(filename_n,"T_YZ%d",t);
			output_face_n=0;
		}
		if(flag_out_f==ON)	sprintf_s(filename_f,"T_XZ_forward%d",t);
		if(flag_out_b==ON)	sprintf_s(filename_b,"T_XZ_backward%d",t);
	}
	else if(output_face==2)
	{
		sprintf_s(filename,"T_XY%d",t);
		if(flag_out_f==ON)	sprintf_s(filename_f,"T_XY_forward%d",t);
		if(flag_out_b==ON)	sprintf_s(filename_b,"T_XY_backward%d",t);
	}

//	ofstream fout(filename);
	ofstream	fout_n(filename_n);
	ofstream	fout_f(filename_f);
	ofstream	fout_b(filename_b);
/*	if(!fout)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}*/

	if(CON->get_dimention()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type==FLUID)
			{
/*				if(PART[i].r[output_face]<cross_section+le && PART[i].r[output_face]>cross_section-le)	
				//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
				{
					double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
					double y=PART[i].r[A_Y]*1.0E+05;
					double z=PART[i].r[A_Z]*1.0E+05;
					double P=T[i];
					fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}*/

				if(CON->get_output_another_face()==ON)
				{
					if(PART[i].r[output_face_n]<cross_section+le && PART[i].r[output_face_n]>cross_section-le)	
					//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
					{
						double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
						double y=PART[i].r[A_Y]*1.0E+05;
						double z=PART[i].r[A_Z]*1.0E+05;
						double P=T[i];
						fout_n << P << "\t" << x << "\t" << y << "\t" << z << endl;
						nn++;
					}
				}

				if(flag_out_f==ON)
				{
					if(PART[i].r[output_face]<cross_section-shold_R+le && PART[i].r[output_face]>cross_section-shold_R-le)	
					//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
					{
						double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
						double y=PART[i].r[A_Y]*1.0E+05;
						double z=PART[i].r[A_Z]*1.0E+05;
						double P=T[i];
						fout_f << P << "\t" << x << "\t" << y << "\t" << z << endl;
						nf++;
					}
				}

				if(flag_out_b==ON)
				{
					if(PART[i].r[output_face]<cross_section+shold_R+le && PART[i].r[output_face]>cross_section+shold_R-le)	
					//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
					{
						double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
						double y=PART[i].r[A_Y]*1.0E+05;
						double z=PART[i].r[A_Z]*1.0E+05;
						double P=T[i];
						fout_b << P << "\t" << x << "\t" << y << "\t" << z << endl;
						nb++;
					}
				}
			}
		}
	}
/*	else if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type!=FLUID)
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;

				//double x=PART[i].r[A_X];//*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				//double y=PART[i].r[A_Y];//*1.0E+05;
				//double z=PART[i].r[A_Z];//*1.0E+05;
				double P=T[i];
				//double P=PART[i].heat_gene_before1;
				fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n++;
			}
		}
	}*/
//	fout.close();
	fout_n.close();
	fout_b.close();
	fout_f.close();


	//sprintf_s(filename,"pressure/pressure%d.fld",t);//フォルダを作成して管理する場合はこちら
	if(output_face==0)
	{
		sprintf_s(filename,"T_YZ%d.fld",t);//他のファイルと同じ階層に生成するならこちら
		if(CON->get_output_another_face()==ON)	sprintf_s(filename_n,"T_XZ%d.fld",t);
		if(flag_out_f==ON)	sprintf_s(filename_f,"T_YZ_forward%d.fld",t);
		if(flag_out_b==ON)	sprintf_s(filename_b,"T_YZ_backward%d.fld",t);
	}
	else if(output_face==1)
	{
	//	sprintf_s(filename,"T_XZ%d.fld",t);
		if(CON->get_output_another_face()==ON)	sprintf_s(filename_n,"T_YZ%d.fld",t);
		if(flag_out_f==ON)	sprintf_s(filename_f,"T_XZ_forward%d.fld",t);
		if(flag_out_b==ON)	sprintf_s(filename_b,"T_XZ_backward%d.fld",t);
	}
	else if(output_face==2)
	{
		sprintf_s(filename,"T_XY%d.fld",t);
		if(flag_out_f==ON)	sprintf_s(filename_f,"T_XY_forward%d.fld",t);
		if(flag_out_b==ON)	sprintf_s(filename_b,"T_XY_backward%d.fld",t);
	}

//	ofstream fout2(filename);
	ofstream fout_n2(filename_n);
	ofstream fout_f2(filename_f);
	ofstream fout_b2(filename_b);
/*	if(!fout2)
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

	if(output_face==0)
	{
		fout2 << "variable 1 file=T_YZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    1 file=T_YZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    2 file=T_YZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    3 file=T_YZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	}
	if(output_face==1)
	{
		fout2 << "variable 1 file=T_XZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    1 file=T_XZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    2 file=T_XZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    3 file=T_XZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	}
	if(output_face==2)
	{
		fout2 << "variable 1 file=T_XY" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    1 file=T_XY" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    2 file=T_XY" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		fout2 << "coord    3 file=T_XY" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	}
	fout2.close();*/

	///////////////////温度分布出力_他断面
	if(CON->get_output_another_face()==ON)
	{
		fout_n2 << "# AVS field file" << endl;
		fout_n2 << "ndim=1" << endl;
		fout_n2 << "dim1=" << nn <<endl;
		fout_n2 << "nspace=3" << endl;
		fout_n2 << "veclen=1" << endl;
		fout_n2 << "data=float" << endl;
		fout_n2 << "field=irregular" << endl;
		fout_n2 << "label=temperature" << endl << endl;
		//fout2n << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2n << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2n << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2n << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		if(output_face_n==0)
		{
			fout_n2 << "variable 1 file=T_YZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_n2 << "coord    1 file=T_YZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_n2 << "coord    2 file=T_YZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_n2 << "coord    3 file=T_YZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		if(output_face_n==1)
		{
			fout_n2 << "variable 1 file=T_XZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_n2 << "coord    1 file=T_XZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_n2 << "coord    2 file=T_XZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_n2 << "coord    3 file=T_XZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		fout_n2.close();
	}

	////////////////粘性分布出力_前方
	if(flag_out_f==ON)
	{
		fout_f2 << "# AVS field file" << endl;
		fout_f2 << "ndim=1" << endl;
		fout_f2 << "dim1=" << nf <<endl;
		fout_f2 << "nspace=3" << endl;
		fout_f2 << "veclen=1" << endl;
		fout_f2 << "data=float" << endl;
		fout_f2 << "field=irregular" << endl;
		fout_f2 << "label=temperature" << endl << endl;
		//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		if(output_face==0)
		{
			fout_f2 << "variable 1 file=T_YZ_forward" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    1 file=T_YZ_forward" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    2 file=T_YZ_forward" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    3 file=T_YZ_forward" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		else if(output_face==1)
		{
			fout_f2 << "variable 1 file=T_XZ_forward" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    1 file=T_XZ_forward" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    2 file=T_XZ_forward" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    3 file=T_XZ_forward" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		else if(output_face==2)
		{
			fout_f2 << "variable 1 file=T_XY_forward" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    1 file=T_XY_forward" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    2 file=T_XY_forward" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_f2 << "coord    3 file=T_XY_forward" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		fout_f2.close();
	}

	////////////温度分布出力＿後方
	if(flag_out_b==ON)
	{
		fout_b2 << "# AVS field file" << endl;
		fout_b2 << "ndim=1" << endl;
		fout_b2 << "dim1=" << nb <<endl;
		fout_b2 << "nspace=3" << endl;
		fout_b2 << "veclen=1" << endl;
		fout_b2 << "data=float" << endl;
		fout_b2 << "field=irregular" << endl;
		fout_b2 << "label=temperature" << endl << endl;
		//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
		if(output_face==0)
		{
			fout_b2 << "variable 1 file=T_YZ_backward" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    1 file=T_YZ_backward" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    2 file=T_YZ_backward" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    3 file=T_YZ_backward" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		else if(output_face==1)
		{
			fout_b2 << "variable 1 file=T_XZ_backward" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    1 file=T_XZ_backward" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    2 file=T_XZ_backward" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    3 file=T_XZ_backward" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		else if(output_face==2)
		{
			fout_b2 << "variable 1 file=T_XY_backward" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    1 file=T_XY_backward" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    2 file=T_XY_backward" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
			fout_b2 << "coord    3 file=T_XY_backward" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
		}
		fout_b2.close();
	}

}

void output_temperature_avs2(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number,double *T,double height)
{

	char filename2[40];

	double le=CON->get_distancebp();
	double r_face=CON->get_speed_face_p();

	double shold_R=4.25*1e-3;
	if(CON->get_tool_angle()>0)	
	{
		double	angle=CON->get_tool_angle();
		shold_R*=cos(angle);
	}

	//出力面の移動
	if(CON->get_process_type()==1||CON->get_process_type()==2)
	{
		double dt=CON->get_dt();
		int dwell_step=CON->get_dwelling_time()/dt;
		int change_step=CON->get_change_step();
		double speed2=CON->get_move_speed2();

		if(t>dwell_step+change_step)	r_face+=speed2*dt*(t-dwell_step-change_step);
	}

	int YZ=0,XZ=1,XY=2;
	int face=XZ,face2=XY,face3=YZ;
	//t=1;//いまはわざと毎ステップ上書き
	
	//sprintf_s(filename,"pressure/pressure%d",t);//フォルダを作成して管理する場合はこちら
	sprintf_s(filename2,"./Temperature/T_XY%d",t);//他のファイルと同じ階層に生成するならこちら

	ofstream	fout2(filename2);

	if(!fout2)
	{
		cout << "cannot open" << filename2 << endl;
		exit(EXIT_FAILURE);
	}

	int n=0,n2=0,n3=0,n4=0;
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID)
		{

			if(PART[i].r[face2]<6.0*1e-3-shold_R*sin(3.0)-le && PART[i].r[face2]>6.0*1e-3-shold_R*sin(3.0)-3*le)	
			//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout2 << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n2++;
			}

		}
	}
	fout2.close();
	sprintf_s(filename2,"./Temperature/T_XY%d.fld",t);//他のファイルと同じ階層に生成するならこちら

	ofstream fout2_fld(filename2);
	if(!fout2_fld)
	{
		cout << "cannot open" << filename2 << endl;
		exit(EXIT_FAILURE);
	}

	fout2_fld << "# AVS field file" << endl;
	fout2_fld << "ndim=1" << endl;
	fout2_fld << "dim1=" << n2 <<endl;
	fout2_fld << "nspace=3" << endl;
	fout2_fld << "veclen=1" << endl;
	fout2_fld << "data=float" << endl;
	fout2_fld << "field=irregular" << endl;
	fout2_fld << "label=temperature" << endl << endl;
	//fout2_fld << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2_fld << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2_fld << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2_fld << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout2_fld << "variable 1 file=T_XY" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld << "coord    1 file=T_XY" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld << "coord    2 file=T_XY" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld << "coord    3 file=T_XY" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld.close();


	char filename3[40];
	char filename4[40];

	sprintf_s(filename3,"./Temperature/T_XZ_back%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename4,"./Temperature/T_YZ%d",t);//他のファイルと同じ階層に生成するならこちら

	ofstream	fout3(filename3);
	ofstream	fout4(filename4);

	if(!fout3)
	{
		cout << "cannot open" << filename3 << endl;
		exit(EXIT_FAILURE);
	}
	if(!fout4)
	{
		cout << "cannot open" << filename4 << endl;
		exit(EXIT_FAILURE);
	}

	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID)
		{
			if(PART[i].r[face]<r_face-shold_R+le && PART[i].r[face]>r_face-shold_R-le)	
			//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout3 << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n3++;
			}
			if(PART[i].r[face3]<le && PART[i].r[face3]>-le)	
			//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout4 << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n4++;
			}
		}
	}
	fout3.close();
	fout4.close();
	sprintf_s(filename3,"./Temperature/T_XZ_back%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename4,"./Temperature/T_YZ%d.fld",t);//他のファイルと同じ階層に生成するならこちら

	
	ofstream fout3_fld(filename3);
	if(!fout3_fld)
	{
		cout << "cannot open" << filename3 << endl;
		exit(EXIT_FAILURE);
	}

	fout3_fld << "# AVS field file" << endl;
	fout3_fld << "ndim=1" << endl;
	fout3_fld << "dim1=" << n3 <<endl;
	fout3_fld << "nspace=3" << endl;
	fout3_fld << "veclen=1" << endl;
	fout3_fld << "data=float" << endl;
	fout3_fld << "field=irregular" << endl;
	fout3_fld << "label=temperature" << endl << endl;
	//fout_fld << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout3_fld << "variable 1 file=T_XZ_back" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld << "coord    1 file=T_XZ_back" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld << "coord    2 file=T_XZ_back" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld << "coord    3 file=T_XZ_back" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld.close();

	ofstream fout4_fld(filename4);
	if(!fout4_fld)
	{
		cout << "cannot open" << filename4 << endl;
		exit(EXIT_FAILURE);
	}
	fout4_fld << "# AVS field file" << endl;
	fout4_fld << "ndim=1" << endl;
	fout4_fld << "dim1=" << n4 <<endl;
	fout4_fld << "nspace=3" << endl;
	fout4_fld << "veclen=1" << endl;
	fout4_fld << "data=float" << endl;
	fout4_fld << "field=irregular" << endl;
	fout4_fld << "label=temperature" << endl << endl;
	//fout4_fld << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout4_fld << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout4_fld << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout4_fld << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout4_fld << "variable 1 file=T_YZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld << "coord    1 file=T_YZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld << "coord    2 file=T_YZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld << "coord    3 file=T_YZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld.close();
	
	/*
	char filename[40];
	char filename2[40];
	char filename3[40];
	char filename4[40];

	double le=CON->get_distancebp();
	double r_face=CON->get_speed_face_p();


	double shold_R=4.25*1e-3;
	if(CON->get_tool_angle()>0)	
	{
		double	angle=CON->get_tool_angle();
		shold_R*=cos(angle);
	}

	//出力面の移動
	if(CON->get_process_type()==1||CON->get_process_type()==2)
	{
		double dt=CON->get_dt();
		int dwell_step=CON->get_dwelling_time()/dt;
		int change_step=CON->get_change_step();
		double speed2=CON->get_move_speed2();

		if(t>dwell_step+change_step)	r_face+=speed2*dt*(t-dwell_step-change_step);
	}




	int YZ=0,XZ=1,XY=2;
	int face=XZ,face2=XY,face3=YZ;
	//t=1;//いまはわざと毎ステップ上書き
	
	//sprintf_s(filename,"pressure/pressure%d",t);//フォルダを作成して管理する場合はこちら
	sprintf_s(filename,"./Temperature/T_XZ%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename2,"./Temperature/T_XY%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename3,"./Temperature/T_XZ_back%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename4,"./Temperature/T_YZ%d",t);//他のファイルと同じ階層に生成するならこちら

	ofstream fout(filename);
	ofstream	fout2(filename2);
	ofstream	fout3(filename3);
	ofstream	fout4(filename4);

	if(!fout)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(!fout2)
	{
		cout << "cannot open" << filename2 << endl;
		exit(EXIT_FAILURE);
	}
	if(!fout3)
	{
		cout << "cannot open" << filename3 << endl;
		exit(EXIT_FAILURE);
	}
	if(!fout4)
	{
		cout << "cannot open" << filename4 << endl;
		exit(EXIT_FAILURE);
	}

	int n=0,n2=0,n3=0,n4=0;
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID)
		{
			if(PART[i].r[face]<r_face+le && PART[i].r[face]>r_face-le)	
			//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n++;
			}

			if(PART[i].r[face2]<6.0*1e-3-shold_R*sin(3.0)-le && PART[i].r[face2]>6.0*1e-3-shold_R*sin(3.0)-3*le)	
			//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout2 << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n2++;
			}

			if(PART[i].r[face]<r_face-shold_R+le && PART[i].r[face]>r_face-shold_R-le)	
			//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout3 << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n3++;
			}
			if(PART[i].r[face3]<r_face+le && PART[i].r[face3]>r_face-le)	
			//if(PART[i].r[A_Y]<0.006+le && PART[i].r[A_Y]>0.006-le)	
			{
				double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
				double y=PART[i].r[A_Y]*1.0E+05;
				double z=PART[i].r[A_Z]*1.0E+05;
				double P=T[i];
				fout4 << P << "\t" << x << "\t" << y << "\t" << z << endl;
				n4++;
			}
		}
	}
	fout.close();
	fout2.close();
	fout3.close();
	fout4.close();
	sprintf_s(filename,"./Temperature/T_XZ%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename2,"./Temperature/T_XY%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename3,"./Temperature/T_XZ_back%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename4,"./Temperature/T_YZ%d.fld",t);//他のファイルと同じ階層に生成するならこちら

	ofstream fout_fld(filename);
	if(!fout_fld)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	fout_fld << "# AVS field file" << endl;
	fout_fld << "ndim=1" << endl;
	fout_fld << "dim1=" << n <<endl;
	fout_fld << "nspace=3" << endl;
	fout_fld << "veclen=1" << endl;
	fout_fld << "data=float" << endl;
	fout_fld << "field=irregular" << endl;
	fout_fld << "label=temperature" << endl << endl;
	//fout_fld << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout_fld << "variable 1 file=T_XZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout_fld << "coord    1 file=T_XZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout_fld << "coord    2 file=T_XZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout_fld << "coord    3 file=T_XZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout_fld.close();


	ofstream fout2_fld(filename2);
	if(!fout2_fld)
	{
		cout << "cannot open" << filename2 << endl;
		exit(EXIT_FAILURE);
	}

	fout2_fld << "# AVS field file" << endl;
	fout2_fld << "ndim=1" << endl;
	fout2_fld << "dim1=" << n2 <<endl;
	fout2_fld << "nspace=3" << endl;
	fout2_fld << "veclen=1" << endl;
	fout2_fld << "data=float" << endl;
	fout2_fld << "field=irregular" << endl;
	fout2_fld << "label=temperature" << endl << endl;
	//fout2_fld << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2_fld << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2_fld << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2_fld << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout2_fld << "variable 1 file=T_XY" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld << "coord    1 file=T_XY" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld << "coord    2 file=T_XY" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld << "coord    3 file=T_XY" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2_fld.close();

	
	ofstream fout3_fld(filename3);
	if(!fout3_fld)
	{
		cout << "cannot open" << filename3 << endl;
		exit(EXIT_FAILURE);
	}

	fout3_fld << "# AVS field file" << endl;
	fout3_fld << "ndim=1" << endl;
	fout3_fld << "dim1=" << n3 <<endl;
	fout3_fld << "nspace=3" << endl;
	fout3_fld << "veclen=1" << endl;
	fout3_fld << "data=float" << endl;
	fout3_fld << "field=irregular" << endl;
	fout3_fld << "label=temperature" << endl << endl;
	//fout_fld << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout_fld << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout3_fld << "variable 1 file=T_XZ_back" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld << "coord    1 file=T_XZ_back" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld << "coord    2 file=T_XZ_back" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld << "coord    3 file=T_XZ_back" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout3_fld.close();

	ofstream fout4_fld(filename4);
	if(!fout4_fld)
	{
		cout << "cannot open" << filename4 << endl;
		exit(EXIT_FAILURE);
	}
	fout4_fld << "# AVS field file" << endl;
	fout4_fld << "ndim=1" << endl;
	fout4_fld << "dim1=" << n4 <<endl;
	fout4_fld << "nspace=3" << endl;
	fout4_fld << "veclen=1" << endl;
	fout4_fld << "data=float" << endl;
	fout4_fld << "field=irregular" << endl;
	fout4_fld << "label=temperature" << endl << endl;
	//fout4_fld << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout4_fld << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout4_fld << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout4_fld << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout4_fld << "variable 1 file=T_YZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld << "coord    1 file=T_YZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld << "coord    2 file=T_YZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld << "coord    3 file=T_YZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4_fld.close();*/

}

///温度場陰的計算関数
void calc_temperature_implicity(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t)
{
	//現在は3Ｄのみ対応
	//解くべき式は	hi(t+1)-hi(t)=Δt*k*(2d/λn0)Σ(hj-hi)w+q*Δt
	//				⇔Δt*k*(2d/λn0)Σ(hj-hi)w-hi(t+1)=-hi(t)-q*Δt
	//				⇔Σ(hj-hi)w-λn0/(2kdΔt)hi(t+1)=-λn0/(2kdΔt)hi(t)-λn0/(2kd)q
	//				⇔kΣ(hj-hi)w-λn0/(2dΔt)hi(t+1)=-λn0/(2dΔt)hi(t)-λn0/(2d)q

	cout<<"温度場計算(陰的)--";

	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;						//ﾗﾌﾟﾗｼｱﾝ用影響半径
	int d=CON->get_dimention();						//次元
	unsigned int timeA=GetTickCount();				//計算開始時刻
	int count=0;
	
	double co=lamda*n0/(2*d*dt);					//計算によく現れる係数
	double co2=lamda*n0/(2*d);

	double V=CON->get_particle_volume();
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
		else cout<<"粒子"<<i<<"の材質IDが不定?"<<endl;
	}
	for(int i=fluid_number;i<particle_number;i++)
	{
		density[i]=CON->get_wall_density();
		Cp[i]=CON->get_wall_Cp();
	}
	double *mass=new double [particle_number];	//各粒子の質量格納
	for(int i=0;i<particle_number;i++) mass[i]=density[i]*V;
		
	
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

   
    double *T=new double [particle_number];//温度
    double *K=new double [particle_number];//各粒子の熱伝導率

   
    ///エンタルピーから各粒子の温度T[i]を求める
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//融解開始点のエンタルピー
		double hs1=hs0+latent_H[i]*mass[i];			//融解終了点のエンタルピー

        if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];			//固体
		else if(hs0<=PART[i].h && PART[i].h<=hs1)
		{
			T[i]=MP[i];	//融点
			//cout<<"溶融中の粒子が存在"<<endl;
		}
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
        T[i]=PART[i].h/mass[i]/Cp[i];//固体
    }
    ///////////////
	
    ///////各粒子の熱伝導率k[i]を求める
    for(int i=0;i<fluid_number;i++)
    {
		double densityCp=density[i]*Cp[i];
		if(PART[i].materialID==1) K[i]=CON->get_k()/densityCp;
		else if(PART[i].materialID==2) K[i]=CON->get_k2()/densityCp;
    }
	for(int i=fluid_number;i<particle_number;i++)
    {
		double densityCp=density[i]*Cp[i];
		K[i]=CON->get_wall_k()/densityCp;//壁の熱伝導率
    }

	//現時点step(k)での温度ラプラシアンを計算　高次近似する際に使用する
	double *T_lap_before1=new double [particle_number];
	for(int i=0;i<particle_number;i++)
	{
		double lam=0;					//各粒子の正確なλ
		double W=0;						//各粒子の正確な粒子数密度
		T_lap_before1[i]=0.0;			//初期化
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k]; 
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);//粒子間の距離
					
			double w=kernel(R,dis);//重み関数
			W+=w;
			lam+=dis*dis*w;
			double dT=T[j]-T[i];
			//T_lap_before1[i]+=dT*w;
			T_lap_before1[i]+=dT*w*K[i]*K[j]/(K[i]+K[j])*2;
		}
		T_lap_before1[i]=T_lap_before1[i]*2*d/(n0*lamda);
	}
   ///////*/

	

	int pn=0;	//未知数
	for(int i=0;i<particle_number;i++) if(PART[i].surface==OFF) pn++; 

	int *ppn = new int[pn];					//行列における第n番目の未知数は粒子番号ppn[n]の粒子に相当
	int *link = new int [particle_number];	//粒子番号iはlink[i]番目の未知数

	double *B   = new double[pn];					//解行列

	count=0;
	//解行列B[]作成 -λn0/(2dΔt)hi(t)
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].surface==OFF)
		{
			ppn[count]=i;
			B[count]=-co*T[i];
			//B[count]=-co*PART[i].h/V;
			
			//B[count]-=PART[i].heat_generation*co2/(density[i]*Cp[i]);	//発熱量[W/m3]
			B[count]-=0.5*(PART[i].heat_generation+PART[i].heat_gene_before1)*co2/(density[i]*Cp[i]);	//発熱量[W/m3] 1ステップ前の発熱量を用いて台形則で近似

			B[count]-=0.5*T_lap_before1[i]*co2;

			link[i]=count;
			count++;
		}
		else link[i]=pn+1;//行列に含まれない粒子にはﾀﾞﾐｰとして(pn+1)を格納
	}

	int number=0;			//係数行列の非ゼロ要素数
	for(int n=0;n<pn;n++)	//pnは断熱条件ならfluid_number,非断熱ならparticle_numberに一致
	{
		int i=ppn[n];		//n番目の未知数の粒子番号
		int num=1;//粒子i周辺の流体粒子数 初期値が1なのは自分自身をｶｳﾝﾄしているから
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k];
			if(link[j]<pn) num++;
		}
		number+=num;
	}///numberがもとまった
	
    double *val = new double [number];
	int *ind = new int [number];//非ゼロ要素の列番号格納配列
	int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
	
	/////////////////////val,ind ,ptrに値を格納
	int index=0;
	for(int n=0;n<pn;n++)
	{   
		
		ptr[n]=index;
		int KK=index;		//matrixの対角成分が格納される場所を記憶
		double AA=0;
	    ind[index]=n;
	    index++;

		int i=ppn[n];

	    for(int k=0;k<PART[i].N2;k++)
	    {
			
	        int j=PART[i].NEI2[k]; 
			
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double w=kernel(R,dis);

			if(link[j]<pn)//粒子jが未知数なら
			{
				val[index]=w*K[i]*K[j]/(K[i]+K[j])*2*0.5;
				AA+=val[index];
				ind[index]=link[j];
				index++;
			}
			else
			{//ノイマン型なら、ここは記述なし？
				//B[n]-=w*PART[j].h/V*K[i];//粒子jが流体でないなら
				//AA+=w*K[i];
			}
	    }
		val[KK]=-AA-co;
	}
	ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない
	////////////////////*/

	//ここで作られる行列は対角成分がすべて負なので、これを正にするために、係数行列と解行列に-1をかける
	for(int n=0;n<pn;n++)
	{
		for(int j=ptr[n];j<ptr[n+1];j++) val[j]*=-1;
		B[n]*=-1;
	}

	double *r=new double[pn];
	double *X=new double[pn];		//行列の答え格納
	double *AP = new double [pn];
	double *P = new double [pn];

	/////////////////////////初期値//////////////////
	if(CON->get_initial_u()==OFF)
	{
		for(int n=0;n<pn;n++) 
		{
			 X[n]=0;
			 r[n]=B[n];
			 P[n]=r[n];
		}
	}
	else if(CON->get_initial_u()==ON) //初期値を与える場合　そんなに効果ない感じ
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//ｎ番目の未知数に該当する粒子番号
			X[n]=T[i];
		}
		for(int n=0;n<pn;n++)
		{
			double AX=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
			r[n]=B[n]-AX;
			P[n]=r[n];
		}
	}
	//////////////////////////////////////////////

	if(CON->get_T_CG()==0)//CG法
	{
		//cout<<CON->get_T_CGep()<<endl;
		CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,CON->get_T_CGep()); //CG法により行列を解く

		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//ｎ番目の未知数に該当する粒子番号
			PART[i].h+=(X[n]-T[i])*V;//エンタルピー更新
	//		laplacian[D][i]=(X[n]-PART[i].u[D])/(dt*vis0);//こっちだとrenewal関数ないで速度更新する形になる
		}
		cout<<"反復回数:"<<count<<"  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	}
	//////////*/

	//cout<<mass[191480]<<" "<<Cp[191480]<<endl;

	if(CON->get_T_CG()==1)//ICCG法
	{
		for(int n=0;n<pn;n++)//ICCG法でとく場合、係数行列をきちんと並び変えしなくてはならない
		{      
			int num=ptr[n+1]-ptr[n];
			for(int j=ptr[n]+1;j<ptr[n+1];j++)
			{
				for(int m=ptr[n];m<j;m++)
				{
					if(ind[j]<ind[m])
					{
						double temp=val[m];
						int tempR=ind[m];
						val[m]=val[j];
						ind[m]=ind[j];
						val[j]=temp;
						ind[j]=tempR;
					}
				}
			}
		}

		iccg(CON,val,ind,ptr,pn,B,number,X,r,P,CON->get_T_CGep(),&count);
		
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//ｎ番目の未知数に該当する粒子番号
			//cout<<i<<endl;
			//PART[i].h+=(X[n]-T[i])*V;//エンタルピー更新
			//PART[i].h=X[n]*V;//エンタルピー更新
			PART[i].h=X[n]*mass[i]*Cp[i];//エンタルピー更新
		}
		cout<<"反復回数:"<<count<<"  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	}
	///////////////////////////*/

	//発熱量を初期化
	for(int i=0;i<particle_number;i++)
	{
		PART[i].heat_gene_before1=PART[i].heat_generation;	//1step前の情報へと格納しなおし
		PART[i].heat_generation=0;							//初期化　発熱量を増やすときは+=の記述にすると、いろんな関数内で発熱が起こったときにも対処できる
	}


	///エンタルピーから各粒子の温度T[i]を求める(温度をプロットするため)
    for(int i=0;i<fluid_number;i++)
    {
		double hs0=mass[i]*Cp[i]*MP[i];		//融解開始点のエンタルピー
		double hs1=hs0+latent_H[i]*mass[i];			//融解終了点のエンタルピー

        if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];			//固体
		else if(hs0<=PART[i].h && PART[i].h<=hs1) T[i]=MP[i];	//融点
		else if(hs1<PART[i].h)											//液体
		{       
			T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];
			if(CON->get_material()==H2O && T[i]>=393) T[i]=393;//現在沸騰は考えないので、もし温度が100度越えすれば100度に戻す(気化熱の実装の方が先決か?）
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
    //////////////*/

	////温度プロット
	double height=CON->get_TplotZ();
    plot_T(CON ,PART,particle_number,T,height);
	if(CON->get_T_AVS()>0)
	{
		//if(t%CON->get_T_AVS()==0 || t==1)output_temperature_avs(CON,PART, t, particle_number, fluid_number,T, height);	//05/08
		if(t%CON->get_T_AVS()==0 || t==1)output_temperature_avs2(CON,PART, t, particle_number, fluid_number,T, height);	//06/03
	}


	if(CON->get_dimention()==3)//XZ平面のTも出力
	{
		
		ofstream fh("T_XZ.dat");
		for(int i=0;i<particle_number;i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fh<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<T[i]<<endl;
		fh.close();
	}

	double Tmax=0;
	for(int i=0;i<particle_number;i++) if(T[i]>Tmax) Tmax=T[i];
	cout<<"最大温度="<<Tmax<<"[K]"<<endl;

	delete [] val;
	delete [] ind;
	delete [] ptr;
    delete [] B;

	delete [] r;
	delete [] X;
	delete [] AP;
	delete [] P;

	delete [] T;
    delete [] K;

	delete [] density;
	delete [] mass;
	delete [] Cp;
	delete [] MP;
	delete [] latent_H;

	delete [] ppn;
	delete [] link;

	delete [] T_lap_before1;
	
}

///温度境界条件設定関数 //境界の温度そのものが実験などで既にもとまっている場合、それを対応する境界に位置する粒子に与える
void set_temperature_boundary(mpsconfig *CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, double *T, double dt,int t)
{
	//cout<<"温度境界条件の読み込み"<<endl;
	////自動で読めるようにすること
	int b_num=16;//測定点の数
	int t_num=469;//時間条件の数
	double time=dt*t;//現在の時刻
	////
	double **bound_T=new double* [t_num];//温度の境界条件
	for(int i=0;i<t_num;i++) bound_T[i]=new double [b_num];
	
	for(int i=0;i<t_num;i++) for(int j=0;j<b_num;j++) bound_T[i][j]=0.0; 
	//データの読み込み
    

	//ファイルからデータを読み込む。
	int x=0;


	ifstream fin("Tbound.prn");//エクセルの元データから、スペース区切りで出力したファイル
	if(!fin) cout<<"cannot open Tbound"<<endl;
	fin.unsetf(ifstream::dec);
	fin.setf(ifstream::skipws);
	for(int i=0;i<t_num;i++)
	{ 
		for(int j=0;j<b_num;j++)
		{ 
			fin>>bound_T[i][j];
		}
	}
	fin.close();	

	for(int i=0;i<t_num;i++) for(int j=0;j<b_num;j++) bound_T[i][j]+=273; //元データは℃なので、Kに修正
	//cout<<"x="<<x<<" boundT "<<bound_T[t_num-1][b_num-1]<<endl;

	if(CON->get_model_number()==26)//IH釜モデル
	{

		
		if(CON->get_dimention()==2)//2次元
		{

			for(int i=0; i<particle_number; i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL) //壁境界部はデータを利用する 座標系で補間した後、時間系で補間する
				{
					double X= abs(PART[i].r[A_X]);
					double Y= PART[i].r[A_Y]+0.04; //model配置時は流体下部と上部のy座標絶対値が一致しているため、わかりやすいように下部の座標が0になるように補正する
					double temp_t1;
					double temp_t2;
					int t1=(int) floor(time);
					int t2=(int) ceil(time);

					int x1=(int) floor(abs(X*100));//測定点は10mmおき
					int x2=(int) ceil(abs(X*100));
					int y1=(int) floor(Y*100);
					int y2=(int) ceil(Y*100);

					if(time<t_num-1)
					{
						if(Y<=0.0)//下壁
						{
							if(abs(X)<=0.08) //下壁のうち、測定値が存在する範囲
							{
								if(t1!=t2)
								{
								
									if(x1!=x2)
									{
										temp_t1= bound_T[t1][x1]+(bound_T[t1][x2]-bound_T[t1][x1])*(X-0.01*x1)/(0.01*(x2-x1)); //一般的な直線の方程式。y-y1=(y2-y1)(x-x1)/(x2-x1)
										temp_t2= bound_T[t2][x1]+(bound_T[t2][x2]-bound_T[t2][x1])*(X-0.01*x1)/(0.01*(x2-x1));
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
									else
									{
										temp_t1= bound_T[t1][x1];
										temp_t2= bound_T[t2][x1];
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
								}
								else
								{
									if(x1!=x2)
									{
										T[i] = bound_T[t1][x1]+(bound_T[t1][x2]-bound_T[t1][x1])*(X-0.01*x1)/(0.01*(x2-x1)); //一般的な直線の方程式。y-y1=(y2-y1)(x-x1)/(x2-x1)
									}
									else
									{
										T[i]= bound_T[t1][x1];
										
									}
								}
							}
							else //下壁、左右方向のoutwall部分。参照する値が無いのでとりあえずX=0.08のときの値を利用
							{
								if(t1!=t2)
								{
									temp_t1= bound_T[t1][8];
									temp_t2= bound_T[t2][8];
									T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
								}
								else
								{
									T[i]= bound_T[t1][8];
								}
							}
						
						}
						else//左右の壁
						{
							if(Y<0.07) //下壁のうち、測定値が存在する範囲
							{
								//y1+=8; //0から8番目の測定値は底面(8が角)なので、8から15までのy座標に応じ参照するように修正する
								//y2+=8;
								if(t1!=t2)
								{
									if(y1!=y2)
									{
										temp_t1= bound_T[t1][y1+8]+(bound_T[t1][y2+8]-bound_T[t1][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1)); //一般的な直線の方程式。y-y1=(y2-y1)(x-x1)/(x2-x1)
										temp_t2= bound_T[t2][y1+8]+(bound_T[t2][y2+8]-bound_T[t2][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1));
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
									else
									{
										temp_t1= bound_T[t1][y1+8];
										temp_t2= bound_T[t2][y1+8];
										T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
									}
								}
								else
								{
									if(y1!=y2)
									{
										T[i]= bound_T[t1][y1+8]+(bound_T[t1][y2+8]-bound_T[t1][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1)); //一般的な直線の方程式。y-y1=(y2-y1)(x-x1)/(x2-x1)
									}
									else
									{
										T[i] = bound_T[t1][y1+8];
									}
								}

							}
							else //左右壁、上方向の測定値が存在しない部分。とりあえずY=0.07のときの値を利用
							{
								if(t1!=t2)
								{
									temp_t1= bound_T[t1][15];
									temp_t2= bound_T[t2][15];
									T[i] = temp_t1+(temp_t2-temp_t1)*(time-t1)/(t2-t1);
								}
								else
								{
									T[i]= bound_T[t1][15];
								}
							}
						
						}
					}
					else//時間条件が測定値より越えている場合、最終時刻のものを利用する
					{
						if(Y<=0.0)//下壁
						{
							if(abs(X)<=0.08) //下壁のうち、測定値が存在する範囲
							{
								if(x1!=x2)
								{
									T[i]= bound_T[t_num-1][x1]+(bound_T[t_num-1][x2]-bound_T[t_num-1][x1])*(X-0.01*x1)/(0.01*(x2-x1)); //一般的な直線の方程式。y-y1=(y2-y1)(x-x1)/(x2-x1)
								}
								else
								{
									T[i]= bound_T[t_num-1][x1];
								}
							}
							else //下壁、左右方向のoutwall部分。参照する値が無いのでとりあえずX=0.08のときの値を利用
							{
								T[i]= bound_T[t_num-1][8];
							}
						}
						else//左右の壁
						{
							if(Y<0.07) //下壁のうち、測定値が存在する範囲
							{
								//y1+=8; //0から8番目の測定値は底面(8が角)なので、8から15までのy座標に応じ参照するように修正する
								//y2+=8;
								if(y1!=y2)
								{
									T[i]= bound_T[t_num-1][y1+8]+(bound_T[t_num-1][y2+8]-bound_T[t_num-1][y1+8])*(Y-0.01*y1)/(0.01*(y2-y1)); //一般的な直線の方程式。y-y1=(y2-y1)(x-x1)/(x2-x1)
								}
								else
								{
									T[i]= bound_T[t_num-1][y1+8];
								}
							}
							else //左右壁、上方向の測定値が存在しない部分。とりあえずY=0.07のときの値を利用
							{
								T[i]= bound_T[t_num-1][15];
							}
						}
					}
				}
			}
		
		}

		if(CON->get_dimention()==3)//3次元
		{
		}



	}

	for(int i=0;i<t_num;i++) delete [] bound_T[i];
    delete [] bound_T;
}

///発熱境界条件設定関数 //境界の温度などが実験などで既にもとまっている場合、それを対応する境界に位置する粒子に与える
void set_Q_boundary(mpsconfig *CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, double dt,int t)
{
	//cout<<"発熱条件の読み込み"<<endl;
	////自動で読めるようにすること
	int b_num=101;//測定点の数
	int t_num=1;//時間条件の数
	double time=dt*t;//現在の時刻
	////
	double **bound_T=new double* [t_num];//温度の境界条件
	for(int i=0;i<t_num;i++) bound_T[i]=new double [b_num];
	
	for(int i=0;i<t_num;i++) for(int j=0;j<b_num;j++) bound_T[i][j]=0.0; 
	//データの読み込み
    

	//ファイルからデータを読み込む。
	int x=0;


	ifstream fin("q0_f.txt");//エクセルの元データから、スペース区切りで出力したファイル
	if(!fin) cout<<"cannot open q0_f"<<endl;
	fin.unsetf(ifstream::dec);
	fin.setf(ifstream::skipws);
	for(int i=0;i<t_num;i++)
	{ 
		for(int j=0;j<b_num;j++)
		{ 
			//fin>>x;
			fin>>bound_T[i][j];
		}
	}
	fin.close();	

	//cout<<"x="<<x<<" boundT "<<bound_T[t_num-1][b_num-1]<<endl;

	if(CON->get_model_number()==26)//IH釜モデル
	{
		if(CON->get_dimention()==2)//2次元
		{

			for(int i=0; i<particle_number; i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL) //壁境界部はデータを利用する 座標系で補間した後、時間系で補間する
				{
					double X= abs(PART[i].r[A_X]);
					double Y= PART[i].r[A_Y]; //model配置時は流体下部と上部のy座標絶対値が一致しているため、わかりやすいように下部の座標が0になるように補正する
					double temp_t1;
					double temp_t2;
					//int t1=(int) floor(time);
					//int t2=(int) ceil(time);
					int t1=0;
					int t2=0;

					int x1=(int) floor(abs(X/0.0008));//測定点は0.8mmおき
					int x2=(int) ceil (abs(X/0.0008));
					int y1=(int) floor(Y/0.0008);
					int y2=(int) ceil(Y/0.0008);

					//if(t1==t2) cout<<"?"<<endl;

					//if(time<t_num-1)
					{
						//if(Y<=-0.05-0.0001)//下壁 INWALLとOUTWALL全部
						if(Y<=-0.058-0.0001)//下壁 OUTWALLの最外部だけになる
						{
							if(abs(X)<=0.08) //下壁のうち、測定値が存在する範囲
							{
								if(x1!=x2)
								{
									PART[i].heat_generation += bound_T[t1][x1]+(bound_T[t1][x2]-bound_T[t1][x1])*(X-0.0008*x1)/(0.0008*(x2-x1)); //一般的な直線の方程式。y-y1=(y2-y1)(x-x1)/(x2-x1)
								}
								else
								{
									PART[i].heat_generation += bound_T[t1][x1];
								}
							}
							else //下壁、左右方向のoutwall部分。参照する値が無いのでとりあえずX=0.08のときの値を利用
							{
								PART[i].heat_generation += bound_T[t1][100];
							}
						
						}
					}
					
				}
				//実際の発熱部分の厚みを考慮し、発熱量密度を変更する
				PART[i].heat_generation/=CON->get_particle_volume()/(0.00055*CON->get_distancebp());
			}
		
		}

		if(CON->get_dimention()==3)//3次元
		{
		}



	}

	for(int i=0;i<t_num;i++) delete [] bound_T[i];
    delete [] bound_T;
}