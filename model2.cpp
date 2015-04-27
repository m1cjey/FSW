#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納

#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
#include<vector>
#define FULL 1
#define HALF 2
#define HALFD 3
#define HALFD_shell 4

void writedata2(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM);
double get_volume(mpsconfig *CON);

//直線を分割するさいの最適な分割数と分割距離の算出関数 他所で使うためfunction.hに移動
void calc_N_and_L(double dis,double le,int *N,double *L);
//円周分割数計算関数
int calc_division_N_circle(double dis,double le);
//半径Rの円の外周
void set_circle_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R);
//半径Rの円内部
void set_circle_in(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &type1,int *number,double le,double R,int edge_startID,int edge_lastID);
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int edge_startID,int edge_lastID);
//半径Rの球作成関数
//長方形作成関数
void set_rectangular(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height);
void set_sphere(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag);
void set_sphere2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag,int *suf_num);
//長方形の辺作成関数
void set_rectangular_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V);
//長方形内部作成関数
void set_rectangular_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V,int edge_startID,int edge_lastID);
//円柱表面作成関数
void set_cylinder_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int circle_start_id,int circle_end_id,int top_flag);
//円錐表面作成関数
void set_circular_cone_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R2,double height,int circle_start_id,int circle_end_id,int top_flag);
//円柱内部設置関数
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int flag);
//ドーナツ作成
void set_doughnut2D(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_big,double R_smal,int edge_startID,int edge_lastID);
//FSWプローブ内部粒子セット関数
void set_hat_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID);
void set_hat_in_2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_mid,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID);
//箱作成関数
void set_box(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &surface,int *number,double le,double Width,double Height,double Depth);
//BOX内作成関数
void set_box_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height,double Depth,int BO_startID,int BO_lastID);
//るつぼ壁内作成関数（るつぼの肉厚部）
void set_crucible_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R_out,double height,double height_out,int flag,int fluid_number);

//分子動力学関数
void MD_2D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN);
void MD_3D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN,double r,double region[3][2]);

//物質合成関数
void make_fusion3D(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<double> &X2,vector<double> &Y2,vector<double> &Z2,vector<int> &surface2,int *number,double le);

void set_initial_placement_using_MD(mpsconfig *CON,int *particle_number)
{
	cout<<"初期粒子配置をMDにより最適化中--";
	unsigned int timeA=GetTickCount();
	int number=0;	//粒子数　最後にはparticle_numberに格納
	int model=CON->get_model_number();
	int Dim=CON->get_dimention();		//解析次元

	double le=CON->get_distancebp();	//初期粒子間距離
	double A=sqrt(3.0)*0.5;				//よく使う係数
	double B=sqrt(2.0/3);				//よく使う係数 高さ方向間隔
	

	vector<double> X;
	vector<double> Y;
	vector<double> Z;
	vector<int> type1;
	vector<int> surface;

	vector<double> X2;					//使い回し用
	vector<double> Y2;
	vector<double> Z2;
	vector<int> surface2;			//ONなら表面　OFFなら内部

	vector<double> X3;					//使い回し用 弾丸２用
	vector<double> Y3;
	vector<double> Z3;
	vector<int> surface3;			//ONなら表面　OFFなら内部

	vector<double> X4;					//使い回し用　弾丸２用
	vector<double> Y4;
	vector<double> Z4;

	vector<double> X5;					//使い回し用　弾丸1用
	vector<double> Y5;
	vector<double> Z5;
	vector<int> surface5;			//ONなら表面　OFFなら内部
	
	vector<double> X6;					//使い回し用　弾丸1用
	vector<double> Y6;
	vector<double> Z6;
	
	ofstream fq("initial_input.dat");
	if(model==2)//////モデル2　正方水滴
	{
		if(Dim==2)
		{
			double origin[3]={0,0,0};//四角形の左下の点の座標
			double Width=CON->get_fluidwidth()*le;				//水平方向長さ
			double Height=CON->get_fluidwidth()*le;				//垂直方向長さ
			set_rectangular(X,Y,Z,&number,le,Width,Height);		//座標の原点は長方形の中心。よってあとで移動すること。
			
			//初期情報書き込み
			for(int i=0;i<number;i++) writedata2(fq,i,0.5*Width+origin[A_X]+X[i],Height+origin[A_Y]+Y[i],Z[i],FRFLUID,1,0,0,0,0,0,0,0,0);

		}
		else if(Dim==3)//立方体の箱
		{
			double Width=CON->get_fluidwidth()*le;		
			double Height=CON->get_fluidwidth()*le;	
			double Depth=CON->get_fluidwidth()*le;
			set_box(X,Y,Z,surface,&number,le,Width,Height,Depth);//最後3つの引数は横、高さ、奥行き。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ
			for(int i=0;i<number;i++) writedata2(fq,  i,X[i],Y[i],Z[i]+0.19,FLUID,1, OFF,0, 0,0, 0, 0, 0, 0);
		}
	}
	else if(model==3)//////モデル３　円形水滴
	{
		double R=CON->get_fluidwidth()*le*0.5;			//作成する円の半径
		double Zg=CON->get_height();					//球の中心高さ
		if(Dim==2)
		{
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
		}
		else if(Dim==3)
		{
			//円作成
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている

			//球作成
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);
		}
		//初期情報書き込み
		//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,Y[i]*40,0,0,0,0);
		for(int i=0;i<number;i++) writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,FLUID,1, OFF,0, 0,0, 0, 0, 0, ON);
		//for(int i=0;i<number;i++) writedata2(fq,  i,X[i],Y[i],Z[i],FLUID,1, OFF,0, 0,Y[i]*40, 0, 0, 0, ON);
	}
	else if(model==14)//14 静電霧化
	{
		double R=CON->get_fluidwidth()*le;
		double height=6*le*A;
		if(Dim==3)
		{
			//円作成
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
			int circle_end_id=number;	//円の粒子idを記憶
			////////

			for(int i=0;i<number;i++)//X2などに円粒子の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			//半球作成
			int flag=HALF;
			set_sphere(X,Y,Z,&number,le,R,flag);//半球作成
			/////////////////////

			for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,0,0,0,0,0,0,1);//ここまでの粒子はすべて流体

			int beforeN=number;
			int top_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,le,R,height,circle_start_id,circle_end_id,top_flag);//円柱表面作成

			for(int i=beforeN;i<number;i++)//X2などに円柱表面の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			for(int i=beforeN;i<number;i++) Z[i]*=-1.0;//Z座標を反転
			for(int i=beforeN;i<number;i++)
			{
				if(Z[i]>-1.5*le) writedata2(fq,i,X[i],Y[i],Z[i],INWALL,1,0,0,0,0,0,0,0,1);//粒子はINWALL
			}
			int num1=beforeN;//ここまでに発生しているOUTWALLはまだ書き込めない。そこで数を記憶
			int num2=number;

			

			beforeN=number;			//beforeNを更新
			int number2=(int) X2.size();//X2などが格納している粒子数
			int n=number2;		//値を保存
			set_cylinder_in(X2,Y2,Z2,&number2,le,R,height,1);//円柱内部 まずはX2などに座標を格納させる
			
			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//円柱内部の座標をX,Y,Zにｺﾋﾟｰ。だたしZは反転する。
				Y.push_back(Y2[i]);
				Z.push_back(-1.0*Z2[i]);
				number++;
			}


			int count=beforeN;
			for(int i=beforeN;i<number;i++) 
			{
				if(Z[i]>-1.5*le) 
				{
					writedata2(fq,count,X[i],Y[i],Z[i],INWALL,1,0,0,0,0,0,0,0,1);//粒子はINWALL
					count++;
				}
			}
			for(int i=beforeN;i<number;i++) 
			{
				if(Z[i]<=-1.5*le) 
				{
					writedata2(fq,count,X[i],Y[i],Z[i],OUTWALL,1,0,0,0,0,0,0,0,0);
					count++;
				}
			}
			for(int i=num1;i<num2;i++)
			{
				if(Z[i]<=-1.5*le) writedata2(fq,i,X[i],Y[i],Z[i],OUTWALL,1,0,0,0,0,0,0,0,1);//粒子はOUTWALL
			}
			
		}
	}
	else if(model==16)//////モデル16　電磁浮遊
	{
		double R=CON->get_fluidwidth()*le*0.5;			//作成する円の半径
		double Zg=CON->get_height();					//球の中心高さ
		if(Dim==2)
		{
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
		}
		else if(Dim==3)
		{
			//円作成
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている

			//球作成
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);
		}
		//初期情報書き込み
		for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zg,FLUID,1,0,0,0,0,0,0,0,1);
	}
	else if(model==20)//////るつぼと球形溶融金属
	{
		double R=CON->get_fluidwidth()*0.001;//3;			//作成する円の半径
		double Zg=CON->get_height();					//球の中心高さ
		double height=0.1;//10;							//るつぼ円筒部の高さ
		double C_R=0.03;//3;			//るつぼの半径
		
		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//温度場関係の変数定義
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		//double T=CON->get_roomT();
		double T=CON->get_initialT();//初期温度
		double V=CON->get_particle_volume();
		double mass=density*V;	//粒子質量

		double h;//エンタルピー
		if(T>=CON->get_MP())
		{
			h=T*mass*Cp+latent_H*mass;
		}
		else h=T*mass*Cp;
									
		double wallmass=V*CON->get_wall_density();	//壁粒子の質量
		double wall_h=T*wallmass*CON->get_wall_Cp();//壁粒子のエンタルピー

		height=height/2;
		//C_R-=le;	//幅1mmのコーティング層を想定				
		if(Dim==3)
		{
		//////////溶融金属
			//球形
			//円作成
			//R=C_R-le;		//るつぼに接した球
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている

			//cout<<"円作成"<<endl;

			//球作成
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);

			//cout<<"球作成"<<endl;

			int fluid_number=number;
			//cout<<"総粒子数="<<fluid_number<<endl;

			if(CON->get_mesher()==0) for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zg+0.13125/*+0.18125-0.002*/,FLUID,materialID,1,0,0,0,0,0,h,1);
			if(CON->get_mesher()==1) for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zg,FLUID,materialID,1,0,0,0,0,0,h,1);

			//cout<<"出力完了="<<fluid_number<<endl;

		//////////るつぼ

			

			////内壁
			//基準円
			int circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			int circle_end_id=number;//円周に配置された粒子数
			//るつぼ側面作成
			int top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//円筒
			//半球殻作成
			int flag2=HALFD_shell;//下半分の球殻
			set_sphere(X,Y,Z,&number,le,C_R,flag2);

			int in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
		
		//外壁
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//基準円
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			circle_end_id=number;//円周に配置された粒子数
			//set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//円内部    vector配列は参照渡ししている
			//るつぼ側面作成
			top_flag=HALF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//円筒
			//半球殻作成
			flag2=HALFD_shell;//下半分の球殻
			set_sphere(X,Y,Z,&number,le,C_R_out,flag2);
			//for(int i=in_number;i<number;i++) Z[i]-=4*le+C_R;

			set_crucible_in(X,Y,Z,&number,le,C_R,C_R_out,height,height_out,1,fluid_number);//るつぼ内部（肉厚部）

			//OUTWALL出力
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		///*////*/
		

		

	      ////////出力
			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125+R,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			
			
		}

	}

	else if(model==26)//////IH釜モデル //円筒状の簡易モデル　本当は厚みが薄いのでoutwallの量に応じ解像度を考慮する必要があるが、測定した温度条件をそのまま境界条件として利用するため気にせずに作成
	{
		double R=CON->get_fluidwidth()*0.001;//3;			//
		double Zg=CON->get_height();					//
		double C_h=0.08;//10;							//釜側面部の高さ
		double C_hout=C_h+4*le;//10;							//釜側面部の高さ
		double C_R=0.07;//3;			//釜底部の半径
		double C_Rout=C_R+4*le;

		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//温度場関係の変数定義
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		//double T=CON->get_roomT();
		double T=CON->get_initialT();//初期温度
		double V=CON->get_particle_volume();
		double mass=density*V;	//粒子質量

		double h;//エンタルピー
		if(T>=CON->get_MP())
		{
			h=T*mass*Cp+latent_H*mass;
		}
		else h=T*mass*Cp;
									
		double wallmass=V*CON->get_wall_density();	//壁粒子の質量
		double wall_h=T*wallmass*CON->get_wall_Cp();//壁粒子のエンタルピー
			

		if(Dim==2)
		{
		}
		if(Dim==3)
		{
		
		//////////釜と流体を同時に作成する

			int fluid_number=0;

			//円作成
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_Rout,0,number);//円内部    vector配列は参照渡ししている
			int circle_end_id=number;	//円の粒子idを記憶
			////////

			for(int i=0;i<number;i++)//X2などに円粒子の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			/////////////////////

			int beforeN=number;
			int top_flag=OFF;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,le,C_Rout,C_hout,circle_start_id,circle_end_id,top_flag);//円柱表面作成

			for(int i=beforeN;i<number;i++)//X2などに円柱表面の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////
			
			beforeN=number;			//beforeNを更新
			int number2=(int) X2.size();//X2などが格納している粒子数
			int n=number2;		//値を保存
			set_cylinder_in(X2,Y2,Z2,&number2,le,C_Rout,C_hout,1);//円柱内部 まずはX2などに座標を格納させる

			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//円柱内部の座標をX,Y,Zにｺﾋﾟｰ。
				Y.push_back(Y2[i]);
				Z.push_back(Z2[i]);
				number++;
			}

			for(int i=0;i<number;i++) Z[i]-=C_hout-C_h;
			
			///////////////////粒子の材質決定

			for(int i=0;i<number;i++) type1.push_back(FLUID);//まずはとりあえずFLUIDを代入

			for(int i=0;i<number;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>C_R-0.4*le)//釜の側面
				{
					if(r<C_R+0.4*le && Z[i]>0.2*le) type1[i]=INWALL;
					else type1[i]=OUTWALL;	
				}
				else//釜の内部および底面
				{
					if(Z[i]<0.2*le && Z[i]>-0.2*le)
					{
						type1[i]=INWALL;
					}
					else if(Z[i]<=-0.2*le) type1[i]=OUTWALL;
				}
			}
			//////////////////////*/


			//粒子座標出力
			int count=0;
			int count_outf=0;//fluidのうち、Z座標で除外される数
			for(int i=0;i<number;i++)
			{
				//if(type1[i]==FLUID)
				if(type1[i]==FLUID)
				{
					if(Z[i]<=0.05)
					{
						count++;
						writedata2(fq,i,X[i],Y[i],Z[i],FLUID,materialID,0,0,0,0,0,0,h,1);
					}
					else count_outf++;
				}
			}

			for(int i=0;i<number;i++)
			{
				if(type1[i]==INWALL)
				{
					count++;
					writedata2(fq,i,X[i],Y[i],Z[i],INWALL,materialID,0,0,0,0,0,0,wall_h,1);
				}
			}
					for(int i=0;i<number;i++)
			{
				if(type1[i]==OUTWALL)
				{
					count++;
					writedata2(fq,i,X[i],Y[i],Z[i],OUTWALL,materialID,0,0,0,0,0,0,wall_h,1);
				}
			}

			number-=count_outf;

	      ////////出力
			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125+R,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			
			
		}

	}
	
	else if(model==21)//るつぼの中に弾丸型の流体：
	{
		double R=CON->get_fluidwidth()*0.001;//3; //流体の半径
		//double R=0.03-le/sqrt(2.0);//3; //流体の半径
		double fluid_h=CON->get_fluid_h()*0.001;//3; //流体円柱部の高さ
		double height=0.1;//10;							//るつぼ円筒部の高さ
		double C_R=0.03;//3;			//るつぼの半径
		//double Zf=0.13125-(0.03-R)+le/sqrt(2.0);//流体のZ方向の位置修正量
		double Zf=0.13125;

		height/=2.0;//上のほうの壁をカット

		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//温度場関係の変数定義
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		double T=CON->get_roomT();
		double V=CON->get_particle_volume();
		double mass=density*V;									//粒子質量
		double h=CON->get_MP()*mass*Cp+latent_H*mass;//エンタルピー
		//double h;								//エンタルピー
		double wallmass=V*CON->get_wall_density();	//壁粒子の質量
		double wall_h=T*wallmass*CON->get_wall_Cp();//壁粒子のエンタルピー
		

		if(Dim==3)
		{
			//流体

			//円作成
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
			int circle_end_id=number;	//円の粒子idを記憶
			////////

			for(int i=0;i<number;i++)//X2などに円粒子の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////


			//半球作成
			int beforeN=number;
			int flag=HALF;
			set_sphere(X,Y,Z,&number,le,R,flag);//半球作成

			for(int i=beforeN;i<number;i++) Z[i]*=-1.0;//Z座標を反転
			/////////////////////

			beforeN=number;
			int top_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,le,R,fluid_h,circle_start_id,circle_end_id,top_flag);//円柱表面作成

			for(int i=beforeN;i<number;i++)//X2などに円柱表面の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////
			
			beforeN=number;			//beforeNを更新
			int number2=(int) X2.size();//X2などが格納している粒子数
			int n=number2;		//値を保存
			set_cylinder_in(X2,Y2,Z2,&number2,le,R,fluid_h,1);//円柱内部 まずはX2などに座標を格納させる

			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//円柱内部の座標をX,Y,Zにｺﾋﾟｰ。
				Y.push_back(Y2[i]);
				Z.push_back(Z2[i]);
				number++;
			}
			
			for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zf,FLUID,1,0,0,0,0,0,0,h,1);//ここまでの粒子はすべて流体
			//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125-le-le/sqrt(2.0),FLUID,1,0,0,0,0,0,0,0,1);//ここまでの粒子はすべて流体

			int fluid_number=number;
			//X.clear(); Y.clear(); Z.clear(); X2.clear(); Y2.clear(); Z2.clear();
			cout<<"流体粒子出力"<<endl;

		//るつぼ
			//number=0;
		////内壁
			//基準円
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			circle_end_id=number;//円周に配置された粒子数
			//るつぼ側面作成
			top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//円筒
			//半球殻作成
			int flag2=HALFD_shell;//下半分の球殻
			set_sphere(X,Y,Z,&number,le,C_R,flag2);

			int in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
		
		//外壁
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//基準円
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			circle_end_id=number;//円周に配置された粒子数
			//set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//円内部    vector配列は参照渡ししている
			//るつぼ側面作成
			top_flag=HALF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//円筒
			//半球殻作成
			flag2=HALFD_shell;//下半分の球殻
			set_sphere(X,Y,Z,&number,le,C_R_out,flag2);
			//for(int i=in_number;i<number;i++) Z[i]-=4*le+C_R;

			int mostout_number=number;

			set_crucible_in(X,Y,Z,&number,le,C_R,C_R_out,height,height_out,1,fluid_number);//るつぼ内部（肉厚部）

			//OUTWALL出力
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		}
	}

	else if(model==23)//壁境界比較モデル
	{
		double R=0.1;//3; //流体の半径
		//double R=0.03-le/sqrt(2.0);//3; //流体の半径
		double fluid_h=0.05;//3; //流体円柱部の高さ
		double height=0.1;//10;							//るつぼ円筒部の高さ
		double C_R=0.1+le;//3;			//るつぼの半径
		//double Zf=0.13125-(0.03-R)+le/sqrt(2.0);//流体のZ方向の位置修正量
		double Zf=0;//0.13125;

		

		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//温度場関係の変数定義
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		double T=CON->get_roomT();
		double V=CON->get_particle_volume();
		double mass=density*V;									//粒子質量
		double h=CON->get_MP()*mass*Cp+latent_H*mass;//エンタルピー
		//double h;								//エンタルピー
		double wallmass=V*CON->get_wall_density();	//壁粒子の質量
		double wall_h=T*wallmass*CON->get_wall_Cp();//壁粒子のエンタルピー
		

		if(Dim==3)
		{
			//流体

			//円作成
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
			int circle_end_id=number;	//円の粒子idを記憶
			////////

			for(int i=0;i<number;i++)//X2などに円粒子の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			int beforeN=number;
			int top_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,le,R,fluid_h,circle_start_id,circle_end_id,top_flag);//円柱表面作成

			for(int i=beforeN;i<number;i++)//X2などに円柱表面の情報をｺﾋﾟｰ
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////
			
			beforeN=number;			//beforeNを更新
			int number2=(int) X2.size();//X2などが格納している粒子数
			int n=number2;		//値を保存
			set_cylinder_in(X2,Y2,Z2,&number2,le,R,fluid_h,1);//円柱内部 まずはX2などに座標を格納させる

			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//円柱内部の座標をX,Y,Zにｺﾋﾟｰ。
				Y.push_back(Y2[i]);
				Z.push_back(Z2[i]);
				number++;
			}

			int in_number=number;
			
			for(int i=0;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i],FLUID,1,0,0,0,0,0,0,h,1);
			}
		
			int fluid_number=number;
			cout<<"流体粒子出力"<<endl;

			////内壁
			//基準円
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			circle_end_id=number;//円周に配置された粒子数
			//るつぼ側面作成
			top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//円筒

			in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]-le,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
			cout<<"内壁"<<endl;
		//外壁
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//基準円
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			circle_end_id=number;//円周に配置された粒子数
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//円内部    vector配列は参照渡ししている
			//るつぼ側面作成
			top_flag=HALF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//円筒

			int mostout_number=number;

			set_cylinder_in(X,Y,Z,&number,le,R,fluid_h,1);//円柱内部 まずはX2などに座標を格納させる

			cout<<"外壁肉厚"<<endl;


			//OUTWALL出力
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]-le,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		}
	}
	else if(model==25)//////渦電流解析簡易モデル
	{
		double R=CON->get_fluidwidth()*0.001;//3;			//作成する円の半径
		double Zg=CON->get_height();					//球の中心高さ
		double height=0.1;//10;							//るつぼ円筒部の高さ
		double C_R=0.03;//3;			//るつぼの半径
		
		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//温度場関係の変数定義
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		//double T=CON->get_roomT();
		double T=CON->get_initialT();//初期温度
		double V=CON->get_particle_volume();
		double mass=density*V;	//粒子質量

		double h;//エンタルピー
		if(T>=CON->get_MP())
		{
			h=T*mass*Cp+latent_H*mass;
		}
		else h=T*mass*Cp;
									
		double wallmass=V*CON->get_wall_density();	//壁粒子の質量
		double wall_h=T*wallmass*CON->get_wall_Cp();//壁粒子のエンタルピー

		height=height/2;
		//C_R-=le;	//幅1mmのコーティング層を想定				
		if(Dim==3)
		{
		//////////溶融金属
			//球形
			//円作成
			//R=C_R-le;		//るつぼに接した球
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている

			//cout<<"円作成"<<endl;

			//球作成
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);

			//cout<<"球作成"<<endl;

			int fluid_number=number;
			//cout<<"総粒子数="<<fluid_number<<endl;

			
			for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FLUID,materialID,1,0,0,0,0,0,h,1);

			//cout<<"出力完了="<<fluid_number<<endl;

		/*/////////るつぼ

			

			////内壁
			//基準円
			int circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			int circle_end_id=number;//円周に配置された粒子数
			//るつぼ側面作成
			int top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//円筒
			//半球殻作成
			int flag2=HALFD_shell;//下半分の球殻
			set_sphere(X,Y,Z,&number,le,C_R,flag2);

			int in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
		
		//外壁
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//基準円
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//円外周   vector配列は参照渡ししている この円を基準に円筒、球殻をつくる
			circle_end_id=number;//円周に配置された粒子数
			//set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//円内部    vector配列は参照渡ししている
			//るつぼ側面作成
			top_flag=HALF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//円筒
			//半球殻作成
			flag2=HALFD_shell;//下半分の球殻
			set_sphere(X,Y,Z,&number,le,C_R_out,flag2);
			//for(int i=in_number;i<number;i++) Z[i]-=4*le+C_R;

			set_crucible_in(X,Y,Z,&number,le,C_R,C_R_out,height,height_out,1,fluid_number);//るつぼ内部（肉厚部）

			//OUTWALL出力
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		///*////*/	
		}
	}
	else if(model==19)//19 FSW
	{
		double probe_R=2.5*1e-3;	//プローブ半径  //CON->get_fluidwidth()*le;
		double probe_R_2=1.0*1e-3;	//プローブ半径(円錐形状の下部)  tool_type=1のときに利用;
		double shold_R=6*1e-3;	//ショルダー半径
		double height=4*1e-3;	//プローブ高さ //6*le*A;
		if(CON->get_tool_type()==2) height=3*1e-3+B*le;//表裏ツールの場合、ツールの長さ*2が流体領域の厚みと同じにならないといけない

		double shold_height=10*B*le;
		double rpm=500;//ツール回転速度
		double rps=rpm/60;
		double w=rps*2*PI;		//角速度
		double U=CON->get_move_speed();//プローブの移動速度[m/sec]		
		double pich=0.7e-3;		//プローブのねじのピッチ0.7[mm]
		double T=CON->get_roomT();	//初期温度
		double vol=get_volume(CON);				//粒子の体積
		double h;								//エンタルピー
		double wallmass=vol*CON->get_wall_density();	//壁粒子の質量
		double wall_h=T*wallmass*CON->get_wall_Cp();//壁粒子のエンタルピー
		double val=0;
		int materialID=1;
		int count2;
		int beforeN;

		//calclation type
		int plunge=0;
		int traverse=1;
		int calc_type=CON->get_process_type();		//0:plunge 1:traverse
		if(Dim==3)
		{
			///プローブ底面作成
			if(CON->get_tool_type()==0 ||CON->get_tool_type()==2)
			{
				set_circle_edge(X,Y,Z,&number,le,probe_R);//円外周
				set_circle_in_using_6_pieces(X,Y,Z,&number,le,probe_R,0,number);//円内部    vector配列は参照渡ししている
			}
			if(CON->get_tool_type()==1)//円錐ツール
			{
				set_circle_edge(X,Y,Z,&number,le,probe_R_2);//円外周
				set_circle_in_using_6_pieces(X,Y,Z,&number,le,probe_R_2,0,number);//円内部    vector配列は参照渡ししている
			}
			int circle_start_id=0;
			int circle_end_id=number;
			int circle_num=number;		//円周に配置された粒子数
			////////
			int top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
			//プローブ側面作成
			if(CON->get_tool_type()==0 ||CON->get_tool_type()==2)
			{
				top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
				set_cylinder_face(X,Y,Z,&number,le,probe_R,height,circle_start_id,circle_end_id,top_flag);//円柱表面作成
			}
			if(CON->get_tool_type()==1)
			{
				top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
				set_circular_cone_face(X,Y,Z,&number,le,probe_R,probe_R_2,height,circle_start_id,circle_end_id,top_flag);//円柱表面作成
			}

			//ショルダー底面作成
			beforeN=number;
			set_doughnut2D(X,Y,Z,&number,le,shold_R,probe_R,number-circle_num,number);//最後4つの引数は大きい半径、小さい半径、すでに作成してある円周のstartID,endID
			for(int i=beforeN;i<number;i++) Z[i]=height;//上の関数で作成した外周のZ座標はゼロのままなので、ここで修正
			int shld_botom_startID=beforeN;
			int shld_botom_endID=number;

			//ショルダー側面作成
			beforeN=number;
			top_flag=OFF;		//円柱の上面は作成しなくていいからフラグをOFF 
			set_cylinder_face(X,Y,Z,&number,le,shold_R,shold_height,shld_botom_startID,shld_botom_endID,top_flag);//円柱表面作成 
			for(int i=beforeN;i<number;i++) Z[i]+=height;

			//ショルダー上面
			count2=0;//ショルダ上面における境界粒子は、上のset_cylinder_face()で最後に付け加えられた粒子群である。まずはその数を数える
			double H=height+shold_height;//ショルダーのてっぺんの高さ
			for(int i=beforeN;i<number;i++) if(Z[i]<H+0.1*le && Z[i]>H-0.1*le) count2++;
			beforeN=number;
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,shold_R,number-count2,number);//円内部
			for(int i=beforeN;i<number;i++) Z[i]=H;

			for(int i=0;i<number;i++) type1.push_back(INWALL);//ここまでの粒子はすべてINWALL

			//内部作成
			beforeN=number;
			if(CON->get_tool_type()==0 ||CON->get_tool_type()==2) set_hat_in(X,Y,Z,&number,le,probe_R,shold_R,height,shold_height,0,number);
			if(CON->get_tool_type()==1) set_hat_in_2(X,Y,Z,&number,le,probe_R_2,probe_R,shold_R,height,shold_height,0,number);

			for(int i=beforeN;i<number;i++) type1.push_back(OUTWALL);//ここまでの粒子はすべてOUTWALL
			
			beforeN=number;

			if(CON->get_tool_type()==2)///
			{
				for(int i=0;i<beforeN;i++)
				{
					if(Z[i]!=0)//Z=0面が裏表コピーの境界面
					{
						X.push_back(X[i]);
						Y.push_back(Y[i]);
						Z.push_back(-Z[i]);
						type1.push_back(type1[i]);
						number++;
					}
					double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
					if(Z[i]==0 && r<probe_R-0.1*le)//境界面の外周部以外
					{
						type1[i]=OUTWALL;
					}
				}
			}
			
			double max_dis_z=0.0;
			if(CON->get_tool_angle()>0)//ツールをx軸まわりに回転させる。進行方向が+y
			{
				double theta=PI/180*CON->get_tool_angle();//回転する角度
				
				for(int i=0;i<number;i++)
				{
					X3.push_back(X[i]);//ツール角度を変えた際の速度条件の設定を簡単にするため、回転前の座標を記憶
					Y3.push_back(Y[i]);
					Z3.push_back(Z[i]+2*1e-3);
					//////////////
					double x=X[i];
					double y=cos(theta)*Y[i]-sin(theta)*Z[i];
					double z=sin(theta)*Y[i]+cos(theta)*Z[i];
					double dis_z=z-Z[i];
					if(max_dis_z<dis_z) max_dis_z=dis_z;
					X[i]=x;
					Y[i]=y;
					Z[i]=z;

				}		
			}

			if(CON->get_tool_type()!=2) for(int i=0;i<number;i++) Z[i]+=4*1e-3;//重心を移動

			if(CON->get_tool_angle()>0) 
			{
				double theta=PI/180*CON->get_tool_angle();//回転する角度
				for(int i=0;i<number;i++)
				{
					//Z[i]-=max_dis_z;//重心を移動
					Z[i]+=0.006*(1-cos(theta));//重心を移動
				}
			}
			
			//////////////////////////////ツール作成終了*/

			int tool_number=number;
			
			//挿入時の解析は下の行をONにし、ツールの位置をあげる
			if(calc_type==plunge) for(int i=0;i<number;i++) Z[i]+=4*1e-3;//重心を移動

			/////////////////////////////流体作成
			int fluid_number=0;
			double width0=18*1e-3;//25*1e-3;//18*1e-3;//流体が占める幅
			double Width=width0+6*le;		//横18mm+壁粒子を左右に4粒子分
			double Height=6*1e-3+4*le*B;	//高さ6mm+壁粒子を下に4粒子分
			if(CON->get_tool_type()==2) Height=6*1e-3;	//表裏ツールの場合、下側の壁を取り払っているので壁粒子の設定の関係で下側に追加する必要がない

			double Depth=18*1e-3+6*le*A;//27*1e-3+6*le*A;	//奥行き27mm+壁粒子を左右に4粒子分
			double depth0=9e-3;				//ツール中心と、手前の壁との距離
			int BOX_startID=0;				//set_box()開始前の粒子数
			set_box(X2,Y2,Z2,surface2,&fluid_number,le,Width,Height,Depth);//最後3つの引数は横、高さ、奥行き。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ
			int BOX_lastID=fluid_number;	//set_box()終了直後の粒子数

			//double Z_mod=le*B+(1-sqrt(2.0)/2)*le;//ツールが流体へめり込むのを防ぐための調整量　この値だけツール以外の物体が下がる
			double Z_mod=le+(1-sqrt(2.0)/2)*le;
			
			for(int i=BOX_startID;i<BOX_lastID;i++)
			{
				X2[i]-=width0*0.5+3*le;		//重心移動
				Y2[i]-=depth0+3*le*A;
				Z2[i]-=4*le*B+Z_mod;
				if(CON->get_tool_type()==2) Z2[i]+=4*le*B+Z_mod-0.5*Height;
				//Z2[i]-=4*le*B+le;//もう1粒子分下げて、ショルダーが流体部分をえぐらないようにする// Bをかけて1層分だと接触しすぎ？
			}/////////////////////////////////////////
			
			////////////////////////////ツールと被る流体粒子を削除する
			beforeN=number;
			make_fusion3D(X,Y,Z,X2,Y2,Z2,surface2,&number,le);
			////////////////////////////////////

			///////////////////box粒子の材質決定

			for(int i=beforeN;i<number;i++) type1.push_back(FRFLUID);//まずはとりあえずFRFLUIDを代入


			if(CON->get_tool_type()==0 ||CON->get_tool_type()==1)
			{
				Width-=6*le;
				Height-=4*le*B;
				Depth-=6*le*A;
				for(int i=beforeN;i<number;i++)
				{
					if(X[i]>(Width*0.5)+0.2*le)//右の壁
					{
						if(X[i]<(Width*0.5)+1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
					else if(X[i]<-(Width*0.5)-0.2*le)//左の壁
					{
						if(X[i]>-(Width*0.5)-1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
				
					if(Y[i]<-(depth0)-0.2*le)//手前の壁
					{
						//if(Y[i]>-(depth0)-1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]>-(depth0)-1.2*le && Z[i]>-1.2*le-Z_mod) 	
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;
						}
						else type1[i]=OUTWALL;
					}
					else if(Y[i]>(Depth-depth0)+0.2*le)//奥の壁
					{
						//if(Y[i]<(Depth-depth0)+1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]<(Depth-depth0)+1.2*le && Z[i]>-1.2*le-Z_mod)
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;// このifがないと、既に左右でoutwallと判断されている一部がinwallになってしまう
						}
						else type1[i]=OUTWALL;
					}

					if(CON->get_tool_type()==0 || CON->get_tool_type()==1)
					{
						if(X[i]<=(Width*0.5)+0.2*le && X[i]>=-(Width*0.5)-0.2*le)
						{
							if(Y[i]>=-(depth0)-0.2*le && Y[i]<=(Depth-depth0)+0.2*le)
							{
								if(Z[i]<-0.2*le-Z_mod)
								{
									if(Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;//下の壁
									else type1[i]=OUTWALL;
								}
							}
						}
					}
				}
			}
			if(CON->get_tool_type()==2)//表裏ツール
			{
				Width-=6*le;
				//Height-=4*le*B;
				Depth-=6*le*A;
				for(int i=beforeN;i<number;i++)
				{
					if(X[i]>(Width*0.5)+0.2*le)//右の壁
					{
						if(X[i]<(Width*0.5)+1.2*le) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
					else if(X[i]<-(Width*0.5)-0.2*le)//左の壁
					{
						if(X[i]>-(Width*0.5)-1.2*le) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
				
					if(Y[i]<-(depth0)-0.2*le)//手前の壁
					{
						//if(Y[i]>-(depth0)-1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]>-(depth0)-1.2*le) 	
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;
						}
						else type1[i]=OUTWALL;
					}
					else if(Y[i]>(Depth-depth0)+0.2*le)//奥の壁
					{
						//if(Y[i]<(Depth-depth0)+1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]<(Depth-depth0)+1.2*le)
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;// このifがないと、既に左右でoutwallと判断されている一部がinwallになってしまう
						}
						else type1[i]=OUTWALL;
					}
					////下側の壁は必要ない
				}
			}
			//////////////////////*/


			//粒子座標出力
			count2=0;
			double Zmax=0;
			for(int i=0;i<number;i++)
			{
				if(type1[i]==FRFLUID)
				{
					//if(X[i]>0) materialID=1;
					//else       materialID=2;
					materialID=1;//とりあえず同材接合
					double density,Cp,latent_H;
					if(materialID==1)
					{
						density=CON->get_density();
						Cp=CON->get_Cp();
						latent_H=CON->get_latent_H();
					}
					else if(materialID==2)
					{
						density=CON->get_density2();
						Cp=CON->get_Cp2();
						latent_H=CON->get_latent_H2();
					}
					double mass=density*vol;									//粒子質量
					//h=T*mass*Cp+latent_H*mass;//エンタルピー

					if(T>=CON->get_MP())
					{
						h=T*mass*Cp+latent_H*mass;
					}
					else h=T*mass*Cp;
					writedata2(fq, count2, X[i],Y[i],Z[i], FLUID,materialID,OFF, val,0,0,0,0,h,1);//ここまでの粒子はすべて流体

					count2++;
					if(Z[i]>=Zmax) Zmax=Z[i];
				}
			}

			//ツール部分の速度条件を与えた上で出力
			if(CON->get_tool_angle()==0)//ツールを回転させない場合
			{
				if(CON->get_tool_type()==0 || CON->get_tool_type()==1)
				{
					for(int i=0;i<number;i++)
					{
					
						if(type1[i]==INWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//ツールの中心が原点でないときは注意
							int AA=1;
							if(r>0.1*le && r<=shold_R+le && Z[i]>0)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								if(Z[i]<Height && r>probe_R-0.3*le) uw=-pich*rps;//INWALLのうち、プローブの側面のみ(底面は除く)に下方向速度追加
								if(calc_type==traverse) v+=U;//y方向にツールを移動
								if(calc_type==plunge) uw-=U;//Z方向にツールを移動(plange phase)
							}
							if(r<=shold_R+le && Z[i]>0) AA=MOVE;
							materialID=1;		
							writedata2(fq, count2, X[i],Y[i],Z[i], INWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
							count2++;
						}
					}
					for(int i=0;i<number;i++)
					{
						if(type1[i]==OUTWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//ツールの中心が原点でないときは注意
							int AA=1;
							if(r>0.1*le && r<=shold_R+le && Z[i]>0)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								uw=-pich*rps;
								if(calc_type==traverse) v+=U;//y方向にツールを移動
								else if(calc_type==plunge) uw-=U;//Z方向にツールを移動(plange phase)
							}
							if(r<=shold_R+le && Z[i]>0) AA=MOVE;
							materialID=1;
							writedata2(fq, count2, X[i],Y[i],Z[i], OUTWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);

							count2++;
						}
					}
				}
				else if(CON->get_tool_type()==2)
				{
					for(int i=0;i<number;i++)
					{
					
						if(type1[i]==INWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//ツールの中心が原点でないときは注意
							int AA=1;
							if(r>0.1*le && r<=shold_R+le)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								//if(Z[i]<Height && r>probe_R-0.3*le) uw=-pich*rps;//INWALLのうち、プローブの側面のみ(底面は除く)に下方向速度追加
								if(calc_type==traverse) v+=U;//y方向にツールを移動
								if(calc_type==plunge) uw-=U;//Z方向にツールを移動(plange phase)
							}
							if(r<=shold_R+le) AA=MOVE;
							materialID=1;		
							writedata2(fq, count2, X[i],Y[i],Z[i], INWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
							count2++;
						}
					}
					for(int i=0;i<number;i++)
					{
						if(type1[i]==OUTWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//ツールの中心が原点でないときは注意
							int AA=1;
							if(r>0.1*le && r<=shold_R+le)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								uw=-pich*rps;
								if(calc_type==traverse) v+=U;//y方向にツールを移動
								else if(calc_type==plunge) uw-=U;//Z方向にツールを移動(plange phase)
							}
							if(r<=shold_R+le) AA=MOVE;
							materialID=1;
							writedata2(fq, count2, X[i],Y[i],Z[i], OUTWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);

							count2++;
						}
					}
				}
			}
			else if(CON->get_tool_angle()>0)//ツールをx軸まわりに回転させる。進行方向が+y
			{
				double theta=PI/180*CON->get_tool_angle();//回転する角度
				
				for(int i=0;i<number;i++)
				{
					if(type1[i]==INWALL)
					{
						double speed=0;
						double u=0;
						double v=0;
						double uw=0;
						int AA=1;

						if(i<tool_number)
						{
							double r=sqrt(X3[i]*X3[i]+Y3[i]*Y3[i]);//ツールの中心が原点でないときは注意
							if(r>0.1*le && r<=shold_R+le && Z3[i]>0)
							{
								speed=r*w;
								u=speed*(-Y3[i]/r);
								v=speed*(X3[i]/r);
								if(Z[i]<Height && r>probe_R-0.3*le) uw=-pich*rps;//INWALLのうち、プローブの側面のみ(底面は除く)に下方向速度追加
								//ここまでの速度条件は、ツールの回転と共に回転するので座標と同様に回転させる

								double u_t=u;
								double v_t=cos(theta)*v-sin(theta)*uw;
								double uw_t=sin(theta)*v+cos(theta)*uw;
			
								u=u_t;
								v=v_t;
								uw=uw_t;

								if(calc_type==traverse) v+=U;//y方向にツールを移動
								if(calc_type==plunge) uw-=U;//Z方向にツールを移動(plange phase)
							}
							if(r<=shold_R+le && Z3[i]>0) AA=MOVE;
						}

						materialID=1;		
						writedata2(fq, count2, X[i],Y[i],Z[i], INWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
						count2++;
					}
				}
				//cout<<"inwall完了"<<endl;

				for(int i=0;i<number;i++)
				{
					if(type1[i]==OUTWALL)
					{
						//cout<<"out="<<i<<endl;
						double speed=0;
						double u=0;
						double v=0;
						double uw=0;
						
						int AA=1;

						if(i<tool_number)
						{
							double r=sqrt(X3[i]*X3[i]+Y3[i]*Y3[i]);//ツールの中心が原点でないときは注意
						
							if(r>0.1*le && r<=shold_R+le && Z3[i]>0)
							{
								speed=r*w;
								u=speed*(-Y3[i]/r);
								v=speed*(X3[i]/r);
								uw=-pich*rps;

								double u_t=u;
								double v_t=cos(theta)*v-sin(theta)*uw;
								double uw_t=sin(theta)*v+cos(theta)*uw;
			
								u=u_t;
								v=v_t;
								uw=uw_t;

								if(calc_type==traverse) v+=U;//y方向にツールを移動
								else if(calc_type==plunge) uw-=U;//Z方向にツールを移動(plange phase)
							}
							if(r<=shold_R+le && Z3[i]>0) AA=MOVE;
						}
							
						materialID=1;
						writedata2(fq, count2, X[i],Y[i],Z[i], OUTWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
						count2++;
					}
						
				}
			}

		}
		else while(1) cout<<"model次元エラー"<<endl;
	}


	else while(1) cout<<"modelエラー"<<endl;
	fq.close();


	//gnuplot用出力
	ofstream fp("plot.dat");
	for(int i=0;i<number;i++) fp<<X[i]<<" "<<Y[i]<<" "<<Z[i]<<endl;
	fp.close();
	/////////////////////
	
	*particle_number=number;
	ofstream fn("particle_number.dat");
	fn<<number<<endl;
	fn.close();
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

void writedata2(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM)
{	
	//ファイル出力
	fp<<id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<P<<"\t"<<h<<"\t"<<val<<"\t"<<type<<"\t"<<materialID<<"\t"<<surface<<"\t"<<toBEM<<"\t"<<endl;
}

/*void writedata(ofstream &fp, int number, double x, double y,double z, int type,double vx,double vy,double vz,double P,double h,int toFEM)
{	
	//fprintf( fp, "%5.10f\t", 0); とかいう数字の放り込みはだめ？
	
	double angle=0;//回転角
	double angle2=0;
	double angle3=0;
	double angle_s=1;
	double anglar_u=0;//角速度
	double anglar_u2=0;//角速度
	double anglar_u3=0;//角速度
	
	fp<<number<<"\t";
	fp<<x<<"\t";
	fp<<y<<"\t";
	fp<<z<<"\t";
	fp<<vx<<"\t";					//速度x成分
	fp<<vy<<"\t";					//速度y成分
	fp<<vz<<"\t";					//速度z成分
	fp<<P<<"\t";					//圧力
	fp<<h<<"\t";					//エンタルピー
	fp<<angle<<"\t";					//回転角
	fp<<type<<"\t";
	fp<<toFEM<<endl;
}*/


//半径Rの円の外周
void set_circle_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R)
{
	//対称性を考えた時、円周を記述する粒子は偶数個でなければならない。

	int N=calc_division_N_circle(2*PI*R,le);//円周の分割数
	
	double L=2*PI*R/N;				//粒子間距離
	double theta=L/R;				//粒子を配置する角度

	for(int n=0;n<N;n++)
	{
		X.push_back(R*cos(theta*n));
		Y.push_back(R*sin(theta*n));
		Z.push_back(0);
	}
	*number=*number+N;
}

//円周分割数計算関数
int calc_division_N_circle(double dis,double le)
{
	//対称性を考えた時、円周を記述する粒子は偶数個でなければならない。だから他の辺分割数とは扱いが少し特殊
	//dis:分割する距離(円周)
	double temp_num=dis/le;		//円外周に設置する『仮の』粒子数。ただし外周がうまくleで割り切れるとは限らない

	int N1=(int)(temp_num/2);
	N1*=2;							
	int N2=N1+2;					//temp_numはN1とN2の間にある。ここでN1,N2は偶数

	double dif1=temp_num-N1;		//各Nとの差
	double dif2=N2-temp_num;
	int N=N1;						//周方向分割数
	if(dif2<dif1) N=N2;				//差の小さい方をNとして採用する。

	return N;
}

//半径Rの円内部
void set_circle_in(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &type1,int *number,double le,double R,int edge_startID,int edge_lastID)
{
	//edge_startIDから(edge_lastID-1)までの粒子が、円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(R/le)+1;  //円を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //円を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-le*0.5;				//少し小さめの半径を設定
	
	//初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii*ii+jj*jj<R2*R2)
			{
				X.push_back(ii);
				Y.push_back(jj);
				Z.push_back(0);
				newN++;
			}
		}
	}

	//分子動力学により位置を最適化
	MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);

	*number=*number+newN;
}

//半径Rの円内部
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int edge_startID,int edge_lastID)
{
	//set_circle_in()と違い、60度分だけ計算し、それを6つｺﾋﾟｰして円を構成する。時間短縮と配置の均等化が目的
	//edge_startIDから(edge_lastID-1)までの粒子が、円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(R/le)+1;  //円を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //円を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-le*0.5;				//少し小さめの半径を設定


	//////////////////////円を6つにわけるための直線6つを生成

	double temp_R_num=R/le;			//半径方向に設置する『仮の』粒子数
	int    R_num=(int)temp_R_num;	//真の粒子数　とりあえず仮の粒子数の整数番とする。ここで、temp_R_num>R_numが成立している。
	double difference=temp_R_num-R_num;	//仮の数と真の数の差
	
	//仮の数と真の数の差が0.5までなら、真の数はNとする。0.5以上ならN+1とする
	if(difference>0.5) R_num++;
	double L=R/R_num;					//粒子間距離

	X.push_back(0);						//中心粒子を追加
	Y.push_back(0);
	Z.push_back(0);
	newN++;

	for(int k=0;k<6;k++)				//6つの直線のloop
	{
		double theta=PI/3*k;			//直線の角度
		for(int n=1;n<R_num;n++)		//中心粒子と最外周粒子はもうあるから、ここでのloopはそれをカウントしない
		{
			double r=L*n;				//中心からの距離
			X.push_back(r*cos(theta));
			Y.push_back(r*sin(theta));
			Z.push_back(0);
			newN++;
		}
	}
	*number=*number+newN;
	newN=0;
	beforeN=*number;
	edge_lastID=beforeN;
	/////////////////直線6つを生成完了

	//内部初期位置 ただし最初の1ピースのみ
	for(int i=0;i<=half_WX;i++)
	{
		for(int j=1;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii*ii+jj*jj<R2*R2)
			{
				if(jj<sqrt(3.0)*ii-le)		//最初のピースの斜め線より低い領域に設置。ただし直線ぎりぎりはまずいので、保険で-leしている
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(0);
					newN++;
				}
			}
		}
	}////////////////////////

	//分子動力学により位置を最適化
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	///内部粒子を周方向に6つｺﾋﾟｰ
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//回転する角度
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*X[i]-sin(theta)*Y[i];//回転後の座標
			double y=sin(theta)*X[i]+cos(theta)*Y[i];

			X.push_back(x);
			Y.push_back(y);
			Z.push_back(0);
			//newN++;
		}
	}///////////////////*/

	*number=*number+newN*6;//newNはひとつのピース内の粒子数を表しているからここでは6倍
}

//半径Rの球作成関数
void set_sphere(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag)
{
	//まずは半球を作る。そのためには半球表面を作成する必要がある。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数
	int beforeN2=*number;		//関数呼び出し時の粒子数。この関数の最後まで記憶しておく

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int half_WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int half_WZ=(int)(R/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.5*le;				//少し小さめの半径を設定

	///////////半球表面
	int Nt;						//球表面の、θ方向の分割数
	double Lt;					//球表面の、θ方向の分割距離
	calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//半円の分割数は偶数・奇数どちらでもよい
	double d_theta=Lt/R;		//弧の長さがLtになる角度

	for(int k=0;k<Nt;k++)//loopはk<Ntで終わらせる。Ntに該当するところはすでに設置済み
	{
		double THETA=k*d_theta;	//θ
		double r=R*sin(THETA);	//その高さにおける円の半径
		double round=2*PI*r;//その高さにおける円周

		int Nf=calc_division_N_circle(round,le);//球表面の、θ方向の分割数
		double Lf=round/Nf;						//球表面の、θ方向の分割距離
		double d_fai=Lf/r;						//弧の長さがLfになる角度
		
		for(int i=0;i<Nf;i++)
		{
			double fai=d_fai*i;
			if(Nt%2==0)
			{
				if(k%2!=0) fai+=0.5*d_fai;//Ntが偶数なら、作成済みの円と接するときは奇数番目。よって奇数をずらす
			}
			else
			{
				if(k%2==0) fai+=0.5*d_fai;//Ntが奇数なら、作成済みの円と接するときは偶数番目。よって奇数をずらす
			}
			double x=r*cos(fai);
			double y=r*sin(fai);
			double z=R*cos(THETA);
			X.push_back(x);
			Y.push_back(y);
			Z.push_back(z);
			
			newN++;
		}
	}
	if(Nt%2!=0)//Ntが奇数のときは、頂上に粒子が置かれなければならない。しかし上のloopはそれが不可能。よってここで追加
	{
		X.push_back(0);
		Y.push_back(0);
		Z.push_back(R);
		newN++;
	}
	//////////////////////////////////////

	*number=*number+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	
	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(ii*ii+jj*jj+kk*kk<R2*R2*0.6*0.6)//半径*0.7内の粒子は位置固定
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;

	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	///////////////////////*/

	

	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				//if(ii*ii+jj*jj+kk*kk<R2*R2)
				if(ii*ii+jj*jj+kk*kk<R2*R2 && ii*ii+jj*jj+kk*kk>=R2*R2*0.6*0.6)
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}///////////////////////*/
	

	//分子動力学により位置を最適化
	//cout<<"位置最適化(sphere)---------";
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-0.2*R; rigion[A_Z][1]=1.2*R;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;
	}
	//cout<<"完了"<<endl;
	///上半球を下半球へｺﾋﾟｰし球を完成
	if(flag==FULL)
	{
		newN=0;					//新しく追加する粒子数
		beforeN=*number;		//この時における粒子数

		for(int k=0;k<beforeN;k++)
		{
			if(Z[k]>0.4*le)
			{
				X.push_back(X[k]);
				Y.push_back(Y[k]);
				Z.push_back(-Z[k]);
				newN++;
			}
		}
		*number=*number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//下半球がほしいときに、つくった上半球を上下反転させる
	{
		for(int k=beforeN2;k<*number;k++) Z[k]*=-1;
	}
	
}

void set_sphere2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag,int *suf_num)
{
	//まずは半球を作る。そのためには半球表面を作成する必要がある。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数
	int beforeN2=*number;		//関数呼び出し時の粒子数。この関数の最後まで記憶しておく

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int half_WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int half_WZ=(int)(R/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.5*le;				//少し小さめの半径を設定

	///////////半球表面
	int Nt;						//球表面の、θ方向の分割数
	double Lt;					//球表面の、θ方向の分割距離
	calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//半円の分割数は偶数・奇数どちらでもよい
	double d_theta=Lt/R;		//弧の長さがLtになる角度

	for(int k=0;k<Nt;k++)//loopはk<Ntで終わらせる。Ntに該当するところはすでに設置済み
	{
		double THETA=k*d_theta;	//θ
		double r=R*sin(THETA);	//その高さにおける円の半径
		double round=2*PI*r;//その高さにおける円周

		int Nf=calc_division_N_circle(round,le);//球表面の、θ方向の分割数
		double Lf=round/Nf;						//球表面の、θ方向の分割距離
		double d_fai=Lf/r;						//弧の長さがLfになる角度
		
		for(int i=0;i<Nf;i++)
		{
			double fai=d_fai*i;
			if(Nt%2==0)
			{
				if(k%2!=0) fai+=0.5*d_fai;//Ntが偶数なら、作成済みの円と接するときは奇数番目。よって奇数をずらす
			}
			else
			{
				if(k%2==0) fai+=0.5*d_fai;//Ntが奇数なら、作成済みの円と接するときは偶数番目。よって奇数をずらす
			}
			double x=r*cos(fai);
			double y=r*sin(fai);
			double z=R*cos(THETA);
			X.push_back(x);
			Y.push_back(y);
			Z.push_back(z);
			
			newN++;
		}
	}
	if(Nt%2!=0)//Ntが奇数のときは、頂上に粒子が置かれなければならない。しかし上のloopはそれが不可能。よってここで追加
	{
		X.push_back(0);
		Y.push_back(0);
		Z.push_back(R);
		newN++;
	}
	//////////////////////////////////////

	*number=*number+newN;
	*suf_num=*suf_num+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	
	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(ii*ii+jj*jj+kk*kk<R2*R2*0.6*0.6)//半径*0.7内の粒子は位置固定
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;

	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	///////////////////////*/

	

	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				//if(ii*ii+jj*jj+kk*kk<R2*R2)
				if(ii*ii+jj*jj+kk*kk<R2*R2 && ii*ii+jj*jj+kk*kk>=R2*R2*0.6*0.6)
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}///////////////////////*/
	

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-0.2*R; rigion[A_Z][1]=1.2*R;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;
	}

	///上半球を下半球へｺﾋﾟｰし球を完成
	if(flag==FULL)
	{
		newN=0;					//新しく追加する粒子数
		beforeN=*number;		//この時における粒子数

		for(int k=0;k<beforeN;k++)
		{
			if(Z[k]>0.4*le)
			{
				X.push_back(X[k]);
				Y.push_back(Y[k]);
				Z.push_back(-Z[k]);
				newN++;
			}
		}
		*number=*number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//下半球がほしいときに、つくった上半球を上下反転させる
	{
		for(int k=beforeN2;k<*number;k++) Z[k]*=-1;
	}
	
}


//直線を分割するさいの最適な分割数と分割距離の算出関数
void calc_N_and_L(double dis,double le,int *N,double *L)
{
	double temp_N=dis/le;			//仮の分割数。leで割り切れたら一番いいけど、そうもいかないときがある
	int Ns=(int) temp_N;				//真の分割数
	double difference=temp_N-Ns;		//仮と真の差
	if(difference>0.5) Ns++;
	*L=dis/Ns;			//粒子の距離
	*N=Ns;
}

//分子動力学関数
void MD_2D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN)
{
	//分子動力学によりnewN個の粒子の位置を最適化　IDがBstartIDからBendIDまでのは境界粒子なので動かさない

	double region[2][2];	//解析領域

	/////////////////////解析領域の決定
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	for(int i=BstartID;i<BendID;i++)
	{
		if(X[i]<region[A_X][0]) region[A_X][0]=X[i];
		else if(X[i]>region[A_X][1]) region[A_X][1]=X[i];

		if(Y[i]<region[A_Y][0]) region[A_Y][0]=Y[i];
		else if(Y[i]>region[A_Y][1]) region[A_Y][1]=Y[i];
	}
	for(int D=0;D<2;D++)
	{
		region[D][0]-=5*le;	//少し領域を広めにとる
		region[D][1]+=5*le;
	}//////////////////////////

	//パラメータ
	double k0=1;
	double r=1.5;
	double dt=0.001;
	
	//力はax^3+bx^2+dの式を採用。文献[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]を参照
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;

	vector<double> Fx(newN);	//各粒子に働くX方向力
	vector<double> Fy(newN);	//各粒子に働くY方向力
	vector<double> Ax(newN,0);	//X方向加速度
	vector<double> Ay(newN,0);	//Y方向加速度
	vector<double> U(newN,0);	//X方向速度
	vector<double> V(newN,0);	//Y方向速度
	vector<double> visX(newN);	//X方向粘性係数
	vector<double> visY(newN);	//Y方向粘性係数

	//計算の高速化のために格子を形成 解析幅がr*leで割り切れるとは限らないので、はみ出したところは切り捨て。なので各軸とも正の方向には余裕を持つこと
	double grid_width=le*((int)(r+1));								//格子の幅。rを含む整数*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X方向の格子の個数
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//各内部粒子を含む格子番号
	vector<int> *MESH=new vector<int>[plane_SIZE];				//各メッシュに格納される粒子ID格納

	for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);	//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);	//Y方向に何個目の格子か
		int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
	{
		int i=beforeN+k;
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
		MESH[number].push_back(i);
		index[k]=number;
	}//////////////////////////////////////////

	
	//計算開始
	for(int t=0;t<100;t++)
	{
		if(t%10==0 &&t>0)//MESHを作り直す
		{
			//まずはMESHを一度破壊する。
			for(int n=0;n<plane_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
			{
				int xn=(int)((X[i]-region[A_X][0])/grid_width);	//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);	//Y方向に何個目の格子か
				int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
			{
				int i=beforeN+k;
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
				MESH[number].push_back(i);
				index[k]=number;
			}
		}////////////

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0;					//初期化
			int i=beforeN+k;					//対応する粒子番号
			double kx=0;						//X方向バネ係数
			double ky=0;
			int G_id=index[k];				//格納する格子番号
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=X[j]-X[i];
						double y=Y[j]-Y[i];
						double dis=sqrt(x*x+y*y);
						if(dis<r*le && dis!=0)			//このloopは自分自身も通過するから、dis!=0は必要
						{
							double F=a*dis*dis*dis+b*dis*dis+d;
							Fx[k]-=F*x/dis;					//Fの値が正のときは斥力なので、-=にする
							Fy[k]-=F*y/dis;
							double K=3*a*dis*dis+2*b*dis;//バネ係数　力の式の微分に相当
							K=sqrt(K*K);					//正の値が欲しい。だから負のときに備えて正に変換
							kx+=K*x*x/(dis*dis);			//kを各方向に分配。ここで、常に正の量が分配されるようにx*x/(dis*dis)となっている
							ky+=K*y*y/(dis*dis);
						}
					}
				}
			}
			visX[k]=1.414*sqrt(kx);//このように各軸方向の粘性係数を決める。文献「物理モデルによる自動メッシュ分割」P6参照。ただし質量は1としている。
			visY[k]=1.414*sqrt(ky);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
		}//各粒子の加速度が求まった。
		
		if(t==0)	//最初のｽﾃｯﾌﾟ時にdtを決定
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//最大加速度が求まった
			dt=sqrt(0.02*le/MaxAccel);
		}

		for(int k=0;k<newN;k++)//速度と位置の更新
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			X[i]+=dt*(U[k]+u)*0.5;
			Y[i]+=dt*(V[k]+v)*0.5;
		}

		//再近接距離がle以下の場合はこれを修正
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//対応する粒子番号
			int G_id=index[k];				//格納する格子番号
			double mindis=le;
			int J=k;						//最近接距離の相手粒子
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=X[j]-X[i];
						double y=Y[j]-Y[i];
						double dis=sqrt(x*x+y*y);
						if(dis<mindis && i!=j)
						{
							mindis=dis;
							J=j;
						}
					}
				}
			}
			if(J!=i && J<beforeN)//leより近接している相手が境界粒子なら
			{
				double L=le-mindis;//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;			
			}
			else if(J!=i && J>=beforeN)//leより近接している相手が内部粒子なら
			{
				double L=0.5*(le-mindis);//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;
				X[J]+=dX/mindis*L;
				Y[J]+=dY/mindis*L;
			}
		}//////////*/
	}/////MD終了

	delete [] index;
	delete [] MESH;
}

//分子動力学関数
void MD_3D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN,double r,double region[3][2])
{
	//分子動力学によりnewN個の粒子の位置を最適化　IDがBstartIDからBendIDまでのは境界粒子なので動かさない
	double k0=1;
	double dt=0.001;
	
	//力はax^3+bx^2+dの式を採用。文献[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]を参照
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;

	//cout<<"F="<<a*le*le*le+b*le*le+d<<" "<<a*1.5*le*1.5*le*1.5*le+b*1.5*le*1.5*le+d<<endl;

	vector<double> Fx(newN);	//各粒子に働くX方向力
	vector<double> Fy(newN);	//各粒子に働くY方向力
	vector<double> Fz(newN);	//各粒子に働くZ方向力
	vector<double> Ax(newN,0);	//X方向加速度
	vector<double> Ay(newN,0);	//Y方向加速度
	vector<double> Az(newN,0);	//Z方向加速度
	vector<double> U(newN,0);	//X方向速度
	vector<double> V(newN,0);	//Y方向速度
	vector<double> W(newN,0);	//Z方向速度
	vector<double> visX(newN);	//X方向粘性係数
	vector<double> visY(newN);	//Y方向粘性係数
	vector<double> visZ(newN);	//Y方向粘性係数
	vector<double> KX(newN);	//X方向バネ係数
	vector<double> KY(newN);	//Y方向バネ係数
	vector<double> KZ(newN);	//Y方向バネ係数

	//計算の高速化のために格子を形成 解析幅がr*leで割り切れるとは限らないので、はみ出したところは切り捨て。なので各軸とも正の方向には余裕を持つこと
	double grid_width=le*((int)(r+1));								//格子の幅。rを含む整数*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X方向の格子の個数
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//各内部粒子を含む格子番号
	vector<int> *MESH=new vector<int>[grid_SIZE];				//各メッシュに格納される粒子ID格納

	for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
	{
		int i=beforeN+k;
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		MESH[number].push_back(i);
		index[k]=number;
	}

	//計算開始
	for(int t=0;t<100;t++)
	{
		if(t%10==0 && t>0)
		{
			//MESHを一度破壊する。
			for(int n=0;n<grid_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
			{
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
			{
				int i=beforeN+k;
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
				MESH[number].push_back(i);
				index[k]=number;
			}
		}

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0, Fz[k]=0;			//初期化
			KX[k]=0;KY[k]=0; KZ[k]=0;			//バネ係数
		}

		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//対応する粒子番号
			int G_id=index[k];				//格納する格子番号
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							if(j>=beforeN && j>i)	//同じ領域内でかつiより大きな番号なら
							{
								int J=j-beforeN;	//newN内での番号
								double x=X[j]-X[i];
								double y=Y[j]-Y[i];
								double z=Z[j]-Z[i];
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le)			//このloopは自分自身も通過するから、dis!=0は必要
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//Fの値が正のときは斥力なので、-=にする
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									Fx[J]+=F*x/dis;					//相手粒子の力もここで計算。符号は反転させる
									Fy[J]+=F*y/dis;
									Fz[J]+=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//バネ係数　力の式の微分に相当
									K=sqrt(K*K);					//正の値が欲しい。だから負のときに備えて正に変換
									KX[k]+=K*x*x/(dis*dis);			//kを各方向に分配。ここで、常に正の量が分配されるようにx*x/(dis*dis)となっている
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
									KX[J]+=K*x*x/(dis*dis);			//kを相手粒子にも分配
									KY[J]+=K*y*y/(dis*dis);
									KZ[J]+=K*z*z/(dis*dis);
								}
							}
							if(j<BendID && j>=BstartID)
							{
								double x=X[j]-X[i];
								double y=Y[j]-Y[i];
								double z=Z[j]-Z[i];
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le && dis>0)			//このloopは自分自身は通過しない、dis!=0は不要
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//Fの値が正のときは斥力なので、-=にする
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//バネ係数　力の式の微分に相当
									K=sqrt(K*K);					//正の値が欲しい。だから負のときに備えて正に変換
									KX[k]+=K*x*x/(dis*dis);			//kを各方向に分配。ここで、常に正の量が分配されるようにx*x/(dis*dis)となっている
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
								}
							}
						}
					}
				}
			}
			//visX[k]=1.414*sqrt(KX[k]);//このように各軸方向の粘性係数を決める。文献「物理モデルによる自動メッシュ分割」P6参照。ただし質量は1としている。
			//visY[k]=1.414*sqrt(KY[k]);
			//visZ[k]=1.414*sqrt(KZ[k]);
			visX[k]=1.414*sqrt(KX[k]);//このように各軸方向の粘性係数を決める。文献「物理モデルによる自動メッシュ分割」P6参照。ただし質量は1としている。
			visY[k]=1.414*sqrt(KY[k]);
			visZ[k]=1.414*sqrt(KZ[k]);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
			Az[k]=(Fz[k]-visZ[k]*W[k]);
		}//各粒子の加速度が求まった。
		
		if(t==0)	//最初のｽﾃｯﾌﾟ時にdtを決定
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k]+Az[k]*Az[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//最大加速度が求まった
			if(MaxAccel!=0)
			{
				dt=sqrt(0.02*le/MaxAccel);
			}
		}

		for(int k=0;k<newN;k++)//速度と位置の更新
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			double w=W[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			W[k]+=dt*Az[k];
			X[i]+=dt*(U[k]+u)*0.5;
			Y[i]+=dt*(V[k]+v)*0.5;
			Z[i]+=dt*(W[k]+w)*0.5;
		}

		//再近接距離がle以下の場合はこれを修正
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//対応する粒子番号
			int G_id=index[k];				//格納する格子番号
			double mindis=le;
			int J=k;						//最近接距離の相手粒子
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							double x=X[j]-X[i];
							double y=Y[j]-Y[i];
							double z=Z[j]-Z[i];
							double dis=sqrt(x*x+y*y+z*z);
							if(dis<mindis && i!=j)
							{
								mindis=dis;
								J=j;
							}
						}
					}
				}
			}
			if(J!=i && J<beforeN)//leより近接している相手が境界粒子なら
			{
				double L=le-mindis;//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				double dZ=Z[J]-Z[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;	
				Z[i]-=dZ/mindis*L;	
			}
			else if(J!=i && J>=beforeN)//leより近接している相手が内部粒子なら
			{
				double L=0.5*(le-mindis);//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				double dZ=Z[J]-Z[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;
				Z[i]-=dZ/mindis*L;
				X[J]+=dX/mindis*L;
				Y[J]+=dY/mindis*L;
				Z[J]+=dZ/mindis*L;
			}
		}//////////*/
	}/////MD終了

	delete [] index;
	delete [] MESH;
}


//長方形の辺作成関数
void set_rectangular_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V)
{
	//水平方向長さLen_H,垂直長さLen_Vの長方形の辺を作成する
	
	int newN=0;					//この関数で追加される粒子数
	
	////////////////////まず水平方向作成

	int Nh;				//水平辺の分割数
	double dL_H;		//水平辺の分割長さ
	calc_N_and_L(Len_H,le,&Nh,&dL_H);

	for(int n=0;n<=Nh;n++)//底辺 ここではY,Z座標はゼロとしておく。
	{
		X.push_back(-Len_H*0.5+dL_H*n);
		Y.push_back(-Len_V*0.5);
		Z.push_back(0);
		newN++;
	}
	for(int n=0;n<=Nh;n++)//上辺
	{
		X.push_back(-Len_H*0.5+dL_H*n);
		Y.push_back(Len_V*0.5);
		Z.push_back(0);
		newN++;
	}////////////////////////////////////////////

	////////////////////次に垂直方向作成

	int Nv;				//水直辺の分割数
	double dL_V;		//水直辺の分割長さ
	calc_N_and_L(Len_V,le,&Nv,&dL_V);

	for(int n=1;n<Nv;n++)//左辺 ここではX座標はゼロとしておく。またn=0に該当する点はすでに設置済み
	{
		X.push_back(-Len_H*0.5);
		Y.push_back(-Len_V*0.5+n*dL_V);
		Z.push_back(0);
		newN++;
	}
	for(int n=1;n<Nv;n++)//右辺 
	{
		X.push_back(Len_H*0.5);
		Y.push_back(-Len_V*0.5+n*dL_V);
		Z.push_back(0);
		newN++;
	}////////////////////////////////////////////

	*number=*number+newN;
}

//長方形作成関数
void set_rectangular_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V,int edge_startID,int edge_lastID)
{
	//水平方向長さLen_H,垂直長さLen_Vの長方形の辺を作成する
	
	int newN=0;					//この関数で追加される粒子数
	int beforeN=*number;			//この関数呼び出し時にすでに存在している粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(Len_H/le)+1;  //長方形を十分含む幅
	int half_WY=(int)(Len_V/(le*A))+1;
	double gap=0.4*le;				//辺ぎりぎりに内部粒子を配置しないよう、隙間を設ける
	
	//初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WX;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii>-Len_H*0.5+gap && ii<Len_H*0.5-gap)
			{
				if(jj>-Len_V*0.5+gap && jj<Len_V*0.5-gap)
				{
					if(ii*ii<Len_H*0.25*Len_H*0.25 && jj*jj<Len_V*0.25*Len_V*0.25)//十分内部は配置を決め打ち
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(0);
						newN++;
					}
				}
			}
		}
	}
	*number=*number+newN;
	newN=0;					
	beforeN=*number;
	edge_lastID=*number;//境界IDを変更　ここまでに設置された粒子が境界粒子になる。

	//初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WX;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii>-Len_H*0.5+gap && ii<Len_H*0.5-gap)
			{
				if(jj>-Len_V*0.5+gap && jj<Len_V*0.5-gap)
				{
					if(ii*ii>=Len_H*0.25*Len_H*0.25 || jj*jj>=Len_V*0.25*Len_V*0.25)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(0);
						newN++;
					}
				}
			}
		}
	}

	//分子動力学により位置を最適化
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	*number=*number+newN;
}

//円柱表面作成関数
void set_cylinder_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int circle_start_id,int circle_end_id,int top_flag)
{
	//半径R,高さheightの円柱の面を作成する。ただしこの関数呼び出し時において、すでに下面の円(Z=0)は作成済みとする
	//top_flag=ONなら円柱上面を作成する。OFFならしないが、側面だけは作成する。
	int beforeN=*number;
	int newN=0;

	int Nv;				//水直の分割数
	double dL_V;		//水直の分割長さ
	double A=sqrt(3.0)/2;		//よく使う係数
	calc_N_and_L(height,le*A,&Nv,&dL_V);

	int Nr=calc_division_N_circle(2*PI*R,le);//円周の分割数
	double Lr=2*PI*R/Nr;				//円周分割距離

	double gap=0.4*le;				//辺ぎりぎりに内部粒子を配置しないよう、隙間を設ける

	///////////////////////////////////側面
	for(int i=0;i<Nr;i++)
	{
		for(int j=1;j<Nv;j++)//j=0,j=Nvは下面、上面に該当するのでここではぬかす
		{
			double jj=j*dL_V;
			double ii=i*Lr;
			if(j%2!=0) ii+=0.5*Lr;//jが奇数ならiiを0.5格子だけずらす
			if(ii<2*PI*R-gap)
			{
				if(jj<height-gap)	
				{
					double theta=2*PI*(ii/(2*PI*R));
					X.push_back(R*cos(theta));
					Y.push_back(R*sin(theta));
					Z.push_back(jj);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	////////////////////////


	//上面作成（下面のコピー）ただしNvが奇数なら上面は下面と半格子ずれなければならない
	beforeN=*number;
	newN=0;
	if(Nv%2==0)	//偶数ならそのままｺﾋﾟｰ
	{
		if(top_flag==ON)
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				X.push_back(X[i]);
				Y.push_back(Y[i]);
				Z.push_back(height);
				newN++;
			}
		}
		if(top_flag==HALF)//内壁部より内側はなし
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R && r<R+4*le-0.1*le)//外周のみ作成
				{
				X.push_back(X[i]);
				Y.push_back(Y[i]);
				Z.push_back(height);
				newN++;
				}
			}
		}
		else if(top_flag==OFF)//上面が必要ないなら
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R-0.1*le)//外周のみ作成
				{
					X.push_back(X[i]);
					Y.push_back(Y[i]);
					Z.push_back(height);
					newN++;
				}
			}
		}
	}
	else
	{
		double d_theta=0.5*Lr/R;//この微小角度だけ回転させる。
		if(top_flag==ON)
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//回転後の座標
				double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
				X.push_back(X2);
				Y.push_back(Y2);
				Z.push_back(height);
				newN++;
			}
		}
		else if(top_flag==HALF)//上面が必要ないなら
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				//if(r>R-0.1*le)//外周のみ作成
				if(r>R && r<R+4*le-0.1*le)//外周のみ作成
				{
					double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//回転後の座標
					double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
					X.push_back(X2);
					Y.push_back(Y2);
					Z.push_back(height);
					newN++;
				}
			}
		}
		else if(top_flag==OFF)//上面が必要ないなら
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R-0.1*le)//外周のみ作成
				{
					double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//回転後の座標
					double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
					X.push_back(X2);
					Y.push_back(Y2);
					Z.push_back(height);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	
	/////////////////////////////

}

//円錐表面作成関数
void set_circular_cone_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R2,double height,int circle_start_id,int circle_end_id,int top_flag)
{
	//上部半径R,下部半径R2,高さheightの円錐の面を作成する。ただしこの関数呼び出し時において、すでに下面の円(Z=0)は作成済みとする
	//top_flag=ONなら円柱上面を作成する。OFFならしないが、側面だけは作成する。
	int beforeN=*number;
	int newN=0;

	int Nv;				//水直の分割数
	double dL_V;		//水直の分割長さ
	double A=sqrt(3.0)/2;		//よく使う係数
	calc_N_and_L(height,le*A,&Nv,&dL_V);

	

	double gap=0.4*le;				//辺ぎりぎりに内部粒子を配置しないよう、隙間を設ける

	///////////////////////////////////側面
	
	for(int j=1;j<=Nv;j++)//j=0,j=Nvは下面、上面に該当するのでここではぬかす
	{
		double jj=j*dL_V;
		double r=jj*(R-R2)/height + R2;
		int Nr=calc_division_N_circle(2*PI*r,le);//円周の分割数
		double Lr=2*PI*r/Nr;				//円周分割距離

		for(int i=0;i<Nr;i++)
		{
			
			double ii=i*Lr;
			if(j%2!=0) ii+=0.5*Lr;//jが奇数ならiiを0.5格子だけずらす
			if(ii<2*PI*R-gap)
			{
				if(jj<height-gap)	
				{
					
					double theta=2*PI*(ii/(2*PI*r));
					X.push_back(r*cos(theta));
					Y.push_back(r*sin(theta));
					Z.push_back(jj);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	////////////////////////


	//上面作成
	beforeN=*number;
	newN=0;
	
		//if(top_flag==ON)
		//{
		//	for(int i=circle_start_id;i<circle_end_id;i++)
		//	{
		//		X.push_back(X[i]);
		//		Y.push_back(Y[i]);
		//		Z.push_back(height);
		//		newN++;
		//	}
		//}
		//if(top_flag==HALF)//内壁部より内側はなし
		//{
		//	for(int i=circle_start_id;i<circle_end_id;i++)
		//	{
		//		double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
		//		if(r>R && r<R+4*le-0.1*le)//外周のみ作成
		//		{
		//		X.push_back(X[i]);
		//		Y.push_back(Y[i]);
		//		Z.push_back(height);
		//		newN++;
		//		}
		//	}
		//}
		//else if(top_flag==OFF)//上面が必要ないなら
		//{
		//	for(int i=circle_start_id;i<circle_end_id;i++)
		//	{
		//		double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
		//		if(r>R-0.1*le)//外周のみ作成
		//		{
		//			X.push_back(X[i]);
		//			Y.push_back(Y[i]);
		//			Z.push_back(height);
		//			newN++;
		//		}
		//	}
		//}
	
	
	*number=*number+newN;
	
	/////////////////////////////

}

//円柱内部設置関数
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int flag)
{
	//半径R,高さheightの円柱内部を作成する。この関数呼び出し時に0<=i<numberの粒子で円柱表面が形成されているとする。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.3*le;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(ii*ii+jj*jj<R2*R2)
				{
					if(kk<height2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;

}

void set_crucible_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R_out,double height,double height_out,int flag,int fluid_number)
{
	//半径R,高さheightの円柱内部を作成する。この関数呼び出し時に0<=i<numberの粒子で円柱表面が形成されているとする。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int WX=(int)(R_out/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R_out/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R+0.3*le;				//少し小さめの半径を設定
	double R2_out=R_out-0.3*le;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=-WZ;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk>=0 && kk<height2)//内壁円柱部
				{
					if(ii*ii+jj*jj<R2_out*R2_out &&ii*ii+jj*jj>R2*R2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<0 && kk>-R2_out)//内壁円柱部
				{
					if(ii*ii+jj*jj+kk*kk<R2_out*R2_out &&ii*ii+jj*jj+kk*kk>R2*R2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				
			}
		}
	}///////////////////////*/

	///分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R_out; rigion[A_X][1]=1.2*R_out;
	rigion[A_Y][0]=-1.2*R_out; rigion[A_Y][1]=1.2*R_out;
	rigion[A_Z][0]=-1.2*R_out;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,fluid_number,beforeN,beforeN,newN,r,rigion);
	
	*number=*number+newN;

}

//ドーナツ作成
void set_doughnut2D(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_big,double R_smal,int edge_startID,int edge_lastID)
{
	//すでに作成された円周を用いて、それより大きな半径のドーナツを作成する。
	//edge_startIDから(edge_lastID-1)までの粒子が、すでに作成されている円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(R_big/le)+1;  //円を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R_big/(le*A))+1;  //円を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R_big-le*0.4;				//少し小さめの半径を設定
	double R1=R_smal+le*0.4;				//少し大きめの半径を設定


	//もうひとつの円周を作成 粒子数は内部で増加することに注意
	set_circle_edge(X,Y,Z,number,le,R_big);

	edge_lastID=*number;//MD2D()における境界粒子はedge_startID<=i<edge_lastID
	beforeN=*number;

	//内部初期位置 ただし最初の1ピースのみ
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii*ii+jj*jj<R2*R2 && ii*ii+jj*jj>R1*R1)
			{
				X.push_back(ii);
				Y.push_back(jj);
				Z.push_back(0);
				newN++;
			}
		}
	}////////////////////////

	//分子動力学により位置を最適化
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	*number=*number+newN;
}

//FSWプローブ内部粒子セット関数
void set_hat_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID)
{
	//帽子型の物質内部を作成する。例えばFSWのツール形状。
	
	//帽子の頭に該当する半径をR_smal,つばに該当する半径をR_big,つばの幅をH_flange,頭の幅をH_hat
	//ここで作成する帽子型の姿勢は、FSWのツールと同じで、頭を下にしてつばが上。
	//プローブ底辺のZ=0とする

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	double height=H_hat+H_flange;

	int WX=(int)(R_big/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R_big/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double gap=0.4*le;					//隙間
	double R_big2=R_big-gap;				//少し小さめの半径を設定
	double R_smal2=R_smal-gap;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定
	

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk<=H_hat+gap)
				{
					if(ii*ii+jj*jj<R_smal2*R_smal2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<height-gap)
				{
					if(ii*ii+jj*jj<R_big2*R_big2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R_big; rigion[A_X][1]=1.2*R_big;
	rigion[A_Y][0]=-1.2*R_big; rigion[A_Y][1]=1.2*R_big;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,bound_startID,bound_endID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//FSWプローブ内部粒子セット関数
void set_hat_in_2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_mid,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID)
{
	//帽子型の物質内部を作成する。例えばFSWのツール形状。
	
	//帽子のてっぺんに該当する半径をR_smal,かぶる部分に該当する半径をR_mid,つばに該当する半径をR_big,つばの幅をH_flange,頭の幅をH_hat
	//ここで作成する帽子型の姿勢は、FSWのツールと同じで、頭を下にしてつばが上。
	//プローブ底辺のZ=0とする

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	double height=H_hat+H_flange;

	int WX=(int)(R_big/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R_big/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double gap=0.4*le;					//隙間
	double R_big2=R_big-gap;				//少し小さめの半径を設定
	double R_mid2=R_mid-gap;				//少し小さめの半径を設定
	double R_smal2=R_smal-gap;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定
	

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk<=H_hat+gap)
				{
					double r=kk*(R_mid2-R_smal2)/H_hat + R_smal2;
					if(ii*ii+jj*jj<r*r)
					//if(ii*ii+jj*jj<R_smal2*R_smal2) double r=jj*(R-R2)/height + R2;
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<height-gap)
				{
					if(ii*ii+jj*jj<R_big2*R_big2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R_big; rigion[A_X][1]=1.2*R_big;
	rigion[A_Y][0]=-1.2*R_big; rigion[A_Y][1]=1.2*R_big;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,bound_startID,bound_endID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//箱作成関数
void set_box(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &surface,int *number,double le,double Width,double Height,double Depth)
{
	//横(Width)×高さ(Height)×奥行き(Depth)の箱を作成する。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ

	int BOX_startID=*number;	//この関数開始時の粒子数
	int beforeN=*number;
	int newN=0;
	double A=sqrt(3.0)/2;		//よく使う係数

	vector<double> X2;					//使い回し用
	vector<double> Y2;
	vector<double> Z2;
	int number2;						//使いまわしよう粒子数

	int Ndepth;					//奥行き方向の分割数
	int Nwidth;					//奥行き方向の分割数
	int Nheight;				//奥行き方向の分割数
	double dL_depth;			//行き方向の分割長さ
	double dL_width;			//行き方向の分割長さ
	double dL_height;			//行き方向の分割長さ

	double gap=0.4*le;
	
	calc_N_and_L(Depth,le,&Ndepth,&dL_depth);//各方向の分割数とその長さが求まる
	calc_N_and_L(Width,le,&Nwidth,&dL_width);
	calc_N_and_L(Height,le,&Nheight,&dL_height);

	//XY平面(上下面)作成//////長方形作成関数を使うのではなく、ここできちんと作成する。理由は、ひとつの長方形の辺を別の長方形が使用するから

	set_rectangular(X,Y,Z,number,le,Width,Depth);		//座標の原点は長方形の中心。よってあとで移動すること。
	for(int i=beforeN;i<*number;i++)					//重心移動
	{
		X[i]+=Width*0.5;
		Y[i]+=Depth*0.5;
	}

	for(int i=beforeN;i<*number;i++)					//上面にｺﾋﾟｰ
	{
		X.push_back(X[i]);
		Y.push_back(Y[i]);
		Z.push_back(Height);
		newN++;
	}
	*number=*number+newN;
	/////////////////////////////////////

	//XZ平面(正面・背面)作成
	number2=0;
	set_rectangular(X2,Y2,Z2,&number2,le,Width,Height);		//座標の原点は長方形の中心。よってあとで移動すること。
	for(int i=0;i<number2;i++)								//重心移動
	{
		X2[i]+=Width*0.5;
		Y2[i]+=Height*0.5;		//set_rectangular()は2D用なので、Zは値がゼロに注意
	}
		//粒子追加
	newN=0;
	beforeN=*number;
	for(int i=0;i<number2;i++)								//正式な座標に移動し、正式配列にｺﾋﾟｰ
	{
		if(Y2[i]>gap && Y2[i]<Height-gap)					//Y2=0,Y2=Heightの粒子はすでに作成済みなので省く
		{
			X.push_back(X2[i]);
			Y.push_back(0.0);			//Y=0すなわち正面
			Z.push_back(Y2[i]);			//Y2から値をもらうことに注意
			newN++;
		}
	}
	*number=*number+newN;				//正面粒子配置完了

	newN=0;
	for(int i=beforeN;i<*number;i++)
	{
		X.push_back(X[i]);
		Y.push_back(Depth);			//Y=Depthすなわち背面
		Z.push_back(Z[i]);
		newN++;
	}
	*number=*number+newN;				//背面粒子配置完了
	////////////////////////////////////////////////////////////


	//YZ平面(側面)作成
	newN=0;
	size_t vector_size=X2.size();
	for(int k=0;k<vector_size;k++)//X2,Y2,Z2をいったん消去
	{
		X2.pop_back();
		Y2.pop_back();
		Z2.pop_back();
	}

	number2=0;
	set_rectangular(X2,Y2,Z2,&number2,le,Depth,Height);		//座標の原点は長方形の中心。よってあとで移動すること。
	for(int i=0;i<number2;i++)								//重心移動
	{
		X2[i]+=Depth*0.5;
		Y2[i]+=Height*0.5;									//set_rectangular()は2D用なので、Zは値がゼロに注意
	}
		//粒子追加
	newN=0;
	beforeN=*number;
	for(int i=0;i<number2;i++)								//正式な座標に移動し、正式配列にｺﾋﾟｰ
	{
		if(X2[i]>gap && X2[i]<Depth-gap)					//辺粒子はすでに4辺とも作成ずみなので省く
		{
			if(Y2[i]>gap && Y2[i]<Height-gap)				
			{
				X.push_back(0.0);			//Z=0すなわち左面
				Y.push_back(X2[i]);			//X2から値をもらうことに注意
				Z.push_back(Y2[i]);			//Y2から値をもらうことに注意
				newN++;
			}
		}
	}
	*number=*number+newN;				//左面粒子配置完了

	newN=0;
	for(int i=beforeN;i<*number;i++)
	{
		X.push_back(Width);
		Y.push_back(Y[i]);			//Y=Depthすなわち背面
		Z.push_back(Z[i]);
		newN++;
	}
	*number=*number+newN;				//背面粒子配置完了
	////////////////////////////////////////////////////////////
	
	for(int i= BOX_startID;i<*number;i++)  surface.push_back(ON);//ここまでは表面
	beforeN=*number;

	set_box_in(X,Y,Z,number,le,Width,Height,Depth,BOX_startID,*number);
	
	for(int i=beforeN;i<*number;i++) surface.push_back(OFF);//内部
}

//長方形作成関数
void set_rectangular(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height)
{
	int BO_startID=*number;
	set_rectangular_edge(X,Y,Z,number,le,Width,Height);	//座標の原点は長方形の中心。よってあとで移動すること。
	int BO_lastID=*number;
	set_rectangular_in(X,Y,Z,number,le,Width,Height,BO_startID,BO_lastID);
}

//BOX内作成関数
void set_box_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height,double Depth,int BO_startID,int BO_lastID)
{
	//横(Width)×高さ(Height)×奥行き(Depth)の箱を作成する。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数
	
	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数

	int WX=(int)(Width/le)+1;		
	int WY=(int)(Depth/(le*A))+1; 
	int WZ=(int)(Height/(le*B))+1;  
	double gap=0.4*le;					//隙間
	
	//内部固定初期位置
	for(int i=1;i<=WX;i++)
	{
		for(int j=1;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk<Height*0.75 && jj<Depth*0.75 && ii<Width*0.75)
				{
					if(kk>Height*0.25 && jj>Depth*0.25 && ii>Width*0.25)
					{
						X.push_back(ii);//十分内部は決め打ち
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}
	*number=*number+newN;
	newN=0;					//この関数で新しく追加する粒子数
	beforeN=*number;
	///////////////////////*/

	//内部流動初期位置
	for(int i=1;i<=WX;i++)
	{
		for(int j=1;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk<Height-gap && jj<Depth-gap && ii<Width-gap)
				{
					if(kk<=Height*0.25 || jj<=Depth*0.25 || ii<=Width*0.25 || kk>=Height*0.75 || jj>=Depth*0.75 || ii>=Width*0.75)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}


	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-5*le; rigion[A_X][1]=1.5*Width;
	rigion[A_Y][0]=-5*le; rigion[A_Y][1]=1.5*Depth;
	rigion[A_Z][0]=-5*le; rigion[A_Z][1]=Height*1.5;

	MD_3D(X,Y,Z,le,BO_startID,BO_lastID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//物質合成関数
void make_fusion3D(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<double> &X2,vector<double> &Y2,vector<double> &Z2,vector<int> &surface2,int *number,double le)
{
	//存在が優先される粒子の座標がX,Y,Zに、存在が優先されない粒子の座標がX2,Y2,Z2に格納されている。
	//また、存在が優先されない粒子が、表面粒子か内部粒子かがsurface2[i]に格納されている。表面粒子は分子動力学で動かさないようにしないといけない

	size_t pri_num=X.size();			//優先物体の構成粒子数
	size_t neg_num=X2.size();			//消される物体の構成粒子数

	//cout<<X2.size()<<" "<<surface2.size()<<endl;

	double r=1.5;
	double region[3][2];	//解析領域
	int *flag=new int[neg_num];			//flagがONなら存在を許される。OFFなら消される
	int newN;
	int beforeN=*number;

	////////////////////解析領域の決定
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	region[A_Z][0]=100; region[A_Z][1]=-100;
	for(int i=0;i<pri_num;i++)
	{
		if(X[i]<region[A_X][0]) region[A_X][0]=X[i];
		else if(X[i]>region[A_X][1]) region[A_X][1]=X[i];
		if(Y[i]<region[A_Y][0]) region[A_Y][0]=Y[i];
		else if(Y[i]>region[A_Y][1]) region[A_Y][1]=Y[i];
		if(Z[i]<region[A_Z][0]) region[A_Z][0]=Z[i];
		else if(Z[i]>region[A_Z][1]) region[A_Z][1]=Z[i];
	}
	for(int i=0;i<neg_num;i++)
	{
		if(X2[i]<region[A_X][0]) region[A_X][0]=X2[i];
		else if(X2[i]>region[A_X][1]) region[A_X][1]=X2[i];
		if(Y2[i]<region[A_Y][0]) region[A_Y][0]=Y2[i];
		else if(Y2[i]>region[A_Y][1]) region[A_Y][1]=Y2[i];
		if(Z2[i]<region[A_Z][0]) region[A_Z][0]=Z2[i];
		else if(Z2[i]>region[A_Z][1]) region[A_Z][1]=Z2[i];
	}

	for(int D=0;D<3;D++)
	{
		region[D][0]-=5*le;//保険の意味で少し広めにとる
		region[D][1]+=5*le;
	}
	//////////////////////////////

	//計算の高速化のために格子を形成 解析幅がr*leで割り切れるとは限らないので、はみ出したところは切り捨て。なので各軸とも正の方向には余裕を持つこと
	double grid_width=le*((int)(r+1));								//格子の幅。rを含む整数*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X方向の格子の個数
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[neg_num];									//消される物体の構成粒子を含む格子番号
	vector<int> *MESH_pri=new vector<int>[grid_SIZE];				//各メッシュに格納される優先粒子ID格納
	//vector<int> *MESH_neg=new vector<int>[grid_SIZE];				//各メッシュに格納される非優先粒子ID格納

	for(int i=0;i<pri_num;i++)	//まずは優先粒子を格子に格納
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		MESH_pri[number].push_back(i);
	}
	for(int i=0;i<neg_num;i++)	//次に非優先粒子を格子に格納
	{
		int xn=(int)((X2[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y2[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z2[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		//MESH_neg[number].push_back(i);
		index[i]=number;
	}

	//flag[i]計算開始
	for(int i=0;i<neg_num;i++)
	{
		int G_id=index[i];				//格納する格子番号
		flag[i]=ON;
		for(int II=G_id-1;II<=G_id+1;II++)
		{       
			for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
			{
				for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
				{
					int M_id=II+JJ+KK;
					for(int L=0;L<MESH_pri[M_id].size();L++)//近隣の優先粒子を探索
					{
						int j=MESH_pri[M_id][L];
						
						double x=X[j]-X2[i];
						double y=Y[j]-Y2[i];
						double z=Z[j]-Z2[i];
						double dis=sqrt(x*x+y*y+z*z);
						if(dis<0.7*le)
						{
							flag[i]=OFF;	
						}
					}
				}
			}
		}
	}///////////////////

	//flag[i]=ONの粒子のみX,Y,Zに追加
	newN=0;
	for(int i=0;i<neg_num;i++)
	{
		if(flag[i]==ON && surface2[i]==ON)//flag=ONかつ表面粒子
		{
			X.push_back(X2[i]);
			Y.push_back(Y2[i]);
			Z.push_back(Z2[i]);
			newN++;
		}
	}
	*number=*number+newN;//分子動力学ではここまでの粒子が固定粒子扱い
	beforeN=*number;

	newN=0;
	for(int i=0;i<neg_num;i++)
	{
		if(flag[i]==ON && surface2[i]==OFF)//flag=ONかつ内部粒子
		{
			X.push_back(X2[i]);
			Y.push_back(Y2[i]);
			Z.push_back(Z2[i]);
			newN++;
		}
	}

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,region);

	*number=*number+newN;

	delete [] index;
	delete [] MESH_pri;
	//delete [] MESH_neg;
	delete [] flag;
	
}

