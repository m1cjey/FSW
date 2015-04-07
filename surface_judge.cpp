#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
//#include"define.h"	//#define 格納
//include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
//#include"PART.h"		//class PART定義
#include"BEMclass.h"	//BEM2D関係のclass 定義
//#include"FEM3Dclass.h"	//FEM3D関係のclass 定義
//#include<omp.h>
//#include<vector>
#include"function.h"

//近隣粒子探索関数
void calc_neighbor_relation(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int fluid_number,int out)
{
	//近隣粒子探索のための準備
	int *INDEX=new int[CON->get_number_of_mesh()];	//各格子に含まれる粒子数を格納
	reload_INDEX(CON,PART,particle_number,INDEX);//格子内の粒子数更新

	int **MESH = new int *[CON->get_number_of_mesh()];
	for(int i=0;i<CON->get_number_of_mesh();i++) MESH[i]=new int [INDEX[i]];
		
	reload_INDEX2(CON,PART,particle_number,MESH);
	////////////////////////*/

		
	unsigned int timeA=GetTickCount();
	if(CON->get_freeon()==1) freeon(CON,PART,particle_number,n0_4,INDEX,MESH,fluid_number,out);//表面判定
	else if(CON->get_freeon()==2) freeon2(CON,PART,particle_number,n0_4,INDEX,MESH,fluid_number,out);//表面判定
	//else if(CON->get_freeon()==4) freeon4(CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,n0,fluid_number,e0,out,t);
	else cout<<"表面判定未解決"<<endl;

	if(CON->get_surface_judge2()==ON && CON->get_adaptive_sw()==OFF) 
	{
		//surface_judge2(CON,PART,fluid_number,particle_number);
		surface_judge2_old(CON,PART,fluid_number,particle_number);
		//surface_judge2_new(CON,PART,fluid_number,particle_number);
	}
	
	delete [] INDEX;
	for(int i=0;i<CON->get_number_of_mesh();i++) delete [] MESH[i];
	delete [] MESH;
	cout<<"粒子依存関係計算終了--time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}


///表面判定関数
void freeon(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,int fluid_number,int out)
{
    ///説明:freeon関数では各粒子の粒子数密度を求めるとともに、得られた密度から表面判定を行う。また、最低粒子間距離もついでにもとめている
	///freeon2より遅いが、そのぶん並列化が容易。
	double le=CON->get_distancebp();
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	
	if(CON->get_T_field()==ON && CON->get_insulate()==1)
	{//非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
		out=particle_number;
	}
	
	//密度を求め、表面判定を行う。表面ならP=0にする
	//omp_set_num_threads(8);//ｽﾚｯﾄﾞ数指定
	#pragma omp parallel for
	for(int i=0;i<out;i++)//OUTWALL以外の粒子。//OUTWALLの粒子数密度などはいらない
	{    
		//printf("%d %d\n",i,omp_get_thread_num());//各iの計算を担当しているｽﾚｯﾄﾞ番号出力
	       
		PART[i].PND=0;//初期化
	    PART[i].PND2=0;
	    PART[i].N=0;
	    PART[i].N2=0;
	    PART[i].N3=0;
		////粒子数密度測定
		double pnd=0;//粒子数密度
		double pnd2=0;//ﾗﾌﾟﾗｼｱﾝ用粒子数密度
		double pnd4=0;//表面判定用
		int N=0;
		int N2=0;
		int N3=0;
		if(PART[i].index>=SIZE && PART[i].index<CON->get_number_of_mesh()-SIZE)
		{       
			for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
			{       
				for(int J=-1*CON->get_X_mesh();J<=CON->get_X_mesh();J+=CON->get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];
				     
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
					
							if(dis<=CON->get_re()*le && j!=i)//勾配・発散
							{       
								double r=CON->get_re()*le;
								double w=kernel(r,dis);
								pnd+=w;
								PART[i].NEI[N]=j;
								N++;
							}
							if(dis<=CON->get_re2()*le && j!=i)
							{       
								double r=CON->get_re2()*le;
								double w=kernel(r,dis);
								pnd2+=w;
								PART[i].NEI2[N2]=j;
								N2++;
							}
							if(dis<=CON->get_re3()*le && j!=i)//表面張力re3
							{       
								PART[i].NEI3[N3]=j;
								N3++;
							}
							if(dis<=CON->get_re4()*le && j!=i)
							{       
							    double r=CON->get_re4()*le;
								pnd4+=kernel(r,dis);
							}
						}
					}
				}
			}
		}
		PART[i].PND=pnd;
		PART[i].PND2=pnd2;
		PART[i].N=N;
		PART[i].N2=N2;
		PART[i].N3=N3;
		if(PART[i].N3>800) cout<<"error in freeon over. Change more than "<<PART[i].N3<<" coods:"<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		////////////////////

		if(PART[i].type==FLUID)
		{
			if(pnd4<n0_4*CON->get_beta())//β以下なら
			{
				PART[i].surface=ON;//表面粒子とする
				//PART[i].P=0;
			}
			else PART[i].surface=OFF;
		}
		else if(PART[i].type==INWALL)
		{
			if(pnd4<n0_4*CON->get_beta())//β以下なら
			{
				PART[i].surface=ON;//壁表面粒子とする
				//PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
	    
	}

	/*///最低粒子間距離をもとめる
	double min0=CON->get_distancebp();//最低粒子間距離
	int type1,type2,surface1,surface2;
	double X1,Y1,Z1;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min0)
			{
				type1=PART[i].type;
				surface1=PART[i].surface;
				type2=PART[j].type;
				surface2=PART[j].surface;
				min0=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}/////最低粒子間距離がもとまった
	cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
	//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
	*mindis=min0;
	///*/
}

///表面判定関数ver.2
void freeon2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,int fluid_number,int out)
{
	///freeon2関数の説明:flag1[i]の導入により高速化。ただし、1CPUなら早いが、ﾏﾙﾁCPUによる並列化は厳しい
	//cout<<"表面判定(freeon2)"<<endl;
	double le=CON->get_distancebp();
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	double d=2;
	if(CON->get_dimention()==3) d=3;

	if(CON->get_T_field()==ON && CON->get_insulate()==1)
	{
		out=particle_number;	//非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
	}

	/////壁境界ポリゴン化のテスト
	int wall_poly_number=0;//壁境界のポリゴン数
	wall_poly_number=3;//読み込みなどでもとめること
	double wall_a[3];//2次元の直線方程式ax+by+c=0　3DFEMを参考にnodeやedgeのような書き方をすること
	double wall_b[3];
	double wall_c[3];
	if(CON->get_wall_poly()==1)
	{
		//壁境界の設定 //メッシュなどから読み込むことを考えると、この位置はまずいかも
		if(CON->get_model_number()==1)
		{
			if(CON->get_dimention()==2)
			{
				wall_a[0]=1; wall_b[0]=0; wall_c[0]=0.011;//左側内壁
				wall_a[1]=1; wall_b[1]=0; wall_c[1]=-0.011;//右側内壁
				wall_a[2]=0; wall_b[2]=1; wall_c[2]=0.021;//下側内壁
			}
		}
	}
	////壁重み関数の計算
	double wallZ_pnd=0;//1E-05x-1.636
	double wallZ_pnd2=0;
	double wallZ_pnd4=0;//4E-06x-1.979
	///////////////
	double *PND4=new double [out];//表面判定用粒子数密度
	int *flag1=new int [out];		//検査フラグ。0:未検査　1:検査済み
	///初期化
	#pragma omp parallel for
	for(int i=0;i<out;i++)
	{
		PART[i].PND=0;//初期化
	    PART[i].PND2=0;
	    PART[i].N=0;
	    PART[i].N2=0;
	    PART[i].N3=0;

		PND4[i]=0;
		flag1[i]=0;
	}
	
	//密度を求め、表面判定を行う。表面ならP=0にする
	for(int i=0;i<out;i++)//OUTWALL以外のCFD粒子。//OUTWALLの粒子数密度などはいらない
	{    
		if(CON->get_wall_poly()==0)
		{
			////粒子数密度測定
			if(PART[i].index>=SIZE && PART[i].index<CON->get_number_of_mesh()-SIZE)
			{       
				for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
				{       
					for(int J=-1*CON->get_X_mesh();J<=CON->get_X_mesh();J+=CON->get_X_mesh())
					{
						for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
						{
							for(int L=0;L<INDEX[I+J+K];L++)
							{       
								int j=MESH[I+J+K][L];
								if(j<out)
								{
									if(flag1[j]==0 && j!=i)//まだ検査してないなら
									{
										double X=PART[j].r[A_X]-PART[i].r[A_X];
										double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
										double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
										double dis=sqrt(X*X+Y*Y+Z*Z);
							
										if(dis<=CON->get_re()*le)//勾配・発散
										{       
											double r=CON->get_re()*le;
											double w=kernel(r,dis);
											PART[i].PND+=w;
											PART[j].PND+=w;
											PART[i].NEI[PART[i].N]=j;
											PART[j].NEI[PART[j].N]=i;
											PART[i].N++;
											PART[j].N++;
										}
										if(dis<=CON->get_re2()*le)
										{       
											double r=CON->get_re2()*le;
											double w=kernel(r,dis);
							
											PART[i].PND2+=w;
											PART[j].PND2+=w;
										
											PART[i].NEI2[PART[i].N2]=j;
											PART[j].NEI2[PART[j].N2]=i;
											PART[i].N2++;
											PART[j].N2++;
										}
										if(dis<=CON->get_re3()*le)//表面張力re3
										{       
											PART[i].NEI3[PART[i].N3]=j;
											PART[j].NEI3[PART[j].N3]=i;
											PART[i].N3++;
											PART[j].N3++;
										}
										if(dis<=CON->get_re4()*le)
										{       
											double r=CON->get_re4()*le;
											double w=kernel(r,dis);
											PND4[i]+=w;
											PND4[j]+=w;
										}
									}
								}
								else
								{
									double X=PART[j].r[A_X]-PART[i].r[A_X];
									double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
									double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
									double dis=sqrt(X*X+Y*Y+Z*Z);
							
									if(dis<=CON->get_re()*le)//勾配・発散
									{       
										double r=CON->get_re()*le;
										double w=kernel(r,dis);
										PART[i].PND+=w;
										PART[i].NEI[PART[i].N]=j;
										PART[i].N++;
									}
									if(dis<=CON->get_re2()*le)
									{       
										double r=CON->get_re2()*le;
										double w=kernel(r,dis);
										PART[i].PND2+=w;
										PART[i].NEI2[PART[i].N2]=j;
										PART[i].N2++;
									}
									if(dis<=CON->get_re3()*le)//表面張力re3
									{       
										PART[i].NEI3[PART[i].N3]=j;
										PART[i].N3++;
									}
									if(dis<=CON->get_re4()*le)
									{       
										double r=CON->get_re4()*le;
										double w=kernel(r,dis);
										PND4[i]+=w;
									}
								}
							}
						}
					}
				}
			}
		}

		if(CON->get_wall_poly()==1)//ポリゴン壁。ここでは流体粒子同士の相互作用だけ計算される
		{
			////粒子数密度測定
			if(PART[i].index>=SIZE && PART[i].index<CON->get_number_of_mesh()-SIZE)
			{       
				for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
				{       
					for(int J=-1*CON->get_X_mesh();J<=CON->get_X_mesh();J+=CON->get_X_mesh())
					{
						for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
						{
							for(int L=0;L<INDEX[I+J+K];L++)
							{       
								int j=MESH[I+J+K][L];
								if(j<out)
								{
									if(flag1[j]==0 && j!=i)//まだ検査してないなら
									{
										double X=PART[j].r[A_X]-PART[i].r[A_X];
										double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
										double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
										double dis=sqrt(X*X+Y*Y+Z*Z);
							
										if(dis<=CON->get_re()*le)//勾配・発散
										{       
											double r=CON->get_re()*le;
											double w=kernel(r,dis);
											PART[i].PND+=w;
											PART[j].PND+=w;
											PART[i].NEI[PART[i].N]=j;
											PART[j].NEI[PART[j].N]=i;
											PART[i].N++;
											PART[j].N++;
										}
										if(dis<=CON->get_re2()*le)
										{       
											double r=CON->get_re2()*le;
											double w=kernel(r,dis);
							
											PART[i].PND2+=w;
											PART[j].PND2+=w;
										
											PART[i].NEI2[PART[i].N2]=j;
											PART[j].NEI2[PART[j].N2]=i;
											PART[i].N2++;
											PART[j].N2++;
										}
										if(dis<=CON->get_re3()*le)//表面張力re3
										{       
											PART[i].NEI3[PART[i].N3]=j;
											PART[j].NEI3[PART[j].N3]=i;
											PART[i].N3++;
											PART[j].N3++;
										}
										if(dis<=CON->get_re4()*le)
										{       
											double r=CON->get_re4()*le;
											double w=kernel(r,dis);
											PND4[i]+=w;
											PND4[j]+=w;
										}
									}
								}
								else
								{
									double X=PART[j].r[A_X]-PART[i].r[A_X];
									double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
									double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
									double dis=sqrt(X*X+Y*Y+Z*Z);
							
									if(dis<=CON->get_re()*le)//勾配・発散
									{       
										double r=CON->get_re()*le;
										double w=kernel(r,dis);
										PART[i].PND+=w;
										PART[i].NEI[PART[i].N]=j;
										PART[i].N++;
									}
									if(dis<=CON->get_re2()*le)
									{       
										double r=CON->get_re2()*le;
										double w=kernel(r,dis);
										PART[i].PND2+=w;
										PART[i].NEI2[PART[i].N2]=j;
										PART[i].N2++;
									}
									if(dis<=CON->get_re3()*le)//表面張力re3
									{       
										PART[i].NEI3[PART[i].N3]=j;
										PART[i].N3++;
									}
									if(dis<=CON->get_re4()*le)
									{       
										double r=CON->get_re4()*le;
										double w=kernel(r,dis);
										PND4[i]+=w;
									}
								}
							}
						}
					}
				}
			}
		}
		flag1[i]=1;//検査終了

		if(CON->get_wall_poly()==1)//ポリゴン壁。ここでは壁からの作用が計算される
		{
			double riw=100*le;//壁境界と粒子の距離
			double dis_w=0;
			for(int p=0; p<3;p++)//各壁境界からの粒子との距離を求める。pの値はポリゴン数を参照できるようにすること
			{
				dis_w=abs(wall_a[p]*PART[i].r[A_X]+wall_b[p]*PART[i].r[A_Y]+wall_c[p])/(sqrt(wall_a[p]*wall_a[p]+wall_b[p]*wall_b[p]));
				riw=dis_w;
					/*/
				if(dis_w<=riw)
				{
					riw=dis_w;

				}
				/*/
				if(riw<=CON->get_re()*le) PART[i].PND+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
				if(riw<=CON->get_re2()*le) PART[i].PND2+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
				if(riw<=CON->get_re4()*le)//pow(riw,-1.979)*pow(10.0,-6)*4;//4E-06x-1.979//= -9.245ln(x) - 55.708
				{
					if(-9.245*log(riw)-55.708>=0) PND4[i]+=-9.245*log(riw)-55.708;
				}
			}
			/*////
			if(riw<=CON->get_re()*le) PART[i].PND+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
			if(riw<=CON->get_re2()*le) PART[i].PND2+=pow(riw,-1.636)*pow(10.0,-5);//1E-05x-1.636
			if(riw<=CON->get_re4()*le)//pow(riw,-1.979)*pow(10.0,-6)*4;//4E-06x-1.979//= -9.245ln(x) - 55.708
			{
				if(-9.245*log(riw)-55.708>=0) PND4[i]+=-9.245*log(riw)-55.708;
			}
			///*/

		}

		if(PART[i].N3>450) cout<<"error in freeon over. Change more than "<<PART[i].N3<<" coods:"<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		////////////////////

		if(PART[i].type==FLUID)
		{
			if(PND4[i]<n0_4*CON->get_beta())//β以下なら
			{
				PART[i].surface=ON;//表面粒子とする
				PART[i].P=0;
			}
			else PART[i].surface=OFF;
		}
		else if(PART[i].type==INWALL)
		{
			if(PND4[i]<n0_4*CON->get_beta())//β以下なら
			{
				PART[i].surface=ON;//壁表面粒子とする
				PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
	}

	/*///最低粒子間距離をもとめる
	double min0=CON->get_distancebp();//最低粒子間距離
	int type1,type2,surface1,surface2;
	type1=0;
	type2=0;
	surface1=0;
	surface2=0;
	double X1,Y1,Z1;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min0)
			{
				type1=PART[i].type;
				surface1=PART[i].surface;
				type2=PART[j].type;
				surface2=PART[j].surface;
				min0=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}/////最低粒子間距離がもとまった
	cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
	//////////
	*mindis=min0;
	

	if(t==1)
	{
		ofstream fouts("mindis.dat");
		fouts.close();
	}

	ofstream fouts2("mindis.dat",ios :: app);
	fouts2<<min0<<endl;
	fouts2.close();
	/////*/

	/*/粒子数密度出力
	plot_PND(CON,PART,fluid_number,PND4,t);
	if(CON->get_PND_interval()>0)
	{
		if(t==1 || t%CON->get_PND_interval()==0) plot_PND_each(CON,PART,fluid_number,PND4,t);
	}
	///*/


	delete [] PND4;
	delete [] flag1;
}

//表面判定関数ver.2
void surface_judge2_old(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//従来の表面判定に加えて、ここでさらにふるいにかける。
	//粒子iの法線ベクトルと、粒子j方向ﾍﾞｸﾄﾙとの角度がある値を超えていたら表面ではないと判断
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	for(int i=0;i<particle_number;i++)
    {
       if(PART[i].type==FLUID ||PART[i].type==INWALL) PART[i].surface==ON;//粒子数密度による表面判定を用いない
	}
	
    //////法線ﾍﾞｸﾄﾙ計算
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)//この表面粒子が本当に表面粒子足りうるか判定する
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位ﾍﾞｸﾄﾙとの内積
					if(inp<-0.5) flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
				}
			}
			PART[i].surface=flag;
		}
	}

	//壁
	//int flag=OFF;
	//if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
	//if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM用のメッシュを作成する際に、すべての壁粒子を出力する場合も、OUTWALLの表面判定をする必要があるため、out=particle_numberとする。

	//if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].surface==ON)//この表面粒子が本当に表面粒子足りうるか判定する
			{
				int flag=ON;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(j>=fluid_number)							//壁粒子のsurface_jusge2は壁粒子だけで行う
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
						double ny=Y/dis;
						double nz=Z/dis;
						double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位ﾍﾞｸﾄﾙとの内積
						if(inp<-0.5) flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
					}
				}
				PART[i].surface=flag;
			}
		}
	}
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
}

//表面判定関数ver.2
void surface_judge2(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//従来の表面判定に加えて、ここでさらにふるいにかける。
	//粒子iの法線ベクトルと、粒子j方向ﾍﾞｸﾄﾙとの角度がある値を超えていたら表面ではないと判断
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	
    //////法線ﾍﾞｸﾄﾙ計算
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
		//if(PART[i].surface==ON)  direct_f2(CON,PART,i,direct);//流体のみが存在する条件下で求める
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}
	double le=CON->get_distancebp();
	//法線ベクトルを出力
	ofstream nn("normal.dat");
	if(CON->get_dimention()==2) for(int i=0;i<particle_number;i++) nn<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<direct[A_X][i]*le<<" "<<direct[A_Y][i]*le<<endl;
	nn.close();
	
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==ON)//この表面粒子が本当に表面粒子足りうるか判定する
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				//if(j<fluid_number)
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(PART[j].type==FLUID || PART[j].type==INWALL)
				//if(PART[j].type==FLUID || (PART[j].type==INWALL && dis<1.5*CON->get_distancebp()))
				//if(PART[j].type==FLUID)// テスト　壁付近の判定がマシになるか？
				{
					
					double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位ﾍﾞｸﾄﾙとの内積
					if(inp<-0.6) 
					{
						flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
						//cout<<"judge2で判定変換 i="<<i<<endl;
					}
				}
			}
			PART[i].surface=flag;
		}
	}
	
	//壁
	int flag=OFF;
	//if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
	//if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM用のメッシュを作成する際に、すべての壁粒子を出力する場合も、OUTWALLの表面判定をする必要があるため、out=particle_numberとする。

	////
	if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL)//この表面粒子が本当に表面粒子足りうるか判定する
			{
				if(PART[i].surface==ON)//この表面粒子が本当に表面粒子足りうるか判定する
				{
					int flag=ON;
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(j>=fluid_number)							//壁粒子のsurface_jusge2は壁粒子だけで行う なぜ？
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
							double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
							double ny=Y/dis;
							double nz=Z/dis;
							double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位ﾍﾞｸﾄﾙとの内積
							if(inp<-0.5) flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
						}
					}
					PART[i].surface=flag;
				}
			}
		}
	}
	////*/

	/*//////
	//孤立して表面粒子なものは、内部粒子と判定する
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			int flag=OFF;					//ONなら表面　OFFなら内部
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==FLUID &&PART[j].surface==ON && dis<CON->get_re()*CON->get_distancebp()) flag=ON;	//周囲に表面粒子がいればそれでよし
			}
			PART[i].surface=flag;
			//if(flag==OFF) cout<<"judge2で孤立表面判定変換 i="<<i<<endl; //これを入れると止まる？

		}
	}//////////////*/

	/*/孤立して表面粒子なINWALLは、内部粒子と判定する
	for(int i=fluid_number;i<particle_number;i++)
	{
		if(PART[i].type==INWALL && PART[i].surface==ON)
		{
			int flag=OFF;					//ONなら表面　OFFなら内部
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==INWALL &&PART[j].surface==ON) flag=ON;	//周囲に表面粒子がいればそれでよし
			}
			PART[i].surface=flag;
			
		}
	}//////////////*/
	
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];

}

////
//表面判定関数
void surface_judge2_new(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{

	//粒子iの法線ベクトルと、粒子j方向ﾍﾞｸﾄﾙとの角度がある値を超えていたら表面ではないと判断
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	for(int i=0;i<particle_number;i++)
    {
       if(PART[i].type==FLUID ||PART[i].type==INWALL) PART[i].surface==ON;//粒子数密度による表面判定を用いない
	}

    //////法線ﾍﾞｸﾄﾙ計算
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}
	
	int d=CON->get_dimention();

	double *L=new double [fluid_number];		//この関数で使用する基準粒子間距離

	for(int i=0;i<fluid_number;i++)  L[i]=PART[i].L;

	/*////
	for(int i=0;i<fluid_number;i++)
	{
		int jnb=PART[i].N;
		if(jnb>4)
		{
			for(int n=0;n<4;n++)
			{
				int temp_k=n;
				int A=PART[i].NEI[n];
				int temp_j=A;
				double X=PART[i].r[A_X]-PART[A].r[A_X];
				double Y=PART[i].r[A_Y]-PART[A].r[A_Y];
				double Z=PART[i].r[A_Z]-PART[A].r[A_Z];
				double mindis=X*X+Y*Y+Z*Z;
					
				for(int k=n+1;k<jnb;k++) 
				{
					int j=PART[i].NEI[k];
					X=PART[j].r[A_X]-PART[i].r[A_X];
					Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=X*X+Y*Y+Z*Z;
					if(dis<mindis)
					{
						mindis=dis; 
						temp_k=k;
						temp_j=j;
					}
				}
				PART[i].NEI[temp_k]=A;
				PART[i].NEI[n]=temp_j;
			}
			double ave_dis=0;
			for(int n=0;n<4;n++)
			{
				int j=PART[i].NEI[n];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				ave_dis+=dis*0.25;
			}
			L[i]=ave_dis;
			//L[i]=PART[i].L;
			//cout<<i<<" "<<L[i]/PART[i].L<<endl;
		}
	}////////////*/

	
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			int flag=ON;					//ONなら表面　OFFなら内部
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(PART[j].type==FLUID || (PART[j].type==INWALL && dis<1.2*PART[i].L))		//壁粒子の場合は接近の定義を狭める
				{
					//if(dis<CON->get_re()*L[i])
					{
					double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位相対距離ﾍﾞｸﾄﾙとの内積
					if(inp<-0.6) flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
					}
				}
			}
			PART[i].surface=flag; 
		}
	}


	//孤立して表面粒子なものは、内部粒子と判定する
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			int flag=OFF;					//ONなら表面　OFFなら内部
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==FLUID &&PART[j].surface==ON && dis<CON->get_re()*L[i]) flag=ON;	//周囲に表面粒子がいればそれでよし
			}
			PART[i].surface=flag;
		}
	}//////////////*/
	
	

	//壁
	//int flag=OFF;
	//if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
	//if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM用のメッシュを作成する際に、すべての壁粒子を出力する場合も、OUTWALLの表面判定をする必要があるため、out=particle_numberとする。

	//if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL)//この表面粒子が本当に表面粒子足りうるか判定する
			{
				int flag=ON;
				double direct2[3];
				for(int D=0;D<d;D++) direct2[D]=direct[D][i]*(-1);	//外向き法線 
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(j<fluid_number && PART[j].surface==OFF)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						if(dis<1.5*CON->get_distancebp())
						{
							X/=dis; Y/=dis; Z/=dis;
							double COS=X*direct2[A_X]+Y*direct2[A_Y]+Z*direct2[A_Z];
							if(COS>sqrt(3.0)*0.5) flag=OFF;
							//flag=OFF;	//内部流体粒子の近くにある壁は内部
						}
					}
				}
				PART[i].surface=flag;
			}
			//else if(PART[i].type==OUTWALL) PART[i].surface=ON;	//OUTWALLはとりあえず表面にしておく
		}
	}

	//孤立して表面粒子なINWALLは、内部粒子と判定する
	for(int i=fluid_number;i<particle_number;i++)
	{
		if(PART[i].type==INWALL && PART[i].surface==ON)
		{
			int flag=OFF;					//ONなら表面　OFFなら内部
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(PART[j].type==INWALL &&PART[j].surface==ON) flag=ON;	//周囲に表面粒子がいればそれでよし
			}
			PART[i].surface=flag;
		}
	}//////////////*/

	//内部粒子と判定されたら、法線ベクトルをゼロにし、かつ圧力をゼロにしておく
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].surface==OFF)
		{
			//for(int D=0;D<d;D++) PART[i].direct[D]=0;		//adaptive2のことも考えると、ここで法線を消さないほうがよい
			PART[i].P=0;
		}
	}

	delete [] L;
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
}
///*/

//表面判定関数ver.3
void surface_judge3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//従来の表面判定に加えて、ここでさらにふるいにかける。
	//粒子iの法線ベクトルと、粒子j方向ﾍﾞｸﾄﾙとの角度がある値を超えていたら表面ではないと判断
	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];
		
	
    //////法線ﾍﾞｸﾄﾙ計算
   // for(int i=0;i<fluid_number;i++)
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)//この表面粒子が本当に表面粒子足りうるか判定する
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位ﾍﾞｸﾄﾙとの内積
					if(inp<-0.5) flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
				}
			}
			PART[i].surface=flag;
		}
	}

	//壁
	int flag=OFF;
	if(CON->get_T_field()==ON && CON->get_insulate()==1) flag=ON; //非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
	if(CON->get_EM_method()==2 && CON->get_output_wall_way()==1) flag=ON;//BEM用のメッシュを作成する際に、すべての壁粒子を出力する場合も、OUTWALLの表面判定をする必要があるため、out=particle_numberとする。

	if(flag==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].surface==ON)//この表面粒子が本当に表面粒子足りうるか判定する
			{
				int flag=ON;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(j>=fluid_number)							//壁粒子のsurface_jusge2は壁粒子だけで行う
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
						double ny=Y/dis;
						double nz=Z/dis;
						double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位ﾍﾞｸﾄﾙとの内積
						if(inp<-0.5) flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
					}
				}
				PART[i].surface=flag;
			}
		}
	}
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];

}

///freeon関数 ver.3 粒子数密度のみ再計算。粒子―粒子関係は変化しないと仮定している
void freeon3(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int out)
{
	double d=2;
	if(CON->get_dimention()==3) d=3;
	double R=CON->get_re()*CON->get_distancebp();
	double R2=CON->get_re2()*CON->get_distancebp();
	for(int i=0;i<out;i++)
	{
		double W=0;
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double w=kernel(R,dis);
			//double w=kernel2(R,dis,d);
			W+=w;
		}
		PART[i].PND=W;
	}
	if(CON->get_re()!=CON->get_re2())
	{
		for(int i=0;i<out;i++)
		{
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R2,dis);
				//double w=kernel2(R2,dis,d);
				W+=w;
			}
			PART[i].PND2=W;
		}
	}
	else for(int i=0;i<out;i++) PART[i].PND2=PART[i].PND;
}

//粒子数密度出力
void plot_PND(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t)
{
	double le=CON->get_distancebp();

	ofstream fp("PND.dat");
	if(CON->get_dimention()==2) for(int i=0;i<fluid_number;i++)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PND4[i]<<endl;
	if(CON->get_dimention()==3) for(int i=0;i<fluid_number;i++)	if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PND4[i]<<endl;
	fp.close();///////////*/
}

//粒子数密度プロット関数(毎ステップ出力)
void plot_PND_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t)
{
	double le=CON->get_distancebp();

	char filename[30];
	sprintf_s(filename,"PND%d.dat", t);

	ofstream fp(filename);
	if(CON->get_dimention()==2) for(int i=0;i<fluid_number;i++)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PND4[i]<<endl;
	if(CON->get_dimention()==3) for(int i=0;i<fluid_number;i++)	if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PND4[i]<<endl;
	fp.close();///////////*/
}