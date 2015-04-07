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

///表面張力計算開始関数
void calc_surface_tension(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double dt,int particle_number,double n0,double **potential,int t)
{
	int d=CON->get_dimention();

	if(CON->get_surface_tension()>0)
	{
		if(CON->get_surface_tension()==1)////重みつき最少二乗法だが、ver4と違い、法線発散ではなく、表面曲面(線)を求めその直接微分により曲率をもとめる
		{
			surface_tension1(CON,PART,fluid_number,potential,particle_number);
		}
		else if(CON->get_surface_tension()==2)////粒子間ポテンシャル
		{
			surface_tension2(CON,PART,fluid_number,potential,particle_number);
		}
		else cout<<"指定された表面張力計算方法が見当たりません"<<endl;
		
		//スムージング
		if(CON->get_smooth()==ON) smoothing(CON,PART,particle_number,fluid_number,potential,n0);

		/*/////表面張力表示
		ofstream vec("vector.dat");
		double le=CON->get_distancebp();
		double times=CON->get_times()*le*le;
		if(d==2) for(int i=0;i<fluid_number;i++) vec<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<potential[A_X][i]*times<<" "<<potential[A_Y][i]*times<<endl;
		else if(d==3) for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) vec<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<potential[A_X][i]*times<<"\t"<<potential[A_Z][i]*times<<endl;
		vec.close();
		////////////////*/

		for(int D=0;D<d;D++) for(int i=0;i<fluid_number;i++) PART[i].potential[D]=potential[D][i];

		//表面張力表示
		plot_ST(CON,PART,fluid_number,potential,t);
		
		if(CON->get_ST_interval()>0)
		{
			//if(t==1 || t%CON->get_ST_interval()==0) plot_ST_each(CON,PART,fluid_number,potential,t);
			if(t==1 || (t-1)%CON->get_EM_interval()==0) plot_ST_each(CON,PART,fluid_number,potential,t);
		}//*/

		/*/表面張力AVSファイル出力関数
		if(CON->get_avs_stension_interval()>0)
		{
			if(CON->get_current_step()==1 || CON->get_current_step()%CON->get_avs_stension_interval()==0) plot_avs_stension(CON,PART,fluid_number);
		}//*/

	}
	else for(int D=0;D<d;D++) for(int i=0;i<fluid_number;i++) potential[D][i]=0;//計算しない場合も初期化だけしておく
}

////WLSMによる、表面曲面(線)解法
void surface_tension1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int particle_number)
{   
	cout<<"WLSMによる表面張力計算ver.1---------";

	unsigned int timeA=GetTickCount();
    double le=CON->get_distancebp();		//平均粒子間距離
    double R=CON->get_re3()*le;				//発散用影響半径
    double mass=CON->get_particle_mass();

	int dim=CON->get_dimention();
	int N=0;							//係数行列の元
	int order=2;						//近似曲面のｵｰﾀﾞｰ。 2=二次 3=3次
	double err=(0.01*le)*(0.01*le);		//エラー値

    double *curv=new double [fluid_number];	//各粒子の曲率格納
	double *W=new double [fluid_number];	//各粒子の重み格納

    double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];//内向き法線ベクトル

	double *modify[DIMENTION];
    for(int D=0;D<DIMENTION;D++) modify[D]=new double [fluid_number];//内向き法線ベクトル

	int *particle_J=new int[particle_number];	//particle_J[i]=ONなら、粒子iは曲率計算に考慮される。通常は流体表面粒子のみON.ただしﾉﾝｽﾘｯﾌﾟ考慮の際は壁表面粒子もON
	int *calced=new int[fluid_number];	//calced[i]=ONなら、粒子iは曲率計算が実行された。OFFなら近隣粒子数不足で実行されていない

	//初期化
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			potential[D][i]=0;
			modify[D][i]=0;	
		}
		curv[i]=0;
		W[i]=0;
		if(PART[i].surface==ON) particle_J[i]=ON;//曲率計算に寄与する
		else particle_J[i]=OFF;
		calced[i]=OFF;
	}
	if(CON->get_non_slip()==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL && PART[i].surface==ON) particle_J[i]=ON;//壁表面粒子も曲率計算に寄与する
			else particle_J[i]=OFF;
		}
	}
	else if(CON->get_non_slip()==OFF) for(int i=fluid_number;i<particle_number;i++) particle_J[i]=OFF;//ﾉﾝｽﾘｯﾌﾟでないなら流体表面粒子だけで曲率計算
		
	//係数行列の大きさの決定
	if(dim==2)
	{
		if(order==2) N=3;
		else if(order==3) N=4;
	}
	else if(dim==3)
	{
		if(order==2) N=6;
		else if(order==3) cout<<"order=3??"<<endl;
	}
	////////////////////////////////

	double *matrix=new double [N*N];	//N×Nの係数行列
	double *Bx=new double [N];	//Nの解行列
	double *By=new double [N];	//Nの解行列
	double *Bz=new double [N];	//Nの解行列
    
	//////法線ﾍﾞｸﾄﾙ計算  
	ofstream fq("direct.dat");
    for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
	
		if(PART[i].type==INWALL && PART[i].surface==ON)
		{
			if(CON->get_wall_direct()==OFF) for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
			else if(CON->get_wall_direct()==2)//2DFEM静電霧化のときはこれ
			{
				if(CON->get_dimention()==2)
				{
					direct[A_X][i]/=sqrt(direct[A_X][i]*direct[A_X][i]);//水平成分のみ残す
					direct[A_Y][i]=0;
				}
				else if(CON->get_dimention()==3)
				{
					double r=sqrt(direct[A_X][i]*direct[A_X][i]+direct[A_Y][i]*direct[A_Y][i]);
					if(r!=0)
					{
						direct[A_X][i]/=r;
						direct[A_Y][i]/=r;
					}
					direct[A_Z][i]=0;//水平成分のみ残す
				}
			}
			else if(CON->get_wall_direct()==3)
			{
				if(CON->get_dimention()==2)
				{
					direct[A_Y][i]/=sqrt(direct[A_Y][i]*direct[A_Y][i]);//垂直成分
					direct[A_X][i]=0;
				}
				else if(CON->get_dimention()==3)
				{
					direct[A_X][i]=0;
					direct[A_Y][i]=0;
					direct[A_Z][i]/=sqrt(direct[A_Z][i]*direct[A_Z][i]);//垂直
				}
			}
		}

		if(CON->get_dimention()==2) fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<direct[A_X][i]*0.0001<<" "<<direct[A_Y][i]*0.0001<<endl;
		else if(CON->get_dimention()==3)
		{
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le)
			{
				fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<direct[A_X][i]*0.0001<<" "<<direct[A_Z][i]*0.0001<<endl;
			}
		}
    }
	fq.close();
    ////////////////*/

	if(dim==2 && order==2)//二次元2次式
	{
		//stand_d[i]はA_XかA_Yのどちらか
		//近似曲線はy=ax2+bx+c ここでのxは残差　yは真値
		///係数行列は
		///  ΣΔx4   ΣΔx3  ΣΔx2  a = ΣΔx2*yj  
		///  ΣΔx3   ΣΔx2  ΣΔx   b = ΣΔx*yj 
		///  ΣΔx2   ΣΔx   Σ1     c = Σyj

		//このときyの2階微分をy'',1階微分をy'とおくなら、曲率kはy''/(1+y'^2)^1.5となる。ここで、y''=2a, y'=2ax+bである(ただし粒子iにおいてはx=0となるから、結局 y'=b)
		
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;//初期化
				for(int n=0;n<N;n++){By[n]=0;}
				int num=1;

				double Cx=PART[i].r[A_X]+direct[A_X][i]*le;//点(Cx,Cy)まわりに、周辺粒子を移動させる
				double Cy=PART[i].r[A_Y]+direct[A_Y][i]*le;

				double SIN=-direct[A_X][i];//このように定義したcosθとsinθをもちいて点(Cx,Cy)まわりに回転させる。
				double COS=-direct[A_Y][i];

				double Xi=COS*PART[i].r[A_X]-SIN*PART[i].r[A_Y]+Cx-Cx*COS+Cy*SIN;//点(Cx,Cy)まわりにθ回転させた粒子iの座標
				double Yi=SIN*PART[i].r[A_X]+COS*PART[i].r[A_Y]+Cy-Cx*SIN-Cy*COS;

				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(particle_J[j]==ON)//ONの粒子のみ曲率計算に寄与
					{
						double product=direct[A_X][i]*direct[A_X][j]+direct[A_Y][i]*direct[A_Y][j];//粒子iと粒子jの法線ベクトルの内積

						if(product>0)//内積が負の粒子は考慮しない。特に流体が尖ってre3>disになるとまずいから
						{
							//double Xj=COS*PART[j].r[A_X]-SIN*PART[j].r[A_Y]+Cx-Cx*COS+Cy*SIN;//点(Cx,Cy)まわりにθ回転させた粒子jの座標
							double Yj=SIN*PART[j].r[A_X]+COS*PART[j].r[A_Y]+Cy-Cx*SIN-Cy*COS;

							double dX=PART[j].r[A_X]-PART[i].r[A_X];//実際の座標系での粒子jと粒子iのX座標の差
							double dY=PART[j].r[A_Y]-PART[i].r[A_Y];//実際の座標系での粒子jと粒子iのY座標の差
	
							double X=COS*dX-SIN*dY;//回転させた座標系での両粒子のX座標の差.けっきょくﾍﾞｸﾄﾙの回転であることに気づくこと.(Xj-Xiと計算してもいいけど)
	
							double Y=Yj;
							double dis=sqrt(X*X);//ここでの距離はx方向の距離
							
							double w=1;
							if(dis>le) w=le*le/(dis*dis);
							
							matrix[0]+=X*X*X*X*w;			
							matrix[1]+=X*X*X*w;
							matrix[2]+=X*X*w;		
								
							matrix[5]+=X*w;		
			
							matrix[8]+=w;			
		
							By[0]+=X*X*Y*w;
							By[1]+=X*Y*w;
							By[2]+=Y*w;
		
							num++;
							
						}
					}
				}
				//if(i==0) cout<<By[0]<<" "<<By[1]<<" "<<By[2]<<endl;
				if(num>2)//周辺表面粒子数がこれより少ないとエラー
				{
					calced[i]=ON;		//計算したという印をつけておく
					matrix[3]=matrix[1];		
					matrix[4]=matrix[2];
					matrix[6]=matrix[2];		
					matrix[7]=matrix[5];
		
					matrix[8]+=1;//自分自身
					//By[2]+=PART[i].r[stand_d[i]];//自分自身
					By[2]+=Yi;//自分自身
					//これ以外はX=0になるので加算する必要はない
					
					//行列をガウスの消去法で解く　解はByに格納される
					gauss(matrix,By,N);

					double a=By[0];
					double b=By[1];
					double c=By[2];
					
					double y2=2*a;					//y''
					double y1=b;					//y' 
					double A=pow((1+y1*y1),1.5);	//(1+y'^2)^1.5
						
					curv[i]=-y2/A;//絶対ﾏｲﾅｽをつける
					//cout<<By[0]<<" "<<By[1]<<" "<<By[2]<<endl;

					//修正量計算
					double Yc=c-Yi;//修正量
					if(fabs(Yc)>0.5*le) Yc=0;
					for(int D=0;D<3;D++) modify[D][i]+=-Yc*direct[D][i];
					W[i]+=1;
					
				}
				else
				{
					calced[i]=OFF;		//計算していないという印をつけておく
					cout<<"粒子数不足？"<<endl;
				}
			}
		}
	}
	else if(dim==2 && order==3)//二次元3次式
	{
		//stand_d[i]はA_XかA_Yのどちらか
		//近似曲線はy=ax3+bx2+cx+d ここでのxは残差　yは真値
		///係数行列は
		///  ΣΔx6   ΣΔx5  ΣΔx4  ΣΔx3  a = ΣΔx3*yj  
		///  ΣΔx5   ΣΔx4  ΣΔx3  ΣΔx2  b = ΣΔx2*yj 
		///  ΣΔx4   ΣΔx3  ΣΔx2  ΣΔx   c = ΣΔx*yj
		///  ΣΔx3   ΣΔx2  ΣΔx1  Σ1     d = ΣΔyj

		//このときyの2階微分をy'',1階微分をy'とおくなら、曲率kはy''/(1+y'^2)^1.5となる。ここで、y''=6ax+2b, y'=3ax2+2bx+cである(ただし粒子iにおいてはx=0となるから、結局 y''=2b, y'=cとなる)
		
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;//初期化
				for(int n=0;n<N;n++){By[n]=0;}
				int num=1;

				double Cx=PART[i].r[A_X]+direct[A_X][i]*le;//点(Cx,Cy)まわりに、周辺粒子を移動させる
				double Cy=PART[i].r[A_Y]+direct[A_Y][i]*le;

				double SIN=-direct[A_X][i];//このように定義したcosθとsinθをもちいて点(Cx,Cy)まわりに回転させる。
				double COS=-direct[A_Y][i];

				double Xi=COS*PART[i].r[A_X]-SIN*PART[i].r[A_Y]+Cx-Cx*COS+Cy*SIN;//点(Cx,Cy)まわりにθ回転させた粒子iの座標
				double Yi=SIN*PART[i].r[A_X]+COS*PART[i].r[A_Y]+Cy-Cx*SIN-Cy*COS;

				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(particle_J[j]==ON)//ONの粒子のみ曲率計算に寄与
					{
						//double Xj=COS*PART[j].r[A_X]-SIN*PART[j].r[A_Y]+Cx-Cx*COS+Cy*SIN;//点(Cx,Cy)まわりにθ回転させた粒子jの座標
						double Yj=SIN*PART[j].r[A_X]+COS*PART[j].r[A_Y]+Cy-Cx*SIN-Cy*COS;

						double dX=PART[j].r[A_X]-PART[i].r[A_X];//実際の座標系での粒子jと粒子iのX座標の差
						double dY=PART[j].r[A_Y]-PART[i].r[A_Y];//実際の座標系での粒子jと粒子iのY座標の差

						double X=COS*dX-SIN*dY;//回転させた座標系での両粒子のX座標の差.けっきょくﾍﾞｸﾄﾙの回転であることに気づくこと.(Xj-Xiと計算してもいいけど)

						double Y=Yj;
						double dis=sqrt(X*X);//ここでの距離はx方向の距離
						
						double w=1;
						if(dis>le) w=le*le/(dis*dis);
						
						matrix[0]+=X*X*X*X*X*X*w;			
						matrix[1]+=X*X*X*X*X*w;	
						matrix[2]+=X*X*X*X*w;
						matrix[3]+=X*X*X*w;
							
						matrix[7]+=X*X*w;	

						matrix[11]+=X*w;
		
						matrix[15]+=w;
								
						By[0]+=X*X*X*Y*w;
						By[1]+=X*X*Y*w;
						By[2]+=X*Y*w;
						By[3]+=Y*w;
	
						num++;
					}
				}
				if(num>4)//周辺表面粒子数がこれより少ないとエラー
				{
					calced[i]=ON;		//計算したという印をつけておく

					matrix[4]=matrix[1];
					matrix[5]=matrix[2];
					matrix[6]=matrix[3];

					matrix[8]=matrix[5];
					matrix[9]=matrix[6];
					matrix[10]=matrix[7];

					matrix[12]=matrix[9];
					matrix[13]=matrix[10];
					matrix[14]=matrix[11];

					matrix[15]+=1;//自分自身
					By[3]+=Yi;//自分自身
					//これ以外はX=0になるので加算する必要はない

					//行列をガウスの消去法で解く　解はByに格納される
					gauss(matrix,By,N);

					double a=By[0];
					double b=By[1];
					double c=By[2];
					double d=By[3];

					double y2=2*b;					//y''
					double y1=c;					//y' 
					double A=pow((1+y1*y1),1.5);	//(1+y'^2)^1.5
						
					curv[i]=-y2/A;//絶対ﾏｲﾅｽをつける

					//修正量計算
					double Yc=d-Yi;//修正量
					for(int D=0;D<3;D++) modify[D][i]+=-Yc*direct[D][i];
					W[i]+=1;
				}
				else cout<<"error in surface_tension ver.6 粒子数<=4"<<endl;
			}
		}
			
	}
	else if(dim==3 && order==2)//3次元1次式
	{
		//近似曲線はz=ax2+by2+cxy+dx+ey+f ここでのx,yは残差　zは真値
		///係数行列は
		///  ΣΔx4      ΣΔx2Δy2  ΣΔx3Δy  ΣΔx3     ΣΔx2Δy  ΣΔx2    a = ΣΔx2*zj  
		///  ΣΔx2Δy2  ΣΔy4      ΣΔxΔy3  ΣΔxΔy2  ΣΔy3     ΣΔy2    b = ΣΔy2*zj
		///  ΣΔx3Δy   ΣΔxΔy3   ΣΔx2Δy2 ΣΔx2Δy  ΣΔxΔy2  ΣΔxΔy  c = ΣΔxΔy*zj
		///  ΣΔx3      ΣΔxΔy2   ΣΔx2Δy  ΣΔx2     ΣΔxΔy   ΣΔx     d = ΣΔx*zj
		///  ΣΔx2Δy   ΣΔy3      ΣΔxΔy2  ΣΔxΔy   ΣΔy2     ΣΔy     e = ΣΔy*zj
		///  ΣΔx2      ΣΔy2      ΣΔxΔy   ΣΔx      ΣΔy      Σ1       f = Σzj

		//このときzの2階微分をZxx,Zyy,1階微分をZx,Zyとおくなら、曲率kは-div{∇Z/(1+|∇Z|^2)^1.5} = -(Zxx+Zyy+Zxx*Zy^2+Zyy*Zx^2)/(1+Zx^2+Zy^2)^1.5となる。ここで、粒子iにおいてはx=y=0となるから、結局 Zxx=2a, Zyy=2b, Zx=d, Zy=eとなる)
		
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;//初期化
				for(int n=0;n<N;n++){Bz[n]=0;}
				int num=1;
				
				double Cx=PART[i].r[A_X]+direct[A_X][i]*le;//点(Cx,Cy,Cz)まわりに、周辺粒子を移動させる
				double Cy=PART[i].r[A_Y]+direct[A_Y][i]*le;
				double Cz=PART[i].r[A_Z]+direct[A_Z][i]*le;

				double normal=sqrt(direct[A_Y][i]*direct[A_Y][i]+direct[A_Z][i]*direct[A_Z][i]);

				double SIN1=-direct[A_Y][i]/normal;//このように定義したcosθとsinθをもちいて点(Cx,Cy,Cz)を通るX軸まわりにまずは回転させる。
				double COS1=-direct[A_Z][i]/normal;

				double X1=PART[i].r[A_X];
				double Y1=COS1*PART[i].r[A_Y]-SIN1*PART[i].r[A_Z]+Cy-COS1*Cy+SIN1*Cz;//点(Cx,Cy,Cz)を通るX軸まわりに回転させた後の座標
				double Z1=SIN1*PART[i].r[A_Y]+COS1*PART[i].r[A_Z]+Cz-SIN1*Cy-COS1*Cz;

				normal=sqrt(direct[A_X][i]*direct[A_X][i]+(SIN1*direct[A_Y][i]+COS1*direct[A_Z][i])*(SIN1*direct[A_Y][i]+COS1*direct[A_Z][i]));

				double SIN2=direct[A_X][i]/normal;//このように定義したcosθとsinθをもちいて,次は点(Cx,Cy,Cz)を通るY軸まわりに回転させる。
				double COS2=-1*(SIN1*direct[A_Y][i]+COS1*direct[A_Z][i])/normal;//これは回転後のnzに相当

				double Xi=COS2*X1+SIN2*Z1+Cx-COS2*Cx-SIN2*Cz;//点(Cx,Cy,Cz)を通るY軸まわりに回転させた後の座標
				double Yi=Y1;
				double Zi=-SIN2*X1+COS2*Z1+Cz+SIN2*Cx-COS2*Cz;

				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(particle_J[j]==ON)//ONの粒子のみ曲率計算に寄与
					{
						double product=direct[A_X][i]*direct[A_X][j]+direct[A_Y][i]*direct[A_Y][j]+direct[A_Z][i]*direct[A_Z][j];//粒子iと粒子jの法線ベクトルの内積
						//if(product>0.5)//内積が負の粒子は考慮しない。特に流体が尖ってre3>disになるとまずいから
						if(product>0)//内積が負の粒子は考慮しない。特に流体が尖ってre3>disになるとまずいから
						{
							
							double Xj1=PART[j].r[A_X];
							double Yj1=COS1*PART[j].r[A_Y]-SIN1*PART[j].r[A_Z]+Cy-COS1*Cy+SIN1*Cz;//点(Cx,Cy,Cz)を通るX軸まわりに回転させた後の座標
							double Zj1=SIN1*PART[j].r[A_Y]+COS1*PART[j].r[A_Z]+Cz-SIN1*Cy-COS1*Cz;

							double Xj=COS2*Xj1+SIN2*Zj1+Cx-COS2*Cx-SIN2*Cz;//点(Cx,Cy,Cz)を通るY軸まわりに回転させた後の座標
							double Yj=Yj1;
							double Zj=-SIN2*Xj1+COS2*Zj1+Cz+SIN2*Cx-COS2*Cz;

							double X=Xj-Xi;
							double Y=Yj-Yi;
							double Z=Zj;
							double dis=sqrt(X*X+Y*Y);//ここでの距離はXY平面上の距離
						
							double w=1;
							if(dis>le) w=le*le/(dis*dis);
						
							matrix[0]+=X*X*X*X*w;			
							matrix[1]+=X*X*Y*Y*w;	
							matrix[2]+=X*X*X*Y*w;
							matrix[3]+=X*X*X*w;
							matrix[4]+=X*X*Y*w;
							matrix[5]+=X*X*w;
							
							matrix[7]+=Y*Y*Y*Y*w;
							matrix[8]+=X*Y*Y*Y*w;
							matrix[9]+=X*Y*Y*w;
							matrix[10]+=Y*Y*Y*w;
							matrix[11]+=Y*Y*w;

							matrix[17]+=X*Y*w;

							matrix[23]+=X*w;

							matrix[29]+=Y*w;

							matrix[35]+=w;
								
							Bz[0]+=X*X*Z*w;
							Bz[1]+=Y*Y*Z*w;
							Bz[2]+=X*Y*Z*w;
							Bz[3]+=X*Z*w;
							Bz[4]+=Y*Z*w;
							Bz[5]+=Z*w;

							num++;
						}
					}
				}
				if(num>4)//周辺表面粒子数がこれより少ないとエラー
				{
					calced[i]=ON;		//計算したという印をつけておく

					matrix[6]=matrix[1];

					matrix[12]=matrix[2];
					matrix[13]=matrix[8];
					matrix[14]=matrix[1];
					matrix[15]=matrix[4];
					matrix[16]=matrix[9];

					matrix[18]=matrix[3];
					matrix[19]=matrix[9];
					matrix[20]=matrix[15];
					matrix[21]=matrix[5];
					matrix[22]=matrix[17];

					matrix[24]=matrix[4];
					matrix[25]=matrix[10];
					matrix[26]=matrix[16];
					matrix[27]=matrix[22];
					matrix[28]=matrix[11];

					matrix[30]=matrix[5];
					matrix[31]=matrix[11];
					matrix[32]=matrix[17];
					matrix[33]=matrix[23];
					matrix[34]=matrix[29];

					matrix[35]+=1;//自分自身
					Bz[5]+=Zi;
					//これ以外はX=0,Y=0になるので加算する必要はない

					//行列をガウスの消去法で解く　解はBzに格納される
					gauss(matrix,Bz,N);

					double a=Bz[0];
					double b=Bz[1];
					double c=Bz[2];
					double d=Bz[3];
					double e=Bz[4];
					double f=Bz[5];
					
					double Zxx=2*a;
					double Zyy=2*b;
					double Zx=d;
					double Zy=e;
					double A=pow((1+Zx*Zx+Zy*Zy),1.5);
					double A2=Zxx*(1+Zy*Zy)+Zyy*(1+Zx*Zx);
					
					curv[i]=-A2/A;
					

					//修正量計算
					double Zc=f-Zi;//修正量
					
					if(2*curv[i]!=curv[i]+curv[i])//エラーなら
					{
						curv[i]=0;//すべてをゼロに
						Zc=0;
					}
					//if(Zc>0.1*le || Zc<-0.1*le) Zc=0;//修正量に上限、下限を設ける
					if(Zc>0.1*le) Zc=0.1*le;
					if(Zc<-0.1*le) Zc=-0.1*le;

					for(int D=0;D<3;D++) modify[D][i]+=-Zc*direct[D][i];
					W[i]+=1;
				}
				else curv[i]=0;
			}
		}
	}

	//位置修正
	if(CON->get_suf_modify()==ON)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(W[i]>0) for(int D=0;D<DIMENTION;D++) PART[i].r[D]+=modify[D][i]/W[i];
		}
	}

	int count;
	count=0;
	if(CON->get_interpolate_curv()==ON)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				if(calced[i]==OFF)//表面粒子なのに曲率が計算されていなかったら
				{
					double W=0;
					for(int k=0;k<PART[i].N3;k++)
					{
						int j=PART[i].NEI3[k];
						if(PART[j].type==FLUID)
						{
							if(calced[j]==ON)//内部粒子のcalced[]はOFFなので、自動的にここでは表面粒子のみ対象に絞られる
							{
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
								double w=(1-dis/R)*(1-dis/R);
								W+=w;
								curv[i]+=curv[j]*w;
							}
						}
					}
					if(W>0) curv[i]/=W;
					else 
					{
						//cout<<"in interpolate_curv() 周辺粒子が不在"<<endl;
						count++;
					}
				}
			}
		}
	}

	if(count>0) cout<<"in interpolate_curv() 周辺が不在の粒子存在 n="<<count<<endl;
	////////////////////////////////////curvのスムージング
    if(CON->get_smth_sumP()!=OFF)
    {
        double *newP=new double [fluid_number];
        for(int n=0;n<CON->get_smth_sumP();n++)
        {
            for(int i=0;i<fluid_number;i++) 
            {  
				if(PART[i].surface==OFF) newP[i]=0;
				else
				{  
					newP[i]=curv[i];
					int num=1;
					double W=1;
					for(int k=0;k<PART[i].N3;k++)
					{
						int j=PART[i].NEI3[k];
						if(PART[j].type==FLUID &&PART[j].surface==ON ) 
						{ 	
							double dif=fabs(curv[i]-curv[j]);//曲率の差
						//	if(dif<400)//曲率の差があまりにも激しいものは平滑化に考慮しない
							{
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
								double w=(1-dis/R)*(1-dis/R);
								W+=w;
								num++;
								newP[i]+=curv[j]*w;
							}
						}
					}
					newP[i]/=W;
				}
            } 
            for(int i=0;i<fluid_number;i++) curv[i]=newP[i];
		}
		delete [] newP;
    }
    ////////////////////*/

	///曲率をﾌｧｲﾙ出力
	ofstream fp2("curv.dat");///曲率出力
	if(dim==2) for(int i=0;i<fluid_number;i++) if(PART[i].surface==ON) fp2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<curv[i]<<endl;
	if(dim==3) for(int i=0;i<fluid_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fp2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<curv[i]<<endl;
	fp2.close();///////////*/

	double *sigma=new double[fluid_number];
	for(int i=0;i<fluid_number;i++) sigma[i]=0;//初期化
	///potential[D][i]に代入
	if(CON->get_SFT()==OFF)
	{
		double L=1.5*le;//表面の厚み
		double density;
		for(int i=0;i<fluid_number;i++)
		{   
			sigma[i]=CON->get_sigma();
		    if(PART[i].surface==ON)
			{
				if(PART[i].materialID==1) density=CON->get_density();
				else if(PART[i].materialID==2) density=CON->get_density2();
				for(int D=0;D<DIMENTION;D++) potential[D][i]+=CON->get_sigma()*curv[i]/(density*L)*direct[D][i];
			}
		}
	}
	else 
	{
		for(int i=0;i<fluid_number;i++)
		{   
			sigma[i]=CON->get_sigma();//ひとまず元の設定値で初期化
		}
		//cout<<"現在は表面張力の温度依存性は計算不可"<<endl;
		double L=1.5*le;//表面の厚み
		double *density=new double[fluid_number];

		calc_physical_property(CON,PART,fluid_number,density,particle_number,1);//密度の温度依存
		calc_physical_property(CON,PART,fluid_number,sigma,particle_number,3);//表面張力係数の温度依存

		for(int i=0;i<fluid_number;i++)
		{   
		    if(PART[i].surface==ON)
			{
				for(int D=0;D<DIMENTION;D++) potential[D][i]+=sigma[i]*curv[i]/(density[i]*L)*direct[D][i];
			}
		}
		
		delete [] density;
		
	}
	
	
	if(CON->get_dir_for_P()==1 ||CON->get_dir_for_P()==3 )//圧力の計算時に、表面張力をディリクレ値として使用する
	{
		//ofstream dir("surface_tension_P.dat");
		for(int i=0;i<fluid_number;i++) PART[i].dir_Pst=sigma[i]*curv[i];			//ディリクレ値をセット
		//for(int i=0;i<fluid_number;i++) PART[i].dir_Pst=sigma[i]*curv[i];			//ディリクレ値をセット
		//for(int i=0;i<fluid_number;i++)  dir<<CON->get_sigma()*curv[i]<<endl;
		//dir.close();
	}

    for(int D=0;D<DIMENTION;D++) delete [] direct[D];
    delete [] curv;
	delete [] W;
	delete [] matrix;
	delete [] Bx;
	delete [] By;
	delete [] Bz;

	for(int D=0;D<DIMENTION;D++) delete [] modify[D];
	delete [] particle_J;
	delete [] calced;
	delete []sigma;

	cout<<"ok  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//粒子間ポテンシャル
void surface_tension2(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int particle_number)
{
	cout<<"粒子間引力にもとづく表面張力計算開始--------";

	unsigned int timeA=GetTickCount();		//計算開始時刻
    double le=CON->get_distancebp();		//平均粒子間距離
    double r=CON->get_re3();				//表面張力用影響半径
	double mass=CON->get_particle_mass();	//粒子の質量[kg]
	double Cst=CON->get_Cst();				//一応最初に計算されたポテンシャル係数
    double C;								//ポテンシャル係数
	double C2=Cst*CON->get_C_times();		//ポテンシャル係数の補正(これをしないと粒子間力が強くなりすぎる)
	double wall_C=CON->get_wall_C();		//壁との親和力係数

	if(CON->get_freeon()==1)
	{
		if(CON->get_SFT()==OFF)//表面張力の温度依存性なし
		{
		    for(int i=0;i<fluid_number;i++)
			{   
				for(int D=0;D<DIMENTION;D++) potential[D][i]=0;
				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
				           
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
						    
					C=C2;
			  
					double p_slope;//p(r)の勾配
					//  else   C=Cst*(1+cos(PI/2))/2; //論文通り
	
					if(PART[j].type!=FLUID) C=C2*wall_C;
					p_slope=C*(dis-le)*(dis-le*r);
					p_slope/=mass;	//仮の位置・速度の計算でmassで割るのでここでは[N]で求める
					potential[A_X][i]-=p_slope*X/dis;
					potential[A_Y][i]-=p_slope*Y/dis;	
					potential[A_Z][i]-=p_slope*Z/dis;		
				}
			}
		}
    }//potential[D][i]が求まった*/
	
	if(CON->get_freeon()!=1)
	{
		if(CON->get_SFT()==OFF)//表面張力の温度依存性なし
		{
			int *check=new int[fluid_number];//0なら未処理　1は処理済み
			for(int i=0;i<fluid_number;i++)
			{
				check[i]=0;//初期化
				for(int D=0;D<DIMENTION;D++) potential[D][i]=0;
			}
		    for(int i=0;i<fluid_number;i++)
			{   
				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
				    if(j<fluid_number)//粒子jが流体粒子なら
					{
						if(check[j]==0)//まだ処理が終わってないなら
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
							    
							C=C2;
					
							double p_slope;//p(r)の勾配
							
							p_slope=C*(dis-le)*(dis-le*r);
							p_slope/=mass;
							potential[A_X][i]-=p_slope*X/dis;//粒子jが粒子iにおよぼす力
							potential[A_Y][i]-=p_slope*Y/dis;
							potential[A_Z][i]-=p_slope*Z/dis;

							potential[A_X][j]+=p_slope*X/dis;//粒子iが粒子jにおよぼす力
							potential[A_Y][j]+=p_slope*Y/dis;
							potential[A_Z][j]+=p_slope*Z/dis;
						}
					}
					else
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
							    
						C=C2*wall_C;
					
						double p_slope;//p(r)の勾配
						
						p_slope=C*(dis-le)*(dis-le*r);
						p_slope/=mass;
						potential[A_X][i]-=p_slope*X/dis;//粒子jが粒子iにおよぼす力
						potential[A_Y][i]-=p_slope*Y/dis;
						potential[A_Z][i]-=p_slope*Z/dis;
					}
				}
				check[i]=1;
			}
			delete [] check;
	    }//potential[D][i]が求まった
	}////*/

    /////表面張力の温度依存性を考える場合
    if(CON->get_SFT()==1)
    {    
        double hs0=mass*CON->get_Cp()*CON->get_MP();//融解開始点のエンタルピー
		double hs1=hs0+CON->get_latent_H()*mass;    //融解終了点のエンタルピー
		double *C1=new double[fluid_number];//粒子iの表面張力係数
		///エンタルピーから各粒子の温度T[i]を求める
		for(int i=0;i<fluid_number;i++)
		{
		    //粒子はすべて液体
		    double T=CON->get_MP()+(PART[i].h-hs1)/mass/CON->get_Cp();
		    C1[i]=-0.0003*T*T+1.0343*T+951.65;//この値はmN/m
		    C1[i]/=1000;//単位をN/mになおす。
		    C1[i]/=CON->get_sigma();//標準値との比率をだす。
		    //cout<<C1[i]<<endl;
		}
		  
		for(int i=0;i<fluid_number;i++)
		{
		    for(int D=0;D<DIMENTION;D++) potential[D][i]=0;
		    for(int k=0;k<PART[i].N3;k++)
		    {
		        int j=PART[i].NEI3[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(dis<r*le && i!=j)
				{   
					double p_slope;//p(r)の勾配
					if(PART[j].type==FRFLUID ||PART[j].type==BOFLUID) C=Cst*C1[i]*C1[j]/(C1[i]+C1[j])*2;
					else if(PART[j].type!=BOFLUID && PART[j].type!=FRFLUID) C=1.0*Cst*C1[i];
					p_slope=C*(dis-le)*(dis-le*r);
					p_slope/=mass;
					potential[A_X][i]-=p_slope*X/dis;
					potential[A_Y][i]-=p_slope*Y/dis;
					potential[A_Z][i]-=p_slope*Z/dis;
				}
			} 
		}//potential[D][i]が求まった */
		delete [] C1;
    }
	/////////////////*/

	
	cout<<"ok  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}


//表面張力スムージング関数
void smoothing(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double *potential[DIMENTION],double n0)
{
	int d=CON->get_dimention();
	double le=CON->get_distancebp();
    double r=CON->get_re()*le;  //ｽﾑｰｼﾞﾝｸﾞの範囲

	///potential[D][i]のスムージング////////
	double *newtension[DIMENTION];
	for(int D=0;D<DIMENTION;D++) newtension[D]=new double [fluid_number];
	
	for(int n=0;n<CON->get_smn();n++)
    {
        for(int i=0;i<fluid_number;i++) 
        {  
            for(int D=0;D<3;D++) newtension[D][i]=potential[D][i];
			int num=1; //自分自身をｶｳﾝﾄするから1
			for(int k=0;k<PART[i].N;k++)
			{       
				int j=PART[i].NEI[k];
	        
				if(PART[j].type==FLUID)
				{
					num++;
					for(int D=0;D<3;D++) newtension[D][i]+=potential[D][j];
				}
			}
			for(int D=0;D<3;D++) newtension[D][i]/=num;
        } 
        for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) potential[D][i]=newtension[D][i];
    }
    /////////////////////////////////////*/
	
	for(int D=0;D<DIMENTION;D++) delete [] newtension[D];
}

double calc_Cst(mpsconfig *CON)
{
	double p;//ポテンシャル
	double sumP=0;//Σp(r)
	double le=CON->get_distancebp();	//初期粒子間距離
	double re=CON->get_re3()*le;		//影響半径
	int calc_type=CON->get_model_set_way();	//モデルのセットタイプ(正方 or MD)

	if(CON->get_dimention()==2)
	{
		for(int k=1;k<10;k++)
		{       
			int count=0;
			for(int i=-10;i<=10;i++)
			{
				for(int j=k;j<k+10;j++)
				{
					double dis=sqrt((double)(i*i+j*j));
					dis*=le;

					if(dis<re)
					{
						p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
						sumP+=p;
						count++;
					}
				}
			}
		}
		return 2*le*CON->get_sigma()/sumP;
	}

    else if(CON->get_dimention()==3) 
	{
		//////正方格子配置の場合
		if(calc_type==0)
		{
			//////１粒子との総和
			if(CON->get_C_type()==0)
			{
				for(int k=1;k<10;k++)
				{       
					int count=0;
					for(int i=-10;i<=10;i++)
					{
						for(int j=-10;j<=10;j++)
						{
							double dis=sqrt((double)(i*i+j*j+k*k));
							dis*=le;

							if(dis<re)
							{
								p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
								sumP+=p;
								count++;
							}
						}
					}
				}
			}
			//////近藤らの方法
			else if(CON->get_C_type()==1)
			{
				for(int k0=(-1);k0>-10;k0--)
				{
					for(int k=1;k<10;k++)
					{       
						int count=0;
						for(int i=-10;i<=10;i++)
						{
							for(int j=-10;j<=10;j++)
							{
								double dis=sqrt((double)(i*i+j*j+(k-k0)*(k-k0)));
								dis*=le;

								if(dis<re)
								{
									p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
									sumP+=p;
									count++;
								}
							}
						}
					}
				}
			}
			return 2*le*le*CON->get_sigma()/sumP;
		}

		//////MDによる配置の場合
		else if(calc_type==1)	
		{
			//必要変数の宣言
			double xa,ya,za,xb,yb,zb;
			double dx=1;
			double dy=sqrt(3.0)/2;
			double dz=sqrt(6.0)/3;
			int size=20;
			int countA=0;
			int countB=0;

			//////１粒子との総和
			if(CON->get_C_type()==0)
			{
				for(int k=1;k<size;k++)
				{       
					int count=0;
					for(int i=-size;i<=size;i++)
					{
						for(int j=-size;j<=size;j++)
						{
							xa=i*dx;
							ya=j*dy;
							za=k*dz;
							if(j%2==1)	xa=i*dx+(dx/2);
							if(k%2==1)	ya=j*dy+(dy*2/3);

							double dis=sqrt((double)(xa*xa+ya*ya+za*za));
							dis*=le;

							if(dis<re)
							{
								p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
								sumP+=p;
								count++;
							}
						}
					}
				}
			}//*/

			//////近藤らの方法
			if(CON->get_C_type()==1)
			{
				for(int k0=(-1);k0>-10;k0--)
				{
					for(int k=0;k<size;k++)
					{       
						int count=0;
						for(int i=-size;i<=size;i++)
						{
							for(int j=-size;j<=size;j++)
							{
								xa=i*dx;
								ya=j*dy;
								za=k*dz;
								zb=k0*dz;
								if(j%2==1)	xa=i*dx+(dx/2);
								if(k%2==1)	ya=j*dy+(dy*2/3);

								double dis=sqrt((double)(xa*xa+ya*ya+(za-zb)*(za-zb)));
								dis*=le;

								if(dis<re)
								{
									p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
									sumP+=p;
									count++;
								}
							}
						}
					}
				}
			}//*/

			//////二分割した影響半径内の総当たり
			if(CON->get_C_type()==2)
			{
				for(int ka=0;ka<=size;ka++)
				{
					for(int ja=-size;ja<=size;ja++)
					{
						for(int ia=-size;ia<=size;ia++)
						{
							xa=ia*dx;
							ya=ja*dy;
							za=((double)ka+0.5)*dz;
							if(ja%2==1)	xa=ia*dx+(dx/2);
							if(ka%2==1)	ya=ja*dy+(dy*2/3);
							double disA=sqrt((double)(xa*xa+ya*ya+za*za));
							disA*=le;

							if(disA<re)
							{
								countA++;
								for(int kb=0;kb<=size;kb++)
								{
									for(int jb=-size;jb<=size;jb++)
									{
										for(int ib=-size;ib<=size;ib++)
										{
											xb=ib*dx;
											yb=jb*dy;
											zb=((double)kb+0.5)*(-dz);
											if(jb%2==1)	xb=ib*dx+(dx/2);
											if(kb%2==1)	yb=jb*dy+(dy*2/3);
											double disB=sqrt((double)(xb*xb+yb*yb+zb*zb));
											disB*=le;
											
											if(disB<re)
											{
												countB++;
												double dis=sqrt((double)((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya)+(zb-za)*(zb-za)));
												dis*=le;
												
												if(dis<re)
												{
													p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3.0;
													sumP+=p;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}//*/
			return 2*le*le*sqrt(3.0)/4.0*CON->get_sigma()/sumP;
			//return 2*le*le*CON->get_sigma()/sumP;

			cout << "countA=" << countA << endl;
			cout << "countB=" << countB << endl;
		}
	}
	else return 0;

	return 0;	//警告「値を返さないコントロールパスがあります」を消すため
}

//表面張力出力（最新ステップ）
void plot_ST(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int t)
{
	int d=CON->get_dimention();
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le;
	//double times=CON->get_times()*le*le;
	double times=CON->get_times()*CON->get_density();

	ofstream st("ST.dat");
	

	if(d==2) for(int i=0;i<fluid_number;i++) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<potential[A_X][i]*times<<" "<<potential[A_Y][i]*times<<endl;
	else if(d==3) for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<potential[A_X][i]*times<<"\t"<<potential[A_Z][i]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
	}

	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	//if(CON->get_legend_F()>0) st<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力

	st.close();
}

//表面張力プロット関数(毎ステップ出力)
void plot_ST_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int t)
{
	int d=CON->get_dimention();
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le;
	//double times=CON->get_times()*le*le;
	double times=CON->get_times()*CON->get_density();

	char filename[20];
	sprintf_s(filename,"ST%d.dat", t);
	ofstream st(filename);

	if(d==2) for(int i=0;i<fluid_number;i++) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<potential[A_X][i]*times<<" "<<potential[A_Y][i]*times<<endl;
	
	else if(d==3) for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<potential[A_X][i]*times<<"\t"<<potential[A_Z][i]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
	}
	
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	//if(CON->get_legend_F()>0) st<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力

	st.close();
}