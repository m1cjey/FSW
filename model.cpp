#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納
#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい


void writedata(ofstream &fp, ofstream &gnu,int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM)
{	
	//ファイル出力
	fp<<id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<P<<"\t"<<h<<"\t"<<val<<"\t"<<type<<"\t"<<materialID<<"\t"<<surface<<"\t"<<toBEM<<"\t"<<endl;
	
	//////gnuplot用に出力
	gnu<<x<<"\t"<<y<<"\t"<<z<<endl;
}


void set_initial_placement(mpsconfig *CON,int *particle_number)
{
	ofstream fp("initial_input.dat");//粒子データ格納ﾌｧｲﾙ
	ofstream gnu("plot.dat");		//gnuplot用初期配置出力ﾌｧｲﾙ
	
    double ii=0;
	double jj=0;
	double kk=0;
	int number=0;	//粒子数
	int model=CON->get_model_number();
	double mass=CON->get_particle_mass();//粒子の質量
	

	//壁重み関数計算用モデル
	if(model==0)
	{
		double T=293;//温度
		double h=T*mass*CON->get_Cp();//エンタルピー
		int materialID=1;
		double val=0;
		if(CON->get_dimention()==2)
		{
			//if(CON->get_wall_poly()==0)
			{
				//
				//////各riw(壁境界と粒子の距離)での重みを調べるための点。
				ii = 0;
				jj = 0;
				writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
				number++;


				//均質な壁
				for(int i=-5;i<=5;i++)
				{
					for(int j=-5;j<=0;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
			}
			
		}
	}
	//壁付長方形液相
	if(model==1)
	{
		double T=293;//温度
		double h=T*mass*CON->get_Cp();//エンタルピー
		int materialID=1;
		double val=0;
		if(CON->get_dimention()==2)
		{
			//////流体
			for(int i=-CON->get_fluidwidth()/2;i<=CON->get_fluidwidth()/2;i++)
			{
				
				for(int j=-CON->get_fluidwidth();j<=CON->get_fluidwidth();j++)
				{
					ii = i*CON->get_distancebp();
					jj = j*CON->get_distancebp();
					writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
					number++;
				}
			}

			if(CON->get_wall_poly()==0)
			{
				//////内壁
				//左
				for(int i=-CON->get_fluidwidth()/2-1;i<=-CON->get_fluidwidth()/2-1;i++)
				{
					for(int j=-CON->get_fluidwidth();j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//右
				for(int i=CON->get_fluidwidth()/2+1;i<=CON->get_fluidwidth()/2+1;i++)
				{
					for(int j=-CON->get_fluidwidth();j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//下
				for(int i=-CON->get_fluidwidth()/2-1;i<=CON->get_fluidwidth()/2+1;i++)
				{
					for(int j=-CON->get_fluidwidth()-1;j<=-CON->get_fluidwidth()-1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//////外壁
				//左
				for(int i=-CON->get_fluidwidth()/2-5;i<=-CON->get_fluidwidth()/2-2;i++)
				{
					for(int j=-CON->get_fluidwidth()-5;j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//右
				for(int i=CON->get_fluidwidth()/2+2;i<=CON->get_fluidwidth()/2+5;i++)
				{
					for(int j=-CON->get_fluidwidth()-5;j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//下
				for(int i=-CON->get_fluidwidth()/2-1;i<=CON->get_fluidwidth()/2+1;i++)
				{
					for(int j=-CON->get_fluidwidth()-5;j<=-CON->get_fluidwidth()-2;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
			}
			
		}
	}


	//立方体液層
	if(model==2)
	{
		double T=293;//温度
		double h=T*mass*CON->get_Cp();//エンタルピー
		int materialID=1;
		double val=0;
		if(CON->get_dimention()==3)
		{
			for(int i=-CON->get_fluidwidth()/2;i<=CON->get_fluidwidth()/2;i++)
			{
				for(int j=-CON->get_fluidwidth()/2;j<=CON->get_fluidwidth()/2;j++)
				{
					for(int k=-CON->get_fluidwidth()/2;k<=CON->get_fluidwidth()/2;k++)		
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						kk = k*CON->get_distancebp();	
						writedata( fp, gnu, number,  ii,  jj, kk+0.19,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
			}
		}
	}
	//////モデル３　円形水滴
	if(model==3)
	{
		double T=293;//温度
		double h=T*mass*CON->get_Cp();//エンタルピー
		int materialID=1;
		double val=0;
		//流体粒子書き込み
		for(int i=-CON->get_fluidwidth()/2;i<=CON->get_fluidwidth()/2;i++)
		{
			for(int j=-CON->get_fluidwidth()/2;j<=CON->get_fluidwidth()/2;j++)
			{
				if(CON->get_dimention()==2)
				{
					if(sqrt((double)(i*i+j*j))<CON->get_fluidwidth()/2)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();	
						writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
				else if(CON->get_dimention()==3)
				{
					for(int k=-CON->get_fluidwidth()/2;k<=CON->get_fluidwidth()/2;k++)
					{
			    		if(sqrt((double)(i*i+j*j+k*k))<CON->get_fluidwidth()/2)
						{
							ii = i*CON->get_distancebp();
							jj = j*CON->get_distancebp();
							kk = k*CON->get_distancebp();	
							writedata( fp, gnu, number,  ii,  jj, kk,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
							number++;
						}
					}
				}
			}
		}
	}
	//////*/
	
	/////////モデル1 静電霧化
	///推奨　distancebp=0.000025 dt=0.000001 or0.0000005  fluidwidth=20
	if(model==14)
	{
		if(CON->get_dimention()==2)
		{
		    double le=CON->get_distancebp();
		    double val=0;
			int materialID=1;
		    int L=10;//解析領域
		    int W=CON->get_fluidwidth();
		    ///水滴
		    for(int i=-W;i<=W;i++)
		    {
		        for(int j=0;j<=W;j++)
				{
					if(i*i+j*j<W*W)
					{
						ii = i*le;
			    		jj = j*le;
			    		double T=293;//温度
			    		double h=T*mass*CON->get_Cp();//エンタルピー
						writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
		    }
		    
		    ////円柱電極
		    for(int i=-W*2;i<=W*2;i++)
		    {
		        for(int j=-1;j>=-W/3;j--)
				{
		            ii = i*le;
					jj = j*le;
					double T=293;//温度
					double h=T*mass*CON->get_Cp();//エンタルピー
					if(i*i<=W*W)
					{
						if(j==-1 ||i==W || i==-(W))
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
							number++;
						}
						else
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
							number++;
						}
					}
			  
				}
			}
		}
		/*else if(CON->get_dimention()==3)
		{
		    double le=CON->get_distancebp();
			double mass=CON->get_particle_mass();
		    int W=CON->get_fluidwidth();
		    ///水滴
			 for(int i=-W;i<=W;i++)
		    {
		        for(int j=-W;j<=W;j++)
					{
					for(int k=0;k<=W;k++)
					{
						if(i*i+j*j+k*k<W*W)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
		    				double T=293;//温度
				 			double h=T*mass*CON->get_Cp();//エンタルピー
							writedata(fp, number, ii, jj,kk, FRFLUID,gnu,0,0,0,0,h,0);
							number++;
						}
					}
				}
			}

			///上球
			for(int i=-W;i<=W;i++)
		    {
		        for(int j=-W;j<=W;j++)
				{
					for(int k=0;k<=W;k++)
					{
						if(i*i+j*j+k*k<3)
						{
							ii = i*le;
			    			jj = j*le;
							kk = k*le+1.2*W*le;
			    			double T=293;//温度
			    			double h=T*mass*CON->get_Cp();//エンタルピー
							//writedata(fp, number, ii, jj,kk, FRFLUID,gnu,0,0,0,0,h,0);
							//number++;
						}
					}
				}
		    }//////
	    
		    ////円柱電極(INWALL)
		    for(int i=-2*W;i<=2*W;i++)
		    {
		        for(int j=-2*W;j<=2*W;j++)
				{
		            for(int k=-1;k>=-5;k--)
					{
		                ii = i*le;
						jj = j*le;
						kk = k*le;
						double T=293;//温度
						double h=T*mass*CON->get_Cp();//エンタルピー
						int r=i*i+j*j;
						//if(r<=(W+1)*(W+1))//1層厚い電極。水滴が横にこぼれるのを防ぐ？
						if(r<=W*W)
						{
							if(k>=-2 || r>=(W-1)*(W-1))
							{
								writedata(fp, number, ii, jj,kk, INWALL,gnu,0,0,0,0,h,0);
								number++;
							}
						}
					}	  
				}
			}
			////円柱電極(OUTWALL)
			for(int i=-2*W;i<=2*W;i++)
			{
			    for(int j=-2*W;j<=2*W;j++)
				{
		         for(int k=-1;k>=-5;k--)
					{
			            ii = i*le;
						jj = j*le;
						kk = k*le;
						double T=293;//温度
						double h=T*mass*CON->get_Cp();//エンタルピー
						int r=i*i+j*j;
						if(r<(W-1)*(W-1))
						{
							if(k<-2)
							{
								writedata(fp, number, ii, jj,kk, OUTWALL,gnu,0,0,0,0,h,0);
								number++;
							}
						}
					}
			  
				}
			}
		}*/
	}
	///////////////////*/

	////モデル19　FSW
	//半径9mm、深さ6mmの液相に、半径2.5mm、深さ4mmのﾌﾟﾛｰﾌﾞ挿入。ショルダー半径は6mm
	if(model==19)
	{
		if(CON->get_dimention()==3)
		{
			double le=CON->get_distancebp();
			double mass=CON->get_particle_mass();
			int W=CON->get_fluidwidth();	//流体のＸ幅
			int YW=2*W;						//流体のＹ幅
			int fH=W*2/3;					//流体高さ
			int WH=W*2/3;					//壁の高さ
			int shoR=W*2/3;					//ショルダー半径
			double proR=2.5e-3;				//ﾌﾟﾛｰﾌﾞ半径
			double proL=4e-3;				//ﾌﾟﾛｰﾌﾞ長さ
			double rpm=500;
			double rps=rpm/60;
			double w=rps*2*PI;				//角速度
			
			double U=CON->get_move_speed();	//プローブの移動速度[m/sec]		
			double pich=0.7e-3;				//プローブのねじのピッチ0.7[mm]
			double T=CON->get_roomT();		//初期温度
			double wall_h=T*(CON->get_wall_density()*le*le*le)*CON->get_wall_Cp();//壁粒子のエンタルピー

			int toBEM=OFF;//このモデルではＦＥＭは使用しないからどちらでもよい
			double val=0;
			int materialID=1;

			//calclation type
			int plunge=0;
			int traverse=1;
			int calc_type=0;
			if(CON->get_process_type()==2)	calc_type=0;	//1:plunge 2:traverse
			else
			{
				calc_type=CON->get_process_type();
			}
			double plunge_H=4*1e-3+1*le;	//plunge解析の際の、ツールを上げる高さ

			//流体 半径9mm、深さ6mm
			for(int i=-W;i<=W;i++)
			{
				for(int j=-W;j<=YW;j++)
				{
				
					for(int k=1;k<=fH;k++)
					{
						int flag=0;
						ii = i*le;
						jj = j*le;
						kk = k*le;
						double r=(ii*ii+jj*jj);
						if(calc_type==traverse)
						{
							if(r<proR*proR && kk>=fH*le-proL) flag=1;//ﾌﾟﾛｰﾌﾞ領域 plunge解析のときはすべて流体で良い
						}
						if(flag==0)
						{
							if(ii>0) materialID=1;
							else materialID=2;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
							writedata( fp, gnu, number,  ii,  jj, kk,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, toBEM);
							number++;
						}
					}
				}	
			}
			materialID=1;
			//壁INWALL側面
			for(int i=-W-2;i<=W+2;i++)
			{
				for(int j=-W-2;j<=YW+2;j++)
				{
					if(i*i>W*W || j<-W || j>YW)
					{
						for(int k=1;k<=fH;k++)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
							writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
							number++;
						}
					}
				}
			}
			//壁INWALL下
			for(int i=-W-2;i<=W+2;i++)
			{
				for(int j=-W-2;j<=YW+2;j++)
				{
					for(int k=0;k>=-1;k--)
					{
						ii = i*le;
						jj = j*le;
						kk = k*le;
						double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
						writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
						number++;
					}
				}
			}
				//プローブINWALL
			for(int i=-W;i<=W;i++)
			{
				for(int j=-W;j<=W;j++)
				{
					if(i*i+j*j<W*W)
					{
						for(int k=1;k<=fH;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							if(r<proR*proR && kk>=fH*le-proL) flag=1;//プローブ
							if(r<(proR-le)*(proR-le) && kk>=fH*le-proL+le) flag=0;
							if(flag==1)
							{
								double speed=0;
								double u=0;double v=0;double uw=0;
								if(r!=0)
								{
									speed=sqrt(r)*w;
									u=speed*(-jj/sqrt(r));
									v=speed*(ii/sqrt(r));
									if(r>=(proR-le)*(proR-le))
									{
										uw=-pich*rps;
									}
								}
								if(calc_type==plunge)
								{
									kk+=plunge_H;
									uw-=U;//Z方向に挿入速度を上乗せ
								}
								if(calc_type==traverse) v+=U;//y方向にツールを移動
								double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
								writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  u, v, uw, 0, wall_h, MOVE);
								number++;
							}
						}
					}
				}
			}
			//ショルダーINWALL
			for(int i=-shoR;i<=shoR;i++)
			{
				for(int j=-shoR;j<=shoR;j++)
				{
					if(i*i+j*j<=shoR*shoR)
					{
						for(int k=fH+1;k<=fH+2;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							
							double speed=0;
							double u=0;double v=0; double uw=0;
							if(r!=0)
							{
								speed=sqrt(r)*w;
								u=speed*(-jj/sqrt(r));
								v=speed*(ii/sqrt(r));
							}
							if(r>=(proR-le)*(proR-le) && r<proR*proR) uw=-pich*rps;//壁における速度発散ゼロ化のため、ショルダー内にもねじ速度導入
							if(calc_type==traverse) v+=U;//y方向にツールを移動
							if(calc_type==plunge)
							{
								kk+=plunge_H;
								uw-=U;//Z方向に挿入速度を上乗せ
							}
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
							writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  u, v, uw, 0, wall_h, MOVE);
							number++;
						}
					}
				}
		    }
		    //OUTWALL側面
			for(int i=-W-4;i<=W+4;i++)
		    {
				for(int j=-W-4;j<=YW+4;j++)
				{
					if(i*i>(W+2)*(W+2) || j<-W-2 || j>YW+2)
					{
						for(int k=-1;k<=fH;k++)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
							writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
							number++;
						}
					}
				}
			}
			//OUTWALL下面
			for(int i=-W-4;i<=W+4;i++)
			{
				for(int j=-W-4;j<=YW+4;j++)
				{
					//if(i*i+j*j<=(W+4)*(W+4))
					{
						for(int k=-2;k>=-3;k--)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
							writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
							number++;
						}
					}
				}
			}
			//プローブOUTWALL
			for(int i=-W;i<=W;i++)
			{
				for(int j=-W;j<=W;j++)
				{	
					if(i*i+j*j<W*W)
					{
						for(int k=1;k<=fH;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							if(r<(proR-le)*(proR-le) && kk>=fH*le-proL+le) flag=1;
							if(flag==1)
							{
								double speed=0;
								double u=0;double v=0; double uw=0;
								if(r!=0)
								{
									speed=sqrt(r)*w;
									u=speed*(-jj/sqrt(r));
									v=speed*(ii/sqrt(r));
									//uw=-pich*rps;
								}
								if(calc_type==traverse) v+=U;//y方向にツールを移動
								if(calc_type==plunge)
								{
									kk+=plunge_H;
									uw-=U;//Z方向に挿入速度を上乗せ
								}
								double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
								writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val, u, v, uw, 0, wall_h, MOVE);
								number++;
							}
						}
					}
				}
		    }
			//ショルダーOUTWALL
			for(int i=-shoR;i<=shoR;i++)
		    {
				for(int j=-shoR;j<=shoR;j++)
				{
					if(i*i+j*j<=shoR*shoR)
					{
						for(int k=fH+3;k<=fH+20;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							
							double speed=0;
							double u=0;double v=0; double uw=0;
							if(r!=0)
							{
								speed=sqrt(r)*w;
								u=speed*(-jj/sqrt(r));
								v=speed*(ii/sqrt(r));
							}
							if(r>=(proR-le)*(proR-le) && r<proR*proR) uw=-pich*rps;//壁における速度発散ゼロ化のため、ショルダー内にもねじ速度導入

							if(calc_type==traverse) v+=U;//y方向にツールを移動
							if(calc_type==plunge)
							{
								kk+=plunge_H;
								uw-=U;//Z方向に挿入速度を上乗せ
							}
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//エンタルピー
							writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val, u, v, uw, 0, wall_h, MOVE);
							number++;
						}
					}
				}
			}
		}
	}

	//ih釜
	if(model==26)
	{
		double le=CON->get_distancebp();
		double R=CON->get_fluidwidth()*0.001;//3;			//
		double Zg=CON->get_height();					//
		double C_h=0.1;//10;							//釜側面部の高さ
		double C_hout=C_h+4*le;//10;							//釜側面部の高さ
		double C_R=0.08;//3;			//釜底部の半径
		double C_Rout=C_R+4*le;
		double real_width=C_R/le;
		double real_height=C_h/le/2;

		int width=(int) real_width;
		int height=(int) real_height;

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
		
		double air_mass=V*1.205;//空気の重さ
		double air_Cp=1006;//空気の比熱
		double air_h=T*air_mass*air_Cp;//固体

		
		double val=0;

		if(CON->get_dimention()==2)
		{
			//////流体
			for(int i=-width;i<=width;i++)
			{
				
				for(int j=-height;j<=height;j++)
				{
					ii = i*CON->get_distancebp();
					jj = j*CON->get_distancebp();
					if(jj<=0.0)//0.04+0.01で50mm。半径80mmの円筒容器に1lのものを入れた際のおおよその高さ
					{
					writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
					//writedata2(fp,i,ii,j,Z[i],FLUID,materialID,1,0,0,0,0,0,h,1);
					number++;
					}
				}
			}

			if(CON->get_wall_poly()==0)
			{
				//////内壁
				//左
				for(int i=-width-1;i<=-width-1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=-height;j<=height+1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						//if(jj<=0.0)//0.04+0.01で50mm。半径80mmの円筒容器に1lのものを入れた際のおおよその高さ
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}

				//右
				for(int i=width+1;i<=width+1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=-height;j<=height+1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						//if(jj<=0.0)//0.04+0.01で50mm。半径80mmの円筒容器に1lのものを入れた際のおおよその高さ
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}

				//下
				for(int i=-width-1;i<=width+1;i++)
				{
					for(int j=-height-1;j<=-height-1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}

				//////
				if(CON->get_airwall()==1)
				{
					//上//表面流動を抑える用
					for(int i=-width;i<=width;i++)
					{
						for(int j=height+1;j<=height+1;j++)
						{
							ii = i*CON->get_distancebp();
							jj = j*CON->get_distancebp();
							writedata( fp, gnu, number,  ii,  jj-0.05+3*CON->get_distancebp(), 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, air_h, 0);//+3*CON->get_distancebp()
							number++;
						}
					}
					////*/
				}

				//////外壁
				//左
				for(int i=-width-5;i<=-width-2;i++)
				{
					for(int j=-height-5;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}
				/*/////
				for(int i=-width-1;i<=-width-1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=height+1;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						if(jj>0.0)//0.04+0.01で50mm。半径80mmの円筒容器に1lのものを入れた際のおおよその高さ
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}
				/////*/

				//右
				for(int i=width+2;i<=width+5;i++)
				{
					for(int j=-height-5;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}
				/*/////
				//右
				for(int i=width+1;i<=width+1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=height+1;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						if(jj>0.0)//0.04+0.01で50mm。半径80mmの円筒容器に1lのものを入れた際のおおよその高さ
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}
				////*/

				//下
				for(int i=-width-1;i<=width+1;i++)
				{
					for(int j=-height-5;j<=-height-2;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}
				if(CON->get_airwall()==ON)
				{
					//////
					//上
					for(int i=-width;i<=width;i++)
					{
						for(int j=height+2;j<=height+5;j++)
						{
							ii = i*CON->get_distancebp();
							jj = j*CON->get_distancebp();
							writedata( fp, gnu, number,  ii,  jj-0.05+3*CON->get_distancebp(), 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, air_h, 0);//+3*CON->get_distancebp()
							number++;
						}
					}
					//////*/
				}
			}
			
		}
	}

	//溶接
	if(model==27)
	{
		double le=CON->get_distancebp();
		double R=CON->get_fluidwidth()*0.001;//3;			//

		double real_width=0.025/le/2;
		double real_height=0.01/le/2;

		int width=(int) real_width;
		int height=(int) real_height;

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
		

		
		double val=0;

		if(CON->get_dimention()==2)
		{
			//////流体
			for(int i=-width;i<=width;i++)
			{
				
				for(int j=-height;j<=height;j++)
				{
					ii = i*CON->get_distancebp();
					jj = j*CON->get_distancebp();
					
					writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
					//writedata2(fp,i,ii,j,Z[i],FLUID,materialID,1,0,0,0,0,0,h,1);
					number++;
					
				}
			}

		}
	}

	
	///////////////////*/

	fp.close();
	gnu.close();
	*particle_number=number;
}