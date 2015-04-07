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
#include"MMM_CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい


///磁気モーメント法計算関数
void Magnetic_Moment_Method(mpsconfig *CON,vector<mpsparticle> &PART,double **F,double n0,double lamda,int fluid_number,int particle_number)
{
	
	double u0=4*PI*1e-7;
	MMM_config MCON;

	cout<<"磁気モーメント法---";

	double le=CON->get_distancebp();
	double R=CON->get_re()*le;	//発散用影響半径
	int d=CON->get_dimention();						//次元
	unsigned int timeA=GetTickCount();				//計算開始時刻
	int count=0;
	int pn=fluid_number*d;							//未知数:粒子数×速度D(次元)成分 
	double co=1/(4*PI*u0);							//計算によく現れる係数
	double RP0=CON->get_RP();						//比透磁率
	double kai=1.0-RP0;								//磁気感受率
	double V=CON->get_particle_volume();
	double S=le;
	if(d==3 && CON->get_model_set_way()==0) S=le*le;
	if(d==3 && CON->get_model_set_way()==1) S=sqrt(3.0)/4*le*le;//断面積

	//行列確保 //未知数の順番はMx0,My0,Mz0,Mx1,My1,Mz1,Mx2,My2,Mz2・・・
	double **matrix=new double*[pn];
	for(int i=0;i<pn;i++) matrix[i]=new double[pn];
	double *B   = new double[pn];					//解行列

	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];//外向き法線ベクトル

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) direct[D][i]=0;				//初期化
        if(PART[i].surface==ON)
		{
			direct_f(CON,PART,i,direct);
			for(int D=0;D<d;D++) direct[D][i]*=-1;//外向きが欲しいから反転する
		}
	}

	for(int i=0;i<pn;i++)
	{
		B[i]=0;
		for(int j=0;j<pn;j++) matrix[i][j]=0;		//初期化
	}
	//////////////////////////////////////////////////////////////////

	//解行列B[]作成
	if(MCON.get_Hf_type()==0)
	{
		double sign[3]={0,0,0};
		sign[d-1]=1;
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<d;D++) B[i*d+D]=MCON.get_Hf_H()*sign[D];//2次元ならA_Y,3次元ならA_Zの方向だけ値が入る。
		}
	}
	else cout<<"外部磁場のタイプが未実装"<<endl;
	/////*/

	
	/////////////////////matrixに値を格納
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<d;D++)
		{
			int n=i*d+D;//対応するmatrix内要素番号
			matrix[n][n]+=1/(RP0-1);
		}
		for(int j=0;j<fluid_number;j++)
		{
			if(j!=i)
			{
				double rA[3]={PART[i].r[A_X]-PART[j].r[A_X],PART[i].r[A_Y]-PART[j].r[A_Y],PART[i].r[A_Z]-PART[j].r[A_Z]};
				double disA=sqrt(rA[A_X]*rA[A_X]+rA[A_Y]*rA[A_Y]+rA[A_Z]*rA[A_Z]);//計算点iと計算点jの距離
	
				double dis3=disA*disA*disA;		//距離の3乗
				
				//if(PART[j].surface==OFF)
				{
					double co2=1/dis3*co*V;	//よく使う係数
					for(int D=0;D<d;D++)
					{
						int n=i*d+D;//対応するmatrix内要素番号
						for(int k=0;k<PART[j].N;k++)
						{
							int j2=PART[j].NEI[k];	//計算点jの周辺粒子
							double rB[3]={PART[j2].r[A_X]-PART[j].r[A_X],PART[j2].r[A_Y]-PART[j].r[A_Y],PART[j2].r[A_Z]-PART[j].r[A_Z]};
							double disB=sqrt(rB[A_X]*rB[A_X]+rB[A_Y]*rB[A_Y]+rB[A_Z]*rB[A_Z]);//計算点jとその近隣粒子との距離
							double w=kernel(R,disB);
							for(int D2=0;D2<d;D2++)
							{
								int m=j2*d+D2;
								int l=j*d+D2;
								matrix[n][m]+=d/n0*co2*rB[D2]/(disB*disB)*w*rA[D];
								matrix[n][l]+=-d/n0*co2*rB[D2]/(disB*disB)*w*rA[D];
							}
						}	
					}
				}
				/*else if(PART[j].surface==ON)
				{
					double co2=1/dis3*co*S;	//よく使う係数
					for(int D=0;D<d;D++)
					{
						int n=i*d+D;//対応するmatrix内要素番号
						for(int D2=0;D2<d;D2++)
						{
							int l=j*d+D2;
							matrix[n][l]+=-co2*direct[D2][j]*rA[D];
						}	
					}
				}///////*/
			}
		}
	}
	cout<<"行列作成"<<endl;

	//行列値出力
	//output_matrix(matrix,B, pn);

	cout<<"未知数:"<<pn<<"に対しガウスの消去法--";
	unsigned int timeB=GetTickCount();
	gauss(matrix,B,pn);//答えはBに格納
	//jacobi(matrix,B,pn);//答えはBに格納
	cout<<"ok time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;

	double *M[DIMENTION];					
	for(int D=0;D<DIMENTION;D++) M[D]=new double [fluid_number];
	double *H[DIMENTION];					
	for(int D=0;D<DIMENTION;D++) H[D]=new double [fluid_number];//粒子位置での磁場H
	
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			M[D][i]=0;
			H[D][i]=0;
		}
		for(int D=0;D<d;D++)
		{
			M[D][i]=B[i*d+D];
			H[D][i]=M[D][i]/kai;
		}
	}
	ofstream fp("M.dat");
	ofstream fq("H.dat");
	double times=1e-11;
	if(d==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<M[A_X][i]*times<<" "<<M[A_Y][i]*times<<endl;
			fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i])<<endl;
		}
	}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le)
			{
				fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]*times<<" "<<M[A_Z][i]*times<<endl;
				fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i])<<endl;
			}
		}
	}
	fp.close();
	fq.close();

	//力を求める
	
	double *Fs[DIMENTION];
    for(int D=0;D<DIMENTION;D++) Fs[D]=new double [particle_number];//単位面積あたりの力
	double *Fv[DIMENTION];
    for(int D=0;D<DIMENTION;D++) Fv[D]=new double [particle_number];//単位体積あたりの力
	double *Hgrad[DIMENTION];			//H勾配ﾍﾞｸﾄﾙ格納
	for(int D=0;D<DIMENTION;D++) Hgrad[D]=new double [fluid_number];
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) Fs[D][i]=0;		//初期化
		if(PART[i].surface==ON)
		{
			double Mn=M[A_X][i]*direct[A_X][i]+M[A_Y][i]*direct[A_Y][i]+M[A_Z][i]*direct[A_Z][i];
			double val=0.5*u0*Mn*Mn;		//応力値
			for(int D=0;D<DIMENTION;D++) Fs[D][i]=val*direct[D][i];
			for(int D=0;D<DIMENTION;D++) F[D][i]=Fs[D][i]*CON->get_distancebp();
		}
	}
	//体積力
	H_gradient1(CON,PART, fluid_number,Hgrad,H);//∇H計算
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) Fv[D][i]=u0*M[D][i]*Hgrad[D][i];
	
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=Fv[D][i]*V;

	////ｽﾑｰｼﾞﾝｸﾞ
	//smoothingF3D(CON,PART,fluid_number,F);

	ofstream fr("Fs.dat");
	times=5e-2/MCON.get_Hf_H();
	if(d==2){ for(int i=0;i<fluid_number;i++) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Y][i]*times<<endl;}
	else if(d==3)
	{
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Fs[A_X][i]*times<<" "<<Fs[A_Z][i]*times<<endl;
	}
	fr.close();
	ofstream ft("Fv.dat");
	times=1e-6;
	if(d==2) for(int i=0;i<fluid_number;i++) ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Fv[A_X][i]*times<<" "<<Fv[A_Y][i]*times<<endl;
	ft.close();

	if(CON->get_dir_for_P()==2 ||CON->get_dir_for_P()==3 )
    {
		ofstream bb("electromagnetic_P.dat");
		for(int i=0;i<fluid_number;i++)
		{
			double fs=0;//表面力
			if(PART[i].surface==ON)
			{
				fs=sqrt(Fs[A_X][i]*Fs[A_X][i]+Fs[A_Y][i]*Fs[A_Y][i]+Fs[A_Z][i]*Fs[A_Z][i]);
				for(int D=0;D<DIMENTION;D++) F[D][i]=0;
			}
			bb<<-fs<<endl;
        }
		bb.close();
	}

	for(int i=0;i<pn;i++) delete [] matrix[i];
	delete [] matrix;
	
    delete [] B;
	for(int D=0;D<DIMENTION;D++) delete [] M[D];
	for(int D=0;D<DIMENTION;D++) delete [] H[D];
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	for(int D=0;D<DIMENTION;D++) delete [] Fs[D];
	for(int D=0;D<DIMENTION;D++) delete [] Fv[D];
	for(int D=0;D<DIMENTION;D++) delete [] Hgrad[D];

	
}

///H勾配計算関数ver.1
void H_gradient1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **Hgrad,double **H)
{
	double le=CON->get_distancebp();//初期粒子間距離
	double r=CON->get_re()*le;
	int d=CON->get_dimention();

	double *HH=new double[fluid_number];//各粒子位置での磁場強度H格納
	for(int i=0;i<fluid_number;i++) HH[i]=sqrt(H[A_X][i]*H[A_X][i]+H[A_Y][i]*H[A_Y][i]+H[A_Z][i]*H[A_Z][i]);
	

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) Hgrad[D][i]=0;//初期化

		double W=0;//粒子数密度　OUTを除いたりするのでPND[i]は微妙

		for(int k=0;k<PART[i].N;k++)
		{       
			int j=PART[i].NEI[k];
			if(PART[j].type==FLUID)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
		
				double w=kernel(r,dis);
				W+=w;
				
				Hgrad[A_X][i]+=(HH[j]-HH[i])*X*w/(dis*dis);
				Hgrad[A_Y][i]+=(HH[j]-HH[i])*Y*w/(dis*dis);
				Hgrad[A_Z][i]+=(HH[j]-HH[i])*Z*w/(dis*dis);
			}

		}
		for(int D=0;D<DIMENTION;D++) if(W!=0) Hgrad[D][i]= Hgrad[D][i]*d/W;
	}///////////////Pgrad[D][i]計算終了

	

	delete [] HH;
}

//行列値出力関数(管理用)
void output_matrix(double **matrix,double *Bmatrix,int node_num)
{
	//node_num=30;					//左上だけ表示したい時
	int *DDN=new int[node_num];		//対角優位　diagonally dominant matrix
	ofstream fp2("matrixMMM.dat");
	for(int n=0;n<node_num;n++)
	{
		double val=0;
		DDN[n]=OFF;
		for(int m=0;m<node_num;m++)
		{
			fp2<<matrix[n][m]<<"\t";
			if(n!=m) val+=sqrt(matrix[n][m]*matrix[n][m]);
		}
		fp2<<endl;
		if(val<sqrt(matrix[n][n]*matrix[n][n])) DDN[n]=ON;//対角優位
	}
	fp2.close();
	ofstream fp4("BmatrixMMM.dat");
	for(int n=0;n<node_num;n++) fp4<<Bmatrix[n]<<endl;
	fp4.close();

	//for(int n=0;n<node_num;n++) if(DDN[n]==ON) cout<<n<<"は対角優位"<<endl;

	delete [] DDN;
}

//jacobiの反復法 解は最終的にBのなかへ
void jacobi(double **matrix,double *B,int N)
{
	double ep=1e-8;//収束判定
	double E=10;		//誤差
	int count=0;
	
	double *X=new double[N];//解
	for(int k=0;k<N;k++) X[k]=0;//初期値
	while(E>ep)
	{
		E=0;
		for(int i=0;i<N;i++)
		{
			double L=0;
			double U=0;
			for(int j=0;j<i;j++) L+=matrix[i][j]*X[j]; 
			for(int j=i+1;j<N;j++) U+=matrix[i][j]*X[j];
			double Xnew=(B[i]-L-U)/matrix[i][i];
			E+=fabs(Xnew-X[i]);
			X[i]=Xnew;
		}
		E/=N;
		count++;
		cout<<count<<" E="<<E<<endl;
	}
	
	for(int k=0;k<N;k++) B[k]=X[k];//Bに答えを格納
	delete [] X;
}