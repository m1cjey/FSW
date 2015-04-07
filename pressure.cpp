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

void P3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,int *jnb,vector<mpsparticle> &PART,int fluid_number,int **nei,double dt,double N0);
void reU3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,vector<mpsparticle> &PART,int fluid_number,double dt,int *jnb);


//陰解析
void negative1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lamda,double N0,double *PND2,double n0,double **Un)
{
	int d=CON->get_dimention();
	double *reU[DIMENTION];//速度修正量
	int negativeP=CON->get_negativeP();			//負の圧力を許可するか、しないか
	for(int D=0;D<DIMENTION;D++) reU[D] = new double [fluid_number];

	/////圧力計算
	pressure(CON,PART,fluid_number,particle_number,out,dt,t,lamda,N0,PND2,n0,CON->get_B_of_P(),negativeP);

	//圧力勾配計算
	calc_Pgradient(CON,PART,particle_number,fluid_number,reU,n0,dt,CON->get_minP());
	///////////////
	
	//速度修正と位置修正
	modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);
	
	for(int D=0;D<DIMENTION;D++) delete [] reU[D];
}

//陰解析を2回（１回目はPND　２回目は速度発散）
void negative1_twice(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lamda,double N0,double *PND2,double n0,double **Un,double n0_4)
{
	int d=CON->get_dimention();
	double le=CON->get_distancebp();
	double limit=CON->get_Pgrad_limit()*le;		//位置修正量の限界値
	int negativeP=CON->get_negativeP();			//負の圧力を許可するか、しないか
	double *reU[DIMENTION];//速度修正量
	for(int D=0;D<DIMENTION;D++) reU[D] = new double [fluid_number];

	//if(t%2==0)
	if(t%CON->get_P_twice()==0)
	{

		//calc_PND_by_minmum_L_main(CON,PART, particle_number);

		/////圧力計算
		//negativeP=OFF;
		pressure(CON,PART,fluid_number,particle_number,out,dt,t,lamda,N0,PND2,n0,0,negativeP);//粒子数密度一定条件で計算
	
		//圧力勾配計算
		int minPsw=OFF;
		calc_Pgradient(CON,PART,particle_number,fluid_number,reU,n0,dt,minPsw);

		//位置のみ修正
		int modified_num=0;						//位置修正量を修正された粒子数
		for(int i=0;i<fluid_number;i++)
		{
			double ABS=0;		//位置修正量の絶対値
			for(int D=0;D<d;D++) ABS+=dt*reU[D][i]*dt*reU[D][i];
			ABS=sqrt(ABS);
			if(ABS>limit)
			{
				for(int D=0;D<d;D++) reU[D][i]*=limit/ABS;		//修正量の大きさを修正
				modified_num++;
			}
		
			for(int D=0;D<d;D++) PART[i].r[D]+=dt*reU[D][i];
		}

		if(t==1 && CON->get_restart()==OFF)
		{
			ofstream fout("Pgrad_modified_num.dat");
	
			fout<<t<<" "<<modified_num<<endl;		
			fout.close();
		}
		else if(modified_num>0)
		{
			ofstream avs("Pgrad_modified_num.dat",ios :: app);
			avs<<t<<" "<<modified_num<<endl;	
			avs.close();
		}

		for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].r[D]+PART[i].r[D]!=2*PART[i].r[D]) cout<<"## i="<<i<<endl;
	
		//速度修正と位置修正
		//modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);
		
		//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].u[D]-=reU[D][i];//速度だけもとに戻す
	
		//粒子位置が変更されたので、freeonを実行しないといけない
		int *surface=new int[particle_number];		//ただし表面の定義は変更してほしくない。そこで値を保存しておく
		for(int i=0;i<particle_number;i++) surface[i]=PART[i].surface;
		calc_neighbor_relation(CON,PART,particle_number,n0_4,fluid_number,out);
		for(int i=0;i<particle_number;i++) PART[i].surface=surface[i];
		delete [] surface;
	}	

	//再度計算
	negativeP=CON->get_negativeP();		
	pressure(CON,PART,fluid_number,particle_number,out,dt,t,lamda,N0,PND2,n0,CON->get_B_of_P(),negativeP);

	//圧力勾配計算
	calc_Pgradient(CON,PART,particle_number,fluid_number,reU,n0,dt,CON->get_minP());
	///////////////
	
	//速度修正と位置修正
	modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);
	
	for(int D=0;D<DIMENTION;D++) delete [] reU[D];
	
}
///圧力計算関数
void pressure(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,double dt,int t,double lamda,double N0,double *PND2,double n0,int B_of_P,int negativeP)
{
	int P_flag=OFF;	//圧力計算を行うか否か
	int pn=0;		//未知数
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==OFF)
		{
			P_flag=ON;//圧力計算を行う
			pn++;
		}
		else if(PART[i].type==INWALL && PART[i].surface==OFF)
		{
			P_flag=ON;//圧力計算を行う
			pn++;
		}
	}
	
	if(P_flag==ON && fluid_number>0)
	{
		if(CON->get_gridless_P()==OFF) calc_P_main(CON,PART,particle_number,lamda,N0,dt,PND2,n0,fluid_number,out,pn,B_of_P);//MPSによる圧力計算
		if(CON->get_gridless_P()==ON) calc_P_main_with_gridless(CON,PART,particle_number,lamda,N0,dt,PND2,n0,fluid_number,out,pn,B_of_P,t);//gridless法による圧力計算
	} 

    ///圧力を強制決定
	//if(CON->get_set_P()!=OFF) set_P(CON,PART,particle_number,fluid_number);
	/////*/
	
	////圧力をsplot
	if(P_flag==ON) plot_P(CON ,PART,particle_number,t,fluid_number);
	////////////////////////////
}

void calc_P_main(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double lamda,double N0,double dt,double *PND2,double n0,int fluid_number,int out,int pn,int B_of_P)
{

	///圧力計算において、使用する重み関数は勾配などに使用するそれとは異なるものを用いる。
	//理由は、勾配に使用する重み関数は単なる重み平均なので何をしようしてもよいが、
	//PPE式の離散化に使用する重み関数は、粒子数密度∝密度となるように設定する必要があるから。
	//もし両者で同一の重み関数を使用したければ、kernel2()の中身をkernelと同一に書き換えること。

	//pn:圧力を解くための連立方程式未知数

	double le=CON->get_distancebp();
	double r2=CON->get_re2()*le;
	unsigned int timeA=GetTickCount();					//計算開始時刻
	double dimention=2;
	if(CON->get_dimention()==3) dimention=3;
	double *density=new double[particle_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=fluid_number;i<particle_number;i++) density[i]=CON->get_density();

	//圧力解析のためにN0とPART[i].PND2をkernel2()用におきかえる
	//set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	//cout<<"N0="<<N0<<endl;
	///////////////////*/

	cout<<"圧力未知数:"<<pn<<" ";

	int *ppn = new int[pn];					//行列における第n番目の未知数は粒子番号ppn[n]の粒子に相当
	int *link = new int [particle_number];	//粒子番号iはlink[i]番目の未知数
	double *B   = new double[pn];			//解行列
	
	int count=0;
	///ppn配列とlink配列作成
	
	double real_lamda=lamda;	//lamdaの値記憶
	if(CON->get_HL_sw()==ON) lamda=2.0*CON->get_dimention()*real_lamda;//高次の離散化の場合、lamda=2dと再定義すればλ/2dの項は消える。ただしそうすると両辺の値が大きくなりすぎるので、reak_lamdaを両辺にかけて小さくしている
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type!=OUTWALL && PART[i].surface==OFF)
		{
			ppn[count]=i;
			if(B_of_P==0)
			{
				B[count]=-density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());//教科書
			}
			else if(B_of_P==1)//速度発散
			{    
				//double div=divergence(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//divがエラーで非実数のときはMPSで計算する
				//cout<<"D="<<div<<" "<<PART[i].PND2<<endl;
				
				//B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i];
				B[count]=div*lamda*PART[i].PND2/(dt*2*CON->get_dimention())*density[i];
			}
			else if(B_of_P==2)//(教科書+速度発散)/2
			{
				double a=CON->get_w_div();double b=1;///各手法の重み
				//double div=divergence(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//divがエラーで非実数のときはMPSで計算する
				
				B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;
				B[count]/=a+b;
			}	
			else if(B_of_P==3)//重み関数の直接微分による速度発散
			{    
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
			}
			else if(B_of_P==4)//重み関数の直接微分による速度発散+PND
			{    
				double a=CON->get_w_div();double b=1;///各手法の重み
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;//教科書
				B[count]/=a+b;
			}
			else if(B_of_P==5)//(ni-nk)+(nk-n0)=div+pnd(nk-n0)
			{    
				//double div=Dndt(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//divがエラーで非実数のときはMPSで計算する
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());
			}
			//if(B[count]>800000) cout<<"B="<<i<<endl;
			link[i]=count;
			count++;
		}
		else link[i]=pn+1;//行列に含まれない粒子にはﾀﾞﾐｰとして(pn+1)を格納
	}
	lamda=real_lamda;//lamdaの値を戻す
	//////*/

	if(CON->get_dir_for_P()!=OFF && B_of_P!=0)//表面粒子の圧力として、ディリクレ値
	{
		int flag=CON->get_dir_for_P();
		if(flag==2 || flag==3)
		{
			if(CON->get_EM_method()==OFF) flag=1;//BEMによる計算を行わないのにflagが2や3なら間違いなので1に戻す
		}
		double *Dirichlet_P=new double [particle_number];
		for(int i=0;i<particle_number;i++) Dirichlet_P[i]=0;

		if(flag==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
		else if(flag==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		else if(flag==3)
		{
			set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		}
		
		for(int i=fluid_number;i<particle_number;i++)//壁表面粒子のDirichlet_P計算
		{
			//cout<<i<<endl;
			if(PART[i].surface==ON)//壁表面粒子なら
			{
				double P=0;
				double W=0;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k];
					if(PART[j].type==FLUID && PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel(r2,dis);
						P+=Dirichlet_P[j]*w;
						W+=w;
					}
				}
				if(W!=0) P/=W;
				Dirichlet_P[i]=P;
			}
		}

		ofstream gg2("Dirichlet_error.dat");
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k];
					int m=link[j];
					if(m<pn)//粒子jが未知数なら
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
				    
						double w=kernel2(r2,dis,dimention);
						if(CON->get_HL_sw()==ON) w*=real_lamda;
						B[m]-=w*Dirichlet_P[i];
						PART[i].P=Dirichlet_P[i];
					}
				}
			}
			else if(PART[i].type==FLUID && Dirichlet_P[i]!=0)
			{
				cout<<"内部粒子なのにﾃﾞｨﾘｸﾚ値が格納されています。値は"<<Dirichlet_P[i]<<endl;
				gg2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		gg2.close();

		//ﾌｧｲﾙ出力
		ofstream gg("Dirichlet_P.dat");
		if(dimention==2) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]<<endl;}
		else if(dimention==3) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;}
		gg.close();///*/

		//ディリクレ値をベクトル表示
		output_dirichlet_vector_files(CON,PART,fluid_number,flag,Dirichlet_P);

		delete [] Dirichlet_P;
	}///////////

	///解行列出力
	ofstream h("Bmatrix_for_P.dat");
	if(CON->get_dimention()==2) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<B[n]<<endl;
		}
	}
	if(CON->get_dimention()==3) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<B[n]<<endl;
		}
	}
	h.close();
	/////////////*/

	
	int number=0;			//係数行列の非ゼロ要素数
	for(int n=0;n<pn;n++)
	{   
	    number++;///自分自身を数にいれる
	    int i=ppn[n];//n番目の未知数は粒子i
	  
	    for(int k=0;k<PART[i].N2;k++)
	    {
	        int j=PART[i].NEI2[k];
			int m=link[j];
			if(m<pn) number++;
	    }
	}
	////numberが求まった
	
    double *val = new double [number];
	int *ind = new int [number];//非ゼロ要素の列番号格納配列
	int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
	
	
	/////////////////////val,ind ,ptrに値を格納
	
	if(CON->get_HL_sw()==OFF)//標準的な離散化
	{
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//値を保存
			ind[index]=n;
			index++;
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				if(m<pn)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					    
					double w=kernel2(r2,dis,dimention);
					val[index]=w;
					ind[index]=m;
					index++;
					W+=w;
				}
				else if(PART[j].surface==ON)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					
					//double w=kernel(r2,dis);
					double w=kernel2(r2,dis,dimention);
					W+=w; //ここでwを計算に組み込むとき、粒子jはﾃﾞｨﾘｸﾚ型、組み込まないと勾配ゼロのﾉｲﾏﾝ型 
				}
			}
			val[kk]=-W;
		}
		ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない
	}
	else if(CON->get_HL_sw()==ON)//高次の離散化 重み関数を変えたらここも書き変えないといけないことに注意
	{
		if(CON->get_dimention()==2) cout<<"2Dは非対応"<<endl;
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//値を保存
			//val[index]=-PART[i].PND;
			ind[index]=n;
			index++;
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				if(m<pn)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					    
					//double w=15*r2*r2*r2/(dis*dis*dis*dis*dis);
					double w=3*r2/(dis*dis*dis);				//R/dis-1
					val[index]=w*real_lamda;
					ind[index]=m;
					index++;
					W+=w;
				}
				else if(PART[j].surface==ON)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					//double w=15*r2*r2*r2/(dis*dis*dis*dis*dis);
					double w=3*r2/(dis*dis*dis);				//R/dis-1
					W+=w; //ここでwを計算に組み込むとき、粒子jはﾃﾞｨﾘｸﾚ型、組み込まないと勾配ゼロのﾉｲﾏﾝ型 
				}
			}
			val[kk]=-W*real_lamda;
		}
		ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない
	}
	////////////////////*/	 

	//ここで作られる行列は対角成分がすべて負なので、これを正にするために、係数行列と解行列に-1をかける
	for(int n=0;n<pn;n++)
	{
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			val[j]*=-1;
			if(ind[j]==n && val[j]<0) cout<<"対角成分が負 n="<<n<<endl;
		}
		B[n]*=-1;
	}
        
	/////////////////////////////////////CG法
    double *r=new double[pn];
	double *X=new double[pn];
	
	double *AP = new double [pn];
	double *P = new double [pn];

	/////////////////////////初期値//////////////////
    if(CON->get_initialP()==OFF)
	{
		for(int n=0;n<pn;n++) 
		{
			 X[n]=0;
			 r[n]=B[n];
			 P[n]=r[n];
		}
	}
	else if(CON->get_initialP()==ON)//初期値として現在の情報を与える。微小に早くなる
	{
		for(int n=0;n<pn;n++) X[n]=PART[ppn[n]].P;
		for(int n=0;n<pn;n++)
		{
			double AX=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
			r[n]=B[n]-AX;
			P[n]=r[n];
		}
	}
	//////////////////////////////////////////////

	//CG法により行列を解く
	if(CON->get_solution()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,CON->get_CGep());//CG法
	else if(CON->get_solution()==1)
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
		iccg(CON,val,ind,ptr,pn,B,number,X,r,P,CON->get_CGep(),&count);
	}
	
	if(CON->get_negativeP()==OFF)//負圧を考慮しない場合
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			PART[i].P=X[n];
			if(PART[i].P<0) PART[i].P=0;
		}
	}
	else if(CON->get_negativeP()==ON)//負圧を考慮する場合
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			//PART[i].P=X[n];
			if(2*X[n]==X[n]+X[n])PART[i].P=X[n];
			
		}
	}
	
		
    delete [] r;
	delete [] X;
	delete [] AP;
	delete [] P;

    //////////////////////////////*/

	//CG_GPU(CON,PART,val,ind,ptr,pn,number,ppn,B,&count);

    delete [] val;
	delete [] ind;
	delete [] ptr;
    delete [] B;
	delete [] ppn;
	delete [] link;

	delete [] density;

	cout<<"反復回数:"<<count<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//圧力解析のためにN0とPART[i].PND2をkernel2()用におきかえる
void set_N0_and_PND2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double *N0,int out)
{
	//0=教科書 1=速度発散 2=0+1 3=速度発散2 4=3+PND 5=(ni-nk)+(nk-n0)

	int SW=CON->get_B_of_P();			//PPE式における解行列 

	////////////////N0の計算
	int size = (int)(CON->get_re2()+1);	//計算領域
	double dis2;						//距離
	int d=CON->get_dimention();			//解析次元
	double R2=CON->get_re2()*CON->get_distancebp();
	double N=0;
	
	if(d==2)
	{
		if(CON->get_model_set_way()==0)//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis2=sqrt((double)(i*i+j*j));
					if(dis2!=0 && dis2<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis2,d);		
				}
			}
		}
		if(CON->get_model_set_way()==1)//初期配置として細密をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
					double dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis,d);			
				}
			}
		}
	}
	if(d==3)
	{
		if(CON->get_model_set_way()==0)//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis2=sqrt((double)(i*i+j*j+k*k));
						if(dis2!=0 && dis2<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis2,d);
					}			
				}
			}
		}
		if(CON->get_model_set_way()==1)//初期配置として細密をとった場合
		{
			for(int i=-2*size;i<=2*size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//kが奇数ならiiとjjをずらす
						double dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis<=CON->get_re2() ) N+=kernel2(CON->get_re2(),dis,d);
					}
				}
			}
		}
	}
	//cout<<"N0="<<N<<endl;
	*N0=N;
	///////////////////*/

	if(SW==0 || SW==2 || SW==4 || SW==5)//粒子数密度を使用する場合は
	{
		///PART[i].PND2の計算しなおし
		for(int i=0;i<particle_number;i++)//outwallはPART[i].N2が0なので、PND2は0になるが、問題はない
		{
			PART[i].PND2=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				    
				double w=kernel2(R2,dis,d);
				PART[i].PND2+=w;
			}
		}
	}///////////*/
}

//ディリクレ値セット関数
void set_Dirichlet_P(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int flag,double *Dirichlet_P)
{
	//fileを開く
	ifstream fin;
	/*/if(flag==1) fin.open("surface_tension_P.dat", ios::in);
	if(flag==2)
	{fin.open("electromagnetic_P.dat", ios::in);
		
		if(!fin) cout<<"cannot open the file for dir_for_P"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
	}
	///////////////////////////////*/
	

	double *val=new double[fluid_number];

	int *err=new int[fluid_number];//err[i]=ONなら、その表面粒子は対応するディリクレ値がゼロということ。主に電磁力がFEM_interval()分の1回しか実行されないことに起因するエラー。
    int modify_sw=OFF;		//modify_sw=ONならエラー発生。対処する。

	if(flag==1)
	{
		for(int i=0;i<fluid_number;i++)
		{
			Dirichlet_P[i]+=PART[i].dir_Pst;	//表面張力項
			val[i]=PART[i].dir_Pst;
			//if(PART[i].surface==1  && abs(val[i])<0.1) cout<<"error tension有りで表面なのにディリクレ値0? cood"<<endl;
		}
	}
	
	if(flag==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			err[i]=OFF;
			val[i]=PART[i].dir_Pem;

			//fin>>val[i];
			Dirichlet_P[i]+=val[i];//ﾌｧｲﾙからﾃﾞｨﾘｸﾚ値読み取り	
			if(flag==2)
			{
				if(PART[i].type==BOFLUID && val[i]==0)
				{
					//cout<<"val["<<i<<"]=0"<<endl;
					err[i]=ON;
					modify_sw=ON;//修正スイッチON
				}
			}
		}
	}

	//if(flag==2) fin.close();

	//エラー修正
	if(modify_sw==ON)
	{
		cout<<endl<<"set_Dirichlet_P()にてエラーを感知。重み補正実行。CON.FEM_intervalの縮小を要求。"<<endl;
		double le=CON->get_distancebp();
		double R=CON->get_re()*le;//修正に使用する影響半径
		ofstream fp("before_val.dat");
		if(CON->get_dimention()==2) {for(int i=0;i<fluid_number;i++) if(PART[i].type==BOFLUID) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<val[i]<<endl;}
		else if(CON->get_dimention()==3) {for(int i=0;i<fluid_number;i++) if(PART[i].type==BOFLUID) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<val[i]<<endl;}
		fp.close();

		ofstream fq("after_val.dat");
		for(int i=0;i<fluid_number;i++)
		{
			if(err[i]==ON)
			{
				double W=0;//重みの総和
				double value=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type==FLUID&& PART[i].surface==ON && err[j]==OFF)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel(R,dis);
						value+=val[j]*w;
						W+=w;
					}
				}
				if(W>0) value/=W;
				val[i]=value;
				Dirichlet_P[i]+=val[i];//エラーだったということは、修正前のval[i]=0を意味する。違うならこの行はそれを考慮すること。
			}
		}
		if(CON->get_dimention()==2) {for(int i=0;i<fluid_number;i++) if(PART[i].type==FLUID&& PART[i].surface) fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<val[i]<<endl;}
		else if(CON->get_dimention()==3) {for(int i=0;i<fluid_number;i++) if(PART[i].type==FLUID&& PART[i].surface) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fq<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<val[i]<<endl;}
		fq.close();
	}


	delete [] val;
	delete [] err;
}

//ディリクレ値表面力出力関数
void output_dirichlet_vector_files(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int flag,double *Dirichlet_P)
{
	//表面力ベクトルを出力する。単位はN/m^2 flag=1なら表面張力、flag=2なら電磁力
	
	double le=CON->get_distancebp();
	int dimention=CON->get_dimention();
	double xp=-100;	//判例を出すX座標
	double yp=-100;	//判例を出すY座標
	double maxP=0;	//絶対値の最大ディリクレ値

	//法線ベクトル
	double *direct[DIMENTION];
	for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];//内向き法線ベクトル
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
		else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
	}

	//maxP求める
	for(int i=0;i<fluid_number;i++) if(fabs(Dirichlet_P[i])>maxP) maxP=fabs(Dirichlet_P[i]);

	double times=4*le/maxP;
	if(CON->get_dir_for_P()==3) times*=0.5;//両方の場合は倍率半分
	

	ofstream gg2("dirichlet_vector.dat");
	if(dimention==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			gg2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]*direct[A_X][i]*times<<" "<<Dirichlet_P[i]*direct[A_Y][i]*times<<endl;
			if(PART[i].r[A_X]>xp) xp=PART[i].r[A_X];
			if(PART[i].r[A_Z]>yp) yp=PART[i].r[A_Y];
		}
		xp+=4*le;
		yp+=4*le;
		gg2<<xp<<"\t"<<yp<<"\t"<<maxP*times<<" "<<0<<" ////凡例の長さは"<<maxP<<endl;
	}
	if(dimention==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le)
			{
				gg2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]*direct[A_X][i]*times<<" "<<Dirichlet_P[i]*direct[A_Z][i]*times<<endl;
				if(PART[i].r[A_X]>xp) xp=PART[i].r[A_X];
				if(PART[i].r[A_Z]>yp) yp=PART[i].r[A_Z];
			}
		}
		xp+=4*le;
		yp+=4*le;
		gg2<<xp<<"\t"<<yp<<"\t"<<maxP*times<<" "<<0<<" ////凡例の長さは"<<maxP<<endl;
	}
	gg2.close();

	for(int D=0;D<DIMENTION;D++) delete [] direct[D]; 
}

///粒子数密度の時間微分計算関数
double Dndt(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0)
{
	///ここで返されるDndtは、重み関数pow(R,D)/pow(r,D)を用いてMPSにより速度の発散を計算した値に等しい
    double R=CON->get_distancebp()*CON->get_re2();	//影響半径
    double Dndt=0;									//導関数の値
	double D=2;										//解析次元。ただし計算の都合上、double型
	if(CON->get_dimention()==3) D=3;
    
    for(int k=0;k<PART[i].N2;k++)
    {    
        int j=PART[i].NEI2[k]; 
        double X=PART[j].r[A_X]-PART[i].r[A_X];
		double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
		double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
		double dis=sqrt(X*X+Y*Y+Z*Z);
		
		double div=(PART[j].u[A_X]-PART[i].u[A_X])*X+(PART[j].u[A_Y]-PART[i].u[A_Y])*Y+(PART[j].u[A_Z]-PART[i].u[A_Z])*Z;
		Dndt-=D*pow(R,D)/pow(dis,D)*div/(dis*dis);	
    }
    return Dndt;
}

///圧力プロット関数（スカラー）
void plot_P(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int t,int fluid_number)
{
	ofstream fout("P.dat");
    double le=CON->get_distancebp();

    if(CON->get_dimention()==2)
    {
		for(int i=0;i<particle_number;i++) fout<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].P<<endl;
    }
    else if(CON->get_dimention()==3)
    {
        for(int i=0;i<particle_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fout<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].P<<endl;
    }
    fout.close();

	ofstream fout2("Pf.dat");
    le=CON->get_distancebp();

    if(CON->get_dimention()==2)
    {
		for(int i=0;i<fluid_number;i++) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].P<<endl;
    }
    else if(CON->get_dimention()==3)
    {
        for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].P<<endl;
    }
    fout.close();

	if(CON->get_P_AVS()>0)//ここでif(CON->get_P_AVS()>0 &&  t%CON->get_P_AVS())とすると0のときバグる
	{
		if(t==1 || t%CON->get_P_AVS()==0) output_pressuer_avs(CON,PART,t,particle_number,fluid_number);
	}
}

//圧力AVSファイル出力関数
void output_pressuer_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	t=1;//いまはわざと毎ステップ上書き

	//sprintf_s(filename,"pressure/pressure%d",t);//フォルダを作成して管理する場合はこちら
	sprintf_s(filename,"pressure%d",t);//他のファイルと同じ階層に生成するならこちら
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
				double P=PART[i].P;
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
			double P=PART[i].P;
			fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
			n++;
		}
	}
	fout.close();
	//sprintf_s(filename,"pressure/pressure%d.fld",t);//フォルダを作成して管理する場合はこちら
	sprintf_s(filename,"pressure%d.fld",t);//他のファイルと同じ階層に生成するならこちら
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
	fout2 << "label=pressure" << endl << endl;
	//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout2 << "variable 1 file=pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    1 file=pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    2 file=pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    3 file=pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2.close();
}

///粒子データ読み取り関数
void calc_Pgradient(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double **reU,double n0,double dt,int minPsw)
{
	if(minPsw==ON) cout<<"圧力勾配計算開始 ver."<<CON->get_Pgrad()<<" minP=ON ---------";
	if(minPsw==OFF) cout<<"圧力勾配計算開始 ver."<<CON->get_Pgrad()<<" minP=OFF ---------";
	double le=CON->get_distancebp();	//初期粒子間距離
	unsigned int timeA=GetTickCount();	//計算開始時刻

	double *direct[DIMENTION];			//内向き法線ﾍﾞｸﾄﾙ格納
	for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];

	double *Pgrad[DIMENTION];			//圧力勾配ﾍﾞｸﾄﾙ格納
	for(int D=0;D<DIMENTION;D++) Pgrad[D]=new double [fluid_number];
	
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			direct[D][i]=0.0;
			Pgrad[D][i]=0.0;
		}
	}

	//////法線ﾍﾞｸﾄﾙ計算開始
	if(CON->get_Pgrad()==2 ||CON->get_Pgrad()==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON ) direct_f(CON,PART,i,direct);
			else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
		}
		
	}/////////////////*/
	
	///////////////////圧力勾配計算
	
	if(CON->get_Pgrad()==3)//表面のみ法線 && minP=0のとき表面粒子からも反発力計算
	{
		P_gradient3(CON,PART,fluid_number,direct,dt,reU,Pgrad,minPsw);
	}
	else if(CON->get_Pgrad()==4)//重みつき最少二乗法(WLSM)
	{
		P_gradient4(CON,PART,dt,fluid_number,reU,Pgrad,minPsw);
	}
	else if(CON->get_Pgrad()==5)//重みつき最少二乗法(WLSM)
	{
		P_gradient5(CON,PART,dt,fluid_number,reU,Pgrad);
	}
	///////////////////////*/

	//速度修正量計算
	double limit_reU=CON->get_Pgrad_limit()*CON->get_distancebp();
	double *density=new double[fluid_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) reU[D][i]=-dt*Pgrad[D][i]/density[i];
	delete [] density;

	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	for(int D=0;D<DIMENTION;D++) delete [] Pgrad[D];

	cout<<"ok  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

///圧力勾配計算関数ver.3
void P_gradient3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *direct[DIMENTION],double dt,double **reU,double **Pgrad,int minPsw)
{
	double le=CON->get_distancebp();//初期粒子間距離
	double r=CON->get_re()*le;
	int d=CON->get_dimention();
	
	double *minP=new double[fluid_number];	//周囲の最少圧力格納

	///minPを求める
	set_minP(CON,PART,fluid_number,minP,minPsw);
	///////////////

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++) Pgrad[D][i]=0;//初期化

		double W=0;//粒子数密度　OUTを除いたりするのでPND[i]は微妙

		if(PART[i].surface==ON)//表面粒子の場合、少し内側の圧力を調べて勾配を計算する。
		{
			double x1=PART[i].r[A_X]+direct[A_X][i]*le;//少し内側の座標
			double y1=PART[i].r[A_Y]+direct[A_Y][i]*le;
			double z1=PART[i].r[A_Z]+direct[A_Z][i]*le;
			double P=0;//少し内側の圧力
			
			for(int k=0;k<PART[i].N3;k++)
			{       
	 			int j=PART[i].NEI3[k];
				if(PART[j].type!=OUTWALL)//BDWALLも省いたほうがよくない？
				{
					double X=PART[j].r[A_X]-x1;
					double Y=PART[j].r[A_Y]-y1;
					double Z=PART[j].r[A_Z]-z1;
					double dis=sqrt(X*X+Y*Y+Z*Z);
					if(dis<r)
					{
						double w=(1-dis/r)*(1-dis/r);//w(0)=∞となる関数は使えない
						W+=w;
						P+=PART[j].P*w;			
					}
				}
			}

			if(W!=0) P/=W; //(x1,y1,z1)での圧力
			Pgrad[A_X][i]=(P-PART[i].P)*direct[A_X][i]/le;
			Pgrad[A_Y][i]=(P-PART[i].P)*direct[A_Y][i]/le;
			Pgrad[A_Z][i]=(P-PART[i].P)*direct[A_Z][i]/le;
		}
		else//内部流体粒子の場合、周囲から圧力勾配を計算
		{
			for(int k=0;k<PART[i].N;k++)
			{       
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					
					double w=kernel(r,dis);
					W+=w;
					if(minP[i]==0 && PART[j].type==FLUID && PART[j].surface==ON)///粒子jの圧力を粒子iの圧力でおきかえる
					{   
						Pgrad[A_X][i]+=(PART[i].P-minP[i])*X*w/(dis*dis);
						Pgrad[A_Y][i]+=(PART[i].P-minP[i])*Y*w/(dis*dis);
						Pgrad[A_Z][i]+=(PART[i].P-minP[i])*Z*w/(dis*dis);
					}
					else
					{   //通常どおり
						Pgrad[A_X][i]+=(PART[j].P-minP[i])*X*w/(dis*dis);
						Pgrad[A_Y][i]+=(PART[j].P-minP[i])*Y*w/(dis*dis);
						Pgrad[A_Z][i]+=(PART[j].P-minP[i])*Z*w/(dis*dis);
					}
				}
			}
			
			for(int D=0;D<DIMENTION;D++) if(W!=0) Pgrad[D][i]= Pgrad[D][i]*d/W;
		}
	}///////////////Pgrad[D][i]計算終了

	///////////////////////////////圧力勾配を表示
	if(CON->get_iteration_count()==1)				//陰解析を複数回行う場合は、毎回出力してしまうのでやめる
	{
		plot_Pgradient(CON,PART,fluid_number,Pgrad);
	}
	////////////////////////*/

	delete [] minP;
}

//最少圧力計算関数
void set_minP(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *minP,int minPsw)
{
	///minPを求める
	for(int i=0;i<fluid_number;i++) minP[i]=PART[i].P;
	if(minPsw==ON)
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int k=0;k<PART[i].N;k++)
			{       
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					if(PART[j].P<minP[i]) minP[i]=PART[j].P;
				}
			}
		}
	}////////////////////*/
}

///圧力勾配表示関数
void plot_Pgradient(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **Pgrad)
{
	int d=CON->get_dimention();
	double times=0;					//表示用倍率
	double le=CON->get_distancebp();
	double density=CON->get_density();

	//通常通り
	times=CON->get_times()*le*le/density*CON->get_Pgrad_times();
	
	ofstream gra("Pgrad.dat");
	ofstream gra2("Pgrad2.dat");//スカラー表示
	if(d==2)//見やすいように符号を反転表示
	{
		for(int i=0;i<fluid_number;i++)
		{
			gra<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<-Pgrad[A_X][i]*times<<" "<<-Pgrad[A_Y][i]*times<<endl;
			gra2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<sqrt(Pgrad[A_X][i]*Pgrad[A_X][i]+Pgrad[A_Y][i]*Pgrad[A_Y][i])<<endl;
		}
	}
	else if(d==3)//見やすいように符号を反転表示
	{
		for(int i=0;i<fluid_number;i++) 
		{
			if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le)
			{
				gra<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<-Pgrad[A_X][i]*times<<" "<<-Pgrad[A_Z][i]*times<<endl;
				gra2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<sqrt(Pgrad[A_X][i]*Pgrad[A_X][i]+Pgrad[A_Y][i]*Pgrad[A_Y][i]+Pgrad[A_Z][i]*Pgrad[A_Z][i])<<endl;
			}
		}
	}
	gra.close();
	gra2.close();
}

///圧力勾配計算関数ver.4 重みつき最少二乗法に基づく勾配(WLSM)
void P_gradient4(mpsconfig *CON,vector<mpsparticle> &PART,double dt,int fluid_number,double *reU[3],double *P_grad[3],int minPsw)
{
	//WLSMに関する注意点
	//・minP=OFFのほうがうまくいく??
	//・2D解析で2次精度の場合、re=2.1だと特に表面粒子の計算において近隣粒子数が少なく、計算結果がおかしいときがある。re=2.5くらいを推奨??。
	//・↑と同じ理由で、孤立粒子やそれに近い粒子の計算は破綻するおそれあり。

	double le=CON->get_distancebp();
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N=0;							//係数行列の元
	int order=CON->get_Pgrad_order();	//近似曲面のｵｰﾀﾞｰ。 1=線形 2=二次

	double *minP=new double[fluid_number];	//周囲の最少圧力格納

	///minPを求める
	set_minP(CON,PART,fluid_number,minP,minPsw);
	///////////////

	//係数行列の大きさの決定
	if(d==2)
	{
		if(order==1) N=2;
		else if(order==2) N=5;
	}
	else if(d==3)
	{
		if(order==1) N=3;
		else if(order==2) N=9;
	}
	////////////////////////////////

	double *matrix=new double [N*N];	//N×Nの係数行列
	double *B=new double [N];	//Nの解行列
	
	if(d==2 && order==1)//二次元
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			for(int n=0;n<N*N;n++) matrix[n]=0;//初期化
			for(int n=0;n<N;n++) B[n]=0;

			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					double X=(PART[j].r[A_X]-PART[i].r[A_X])/le;// leで割るのは打ち切り誤差防止
					double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/le;
					double dis=sqrt(X*X+Y*Y);
					
					double w=1;
					if(dis>1) w=1/(dis*dis*dis*dis);
					
					matrix[0]+=X*X*w;			//ΣXjwj
					matrix[1]+=X*Y*w;		//ΣXjYjwj
					matrix[3]+=Y*Y*w;			//ΣYjwj
				
					B[0]+=(PART[j].P-minP[i])*X*w;//ΣfjXjwj
					B[1]+=(PART[j].P-minP[i])*Y*w;//ΣfjYjwj
				}
			}
			matrix[2]=matrix[1];		//ΣXjYjwj
			for(int n=0;n<N;n++) B[n]/=le;//打ち切り誤差防止

			double determinant=matrix[0]*matrix[3]-matrix[1]*matrix[2];//行列式
			
			P_grad[A_X][i]=(B[0]*matrix[3]-matrix[1]*B[1])/determinant;
			P_grad[A_Y][i]=(B[1]*matrix[0]-matrix[2]*B[0])/determinant;

			int flag=OFF;//OFFなら問題ない。ONならMPSでやり直す
			for(int D=0;D<2;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//非実数ならON
			if(flag==ON)
			{
				//cout<<"エラーのためMPSにより離散化 "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				double grad[3]={0,0,0};
				double n0=1;//ここでは無意味な変数なので適当にいれとく
				P_gradient_MPS(CON,PART,i,n0,grad);//
				for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
			}
			
		}
	}
	else if(d==2 && order==2)//二次元2次式
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;

			for(int n=0;n<N*N;n++) matrix[n]=0;//初期化
			for(int n=0;n<N;n++) B[n]=0;
		
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				//if(PART[j].type!=OUTWALL)
				if(PART[j].type==FLUID)
				{
					double X=(PART[j].r[A_X]-PART[i].r[A_X]);
					double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
					double dis=sqrt(X*X+Y*Y);
					double dP=PART[j].P-minP[i];
						
					double w=1;
					
					if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
						
					matrix[0]+=X*X*w;			//ΣXjwj
					matrix[1]+=X*Y*w;		//ΣXjYjwj
					matrix[2]+=X*X*X*w;			//ΣXj^3wj
					matrix[3]+=X*X*Y*w;			//ΣXj^2Yjwj
					matrix[4]+=X*Y*Y*w;			//ΣXjYj^2wj
	
					matrix[6]+=Y*Y*w;			//ΣYj^2wj
					matrix[9]+=Y*Y*Y*w;			//ΣYj^3wj
	
					matrix[12]+=X*X*X*X*w;			//ΣXj^4wj
					matrix[13]+=X*X*X*Y*w;			//ΣXj^3Yjwj
					matrix[14]+=X*X*Y*Y*w;			//ΣXj^2Yj^2wj
		
					matrix[19]+=X*Y*Y*Y*w;			//ΣXjYj^3wj
	
					matrix[24]+=Y*Y*Y*Y*w;			//ΣYj^4wj

				
					B[0]+=dP*X*w;//ΣdPjXjwj
					B[1]+=dP*Y*w;//ΣdPjYjwj
					B[2]+=dP*X*X*w;//ΣdPjXj^2wj
					B[3]+=dP*X*Y*w;//ΣdPjXjYjwj
					B[4]+=dP*Y*Y*w;//ΣdPjYj^2wj
					//if(i==551) cout<<endl<<dP*Y*w<<" "<<dP<<" "<<Y<<" "<<w<<endl;
				}
			}
			//if(i==551) cout<<"B[1]="<<B[1]<<endl;
			
			matrix[5]=matrix[1];		//ΣXjYjwj
			matrix[7]=matrix[3];		//ΣXj^2Yjwj
			matrix[8]=matrix[4];		//ΣXjYj^2wj
			matrix[10]=matrix[2];		//ΣXj^3Yjwj
			matrix[11]=matrix[3];		//ΣXj^2Yjwj
			matrix[15]=matrix[3];		//ΣXj^2Yjwj
			matrix[16]=matrix[4];		//ΣXjYj^2wj
			matrix[17]=matrix[13];		//ΣXj^3Yjwj
			matrix[18]=matrix[14];		//ΣXj^2Yj^2wj
			matrix[20]=matrix[4];		//ΣXjYj^2wj
			matrix[21]=matrix[9];		//ΣYj^3wj
			matrix[22]=matrix[14];		//ΣXj^2Yj^2wj
			matrix[23]=matrix[19];		//ΣXjYj^3wj

			/*/丸め誤差防止
			for(int n=0;n<N;n++)
			{
				for(int m=0;m<N;m++) matrix[n*N+m]*=1e14;
				B[n]*=1e14;
			}//*/

			double dPdx=0;
			double dPdy=0;
			
			return_X_for5N(matrix,N,B,B,&dPdx,&dPdy);//5次連立方程式の、解１，２を返す関数(ここではBを２つ渡す)
		
			P_grad[A_X][i]=dPdx;
			P_grad[A_Y][i]=dPdy;
		}
			
	}
	else if(d==3 && order==1)//3次元1次式
	{
		///係数行列は
		///   ΣΔx2    ΣΔxΔy  ΣΔxΔz  a = ΣΔxΔf  
		///  ΣΔxΔy    ΣΔy2   ΣΔyΔz  b = ΣΔyΔf 
		///  ΣΔxΔz   ΣΔyΔz   ΣΔz2   c = ΣΔzΔf 
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			for(int n=0;n<N*N;n++) matrix[n]=0;//初期化
			for(int n=0;n<N;n++) B[n]=0;

			//if(PART[i].N>4)//周辺粒子数が極端に少ないと行列式がゼロになって解けない
			if(PART[i].N>=3)
			{
				cacl_WLSM_P_D3_order1(CON,PART,matrix,B,P_grad, i,N,minP);//1次近似 本来は3個の近隣粒子で計算可能だけど、保険をかけてここでは5以上の計算点が必要としている。
				int flag=OFF;//OFFなら問題ない。ONならMPSでやり直す
				for(int D=0;D<DIMENTION;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//非実数ならON
				if(flag==ON)
				{
					//cout<<"エラーのためMPSにより離散化 "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					double grad[3]={0,0,0};
					double n0=1;//ここでは無意味な変数なので適当にいれとく
					P_gradient_MPS(CON,PART,i,n0,grad);//
					for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
				}
			}
			else//周辺粒子の重み平均で. 粒子数が少ないので、MPSはちがうかと
			{
				double W=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					double X[3];
					for(int D=0;D<3;D++) X[D]=PART[j].r[D]-PART[i].r[D];
					double dis=sqrt(X[A_X]*X[A_X]+X[A_Y]*X[A_Y]+X[A_Z]*X[A_Z]);
					double dP=PART[j].P-PART[i].P;
					double w=kernel(R,dis);
					W+=w;
					for(int D=0;D<3;D++) P_grad[D][i]+=dP*X[D]/dis*w;
				}
				if(W!=0) for(int D=0;D<3;D++) P_grad[D][i]/=W;
				cout<<"警告 近隣粒子数が4以下("<<PART[i].N<<")です"<<endl;
			}
		}
	}
	else if(d==3 && order==2)//3次元2次式
	{
		//P=Pi+aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+iΔzΔxとおくと、
		///係数行列は
		///   ΣΔx2      ΣΔxΔy    ΣΔxΔz    ΣΔx3       ΣΔxΔy2    ΣΔxΔz2    ΣΔx2Δy     ΣΔxΔyΔz  ΣΔx2Δz     a = ΣΔxΔP  
		///   ΣΔxΔy    ΣΔy2      ΣΔyΔz    ΣΔx2Δy    ΣΔy3       ΣΔyΔz2    ΣΔxΔy2     ΣΔy2Δz    ΣΔxΔyΔz   b = ΣΔyΔP
		///   ΣΔxΔz    ΣΔyΔz    ΣΔz2      ΣΔx2Δz    ΣΔy2Δz    ΣΔz3       ΣΔxΔyΔz   ΣΔyΔz2    ΣΔxΔz2     c = ΣΔzΔP
		///   ΣΔx3      ΣΔx2Δy   ΣΔx2Δz   ΣΔx4       ΣΔx2Δy2   ΣΔx2Δz2   ΣΔx3Δy     ΣΔx2ΔyΔz ΣΔx3Δz     d = ΣΔx2ΔP
		///   ΣΔxΔy2   ΣΔy3      ΣΔy2Δz   ΣΔx2Δy2   ΣΔy4       ΣΔy2Δz2   ΣΔxΔy3     ΣΔy3Δz    ΣΔxΔy2Δz  e = ΣΔy2ΔP
		///   ΣΔxΔz2   ΣΔyΔz2   ΣΔz3      ΣΔx2Δz2   ΣΔy2Δz2   ΣΔz4       ΣΔxΔyΔz2  ΣΔyΔz3    ΣΔxΔz3     f = ΣΔz2ΔP
		///   ΣΔx2Δy   ΣΔxΔy2   ΣΔxΔyΔz ΣΔx3Δy    ΣΔxΔy3    ΣΔxΔyΔz2 ΣΔx2Δy2    ΣΔxΔy2Δz ΣΔx2ΔyΔz  g = ΣΔxΔyΔP
		///   ΣΔxΔyΔz ΣΔy2Δz   ΣΔyΔz2   ΣΔx2ΔyΔz ΣΔy3Δz    ΣΔyΔz3    ΣΔxΔy2Δz  ΣΔy2Δz2   ΣΔxΔyΔz2  h = ΣΔyΔzΔP
		///   ΣΔx2Δz   ΣΔxΔyΔz ΣΔxΔz2   ΣΔx3Δz    ΣΔxΔy2Δz  ΣΔxΔz3   ΣΔx2Δy Δz ΣΔxΔyΔz2 ΣΔx2Δz2    g = ΣΔxΔzΔP

		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			for(int n=0;n<N*N;n++) matrix[n]=0;//初期化
			for(int n=0;n<N;n++) B[n]=0;

			if(PART[i].N>8) cacl_WLSM_P_D3_order2(CON,PART,matrix,B,P_grad, i,N,minP);//2次近似 本来は6個の近隣粒子で計算可能だけど、保険をかけてここでは9以上の計算点が必要としている。
			else if(PART[i].N>4) cacl_WLSM_P_D3_order1(CON,PART,matrix,B,P_grad, i,3,minP);//1次近似
			int flag=OFF;//OFFなら問題ない。ONならMPSでやり直す
			for(int D=0;D<DIMENTION;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//非実数ならON
			if(flag==ON)
			{
				cout<<"エラーのためMPSにより離散化 "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				double grad[3]={0,0,0};
				double n0=1;//ここでは無意味な変数なので適当にいれとく
				P_gradient_MPS(CON,PART,i,n0,grad);//
				for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
			}
			if(PART[i].N<=4)
			{
				double W=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					double X[3];
					for(int D=0;D<3;D++) X[D]=PART[j].r[D]-PART[i].r[D];
					double dis=sqrt(X[A_X]*X[A_X]+X[A_Y]*X[A_Y]+X[A_Z]*X[A_Z]);
					double dP=PART[j].P-PART[i].P;
					double w=kernel(R,dis);
					W+=w;
					for(int D=0;D<3;D++) P_grad[D][i]+=dP*X[D]/dis*w;
				}
				if(W!=0) for(int D=0;D<3;D++) P_grad[D][i]/=W;
				cout<<"警告 近隣粒子数が4以下("<<PART[i].N<<")です"<<endl;
				//for(int D=0;D<3;D++) P_grad[D][i]=0;
			}
		}
	}

	/*/速度修正量計算
	double *density=new double[fluid_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) reU[D][i]=-dt*P_grad[D][i]/density[i];
	delete [] density;*/

	///////////////////////////////圧力勾配を表示
	if(CON->get_iteration_count()==1)				//陰解析を複数回行う場合は、毎回出力してしまうのでやめる
	{
		plot_Pgradient(CON,PART,fluid_number,P_grad);
	}
	////////////////////////*/

	delete [] matrix;
	delete [] B;
	delete [] minP;
}

///圧力勾配計算関数ver.5 自身を通ると仮定しないWLSM
void P_gradient5(mpsconfig *CON,vector<mpsparticle> &PART,double dt,int fluid_number,double *reU[3],double *P_grad[3])
{
	//P_gradient4()より未知数が１つ多い場合　近似局面が自身の値を通ると仮定しない。
	//その都合上、minPは考慮しない。

	double le=CON->get_distancebp();
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N0=0;							//係数行列の元
	int order0=CON->get_Pgrad_order();	//近似曲面のｵｰﾀﾞｰ。 1=線形 2=二次
	double threshold=CON->get_threshold();	//圧力誤差の閾値(アダプティブをする際に使用)

	//係数行列の大きさの決定
	if(d==2)
	{
		if(order0==1) N0=3;
		else if(order0==2) N0=6;
		else if(order0==3) N0=10;
	}
	else if(d==3)
	{
		if(order0==1) N0=4;
		else if(order0==2) N0=10;
		else if(order0==3) N0=20;
	}
	////////////////////////////////

	double *matrix=new double [N0*N0];	//N×Nの係数行列
	double *B=new double [N0];	//Nの解行列
	double **MAT= new double*[N0];		//N×Nの係数行列(配列は2次元)
	for(int n=0;n<N0;n++) MAT[n]=new double[N0];
	double *base=new double[N0];			//基底ベクトル格納
	ofstream fg("error_in_P.dat");
	
	int N=N0;
	int order=order0;
	if(d==2 )//二次元
	{
		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			int jnb=0;		//周辺粒子数。ただしOUTWALLなどは除く
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL) jnb++;
			}

			order=order0;
			if(order==3) {if(jnb<10) {order=2; N=6;}}	//周辺粒子数が少なかったら2次近似にする
			if(order==2) {if(jnb<6) {order=1; N=3;}}	//周辺粒子数が少なかったら1次近似にする

			if(jnb>=2)		//周辺粒子数がこれより少ない時は、圧力勾配はゼロ
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;//初期化
				for(int n=0;n<N0;n++) B[n]=0;
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//初期化

				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type!=OUTWALL)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// leで割るのは打ち切り誤差防止
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
						double dP=(PART[j].P)/PART[i].L;
						double dis=sqrt(X*X+Y*Y);
					
						double w=1;
						if(dis>1) w=1/(dis*dis*dis*dis);
						//if(dis>1) w=kernel_in_WLSM( dis, CON->get_re());

						if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//基底ベクトル
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y; base[9]=1;}

						//行列作成
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B[n]+=base[n]*w*dP;
						}
					}
				}

				MAT[N-1][N-1]+=1;		//一番右下の配列に自分自身の寄与を加える
				B[N-1]+=PART[i].P/PART[i].L;	//粒子iの寄与。
				//これ以外はX=0になるので加算する必要はない
				
				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//値をmatrixに転送

				gauss(matrix,B,N);//ガウスの消去法で解く
				
				P_grad[A_X][i]=B[0];	//1次近似でも2次近似でも、未知数の順番的にこうなるようにしてある
				P_grad[A_Y][i]=B[1];
				
				int flag=OFF;//OFFなら問題ない。ONならMPSでやり直す
				for(int D=0;D<2;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//非実数ならON
				if(flag==ON)
				{
					cout<<"エラーのためMPSにより離散化 "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					double grad[3]={0,0,0};
					double n0=1;//ここでは無意味な変数なので適当にいれとく
					P_gradient_MPS(CON,PART,i,n0,grad);//
					for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
				}
				/*double Q=0;
				int count_for_Q=0;
				if(flag==OFF)
				{
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type!=OUTWALL)
						{
							double X=(PART[j].r[A_X]-PART[i].r[A_X]);
							double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
							double dis=sqrt(X*X+Y*Y);
					
							double w=1;
							//if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
							if(dis>PART[i].L) w=kernel_in_WLSM( dis, R);
							Q+=(a*X+b*Y+PART[i].P-PART[j].P)*(a*X+b*Y+PART[i].P-PART[j].P)*w;
							count_for_Q++;
						}
					}
					if(count_for_Q!=0) Q/=count_for_Q;
					double P2=PART[i].P*PART[i].P;
					if(P2>1e-4) Q/=P2;
					//if(Q>threshold && PART[i].surface==OFF) PART[i].division_flag=1;		//誤差が大きいと判断して、分割する。
					if(Q>threshold) PART[i].division_flag=1;		//誤差が大きいと判断して、分割する。
				}
				//fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				//if(PART[i].division_flag==1) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				if(Q>threshold) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;*/
			}
		}
	}
	else if(d==3 )//3次元1次式
	{
		///係数行列は
		///   ΣΔx2    ΣΔxΔy  ΣΔxΔz  a = ΣΔxΔf  
		///  ΣΔxΔy    ΣΔy2   ΣΔyΔz  b = ΣΔyΔf 
		///  ΣΔxΔz   ΣΔyΔz   ΣΔz2   c = ΣΔzΔf 

		//２次近似なら
		//P=Pi+aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+iΔzΔxとおくと、
		///係数行列は
		///   ΣΔx2      ΣΔxΔy    ΣΔxΔz    ΣΔx3       ΣΔxΔy2    ΣΔxΔz2    ΣΔx2Δy     ΣΔxΔyΔz  ΣΔx2Δz     a = ΣΔxΔP  
		///   ΣΔxΔy    ΣΔy2      ΣΔyΔz    ΣΔx2Δy    ΣΔy3       ΣΔyΔz2    ΣΔxΔy2     ΣΔy2Δz    ΣΔxΔyΔz   b = ΣΔyΔP
		///   ΣΔxΔz    ΣΔyΔz    ΣΔz2      ΣΔx2Δz    ΣΔy2Δz    ΣΔz3       ΣΔxΔyΔz   ΣΔyΔz2    ΣΔxΔz2     c = ΣΔzΔP
		///   ΣΔx3      ΣΔx2Δy   ΣΔx2Δz   ΣΔx4       ΣΔx2Δy2   ΣΔx2Δz2   ΣΔx3Δy     ΣΔx2ΔyΔz ΣΔx3Δz     d = ΣΔx2ΔP
		///   ΣΔxΔy2   ΣΔy3      ΣΔy2Δz   ΣΔx2Δy2   ΣΔy4       ΣΔy2Δz2   ΣΔxΔy3     ΣΔy3Δz    ΣΔxΔy2Δz  e = ΣΔy2ΔP
		///   ΣΔxΔz2   ΣΔyΔz2   ΣΔz3      ΣΔx2Δz2   ΣΔy2Δz2   ΣΔz4       ΣΔxΔyΔz2  ΣΔyΔz3    ΣΔxΔz3     f = ΣΔz2ΔP
		///   ΣΔx2Δy   ΣΔxΔy2   ΣΔxΔyΔz ΣΔx3Δy    ΣΔxΔy3    ΣΔxΔyΔz2 ΣΔx2Δy2    ΣΔxΔy2Δz ΣΔx2ΔyΔz  g = ΣΔxΔyΔP
		///   ΣΔxΔyΔz ΣΔy2Δz   ΣΔyΔz2   ΣΔx2ΔyΔz ΣΔy3Δz    ΣΔyΔz3    ΣΔxΔy2Δz  ΣΔy2Δz2   ΣΔxΔyΔz2  h = ΣΔyΔzΔP
		///   ΣΔx2Δz   ΣΔxΔyΔz ΣΔxΔz2   ΣΔx3Δz    ΣΔxΔy2Δz  ΣΔxΔz3   ΣΔx2Δy Δz ΣΔxΔyΔz2 ΣΔx2Δz2    g = ΣΔxΔzΔP

		for(int i=0;i<fluid_number;i++)
		{
			for(int D=0;D<DIMENTION;D++) P_grad[D][i]=0;
			int jnb=0;		//周辺粒子数。ただしOUTWALLなどは除く
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL) jnb++;
			}

			order=order0;
			if(order==3) if(jnb<20) {order=2; N=10;}	//周辺粒子数が少なかったら2次近似にする
			if(order==2) if(jnb<10) {order=1; N=4;}	//周辺粒子数が少なかったら1次近似にする

			if(jnb>=3)
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;//初期化
				for(int n=0;n<N0;n++) B[n]=0;
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//初期化
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type!=OUTWALL)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
						double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double dP=(PART[j].P);//Δfjに相当
					
						double w=1;
						if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
						//if(dis>PART[i].L) w=kernel_in_WLSM( dis, R);
						
						if(order==1) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=1;}	//基底ベクトル
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=X*X*X; base[10]=Y*Y*Y; base[11]=Z*Z*Z; base[12]=X*X*Y; base[13]=X*Y*Y; base[14]=Y*Y*Z; base[15]=Y*Z*Z; base[16]=X*X*Z; base[17]=X*Z*Z; base[18]=X*Y*Z; base[19]=1;}
				
						//行列作成
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B[n]+=base[n]*w*dP;
						}
					}
				}
				MAT[N-1][N-1]+=1;		//一番右下の配列に自分自身の寄与を加える
				B[N-1]+=PART[i].P;	//粒子iの寄与。
				
				//これ以外はX=0になるので加算する必要はない*/

				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//値をmatrixに転送
				
				gauss(matrix,B,N);//ガウスの消去法で解く
				double Px=B[0];//X方向微分
				double Py=B[1];//Y方向微分
				double Pz=B[2];//Z方向微分

				P_grad[A_X][i]=Px;
				P_grad[A_Y][i]=Py;
				P_grad[A_Z][i]=Pz;
				
			}
			else//周辺粒子の重み平均で. 粒子数が少ないので、MPSはちがうかと
			{
				double W=0;
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					double X[3];
					for(int D=0;D<3;D++) X[D]=PART[j].r[D]-PART[i].r[D];
					double dis=sqrt(X[A_X]*X[A_X]+X[A_Y]*X[A_Y]+X[A_Z]*X[A_Z]);
					double dP=PART[j].P-PART[i].P;
					double w=kernel(R,dis);
					W+=w;
					for(int D=0;D<3;D++) P_grad[D][i]+=dP*X[D]/dis*w;
				}
				if(W!=0) for(int D=0;D<3;D++) P_grad[D][i]/=W;
				cout<<"警告 近隣粒子数が4以下("<<PART[i].N<<")です"<<endl;
			}

			int flag=OFF;//OFFなら問題ない。ONならMPSでやり直す
			for(int D=0;D<3;D++) if(P_grad[D][i]*2!=P_grad[D][i]+P_grad[D][i]) flag=ON;//非実数ならON
			if(flag==ON)
			{
				cout<<"エラーのためMPSにより離散化 "<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				double grad[3]={0,0,0};
				double n0=1;//ここでは無意味な変数なので適当にいれとく
				P_gradient_MPS(CON,PART,i,n0,grad);//
				for(int D=0;D<DIMENTION;D++) P_grad[D][i]=grad[D];
			}
			double Q=0;
			int count_for_Q=0;
			if(order==1)
			{
				if(flag==OFF && PART[i].N>=3)
				{
					//double a,b,c,d,e;
					double a,b,c;
					a=B[0];//X方向微分
					b=B[1];//Y方向微分
					c=B[2];//Z方向微分
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type!=OUTWALL)
						{
							double X=(PART[j].r[A_X]-PART[i].r[A_X]);
							double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
					
							double w=1;
							//if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
							if(dis>PART[i].L) w=kernel_in_WLSM( dis, R);
							double error=a*X+b*Y+c*Z+PART[i].P-PART[j].P;
							Q+=error*error*w;
							count_for_Q++;
						}
					}
					if(count_for_Q!=0) Q/=count_for_Q;
					double P2=PART[i].P*PART[i].P;
					if(P2>1e-4) Q/=P2;
				
					//if(Q>threshold && PART[i].surface==OFF) PART[i].division_flag=1;		//誤差が大きいと判断して、分割する。
					if(Q>threshold) 
					{
						PART[i].division_flag=1;		//誤差が大きいと判断して、分割する。
						//cout<<"圧力誤差が大きいため分割(threshold)"<<endl;
					}
				}
				//if(PART[i].division_flag==1) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				if(fabs(PART[i].r[A_Y])<le) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Q<<endl;
			}
			else if(order==2)
			{
				if(flag==OFF)
				{
					double a,b,c,d,e,f,g,h,ii;
					a=B[0];//X方向微分
					b=B[1];//Y方向微分
					c=B[2];//Z方向微分
					d=B[3];//X方向2階微分
					e=B[4];//Y方向2階微分
					f=B[5];//Z方向2階微分
					g=B[6];//XY方向微分
					h=B[7];//YZ方向微分
					ii=B[8];//ZX方向微分
					
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(PART[j].type!=OUTWALL)
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
					
							double w=1;
							if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
							double error=a*X+b*Y+c*Z+d*X*X+e*Y*Y+f*Z*Z+g*X*Y+h*Y*Z+ii*X*Z+PART[i].P-PART[j].P;
							Q+=error*error*w;
							count_for_Q++;
						}
					}
					if(count_for_Q!=0) Q/=count_for_Q;
					double P2=PART[i].P*PART[i].P;
					if(P2>1e-4) Q/=P2;
				
					//if(Q>threshold && PART[i].surface==OFF) PART[i].division_flag=1;		//誤差が大きいと判断して、分割する。	
					if(Q>threshold)
					{
						PART[i].division_flag=1;		//誤差が大きいと判断して、分割する。
						//cout<<"圧力誤差が大きいため分割(threshold)"<<endl;
					}
				}
				//if(PART[i].division_flag==1) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<Q<<endl;
				if(fabs(PART[i].r[A_Y])<le) fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<Q<<endl;
			}
		}
	}

	///////////////////////////////圧力勾配を表示
	if(CON->get_iteration_count()==1)				//陰解析を複数回行う場合は、毎回出力してしまうのでやめる
	{
		plot_Pgradient(CON,PART,fluid_number,P_grad);
	}
	////////////////////////*/
	fg.close();
	delete [] matrix;
	delete [] B;
	for(int n=0;n<N0;n++) delete [] MAT[n];
    delete [] MAT;
	delete [] base;
}


//Pgradient4における、3次元1次近似を行う関数
void cacl_WLSM_P_D3_order1(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP)
{
	double le=CON->get_distancebp();
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		if(PART[j].type!=OUTWALL)
		{
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double dP=(PART[j].P-minP[i]);//Δfjに相当
					
			double w=1;
			if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
						
			matrix[0]+=X*X*w;			//ΣXj^2wj
			matrix[1]+=X*Y*w;		//ΣXjYjwj
			matrix[2]+=X*Z*w;		//ΣXjZjwj
						
			matrix[4]+=Y*Y*w;			//ΣYj^2wj
			matrix[5]+=Y*Z*w;		//ΣYjZjwj
	
			matrix[8]+=Z*Z*w;			//ΣZj^2wj
					
			B[0]+=dP*X*w;//ΣfjXjwj
			B[1]+=dP*Y*w;//ΣfjYjwj
			B[2]+=dP*Z*w;//ΣfjZjwj
		}
	}
	matrix[3]=matrix[1];		//ΣXjYjwj
	matrix[6]=matrix[2];		//ΣXjZjwj
	matrix[7]=matrix[5];		//ΣYjZjwj

	//計算するかしないかを判定
	if(B[0]==0 && B[1]==0 && B[2]==0)
	{
		//解行列がすべてゼロのときは解けなくなるのでガウスの消去法にもちこんではいけない。
		for(int D=0;D<3;D++) P_grad[D][i]=0;//解行列がゼロということは勾配がゼロということだからそれを代入
	}
	else
	{
		gauss(matrix,B,N);//ガウスの消去法で解く
		double Px=B[0];//X方向微分
		double Py=B[1];//Y方向微分
		double Pz=B[2];//Z方向微分

		P_grad[A_X][i]=Px;
		P_grad[A_Y][i]=Py;
		P_grad[A_Z][i]=Pz;
	}
	
	///*/

	/*double determinant=(matrix[0]*matrix[4]*matrix[8]-matrix[0]*matrix[5]*matrix[7]-matrix[1]*matrix[3]*matrix[8]+matrix[1]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*matrix[7]-matrix[2]*matrix[4]*matrix[6]);//行列式
	//if(i==48155)cout<<"D="<<determinant<<endl;		
	P_grad[A_X][i]=(B[0]*matrix[4]*matrix[8]-B[0]*matrix[5]*matrix[7]-matrix[1]*B[1]*matrix[8]+matrix[1]*matrix[5]*B[2]+matrix[2]*B[1]*matrix[7]-matrix[2]*matrix[4]*B[2])/determinant;
	P_grad[A_Y][i]=(matrix[0]*B[1]*matrix[8]-matrix[0]*matrix[5]*B[2]-B[0]*matrix[3]*matrix[8]+B[0]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*B[2]-matrix[2]*B[1]*matrix[6])/determinant;
	P_grad[A_Z][i]=(matrix[0]*matrix[4]*B[2]-matrix[0]*B[1]*matrix[7]-matrix[1]*matrix[3]*B[2]+matrix[1]*B[1]*matrix[6]+B[0]*matrix[3]*matrix[7]-B[0]*matrix[4]*matrix[6])/determinant;	
	///*/
}

//Pgradient4における、3次元2次近似を行う関数
void cacl_WLSM_P_D3_order2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP)
{
	//P=Pi+aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+iΔzΔxとおくと、
		///係数行列は
		///   ΣΔx2      ΣΔxΔy    ΣΔxΔz    ΣΔx3       ΣΔxΔy2    ΣΔxΔz2    ΣΔx2Δy     ΣΔxΔyΔz  ΣΔx2Δz     a = ΣΔxΔP  
		///   ΣΔxΔy    ΣΔy2      ΣΔyΔz    ΣΔx2Δy    ΣΔy3       ΣΔyΔz2    ΣΔxΔy2     ΣΔy2Δz    ΣΔxΔyΔz   b = ΣΔyΔP
		///   ΣΔxΔz    ΣΔyΔz    ΣΔz2      ΣΔx2Δz    ΣΔy2Δz    ΣΔz3       ΣΔxΔyΔz   ΣΔyΔz2    ΣΔxΔz2     c = ΣΔzΔP
		///   ΣΔx3      ΣΔx2Δy   ΣΔx2Δz   ΣΔx4       ΣΔx2Δy2   ΣΔx2Δz2   ΣΔx3Δy     ΣΔx2ΔyΔz ΣΔx3Δz     d = ΣΔx2ΔP
		///   ΣΔxΔy2   ΣΔy3      ΣΔy2Δz   ΣΔx2Δy2   ΣΔy4       ΣΔy2Δz2   ΣΔxΔy3     ΣΔy3Δz    ΣΔxΔy2Δz  e = ΣΔy2ΔP
		///   ΣΔxΔz2   ΣΔyΔz2   ΣΔz3      ΣΔx2Δz2   ΣΔy2Δz2   ΣΔz4       ΣΔxΔyΔz2  ΣΔyΔz3    ΣΔxΔz3     f = ΣΔz2ΔP
		///   ΣΔx2Δy   ΣΔxΔy2   ΣΔxΔyΔz ΣΔx3Δy    ΣΔxΔy3    ΣΔxΔyΔz2 ΣΔx2Δy2    ΣΔxΔy2Δz ΣΔx2ΔyΔz  g = ΣΔxΔyΔP
		///   ΣΔxΔyΔz ΣΔy2Δz   ΣΔyΔz2   ΣΔx2ΔyΔz ΣΔy3Δz    ΣΔyΔz3    ΣΔxΔy2Δz  ΣΔy2Δz2   ΣΔxΔyΔz2  h = ΣΔyΔzΔP
		///   ΣΔx2Δz   ΣΔxΔyΔz ΣΔxΔz2   ΣΔx3Δz    ΣΔxΔy2Δz  ΣΔxΔz3   ΣΔx2Δy Δz ΣΔxΔyΔz2 ΣΔx2Δz2    g = ΣΔxΔzΔP

	double le=CON->get_distancebp();
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。
	double B_val[3];						//解行列を3要素分だけ保存しておく。あとで1次近似でやりなおす場合に必要
	double matrix_val[9];					//係数行列を9要素分だけ保存しておく。あとで1次近似でやりなおす場合に必要

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		if(PART[j].type!=OUTWALL)
		{
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double dP=(PART[j].P-minP[i]);//Δfjに相当
						
			double w=1;
			//if(dis>le) w=le*le/(dis*dis);
			//if(dis>le) w=le/(dis);
			if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
			weight[k]=w;

			matrix[0]+=X*X*w;		//ΣXj^2wj
			matrix[1]+=X*Y*w;		//ΣXjYjwj
			matrix[2]+=X*Z*w;		//ΣXjZjwj
			matrix[3]+=X*X*X*w;		//ΣXj^3wj
			matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
			matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
			matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
			matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
			matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj
	
			matrix[10]+=Y*Y*w;		
			matrix[11]+=Y*Z*w;		
			matrix[12]+=X*X*Y*w;		
			matrix[13]+=Y*Y*Y*w;		
			matrix[14]+=Y*Z*Z*w;
			matrix[15]+=X*Y*Y*w;
			matrix[16]+=Y*Y*Z*w;
						
			matrix[20]+=Z*Z*w;			
			matrix[23]+=Z*Z*Z*w;		
					
			matrix[30]+=X*X*X*X*w;
			matrix[31]+=X*X*Y*Y*w;
			matrix[32]+=X*X*Z*Z*w;	
			matrix[33]+=X*X*X*Y*w;	
			matrix[34]+=X*X*Y*Z*w;	
			matrix[35]+=X*X*X*Z*w;	
					
			matrix[40]+=Y*Y*Y*Y*w;
			matrix[41]+=Y*Y*Z*Z*w;
			matrix[42]+=X*Y*Y*Y*w;
			matrix[43]+=Y*Y*Y*Z*w;
			matrix[44]+=X*Y*Y*Z*w;

			matrix[50]+=Z*Z*Z*Z*w;	//6行目
			matrix[51]+=X*Y*Z*Z*w;
			matrix[52]+=Y*Z*Z*Z*w;
			matrix[53]+=X*Z*Z*Z*w;

			//7～9行目はすべて既存の要素から転用が可能


			B[0]+=dP*X*w;		//a
			B[1]+=dP*Y*w;		//b
			B[2]+=dP*Z*w;		//c
			B[3]+=dP*X*X*w;		//d
			B[4]+=dP*Y*Y*w;		//e
			B[5]+=dP*Z*Z*w;		//f
			B[6]+=dP*X*Y*w;		//g
			B[7]+=dP*Y*Z*w;		//h
			B[8]+=dP*X*Z*w;		//i
		}
	}
	matrix[9]=matrix[1];		//ΣXjYjwj
	matrix[17]=matrix[7];

	matrix[18]=matrix[2];
	matrix[19]=matrix[11];
	matrix[21]=matrix[8];
	matrix[22]=matrix[16];
	matrix[24]=matrix[7];
	matrix[25]=matrix[14];
	matrix[26]=matrix[5];

	matrix[27]=matrix[3];
	matrix[28]=matrix[12];
	matrix[29]=matrix[21];

	matrix[36]=matrix[4];
	matrix[37]=matrix[13];
	matrix[38]=matrix[22];
	matrix[39]=matrix[31];

	matrix[45]=matrix[5];
	matrix[46]=matrix[14];
	matrix[47]=matrix[23];
	matrix[48]=matrix[32];
	matrix[49]=matrix[41];

	matrix[54]=matrix[6];
	matrix[55]=matrix[15];
	matrix[56]=matrix[24];
	matrix[57]=matrix[33];
	matrix[58]=matrix[42];
	matrix[59]=matrix[51];
	matrix[60]=matrix[31];
	matrix[61]=matrix[44];
	matrix[62]=matrix[34];

	matrix[63]=matrix[7];
	matrix[64]=matrix[16];
	matrix[65]=matrix[25];
	matrix[66]=matrix[34];
	matrix[67]=matrix[43];
	matrix[68]=matrix[52];
	matrix[69]=matrix[61];
	matrix[70]=matrix[41];
	matrix[71]=matrix[51];

	matrix[72]=matrix[8];
	matrix[73]=matrix[17];
	matrix[74]=matrix[26];
	matrix[75]=matrix[35];
	matrix[76]=matrix[44];
	matrix[77]=matrix[53];
	matrix[78]=matrix[62];
	matrix[79]=matrix[71];
	matrix[80]=matrix[32];
		
	for(int L=0;L<3;L++) B_val[L]=B[L];//解行列の値を保存
	for(int L=0;L<3;L++)
	{
		for(int M=0;M<3;M++)
		{
			matrix_val[L*3+M]=matrix[L*9+M];//解行列の値を保存
		}
	}
		

	//行列をガウスの消去法で解く　解はBに格納される
	gauss(matrix,B,N);

	//誤差を調査
	double Q=0;
	double W=0;//重みの総和
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		if(PART[j].type!=OUTWALL)
		{
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
			double dP=(PART[j].P-minP[i]);//Δfjに相当
			double w=weight[k];
			W+=w;
			double err=B[0]*X+B[1]*Y+B[2]*Z+B[3]*X*X+B[4]*Y*Y+B[5]*Z*Z+B[6]*X*Y+B[7]*Y*Z+B[8]*X*Z-dP;
			Q+=err*err*w;
		}
	}
	if(W>0) Q/=W;

	if(Q>100)//誤差が100以上なら(ここ見直すこと)
	{
		//cout<<i<<endl;
		for(int n=0;n<3;n++) B[n]=B_val[n];	//Px=B[0]のように、B[n]に答えを格納させるようにしてるので、ここではBを使う
		N=3;								//N=3にしておかないといけない
		gauss(matrix_val,B,N);				//matrix_valには今回使用した行列のうち[3][3]だけが１次元に格納されている。
	}///*/
	
	double Px=B[0];//X方向微分
	double Py=B[1];//Y方向微分
	double Pz=B[2];//Z方向微分
	
	P_grad[A_X][i]=Px;
	P_grad[A_Y][i]=Py;
	P_grad[A_Z][i]=Pz;/////////*/

	delete [] weight;
}

///圧力勾配計算関数
void P_gradient_MPS(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0,double *P_grad)
{
    double R=CON->get_re()*CON->get_distancebp();
	double W=0;//粒子数密度　OUTを除いたりするのでPND[i]は微妙
	double minP=PART[i].P;//周囲の最少圧力
	int d=CON->get_dimention();
	int J=i;
	for(int k=0;k<PART[i].N;k++)
	{       
	    int j=PART[i].NEI[k];
	    if(PART[j].type!=OUTWALL)
	    {
	        if(PART[j].P<minP)
			{
				minP=PART[j].P;
				J=j;
			}
	    }
	}///minP[i]が求まった
    //minP=PART[i].P;
	////////
	for(int k=0;k<PART[i].N;k++)
	{       
	    int j=PART[i].NEI[k];
	    if(PART[j].type!=OUTWALL)
	    {
	        double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
		
			double w=kernel(R,dis);
			W+=w;
		
			P_grad[A_X]+=(PART[j].P-minP)*X*w/(dis*dis);
			P_grad[A_Y]+=(PART[j].P-minP)*Y*w/(dis*dis);
			P_grad[A_Z]+=(PART[j].P-minP)*Z*w/(dis*dis);
		
		
			//P_grad[A_X]+=(PART[i].P-minP)*X*w/(dis*dis)/2;
			//P_grad[A_Y]+=(PART[i].P-minP)*Y*w/(dis*dis)/2;
			//P_grad[A_Z]+=(PART[i].P-minP)*Z*w/(dis*dis)/2;
		
	    }
	}////*/
	
	for(int D=0;D<DIMENTION;D++) if(W!=0) P_grad[D]= P_grad[D]*d/W;
	//for(int D=0;D<DIMENTION;D++) if(W!=0) P_grad[D]= P_grad[D]*d/n0;
}

//圧力勾配による速度・位置修正関数
void modify_u_and_x_after_Pcalc(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **reU,double dt,double **Un)
{
	//速度修正による位置修正
	int sw=CON->get_temporary_r();		//ONなら陽解析後に仮の位置を計算している(通常の解析)
	int d=CON->get_dimention();

	//速度更新
	for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].u[D]+=reU[D][i];

	//表面粒子の速度を、上のように求めるのではなく、内部粒子からの近似で求める場合
	if(CON->get_interpolate_surface_u()==ON) reset_surface_velocity_by_WLSM(CON,PART, fluid_number);


	//表面だったら速度の更新をやめる
	//if(CON->get_fix_surface()==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].surface==ON) PART[i].u[D]=0;

	//位置更新
	//if(CON->get_Position_by_WLSM()==ON) update_position_by_WLSM(CON,PART, fluid_number,dt,Un);//最小二乗法による位置更新
	//if(CON->get_Position_by_WLSM()==ON) update_position_by_WLSM_2(CON,PART, fluid_number,dt,Un);//最小二乗法による位置更新
	//else
	{
		if(sw==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*reU[D][i]*0.5;//台形則
		else if(sw==OFF) 
		{
			//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*PART[i].u[D];
			for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*(PART[i].u[D]+Un[D][i])*0.5;//台形則
		}///*/
	}

	//表面だったら位置の更新をやめる
	//if(CON->get_fix_surface()==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].surface==ON) PART[i].r[D]-=dt*reU[D][i]*0.5;//台形則
}

//陰解析ver,3
void negative3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double lamda,double N0,vector<point3D> &NODE,vector<element3D> &ELEM,int out)
{
	cout<<"FEMによる陰解析実行"<<endl;

	if(CON->get_dimention()==2) cout<<"2次元は非対応"<<endl;

	point3D NODE0;
	element3D ELEM0;
	NODE.push_back(NODE0);//0番目の要素だけ生成しておく
	ELEM.push_back(ELEM0);//0番目の要素だけ生成しておく

	//INWALLの数を数える
	int inwall_number=0;
	for(int i=fluid_number;i<particle_number;i++) if(PART[i].type==INWALL && PART[i].surface==OFF) inwall_number++;//BDWALLは計算には組み込まない

	int node=fluid_number+inwall_number;		//節点数=流体粒子数
    int KTJ=node;							//最大節点数
    int KTE=12*KTJ;								//最大要素数　3次元式もとめよ
	int nelm=0;									//現在の要素数
    double err=1.0e-14;							//誤差判定のしきい値
	double femfactor=1;//1000;							//打ち切り誤差対策係数
	
	///ﾃﾞﾛｰﾆを行うか否かの判定
	int delaun_flag;							//ﾃﾞﾛｰﾆ分割を行うか行わないか
	if(t==1 || t%100==0) delaun_flag=ON;
	else delaun_flag=OFF;
	///////////////////////
	
	if(delaun_flag==OFF)
	{
		nelm=NODE[0].material;
		node=NODE[0].particleID;
		//nelm=(int)ELEM.size();
		//node=(int)NODE.size()-(8+1);//スーパーボックのぶんと0番目の要素は受け取ってはいけない
	}


    if(delaun_flag==OFF)	//ﾃﾞﾛｰﾆはしなくても節点の座標は変更する
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;//節点iは粒子p
			for(int D=0;D<3;D++) NODE[i].r[D]=PART[p].r[D];
		}
	}
	else if(delaun_flag==ON)
	{
		/////////////input
		int num=1;
		for(int i=0;i<fluid_number;i++)
		{
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].material=FLUID;
			NODE[num].particleID=i;									//対応する粒子番号格納
			if(PART[i].surface==ON) NODE[num].boundary_condition=1;	//固定境界条件
			else NODE[num].boundary_condition=0;
			num++;
		}
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL && PART[i].surface==OFF)//流体と接している壁なら
			{
				NODE.push_back(NODE0);
				for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
				NODE[num].material=FLUID;		//面倒なのでここではWATERと定義
				NODE[num].particleID=i;			//対応する粒子番号格納
				NODE[num].boundary_condition=0;	//自然境界条件
				num++;
			}
		}//////////////*/
    
		cout<<"節点数="<<node<<"  最大節点数="<<KTJ<<"   最大要素数="<<KTE<<endl;

		///femfacterによる拡大
		for(int i=1;i<=node;i++) for(int D=0;D<3;D++) NODE[i].r[D]*=femfactor;

		for(int i=KTJ+8;i>node;i--) NODE.push_back(NODE0);//スーパーボック用の節点を確保
		for(int i=1;i<=KTE;i++) ELEM.push_back(ELEM0);
	
	    /////////////節点座標の正規化
		double xmin=NODE[1].r[A_X];
		double ymin=NODE[1].r[A_Y];
		double zmin=NODE[1].r[A_Z];
		double xmax=xmin;
		double ymax=ymin;
		double zmax=zmin;

		///座標の最大、最小値を求める
		for(int i=2;i<=node;i++)
		{
			if(NODE[i].r[A_X]<xmin) xmin=NODE[i].r[A_X];
			else if(NODE[i].r[A_X]>xmax) xmax=NODE[i].r[A_X];
	
			if(NODE[i].r[A_Y]<ymin) ymin=NODE[i].r[A_Y];
			else if(NODE[i].r[A_Y]>ymax) ymax=NODE[i].r[A_Y];
		
			if(NODE[i].r[A_Z]<zmin) zmin=NODE[i].r[A_Z];
			else if(NODE[i].r[A_Z]>zmax) zmax=NODE[i].r[A_Z];
		}
		////

		double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X軸方向の寸法
		double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y軸方向の寸法
		double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z軸方向の寸法
		double rmax=rax;							///最大寸法
		if(ray>rmax) rmax=ray;
		if(raz>rmax) rmax=raz;						//ここはelseにしたらダメ
	
		///座標変換
		double rrm=1.000000/rmax;///こういう書き方をすることで、数値誤差を減らせる・・？
		for(int i=1;i<=node;i++)
		{   //   A/Bという計算をしたとき、Ａの値によって微妙に1/Bという倍率がちがってくるのではないかと考えて、下のような書き方にしている
		    NODE[i].r[A_X]=(NODE[i].r[A_X]-xmin)*rrm;
			NODE[i].r[A_Y]=(NODE[i].r[A_Y]-ymin)*rrm;
			NODE[i].r[A_Z]=(NODE[i].r[A_Z]-zmin)*rrm;
		}
		rax*=rrm;
		ray*=rrm;
		raz*=rrm;
		/////

		///ﾃﾞﾛｰﾆ分割
		int FINE_sw=OFF;
		delaun3D(CON,NODE,ELEM,KTJ,KTE,rax,ray,raz,&node,&nelm,FINE_sw,rrm);
	
		////座標を元に戻す
		for(int i=1;i<=node;i++)
		{
			NODE[i].r[A_X]=rmax*NODE[i].r[A_X]+xmin;
			NODE[i].r[A_Y]=rmax*NODE[i].r[A_Y]+ymin;
			NODE[i].r[A_Z]=rmax*NODE[i].r[A_Z]+zmin;
		}
    
		for(int i=1;i<=node;i++) for(int D=0;D<3;D++) NODE[i].r[D]/=femfactor;
		/////ﾒｯｼｭ生成完了

		NODE[0].material=nelm;//要素数などを保存
		NODE[0].particleID=node;
	}
	/////////////////////////

	cout<<"要素数="<<nelm<<" 節点数＝"<<node<<endl;

	double *P=new double[node+1];
	for(int i=0;i<=node;i++) P[i]=1;

	 /////節点-要素関係
    int *jnb=new int[node+1];///各節点に隣接する要素数格納
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//各節点の周辺要素番号格納
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);
	
	//圧力解析のためにN0とPART[i].PND2をkernel2()用におきかえる
	set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	cout<<"N0="<<N0<<endl;
	///////////////////*/

	///////////有限要素法計算開始
	P3D(CON,NODE,ELEM,node,nelm,P,jnb,PART,fluid_number,nei,dt,N0);

	reU3D(CON,NODE,ELEM,node,nelm,P,PART,fluid_number,dt,jnb);

    delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
	delete [] P;

}

///圧力計算関数（FEM)
void P3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,int *jnb,vector<mpsparticle> &PART,int fluid_number,int **nei,double dt,double N0)
{
	cout<<"圧力計算開始--";
    double P0=0;					//表面圧力
	double le=CON->get_distancebp();//粒子間距離
	double density=CON->get_density();
	double CO=density/(dt*dt)/N0;		//計算に必要な係数
	unsigned timeA=GetTickCount();	//計算開始時刻

    int NN=0;					//ディリクレ型境界節点数
    int *dn=new int [node+1];	//各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
    double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
	
	double *Dirichlet_P=new double [fluid_number];//各粒子のディリクレ値

	if(CON->get_dir_for_P()==OFF) 
	{
		for(int i=0;i<fluid_number;i++) Dirichlet_P[i]=P0;
	}
	else //表面粒子の圧力として、ディリクレ値
	{
		int flag=CON->get_dir_for_P();
		
		
		if(CON->get_dir_for_P()==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
		else if(CON->get_dir_for_P()==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		else if(CON->get_dir_for_P()==3)
		{
			set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		}

		//ﾌｧｲﾙ出力
		ofstream gg("Dirichlet_P.dat");
		for(int i=0;i<fluid_number;i++) if(PART[i].surface==ON ) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;
		gg.close();///*/

	}///////////

    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			int PN=NODE[i].particleID;//MPSでの粒子番号
	        if(PN<fluid_number) PHAT[NN]=Dirichlet_P[PN];
			else PHAT[NN]=P0;
	        P[i]=PHAT[NN];
	        NN++;
		}
		else dn[i]=node+1;
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数="<<NN<<" ";
	    
    int pn=node-NN;///未知数
    int *ppn=new int [pn];		//行列でn番目の未知数は節点番号ppn[n]に相当
    int *npp=new int [node+1];	//i番目の節点は行列のnpp[i]番目に相当。ﾃﾞｨﾘｸﾚ型の場合は行列に入っていないのでpn+1を格納
    int num=0; 
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)	//未知数なら
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////行列の最大幅計算
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//メモリをきっちりもとめる。だけど若干の計算コストになる。
    //////
	
    ////配列確保
    double **G=new double *[pn+1];						///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1];							///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1];							///各行の、非ゼロ要素数
    
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    ////


    /////////全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double U[4+1];//各点の速度
    double V[4+1];
    double W[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	double dN[4+1];	//粒子数密度の差
    for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==FLUID)
		{
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				int particle_num=NODE[N[j]].particleID;//MPSにおける粒子番号
				U[j]=PART[particle_num].u[A_X];
				V[j]=PART[particle_num].u[A_Y];
				W[j]=PART[particle_num].u[A_Z];
				dN[j]=PART[particle_num].PND2-N0;
			}
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
			
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
		
			double delta=ELEM[je].volume/6;//本当の体積
			double co1=delta/20.0;			//係数

			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
			
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]))*delta6;//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]))*delta6;//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]))*delta6;//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
		
			double udiv=0;//速度発散
			udiv+=c[1]*U[1]+c[2]*U[2]+c[3]*U[3]+c[4]*U[4];
			udiv+=d[1]*V[1]+d[2]*V[2]+d[3]*V[3]+d[4]*V[4];
			udiv+=e[1]*W[1]+e[2]*W[2]+e[3]*W[3]+e[4]*W[4];
			udiv*=density/dt;
			////要素ﾏﾄﾘｸｽ作成開始
			for(int i=1;i<=4;i++)
			{
				if(NODE[N[i]].boundary_condition==0)///未知なら
				{   
					int I=npp[N[i]]+1;///節点N[i]は行列のI番目
					for(int j=1;j<=4;j++)
					{					
						int J=npp[N[j]]+1;///節点N[j]は行列のJ番目
						if(NODE[N[j]].boundary_condition==0)///未知なら
						{   
							int flag=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
									flag=1;
								}
							}
							if(flag==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
				    
								G[I][H]+=(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
								ROW[I][H]=J;
							}
						}
						else //N[j]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int n=dn[N[j]];
							B[I-1]-=(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta*PHAT[n];
						}
					}
					B[I-1]+=-udiv/4*delta;
					/*double A=0;				//粒子数密度の差を要素で体積積分した値
					for(int j=1;j<=4;j++)
					{
						if(j==i) A+=co1*2*dN[j];
						else A+=co1*dN[j];
					}
					B[I-1]+=CO*A;*/
				}
			}   
		}
	}
    ///////////////////////*/
    
	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//非ゼロ要素の列番号格納配列
    int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
    /////////////////////val,ind ,ptrに値を格納
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    
	///////////////////////行列計算開始
	double *XX=new double [pn];//行列の答え格納
	ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		P[i]=XX[n];
	}
	delete [] XX;
	////////////////////////////
    
    ///////圧力をﾌｧｲﾙ出力
	ofstream fp("P.dat");
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].r[A_Y]<0.5*le && NODE[i].r[A_Y]>-0.5*le) fp<<NODE[i].r[A_X]<<"\t"<<NODE[i].r[A_Z]<<"\t"<<P[i]<<endl;
    }
	fp.close(); 
    
    ////////////////////////

	delete [] Dirichlet_P;
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

}

///速度修正関数(FEM)
void reU3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *P,vector<mpsparticle> &PART,int fluid_number,double dt,int *jnb)
{
	unsigned timeA=GetTickCount();		//計算開始時刻
    
    int N[4+1];							//要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	double density=CON->get_density();
	double A=dt/density;				//よく使うので係数化

	double *newU[3];
	for(int D=0;D<3;D++) newU[D]=new double [node+1];	//各節点の新しい速度格納
	double *reU[3];
	for(int D=0;D<3;D++) reU[D]=new double [nelm+1];	//各要素の速度修正量格納

	int *num=new int[node+1];
	for(int i=1;i<=node;i++) num[i]=0;
	///newUに値を格納
	for(int i=1;i<=node;i++)
	{
		int particle_id=NODE[i].particleID;	//節点iはparticle_idの粒子
		for(int D=0;D<3;D++)  newU[D][i]=0;//初期化
	}////////////

    ///各節点の速度を更新する
    for(int je=1;je<=nelm;je++)
    {
		for(int D=0;D<3;D++) reU[D][je]=0;
		
        for(int j=1;j<=4;j++)
		{
			N[j]=ELEM[je].node[j];
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
	
		double delta6=ELEM[je].volume;//体積の6倍  体積は圧力求めるときに計算してある
	
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
		
		///係数c,d,e計算
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]))*delta6;//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]))*delta6;//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]))*delta6;//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}   
		}
		/////////
    
		reU[A_X][je]=c[1]*P[N[1]]+c[2]*P[N[2]]+c[3]*P[N[3]]+c[4]*P[N[4]];
        reU[A_Y][je]=d[1]*P[N[1]]+d[2]*P[N[2]]+d[3]*P[N[3]]+d[4]*P[N[4]];
		reU[A_Z][je]=e[1]*P[N[1]]+e[2]*P[N[2]]+e[3]*P[N[3]]+e[4]*P[N[4]];

		for(int j=1;j<=4;j++)
		{
			newU[A_X][N[j]]-=A*reU[A_X][je];
			newU[A_Y][N[j]]-=A*reU[A_Y][je];
			newU[A_Z][N[j]]-=A*reU[A_Z][je];
			
			num[N[j]]=num[N[j]]+1;//jnbと一致するのでは？
		}////*/

    }///速度更新終了
    
	//速度更新を粒子に変換
	for(int i=1;i<=node;i++)
	{
		int id=NODE[i].particleID;
		if(id<fluid_number)
		{	
			for(int D=0;D<3;D++)
			{
				if(num[i]<=1)cout<<"num<=1 ??"<<endl;
				PART[id].u[D]+=newU[D][i]/num[i];
				if(CON->get_temporary_r()==OFF) PART[id].r[D]+=PART[id].u[D]*dt;
				if(CON->get_temporary_r()==ON)  PART[id].r[D]+=newU[D][i]/num[i]*dt;
			}
		}
	}

	for(int D=0;D<3;D++) delete [] newU[D];
	for(int D=0;D<3;D++) delete [] reU[D];
	delete [] num;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//陰解析ver,2
void negativeP_iterative(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double lamda,double n0,double N0,int out,double **Un)
{
	cout<<"反復法による陰解析実行"<<endl;

	int d=CON->get_dimention();						//次元
	double le=CON->get_distancebp();				//初期粒子間距離
	double r2=CON->get_re2()*CON->get_distancebp();
	double *density=new double[particle_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=fluid_number;i<particle_number;i++) density[i]=CON->get_density();


	int count2=0;									//反復回数
	
	double error=0;
	unsigned timeA=GetTickCount();					//計算開始時刻


	int calctype=0;									//計算タイプ 0:速度のみ更新 1:位置も更新
	int pn=0;										//圧力を解くための連立方程式未知数
	
	double div_error=0;								//発散誤差
	int numf=0;										//内部流体粒子数
	int divflag=2;									//発散計算ﾌﾗｸﾞ 1:MPS 2:WLSM 3:Dndt
	int Btype=CON->get_B_of_P();					//解行列の種類
	
	if(Btype==3) divflag=3;

	for(int i=0;i<particle_number;i++) if(PART[i].surface==OFF) if(PART[i].type==FLUID || PART[i].type==INWALL) pn++;//方程式の数

	double *reU[DIMENTION];//速度修正量
	for(int D=0;D<DIMENTION;D++) reU[D] = new double [fluid_number];

	double *old_U[DIMENTION];//仮の速度記憶
	for(int D=0;D<DIMENTION;D++) old_U[D] = new double [fluid_number];

	double *old_r[DIMENTION];//仮の速度記憶
	for(int D=0;D<DIMENTION;D++) old_r[D] = new double [fluid_number];

	double *old_reU[DIMENTION];//1反復前の速度修正量記憶
	for(int D=0;D<DIMENTION;D++) old_reU[D] = new double [fluid_number];

	double *direct[DIMENTION];			//内向き法線ﾍﾞｸﾄﾙ格納
	for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];

	double *P_grad[DIMENTION];			//圧力勾配ﾍﾞｸﾄﾙ格納
	for(int D=0;D<DIMENTION;D++) P_grad[D]=new double [fluid_number];

	double *udiv=new double [particle_number];	//各粒子の速度発散(壁粒子も含む)

	//圧力解析のためにN0とPART[i].PND2をkernel2()用におきかえる
	set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	cout<<"N0="<<N0<<endl;
	///////////////////*/

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<d;D++)
		{
			old_U[D][i]=PART[i].u[D];//仮の速度を記憶
			old_reU[D][i]=0;//初期化
			old_r[D][i]=PART[i].r[D];
		}
	}

	//////法線ﾍﾞｸﾄﾙ計算開始
	if(CON->get_Pgrad()==2 ||CON->get_Pgrad()==3)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].type==BOFLUID ) direct_f(CON,PART,i,direct);
			else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
		}
		
	}/////////////////*/

	///各粒子の速度発散計算
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type!=OUTWALL)
		{
			if(divflag==1) udiv[i]=divergence(CON,PART,i,n0);
			else if(divflag==2)
			{
				udiv[i]=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*udiv[i]!=udiv[i]+udiv[i]) udiv[i]=divergence(CON,PART,i,n0);//udiv[i]がエラーで非実数のときはMPSで計算する
			}
			else if(divflag==3) udiv[i]=Dndt(CON,PART,i,n0);
		}
	}///////

	//最初の速度発散分布出力
	ofstream f0("div0.dat");
	if(d==2) for(int i=0;i<fluid_number;i++) f0<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<udiv[i]<<endl;
	if(d==3) for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) f0<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<udiv[i]<<endl;
	f0.close();////////////*/

	///速度の発散誤差を計算
	div_error=0;
	numf=0;//内部流体粒子数
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==OFF)
		{
			div_error+=udiv[i]*udiv[i];
			numf++;
		}
	}
	if(numf!=0) 
	{
		div_error=div_error/numf;
		div_error=sqrt(div_error);	//標準偏差
	}


	int *ppn = new int[pn];					//行列における第n番目の未知数は粒子番号ppn[n]の粒子に相当
	int *link = new int [particle_number];	//粒子番号iはlink[i]番目の未知数
	double *B   = new double[pn];			//解行列

	int count=0;
	///ppn配列とlink配列作成
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].surface==OFF)
		{
			if(PART[i].type==FLUID || PART[i].type==INWALL)
			{
				ppn[count]=i;
				link[i]=count;
				count++;
			}
			else link[i]=pn+1;//行列に含まれない粒子にはﾀﾞﾐｰとして(pn+1)を格納
		}
		else link[i]=pn+1;//行列に含まれない粒子にはﾀﾞﾐｰとして(pn+1)を格納
	}
	//////*/

	int number=0;			//係数行列の非ゼロ要素数
	for(int n=0;n<pn;n++)
	{   
		number++;///自分自身を数にいれる
		int i=ppn[n];//n番目の未知数は粒子i
	  
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k];
			int m=link[j];
			if(m<pn) number++;
		}
	}
	////numberが求まった

	double *val = new double [number];
	int *ind = new int [number];//非ゼロ要素の列番号格納配列
	int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
	
	
	/////////////////////val,ind ,ptrに値を格納
	double real_lamda=lamda;	//lamdaの値記憶
	if(CON->get_HL_sw()==ON) lamda=2.0*CON->get_dimention()*real_lamda;//高次の離散化の場合、lamda=2dと再定義すればλ/2dの項は消える。ただしそうすると両辺の値が大きくなりすぎるので、reak_lamdaを両辺にかけて小さくしている
		
	

	
	///////////////////////////////////反復開始
	int flag=0;
	while(flag==0)
	{
		unsigned int timeC=GetTickCount();	//計算開始時刻
	
		cout<<"圧力未知数:"<<pn<<" に対して";

		int remake_SW=OFF;			//係数行列を再度作成するか、しないか
		if(calctype==0) if(count2==0) remake_SW=ON;	//速度のみ修正する場合は最初だけ、
		if(calctype==1) remake_SW=ON;//位置も修正する場合は毎回作りなおす

		if(CON->get_HL_sw()==OFF && remake_SW==ON)//標準的な離散化
		{
			int index=0;
			for(int n=0;n<pn;n++)
			{
				ptr[n]=index;
				int i=ppn[n];
				int kk=index;//値を保存
				ind[index]=n;
				index++;
				double W=0;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					if(m<pn)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						    
						double w=kernel2(r2,dis,d);
						val[index]=w;
						ind[index]=m;
						index++;
						W+=w;
					}
					else if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel2(r2,dis,d);
						W+=w; //ここでwを計算に組み込むとき、粒子jはﾃﾞｨﾘｸﾚ型、組み込まないと勾配ゼロのﾉｲﾏﾝ型 
					}
				}
				val[kk]=-W;
			}
			ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない
		}
		////////////////////*/ 

		//if(count2==CON->get_iteration_count()) Btype=4;		//一番最後の反復だけは密度型

		///////////////////////////////////////////解行列Bの計算
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];//行列のn番目の粒子はi
			if(Btype==0) B[n]=-density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());//教科書
			else if(Btype==1)//速度発散
			{    
				double div=udiv[i];
				B[n]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i];
			}
			else if(Btype==2)//(教科書+速度発散)/2
			{
				double a=CON->get_w_div();double b=1;///各手法の重み
				double div=udiv[i];
				
				B[n]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i]*a;
				B[n]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;
				B[n]/=a+b;
			}
			else if(Btype==3)//重み関数の直接微分による速度発散
			{    
				double div=udiv[i];
				B[n]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
			}
			else if(Btype==4)//重み関数の直接微分による速度発散+PND
			{    
				double a=CON->get_w_div();double b=1;///各手法の重み
				double div=udiv[i];
				B[n]=-div*lamda/(dt*2*CON->get_dimention())*density[i]*a;
				B[n]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;//教科書
				B[n]/=a+b;
			}
			else if(Btype==5)//(ni-nk)+(nk-n0)=div+pnd(nk-n0)
			{    
				double div=udiv[i];
				B[n]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
				B[n]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());
			}
			else cout<<"解行列未解決 check B_of_P"<<endl;
		}
		////解行列Bがもとまった//*/
	
	
		if(CON->get_dir_for_P()!=OFF)//表面粒子の圧力として、ディリクレ値
		{
			int flag=CON->get_dir_for_P();
			if(flag==2 || flag==3)
			{
				if(CON->get_EM_method()==OFF) flag=1;//EMによる計算を行わないのにflagが2や3なら間違いなので1に戻す
			}
			double *Dirichlet_P=new double [particle_number];
			for(int i=0;i<particle_number;i++) Dirichlet_P[i]=0;

			if(flag==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			else if(flag==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
			else if(flag==3)
			{
				set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
				set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
			}

			for(int i=fluid_number;i<particle_number;i++)//BDWALLのDirichlet_P計算
			{
				if(PART[i].surface==ON)//壁表面粒子なら
				{
					double P=0;
					double W=0;
					for(int k=0;k<PART[i].N2;k++)
					{
						int j=PART[i].NEI2[k];
						if(PART[j].type==FLUID && PART[j].surface==ON)
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
							double w=kernel(r2,dis);
							P+=Dirichlet_P[j]*w;
							W+=w;
						}
					}
					if(W!=0) P/=W;
					Dirichlet_P[i]=P;
				}
			}

			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].surface==ON)
				{
					for(int k=0;k<PART[i].N2;k++)
					{
						int j=PART[i].NEI2[k];
						int m=link[j];
						if(m<pn)//粒子jが未知数なら
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
				    
							double w=kernel2(r2,dis,d);
							if(CON->get_HL_sw()==ON) w*=real_lamda;
							B[m]-=w*Dirichlet_P[i];
							PART[i].P=Dirichlet_P[i];
						}
					}
				}
			
			}

			//ﾌｧｲﾙ出力
			if(count2==0)
			{
				ofstream gg("Dirichlet_P.dat");
				if(d==2) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]<<endl;}
				else if(d==3) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;}
				gg.close();///*/
	
				
			}
			delete [] Dirichlet_P;
		}///////////
	
	
		/////////////////////////////////////CG法
		double *r=new double[pn];
		double *X=new double[pn];
		double *AP = new double [pn];
		double *P = new double [pn];

		/////////////////////////初期値//////////////////
		if(CON->get_initialP()==OFF)
		{
			for(int n=0;n<pn;n++) 
			{
				 X[n]=0;
				 r[n]=B[n];
				 P[n]=r[n];
			}
		}
		else if(CON->get_initialP()==ON)//初期値として現在の情報を与える。微小に早くなる
		{
			for(int n=0;n<pn;n++) X[n]=PART[ppn[n]].P;
			for(int n=0;n<pn;n++)
			{
				double AX=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
				r[n]=B[n]-AX;
				P[n]=r[n];
			}
		}
		//////////////////////////////////////////////

		//CG法により行列を解く
		if(CON->get_solution()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,CON->get_CGep());
		else if(CON->get_solution()==1)
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
			iccg(CON,val,ind,ptr,pn,B,number,X,r,P,CON->get_CGep(),&count);
		}
		
		if(CON->get_negativeP()==OFF)//負圧を考慮しない場合
		{
			for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				PART[i].P=X[n];
				if(PART[i].P<0) PART[i].P=0;
			}
		}
		else if(CON->get_negativeP()==ON)//負圧を考慮する場合
		{
			for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				PART[i].P=X[n];
			}
		}
		
		delete [] r;
		delete [] X;
		delete [] AP;
		delete [] P;

		//////////////////////////////*/

		cout<<"反復回数:"<<count<<"  time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;

		if(count2==0) plot_P(CON ,PART,particle_number,t,fluid_number);//最初のみ圧力出力

		///////////////////////////////////////////////////圧力勾配計算
		cout<<"圧力勾配計算開始 ver."<<CON->get_Pgrad()<<" ---------";
	
		unsigned int timeB=GetTickCount();	//計算開始時刻

		int minPsw=CON->get_minP();

		if(CON->get_Pgrad()==3)//表面のみ法線 && minP=0のとき表面粒子からも反発力計算
		{
			P_gradient3(CON,PART,fluid_number,direct,dt,reU,P_grad,minPsw);
		}
		else if(CON->get_Pgrad()==4)//重みつき最少二乗法(WLSM)
		{
			P_gradient4(CON,PART,dt,fluid_number,reU,P_grad,minPsw);
		}
		else if(CON->get_Pgrad()==5)//自身を通らない重みつき最少二乗法(WLSM)
		{
			P_gradient5(CON,PART,dt,fluid_number,reU,P_grad);
		}
	
		cout<<"ok  time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;
		//////////////////////////////////////////////////////////////////////////////*/

		//速度修正量を計算
		for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) reU[D][i]=-dt*P_grad[D][i]/density[i];
		
		///速度修正量の修正
		int mode=0;
		if(mode==1) modify_reU(CON,PART,fluid_number,particle_number,t,dt,n0,udiv,reU);


		//速度修正
		if(calctype==0)//速度のみ修正
		{
			for(int i=0;i<fluid_number;i++)
			{
				//double divreU=div_of_reU(CON,PART,reU,i,fluid_number);
				double w=1;//0.8;
				for(int D=0;D<d;D++)
				{
					PART[i].u[D]+=w*reU[D][i];//なんか0.6～0.9くらいの倍率かけたら早く収束する
					//PART[i].u[D]+=w*reU[D][i]+(1-w)*old_reU[D][i];
					old_reU[D][i]=reU[D][i];//更新
					
				}
			}
		}
		else if(calctype==1)//位置も修正。その場合は重み関数を再計算。
		{
			int sw=CON->get_temporary_r();		//ONなら陽解析後に仮の位置を計算している(通常の解析)
			for(int i=0;i<fluid_number;i++)
			{
				for(int D=0;D<d;D++)
				{
					PART[i].u[D]+=reU[D][i];
					//PART[i].r[D]+=reU[D][i]*dt;
					if(sw==ON) PART[i].r[D]+=reU[D][i]*dt*0.5;
					else 
					{
						PART[i].r[D]=old_r[D][i]+dt*(PART[i].u[D]+old_U[D][i])*0.5;
						//PART[i].r[D]+=dt*(PART[i].u[D]+old_U[D][i])*0.5;//台形則 old_Uはこの関数に突入した際の速度。反復の間ずっと記憶してる
					}
				}
			}

			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type!=OUTWALL)
				{
					double W=0;
					for(int k=0;k<PART[i].N2;k++)
					{
						int j=PART[i].NEI2[k];
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel2(r2,dis,d);
						W+=w;
					}
					PART[i].PND2=W;
				}
			}
		}////////////////////////////

		///速度の発散誤差を計算
		div_error=0;
		numf=0;//内部流体粒子数
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type!=OUTWALL && PART[i].surface==OFF)
			{
				if(divflag==1) udiv[i]=divergence(CON,PART,i,n0);
				else if(divflag==2)
				{
					udiv[i]=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
					if(2*udiv[i]!=udiv[i]+udiv[i]) udiv[i]=divergence(CON,PART,i,n0);//udiv[i]がエラーで非実数のときはMPSで計算する
				}
				else if(divflag==3) udiv[i]=Dndt(CON,PART,i,n0);
			}
			if(PART[i].type==FLUID && PART[i].surface==OFF)
			{
				div_error+=udiv[i]*udiv[i];
				numf++;
			}
		}
		if(numf!=0) 
		{
			div_error=div_error/numf;
			div_error=sqrt(div_error);	//標準偏差
		}///////*/


		if(count2>0 && div_error>=error) flag=1;//2回目以降の計算で誤差が増えてしまったら計算中止
		
		cout<<div_error<<endl;

		error=div_error;//誤差更新（どんどん小さくなる）
		count2++;

		if(count2>CON->get_iteration_count()) flag=1;//脱出
		if(mode==1) flag=1;
		
	}
	lamda=real_lamda;//lamdaの値を戻す

	//速度修正と位置修正
	if(calctype==0) modify_u_and_x_after_Pcalc(CON,PART, fluid_number,reU,dt,Un);

	//最終的な速度発散分布出力
	ofstream ff("div.dat");
	if(d==2) for(int i=0;i<fluid_number;i++) ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<udiv[i]<<endl;
	if(d==3) for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) ff<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<udiv[i]<<endl;
	ff.close();

	for(int D=0;D<DIMENTION;D++) delete [] reU[D];
	for(int D=0;D<DIMENTION;D++) delete [] old_U[D];
	for(int D=0;D<DIMENTION;D++) delete [] old_r[D];
	for(int D=0;D<DIMENTION;D++) delete [] old_reU[D];
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	for(int D=0;D<DIMENTION;D++) delete [] P_grad[D];

	delete [] B;
	delete [] ppn;
	delete [] link;

	delete [] val;
	delete [] ind;
	delete [] ptr;

	delete [] udiv;

	delete [] density;

	cout<<"ok count="<<count2-1<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//速度修正量の修正関数
void modify_reU(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double n0,double *udiv, double **reU)
{
	cout<<"reUの修正開始--";
	///ここでもらうn0は拡散のn0なのでは？要チェック
	double le=CON->get_distancebp();
	double R=CON->get_re()*le;
	int    d=CON->get_dimention();
	double dir_val=1;				//ディリクレ値　現在は1を使用 たとえば表面粒子だけは個人勝手に決めそれをﾃﾞｨﾘｸﾚ値にしたたらどうなる？無理か？

	int pn=0;	//未知数
	for(int i=0;i<fluid_number;i++) if(PART[i].surface==OFF) pn++;//現在のところ、表面粒子は固定境界としている


	int *ppn = new int[pn];					//行列における第n番目の未知数は粒子番号ppn[n]の粒子に相当
	int *link = new int [particle_number];		//粒子番号iはlink[i]番目の未知数
	double *B   = new double[pn];			//解行列

	int count=0;
	///ppn配列とlink配列,および解行列Ｂ作成
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==OFF)
		{
			ppn[count]=i;
			link[i]=count;
			B[count]=-udiv[i];
			count++;
		}
		else link[i]=pn+1;//行列に含まれない粒子にはﾀﾞﾐｰとして(pn+1)を格納
	}	
	//////*/

	int number=0;			//係数行列の非ゼロ要素数
	for(int n=0;n<pn;n++)
	{   
		number++;			//自分自身を数にいれる
		int i=ppn[n];		//n番目の未知数は粒子i
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			int m=link[j];
			if(m<pn) number++;
		}
	}
	////numberが求まった

	double *val = new double [number];
	int *ind = new int [number];//非ゼロ要素の列番号格納配列
	int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
	
	
	/////////////////////val,ind ,ptrに値を格納
	
	int index=0;
	for(int n=0;n<pn;n++)
	{
	    ptr[n]=index;
	    int i=ppn[n];
		int kk=index;//値を保存
	    val[index]=-PART[i].PND;
	    ind[index]=n;
	    index++;
		double W=0;
		double d_val=0;//対角線上の値
		//n0=PART[i].PND;
	    for(int k=0;k<PART[i].N;k++)
	    {
	        int j=PART[i].NEI[k];
			
			int m=link[j];
			if(m<pn)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R,dis);
				val[index]=d/n0*w*(reU[A_X][j]*X+reU[A_Y][j]*Y+reU[A_Z][j]*Z)/(dis*dis);
				ind[index]=m;
				index++;
				W+=w;
	
				d_val+=-d/n0*w*(reU[A_X][i]*X+reU[A_Y][i]*Y+reU[A_Z][i]*Z)/(dis*dis);
			}
			else if(PART[j].type==FLUID && PART[j].surface==ON)
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R,dis);
					
				B[n]-=d/n0*w*(reU[A_X][j]*X+reU[A_Y][j]*Y+reU[A_Z][j]*Z)/(dis*dis)*dir_val;
				W+=w; 
				d_val+=-d/n0*w*(reU[A_X][i]*X+reU[A_Y][i]*Y+reU[A_Z][i]*Z)/(dis*dis);
			}
			else	//壁粒子
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R,dis);
					
				W+=w; 
				d_val+=-d/n0*w*(reU[A_X][i]*X+reU[A_Y][i]*Y+reU[A_Z][i]*Z)/(dis*dis);
			}//*/
			
		}
		val[kk]=d_val;
		//val[kk]=-W;
	}
	ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない
	cout<<"行列作成 ";
	////////////////////*/
	

	/////////////////////////////////////GCR法
	double alp;
	double E=1;//誤差
	unsigned int timeC=GetTickCount();
	
	double *r=new double[pn];
	double *XX=new double[pn];
	for(int n=0;n<pn;n++) XX[n]=0;
	
	double *AP = new double [pn];
	double *Ar = new double [pn];
	double *P = new double [pn];
	double *nextP = new double [pn];

	int max=500;
	double **oldP=new double *[max];
	for(int m=0;m<max;m++) oldP[m]=new double [pn];
	double **oldAP=new double *[max];
	for(int m=0;m<max;m++) oldAP[m]=new double [pn];
	double *APjAPj=new double [max];

	/////////////////////////初期値//////////////////
	
	for(int n=0;n<pn;n++) XX[n]=dir_val;
	/*for(int n=0;n<pn;n++)
	{
		double AX=0;
		for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*XX[ind[j]];
		r[n]=B[n]-AX;
		P[n]=r[n];
	}/////////*/
		
	//////////////////////////////////////////////

	double EP=1e-2;//CON->get_CGep();//収束判定(convergence test)

	
	cout<<"GCR法スタート------";//このﾌﾟﾛｸﾞﾗﾑは『新数値計算』福井 義成のP46を参考にしている
	for(int n=0;n<pn;n++)
	{
		double AX=0;
		for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*XX[ind[j]];
		r[n]=B[n]-AX;
		P[n]=r[n];
	}/////////*/
	
	count=0;
	
	while(E>EP && count<max)
	{
		for(int n=0;n<pn;n++) oldP[count][n]=P[n];//探索方向ﾍﾞｸﾄﾙの保存

		//////////////alpを求める alp=(ri,APi)/(APi,APi)
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
		}
		for(int n=0;n<pn;n++) oldAP[count][n]=AP[n];//値を保存


		double rAP=0;
		for(int n=0;n<pn;n++) rAP+=r[n]*AP[n];
		double APAP=0;
		for(int n=0;n<pn;n++) APAP+=AP[n]*AP[n];
		alp=rAP/APAP;
			//cout<<"alp="<<alp<<" "<<rAP<<" "<<APAP<<endl;
		//////////////////////

		//(AP,AP)を保存
		APjAPj[count]=APAP;
		////////////
		
		//////////////// 解更新　X(k+1)=X(k)+alp*P
		for(int n=0;n<pn;n++) XX[n]+=alp*P[n];
		//////////////////////////////
			
		//////////////// r=r-alp*AP
		for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
			
		//////////////////誤差
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
			
		////////////////////////とりあえず次の探索方向として残差ｒを格納
		for(int n=0;n<pn;n++) nextP[n]=r[n];
		///////////////////////////////////

		////////////////////////ΣβP 過去の値を使って計算
		double betaP=0;//ΣβP
		for(int n=0;n<pn;n++)
		{      
			Ar[n]=0;		
			for(int j=ptr[n];j<ptr[n+1];j++) Ar[n]+=val[j]*r[ind[j]];
		}

		for(int j=0;j<=count;j++)
		{
			double ArAPj=0;
			for(int n=0;n<pn;n++) ArAPj+=Ar[n]*oldAP[j][n];
			double betaj=-ArAPj/APjAPj[j];
			
			for(int n=0;n<pn;n++) nextP[n]+=betaj*oldP[j][n];
		}//////////*/
			
		///////////////////// Pの更新
		for(int n=0;n<pn;n++) P[n]=nextP[n];
		
		count++;
	}
	
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		for(int D=0;D<d;D++) reU[D][i]*=XX[n];

	}

	ofstream fx("kreU.dat");
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fx<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<XX[n]<<endl;

	}
	fx.close();
		
	///配列削除
	delete [] r;
	delete [] XX;
	delete [] AP;
	delete [] Ar;
	delete [] P;
	delete [] nextP;

	//////////////////////////////*/

	cout<<"反復回数:"<<count<<"  time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;


	delete [] ppn;
	delete [] link;
	delete [] B;

	delete [] val;
	delete [] ind;
	delete [] ptr;

	
	for(int m=0;m<max;m++) delete [] oldP[m];
	delete [] oldP;
	for(int m=0;m<max;m++) delete [] oldAP[m];
	delete [] oldAP;
	delete [] APjAPj;
	//////////

	
}

void calc_P_main_with_gridless(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double lamda,double N0,double dt,double *PND2,double n0,int fluid_number,int out,int pn,int B_of_P,int t)
{
	///圧力計算において、使用する重み関数は勾配などに使用するそれとは異なるものを用いる。
	//理由は、勾配に使用する重み関数は単なる重み平均なので何をしようしてもよいが、
	//PPE式の離散化に使用する重み関数は、粒子数密度∝密度となるように設定する必要があるから。
	//もし両者で同一の重み関数を使用したければ、kernel2()の中身をkernelと同一に書き換えること。

	//pn:圧力を解くための連立方程式未知数
	cout<<"gridlessによる圧力計算"<<endl;

	double le=CON->get_distancebp();
	double r2=CON->get_re2()*le;
	unsigned int timeA=GetTickCount();					//計算開始時刻
	double dimention=2;
	if(CON->get_dimention()==3) dimention=3;
	double *density=new double[particle_number];
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	for(int i=fluid_number;i<particle_number;i++) density[i]=CON->get_density();
	//double lamda=calclamda(CON); //ﾗﾌﾟﾗｼｱﾝ用λ

	//圧力解析のためにN0とPART[i].PND2をkernel2()用におきかえる
	//set_N0_and_PND2(CON,PART,particle_number,fluid_number,&N0,out);
	//cout<<"N0="<<N0<<endl;
	///////////////////*/

	//この関数に入った段階でpnが決まっている仕様になっているが、gridless法で行う場合は、周辺粒子数が少ない内部粒子の演算を省きたい。
	//そうするとpnが変化することになる。そこで、再度pnを求める。

	int min_neiber=4;				//周辺粒子数がこれより少ない場合は圧力を計算しない
	if(dimention==3) min_neiber=6;
	pn=0;
	for(int i=0;i<particle_number;i++) if(PART[i].type!=OUTWALL && PART[i].surface==OFF && PART[i].N2>=min_neiber) pn++;

	cout<<"圧力未知数:"<<pn<<" ";

	int *ppn = new int[pn];					//行列における第n番目の未知数は粒子番号ppn[n]の粒子に相当
	int *link = new int [particle_number];	//粒子番号iはlink[i]番目の未知数
	double *B   = new double[pn];			//解行列
	
	int count=0;

	///ppn配列とlink配列作成
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type!=OUTWALL && PART[i].surface==OFF && PART[i].N2>=min_neiber)
		{
			ppn[count]=i;
			if(B_of_P==0)
			{
				//B[count]=-density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());//教科書
				B[count]=-density[i]*(PART[i].PND2-N0)/(dt*dt*N0);//教科書
			}
			else if(B_of_P==1)//速度発散
			{    
				//double div=divergence(CON,PART,i,n0);
				double div;
				if(CON->get_divU_method()==1) div=divergence(CON,PART,i,n0);
				else if(CON->get_divU_method()==2) div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());//最小二乗法（自身を通る曲面）
				else if(CON->get_divU_method()==3) div=divergence3(CON,PART,i,CON->get_interpolate_surface_u());//最小二乗法（自身を通るとは限らない曲面）
				else if(CON->get_divU_method()==4) div=divergence4(CON,PART,i);//入部ら
				if(2*div!=div+div)
				{
					div=divergence(CON,PART,i,n0);//divがエラーで非実数のときはMPSで計算する
					//cout<<"div_err i="<<i<<" "<<div<<endl;
				}
				
				//B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i];
				//B[count]=div*lamda*PART[i].PND2/(dt*2*CON->get_dimention())*density[i];
				B[count]=div/(dt)*density[i];
			}
			else if(B_of_P==2)//(教科書+速度発散)/2
			{
				double a=CON->get_w_div();double b=1;///各手法の重み
				//double div=divergence(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//divがエラーで非実数のときはMPSで計算する
				
				B[count]=div*lamda*N0/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;
				B[count]/=a+b;
			}	
			else if(B_of_P==3)//重み関数の直接微分による速度発散
			{    
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
			}
			else if(B_of_P==4)//重み関数の直接微分による速度発散+PND
			{    
				double a=CON->get_w_div();double b=1;///各手法の重み
				double div=Dndt(CON,PART,i,n0);
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i]*a;
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention())*b;//教科書
				B[count]/=a+b;
			}
			else if(B_of_P==5)//(ni-nk)+(nk-n0)=div+pnd(nk-n0)
			{    
				//double div=Dndt(CON,PART,i,n0);
				double div=divergence2(CON,PART,i,CON->get_interpolate_surface_u());
				if(2*div!=div+div) div=divergence(CON,PART,i,n0);//divがエラーで非実数のときはMPSで計算する
				B[count]=-div*lamda/(dt*2*CON->get_dimention())*density[i];
				B[count]-=density[i]*lamda*(PART[i].PND2-N0)/(dt*dt*2*CON->get_dimention());
			}
			link[i]=count;
			count++;
		}
		else
		{
			//if(PART[i].N2<=6) cout<<i<<" "<<PART[i].N2<<endl;
			link[i]=pn+1;//行列に含まれない粒子にはﾀﾞﾐｰとして(pn+1)を格納
		}
	}
	//////*/

	if(CON->get_dir_for_P()!=OFF && B_of_P!=0)//表面粒子の圧力として、ディリクレ値
	{
		int flag=CON->get_dir_for_P();
		if(flag==2 || flag==3)
		{
			if(CON->get_EM_method()==OFF) flag=1;//BEMによる計算を行わないのにflagが2や3なら間違いなので1に戻す
		}
		double *Dirichlet_P=new double [particle_number];
		for(int i=0;i<particle_number;i++) Dirichlet_P[i]=0;

		if(flag==1) set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
		else if(flag==2) set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		else if(flag==3)
		{
			set_Dirichlet_P(CON,PART,fluid_number,1,Dirichlet_P);
			set_Dirichlet_P(CON,PART,fluid_number,2,Dirichlet_P);
		}
		
		for(int i=fluid_number;i<particle_number;i++)//壁表面粒子のDirichlet_P計算
		{
			if(PART[i].surface==ON)//壁表面粒子なら
			{
				double P=0;
				double W=0;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k];
					if(PART[j].type==FLUID && PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						double w=kernel(r2,dis);
						P+=Dirichlet_P[j]*w;
						W+=w;
					}
				}
				if(W!=0) P/=W;
				Dirichlet_P[i]=P;
			}
		}

		//表面粒子のPART[i].Pにディリクレ境界条件を代入
		ofstream gg2("Dirichlet_error.dat");
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].surface==ON) PART[i].P=Dirichlet_P[i];
			else if(PART[i].type==FLUID && Dirichlet_P[i]!=0)
			{
				cout<<"内部粒子なのにﾃﾞｨﾘｸﾚ値が格納されています。i="<<i<<" 値は"<<Dirichlet_P[i]<<endl;
				gg2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		gg2.close();

		//ﾌｧｲﾙ出力
		ofstream gg("Dirichlet_P.dat");
		if(dimention==2) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<Dirichlet_P[i]<<endl;}
		else if(dimention==3) {for(int i=0;i<particle_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) gg<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<Dirichlet_P[i]<<endl;}
		gg.close();////

		//ディリクレ値をベクトル表示
		output_dirichlet_vector_files(CON,PART,fluid_number,flag,Dirichlet_P);

		delete [] Dirichlet_P;
	}//////////*/

	///解行列出力
	ofstream h("Bmatrix_for_P.dat");
	if(CON->get_dimention()==2) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<B[n]<<endl;
		}
	}
	if(CON->get_dimention()==3) 
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			//cout<<B[n]<<endl;
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le) h<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<B[n]<<endl;
		}
	}
	h.close();
	/////////////*/

	
	int number=0;			//係数行列の非ゼロ要素数
	for(int n=0;n<pn;n++)
	{   
	    number++;///自分自身を数にいれる
	    int i=ppn[n];//n番目の未知数は粒子i
	  
	    for(int k=0;k<PART[i].N2;k++)
	    {
	        int j=PART[i].NEI2[k];
			int m=link[j];
			if(m<pn) number++;
	    }
	}
	////numberが求まった
	
    double *val = new double [number];
	int *ind = new int [number];//非ゼロ要素の列番号格納配列
	int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
	
	/////////////////////val,ind ,ptrに値を格納
	unsigned int timeC=GetTickCount();
	ofstream ft("gridless_to_MPS.dat");				//離散化手法がgridlessからMPSに変更になった粒子の座標を出力
	if(CON->get_dimention()==2)//3次元
	{
		int N=5;
		double *matrix=new double [N*N];	//N×Nの係数行列
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//値を保存
			ind[index]=n;
			index++;
			double valII=0;						//対角成分の値
			for(int k=0;k<N*N;k++) matrix[k]=0;//初期化
			double B0=B[n];			//B[n]の値を保存
			int index0=index;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				
				int m=link[j];
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				double w=1;
				//if(dis>le) w=le*le/(dis*dis);
				if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
				//if(dis>1) w=1/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
				//if(PART[j].surface==ON) w=10;	//2回、足すのと同じ事
				//if(PART[j].surface==ON) w=1;//表面はきちんと通って欲しいから、重みを大きくとる

				matrix[0]+=X*X*w;			//ΣXjwj
				matrix[1]+=X*Y*w;		//ΣXjYjwj
				matrix[2]+=X*X*X*w;			//ΣXj^3wj
				matrix[3]+=X*X*Y*w;			//ΣXj^2Yjwj
				matrix[4]+=X*Y*Y*w;			//ΣXjYj^2wj
	
				matrix[6]+=Y*Y*w;			//ΣYj^2wj
				matrix[9]+=Y*Y*Y*w;			//ΣYj^3wj
	
				matrix[12]+=X*X*X*X*w;			//ΣXj^4wj
				matrix[13]+=X*X*X*Y*w;			//ΣXj^3Yjwj
				matrix[14]+=X*X*Y*Y*w;			//ΣXj^2Yj^2wj
		
				matrix[19]+=X*Y*Y*Y*w;			//ΣXjYj^3wj
	
				matrix[24]+=Y*Y*Y*Y*w;			//ΣYj^4wj
	
				//B[0]+=dP*X*w;//ΣdPjXjwj
				//B[1]+=dP*Y*w;//ΣdPjYjwj
				//B[2]+=dP*X*X*w;//ΣdPjXj^2wj
				//B[3]+=dP*X*Y*w;//ΣdPjXjYjwj
				//B[4]+=dP*Y*Y*w;//ΣdPjYj^2wj
				
				/*if(PART[j].surface==ON)
				{
					X*=2; Y*=2; Z*=2; dis*=2;
					
					//if(dis>le) w=le*le/(dis*dis);
					if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる

					matrix[0]+=X*X*w;			//ΣXjwj
					matrix[1]+=X*Y*w;		//ΣXjYjwj
					matrix[2]+=X*X*X*w;			//ΣXj^3wj
					matrix[3]+=X*X*Y*w;			//ΣXj^2Yjwj
					matrix[4]+=X*Y*Y*w;			//ΣXjYj^2wj
	
					matrix[6]+=Y*Y*w;			//ΣYj^2wj
					matrix[9]+=Y*Y*Y*w;			//ΣYj^3wj
	
					matrix[12]+=X*X*X*X*w;			//ΣXj^4wj
					matrix[13]+=X*X*X*Y*w;			//ΣXj^3Yjwj
					matrix[14]+=X*X*Y*Y*w;			//ΣXj^2Yj^2wj
		
					matrix[19]+=X*Y*Y*Y*w;			//ΣXjYj^3wj
	
					matrix[24]+=Y*Y*Y*Y*w;			//ΣYj^4wj
				}*/
			}

			matrix[5]=matrix[1];		//ΣXjYjwj
			matrix[7]=matrix[3];		//ΣXj^2Yjwj
			matrix[8]=matrix[4];		//ΣXjYj^2wj
			matrix[10]=matrix[2];		//ΣXj^3Yjwj
			matrix[11]=matrix[3];		//ΣXj^2Yjwj
			matrix[15]=matrix[3];		//ΣXj^2Yjwj
			matrix[16]=matrix[4];		//ΣXjYj^2wj
			matrix[17]=matrix[13];		//ΣXj^3Yjwj
			matrix[18]=matrix[14];		//ΣXj^2Yj^2wj
			matrix[20]=matrix[4];		//ΣXjYj^2wj
			matrix[21]=matrix[9];		//ΣYj^3wj
			matrix[22]=matrix[14];		//ΣXj^2Yj^2wj
			matrix[23]=matrix[19];		//ΣXjYj^3wj

			//matrixの逆行列を求める 求まった逆行列はmatrixの中を上書きして格納される

			/*if(n==110)
			{
				cout<<endl;
				for(int k=0;k<N;k++)
				{
					for(int k2=0;k2<N;k2++)
					{
						cout<<matrix[k*N+k2]<<" ";
					}
					cout<<endl;
				}
			}*/
			
			calc_inverse_matrix(CON,PART, N, matrix);

			int flag=ON;

			for(int k=0;k<PART[i].N2;k++)
			{
				//if(flag==ON)
				{
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
						
					double w=1;
					if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
					//if(dis>1) w=1/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
					//if(PART[j].surface==ON) w=1;//表面はきちんと通って欲しいから、重みを大きくとる
				
					//変数cの行×jに関する解行列
					double valX=matrix[10]*X+matrix[11]*Y+matrix[12]*X*X+matrix[13]*X*Y+matrix[14]*Y*Y;
				
					//変数eの行×jに関する解行列
					double valY=matrix[20]*X+matrix[21]*Y+matrix[22]*X*X+matrix[23]*X*Y+matrix[24]*Y*Y;

					if(valX+valY<0) flag=OFF;

				//	if(flag==ON)
					{
						if(m<pn)
						{
							val[index]=2*w*(valX+valY);		//2階微分値はd,e,fを2倍したものである点に注意
							ind[index]=m;
							index++;
							valII+=2*w*(valX+valY);
						}
						else if(PART[j].surface==ON)
						{
							valII+=2*w*(valX+valY);//ここで値を計算に組み込むとき、粒子jはﾃﾞｨﾘｸﾚ型、組み込まないと勾配ゼロのﾉｲﾏﾝ型
							B[n]-=2*w*(valX+valY)*PART[j].P;		//すでに表面粒子にはディリクレ型の圧力が代入されている。
						}
					}
				}
			}
			if(valII<0) flag=OFF;
			if(flag==OFF)
			{
				valII=0;
				B[n]=B0;
				index=index0;
				ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				for(int k=0;k<PART[i].N2;k++)
				{
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					if(m<pn)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					    
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						val[index]=val2;
						ind[index]=m;
						index++;
						valII+=val2;
					}
					else if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						//double w=kernel(r2,dis);
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						valII+=val2; //ここでwを計算に組み込むとき、粒子jはﾃﾞｨﾘｸﾚ型、組み込まないと勾配ゼロのﾉｲﾏﾝ型
						B[n]-=val2*PART[j].P;		//すでに表面粒子にはディリクレ型の圧力が代入されている。
					}
				}

			}
			val[kk]=-valII;
		}
		ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない

		delete [] matrix;
	}
	else if(CON->get_dimention()==3)//3次元
	{
		int N=9;
		double *matrix=new double [N*N];	//N×Nの係数行列
		int index=0;
		for(int n=0;n<pn;n++)
		{
			//cout<<n<<"/"<<pn<<endl;
			ptr[n]=index;
			int i=ppn[n];
			int kk=index;//値を保存
			ind[index]=n;
			index++;
			double valII=0;						//対角成分の値
			for(int k=0;k<N*N;k++) matrix[k]=0;//初期化
			double B0=B[n];			//B[n]の値を保存
			int index0=index;		//indexの値を保存
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
				//double dP=(PART[j].P-minP[i]);//Δfjに相当

				double w=1;
				//if(dis>le) w=le*le/(dis*dis);
				//if(dis>le) w=le/(dis);
				if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
				//if(PART[j].surface==ON) w=1;//表面はきちんと通って欲しいから、重みを大きくとる

				matrix[0]+=X*X*w;		//ΣXj^2wj
				matrix[1]+=X*Y*w;		//ΣXjYjwj
				matrix[2]+=X*Z*w;		//ΣXjZjwj
				matrix[3]+=X*X*X*w;		//ΣXj^3wj
				matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
				matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
				matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
				matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
				matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj
	
				matrix[10]+=Y*Y*w;		
				matrix[11]+=Y*Z*w;		
				matrix[12]+=X*X*Y*w;		
				matrix[13]+=Y*Y*Y*w;		
				matrix[14]+=Y*Z*Z*w;
				matrix[15]+=X*Y*Y*w;
				matrix[16]+=Y*Y*Z*w;
							
				matrix[20]+=Z*Z*w;			
				matrix[23]+=Z*Z*Z*w;		
					
				matrix[30]+=X*X*X*X*w;
				matrix[31]+=X*X*Y*Y*w;
				matrix[32]+=X*X*Z*Z*w;	
				matrix[33]+=X*X*X*Y*w;	
				matrix[34]+=X*X*Y*Z*w;	
				matrix[35]+=X*X*X*Z*w;	
					
				matrix[40]+=Y*Y*Y*Y*w;
				matrix[41]+=Y*Y*Z*Z*w;
				matrix[42]+=X*Y*Y*Y*w;
				matrix[43]+=Y*Y*Y*Z*w;
				matrix[44]+=X*Y*Y*Z*w;

				matrix[50]+=Z*Z*Z*Z*w;	//6行目
				matrix[51]+=X*Y*Z*Z*w;
				matrix[52]+=Y*Z*Z*Z*w;
				matrix[53]+=X*Z*Z*Z*w;
	
				//7～9行目はすべて既存の要素から転用が可能


				/*B[0]+=dP*X*w;		//a
				B[1]+=dP*Y*w;		//b
				B[2]+=dP*Z*w;		//c
				B[3]+=dP*X*X*w;		//d	(Xの2階微分)
				B[4]+=dP*Y*Y*w;		//e(Yの2階微分)
				B[5]+=dP*Z*Z*w;		//f(Zの2階微分)
				B[6]+=dP*X*Y*w;		//g
				B[7]+=dP*Y*Z*w;		//h
				B[8]+=dP*X*Z*w;		//i*/
			}

			matrix[9]=matrix[1];		//ΣXjYjwj
			matrix[17]=matrix[7];

			matrix[18]=matrix[2];
			matrix[19]=matrix[11];
			matrix[21]=matrix[8];
			matrix[22]=matrix[16];
			matrix[24]=matrix[7];
			matrix[25]=matrix[14];
			matrix[26]=matrix[5];

			matrix[27]=matrix[3];
			matrix[28]=matrix[12];
			matrix[29]=matrix[21];

			matrix[36]=matrix[4];
			matrix[37]=matrix[13];
			matrix[38]=matrix[22];
			matrix[39]=matrix[31];

			matrix[45]=matrix[5];
			matrix[46]=matrix[14];
			matrix[47]=matrix[23];
			matrix[48]=matrix[32];
			matrix[49]=matrix[41];

			matrix[54]=matrix[6];
			matrix[55]=matrix[15];
			matrix[56]=matrix[24];
			matrix[57]=matrix[33];
			matrix[58]=matrix[42];
			matrix[59]=matrix[51];
			matrix[60]=matrix[31];
			matrix[61]=matrix[44];
			matrix[62]=matrix[34];

			matrix[63]=matrix[7];
			matrix[64]=matrix[16];
			matrix[65]=matrix[25];
			matrix[66]=matrix[34];
			matrix[67]=matrix[43];
			matrix[68]=matrix[52];
			matrix[69]=matrix[61];
			matrix[70]=matrix[41];
			matrix[71]=matrix[51];

			matrix[72]=matrix[8];
			matrix[73]=matrix[17];
			matrix[74]=matrix[26];
			matrix[75]=matrix[35];
			matrix[76]=matrix[44];
			matrix[77]=matrix[53];
			matrix[78]=matrix[62];
			matrix[79]=matrix[71];
			matrix[80]=matrix[32];

			//matrixの逆行列を求める 求まった逆行列はmatrixの中を上書きして格納される
			
			calc_inverse_matrix(CON,PART, N, matrix);

			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k]; 
				int m=link[j];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
						
				double w=1;
				//if(dis>le) w=le*le/(dis*dis);
				//if(dis>le) w=le/(dis);
				if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
				//if(PART[j].surface==ON) w=1;//表面はきちんと通って欲しいから、重みを大きくとる

				
				//変数dの行×jに関する解行列
				double valX=matrix[27]*X+matrix[28]*Y+matrix[29]*Z+matrix[30]*X*X+matrix[31]*Y*Y+matrix[32]*Z*Z+matrix[33]*X*Y+matrix[34]*Y*Z+matrix[35]*X*Z;
				
				//変数eの行×jに関する解行列
				double valY=matrix[36]*X+matrix[37]*Y+matrix[38]*Z+matrix[39]*X*X+matrix[40]*Y*Y+matrix[41]*Z*Z+matrix[42]*X*Y+matrix[43]*Y*Z+matrix[44]*X*Z;
				
				double valZ=matrix[45]*X+matrix[46]*Y+matrix[47]*Z+matrix[48]*X*X+matrix[49]*Y*Y+matrix[50]*Z*Z+matrix[51]*X*Y+matrix[52]*Y*Z+matrix[53]*X*Z;
				
				//if(valX+valY+valZ<0) w=0;					//普通は、正になるべきものだと思う。負なら考慮しない
				if(m<pn)
				{
					val[index]=2*w*(valX+valY+valZ);		//2階微分値はd,e,fを2倍したものである点に注意
					ind[index]=m;
					index++;
					valII+=2*w*(valX+valY+valZ);
				}
				else if(PART[j].surface==ON)
				{
					valII+=2*w*(valX+valY+valZ);//ここで値を計算に組み込むとき、粒子jはﾃﾞｨﾘｸﾚ型、組み込まないと勾配ゼロのﾉｲﾏﾝ型
					B[n]-=2*w*(valX+valY+valZ)*PART[j].P;		//すでに表面粒子にはディリクレ型の圧力が代入されている。
				}
			}
			int flag=ON;
			if(valII<0) flag=OFF;
			else if(valII+valII!=2*valII) flag=OFF;
			if(flag==OFF)
			{
				valII=0;
				B[n]=B0;
				index=index0;
				ft<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<i<<endl;
				for(int k=0;k<PART[i].N2;k++)
				{
				
					int j=PART[i].NEI2[k]; 
					int m=link[j];
					if(m<pn)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					    
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						val[index]=val2;
						ind[index]=m;
						index++;
						valII+=val2;
					}
					else if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						//double w=kernel(r2,dis);
						double w=kernel2(r2,dis,dimention);
						double val2=2.0*CON->get_dimention()/(lamda*PART[i].PND2)*w;
						valII+=val2; //ここでwを計算に組み込むとき、粒子jはﾃﾞｨﾘｸﾚ型、組み込まないと勾配ゼロのﾉｲﾏﾝ型
						B[n]-=val2*PART[j].P;		//すでに表面粒子にはディリクレ型の圧力が代入されている。
					}
				}

			}
			val[kk]=-valII;
		}
		ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない

		delete [] matrix;
	}
	ft.close();

	cout<<"行列作成-time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
        
	/////////////////////////////////////CG法
    double *r=new double[pn];
	double *X=new double[pn];
	
	double *AP = new double [pn];
	double *P = new double [pn];
	/////////////////////////初期値//////////////////
    if(CON->get_initialP()==OFF)
	{
		for(int n=0;n<pn;n++) 
		{
			 X[n]=0;
			 r[n]=B[n];
			 P[n]=r[n];
		}
	}
	else if(CON->get_initialP()==ON)//初期値として現在の情報を与える。微小に早くなる
	{
		for(int n=0;n<pn;n++) X[n]=PART[ppn[n]].P;
		for(int n=0;n<pn;n++)
		{
			double AX=0;
			for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
			r[n]=B[n]-AX;
			P[n]=r[n];
		}
	}
	//////////////////////////////////////////////

	//BiCGStab2_method_with_D_scale_for_sparse(CON,r, pn,X,&count,1e-12,val,ind,ptr);//なんか、非対称な行列は、収束判定をかなり小さくしないといけないみたい。
	BiCGStab2_method_with_D_scale_for_sparse(CON,r, pn,X,&count,1e-12,val,ind,ptr);//なんか、非対称な行列は、収束判定をかなり小さくしないといけないみたい。
	//MRTR(CON,r, pn,X,&count,1e-10,val,ind,ptr);
	
	if(CON->get_negativeP()==OFF)//負圧を考慮しない場合
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			PART[i].P=X[n];
			if(PART[i].P<0) PART[i].P=0;
		}
	}
	else if(CON->get_negativeP()==ON)//負圧を考慮する場合
	{
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			if(2*X[n]==X[n]+X[n])PART[i].P=X[n];
			//else cout<<"!!! i="<<i<<endl;
		}
	}
	
		
    delete [] r;
	delete [] X;
	delete [] AP;
	delete [] P;

    //////////////////////////////*/

	//CG_GPU(CON,PART,val,ind,ptr,pn,number,ppn,B,&count);

    delete [] val;
	delete [] ind;
	delete [] ptr;
    delete [] B;
	delete [] ppn;
	delete [] link;

	delete [] density;

	cout<<"反復回数:"<<count<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//最低粒子間距離より粒子数密度を推算する。
void calc_PND_by_minmum_L_main(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number)
{
	int dimention=CON->get_dimention();
	double Re=CON->get_re2();
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].N>0)
		{
			double minL=100;//最低粒子間距離
			double minL2=110;
			double minL3=110;
			double minL4=110;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(dis<minL) minL=dis;
				else if(dis<minL2) minL2=dis;
				else if(dis<minL3) minL3=dis;
				else if(dis<minL4) minL4=dis;
			}
			//minLが求まった
			//cout<<minL<<" "<<minL2<<endl;
			double ave=(minL+minL2+minL3+minL4)/4;
			//cout<<i<<" "<<PART[i].PND2<<" "<<ave<<" "<<CON->get_distancebp()<<endl;
			double before=PART[i].PND2;
			PART[i].PND2=calc_PND_by_minmum_L(CON, dimention, Re, ave);
			//if(i==586) cout<<i<<" "<<PART[i].PND2<<" "<<before<<endl;
		}
	}
}

//最低粒子間距離より粒子数密度を推算する。
double calc_PND_by_minmum_L(mpsconfig *CON,int dimention,double Re,double L)
{
	//int size = (int)(r+1);//計算領域
	int size = (int)(Re+2);//計算領域
	int calc_type=CON->get_model_set_way();
	double le=CON->get_distancebp();
	double dis;//距離
	double pnd=0;
	if(dimention==2)
	{
		if(calc_type==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis*L<=Re*le )
					{
						pnd+=kernel(Re*le,dis*L);
					}			
				}
			}
		}
		if(calc_type==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
					dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis*L<=Re*le )
					{
						pnd+=kernel(Re*le,dis*L);
					}			
				}
			}
		}
	}
	else if(dimention==3)
	{
		if(calc_type==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis=sqrt((double)(i*i+j*j+k*k));
						if(dis!=0 && dis*L<=Re*le )
						{
							pnd+=kernel(Re*le,dis*L);
						}
					}			
				}
			}
		}
		if(calc_type==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-2*size;i<=2*size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//kが奇数ならiiとjjをずらす
						dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis*L<=Re*le )
						{
							pnd+=kernel(Re*le,dis*L);
						}
					}			
				}
			}
		}
	}
	return pnd;  
}

//逆行列を求める関数
void calc_inverse_matrix(mpsconfig *CON,vector<mpsparticle> &PART,int N, double *matrix)
{
	//N:未知数
	double buf;
	int i,j,k;


	double **inv_a=new  double*[N];					//逆行列格納
	for(int n=0;n<N;n++) inv_a[n]=new double [N]; 
	double **a=new  double*[N];						//matrixの値をコピー
	for(int n=0;n<N;n++) a[n]=new double [N]; 

	int count=0;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			a[i][j]=matrix[count];
			count++;
		}
	}
	//単位行列を作る
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			if(i==j) inv_a[i][j]=1;
			else		inv_a[i][j]=0;
		}
	}
	//掃き出し法
	for(i=0;i<N;i++)
	{
		buf=1/a[i][i];
		for(j=0;j<N;j++)
		{
			a[i][j]*=buf;
			inv_a[i][j]*=buf;
		}
		for(j=0;j<N;j++)
		{
			if(i!=j)
			{
				buf=a[j][i];
				for(k=0;k<N;k++)
				{
					a[j][k]-=a[i][k]*buf;
					inv_a[j][k]-=inv_a[i][k]*buf;
				}
			}
		}
	}

	count=0;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			matrix[count]=inv_a[i][j];
			count++;
		}
	}


	for(int n=0;n<N;n++) delete [] inv_a[n];
	delete [] inv_a;
	for(int n=0;n<N;n++) delete [] a[n];
	delete [] a;
}

//MRTR法//
void MRTR(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double *val,int *ind,int *ptr)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//誤差
	double BB=0;
	double rr=0;
		
	double zeta;
	double zeta_old;
	double eta;
	double v;

	double Ar_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	double Ar_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	double y_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	double Ar_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

    //double *r= new double[pn];
	double *Ar = new double [pn];  //   _ _
	
	double *P = new double [pn];
	double *y = new double [pn];

	zeta=0;
	zeta_old=0;
	eta=0.0;
	//v=complex<double>(1.0,0.0);//1回目の反復において、0でないことに注意
	v=1.0;

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0.0;
		//r[n]=B[n];
		P[n]=r[n];
		y[n]=0.0;
	}

	for(int n=0;n<pn;n++) BB+=r[n]*r[n];
	//BB=sqrt(BB);
	//cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;
	
	cout<<"MRTR法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		//while(E>EP)// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>EP && count<=pn)// EP=CON->get_CGep();//収束判定(convergence test)
		{
			//if(count==pn) cout<<"count=pn"<<endl;
			count++;

			for(int n=0;n<pn;n++)//Ar,Ar_c計算
			{    
				Ar[n]=0.0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
				}
			}

			Ar_r=0.0;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
			Ar_Ar=0.0;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			y_Ar=0.0;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			Ar_y=0.0;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				Ar_r+=Ar[n]*r[n];
				Ar_Ar+=Ar[n]*Ar[n];
				y_Ar+=y[n]*Ar[n];
				Ar_y+=Ar[n]*y[n];
			}

			zeta=v*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);
			eta=-y_Ar*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);

			//////vの計算//////////////
			v=0.0;
			for(int n=0;n<pn;n++) v+=zeta*Ar[n]*r[n];

			//////pの計算//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////Xの計算////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////yの計算////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

			//////rの計算////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//誤差評価
			rr=0;
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
	
			
		}
	}
	*countN=count;//反復回数を渡す
	//cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	//delete [] r;
	delete [] Ar;
	
	delete [] P;
	delete [] y;
	c.close();
}

//疎行列用の対角スケーリングつきbicgstab2法
void BiCGStab2_method_with_D_scale_for_sparse(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double *val,int *ind,int *ptr)
{
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double pk=1;
	double E=1;//誤差
	double alp,beta,rr,w,ita;

	double *P=new double [pn];	//探索ベクトル
	double *Pj=new double [pn];
	double *rj=new double [pn];
	double *AP=new double [pn];
	double *AtPj=new double [pn];
	double *e=new double [pn];
	double *Ae=new double [pn];
	double *y=new double [pn];
	double *u=new double [pn];
	double *W=new double [pn];
	double *Z=new double [pn];
	double *e_pre=new double [pn];
	double *D=new double [pn];			//D^(-1)
	beta=0;
	w=0;
	ita=0;

	for(int n=0;n<pn;n++)
	{
		P[n]=r[n];//初期化
		rj[n]=r[n];
		Pj[n]=rj[n];
		AP[n]=0;
		y[n]=0;
		u[n]=0;
		W[n]=0;
		Z[n]=0;
		e[n]=0;
		e_pre[n]=0;//1ステップ前のe[]
		for(int m=ptr[n];m<ptr[n+1];m++)
		{
			if(ind[m]==n)
			{
				if(val[m]==0) cout<<"RR"<<endl;
				D[n]=1.0/val[m];
			}
		}
		
	}

	ofstream c("convergenceP.dat");
	double rr0=0;
	for(int n=0;n<pn;n++) rr0+=r[n]*r[n];
	cout<<"BiCGstab2法:未知数="<<pn<<" ---";
	while(E>EP)// EP=CON->get_CGep();//収束判定(convergence test)
	{
		count++;
		for(int n=0;n<pn;n++)
		{
			//P[n]=r[n]+beta*(P[n]-w*AP[n]);
			P[n]=r[n]+beta*(P[n]-u[n]);
		} 
		////pk(rとrjの内積)を求める
		pk=0;
		for(int n=0;n<pn;n++) pk+=r[n]*rj[n];
		//////////////alpを求める
	//	#pragma omp parallel for
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]]*D[ind[m]];	//new
			/*for(int m=0;m<pn;m++)
			{
				AP[n]+=A[n][m]*P[m]*D[m];
			}*/
		}
		double APrj=0;
		for(int n=0;n<pn;n++)  APrj+=rj[n]*AP[n];
		alp=pk/APrj;
		//cout<<"alp="<<alp<<" APrj="<<APrj<<endl;
		//////////////////////

		for(int n=0;n<pn;n++) y[n]=e[n]-r[n]-alp*W[n]+alp*AP[n];

		for(int n=0;n<pn;n++)
		{
			e_pre[n]=e[n];//変更前の値を記憶
			e[n]=r[n]-alp*AP[n];
		}

	//	#pragma omp parallel for
		for(int n=0;n<pn;n++)
		{
			Ae[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) Ae[n]+=val[m]*e[ind[m]]*D[ind[m]];//new
			//for(int m=0;m<pn;m++) Ae[n]+=A[n][m]*e[m]*D[m];
		}

		if(count%2!=0)
		{
			double e_dot_ae  = 0;
			double ae_dot_ae = 0;
			for(int n=0;n<pn;n++)
			{
				e_dot_ae+=e[n]*Ae[n];
				ae_dot_ae+=Ae[n]*Ae[n];
			}
			w = e_dot_ae / ae_dot_ae;
			ita=0;
		}
		else
		{
			double e_dot_ae  = 0;
			double y_dot_y=0;
			double y_dot_e=0;
			double Ae_dot_y=0;
			double Ae_dot_Ae=0;
			for(int n=0;n<pn;n++) e_dot_ae+=e[n]*Ae[n];
			for(int n=0;n<pn;n++) y_dot_y+=y[n]*y[n];
			for(int n=0;n<pn;n++) y_dot_e+=y[n]*e[n];
			for(int n=0;n<pn;n++) Ae_dot_y+=Ae[n]*y[n];
			for(int n=0;n<pn;n++) Ae_dot_Ae+=Ae[n]*Ae[n];

			double co=Ae_dot_Ae*y_dot_y-Ae_dot_y*Ae_dot_y;

			w=y_dot_y*e_dot_ae-y_dot_e*Ae_dot_y;
			w/=co;

			ita=Ae_dot_Ae*y_dot_e-Ae_dot_y*e_dot_ae;
			ita/=co;

		}
		

		for(int n=0;n<pn;n++) 
		{
			u[n]=w*AP[n]+ita*(e_pre[n]-r[n]+beta*u[n]);
			Z[n]=w*r[n]+ita*Z[n]-alp*u[n];
			X[n]+=alp*P[n]*D[n]+Z[n]*D[n];
			r[n]=e[n]-ita*y[n]-w*Ae[n];
		}

		///beta
		beta=0;
		for(int n=0;n<pn;n++) beta+=r[n]*rj[n];
		beta/=w*APrj;
		
		//cout<<"beta="<<beta<<endl;
		//////////////////////

		//W[n]
		for(int n=0;n<pn;n++) W[n]=Ae[n]+beta*AP[n];
		
		//////////////////誤差
		rr=0;
		for(int n=0;n<pn;n++) rr+=r[n]*r[n];
		E=sqrt(rr/rr0);
		c<<count<<" "<<E<<endl;
		//cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
	}
	
	*countN=count;//反復回数を渡す
	delete [] Pj;
	delete [] rj;
	delete [] AP;
	delete [] AtPj;
	delete [] P;
	delete [] e;
	delete [] Ae;
	delete [] y;
	delete [] u;
	delete [] W;
	delete [] Z;
	delete [] e_pre;

	delete [] D;
	c.close();
}

//最小二乗法による表面速度定義関数
void reset_surface_velocity_by_WLSM(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number)
{
	//単純に内部計算点に関して圧力の方程式を解いただけでは、表面計算点の速度発散はゼロにならない。よって、表面の速度は汚い。そこで、これを内部粒子に関する最小二乗法により補完する

	double le=CON->get_distancebp();
	
	int d=CON->get_dimention();
	int N0=0;					//係数行列の元
	int order=1;//CON->get_divU_order();				//近似曲面のｵｰﾀﾞｰ。 1=線形 2=二次
    
	//係数行列の大きさの決定
	if(d==2)
	{
		if(order==1) N0=3;
		else if(order==2) N0=6;
		else if(order==3) N0=10;
		N0=9+1;	//最大値
	}
	else if(d==3)
	{
		if(order==1) N0=4;
		else if(order==2) N0=10;
		else if(order==3) N0=20;
		N0=20;
	}
	////////////////////////////////

	double *matrix=new double [N0*N0];	//N×Nの係数行列
	double *matrix2=new double [N0*N0];	//matrixのバックアップ
	double *B1=new double [N0];			//Nの解行列
	double *B2=new double [N0];			//Nの解行列
	double *B3=new double [N0];			//Nの解行列
	double **MAT= new double*[N0];		//N×Nの係数行列(配列は2次元)
	for(int n=0;n<N0;n++) MAT[n]=new double[N0];
	double *base=new double[N0];			//基底ベクトル格納

	if(d==2)				//二次元
	{
		for(int i=0;i<fluid_number;i++)
		{
			int N=N0;
			int jnb=0;			//周辺内部粒子数
			int J=i;
			double dis0=100;
			if(PART[i].surface==ON) 
			{
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					if(PART[j].type==FLUID && PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
						double dis=sqrt(X*X+Y*Y);
						if(dis<dis0)
						{
							dis0=dis; J=j;
						}
					}
				}
				jnb=0;
				//J=i;
				for(int k=0;k<PART[J].N;k++)
				{
					int j=PART[J].NEI[k]; 
					if(PART[j].surface==OFF) jnb++;
					//if(PART[j].type==FLUID && PART[j].surface==OFF) jnb++;
				}
			}


			if(jnb>=3 &&  PART[i].surface==ON) //周辺粒子がこれより少ない場合は速度補間を行わない
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;	//初期化
				for(int n=0;n<N0;n++) {B1[n]=0;B2[n]=0;}
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//初期化

				//if(PART[i].N>=15) {order=3; N=10;}
				if(jnb>=15) {order=2; N=6;}		//表面だから２次近似で十分な気はする
				else if(jnb>=6) {order=2; N=6;}
				else {order=1; N=3;}
				double R=PART[J].L*CON->get_re3();
				for(int k=0;k<PART[J].N3;k++)
				{
					int j=PART[J].NEI3[k];
					//if(PART[j].type==FLUID && PART[j].surface==OFF)
					if(PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[J].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[J].r[A_Y]);
						double U=(PART[j].u[A_X]);
						double V=(PART[j].u[A_Y]);
						double dis=sqrt(X*X+Y*Y);
					
						double w=kernel_in_WLSM(dis,R);	//スプライン
						if(dis>R) w=0;

						if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//基底ベクトル
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y;  base[5]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y;  base[9]=1;}

						//行列作成
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B1[n]+=base[n]*w*U;
							B2[n]+=base[n]*w*V;
						}
					}
				}
				if(J!=i)		//粒子Jの寄与。もし粒子Jがiに等しいなら考慮する必要はない
				{
					MAT[N-1][N-1]+=1;		//一番右下の配列に自分自身の寄与を加える
					B1[N-1]+=PART[J].u[A_X];	//粒子iの寄与。
					B2[N-1]+=PART[J].u[A_Y];	//粒子iの寄与。
					//これ以外はX=0になるので加算する必要はない*/
				}
				
				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//値をmatrixに転送

				for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//値を保存

				gauss(matrix,B1,N);//ガウスの消去法で解く

				for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
				gauss(matrix,B2,N);//ガウスの消去法で解く

				//速度の値を補間
				double X=(PART[i].r[A_X]-PART[J].r[A_X]);
				double Y=(PART[i].r[A_Y]-PART[J].r[A_Y]);
				if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//基底ベクトル
				else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y;  base[5]=1;}
				else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y;  base[9]=1;}

				for(int D=0;D<d;D++) PART[i].u[D]=0;
				for(int n=0;n<N;n++)
				{
					PART[i].u[A_X]+=B1[n]*base[n];
					PART[i].u[A_Y]+=B2[n]*base[n];
				}
			}
		}
	}
	else if(d==3)				//3次元
	{
		for(int i=0;i<fluid_number;i++)
		{
			int N=N0;
			int jnb=0;			//周辺内部粒子数
			int J=i;			//最近接内部流体粒子
			double dis0=100;
			if(PART[i].surface==ON) 
			{
				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(PART[j].type==FLUID && PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[i].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
						double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
						double dis=sqrt(X*X+Y*Y+Z*Z);
						if(dis<dis0) {dis0=dis; J=j;}
					}
				}
				jnb=0;
				for(int k=0;k<PART[J].N;k++)
				{
					int j=PART[J].NEI[k]; 
					if(PART[j].surface==OFF)jnb++;
					//if(PART[j].type==FLUID && PART[j].surface==OFF) jnb++;
				}
			}


			if(jnb>=4 &&  PART[i].surface==ON) //周辺粒子がこれより少ない場合は速度補間を行わない
			{
				for(int n=0;n<N0*N0;n++) matrix[n]=0;	//初期化
				for(int n=0;n<N0;n++) {B1[n]=0;B2[n]=0; B3[n]=0;}
				for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//初期化

				if(jnb>=25) {order=3; N=20;}
				else if(jnb>=10) {order=2; N=10;}
				else {order=1; N=4;}
				double R=PART[J].L*CON->get_re3();

				for(int k=0;k<PART[J].N3;k++)
				{
					int j=PART[J].NEI3[k];
					//if(PART[j].type==FLUID && PART[j].surface==OFF)
					if(PART[j].surface==OFF)
					{
						double X=(PART[j].r[A_X]-PART[J].r[A_X]);
						double Y=(PART[j].r[A_Y]-PART[J].r[A_Y]);
						double Z=(PART[j].r[A_Z]-PART[J].r[A_Z]);
						double U=(PART[j].u[A_X]);
						double V=(PART[j].u[A_Y]);
						double W=(PART[j].u[A_Z]);
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						double w=kernel_in_WLSM(dis,R);	//スプライン
						if(dis>R) w=0;								//スプライン

						if(order==1) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=1;}	//基底ベクトル
						else if(order==2) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=1;}
						else if(order==3) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=X*X*X; base[10]=Y*Y*Y; base[11]=Z*Z*Z; base[12]=X*X*Y; base[13]=X*Y*Y; base[14]=Y*Y*Z; base[15]=Y*Z*Z; base[16]=X*X*Z; base[17]=X*Z*Z; base[18]=X*Y*Z; base[19]=1;}
				
						//行列作成
						for(int n=0;n<N;n++)
						{
							for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
							B1[n]+=base[n]*w*U;
							B2[n]+=base[n]*w*V;
							B3[n]+=base[n]*w*W;
						}
					}
				}
				if(J!=i)		//粒子Jの寄与。もし粒子Jがiに等しいなら考慮する必要はない
				{
					MAT[N-1][N-1]+=1;		//一番右下の配列に自分自身の寄与を加える
					B1[N-1]+=PART[J].u[A_X];	//粒子iの寄与。
					B2[N-1]+=PART[J].u[A_Y];	//粒子iの寄与。
					B3[N-1]+=PART[J].u[A_Z];	//粒子iの寄与。
					//これ以外はX=0になるので加算する必要はない*/
				}

				for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//値をmatrixに転送

				for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//値を保存

				gauss(matrix,B1,N);//ガウスの消去法で解く

				for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
				gauss(matrix,B2,N);//ガウスの消去法で解く

				for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
				gauss(matrix,B3,N);//ガウスの消去法で解く

				//速度の値を補間
				double X=(PART[i].r[A_X]-PART[J].r[A_X]);
				double Y=(PART[i].r[A_Y]-PART[J].r[A_Y]);
				double Z=(PART[i].r[A_Z]-PART[J].r[A_Z]);
				if(order==1) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=1;}	//基底ベクトル
				else if(order==2) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=1;}
				else if(order==3) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=X*X*X; base[10]=Y*Y*Y; base[11]=Z*Z*Z; base[12]=X*X*Y; base[13]=X*Y*Y; base[14]=Y*Y*Z; base[15]=Y*Z*Z; base[16]=X*X*Z; base[17]=X*Z*Z; base[18]=X*Y*Z; base[19]=1;}
				
				for(int D=0;D<d;D++) PART[i].u[D]=0;
				for(int n=0;n<N;n++)
				{
					PART[i].u[A_X]+=B1[n]*base[n];
					PART[i].u[A_Y]+=B2[n]*base[n];
					PART[i].u[A_Z]+=B3[n]*base[n];
				}
			}
		}
	}
	
	delete [] matrix;
	delete [] matrix2;
	delete [] B1;
	delete [] B2;
	delete [] B3;
	for(int n=0;n<N0;n++) delete [] MAT[n];
	delete [] MAT;
	delete [] base;
}

