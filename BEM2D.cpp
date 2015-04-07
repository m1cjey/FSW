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

void BEM2D(mpsconfig *CON,vector<point2D> &s_NODE,vector<element2D> &s_ELEM,int *s_node_num,int *s_elem_num,vector<mpsparticle> &PART,int fluid_number,int particle_number,double **F,int t)
{

	//静的計算点の配置を得る
	if(t==1) set_BEM_static_model(CON,s_NODE,s_ELEM,s_node_num,s_elem_num);
	//cout<<"固定節点数="<<s_NODE.size()<<endl;

	//動的計算点の配置を得る
	vector<point2D> dy_NODE;
	vector<element2D> dy_ELEM;
	int d_node_num,d_elem_num;
	set_BEM2D_dynaic_model_for_CONSTANT(CON,PART,dy_NODE,dy_ELEM,&d_node_num,&d_elem_num, particle_number, fluid_number);

	//両者を合体
	vector<point2D> NODE;
	vector<element2D> ELEM;
	vector<REGION> region;
	
	couple_NODE_and_ELEM(CON,NODE,ELEM,s_NODE,s_ELEM,dy_NODE,dy_ELEM,region);

	//main関数実行
	BEM_main_for_CONSTANT(CON,NODE,ELEM,PART, fluid_number,F,region);

}

//一定要素用のBEM関数
void BEM_main_for_CONSTANT(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,vector<mpsparticle> &PART,int fluid_number,double **F,vector<REGION> &region)
{
	int node_num=int (NODE.size());		//節点数
	int elemnum=int (ELEM.size());		//要素数
	int region_num=(int)(region.size());//領域の数
	int uk=0;								//未知数
	int count;

	cout<<"領域数＝"<<region_num<<endl;
	
	//未知数計算
	for(int i=0;i<node_num;i++)
	{
		if(NODE[i].boundary_condition==BOTH) uk+=2;	//BOTHは多媒質境界上の節点であり、ﾎﾟﾃﾝｼｬﾙ、法線微分の両方が未知であることを示す
		else uk++;									//こちらは通常。DiricだろうとNeumnだろうと、片方は未知なので++;
	}
	cout<<"未知数="<<uk<<endl;

	int *Nid=new int [uk];							//i番目の未知数はi番目の計算点であることを示す
	int *bd_type=new int [uk];
	int *BOTH_column=new int [node_num];
	count=0;
	
	for(int i=0;i<elemnum;i++)
	{
		if(ELEM[i].boundary_condition==BOTH)
		{
			Nid[i]=i;		//i番目の未知数はi番目の計算点であることを示す
			bd_type[i]=Neumn;//BOTHの場合、未知数はﾎﾟﾃﾝｼｬﾙ、法線微分の順に定義する。よってbd_typeはこのようにしておく
			Nid[elemnum+count]=i;
			bd_type[elemnum+count]=Diric;
			BOTH_column[i]=elemnum+count;//i番目の計算点がBOTHのとき、法線微分を表す未知数はBOTH_column[i]番目の未知数である。
			count++;
		}
		else
		{
			Nid[i]=i;				//i番目の未知数はi番目の計算点であることを示す
			bd_type[i]=ELEM[i].boundary_condition;
			BOTH_column[i]=-1;		//BOTHでないならダミーとして－１を格納
		}
	}
	
	
	int gauss_N=4;				//ガウス積分における評価点数 //3,4,8
	int NGS[32];				//たとえばGauss積分点が4のときはfor(int i=NGS[4];i<NGS[4]+4;i++) ep1[i]=~~~~という使いかたをする
	double ep1[32];				//ガウス積分における局所座標格納
	double ep_ln[32];			//自然対数に関するガウス積分における局所座標格納
	double w1[32];				//ガウス積分における重み格納
	double w_ln[32];			//自然対数に関するガウス積分における重み格納

	//ガウス積分の準備
	set_gauss_parameter(NGS, ep1,w1,ep_ln,w_ln);

	//matrix作成
	double **matrixC=new double*[uk];
	for(int n=0;n<uk;n++) matrixC[n]=new double[uk];
	double **matrixC2=new double*[uk];
	for(int n=0;n<uk;n++) matrixC2[n]=new double[uk];
	double **H=new double*[uk];
	for(int n=0;n<uk;n++) H[n]=new double[uk];
	double **G=new double*[uk];
	for(int n=0;n<uk;n++) G[n]=new double[uk];
	double *Bmatrix=new double[uk];//解行列
	
	for(int n=0;n<uk;n++)
	{
		for(int m=0;m<uk;m++)
		{
			matrixC[n][m]=0;				//初期化
			matrixC2[n][m]=0;
			H[n][m]=0;
			G[n][m]=0;
		}
		Bmatrix[n]=0;
	}
	
	for(int ID=0;ID<1;ID++)				//まずは第１領域のみ計算
	{
		for(int n1=region[ID].start;n1<region[ID].end;n1++)
		{	
			int n=n1;
			int n2=BOTH_column[n1];
				
			double SR=0.5*ELEM[n1].L;//要素長さの半分
					
			for(int m1=region[ID].start;m1<region[ID].end;m1++)
			{
				int m=m1;
				int m2=BOTH_column[m1];
					
				int type=0;//type=0なら通常　1なら、その要素は節点nを含む特異積分
				if(n1==m1) type=1;
					
				if(type==0)
				{
					double w;
					double ep;		//評価点の局所座標
					for(int i=NGS[gauss_N];i<NGS[gauss_N]+gauss_N;i++)
					{
						ep=ep1[i];
						w=w1[i];
							
						set_GHmatrix_for_constant(ELEM, ep,  H, G, w, type,n, m,n1,m1);
					}
					if(m2!=-1)
					{
						G[n][m2]=G[n][m];
						G[n][m]=0;
					}
				}
			}
			H[n][n]=PI;
			G[n][n]=2*SR*(1-log(SR));
				

			if(n2!=-1)
			{
				G[n][n2]=G[n][n];
				G[n][n]=0;
			}
		}
	}
	//法線ベクトルを反転
	if(region_num>1) for(int i=0;i<elemnum;i++) for(int D=0;D<2;D++) ELEM[i].direct[D]*=-1.0;
	
	for(int ID=1;ID<region_num;ID++)				//第2領域以降を計算
	{
		for(int n1=region[ID].start;n1<region[ID].end;n1++)
		{	
			int n=BOTH_column[n1];					//計算点n1の、法線微分の未知変数はBOTH_column[n1]番目の未知数
			
			double SR=0.5*ELEM[n1].L;//要素長さの半分
			for(int m1=region[ID].start;m1<region[ID].end;m1++)
			{
				int m=m1;
				int m2=BOTH_column[m1];
				int type=0;//type=0なら通常　1なら、その要素は節点nを含む特異積分
				if(n1==m1) type=1;
					
				if(type==0)
				{
					double w;
					double ep;		//評価点の局所座標
						
					for(int i=NGS[gauss_N];i<NGS[gauss_N]+gauss_N;i++)
					{
						ep=ep1[i];
						w=w1[i];
							
						set_GHmatrix_for_constant(ELEM, ep,  H, G, w, type,n, m,n1,m1);
					}	
					G[n][m2]=-G[n][m]/80;
					//G[n][m2]=-G[n][m]/1;
					G[n][m]=0;
				}
			}
			H[n][n1]=PI;
			G[n][n]=-2*SR*(1-log(SR));
		}
	}
	//法線ベクトルを反転
	if(region_num>1) for(int i=0;i<elemnum;i++) for(int D=0;D<2;D++) ELEM[i].direct[D]*=-1.0;

	//実際に解く係数行列作成
	for(int n=0;n<elemnum;n++)	
	{
		if(ELEM[n].boundary_condition==Diric)
		{
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n]=-G[m][n];
				matrixC2[m][n]=-H[m][n];
			}
		}
		else if(ELEM[n].boundary_condition==Neumn)
		{
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n]=H[m][n];
				matrixC2[m][n]=G[m][n];
			}
		}
		else if(ELEM[n].boundary_condition==BOTH)
		{
			int n2=BOTH_column[n];
			
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n2]=-G[m][n2];
				if(H[m][n2]!=0) cout<<"H "<<m<<" "<<n2<<" "<<H[m][n2]<<endl;
			}
				
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n]=H[m][n];
				if(G[m][n]!=0) cout<<"G "<<m<<" "<<n<<" "<<G[m][n]<<endl;
			}
		}
	}

	//解行列作成
	for(int n=0;n<uk;n++)
	{
		double val=0;//解行列の値
		for(int m=0;m<uk;m++)
		{
			int elem=Nid[m];			//未知数nに対応する要素番号
			int node=ELEM[elem].node[0];	//要素elemに対応するのはELEM[elem].node[0]に格納している
			if(bd_type[m]==Diric) val+=matrixC2[n][m]*NODE[node].potential;
			else if(bd_type[m]==Neumn) val+=matrixC2[n][m]*NODE[node].slop1;
		}
		Bmatrix[n]=val;
	}
	/////matrix作成完了

	cout<<"行列作成完了 ガウスの消去法---";


	gauss(matrixC,Bmatrix,uk);//答えはBmatrixに格納



	cout<<"ok"<<endl;

	//答え格納
	for(int n=0;n<uk;n++)
	{
		int elem=Nid[n];		//未知数nに対応する要素番号
		int node=ELEM[elem].node[0];
		if(bd_type[n]==Neumn) NODE[node].potential=Bmatrix[n];
		else if(bd_type[n]==Diric) NODE[node].slop1=Bmatrix[n]; 
	}

	//答え確認
	ofstream fp3("V.dat");
	ofstream fp4("slop.dat");
	for(int n=0;n<elemnum;n++)
	{
		int node=ELEM[n].node[0];
		fp3<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Y]<<" "<<NODE[node].potential<<endl;
		fp4<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Y]<<" "<<NODE[node].slop1<<endl;
	}
	fp3.close();
	fp4.close();

	

	ofstream fp5("F.dat");
	double ep0=8.854e-12;
	double times=1e-3;
	double le=CON->get_distancebp();
	for(int n=0;n<elemnum;n++)
	{
		if(ELEM[n].material==FLUID)
		{
			int node=ELEM[n].node[0];//対応する節点番号
			int i=NODE[node].particle;//対応する粒子番号
			double E=NODE[node].slop1;
			double Fs=0.5*ep0*E*E;	//応力
			double Fn=Fs*ELEM[n].L;	//力[N]
			for(int D=0;D<2;D++) F[D][i]=-Fn*ELEM[n].direct[D];
			fp5<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<F[A_X][i]*times/le<<" "<<F[A_Y][i]*times/le<<endl;
			//fp5<<NODE[node].r[A_X]<<" "<<NODE[node].r[A_Y]<<" "<<F[A_X][i]*times/le<<" "<<F[A_Y][i]*times/le<<endl;
		//	fp5<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Y]<<" "<<Fs<<endl;
		}
	}
	fp5.close();

	////ｽﾑｰｼﾞﾝｸﾞ
	//smoothingF3D(CON,PART,fluid_number,F,t);

	if(CON->get_dir_for_P()==2 ||CON->get_dir_for_P()==3 )
    {
		ofstream bb("electromagnetic_P.dat");
		for(int i=0;i<fluid_number;i++)
		{
			double fs=0;//表面力
			if(PART[i].surface==ON)//内部流体の場合はfs=0とする
			{
				double Fn=sqrt(F[A_X][i]*F[A_X][i]+F[A_Y][i]*F[A_Y][i]);
				fs=Fn/le;
				for(int D=0;D<3;D++) F[D][i]=0;//圧力デﾞｨﾘｸﾚとして電磁力を使用するのでここでは初期化
			}
			bb<<-fs<<endl;
        }
		bb.close();
	}


	for(int n=0;n<node_num;n++) delete[] matrixC[n];
	delete [] matrixC;
	for(int n=0;n<node_num;n++) delete[] matrixC2[n];
	delete [] matrixC2;
	for(int n=0;n<node_num;n++) delete[] H[n];
	delete [] H;
	for(int n=0;n<node_num;n++) delete[] G[n];
	delete [] G;
	delete [] Bmatrix;

	delete [] Nid;
	delete [] bd_type;
	delete [] BOTH_column;
}

//Gauss積分の準備関数
void set_gauss_parameter(int *NGS,double *ep1,double *w1,double *ep_ln,double *w_ln)
{
	//ガウス積分の準備
	for(int i=0;i<32;i++) NGS[i]=0;
	NGS[3]=0; NGS[4]=3; NGS[8]=7;	//たとえばGauss積分点が4のときはfor(int i=NGS[4];i<NGS[4]+4;i++) ep1[i]=~~~~という使いかたをする

	//gauss_N=3のとき
	ep1[0]=-0.774596669241483;
	ep1[1]=0;
	ep1[2]=0.774596669241483;
	w1[0]=0.555555555555556;
	w1[1]=0.888888888888889;
	w1[2]=0.555555555555556;
	

	ep_ln[0]=0.6389079308732540e-1;
	ep_ln[1]=0.3689970637156187;
	ep_ln[2]=0.7668803039389414;
	w_ln[0]=0.5134045522323633;
	w_ln[1]=0.3919800412014875;
	w_ln[2]=0.9461540656614912e-1;

	//gauss_N=4のとき
	ep1[3]=-0.861136311594053;
	ep1[4]=0.861136311594053;
	ep1[5]=-0.339981043584856;
	ep1[6]=0.339981043584856;
	w1[3]=0.347854845137454;
	w1[4]=0.347854845137454;
	w1[5]=0.652145154862546;
	w1[6]=0.652145154862546;

	ep_ln[3]=0.4144848019938322e-1;
	ep_ln[4]=0.2452749143206022;
	ep_ln[5]=0.5561654535602758;
	ep_ln[6]=0.8489823945329851;
	w_ln[3]=0.3834640681451351;
	w_ln[4]=0.3868753177747626;
	w_ln[5]=0.1904351269501424;
	w_ln[6]=0.3922548712995983e-1;

	//gauss_N=8のとき
	ep1[7]=-0.960289856497536;
	ep1[8]=0.960289856497536;
	ep1[9]=-0.796666477413627;
	ep1[10]=0.796666477413627;
	ep1[11]=-0.525532409916329;
	ep1[12]=0.525532409916329;
	ep1[13]=-0.183434642495650;
	ep1[14]=0.183434642495650;
	w1[7]=0.101228536290376;
	w1[8]=0.101228536290376;
	w1[9]=0.222381034453374;
	w1[10]=0.222381034453374;
	w1[11]=0.313706645877887;
	w1[12]=0.313706645877887;
	w1[13]=0.362683783378362;
	w1[14]=0.362683783378362;

	ep_ln[7]=0.1332024416089246e-1;
	ep_ln[8]=0.7975042901389493e-1;
	ep_ln[9]=0.1978710293261880;
	ep_ln[10]=0.3541539943519094;
	ep_ln[11]=0.5294585752349172;
	ep_ln[12]=0.7018145299390999;
	ep_ln[13]=0.8493793204411066;
	ep_ln[14]=0.9533264500563597;
	w_ln[7]=0.1644166047280028;
	w_ln[8]=0.2375256100233060;
	w_ln[9]=0.2268419844319191;
	w_ln[10]=0.1757540790060702;
	w_ln[11]=0.1129240302467590;
	w_ln[12]=0.5787221071778207e-1;
	w_ln[13]=0.2097907374213297e-1;
	w_ln[14]=0.3686407104027619e-1;
}

//一定要素のためのHマトリクス計算関数
void set_GHmatrix_for_constant(vector<element2D> &ELEM,double ep,double **H,double **G,double w,int type,int n,int m,int n1,int m1)
{
	//n1:ソース点要素番号 n:n1に該当する行列要素番号
	//m1:相手要素番号 m:m1に該当する行列要素番号
	//type=0なら通常 1なら計算しない。
	//L:ソース要素長さ
	//Xms,Yms:相手要素の中点

	double Xs=ELEM[n1].r[A_X];//ソース点座標(要素の中点)
	double Ys=ELEM[n1].r[A_Y];

	double Xms=ELEM[m1].r[A_X];//相手要素の中点
	double Yms=ELEM[m1].r[A_Y];

	double SL=0.5*ELEM[m1].L;//相手要素長さの半分
	double nx=ELEM[m1].direct[A_X];					//相手要素の外向き単位法線ベクトル
	double ny=ELEM[m1].direct[A_Y];
	double AX=-ny*SL;		//要素の中点から終点へ向かうベクトルのX成分(nyの定義を考えればこう記述できることがわかる)
	double AY=nx*SL;		//要素の中点から終点へ向かうベクトルのY成分(nxの定義を考えればこう記述できることがわかる)
	
	//if(type==0)//通常
	{
		double Gx=Xms+AX*ep;//ガウス積分評価点のX座標格納 
		double Gy=Yms+AY*ep;//ガウス積分評価点のY座標格納

		double RA=sqrt((Xs-Gx)*(Xs-Gx)+(Ys-Gy)*(Ys-Gy));//ソース点と評価点の距離
		double RD1=(Gx-Xs)/RA;		//ソース点から評価点へ向かう単位ベクトル
		double RD2=(Gy-Ys)/RA;
		double RDN=RD1*nx+RD2*ny;	//上記ベクトルと法線ベクトルとの内積

		H[n][m]+=-RDN*w*SL/RA;
		G[n][m]+=log(1/RA)*w*SL;
	}
}

//ガウスの消去法 解は最終的にBのなかへ
void gauss(double **matrix,double *B,int N)
{
	for(int k=0;k<N;k++)
	{
		double akk=matrix[k][k];//対角成分

		for(int i=0;i<N;i++)
		{
			if(i!=k)
			{
				double A=matrix[i][k]/akk;
				//for(int j=0;j<N;j++)
				for(int j=k;j<N;j++)
				{					
					matrix[i][j]-=A*matrix[k][j];
				}
				B[i]-=A*B[k];
			}
		}
	}
	for(int k=0;k<N;k++)
	{
		B[k]/=matrix[k][k];
		if(fabs(matrix[k][k])<1e-4) cout<<"要ピボット matrix["<<k<<"]["<<k<<"]="<<matrix[k][k]<<endl;//これがでたらピボット必要
	}
}

//電磁力スムージング関数
void smoothingF3D(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *F[3],int t)
{
	double le=CON->get_distancebp();
    double ep=8.854e-12;///真空の誘電率。
    double *newF[3];
    for(int D=0;D<3;D++) newF[D]=new double [fluid_number];

	if(CON->get_FEM_smn()>0)
	{
		for(int n=0;n<CON->get_FEM_smn();n++)
		{
		    for(int i=0;i<fluid_number;i++) 
		    {  
		        for(int D=0;D<3;D++) newF[D][i]=F[D][i];
				int num=1; //自分自身をｶｳﾝﾄするから1
				for(int k=0;k<PART[i].N;k++)
				{       
					int j=PART[i].NEI[k];
					if(PART[j].type==FLUID)
					{
						num++;
						for(int D=0;D<3;D++) newF[D][i]+=F[D][j];
					}
				}
				for(int D=0;D<3;D++) newF[D][i]/=num;
		    } 
		    for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) F[D][i]=newF[D][i];
		}
	}
	else if(CON->get_FEM_smn()<0)//表面のみでスムージング
	{
		int N=-1*CON->get_FEM_smn();
		for(int n=0;n<N;n++)
		{
			for(int i=0;i<fluid_number;i++) 
			{  
			    for(int D=0;D<3;D++) newF[D][i]=F[D][i];
				if(PART[i].surface==ON)
				{
					int num=1; //自分自身をｶｳﾝﾄするから1
					for(int k=0;k<PART[i].N;k++)
					{       
						int j=PART[i].NEI[k];
						if(PART[j].surface==ON && PART[j].type==FLUID)
						{
							num++;
							for(int D=0;D<3;D++) newF[D][i]+=F[D][j];
						}
					}
					for(int D=0;D<3;D++) newF[D][i]/=num;
				}
			} 
			for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) F[D][i]=newF[D][i];
		}
	}

    for(int D=0;D<3;D++) delete [] newF[D];
    /////////////////////////////////////*/
    
    
    ////電磁力をプロット
	ofstream fp("F.dat");
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	double mass=CON->get_particle_mass();
	double times=CON->get_times()/mass*le*le*CON->get_FEMtimes();
    for(int i=0;i<fluid_number;i++)
    {
		//if(PART[i].r[A_Y]>-le*0.5&& PART[i].r[A_Y]<+le*0.5)
		{
			fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<F[A_X][i]*times<<"\t"<<F[A_Y][i]*times<<endl;
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
		}
    }
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON->get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力
	fp.close();/////

	int BOnum=0;
	for(int i=0;i<fluid_number;i++) if(PART[i].surface==ON) BOnum++;

	ofstream fout2("Lorentz.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	fout2 << "dim1=" << fluid_number <<endl;
	//fout2 << "dim1=" << BOnum <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Lorentz filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Lorentz filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Lorentz filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Lorentz filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Lorentz filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Lorentz filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("Lorentz");
	fout<<"e-x e-y e-z x y z"<<endl;
	for(int i=0;i<fluid_number;i++)
    {
		//if(PART[i].surface==ON)
		{
			fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		}
	}
	fout.close();

	int flag=0;
	if(CON->get_EM_interval()>1 &&t==1) flag=ON;
	//else if(t==1 || t%10==0) flag=ON;


	///////////////////////////////
	if(flag==ON)
	{
		char filename[25];
		sprintf_s(filename,"Lorentz%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		fout2 << "dim1=" << fluid_number <<endl;
		//fout2 << "dim1=" << BOnum <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Lorentz"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Lorentz"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Lorentz"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Lorentz"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Lorentz"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Lorentz"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[25];
		sprintf_s(filename2,"Lorentz%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=0;i<fluid_number;i++)
		{
			//if(PART[i].surface==ON)
			{
				fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		fout.close();
	}

	int BOnumsec=0;
	for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]>sin(PI/24)*PART[i].r[A_X]-0.5*le && PART[i].r[A_Y]<sin(PI/24)*PART[i].r[A_X]+0.5*le) BOnumsec++;

	ofstream foutz("Lorentzsec.fld");
	foutz << "# AVS field file" << endl;
	foutz << "ndim=1" << endl;
	//foutz << "dim1=" << fluid_number <<endl;
	foutz << "dim1=" << BOnumsec <<endl;
	foutz << "nspace=3" << endl;
	foutz << "veclen=3" << endl;
	foutz << "data=float" << endl;
	foutz << "field=irregular" << endl;
	foutz << "label=e-x e-y e-z" << endl << endl;
	foutz << "variable 1 file=./Lorentzsec filetype=ascii skip=1 offset=0 stride=6" << endl;
	foutz << "variable 2 file=./Lorentzsec filetype=ascii skip=1 offset=1 stride=6" << endl;
	foutz << "variable 3 file=./Lorentzsec filetype=ascii skip=1 offset=2 stride=6" << endl;
	foutz << "coord    1 file=./Lorentzsec filetype=ascii skip=1 offset=3 stride=6" << endl;
	foutz << "coord    2 file=./Lorentzsec filetype=ascii skip=1 offset=4 stride=6" << endl;
	foutz << "coord    3 file=./Lorentzsec filetype=ascii skip=1 offset=5 stride=6" << endl;
	foutz.close();

	ofstream fout2z("Lorentzsec");
	fout2z<<"e-x e-y e-z x y z"<<endl;
	for(int i=0;i<fluid_number;i++)
    {
		//if(PART[i].surface==ON)
		//if(PART[i].r[A_Z]>0.14125-0.5*le && PART[i].r[A_Z]<0.14125+0.5*le)
		if(PART[i].r[A_Y]>sin(PI/24)*PART[i].r[A_X]-0.5*le && PART[i].r[A_Y]<sin(PI/24)*PART[i].r[A_X]+0.5*le)
		{
			fout2z<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		}
	}
	fout2z.close();
}