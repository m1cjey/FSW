#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
//#include"define.h"	//#define 格納
//#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
//#include"PART.h"		//class PART定義
#include"BEMclass.h"	//BEM2D関係のclass 定義
//#include"FEM3Dclass.h"	//FEM3D関係のclass 定義
#include<omp.h>
#include<vector>
#include"function.h"

void set_BEM3D_static_model(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,int *s_node_num,int *s_elem_num);
void set_gauss_parameter3D(int *NGS,double ep1[14][3],double *w1,double *ep_ln,double *w_ln);
void set_GHmatrix3D_for_constant(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,double ep[3],double **H,double **G,double w,int type,int n,int m,int n1,int m1);
double calc_Gnn_for_CONSTANT(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,int n1);


void BEM3D(mpsconfig *CON,vector<BEMpoint3D> &s_NODE,vector<BEMelement3D> &s_ELEM,int *s_node_num,int *s_elem_num,vector<mpsparticle> &PART,int fluid_number,int particle_number,double **F,int t)
{
	// t=1のときに、固定要素を作成してs_NODE,s_ELEMに格納する。節点数と要素数もs_node_numとs_elem_numに格納。以後、固定要素に関してはこれらを使い続ける

	//静的計算点の配置を得る
	if(t==1) set_BEM3D_static_model(CON,s_NODE,s_ELEM,s_node_num,s_elem_num);
	cout<<"固定節点数="<<s_NODE.size()<<"固定要素数="<<s_ELEM.size()<<endl;
	cout<<*s_elem_num<<endl;

	//動的計算点の配置を得る
	vector<BEMpoint3D> dy_NODE;
	vector<BEMelement3D> dy_ELEM;
	int d_node_num,d_elem_num;
	set_BEM3D_dynaic_model_for_CONSTANT(CON,PART,dy_NODE,dy_ELEM,&d_node_num,&d_elem_num, particle_number, fluid_number);

	//両者を合体
	vector<BEMpoint3D> NODE;
	vector<BEMelement3D> ELEM;
	vector<REGION> region;
	
	couple3D_NODE_and_ELEM(CON,NODE,ELEM,s_NODE,s_ELEM,dy_NODE,dy_ELEM,region);

	//main関数実行
	BEM3D_main_for_CONSTANT(CON,NODE,ELEM,PART, fluid_number,F,region);

}

//一定要素用のBEM関数
void BEM3D_main_for_CONSTANT(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,vector<mpsparticle> &PART,int fluid_number,double **F,vector<REGION> &region)
{
	int node_num=int (NODE.size());		//節点数
	int elemnum=int (ELEM.size());		//要素数
	int region_num=(int)(region.size());//領域の数
	int uk=0;								//未知数
	int count;

	cout<<"領域数＝"<<region_num<<endl;
	
	//未知数計算
	for(int i=0;i<elemnum;i++)
	{
		if(ELEM[i].boundary_condition==BOTH) uk+=2;	//BOTHは多媒質境界上の節点であり、ﾎﾟﾃﾝｼｬﾙ、法線微分の両方が未知であることを示す
		else uk++;									//こちらは通常。DiricだろうとNeumnだろうと、片方は未知なので++;
	}
	cout<<"未知数="<<uk<<endl;

	int *Nid=new int [uk];							//i番目の未知数はi番目の計算点であることを示す
	int *bd_type=new int [uk];
	int *BOTH_column=new int [elemnum];
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
			//cout<<i<<endl;
			Nid[i]=i;				//i番目の未知数はi番目の計算点であることを示す
			bd_type[i]=ELEM[i].boundary_condition;
			BOTH_column[i]=-1;		//BOTHでないならダミーとして－１を格納
			//cout<<i<<endl;
		}
	}
	
	int gauss_N=4;				//ガウス積分における評価点数 //3,4,8
	int NGS[14];				//たとえばGauss積分点が4のときはfor(int i=NGS[4];i<NGS[4]+4;i++) ep1[i]=~~~~という使いかたをする
	double ep1[14][3];			//ガウス積分における局所座標格納
	double ep_ln[14];			//自然対数に関するガウス積分における局所座標格納
	double w1[14];				//ガウス積分における重み格納
	double w_ln[14];			//自然対数に関するガウス積分における重み格納

	//ガウス積分の準備
	set_gauss_parameter3D(NGS, ep1,w1,ep_ln,w_ln);

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
	cout<<"行列作成開始"<<endl;
	for(int ID=0;ID<1;ID++)				//まずは第１領域のみ計算
	{
		for(int n1=region[ID].start;n1<region[ID].end;n1++)
		{	
			int n=n1;
			int n2=BOTH_column[n1];
			if(n==2036) cout<<n1<<" "<<ELEM[n1].material<<endl;
				
			double SR=0.5*ELEM[n1].S;//要素面積
					
			for(int m1=region[ID].start;m1<region[ID].end;m1++)
			{
				int m=m1;
				int m2=BOTH_column[m1];
					
				int type=0;//type=0なら通常　1なら、その要素は節点nを含む特異積分
				if(n1==m1) type=1;
					
				if(type==0)//n=mの場合も数値積分する
				{
					double w;
					double ep[3];		//評価点の局所座標
					for(int i=NGS[gauss_N];i<NGS[gauss_N]+gauss_N;i++)
					{
						for(int j=0;j<3;j++) ep[j]=ep1[i][j];
						w=w1[i];
							
						set_GHmatrix3D_for_constant(ELEM, NODE,ep,  H, G, w, type,n, m,n1,m1);
					}
					if(m2!=-1)//BOTHなら
					{
						G[n][m2]=G[n][m];
						G[n][m]=0;
					}
				}
			}
			H[n][n]=0.5;		//Hのみ対角項を書き換える
			G[n][n]=calc_Gnn_for_CONSTANT(ELEM,NODE, n1);
				

			if(n2!=-1)
			{
				G[n][n2]=G[n][n];
				G[n][n]=0;
			}
		}
	}
	//法線ベクトルを反転
	if(region_num>1) for(int i=0;i<elemnum;i++) for(int D=0;D<3;D++) ELEM[i].direct[D]*=-1.0;
	
	for(int ID=1;ID<region_num;ID++)				//第2領域以降を計算
	{
		for(int n1=region[ID].start;n1<region[ID].end;n1++)
		{	
			int n=BOTH_column[n1];					//計算点n1の、法線微分の未知変数はBOTH_column[n1]番目の未知数
			
			double SR=0.5*ELEM[n1].S;//要素面積
			for(int m1=region[ID].start;m1<region[ID].end;m1++)
			{
				int m=m1;
				int m2=BOTH_column[m1];
				int type=0;//type=0なら通常　1なら、その要素は節点nを含む特異積分
				if(n1==m1) type=1;
					
				if(type==0)//n=mの場合も数値積分する
				{
					double w;
					double ep[3];		//評価点の局所座標
						
					for(int i=NGS[gauss_N];i<NGS[gauss_N]+gauss_N;i++)
					{
						for(int j=0;j<3;j++) ep[j]=ep1[i][j];
						w=w1[i];
							
						set_GHmatrix3D_for_constant(ELEM, NODE,ep,  H, G, w, type,n, m,n1,m1);
					}
					G[n][m2]=-G[n][m]/80;
					//G[n][m2]=-G[n][m]/1;
					G[n][m]=0;
				}
			}
			H[n][n1]=0.5;
			G[n][n]=calc_Gnn_for_CONSTANT(ELEM,NODE, n1);
		}
	}
	//法線ベクトルを反転
	if(region_num>1) for(int i=0;i<elemnum;i++) for(int D=0;D<3;D++) ELEM[i].direct[D]*=-1.0;
	cout<<"H,G作成完了"<<endl;

	//実際に解く係数行列作成

	for(int n=0;n<elemnum;n++)	
	{
		//if(n==2036) cout<<ELEM[n].boundary_condition<<endl;
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
				//if(H[m][n2]!=0) cout<<"H "<<m<<" "<<n2<<" "<<H[m][n2]<<endl;
			}
				
			for(int m=0;m<uk;m++)
			{
				matrixC[m][n]=H[m][n];
				if(G[m][n]!=0) cout<<"G "<<m<<" "<<n<<" "<<G[m][n]<<endl;
			}
		}
		//cout<<matrixC[2036][2036]<<endl;
	}
	cout<<"係数行列作成完了"<<endl;
	cout<<matrixC[2036][2036]<<endl;
	cout<<G[2036][2036]<<endl;
	cout<<H[2036][2036]<<endl;

	//これより下でエラー。どこ？？？
	//解行列作成
	for(int n=0;n<uk;n++)
	{
		double val=0;//解行列の値
		for(int m=0;m<uk;m++)
		{
			int elem=Nid[m];			//未知数nに対応する要素番号
			if(bd_type[m]==Diric)
			{
				int node=ELEM[elem].node[0];	//要素がﾃﾞｨﾘｸﾚ型なら、3頂点もﾃﾞｨﾘｸﾚ型で、そのﾃﾞｨﾘｸﾚ値も等しい。なのでnode[0]のﾃﾞｨﾘｸﾚ値を使う
				val+=matrixC2[n][m]*NODE[node].potential;
			}
			else if(bd_type[m]==Neumn)
			{
				int node=ELEM[elem].node[0];	//要素がノイマン型のなら、ノイマン型の節点を探してその微分値を使用
				if(NODE[node].boundary_condition==Neumn) node=ELEM[elem].node[1];
				if(NODE[node].boundary_condition==Neumn) node=ELEM[elem].node[2];

				val+=matrixC2[n][m]*NODE[node].slop1;
			}
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
	double le=CON->get_distancebp();
	for(int n=0;n<elemnum;n++)
	{
		int node=ELEM[n].node[0];
		double Ys=ELEM[n].r[A_Y];
		if(Ys>0 && Ys<5*le)
		{
			fp3<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Z]<<" "<<NODE[node].potential<<endl;
			fp4<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Z]<<" "<<NODE[node].slop1<<endl;
		}
	}
	fp3.close();
	fp4.close();

	
/*
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
	smoothingF3D(CON,PART,fluid_number,F);

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

*/
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

//三角形領域にたいするGauss積分の準備関数
void set_gauss_parameter3D(int *NGS,double ep1[14][3],double *w1,double *ep_ln,double *w_ln)
{
	//ガウス積分の準備 値は[境界要素法-基礎と応用- :J.T.カチカデーリス]のp222より引用
	for(int i=0;i<32;i++) NGS[i]=0;
	NGS[3]=0; NGS[4]=3; NGS[7]=7;	//たとえばGauss積分点が4のときはfor(int i=NGS[4];i<NGS[4]+4;i++) ep1[i][0]=~~~~という使いかたをする

	//gauss_N=3のとき
	ep1[0][0]=0.5;		//評価点1を表す面積座標
	ep1[0][1]=0.5;
	ep1[0][2]=0;
	ep1[1][0]=0;		//評価点2を表す面積座標
	ep1[1][1]=0.5;
	ep1[1][2]=0.5;
	ep1[2][0]=0.5;		//評価点3を表す面積座標
	ep1[2][1]=0;
	ep1[2][2]=0.5;
	w1[0]=1.0/3;		//評価点1での重み
	w1[1]=1.0/3;		//評価点2での重み
	w1[2]=1.0/3;		//評価点3での重み

	//gauss_N=4のとき
	ep1[3][0]=1.0/3;	//評価点1を表す面積座標
	ep1[3][1]=1.0/3;
	ep1[3][2]=1.0/3;
	ep1[4][0]=0.6;		//評価点2を表す面積座標
	ep1[4][1]=0.2;
	ep1[4][2]=0.2;
	ep1[5][0]=0.2;		//評価点3を表す面積座標
	ep1[5][1]=0.6;
	ep1[5][2]=0.2;
	ep1[6][0]=0.2;		//評価点4を表す面積座標
	ep1[6][1]=0.2;
	ep1[6][2]=0.6;
	w1[3]=-0.5625;		//評価点1での重み
	w1[4]=25.0/48.0;	//評価点2での重み
	w1[5]=25.0/48.0;	//評価点3での重み
	w1[6]=25.0/48.0;	//評価点4での重み

	//gauss_N=7のとき
	ep1[7][0]=1.0/3;	//評価点1を表す面積座標
	ep1[7][1]=1.0/3;
	ep1[7][2]=1.0/3;
	ep1[8][0]=0.79742699;		//評価点2を表す面積座標
	ep1[8][1]=0.10128651;
	ep1[8][2]=0.10128651;
	ep1[9][0]=0.10128651;		//評価点3を表す面積座標
	ep1[9][1]=0.79742699;
	ep1[9][2]=0.10128651;
	ep1[10][0]=0.10128651;		//評価点4を表す面積座標
	ep1[10][1]=0.10128651;
	ep1[10][2]=0.79742699;
	ep1[11][0]=0.05971587;		//評価点5を表す面積座標
	ep1[11][1]=0.47014206;
	ep1[11][2]=0.47014206;
	ep1[12][0]=0.47014206;		//評価点6を表す面積座標
	ep1[12][1]=0.05971587;
	ep1[12][2]=0.47014206;
	ep1[13][0]=0.47014206;		//評価点7を表す面積座標
	ep1[13][1]=0.47014206;
	ep1[13][2]=0.05971587;
	w1[7]=0.225;		//評価点1での重み
	w1[8]=0.12593918;	//評価点2での重み
	w1[9]=0.12593918;	//評価点3での重み
	w1[10]=0.12593918;	//評価点4での重み
	w1[11]=0.13239415;	//評価点5での重み
	w1[12]=0.13239415;	//評価点6での重み
	w1[13]=0.13239415;	//評価点7での重み

}

//一定要素のためのHマトリクス計算関数
void set_GHmatrix3D_for_constant(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,double ep[3],double **H,double **G,double w,int type,int n,int m,int n1,int m1)
{
	//n1:ソース点要素番号 n:n1に該当する行列要素番号
	//m1:相手要素番号 m:m1に該当する行列要素番号
	//type=0なら通常 1なら計算しない。
	
	//Xms,Yms:相手要素の中点
	//ep[i]:評価点の面積座標
	//w 評価点における重み

	double S=ELEM[n1].S;		//要素面積

	double Xs[3]={ELEM[n1].r[A_X],ELEM[n1].r[A_Y],ELEM[n1].r[A_Z]};//ソース点座標(要素の中点)

	double Gc[3];				//相手要素のGauss積分評価点座標
	for(int D=0;D<3;D++) Gc[D]=ep[0]*NODE[ELEM[m1].node[0]].r[D]+ep[1]*NODE[ELEM[m1].node[1]].r[D]+ep[2]*NODE[ELEM[m1].node[2]].r[D];

	double RA=0;				//ソース点と評価点の距離
	for(int D=0;D<3;D++) RA+=(Xs[D]-Gc[D])*(Xs[D]-Gc[D]);
	RA=sqrt(RA);

	double RD[3];				//ソース点から評価点へ向かう単位ベクトル
	for(int D=0;D<3;D++) RD[D]=(Gc[D]-Xs[D])/RA;

	double RDN=0;				//上記ベクトルと法線ベクトルとの内積
	for(int D=0;D<3;D++) RDN+=RD[D]*ELEM[n1].direct[D];

	if(type==0)
	{
		H[n][m]+=-1.0/(4*PI*RA*RA)*RDN*w*S;
		G[n][m]+=1.0/(4*PI*RA)*w*S;
	}
}

//一定要素を使用するさいのGの得意積分計算関数
double calc_Gnn_for_CONSTANT(vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &NODE,int n1)
{
	//値は[境界要素解析-理論と応用:田中正隆]のp93を参照した
	//n1:ソース点要素番号 n:n1に該当する行列要素番号
	//m1:相手要素番号 m:m1に該当する行列要素番号
	double val=0;
	int ia=ELEM[n1].node[0];
	int ib=ELEM[n1].node[1];
	int ic=ELEM[n1].node[2];

	double Gp[3];				//三角形の重心
	for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+NODE[ib].r[D]+NODE[ic].r[D])/3;

	double H[3]={0,0,0};//重心と各頂点との距離
	for(int D=0;D<3;D++)
	{
		H[0]+=(Gp[D]-NODE[ia].r[D])*(Gp[D]-NODE[ia].r[D]);//重心とiaとの距離
		H[1]+=(Gp[D]-NODE[ib].r[D])*(Gp[D]-NODE[ib].r[D]);//重心とibとの距離
		H[2]+=(Gp[D]-NODE[ic].r[D])*(Gp[D]-NODE[ic].r[D]);//重心とicとの距離
	}
	for(int j=0;j<3;j++) H[j]=sqrt(H[j]);

	double r12,r23,r31;		//三角形の辺の長さ
	r12=0; r23=0; r31=0;
	for(int D=0;D<3;D++)
	{
		r12+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ib].r[D]-NODE[ia].r[D]);//ia-ibの距離
		r23+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ib].r[D]-NODE[ic].r[D]);//ib-icの距離
		r31+=(NODE[ic].r[D]-NODE[ia].r[D])*(NODE[ic].r[D]-NODE[ia].r[D]);//ic-iaの距離
	}
	r12=sqrt(r12);	r23=sqrt(r23);	r31=sqrt(r31);


	double theta[3];//定義は教科書参照
	//余弦定理よりθを求める
	theta[0]=(H[1]*H[1]+H[2]*H[2]-r23*r23)/(2*H[1]*H[2]);		//この時点ではcosθが格納されていることに注意
	theta[1]=(H[0]*H[0]+H[2]*H[2]-r31*r31)/(2*H[0]*H[2]);
	theta[2]=(H[1]*H[1]+H[0]*H[0]-r12*r12)/(2*H[1]*H[0]);
	for(int j=0;j<3;j++) theta[j]=acos(theta[j]);

	double alpha[3];//定義は教科書参照
	for(int j=0;j<3;j++) alpha[j]=0;
	for(int D=0;D<3;D++)
	{
		alpha[0]+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ic].r[D]-NODE[ia].r[D]);
		alpha[1]+=(NODE[ia].r[D]-NODE[ib].r[D])*(NODE[ic].r[D]-NODE[ib].r[D]);
		alpha[2]+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ia].r[D]-NODE[ic].r[D]);
	}
	alpha[0]/=r12*r31;		//この段階ではcos(2α)が格納されている
	alpha[1]/=r12*r23;
	alpha[2]/=r23*r31;
	for(int j=0;j<3;j++)
	{
		alpha[j]=acos(alpha[j]);//この段階では2αが格納されている
		alpha[j]*=0.5;
	}//alphaが求まった

	double SS=ELEM[n1].S;//面積

	double termA,termB,termC;
	termA=0;termB=0;termC=0;
	termA=log(tan((theta[0]+alpha[1])*0.5)/(tan(alpha[1]*0.5)))/r23;
	termB=log(tan((theta[1]+alpha[2])*0.5)/(tan(alpha[2]*0.5)))/r31;
	termC=log(tan((theta[2]+alpha[0])*0.5)/(tan(alpha[0]*0.5)))/r12;

	val=2*SS/3*(termA+termB+termC);

	return val;

}
