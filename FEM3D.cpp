#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
//#include"define.h"	//#define 格納
//#include"CONFIG.h"	//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
//#include"PART.h"		//class PART定義
#include"BEMclass.h"	//BEM2D関係のclass 定義
//#include"FEM3Dclass.h"
//#include<omp.h>
//#include<vector>
#include<complex>

#include"function.h"

#define FULL 3
#define REMESH 4
#define FULL_INPORT 5

int make_edge_element(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,vector<edge3D> &EDGE,int *branch_num,int **nei2,int KTE,vector<edge3D> &static_EDGE,int t,int node_sta);
void calc_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,int *jnb,int **nei,int *branch_num,double **current);
void denryu_side(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *T,double **current);
void VOLT3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei,double *RP);
void potential_calculation(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei);
void inport_J0_density(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void check_J0(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void check_J0Je(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void Avector3D_node_eddy2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t);
void Avector3D_node_eddy2_jw(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_node_eddy2_jw2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_node_eddy2_jw3(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw);
void Avector3D_node_eddy2_jw5(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw);
void Avector3D_node_eddy2_jw5_ver2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta);
void Avector3D_node_eddy2_jw4(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP);
void Avector3D_edge_eddy(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_edge_eddy_jw(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void Avector3D_edge_eddy_jw2(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,int node_sta,vector<edge3D> &static_EDG,double *Am,double *phi);
void Avector3D_edge_eddy_jw_with_parabolic_node_element(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega);
void calc_transitional_EM_field(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,int *jnb,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int **nei,int particle_node,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta,vector<edge3D> &static_EDGE);
void Bflux3D_node(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **B,int t,int flag);
void Bflux3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *A,double **B,int t,int flag);
///境界条件適用関数(３Ｄ辺要素用)
void set_boundary_condition3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A);
void NODE_F3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpsparticle> &PART,double **F,int fluid_number,int t);
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp);
void calc_eddy_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **old_A,double dt,double *V,double **Je,int t,double *sigma);
void calc_node_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **AR,double **AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t, double TIME,double *sigma, double omega);
void calc_edge_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,double *AR,double *AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t, double TIME,double *sigma, double omega);
int poly3D2(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void laplacian_smoothing(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int node0,int nelm0);
void laplacian_smoothing2(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int *jnb, int **nei,int node0);
void node_sorting(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei);
void node_sorting2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei);
void arrange_matrix_complex(int pn,int *NUM,int **ROW,complex<double>**G);
void check_matrix_symmetry_complex(int pn,int *NUM,int **ROW,complex<double> **G);
void diagonal_distribution(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int pn,int *NUM,int **ROW,double **G);
void COCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void parallel_ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void DS_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void parallel_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);
void cs_ICMRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X);

void modify_node_info(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int**nei,int *newID);
void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void calc_jw_field_node(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm);
void calc_jw_field_edge(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm,int nedge,double *Am, double *phi);



int FEM3D_calculation(mpsconfig *CON,int *static_node,int *static_nelm,int *static_nedge,vector<point3D> &static_NODE, vector<element3D> &static_ELEM,vector<edge3D> &static_EDGE,int particle_node,double **F,int t,double TIME,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,mpsconfig& CONF)
{
	/*//////////////////

	//COCGがまともかどうかを確認するため、簡単な複素数の方程式をを解かせる
	 
	cout<<"test.cocg法スタート"<<endl;

	int pn3=2;
	int number3=4;
	complex<double> Re=(2,5.8);
	complex<double> Im;
	complex<double> z;
	double x=1; double y=2;
	complex<double> *val3 = new complex<double> [number3];
    int *ind3 = new int [number3];//非ゼロ要素の列番号格納配列
    int *ptr3 = new int [pn3+1];//各行の要素がvalの何番目からはじまるのかを格納
	complex<double> *B3=new complex<double> [pn3];//解行列
	complex<double> *XX3=new complex<double> [pn3];//行列の答え格納

	double aaa=0;
	aaa=cos(2*PI*1000000);
	cout<<"aaa="<<aaa<<endl;
	Im=complex<double> (0,1);
	cout<<"Re="<<Re<<endl;
	Im*=2;
	cout<<"Im="<<Im<<endl;
	y*=4;
	z=complex<double> (1,1);
	complex<double> xx;
	xx=complex<double> (2,0);
	z*=xx;
	//z+=complex<double> (0,y);
	//////
	
	cout<<"z="<<z<<endl;

	//////
	//3x+(1+i)y=-i;  (1+i)x+(2-i)y=4;  =>   x=-i; y=1+i; 
	val3[0]=complex<double> (3.0,0.0);
	val3[1]=complex<double> (1.0,1.0);
	val3[2]=complex<double> (1.0,1.0);
	val3[3]=complex<double> (2.0,-1.0);

	ind3[0]=0;
	ind3[1]=1;
	ind3[2]=0;
	ind3[3]=1;

	ptr3[0]=0;
	ptr3[1]=2;
	ptr3[2]=4;

	B3[0]=complex<double> (0.0,-1.0);
	B3[1]=complex<double> (4.0,0.0);
	////////

	//COCG(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);
	//cs_MRTR(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);
	//ICCOCG(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);
	cs_ICMRTR(CON,val3,ind3,ptr3,pn3,B3,number3,XX3);

	for(int i=0;i<pn3;i++) cout<<"xx="<<XX3[i]<<endl;

	delete [] val3;
    delete [] ind3;
    delete [] ptr3;
	delete [] B3;
	delete [] XX3;

	///////*/	

	//jω法で計算した後のステップの間は、周波数特性を利用して必要な物性値を求め、関数を抜ける //
	
	if(CON->get_m_A()==1)
	{
		if(CON->get_jw_Faverage()==OFF)//流体計算の時間刻み幅ごとに、周波数特性をもとに電磁力を計算する//現在廃止
		{   
			/*
			int dh=(int) (1.0/(CON->get_dt()*CON->get_Hz()));//１周期が何分割されているか。dtは1/fをdh等分できるように設定すること
			if(t%(CON->get_jw_interval()*dh)!=1)
			{
				cout<<"周波数特性からAを計算"<<endl;
				calc_jw_field(CON,NODE_jw,ELEM_jw,dt,TIME,PART,fluid_number,F,t,node,nelm);
				return 0;
			}
			*/
		}
		else if(CON->get_jw_Faverage()==ON)//磁束や強制電流の1周期に働く電磁力の平均値を利用する
		{/*
			if(t%(CON->get_jw_interval())!=1)
			{
				//平均値を求める処理はFEMを行ったステップにてすでに終えているので、前のステップの電磁力を読み込む
				cout<<"電磁力の平均値の読み込み"<<endl;
				ifstream f("F_FEM.dat");
				if(!f) cout<<"cannot open F_FEM.dat"<<endl;
				f.unsetf(ifstream::dec);
				f.setf(ifstream::skipws);
				for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) f>>F[D][i];			
				f.close();
				return 0;
			}
		*/
		}
	}

	vector<point3D> NODE;
	vector<element3D> ELEM;
	vector<edge3D> EDGE;

	int node=0;					//全節点数
	int nelm=0;					//全要素数
	int nedge=0;				//全辺数
	int KTJ;					//最大節点数
	int KTE;					//最大要素数　3次元式もとめよ
	double err=1.0e-14;			//誤差判定のしきい値
	

	
	if(CON->get_mesher()==0) //magnetで作ったメッシュをもとにるつぼ内部のメッシュ分割
	{
		int delaun_flag;			//デローニ分割を行うか、行わないか
		int node0=0;

		if(CON->get_mesh_input()==0)							//MPSTOFEMによる節点、要素生成
		{
			if(CON->get_remesh_sw()==OFF) delaun_flag=FULL;		//remesh領域を想定せず、常にすべてをデローニ分割
			else if(CON->get_remesh_sw()==ON)
			{
				if(t==1) delaun_flag=FULL;	//全モデルをデローニ分割
				else delaun_flag=REMESH;	//remesh領域のみデローニ分割
			}
		}
		else if(CON->get_mesh_input()==1)						//Magnetより読み込み
		{
			if(CON->get_remesh_sw()==OFF) delaun_flag=FULL_INPORT;		//常にMagnetの要素を読み込み解析 管理者用？
			else if(CON->get_remesh_sw()==ON)
			{
				if(t==1) delaun_flag=FULL_INPORT;	//全モデルをMagnetファイルより読み込み
				else delaun_flag=REMESH;	//remesh領域のみデローニ分割
			}
		}

		if(delaun_flag==FULL) cout<<"FULL デローニ分割実行"<<endl;
		else if(delaun_flag==REMESH) cout<<"remesh領域のみデローニ分割実行"<<endl;
		else if(delaun_flag==FULL_INPORT) cout<<"Magnet生成ファイルより要素情報等読み込み"<<endl;

		if(delaun_flag==FULL)
		{
			MPS_TO_FEM3Dmain(CON,&node,NODE,PART,  fluid_number,  particle_number);//粒子配置より節点配置を入手
			KTJ=node;	
			if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
			KTE=12*KTJ;
		}
		else if(delaun_flag==REMESH)
		{
			node=(int) static_NODE.size()-1;	//静的節点数
			nelm=(int) static_ELEM.size()-1;
			KTJ=node+particle_number;				//このあと動的節点(流体)を格納しないといけないから、KTJを増加
			if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
			KTE=12*KTJ;
		}
		else if(delaun_flag==FULL_INPORT)
		{
			ifstream fin("input_from_MAGNET.dat");
			if(!fin) cout<<"cannot open input_from_MAGNET.dat"<<endl;
			fin.unsetf(ifstream::dec);
			fin.setf(ifstream::skipws);

			fin>>node;			//節点数読み込み
			fin>>nelm;			//要素数読み込み
			KTJ=node+particle_number;				//このあと動的節点(流体)を格納しないといけないから、KTJを増加
			if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
			KTE=12*KTJ;

			point3D NODE0;
			element3D ELEM0;
			edge3D EDGE0;
			int ID;

			for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);
			for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);	//配列を確保
			for(int i=0;i<KTE;i++) EDGE.push_back(EDGE0);	//配列を確保

			//節点情報読み込み
			for(int i=1;i<=node;i++)
			{
				fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
				for(int D=0;D<3;D++) fin>>NODE[i].r[D];
				fin>>NODE[i].material;
				if(NODE[i].material==21) NODE[i].material=CRUCIBLE;
				NODE[i].boundary_condition=0;			//とりあえずゼロを格納
				NODE[i].particleID=-1;					//対応する粒子は存在しない
				NODE[i].remesh=OFF;						//non-remesh領域の節点である。remesh領域との境界に位置する節点に関しては後に処理を施す
				NODE[i].BD_node=OFF;	
			}
			//要素-節点情報読み込み
			for(int i=1;i<=nelm;i++)
			{
				fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
				for(int j=1;j<=4;j++) fin>>ELEM[i].node[j];
			}
			//要素-要素情報読み込み
			for(int i=1;i<=nelm;i++)
			{
				fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
				for(int j=1;j<=4;j++) fin>>ELEM[i].elm[j];
			}
			//要素材質情報読み込み
			for(int i=1;i<=nelm;i++)
			{
			fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
				fin>>ELEM[i].material;
				//////
				if(ELEM[i].material==21) ELEM[i].material=CRUCIBLE;
				/////
				ELEM[i].map=0;			//初期化
				for(int D=0;D<3;D++)ELEM[i].r[D]=0;
				ELEM[i].RR=0;
				ELEM[i].volume=0;		
			}

			fin.close();

			//境界条件設定
			int *BD_flag=new int[node+1];			//BD_flag=ONなら境界節点(解析境界かもしれないし、remesh境界かもしれない)
			for(int i=0;i<=node;i++) BD_flag[i]=OFF;//初期化
			for(int i=1;i<=nelm;i++)
			{
				for(int j=1;j<=4;j++)
				{
					if(ELEM[i].elm[j]==0)
					{
						int ia=ELEM[i].node[j%4+1];
						int ib=ELEM[i].node[4-(j-1)/2*2];
						int ic=ELEM[i].node[3-(j/2%2)*2];
						BD_flag[ia]=ON;						//境界節点という印
						BD_flag[ib]=ON;
						BD_flag[ic]=ON;
					}
				}
			}
			vector <int> BD_NODE_ID;						//境界節点番号格納
			for(int i=1;i<=node;i++)
			{
				if(BD_flag[i]==ON) BD_NODE_ID.push_back(i);
			}
			int BD_num=(int) BD_NODE_ID.size();
			ofstream fs("remesh.dat");
			if(CON->get_region_shape()==0)				//解析領域が直方体なら
			{
				double Xmax=0;
				double Ymax=0;
				double Zmax=0;

				

				for(int i=0;i<BD_num;i++)
				{
					int n=BD_NODE_ID[i];	//境界節点番号
					if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z面
					else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z面
					else if(NODE[n].r[A_X]<CON->get_XL()+err) NODE[n].boundary_condition=1;		//-X面
					else if(NODE[n].r[A_X]>CON->get_XR()-err) NODE[n].boundary_condition=1;		//+X面
					else if(NODE[n].r[A_Y]<CON->get_YD()+err) NODE[n].boundary_condition=1;		//-Y面
					else if(NODE[n].r[A_Y]>CON->get_YU()-err) NODE[n].boundary_condition=1;		//+Y面
					else
					{
						NODE[n].boundary_condition=0;					//解析境界節点ではなく、remesh領域との境界節点なので、境界条件はゼロ
						NODE[n].remesh=ON;
						NODE[n].BD_node=ON;							//remesh領域の境界部であるというしるし
						fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
					}
					if(NODE[n].r[A_Z]>Zmax) Zmax=NODE[n].r[A_Z];
					if(NODE[n].r[A_Y]>Ymax) Ymax=NODE[n].r[A_Y];
					if(NODE[n].r[A_X]>Xmax) Xmax=NODE[n].r[A_X];
				}
				cout<<Xmax<<" "<<Ymax<<" "<<Zmax<<endl;
			}
			else if(CON->get_region_shape()==1)				//解析領域が円筒なら
			{
				for(int i=0;i<BD_num;i++)
				{
					int n=BD_NODE_ID[i];	//境界節点番号
					double R=sqrt(NODE[n].r[A_X]*NODE[n].r[A_X]+NODE[n].r[A_Y]*NODE[n].r[A_Y]);
					if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z面
					else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z面
					else if(R>CON->get_RU()*0.99) NODE[n].boundary_condition=1;	
					else
					{
						NODE[n].boundary_condition=0;											//解析境界節点ではなく、remesh領域との境界節点なので、境界条件はゼロ
						NODE[n].remesh=ON;
						NODE[n].BD_node=ON;							//remesh領域の境界部であるというしるし
						fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
					}
				}
			}
			delete [] BD_flag;
			fs.close();

			//NODE,ELEMのうち、動かない要素、節点だけをstatic_NODE,staticELEMに格納する
			static_ELEM.clear();
			static_NODE.clear();
			//static_EDGE.clear();
			memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);
			
			if(CON->get_remesh_sw()==ON)
			{
				node=(int) static_NODE.size()-1;
				nelm=(int) static_ELEM.size()-1;
				delaun_flag=REMESH;				//delaun_flagをREMESHにすることで、下のif文に入って動的要素を生成する
				NODE.clear();
				ELEM.clear();
			}
			else cout<<"要素数="<<nelm<<" 節点数＝"<<node<<endl;

		}
		

		if(delaun_flag==FULL)									//すべてをデローニ分割
		{
			point3D NODE0;
			element3D ELEM0;
			for(int i=KTJ+8;i>node;i--) NODE.push_back(NODE0);
			for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);	//配列を確保
			int FINE_sw=CON->get_fine();					//自動再分割スイッチ
			delaun3D_main(CON,NODE,ELEM, KTJ, KTE,&node,&nelm, FINE_sw);

			cout<<"要素数="<<nelm<<" 節点数＝"<<node<<endl;
			///メッシュ生成を確認
			double *val=new double[KTJ+1];
			for(int i=1;i<=node;i++) val[i]=NODE[i].material;
			//data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
			//data_avs2(CON,node,nelm,NODE,ELEM,KTJ,val,t);//断面図
			data_avs3(node,nelm,NODE,ELEM,CON,t);//材質

			if(CON->get_remesh_sw()==ON)				//remesh領域を想定するなら、静的節点・要素情報を記憶しておく
			{
				//NODE,ELEMのうち、動かない要素、節点だけをstatic_NODE,staticELEMに格納する
				static_ELEM.clear();
				static_NODE.clear();
				memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);

				/*/チェック
				int snode=(int) static_NODE.size()-1;
				int snelm=(int) static_ELEM.size()-1;
				for(int i=1;i<=snode;i++) val[i]=static_NODE[i].remesh;
				data_avs2(CON,snode,snelm,static_NODE,static_ELEM,KTJ,val);//断面図*/
			}

			delete [] val;
		}
		else if(delaun_flag==REMESH)									//remesh領域をデローニ分割
		{
			point3D NODE0;
			element3D ELEM0;
			edge3D EDGE0;
			for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);			//節点配列確保 +8はスーパーボックス
			for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);				//配列を確保
			for(int i=0;i<KTE;i++) EDGE.push_back(EDGE0);				//配列を確保

			for(int i=1;i<=node;i++)										//この時点でnodeには静的節点数が格納されている
			{
				for(int D=0;D<3;D++) NODE[i].r[D]=static_NODE[i].r[D];
				NODE[i].boundary_condition=static_NODE[i].boundary_condition;
				NODE[i].material=static_NODE[i].material;
				NODE[i].particleID=static_NODE[i].particleID;
				NODE[i].remesh=static_NODE[i].remesh;
				NODE[i].BD_node=static_NODE[i].BD_node;
			}

			//static_ELEMから情報をcopy
			for(int i=1;i<=nelm;i++)										//この時点でnelmには静的要素数が格納されている
			{
				for(int D=0;D<3;D++) ELEM[i].r[D]=static_ELEM[i].r[D];
				for(int j=1;j<=4;j++)
				{
					ELEM[i].node[j]=static_ELEM[i].node[j];
					ELEM[i].elm[j]=static_ELEM[i].elm[j];
					//if(ELEM[i].elm[j]==0) cout<<i<<endl;
				}
				ELEM[i].map=static_ELEM[i].map;
				ELEM[i].material=static_ELEM[i].material;
				ELEM[i].RR=static_ELEM[i].RR;
				ELEM[i].volume=static_ELEM[i].volume;
				//辺は？？？
			}

			//静的要素はremesh領域に接するところがELEM[i].elm=0となっている.そんな要素を探す
			vector<int> BD_static_ELEM;		//動的要素に接する静的要素番号格納
			vector<int> BD_static_ELEM_elm;	//静的要素が動的要素に接する面番号格納
			for(int i=1;i<=nelm;i++)
			{
				int flag=OFF;
				int J;
				for(int j=1;j<=4;j++) if(NODE[ELEM[i].node[j]].remesh==ON) flag=ON; 
				if(flag==ON)
				{
					
					flag=OFF;
					int num=0;
					for(int j=1;j<=4;j++)
					{
						if(ELEM[i].elm[j]==0)
						{
							flag=ON;
							J=j;
							num++;
						}
					}
					//if(num>1) cout<<i<<" 静的要素が複数面で動的要素と接しています mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].boundary_condition<<NODE[ELEM[i].node[2]].boundary_condition<<NODE[ELEM[i].node[3]].boundary_condition<<NODE[ELEM[i].node[4]].boundary_condition<<endl;
					if(num>1) cout<<i<<" 静的要素が複数面で動的要素と接しています mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].r[A_Z]<<endl;
				}
				if(flag==ON)
				{
					BD_static_ELEM.push_back(i);
					BD_static_ELEM_elm.push_back(J);//要素iは第J面で動的要素に接している
				}
			}
			cout<<"KK="<<BD_static_ELEM.size()<<endl;

			double Pn[3];				//起点の座標

			for(int k=0;k<1;k++)
			{
				int i=BD_static_ELEM[k];
				int j=BD_static_ELEM_elm[k];
				int ia=ELEM[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
				int ib=ELEM[i].node[4-(j-1)/2*2];
				int ic=ELEM[i].node[3-(j/2%2)*2];

				double iaic[3];//ia→icのﾍﾞｸﾄﾙ成分格納
				double iaib[3];//ia→ibのﾍﾞｸﾄﾙ成分格納
				double iaicL=0;	//ia→icの長さ
				for(int D=0;D<3;D++)
				{
					iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
					iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
					iaicL+=iaic[D]*iaic[D];
				}
				iaicL=sqrt(iaicL);
				///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのﾍﾞｸﾄﾙとなる
				double S[3];//上記のﾍﾞｸﾄﾙ成分格納
				S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
				S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
				S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
				double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
				////面積Sがもとまった

				double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//法線ﾍﾞｸﾄﾙ

				double Gp[3];								//表面三角形の重心
				for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

				for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//起点の座標

				//////////////////////
				///起点を強制的に座標指定
				
				Pn[A_X]=0;
				Pn[A_Y]=0;
				Pn[A_Z]=0.22;//0.10325;	//0.225;//	//るつぼ底部：0.10125

				//流体の最大高さを求め、それの少し上に設置

				/*//流体節点の最大高さ、最小高さ
				double Zmax=0;
				double Zmin=10;
				for(int i=1;i<=node;i++)
				{
					if(NODE[i].material==FLUID)
					{
						if(NODE[i].r[A_Z]>Zmax) Zmax=NODE[i].r[A_Z];
						if(NODE[i].r[A_Z]<Zmin) Zmin=NODE[i].r[A_Z];
					}
				}

				Pn[A_X]=0;
				Pn[A_Y]=0;
				Pn[A_Z]=Zmax-CON->get_distancebp();
				*///

				node++;
				for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
				NODE[node].boundary_condition=0;
				NODE[node].material=AIR;
				NODE[node].particleID=-1;
				NODE[node].remesh=ON;

			}
			cout<<"起点解決"<<endl;
			int countBD=0;
			for(int i=1;i<=node;i++) if(NODE[i].BD_node==ON) countBD++;
			cout<<"境界の節点数="<<countBD<<endl;

			int *imen[4];
			for(int D=0;D<4;D++) imen[D]=new int [BD_static_ELEM.size()+1];
			int *jmen =new int [BD_static_ELEM.size()+1];
			int *kmen =new int [BD_static_ELEM.size()+1];
			double *vol =new double [BD_static_ELEM.size()+1];
			int ip=node;

			/*
			if(BD_static_ELEM.size()>100000) cout<<"imenなどのメモリをあげてください"<<endl;
			int ip=node;			
			int imen[100000][3+1];	//多面体表面三角形の節点番号格納
			int jmen[100000];		//多面体表面三角形に隣接する四面体番号格納
			int kmen[100000];		//多面体表面三角形に隣接する四面体の隣接面番号 (相手は第何面で自分と接しているか)
			double vol[100000];		//多面体の体積の６倍
			*/

			for(int k=1;k<=BD_static_ELEM.size();k++)
			{
				int kelm=BD_static_ELEM[k-1];//remesh領域に接する静的要素
				int j=BD_static_ELEM_elm[k-1];
				int ia=ELEM[kelm].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
				int ib=ELEM[kelm].node[4-(j-1)/2*2];
				int ic=ELEM[kelm].node[3-(j/2%2)*2];

				/*
				imen[k][1]=ic;
				imen[k][2]=ib;
				imen[k][3]=ia;
				*/
				imen[1][k]=ic;
				imen[2][k]=ib;
				imen[3][k]=ia;

				jmen[k]=kelm;
				kmen[k]=j;//jelmがielmに接する面番号
				//vol[k]=volume3D(NODE,ia,ib,ic,ip);
				vol[k]=volume3D(NODE,ic,ib,ia,ip);	//節点の順番に注意
				//if(vol[k]<0) cout<<"体積負の要素あり "<<vol[k]<<endl;
			}

			int ibound=(int) BD_static_ELEM.size();//表面の数を表す
			int nelm0=nelm;//変更前の要素数を記憶
			
			for(int i=1;i<=ibound;i++)//要素情報生成
			{   
				nelm++;
				int ielm=nelm0+i;
				int ia=imen[1][i];
				int ib=imen[2][i];
				int ic=imen[3][i];
				ELEM[ielm].node[1]=ia;
				ELEM[ielm].node[2]=ib;
				ELEM[ielm].node[3]=ic;
				ELEM[ielm].node[4]=ip;//新点は４番目と定義	
				ELEM[ielm].elm[4]=jmen[i];
				if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
				ELEM[ielm].volume=vol[i];
				
				sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//外接球の中心（ボロノイ点)と半径の二乗を計算
			
				ELEM[ielm].material=AIR;//通常のpoly関数との相違点。ここで材質を空気と決定する
			}
			///////////////////

			//要素-要素関係修正/////////上の処理で第4面で接する要素番号はわかっているので、残りを求める
			//						ここで、1～3面は多面体を構成する要素との境界面であることに注意
			int ix=0;
			
			for(int i=1;i<=ibound;i++)
			{
				int ielm=nelm0+i;
				for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]はすでにもとまったから、それ以外をもとめる
				{
					///ELEM[ielm].node[4]=ipである
					int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3のとき、2,3,1の順
					int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3のとき、3,1,2の順
					int flag=0;
					for(int k=1;k<=ix;k++)
					{
						if(flag==0)
						{
							int ja=imen[1][k];
							int jb=imen[2][k];
							if(ia==ja && ib==jb)//節点が一致したら
							{
								ELEM[ielm].elm[j]=jmen[k];//あらかじめﾘｽﾄしてあった情報を格納
								ELEM[jmen[k]].elm[kmen[k]]=ielm;
								imen[1][k]=imen[1][ix];		//k番目の情報はもう不要。なので配列の一番最後の情報をk番目にもってきて、それまでの情報は破棄する
								imen[2][k]=imen[2][ix];
								jmen[k]=jmen[ix];
								kmen[k]=kmen[ix];
								ix--;						//待ち辺数減少
								flag=1;						//ELEM[ielm].elm[j]はもとまったので、下のネストに入る必要はないのでflag=1
							}
						}
					}
					if(flag==0)
					{
						ix++;			//ここでのixは、[隣接関係を満たす要素]をまっている[辺]の数を表す。
						imen[1][ix]=ib;	//自分の節点の並びを記憶させ、別の要素がこの並びを満たすのを待つ。ibとiaの並びを逆にしてあることに注意
						imen[2][ix]=ia;
						jmen[ix]=ielm;
						kmen[ix]=j;
					}
				}	
			}///要素-要素関係修正完了

			cout<<"要素作成完了"<<endl;

			fill3D(NODE,ELEM,nelm);

			node0=node;
			cout<<"node0="<<node0<<endl;

			/////////
			for(int D=0;D<4;D++) delete [] imen[D];
			delete [] jmen;
			delete [] kmen;
			delete [] vol;


			/*/境界面から層を形成
			for(int m=1;m<=1;m++)
			{
				for(int k=1;k<BD_static_ELEM.size();k++)
				{
					int i=BD_static_ELEM[k];
					int j=BD_static_ELEM_elm[k];
					int ia=ELEM[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];

					double iaic[3];//ia→icのﾍﾞｸﾄﾙ成分格納
					double iaib[3];//ia→ibのﾍﾞｸﾄﾙ成分格納
					double iaicL=0;	//ia→icの長さ
					for(int D=0;D<3;D++)
					{
						iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
						iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
						iaicL+=iaic[D]*iaic[D];
					}
					iaicL=sqrt(iaicL);
					///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのﾍﾞｸﾄﾙとなる
					double S[3];//上記のﾍﾞｸﾄﾙ成分格納
					S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
					S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
					S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
					double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
					////面積Sがもとまった

					double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//法線ﾍﾞｸﾄﾙ

					double Gp[3];								//表面三角形の重心
					for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

					//for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//起点の座標 オリジナル
					for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*0.0001*m;		//起点の座標

					node++;
					for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
					NODE[node].boundary_condition=0;
					NODE[node].material=AIR;
					NODE[node].particleID=-1;
					NODE[node].remesh=ON;
					NODE[node].BD_node=OFF;
				}
			}
			*/


			double u0=PI*4E-7;			//空気の透磁率
			double skin_depth=sqrt(1.0/(PI*CON->get_Hz()*CON->get_ele_conduc()*u0));//表皮深さ
			cout<<"流体表皮深さ="<<skin_depth<<endl;
			
			ofstream fp("rr.dat");
			/////流体粒子を動的節点として格納
			
			if(CON->get_thinout_fluid()==0)
			{
				/////
				for(int i=0;i<fluid_number;i++)
				{
					//if(PART[i].surface==ON  || i%4==0)
					node++;
					for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
					NODE[node].boundary_condition=0;
					NODE[node].material=FLUID;
					NODE[node].particleID=i;
					NODE[node].remesh=ON;
					NODE[node].BD_node=OFF;
					fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				}
				
			}
			
			if(CON->get_thinout_fluid()>0)
			{
				double Xf=0.0;
				double Yf=0.0;
				double Zf=0.0;
				for(int i=0;i<fluid_number;i++)
				{
					Xf+=PART[i].r[A_X];
					Yf+=PART[i].r[A_Y];
					Zf+=PART[i].r[A_Z];
				}
				Xf/=fluid_number;
				Yf/=fluid_number;
				Zf/=fluid_number;
				cout<<"流体重心=("<<Xf<<","<<Yf<<","<<Zf<<")"<<endl;
				
				for(int i=0;i<fluid_number;i++)
				{
					double rx=PART[i].r[A_X]-Xf;
					double ry=PART[i].r[A_Y]-Yf;
					double rz=PART[i].r[A_Z]-Zf;
					if(PART[i].surface==ON || sqrt(rx*rx+ry*ry+rz*rz)>CON->get_fluidwidth()*0.001-1.1*skin_depth)
					{
						node++;
						for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
						NODE[node].boundary_condition=0;
						NODE[node].material=FLUID;
						NODE[node].particleID=i;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;
						fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					}
					else if(i%CON->get_thinout_fluid()==0)
					{
						node++;
						for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
						NODE[node].boundary_condition=0;
						NODE[node].material=FLUID;
						NODE[node].particleID=i;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;
						fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
					}
				}
			}
			fp.close();
			/////*/



			/*///////
			//表面の精度向上のため、流体節点を追加
			double Rz=0.13125-0.003;
			//まずは半球を作る。そのためには半球表面を作成する必要がある。
			for(int Ri=1;Ri<=5;Ri++)
			{
				double R=0.025-5*0.0001+Ri*0.0002;
				double le=CON->get_distancebp();

				double A=sqrt(3.0)/2;				//よく使う係数
				double B=sqrt(2.0/3);						////よく使う係数
				int half_WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
				int half_WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
				int half_WZ=(int)(R/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
				double R2=R-0.5*le;				//少し小さめの半径を設定

				///////////半球表面
				int Nt;						//球表面の、θ方向の分割数
				double Lt;					//球表面の、θ方向の分割距離
				//calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//半円の分割数は偶数・奇数どちらでもよい
				double temp_N=PI/2*R/(le*A);			//仮の分割数。leで割り切れたら一番いいけど、そうもいかないときがある
				int Ns=(int) temp_N;				//真の分割数
				double difference=temp_N-Ns;		//仮と真の差
				if(difference>0.5) Ns++;
				Lt=PI/2*R/Ns;			//粒子の距離
				Nt=Ns;

				double d_theta=Lt/R;		//弧の長さがLtになる角度

				for(int k=0;k<=Nt;k++)//loopはk<Ntで終わらせる。Ntに該当するところはすでに設置済み
				{
					double THETA=k*d_theta;	//θ
					double r=R*sin(THETA);	//その高さにおける円の半径
					double round=2*PI*r;//その高さにおける円周

					//int Nf=calc_division_N_circle(round,le);//球表面の、θ方向の分割数
					//対称性を考えた時、円周を記述する粒子は偶数個でなければならない。だから他の辺分割数とは扱いが少し特殊
					//dis:分割する距離(円周)
					double temp_num=round/le;		//円外周に設置する『仮の』粒子数。ただし外周がうまくleで割り切れるとは限らない

					int N1=(int)(temp_num/2);
					N1*=2;							
					int N2=N1+2;					//temp_numはN1とN2の間にある。ここでN1,N2は偶数

					double dif1=temp_num-N1;		//各Nとの差
					double dif2=N2-temp_num;
					int N=N1;						//周方向分割数
					if(dif2<dif1) N=N2;				//差の小さい方をNとして採用する。

					int Nf=N;
					
					
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
						double xf=r*cos(fai);
						double yf=r*sin(fai);
						double zf=R*cos(THETA);
						
						node++;
						NODE[node].r[A_X]=xf;
						NODE[node].r[A_Y]=yf;
						NODE[node].r[A_Z]=zf+Rz;
						NODE[node].boundary_condition=0;
						//NODE[node].material=FLUID;
						if(Ri<=3) NODE[node].material=FLUID;
						if(Ri>3) NODE[node].material=AIR;
						NODE[node].particleID=0;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;

						if(k!=Nt)
						{
							node++;
							NODE[node].r[A_X]=xf;
							NODE[node].r[A_Y]=yf;
							NODE[node].r[A_Z]=-zf+Rz;
							NODE[node].boundary_condition=0;
							//NODE[node].material=FLUID;
							if(Ri<=3) NODE[node].material=FLUID;
							if(Ri>3) NODE[node].material=AIR;
							NODE[node].particleID=0;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;
						}
					}
				}
				if(Nt%2!=0)//Ntが奇数のときは、頂上に粒子が置かれなければならない。しかし上のloopはそれが不可能。よってここで追加
				{
					node++;
					NODE[node].r[A_X]=0;
					NODE[node].r[A_Y]=0;
					NODE[node].r[A_Z]=R+Rz;
					NODE[node].boundary_condition=0;
					//NODE[node].material=FLUID;
					if(Ri<=3) NODE[node].material=FLUID;
					if(Ri>3) NODE[node].material=AIR;
					NODE[node].particleID=0;
					NODE[node].remesh=ON;
					NODE[node].BD_node=OFF;

					
						node++;
						NODE[node].r[A_X]=0;
						NODE[node].r[A_Y]=0;
						NODE[node].r[A_Z]=-R+Rz;
						NODE[node].boundary_condition=0;
						//NODE[node].material=FLUID;
						if(Ri<=3) NODE[node].material=FLUID;
						if(Ri>3) NODE[node].material=AIR;
						NODE[node].particleID=0;
						NODE[node].remesh=ON;
						NODE[node].BD_node=OFF;
				
				}
				//////////////////////////////////////
			}

			//*////
			unsigned int timeD=GetTickCount();
			
			/*
			int checkmax=0;
			int checkmin=node;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].BD_node==ON)
				{
					if(i>checkmax) checkmax=i;
					if(i<checkmin) checkmin=i;
				}
			}
			cout<<"境界節点番号の最大値、最小値＝"<<checkmax<<","<<checkmin<<endl;

			int checkbig=0;
			int checksmall=0;
			double check=node/2.0;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].BD_node==ON)
				{
					if(i>check) checkbig++;
					if(i<check) checksmall++;
				}
			}
			cout<<"境界節点番号：総数の半分と比較した大小＝"<<checkbig<<","<<checksmall<<endl;

			*/
			
			cout<<"node="<<node<<" デローニ分割"<<endl;

			for(int i=1;i<=nelm;i++)
			{
				ELEM[i].remesh=OFF; //初期化
				ELEM[i].volume=volume3D(NODE,ELEM[i].node[1],ELEM[i].node[2],ELEM[i].node[3],ELEM[i].node[4]);
				if(ELEM[i].volume<0) cout<<"体積負の要素番号="<<i<<" "<<"体積="<<ELEM[i].volume<<endl;

				int R1=NODE[ELEM[i].node[1]].remesh;
				int R2=NODE[ELEM[i].node[2]].remesh;
				int R3=NODE[ELEM[i].node[3]].remesh;
				int R4=NODE[ELEM[i].node[4]].remesh;
				
				if(R1==ON && R2==ON && R3==ON && R4==ON) ELEM[i].remesh=ON;//NODE.remeshは、ほかの関数の都合上リメッシュ境界でONになっている

				sphere3D(NODE,ELEM,ELEM[i].node[1],ELEM[i].node[2],ELEM[i].node[3],ELEM[i].node[4],i);
			}

			for(int i=1;i<KTE;i++) ELEM[i].map=0;//初期化
			err=1e-14;
			///順次節点を導入していく
			int *kv=new int[KTE];//新節点を外接球にふくむ要素群
			int *istack=new int[KTE];//一時配列

			//////////ふつうの順番に流体を追加
			if(CON->get_defer_f()==OFF)
			{
				for(int i=node0+1;i<=node;i++)
				{   
					//cout<<i<<endl;
					int ip=i;
					double xp=NODE[ip].r[A_X];//導入する節点の座標
					double yp=NODE[ip].r[A_Y];
					double zp=NODE[ip].r[A_Z];
			
					///新節点を含む要素の探索
					int loc=locate3D(NODE,ELEM,nelm,xp,yp,zp);
				
					//////////外接球内に新節点を含む要素の抽出
					int iv=0;
					int msk=0;
			
					iv++;//外接球内に新節点を含む要素数
					kv[iv]=loc;
					ELEM[loc].map=1;//mapが1の要素は、外接球に節点iを含むかどうかを検査済みということ
					msk++;
					istack[msk]=loc;
					if(loc==0)cout<<"loc==0"<<endl;
					
					if(CON->get_CDT_sw()==OFF)//通常のデローニ
					{
						while(msk!=0)
						{   
							int isk=istack[msk];//いま注目している要素の番号
							msk--;
							for(int j=1;j<=4;j++)
							{
								int jelm=ELEM[isk].elm[j];//iskと接する要素
								if(jelm!=0)//それが表面でないなら
								{
									if(ELEM[jelm].map!=1) //まだ検査してないなら
									{   
					         
										double rad=ELEM[jelm].RR*(1.000000+err);//外接球半径の２乗
										double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
									
										if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
										{
											dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
											if(dst<rad)
											{
												dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
												if(dst<rad)//外接球内に含む
												{
													iv++;//外接球内に新節点を含む要素数を＋１
													kv[iv]=jelm;//リストにいれる
													ELEM[jelm].map=1;//外接球に新点含むというしるし
													msk++;
													istack[msk]=jelm;
												}
											}
										}
									}
								}
							}
						}
					}//外接球内に新節点を含む要素数ivと、その要素番号kv[iv]がもとまった

					if(CON->get_CDT_sw()==ON)//制約付きデローニ。remesh境界を壊さないようにしたい
					{
						while(msk!=0)
						{   
							int isk=istack[msk];//いま注目している要素の番号
							msk--;
							for(int j=1;j<=4;j++)
							{
								int jelm=ELEM[isk].elm[j];//iskと接する要素
								if(jelm!=0)//それが表面でないなら
								{
									if(ELEM[jelm].map!=1) //まだ検査してないなら
									{   
										double rad=ELEM[jelm].RR*(1.000000+err);//外接球半径の２乗
										double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
								
										if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
										{
											dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
											if(dst<rad)
											{
												dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
												if(dst<rad)//外接球内に含む
												{
													if(ELEM[jelm].remesh==ON)//iskと接する要素を構成する節点がすべてremesh節点、すなわちremesh領域内の要素なら
													{
														iv++;//外接球内に新節点を含む要素数を＋１
														kv[iv]=jelm;//リストにいれる
														ELEM[jelm].map=1;//外接球に新点含むというしるし
														msk++;
														istack[msk]=jelm;
													}
												}
											}
										}
									}
								}
							}
						}
					}//外接球内に新節点を含む要素数ivと、その要素番号kv[iv]がもとまった
					
					/////////////////
					
					////得られた多面体を四面体に分割する
					poly3D2(NODE,ELEM,&iv,kv,ip,&nelm,CON); 
				}
			}

			/////////多面体形成がうまくいかなかった流体節点の追加を後回しにする
			if(CON->get_defer_f()==ON)
			{
				int *flag=new int [node+1];
				for(int i=1;i<=node;i++) flag[i]=ON;
				int *jnb3=new int[node+1];///各節点に隣接する要素数格納 流体節点を再投入するかの条件(jnb=0)
				int num_ON=node;//追加できていない節点の数
				int count=0;
				for(int i=1;i<=node0;i++)//リメッシュ節点（起点を除く）のみが新しく追加する節点
				{
					flag[i]=OFF;
					num_ON--;
				}
				
				while(num_ON!=0)
				{
					if(count>0) cout<<"要素形成に失敗した流体節点の再配置　num="<<num_ON<<endl;
					//if(count>0) cout<<"流体節点の再配置開始 count="<<count+1<<endl;
					count++;
					for(int i=node0+1;i<=node;i++)
					{
						//cout<<i<<endl;
						if(flag[i]==ON)
						{
							int ip=i;
							double xp=NODE[ip].r[A_X];//導入する節点の座標
							double yp=NODE[ip].r[A_Y];
							double zp=NODE[ip].r[A_Z];
						
							///新節点を含む要素の探索
							//unsigned int timeA=GetTickCount();
							int loc=locate3D(NODE,ELEM,nelm,xp,yp,zp);
							//cout<<"locate完了－－time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
						
							//////////外接球内に新節点を含む要素の抽出
							int iv=0;
							int msk=0;
					
							iv++;//外接球内に新節点を含む要素数
							kv[iv]=loc;
							ELEM[loc].map=1;//mapが1の要素は、外接球に節点iを含むかどうかを検査済みということ
							msk++;
							istack[msk]=loc;
							if(loc==0)cout<<"loc==0"<<endl;
							
							if(CON->get_CDT_sw()==OFF)//通常のデローニ
							{
								while(msk!=0)
								{   
									int isk=istack[msk];//いま注目している要素の番号
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//iskと接する要素
										if(jelm!=0)//それが表面でないなら
										{
											if(ELEM[jelm].map!=1) //まだ検査してないなら
											{   
							         
												//double rad=ELEM[jelm].RR*(1.000000+err);//外接球半径の２乗
												double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//外接球半径の２乗
												double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
											
												if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
												{
													dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
													if(dst<rad)
													{
														dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
														if(dst<rad)//外接球内に含む
														{
															iv++;//外接球内に新節点を含む要素数を＋１
															kv[iv]=jelm;//リストにいれる
															ELEM[jelm].map=1;//外接球に新点含むというしるし
															msk++;
															istack[msk]=jelm;
														}
													}
												}
											}
										}
									}
								}
							}//外接球内に新節点を含む要素数ivと、その要素番号kv[iv]がもとまった

							if(CON->get_CDT_sw()==ON)//制約付きデローニ。remesh境界を壊さないようにしたい
							{
								while(msk!=0)
								{   
									int isk=istack[msk];//いま注目している要素の番号
									msk--;
									for(int j=1;j<=4;j++)
									{
										int jelm=ELEM[isk].elm[j];//iskと接する要素
										if(ELEM[jelm].remesh==ON)//iskと接する要素を構成する節点がすべてremesh節点、すなわちremesh領域内の要素なら
										{
											if(jelm!=0)//それが表面でないなら
											{
												if(ELEM[jelm].map!=1) //まだ検査してないなら
												{   
													//double rad=ELEM[jelm].RR*(1.000000+err);//外接球半径の２乗
													double rad=ELEM[jelm].RR*(1.000000+1e-12)*(1.000000+1e-12);//外接球半径の２乗
													double dst=ELEM[jelm].r[A_X]*ELEM[jelm].r[A_X]-2.000000*ELEM[jelm].r[A_X]*xp+xp*xp;///外接球中心と新節点の距離
											
													if(dst<rad)///ここでdst>radなら絶対に外接球に含まない
													{
														dst+=ELEM[jelm].r[A_Y]*ELEM[jelm].r[A_Y]-2.000000*ELEM[jelm].r[A_Y]*yp+yp*yp;
														if(dst<rad)
														{
															dst+=ELEM[jelm].r[A_Z]*ELEM[jelm].r[A_Z]-2.000000*ELEM[jelm].r[A_Z]*zp+zp*zp;
															if(dst<rad)//外接球内に含む
															{
																iv++;//外接球内に新節点を含む要素数を＋１
																kv[iv]=jelm;//リストにいれる
																ELEM[jelm].map=1;//外接球に新点含むというしるし
																msk++;
																istack[msk]=jelm;
															}
														}
													}
												}
											}
										}
									}
								}
							}//外接球内に新節点を含む要素数ivと、その要素番号kv[iv]がもとまった


							
							/////////////////
							
							////得られた多面体を四面体に分割する
							//unsigned int timeB=GetTickCount();
							poly3D2(NODE,ELEM,&iv,kv,ip,&nelm,CON);
							//cout<<"poly完了－－time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;

							
						}
					}
				
					set_jnb3D(NODE,ELEM,node,nelm,jnb3);
					for(int i=node0+1;i<=node;i++)
					{
						if(flag[i]==ON) num_ON--;
						flag[i]=OFF;
						
						if(jnb3[i]==0)
						{
							flag[i]=ON;
							num_ON++;
						}
						/*
						if(flag[i]==ON)
						{
							if(jnb3[i]>0)
							{
								num_ON--;
								flag[i]=OFF;
							}
						}
						*/
					}
				}
				delete [] jnb3;
				delete [] flag;
			}
			
			delete [] kv;
			delete [] istack;//*/

			cout<<"デローニ分割完了－－time="<<(GetTickCount()-timeD)*0.001<<"[sec]"<<endl;

			/////////

			////材質判定
			//set_material(CON,NODE,ELEM,node,nelm);
			
			int countmate=0;
			////材質判定
			for(int i=1;i<=nelm;i++)
			{
				if(ELEM[i].material==AIR)
				{
					
					int M1=NODE[ELEM[i].node[1]].material;
					int M2=NODE[ELEM[i].node[2]].material;
					int M3=NODE[ELEM[i].node[3]].material;
					int M4=NODE[ELEM[i].node[4]].material;

					int B1=NODE[ELEM[i].node[1]].BD_node;
					int B2=NODE[ELEM[i].node[2]].BD_node;
					int B3=NODE[ELEM[i].node[3]].BD_node;
					int B4=NODE[ELEM[i].node[4]].BD_node;

					int F1=0;
					int F2=0;
					int F3=0;
					int F4=0;

					if(M1==FLUID || B1==ON) F1=ON;
					if(M2==FLUID || B2==ON) F2=ON;
					if(M3==FLUID || B3==ON) F3=ON;
					if(M4==FLUID || B4==ON) F4=ON;
			
					///4頂点すべてが同じ材質なら要素もそれにならう。
					///ひとつでも異なっていたら空気と定義
					if(M1==M2 && M2==M3 && M3==M4)
					{
						if(M1==FLUID) ELEM[i].material=FLUID;
					}

					////リメッシュ境界節点と流体節点のみで構成される要素を流体とする
					if(CON->get_BD_fluid()==ON)
					{
						if(F1==F2 && F2==F3 && F3==F4)
						{
							if(F1==ON) ELEM[i].material=FLUID;
						}
					}
				}
			}
					//*//////
			for(int i=1;i<=nelm;i++)
			{
				//最大辺長さが一定値を超えた場合、流体要素ではないと判断し空気要素に変える
				if( ELEM[i].material==FLUID)
				{
					if(CON->get_material_judge()>0)
					{
						int ia=ELEM[i].node[1];
						int ib=ELEM[i].node[2];
						int ic=ELEM[i].node[3];
						int id=ELEM[i].node[4];
							
						double iaib=0;//点ia,ibの距離
						double iaic=0;//点ia,icの距離
						double iaid=0;//点ia,idの距離
						double ibic=0;//点ib,icの距離
						double ibid=0;//点ib,idの距離
						double icid=0;//点ic,idの距離

						for(int D=0;D<3;D++)
						{
							iaib+=(NODE[ib].r[D]-NODE[ia].r[D])*(NODE[ib].r[D]-NODE[ia].r[D]);
							iaic+=(NODE[ia].r[D]-NODE[ic].r[D])*(NODE[ia].r[D]-NODE[ic].r[D]);
							iaid+=(NODE[ia].r[D]-NODE[id].r[D])*(NODE[ia].r[D]-NODE[id].r[D]);
							ibic+=(NODE[ib].r[D]-NODE[ic].r[D])*(NODE[ib].r[D]-NODE[ic].r[D]);
							ibid+=(NODE[ib].r[D]-NODE[id].r[D])*(NODE[ib].r[D]-NODE[id].r[D]);
							icid+=(NODE[ic].r[D]-NODE[id].r[D])*(NODE[ic].r[D]-NODE[id].r[D]);
						}
						iaib=sqrt(iaib);
						iaic=sqrt(iaic);
						iaid=sqrt(iaid);
						ibic=sqrt(ibic);
						ibid=sqrt(ibid);
						icid=sqrt(icid);

						double maxL=iaib;//最大辺長さ
						if(iaic>maxL) maxL=iaic;
						if(iaid>maxL) maxL=iaid;
						if(ibic>maxL) maxL=ibic;
						if(ibid>maxL) maxL=ibid;
						if(icid>maxL) maxL=icid;

						if(maxL>CON->get_material_judge()*CON->get_distancebp())
						{
							countmate++;
							ELEM[i].material=AIR;
						}
					}
				}
			}
			if(countmate>0) cout<<"最大辺長さが"<<CON->get_material_judge()<<"leを超え、空気に変換された流体要素の数="<<countmate<<endl;


			/////メッシュの細分化
			FINE3D(NODE,ELEM,KTJ,KTE,&node,&nelm,CON,1,nelm0);
			//*/
		}
		

		////要素がうまく生成されているかチェックする
		fill3D(NODE,ELEM,nelm);

		cout<<"要素数="<<nelm<<" 節点数＝"<<node<<endl;
		/*//メッシュ生成を確認
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
		//for(int i=1;i<=dnode;i++) val[i]=dy_NODE[i].material;
		//data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
		//data_avs2(CON,dnode,dnelm,dy_NODE,dy_ELEM,KTJ,val);//断面図
		
		if(TIME==0)
		{
			data_avs2(CON,node,nelm,NODE,ELEM,KTJ,val,t);//断面図
			data_avs3(node,nelm,NODE,ELEM,CON,t);//材質
		}
		delete [] val;
		/////*/
		
	}
	else if(CON->get_mesher()==1)	//TetGenによるメッシュ生成
	{
		if(CON->get_dimention()==3)
		{
			//delaun3D関係の変数は全て不要

			vector<int> TRANS;		//iは節点番号  節点iはTRANS[i]番目の粒子に相当
			TRANS.push_back(-1);	//最初のデータは使わないので詰めておく

			//TetGen用配列作成
			vector<tetgen_node> NODEall;
			vector<tetgen_facet> FACEall;
			vector<tetgen_element> ELEMall;

			///TetGenによるメッシュ生成////////////////////////////////////////////////////

			tetgen_function TETFUNC;

			if(CON->get_mesh_input()==0)//すべて自分で用意
			{
				TETFUNC.call_TetGen(CONF,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
			}

			//外部からの読込は現在対応していない
			if(CON->get_mesh_input()==1)//magnetのメッシュを読込
			{
				int delaun_flag;			//デローニ分割を行うか、行わないか
				int node0=0;

				if(CON->get_mesh_input()==0)							//MPSTOFEMによる節点、要素生成
				{
					if(CON->get_remesh_sw()==OFF) delaun_flag=FULL;		//remesh領域を想定せず、常にすべてをデローニ分割
					else if(CON->get_remesh_sw()==ON)
					{
						if(t==1) delaun_flag=FULL;	//全モデルをデローニ分割
						else delaun_flag=REMESH;	//remesh領域のみデローニ分割
					}
				}
				else if(CON->get_mesh_input()==1)						//Magnetより読み込み
				{
					if(CON->get_remesh_sw()==OFF) delaun_flag=FULL_INPORT;		//常にMagnetの要素を読み込み解析 管理者用？
					else if(CON->get_remesh_sw()==ON)
					{
						if(t==1) delaun_flag=FULL_INPORT;	//全モデルをMagnetファイルより読み込み
						else delaun_flag=REMESH;	//remesh領域のみデローニ分割
					}
				}

				if(delaun_flag==FULL) cout<<"FULL デローニ分割実行"<<endl;
				else if(delaun_flag==REMESH) cout<<"remesh領域のみデローニ分割実行(tetgen)"<<endl;
				else if(delaun_flag==FULL_INPORT) cout<<"Magnet生成ファイルより要素情報等読み込み(tetgen)"<<endl;

				if(delaun_flag==FULL)
				{
					MPS_TO_FEM3Dmain(CON,&node,NODE,PART,  fluid_number,  particle_number);//粒子配置より節点配置を入手
					KTJ=node;	
					if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
					KTE=12*KTJ;
				}
				else if(delaun_flag==REMESH)
				{
					node=(int) static_NODE.size()-1;	//静的節点数
					nelm=(int) static_ELEM.size()-1;
					KTJ=node+particle_number;				//このあと動的節点(流体)を格納しないといけないから、KTJを増加
					if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
					KTE=12*KTJ;
				}
				else if(delaun_flag==FULL_INPORT)
				{
					ifstream fin("input_from_MAGNET.dat");
					if(!fin) cout<<"cannot open input_from_MAGNET.dat"<<endl;
					fin.unsetf(ifstream::dec);
					fin.setf(ifstream::skipws);

					fin>>node;			//節点数読み込み
					fin>>nelm;			//要素数読み込み
					KTJ=node+particle_number;				//このあと動的節点(流体)を格納しないといけないから、KTJを増加
					if(CON->get_fine()!=OFF) KTJ+=CON->get_add_points();
					KTE=12*KTJ;

					point3D NODE0;
					element3D ELEM0;
					edge3D EDGE0;
					int ID;

					for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);
					for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);	//配列を確保
					for(int i=0;i<KTE;i++) EDGE.push_back(EDGE0);	//配列を確保
					
					//節点情報読込
					for(int i=1;i<=node;i++)
					{
						fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
						for(int D=0;D<3;D++) fin>>NODE[i].r[D];
						fin>>NODE[i].material;
						if(NODE[i].material==21) NODE[i].material=CRUCIBLE;
						NODE[i].boundary_condition=0;			//とりあえずゼロを格納
						NODE[i].particleID=-1;					//対応する粒子は存在しない
						NODE[i].remesh=OFF;						//non-remesh領域の節点である。remesh領域との境界に位置する節点に関しては後に処理を施す
						NODE[i].BD_node=OFF;		
					}

					/*//////節点情報をstatic.nodeに出力
					ofstream fout("static.node");
	
					fout<<"#node"<<endl;
					fout<<node<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=node;i++)
					{
						fout<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].material<<" "<<NODE[i].material<<endl;
					}

					fout.close();
					//////*/

					
					//要素-節点情報読み込み
					for(int i=1;i<=nelm;i++)
					{
						fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
						for(int j=1;j<=4;j++) fin>>ELEM[i].node[j];
					}
					//要素-要素情報読み込み
					for(int i=1;i<=nelm;i++)
					{
						fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
						for(int j=1;j<=4;j++) fin>>ELEM[i].elm[j];
					}
					//要素材質情報読み込み
					for(int i=1;i<=nelm;i++)
					{
						fin>>ID;								//どうでもいい情報なので何にも使用しない。ただしファイルの左端にはこのように節点番号をつけておくように。そのほうがチェックが楽
						fin>>ELEM[i].material;
						//////
						if(ELEM[i].material==21) ELEM[i].material=CRUCIBLE;
						/////
						ELEM[i].map=0;			//初期化
						for(int D=0;D<3;D++)ELEM[i].r[D]=0;
						ELEM[i].RR=0;
						ELEM[i].volume=0;		
					}
					/*//////要素情報をstatic.eleに出力
					ofstream fout2("static.ele");
	
					fout2<<nelm<<" "<<"4"<<" "<<"1"<<endl;
					for(int i=1;i<=nelm;i++)
					{
						fout2<<i<<" "<<ELEM[i].node[1]<<" "<<ELEM[i].node[2]<<" "<<ELEM[i].node[3]<<" "<<ELEM[i].node[4]<<" "<<ELEM[i].material<<endl;
					}

					fout2.close();

					fin.close();
					/////////////*/
					/*///////////流体節点をfluid.nodeに出力
					ofstream fp("rr.dat");
					ofstream fl("fluid.node");
					/////
					fl<<"#node"<<endl;
					fl<<fluid_number<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=fluid_number;i++)
					{
						fl<<i<<" "<<PART[i-1].r[A_X]<<" "<<PART[i-1].r[A_Y]<<" "<<PART[i-1].r[A_Z]<<" "<<WATER<<" "<<WATER<<endl;
					}
					

					fp.close();
					fl.close();
					///////*/

					//////static.node,eleからメッシュ再構成
					//TETFUNC.call_TetGen(CONF,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
					////////

					//境界条件設定
					int *BD_flag=new int[node+1];			//BD_flag=ONなら境界節点(解析境界かもしれないし、remesh境界かもしれない)
					for(int i=0;i<=node;i++) BD_flag[i]=OFF;//初期化
					for(int i=1;i<=nelm;i++)
					{
						for(int j=1;j<=4;j++)
						{
							if(ELEM[i].elm[j]==0)
							{
								int ia=ELEM[i].node[j%4+1];
								int ib=ELEM[i].node[4-(j-1)/2*2];
								int ic=ELEM[i].node[3-(j/2%2)*2];
								BD_flag[ia]=ON;						//境界節点という印
								BD_flag[ib]=ON;
								BD_flag[ic]=ON;
							}
						}
					}
					vector <int> BD_NODE_ID;						//境界節点番号格納
					for(int i=1;i<=node;i++)
					{
						if(BD_flag[i]==ON) BD_NODE_ID.push_back(i);
					}
					int BD_num=(int) BD_NODE_ID.size();
					ofstream fs("remesh.dat");
					if(CON->get_region_shape()==0)				//解析領域が直方体なら
					{
						double Xmax=0;
						double Ymax=0;
						double Zmax=0;

				

						for(int i=0;i<BD_num;i++)
						{
							int n=BD_NODE_ID[i];	//境界節点番号
							if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z面
							else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z面
							else if(NODE[n].r[A_X]<CON->get_XL()+err) NODE[n].boundary_condition=1;		//-X面
							else if(NODE[n].r[A_X]>CON->get_XR()-err) NODE[n].boundary_condition=1;		//+X面
							else if(NODE[n].r[A_Y]<CON->get_YD()+err) NODE[n].boundary_condition=1;		//-Y面
							else if(NODE[n].r[A_Y]>CON->get_YU()-err) NODE[n].boundary_condition=1;		//+Y面
							else
							{
								NODE[n].boundary_condition=0;					//解析境界節点ではなく、remesh領域との境界節点なので、境界条件はゼロ
								NODE[n].remesh=ON;
								NODE[n].BD_node=ON;							//remesh領域の境界部であるというしるし
								fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
							}
							if(NODE[n].r[A_Z]>Zmax) Zmax=NODE[n].r[A_Z];
							if(NODE[n].r[A_Y]>Ymax) Ymax=NODE[n].r[A_Y];
							if(NODE[n].r[A_X]>Xmax) Xmax=NODE[n].r[A_X];
						}
						cout<<Xmax<<" "<<Ymax<<" "<<Zmax<<endl;
					}
					else if(CON->get_region_shape()==1)				//解析領域が円筒なら
					{
						for(int i=0;i<BD_num;i++)
						{
							int n=BD_NODE_ID[i];	//境界節点番号
							double R=sqrt(NODE[n].r[A_X]*NODE[n].r[A_X]+NODE[n].r[A_Y]*NODE[n].r[A_Y]);
							if(NODE[n].r[A_Z]<CON->get_ZD()+err) NODE[n].boundary_condition=1;		//-Z面
							else if(NODE[n].r[A_Z]>CON->get_ZU()-err) NODE[n].boundary_condition=1;		//+Z面
							else if(R>CON->get_RU()*0.99) NODE[n].boundary_condition=1;	
							else
							{
								NODE[n].boundary_condition=0;											//解析境界節点ではなく、remesh領域との境界節点なので、境界条件はゼロ
								NODE[n].remesh=ON;
								NODE[n].BD_node=ON;							//remesh領域の境界部であるというしるし
								fs<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<NODE[n].r[A_Z]<<endl;
							}
						}
					}
					delete [] BD_flag;
					fs.close();

					//NODE,ELEMのうち、動かない要素、節点だけをstatic_NODE,staticELEMに格納する
					static_ELEM.clear();
					static_NODE.clear();
					memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);
			
					if(CON->get_remesh_sw()==ON)
					{
						node=(int) static_NODE.size()-1;
						nelm=(int) static_ELEM.size()-1;
						delaun_flag=REMESH;				//delaun_flagをREMESHにすることで、下のif文に入って動的要素を生成する
						NODE.clear();
						ELEM.clear();
					}
					else cout<<"要素数="<<nelm<<" 節点数＝"<<node<<endl;

				}
		
				if(delaun_flag==REMESH)									//remesh領域をデローニ分割
				{
					point3D NODE0;
					element3D ELEM0;
					for(int i=0;i<=KTJ+8;i++) NODE.push_back(NODE0);			//節点配列確保 +8はスーパーボックス
					for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);				//配列を確保

					for(int i=1;i<=node;i++)										//この時点でnodeには静的節点数が格納されている
					{
						for(int D=0;D<3;D++) NODE[i].r[D]=static_NODE[i].r[D];
						NODE[i].boundary_condition=static_NODE[i].boundary_condition;
						NODE[i].material=static_NODE[i].material;
						NODE[i].particleID=static_NODE[i].particleID;
						NODE[i].remesh=static_NODE[i].remesh;
						NODE[i].BD_node=static_NODE[i].BD_node;
					}

					//static_ELEMから情報をcopy
					for(int i=1;i<=nelm;i++)										//この時点でnelmには静的要素数が格納されている
					{
						for(int D=0;D<3;D++) ELEM[i].r[D]=static_ELEM[i].r[D];
						for(int j=1;j<=4;j++)
						{
							ELEM[i].node[j]=static_ELEM[i].node[j];
							ELEM[i].elm[j]=static_ELEM[i].elm[j];
							//if(ELEM[i].elm[j]==0) cout<<i<<endl;
						}
						ELEM[i].map=static_ELEM[i].map;
						ELEM[i].material=static_ELEM[i].material;
						ELEM[i].RR=static_ELEM[i].RR;
						ELEM[i].volume=static_ELEM[i].volume;
						//辺は？？？
					}

					//静的要素はremesh領域に接するところがELEM[i].elm=0となっている.そんな要素を探す
					vector<int> BD_static_ELEM;		//動的要素に接する静的要素番号格納
					vector<int> BD_static_ELEM_elm;	//静的要素が動的要素に接する面番号格納
					for(int i=1;i<=nelm;i++)
					{
						int flag=OFF;
						int J;
						for(int j=1;j<=4;j++) if(NODE[ELEM[i].node[j]].remesh==ON) flag=ON; 
						if(flag==ON)
						{
					
							flag=OFF;
							int num=0;
							for(int j=1;j<=4;j++)
							{
								if(ELEM[i].elm[j]==0)
								{
									flag=ON;
									J=j;
									num++;
								}
							}
							//if(num>1) cout<<i<<" 静的要素が複数面で動的要素と接しています mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].boundary_condition<<NODE[ELEM[i].node[2]].boundary_condition<<NODE[ELEM[i].node[3]].boundary_condition<<NODE[ELEM[i].node[4]].boundary_condition<<endl;
							if(num>1) cout<<i<<" 静的要素が複数面で動的要素と接しています mate="<<ELEM[i].material<<"BD "<<NODE[ELEM[i].node[1]].r[A_Z]<<endl;
						}
						if(flag==ON)
						{
							BD_static_ELEM.push_back(i);
							BD_static_ELEM_elm.push_back(J);//要素iは第J面で動的要素に接している
						}
					}
					cout<<"KK="<<BD_static_ELEM.size()<<endl;

					double Pn[3];				//起点の座標

					for(int k=0;k<1;k++)
					{
						int i=BD_static_ELEM[k];
						int j=BD_static_ELEM_elm[k];
						int ia=ELEM[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
						int ib=ELEM[i].node[4-(j-1)/2*2];
						int ic=ELEM[i].node[3-(j/2%2)*2];

						double iaic[3];//ia→icのﾍﾞｸﾄﾙ成分格納
						double iaib[3];//ia→ibのﾍﾞｸﾄﾙ成分格納
						double iaicL=0;	//ia→icの長さ
						for(int D=0;D<3;D++)
						{
							iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
							iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
							iaicL+=iaic[D]*iaic[D];
						}
						iaicL=sqrt(iaicL);
						///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのﾍﾞｸﾄﾙとなる
						double S[3];//上記のﾍﾞｸﾄﾙ成分格納
						S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
						S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
						S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
						double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
						////面積Sがもとまった

						double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//法線ﾍﾞｸﾄﾙ

						double Gp[3];								//表面三角形の重心
						for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

						for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//起点の座標

						//////////////////////
						///起点を強制的に座標指定
				
						Pn[A_X]=0;
						Pn[A_Y]=0;
						Pn[A_Z]=0.225;//0.10325;		//るつぼ底部：0.10125

						//流体の最大高さを求め、それの少し上に設置

						/*//流体節点の最大高さ、最小高さ
						double Zmax=0;
						double Zmin=10;
						for(int i=1;i<=node;i++)
						{
							if(NODE[i].material==FLUID)
							{
								if(NODE[i].r[A_Z]>Zmax) Zmax=NODE[i].r[A_Z];
								if(NODE[i].r[A_Z]<Zmin) Zmin=NODE[i].r[A_Z];
							}
						}

						Pn[A_X]=0;
						Pn[A_Y]=0;
						Pn[A_Z]=Zmax-CON->get_distancebp();
						*///

						node++;
						for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
						NODE[node].boundary_condition=0;
						NODE[node].material=AIR;
						NODE[node].particleID=-1;
						NODE[node].remesh=ON;

					}
					cout<<"起点解決"<<endl;
					int countBD=0;
					for(int i=1;i<=node;i++) if(NODE[i].BD_node==ON) countBD++;
					cout<<"境界の節点数="<<countBD<<endl;

					int *imen[4];
					for(int D=0;D<4;D++) imen[D]=new int [BD_static_ELEM.size()+1];
					int *jmen =new int [BD_static_ELEM.size()+1];
					int *kmen =new int [BD_static_ELEM.size()+1];
					double *vol =new double [BD_static_ELEM.size()+1];
					int ip=node;

					/*
					if(BD_static_ELEM.size()>100000) cout<<"imenなどのメモリをあげてください"<<endl;
					int ip=node;			
					int imen[100000][3+1];	//多面体表面三角形の節点番号格納
					int jmen[100000];		//多面体表面三角形に隣接する四面体番号格納
					int kmen[100000];		//多面体表面三角形に隣接する四面体の隣接面番号 (相手は第何面で自分と接しているか)
					double vol[100000];		//多面体の体積の６倍
					*/

					for(int k=1;k<=BD_static_ELEM.size();k++)
					{
						int kelm=BD_static_ELEM[k-1];//remesh領域に接する静的要素
						int j=BD_static_ELEM_elm[k-1];
						int ia=ELEM[kelm].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
						int ib=ELEM[kelm].node[4-(j-1)/2*2];
						int ic=ELEM[kelm].node[3-(j/2%2)*2];

						/*
						imen[k][1]=ic;
						imen[k][2]=ib;
						imen[k][3]=ia;
						*/
						imen[1][k]=ic;
						imen[2][k]=ib;
						imen[3][k]=ia;

						jmen[k]=kelm;
						kmen[k]=j;//jelmがielmに接する面番号
						//vol[k]=volume3D(NODE,ia,ib,ic,ip);
						vol[k]=volume3D(NODE,ic,ib,ia,ip);	//節点の順番に注意
						//if(vol[k]<0) cout<<"体積負の要素あり "<<vol[k]<<endl;
					}

					int ibound=(int) BD_static_ELEM.size();//表面の数を表す
					int nelm0=nelm;//変更前の要素数を記憶
			
					for(int i=1;i<=ibound;i++)//要素情報生成
					{   
						nelm++;
						int ielm=nelm0+i;
						int ia=imen[1][i];
						int ib=imen[2][i];
						int ic=imen[3][i];
						ELEM[ielm].node[1]=ia;
						ELEM[ielm].node[2]=ib;
						ELEM[ielm].node[3]=ic;
						ELEM[ielm].node[4]=ip;//新点は４番目と定義	
						ELEM[ielm].elm[4]=jmen[i];
						if(jmen[i]!=0) ELEM[jmen[i]].elm[kmen[i]]=ielm;
						ELEM[ielm].volume=vol[i];
				
						sphere3D(NODE,ELEM,ia,ib,ic,ip,ielm);//外接球の中心（ボロノイ点)と半径の二乗を計算
			
						ELEM[ielm].material=AIR;//通常のpoly関数との相違点。ここで材質を空気と決定する
					}
					///////////////////

					//要素-要素関係修正/////////上の処理で第4面で接する要素番号はわかっているので、残りを求める
					//						ここで、1～3面は多面体を構成する要素との境界面であることに注意
					int ix=0;
			
					for(int i=1;i<=ibound;i++)
					{
						int ielm=nelm0+i;
						for(int j=1;j<=3;j++)//ELEM[ielm].elm[4]はすでにもとまったから、それ以外をもとめる
						{
							///ELEM[ielm].node[4]=ipである
							int ia=ELEM[ielm].node[j%3+1];		//j=1,2,3のとき、2,3,1の順
							int ib=ELEM[ielm].node[(j%3+1)%3+1];//j=1,2,3のとき、3,1,2の順
							int flag=0;
							for(int k=1;k<=ix;k++)
							{
								if(flag==0)
								{
									int ja=imen[1][k];
									int jb=imen[2][k];
									if(ia==ja && ib==jb)//節点が一致したら
									{
										ELEM[ielm].elm[j]=jmen[k];//あらかじめﾘｽﾄしてあった情報を格納
										ELEM[jmen[k]].elm[kmen[k]]=ielm;
										imen[1][k]=imen[1][ix];		//k番目の情報はもう不要。なので配列の一番最後の情報をk番目にもってきて、それまでの情報は破棄する
										imen[2][k]=imen[2][ix];
										jmen[k]=jmen[ix];
										kmen[k]=kmen[ix];
										ix--;						//待ち辺数減少
										flag=1;						//ELEM[ielm].elm[j]はもとまったので、下のネストに入る必要はないのでflag=1
									}
								}
							}
							if(flag==0)
							{
								ix++;			//ここでのixは、[隣接関係を満たす要素]をまっている[辺]の数を表す。
								imen[1][ix]=ib;	//自分の節点の並びを記憶させ、別の要素がこの並びを満たすのを待つ。ibとiaの並びを逆にしてあることに注意
								imen[2][ix]=ia;
								jmen[ix]=ielm;
								kmen[ix]=j;
							}
						}	
					}///要素-要素関係修正完了

					cout<<"要素作成完了"<<endl;

					fill3D(NODE,ELEM,nelm);

					node0=node;
					cout<<"node0="<<node0<<endl;

					///////節点情報をstatic.nodeに出力
					ofstream fout("static.node");
	
					fout<<"#node"<<endl;
					fout<<node<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=node;i++)
					{
						fout<<i<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].material<<" "<<NODE[i].material<<endl;
					}

					fout.close();
					//////*/

					
					///////要素情報をstatic.eleに出力
					ofstream fout2("static.ele");
	
					fout2<<nelm<<" "<<"4"<<" "<<"1"<<endl;
					for(int i=1;i<=nelm;i++)
					{
						fout2<<i<<" "<<ELEM[i].node[1]<<" "<<ELEM[i].node[2]<<" "<<ELEM[i].node[3]<<" "<<ELEM[i].node[4]<<" "<<ELEM[i].material<<endl;
					}

					fout2.close();

					/////////////*/

					////////////流体節点をfluid.nodeに出力
					ofstream fl("static-a.node");
					/////
					fl<<"#node"<<endl;
					fl<<fluid_number<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
					for(int i=1;i<=fluid_number;i++)
					{
						fl<<i<<" "<<PART[i-1].r[A_X]<<" "<<PART[i-1].r[A_Y]<<" "<<PART[i-1].r[A_Z]<<" "<<AIR<<" "<<AIR<<endl;
					}
					
					fl.close();
					////////

					//////static.node,eleからメッシュ再構成
					TETFUNC.call_TetGen(CONF,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
					////////

					/////////
					for(int D=0;D<4;D++) delete [] imen[D];
					delete [] jmen;
					delete [] kmen;
					delete [] vol;


					/*/境界面から層を形成
					for(int m=1;m<=1;m++)
					{
						for(int k=1;k<BD_static_ELEM.size();k++)
						{
							int i=BD_static_ELEM[k];
							int j=BD_static_ELEM_elm[k];
							int ia=ELEM[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
							int ib=ELEM[i].node[4-(j-1)/2*2];
							int ic=ELEM[i].node[3-(j/2%2)*2];

							double iaic[3];//ia→icのﾍﾞｸﾄﾙ成分格納
							double iaib[3];//ia→ibのﾍﾞｸﾄﾙ成分格納
							double iaicL=0;	//ia→icの長さ
							for(int D=0;D<3;D++)
							{
								iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
								iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
								iaicL+=iaic[D]*iaic[D];
							}
							iaicL=sqrt(iaicL);
							///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのﾍﾞｸﾄﾙとなる
							double S[3];//上記のﾍﾞｸﾄﾙ成分格納
							S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
							S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
							S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
							double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
							////面積Sがもとまった

							double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//法線ﾍﾞｸﾄﾙ

							double Gp[3];								//表面三角形の重心
							for(int D=0;D<3;D++) Gp[D]=(NODE[ia].r[D]+ NODE[ib].r[D]+NODE[ic].r[D])/3;

							//for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*iaicL;		//起点の座標 オリジナル
							for(int D=0;D<3;D++) Pn[D]=Gp[D]+n[D]*0.0001*m;		//起点の座標

							node++;
							for(int D=0;D<3;D++) NODE[node].r[D]=Pn[D];
							NODE[node].boundary_condition=0;
							NODE[node].material=AIR;
							NODE[node].particleID=-1;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;
						}
					}
					*/


					double u0=PI*4E-7;			//空気の透磁率
					double skin_depth=sqrt(1.0/(PI*CON->get_Hz()*CON->get_ele_conduc()*u0));//表皮深さ
					cout<<"流体表皮深さ="<<skin_depth<<endl;
			
					/////流体粒子を動的節点として格納
					ofstream fp("rr.dat");
					if(CON->get_thinout_fluid()==0)
					{
						/////
						for(int i=0;i<fluid_number;i++)
						{
							//if(PART[i].surface==ON  || i%4==0)
							node++;
							for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
							NODE[node].boundary_condition=0;
							NODE[node].material=FLUID;
							NODE[node].particleID=i;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;
							fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
						}
						////*/
					}
					if(CON->get_thinout_fluid()>0)
					{
						double Xf=0.0;
						double Yf=0.0;
						double Zf=0.0;
						for(int i=0;i<fluid_number;i++)
						{
							Xf+=PART[i].r[A_X];
							Yf+=PART[i].r[A_Y];
							Zf+=PART[i].r[A_Z];
						}
						Xf/=fluid_number;
						Yf/=fluid_number;
						Zf/=fluid_number;
						cout<<"流体重心=("<<Xf<<","<<Yf<<","<<Zf<<")"<<endl;
				
						for(int i=0;i<fluid_number;i++)
						{
							double rx=PART[i].r[A_X]-Xf;
							double ry=PART[i].r[A_Y]-Yf;
							double rz=PART[i].r[A_Z]-Zf;
							if(PART[i].surface==ON || sqrt(rx*rx+ry*ry+rz*rz)>CON->get_fluidwidth()*0.001-1.1*skin_depth)
							{
								node++;
								for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
								NODE[node].boundary_condition=0;
								NODE[node].material=FLUID;
								NODE[node].particleID=i;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;
								fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
							}
							else if(i%CON->get_thinout_fluid()==0)
							{
								node++;
								for(int D=0;D<3;D++) NODE[node].r[D]=PART[i].r[D];
								NODE[node].boundary_condition=0;
								NODE[node].material=FLUID;
								NODE[node].particleID=i;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;
								fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
							}
						}
					}
					fp.close();

					/*///////
					//表面の精度向上のため、流体節点を追加
					double Rz=0.13125-0.003;
					//まずは半球を作る。そのためには半球表面を作成する必要がある。
					for(int Ri=1;Ri<=5;Ri++)
					{
						double R=0.025-5*0.0001+Ri*0.0002;
						double le=CON->get_distancebp();

						double A=sqrt(3.0)/2;				//よく使う係数
						double B=sqrt(2.0/3);						////よく使う係数
						int half_WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
						int half_WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
						int half_WZ=(int)(R/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
						double R2=R-0.5*le;				//少し小さめの半径を設定

						///////////半球表面
						int Nt;						//球表面の、θ方向の分割数
						double Lt;					//球表面の、θ方向の分割距離
						//calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//半円の分割数は偶数・奇数どちらでもよい
						double temp_N=PI/2*R/(le*A);			//仮の分割数。leで割り切れたら一番いいけど、そうもいかないときがある
						int Ns=(int) temp_N;				//真の分割数
						double difference=temp_N-Ns;		//仮と真の差
						if(difference>0.5) Ns++;
						Lt=PI/2*R/Ns;			//粒子の距離
						Nt=Ns;

						double d_theta=Lt/R;		//弧の長さがLtになる角度

						for(int k=0;k<=Nt;k++)//loopはk<Ntで終わらせる。Ntに該当するところはすでに設置済み
						{
							double THETA=k*d_theta;	//θ
							double r=R*sin(THETA);	//その高さにおける円の半径
							double round=2*PI*r;//その高さにおける円周

							//int Nf=calc_division_N_circle(round,le);//球表面の、θ方向の分割数
							//対称性を考えた時、円周を記述する粒子は偶数個でなければならない。だから他の辺分割数とは扱いが少し特殊
							//dis:分割する距離(円周)
							double temp_num=round/le;		//円外周に設置する『仮の』粒子数。ただし外周がうまくleで割り切れるとは限らない

							int N1=(int)(temp_num/2);
							N1*=2;							
							int N2=N1+2;					//temp_numはN1とN2の間にある。ここでN1,N2は偶数

							double dif1=temp_num-N1;		//各Nとの差
							double dif2=N2-temp_num;
							int N=N1;						//周方向分割数
							if(dif2<dif1) N=N2;				//差の小さい方をNとして採用する。

							int Nf=N;
					
					
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
								double xf=r*cos(fai);
								double yf=r*sin(fai);
								double zf=R*cos(THETA);
						
								node++;
								NODE[node].r[A_X]=xf;
								NODE[node].r[A_Y]=yf;
								NODE[node].r[A_Z]=zf+Rz;
								NODE[node].boundary_condition=0;
								//NODE[node].material=FLUID;
								if(Ri<=3) NODE[node].material=FLUID;
								if(Ri>3) NODE[node].material=AIR;
								NODE[node].particleID=0;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;

								if(k!=Nt)
								{
									node++;
									NODE[node].r[A_X]=xf;
									NODE[node].r[A_Y]=yf;
									NODE[node].r[A_Z]=-zf+Rz;
									NODE[node].boundary_condition=0;
									//NODE[node].material=FLUID;
									if(Ri<=3) NODE[node].material=FLUID;
									if(Ri>3) NODE[node].material=AIR;
									NODE[node].particleID=0;
									NODE[node].remesh=ON;
									NODE[node].BD_node=OFF;
								}
							}
						}
						if(Nt%2!=0)//Ntが奇数のときは、頂上に粒子が置かれなければならない。しかし上のloopはそれが不可能。よってここで追加
						{
							node++;
							NODE[node].r[A_X]=0;
							NODE[node].r[A_Y]=0;
							NODE[node].r[A_Z]=R+Rz;
							NODE[node].boundary_condition=0;
							//NODE[node].material=FLUID;
							if(Ri<=3) NODE[node].material=FLUID;
							if(Ri>3) NODE[node].material=AIR;
							NODE[node].particleID=0;
							NODE[node].remesh=ON;
							NODE[node].BD_node=OFF;

					
								node++;
								NODE[node].r[A_X]=0;
								NODE[node].r[A_Y]=0;
								NODE[node].r[A_Z]=-R+Rz;
								NODE[node].boundary_condition=0;
								//NODE[node].material=FLUID;
								if(Ri<=3) NODE[node].material=FLUID;
								if(Ri>3) NODE[node].material=AIR;
								NODE[node].particleID=0;
								NODE[node].remesh=ON;
								NODE[node].BD_node=OFF;
				
						}
						//////////////////////////////////////
					}

					//*////
					unsigned int timeD=GetTickCount();
			
					/*
					int checkmax=0;
					int checkmin=node;
					for(int i=1;i<=node;i++)
					{
						if(NODE[i].BD_node==ON)
						{
							if(i>checkmax) checkmax=i;
							if(i<checkmin) checkmin=i;
						}
					}
					cout<<"境界節点番号の最大値、最小値＝"<<checkmax<<","<<checkmin<<endl;

					int checkbig=0;
					int checksmall=0;
					double check=node/2.0;
					for(int i=1;i<=node;i++)
					{
						if(NODE[i].BD_node==ON)
						{
							if(i>check) checkbig++;
							if(i<check) checksmall++;
						}
					}
					cout<<"境界節点番号：総数の半分と比較した大小＝"<<checkbig<<","<<checksmall<<endl;

					*/
			
					cout<<"node="<<node<<" デローニ分割"<<endl;
				}
			}
			
			///////////////////////////////////////////////////////////////////////////////

			int N=(int)TRANS.size()-1;			//FEM節点に含まれる粒子数(0番目を除くTRANS[]の長さ)
			node=(int)NODEall.size();		//節点数
			nelm=(int)ELEMall.size();		//要素数
			KTJ=node;						//最大節点数
			KTE=nelm;						//最大要素数
			//int *depth=new int [KTE+1];			//各要素の深さ格納

			cout<<"tetgen後のnode="<<node<<endl;
			cout<<"tetgen後のnelm="<<nelm<<endl;

			///配列確保
			//point3D *NODE2=new point3D [NODEall.size()+1];
			//element3D *ELEM2=new element3D [ELEMall.size()+1];
			
			//NODE.clear();
			//ELEM.clear();
			point3D NODE0;
			element3D ELEM0;
			edge3D EDGE0;
			//int ID;
			
			for(int i=0;i<=node+1;i++) NODE.push_back(NODE0);
			for(int i=0;i<=nelm+1;i++) ELEM.push_back(ELEM0);	//配列を確保
			for(int i=0;i<2*KTE;i++) EDGE.push_back(EDGE0);	//配列を確保
			//NODE.resize(NODEall.size()+1);
			//ELEM.resize(ELEMall.size()+1);
			
			for(int i=0;i<node;i++) if(NODEall[i].attribute==WATER) NODEall[i].attribute=FLUID;//tetgenでは流体をWATERとしているので、FLUIDに変更
			for(int i=0;i<nelm;i++) if(ELEMall[i].attribute==WATER) ELEMall[i].attribute=FLUID;//tetgenでは流体をWATERとしているので、FLUIDに変更

			for(int i=0;i<node;i++) if(NODEall[i].attribute==CONDUCT) NODEall[i].attribute=CRUCIBLE;
			for(int i=0;i<nelm;i++) if(ELEMall[i].attribute==CONDUCT) ELEMall[i].attribute=CRUCIBLE;
			
			//TetGenの節点・要素データを取得	節点番号を1つずらす
			//節点データ	
			for(int i=0;i<node;i++)
			{
				NODE[i+1].r[A_X]=NODEall[i].r[A_X];
				NODE[i+1].r[A_Y]=NODEall[i].r[A_Y];
				NODE[i+1].r[A_Z]=NODEall[i].r[A_Z];
				NODE[i+1].material=(int)NODEall[i].attribute;
			}
			//要素データ
			for(int i=0;i<nelm;i++)
			{
				ELEM[i+1].material=ELEMall[i].attribute;								//材質
				for(int n=0;n<4;n++)	ELEM[i+1].node[n+1]=ELEMall[i].node[n]+1;		//構成節点
				for(int n=0;n<4;n++)													//要素-要素関係 
				{
					if(ELEMall[i].nei_elem[n]==-1)	ELEM[i+1].elm[n+1]=0;	//近隣要素なしの場合TetGenは-1を返してくるため0に修正する
					else							ELEM[i+1].elm[n+1]=ELEMall[i].nei_elem[n]+1;
				}

			}

			//体積および外接球パラメータ計算
			for(int i=1;i<=nelm;i++)
			{
				//4つの節点
				int ia=ELEM[i].node[1];
				int ib=ELEM[i].node[2];
				int ic=ELEM[i].node[3];
				int ip=ELEM[i].node[4];

				//要素体積計算
				ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);
				
				///体積チェック
				if(ELEM[i].volume<0)			{cout<<"i="<<i<<" 体積が負です"<<endl;				cout<<"volume="<<ELEM[i].volume<<endl;}
				else if(ELEM[i].volume==0)		{cout<<"i="<<i<<" 体積が0です"<<endl;				cout<<"volume="<<ELEM[i].volume<<endl;}
				//else if(ELEM[i].volume<1e-20)	{cout<<"i="<<i<<" 体積が小さすぎます 1e-20"<<endl;	cout<<"volume="<<ELEM[i].volume<<endl;}
				//*/
				//if(ELEM[i].volume==0)		{cout<<"i="<<i<<" 体積が0です"<<endl;	ELEM[i].volume=1e-50;}
				//if(ELEM[i].volume<min_volume && ELEM[i].volume!=0)	min_volume=ELEM[i].volume;

				//外接球中心座標および半径計算
				sphere3D(NODE,ELEM,ia,ib,ic,ip,i);
			}

			////モデル毎の設定(境界条件など)
			//静電霧化
			if(CON->get_model_number()==14)
			{
				for(int i=0;i<node;i++)
				{
					//iが1ずれていることに注意
					if(NODEall[i].boundary==ELECTRODE1)			NODE[i+1].boundary_condition=1;	//円柱電極および土台
					else if(NODEall[i].boundary==ELECTRODE2)	NODE[i+1].boundary_condition=2;	//平板電極
					else										NODE[i+1].boundary_condition=0;	//その他の節点は未知数
				}
			}
			if(CON->get_model_number()==25)
			{
				for(int i=1;i<=node;i++)
				{
							
					if(NODE[i].r[A_Z]<CON->get_ZD()+err) NODE[i].boundary_condition=1;		//-Z面
					else if(NODE[i].r[A_Z]>CON->get_ZU()-err) NODE[i].boundary_condition=1;		//+Z面
					else if(NODE[i].r[A_X]<CON->get_XL()+err) NODE[i].boundary_condition=1;		//-X面
					else if(NODE[i].r[A_X]>CON->get_XR()-err) NODE[i].boundary_condition=1;		//+X面
					else if(NODE[i].r[A_Y]<CON->get_YD()+err) NODE[i].boundary_condition=1;		//-Y面
					else if(NODE[i].r[A_Y]>CON->get_YU()-err) NODE[i].boundary_condition=1;		//+Y面
					else
					{
						NODE[i].boundary_condition=0;					//解析境界節点ではない
					}
				}
			}
	//cout<<"最小体積="<<min_volume<<endl;
		}
	}
	/////
	int countf=0;
	int countco=0;
	int countcru=0;
	int counta=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].material==FLUID) countf++;
		if(NODE[i].material==COIL) countco++;
		if(NODE[i].material==CRUCIBLE) countcru++;
		if(NODE[i].material==AIR) counta++;
	}

	if(countf>=0) cout<<"流体節点="<<countf<<endl;
	if(countco>=0) cout<<"コイル節点="<<countco<<endl;
	if(countcru>=0) cout<<"るつぼ節点="<<countcru<<endl;
	if(counta>=0) cout<<"空気節点="<<counta<<endl;
	///*/

	/*///るつぼ最大高さ
	double hmax=0;
	double hmin=10000;
	for(int i=1;i<=node;i++)
	{
		
		if(NODE[i].material==CRUCIBLE)
		{
			if(NODE[i].r[A_Z]>hmax) hmax=NODE[i].r[A_Z];
			if(NODE[i].r[A_Z]<hmin) hmin=NODE[i].r[A_Z];
		}
	}
	cout<<"るつぼZ座標 最大="<<hmax<<" 最小="<<hmin<<endl;
	///*/

	
	
	
	
	/////節点-要素関係
    int *jnb=new int[node+1];///各節点に隣接する要素数格納
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//各節点の周辺要素番号格納
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	//メッシュスムージング
	//if(CON->get_mesh_smth()>0) for(int i=0;i<CON->get_mesh_smth();i++) laplacian_smoothing2(NODE,ELEM,&node,&nelm,CON,jnb,nei,node0);
	
	/////メッシュ生成を確認
	//if(TIME==0)
	//if(t==1 || (t-1)%(CON->get_EM_interval()*CON->get_mesh_output_interval())==0) 
	if(t==1)
	{
		data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);//断面図。動画ファイルがつくられる
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
		data_avs2node(CON,node,nelm,NODE,ELEM,val,t);
		data_avs3(node,nelm,NODE,ELEM,CON,t);//材質
		delete [] val;
	}
	//*/

	countf=0;
	countco=0;
	countcru=0;
	counta=0;
	for(int i=1;i<=node;i++)
	{
		if(jnb[i]==0)
		{
			//cout<<"jnb=0 i="<<i<<" material="<<NODE[i].material<<endl;
			NODE[i].boundary_condition=1;//境界条件をﾃﾞｨﾘｸﾚ型にすることで、ICCGに参加させない
			if(NODE[i].material==FLUID)
			{
				NODE[i].particleID=-1;//流体節点が消失する場合、必要な値が計算されないため、粒子にフィードバックできない(してはいけない)
				countf++;
			}
			if(NODE[i].material==COIL) countco++;
			if(NODE[i].material==CRUCIBLE) countcru++;
			if(NODE[i].material==AIR) counta++;
		}
	}

	if(countf>0) cout<<"要素をつくらない流体節点あり 数は"<<countf<<endl;
	if(countco>0) cout<<"要素をつくらないコイル節点あり 数は"<<countco<<endl;
	if(countcru>0) cout<<"要素をつくらないるつぼ節点あり 数は"<<countcru<<endl;
	if(counta>0) cout<<"要素をつくらない空気節点あり 数は"<<counta<<endl;

	int node_sta=(int) static_NODE.size()-1;
	cout<<"static_NODE.size()="<<node_sta<<endl;
	//流体の体積調査
	double volume_f=0;//流体要素の体積の合計
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID)
		{
			int ia=ELEM[i].node[1];
			int ib=ELEM[i].node[2];
			int ic=ELEM[i].node[3];
			int ip=ELEM[i].node[4];
		
			double volme=volume3D(NODE,ia,ib,ic,ip);
			//if(volme<=0) cout<<"体積がゼロか負"<<endl;
			volume_f+=volme;
			
			//ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);//体積の6倍であることに注意
			//sphere3D(NODE,ELEM,ia,ib,ic,ip,i);//外接球の中心（ボロノイ点)と半径の二乗を計算
		}
	}
	volume_f/=6;
	cout<<"流体体積="<<volume_f<<endl;
	
	if(t==1)
	{
		ofstream fv("volume.dat");
		fv.close();
	}
	ofstream vo("volume.dat",ios :: app);
	vo<<t<<" "<<volume_f<<endl;
	vo.close();

	//////節点番号の並び替え
	if(CON->get_node_sort()==1) node_sorting(CON,node,nelm,NODE,ELEM,jnb,nei);
	else if(CON->get_node_sort()==2) node_sorting2(CON,node,nelm,NODE,ELEM,jnb,nei);
	//////
	

	/////節点-要素関係 節点番号が変更されたのでもう一度求める
    int *jnb3=new int[node+1];///各節点に隣接する要素数格納
    set_jnb3D(NODE,ELEM,node,nelm,jnb3);
    int **nei3=new int* [node+1];
    for(int i=1;i<=node;i++) nei3[i]=new int [jnb[i]+1];//各節点の周辺要素番号格納
    set_nei3D(NODE,ELEM,node,nelm,jnb3,nei3);


	//楕円方程式
	if(CON->get_EM_calc_type()==1 || CON->get_EM_calc_type()==4) potential_calculation(CON,NODE,ELEM, node, nelm,jnb3, TIME,PART, fluid_number,nei3);
	if(CON->get_EM_calc_type()==3) calc_transitional_EM_field(CON, node, nelm,nedge,NODE,ELEM,EDGE,jnb3, dt, TIME,PART, fluid_number,F,t,nei3,particle_node,NODE_jw,ELEM_jw,node_sta,static_EDGE);
		
	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
	delete [] nei;

	delete [] jnb3;
    for(int i=1;i<=node;i++) delete [] nei3[i];
    delete [] nei3;
	return 0;
}

//渦電流を考慮にいれた節点要素のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算関数ver.2 ｸｰﾛﾝｹﾞｰｼﾞを適用して左辺をﾗﾌﾟﾗｼｱﾝとする
void Avector3D_node_eddy2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t)
{   
	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	double j0x,j0y,j0z;					//電流密度
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num;///未知数
    int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
    int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=node;i++)//行列はA1x,A1y,A1z,A1φ,A2x,A2y・・・の順に格納
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			npp[i]=num;
			num+=3;

			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-3;
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-3;
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///行列の最大幅計算  もっとも幅の大きいのはφの要素だろうと仮定している。確実ではないことに注意
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=4*mat_w;//X,Y,Z,φ成分とあるので×4 //左辺をﾗﾌﾟﾗｼｱﾝにしても、φに関しては４成分のままに注意

	////配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    //*////

	////行列の幅計算 各行ごとに幅を求め、メモリを削減する
	int mat_w=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=4*mat_w;//
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////配列確保
	double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

	int width_max=0;
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
		width_max+=width_mat[i];
    }

	cout<<"合計幅="<<width_max<<endl;

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;

	delete [] width_mat;	
	delete [] width_node;	

	//*/
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
   
    
    /////////全体行列を作成する
	cout<<"全体行列作成開始"<<endl;
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			j0x=current[A_X][je];
			j0y=current[A_Y][je];
			j0z=current[A_Z][je];
		}
		else
		{
			j0x=0;
			j0y=0;
			j0z=0;
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=delta/4*j0x;//x,y,zは１つずつずれる
				B[I]+=delta/4*j0y;//
				B[I+1]+=delta/4*j0z;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flag=1;
							}
						}
						if(flag==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=J;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y方向
				I=npp[N[n]]+1+1;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
				/////*/

				/////Z方向
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						J3=J2+1;		//Aizの項
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J3;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	//////渦電流項計算
	int J4,flag4;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sigma;//σV/(20dt)　何度も計算することになるので係数化する	
			//double co2=delta6*sigma;//   σ/(36V)
			//double co3=sigma*dt*delta6;

			double co=delta/20.0/dt*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			double co2=delta6*sig[je];//   σ/(36V)
			double co3=sig[je]*dt*delta6;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J4=J+3;			//φの項
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aixの項
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==J) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag=1;
								}
							
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flag4=1;
								}///
							}
							if(flag==0)//相当するＪが存在しなかったら作る(ここではそんなはずないけど)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==J) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=J4;
							}

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							if(I==J) Sx+=2*Ax;
							else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Y方向

					I=npp[N[n]]+1+1;/////Y方向
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J4=J+3;			//φの項
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==J2) G[I][h]+=co*2;//Aiyの項
									else G[I][h]+=co;
									flag2=1;
								}

								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(I==J2) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J2;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=J4;
							}//*/

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							if(I==J2) Sy+=2*Ay;
							else Sy+=Ay;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sy;//支配方程式のX成分から得られるBの値

					/////Z方向
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==J3) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    if(I==J3) G[I][H]+=co*2;
								else G[I][H]+=co;
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=J4;
							}//*/
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							if(I==J3) Sz+=2*Az;
							else Sz+=Az;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sz;//支配方程式のX成分から得られるBの値

					/////φ
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aixの項
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flag=1;
								}
								//Aiyの項
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flag2=1;
								}
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=J4;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}//////*/

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	check_matrix_symmetry(pn,NUM,ROW,G);

	//対角成分の値チェック
	cout<<"対角成分の値チェック"<<endl;
	int count_plus=0;
	int count_minus=0;
	int count_zero=0;
	for(int i=1;i<=pn;i++)
	{
		int flag=OFF;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			if(i==J)
			{
				if(G[i][j]==0)
				{
					count_zero++;
					if(ppn[i-1]==-6) cout<<"対角成分が零(φi） i="<<i<<endl;
					else if(ppn[i-1]==-7) cout<<"対角成分が零(φr） i="<<i<<endl;
					else cout<<"対角成分が零(A） i="<<i<<endl;
				}
				if(G[i][j]>0) count_plus++;
				if(G[i][j]<0) count_minus++;
				flag=ON;

			}
		}
		if(flag==OFF) count_zero++;
	}
	cout<<"正の対角項の数="<<count_plus<<endl;
	cout<<"負の対角項の数="<<count_minus<<endl;
	cout<<"零の対角項の数="<<count_zero<<endl;
	///

	//対角成分の分布
	diagonal_distribution(CON,node,nelm,NODE,ELEM,pn,NUM,ROW,G);

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

	cout<<"非零数="<<number<<endl;


	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	//*///

	//行列の視覚化
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	

	//CG法実行
	double *XX=new double [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) A[A_X][i]=XX[n];
		else if(i==-1) A[A_Y][ppn[n-1]]=XX[n];
		else if(i==-2) A[A_Z][ppn[n-2]]=XX[n];
		else if(i==-3) V[ppn[n-3]]=XX[n];//電位φ
	}	
	
	delete [] XX;

	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///渦電流と強制電流を足した値を出力させる
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_Aの更新 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_Aを粒子に渡す
	cout<<"old_Aを粒子に記憶";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_Aのﾌｧｲﾙ出力
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

////////////
////////////
////////////
//渦電流を考慮にいれた節点要素のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算関数 jω法 節点ごとに実部、虚部をそれぞれ格納
void Avector3D_node_eddy2_jw(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{   
	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始（ｊω法）"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=(3*(node-NN)+conducter_num)*2;///未知数 複素数なので各成分に実部と虚部がある
    int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
    int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=node;i++)//行列はRe(A1x),Im(A1x),Re(A1y),Im(A1y),Re(A1z),Re(φ1),Im(φ1),Re(A2x),・・・の順に格納
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			ppn[num+3]=-3;
			ppn[num+4]=-4;
			ppn[num+5]=-5;
			npp[i]=num;
			num+=6;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-6;
						ppn[num+1]=-7;
						num+=2;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-6;
						ppn[num+1]=-7;
						num+=2;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///行列の最大幅計算  もっとも幅の大きいのはφの要素だろうと仮定している。確実ではないことに注意
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=8*mat_w;//X,Y,Z,φ成分にそれぞれRe,Imがあるので×8
    *//////

	/*///配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
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
    //*/
    

	////行列の幅計算 各行ごとに幅を求め、メモリを削減する
	int mat_w=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=8*mat_w;//X,Y,Z,φ成分にそれぞれRe,Imがあるので×8
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=8;
				for(int j=1;j<=8;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=8;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=6;
				for(int j=1;j<=6;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=6;
			}
		}
	}	
    //////

	////配列確保
	double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;

	delete [] width_mat;	
	delete [] width_node;	

    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    ////
    
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	//int J,J2,J3,flag,flag2,flag3;
	int Jxr,Jxi,Jyr,Jyi,Jzr,Jzi,flagxr,flagxi,flagyr,flagyi,flagzr,flagzi;
    for(int je=1;je<=nelm;je++)
    {   
		//cout<<je<<endl;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0);
			j0y=complex<double> (current[A_Y][je],0);
			j0z=complex<double> (current[A_Z][je],0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向、実部

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=delta/4*j0x.real();//Rex
				B[I]+=delta/4*j0x.imag();//Imx
				B[I+1]+=delta/4*j0y.real();//Rey
				B[I+2]+=delta/4*j0y.imag();//Imy
				B[I+3]+=delta/4*j0z.real();//Rez
				B[I+4]+=delta/4*j0z.imag();//Imz
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxr=npp[N[m]]+1;	//Re(Aix)の項
						flagxr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flagxr=1;
							}
						}
						if(flagxr==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=Jxr;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////X方向、虚部

				I=npp[N[n]]+2;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxi=npp[N[m]]+2;	//Im(Aix)の項
						flagxi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flagxi=1;
							}
						}
						if(flagxi==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=Jxi;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}


				/////Y方向,実部
				I=npp[N[n]]+3;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyr=npp[N[m]]+3;			//Re(Aiy)の項
						flagyr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flagyr=1;
							}
						}
						if(flagyr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyr;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y方向,虚部
				I=npp[N[n]]+4;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyi=npp[N[m]]+4;			//Im(Aiy)の項
						flagyi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flagyi=1;
							}
						}
						if(flagyi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyi;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
				/////*/

				/////Z方向,実部
				I=npp[N[n]]+5;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzr=npp[N[m]]+5;	//Aixの項
						flagzr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flagzr=1;
							}
						}
						if(flagzr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzr;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Z方向,虚部
				I=npp[N[n]]+6;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzi=npp[N[m]]+6;	//Aixの項
						flagzi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flagzi=1;
							}
						}
						if(flagzi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzi;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	//////渦電流項計算
	//int J4,flag4;
	int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==2)
		{
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			//double co2=delta6*sig[je];//   σ/(36V)
			//double co3=sig[je]*dt*delta6;

			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=(delta/20.0)*omega*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			double co2=delta6*sig[je];//   σ/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向,実部

					int I=npp[N[n]]+1;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)の項
							Jxi=npp[N[m]]+2;	//Im(Aix)の項
							Jvr=npp[N[m]]+7;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)の項
								if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jxi)
									{
										//cout<<"対角項に足されている？"<<endl;
										G[I][h]-=co*2;//起こらないはず
									}
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagxi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvr=1;
								}///
							}
							if(flagxi==0)
							{  
								
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jxi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////X方向,虚部

					I=npp[N[n]]+2;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1;	//Re(Aix)の項
							//Jxi=npp[N[m]]+2;	//Im(Aix)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+8;	//Imφの項
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)の項
								if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jxr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagxr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvi=1;
								}///
							}
							if(flagxr==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jxr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Y方向,実部

					I=npp[N[n]]+3;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jyr=npp[N[m]]+3;	//Re(Aiy)の項
							Jyi=npp[N[m]]+4;	//Im(Aiy)の項
							Jvr=npp[N[m]]+7;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiy)の項
								if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jyi) G[I][h]-=co*2;//起こらないはず
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagyi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvr=1;
								}///
							}
							if(flagyi==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jyi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Ｙ方向,虚部

					I=npp[N[n]]+4;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jyr=npp[N[m]]+3;	//Re(Aiy)の項
							//Jyi=npp[N[m]]+4;	//Im(Aiy)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+8;	//Imφの項
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiy)の項
								if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jyr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagyr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvi=1;
								}///
							}
							if(flagyr==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jyr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値//ｊω法だとBが増えない

					/////Z方向,実部

					I=npp[N[n]]+5;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jzr=npp[N[m]]+5;	//Re(Aiz)の項
							Jzi=npp[N[m]]+6;	//Im(Aiz)の項
							Jvr=npp[N[m]]+7;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiz)の項
								if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jzi) G[I][h]-=co*2;//起こらないはず
									else G[I][h]-=co;//G=co(-Ai+jAr)
									flagzi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvr=1;
								}///
							}
							if(flagzi==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jzi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Ｚ方向,虚部

					I=npp[N[n]]+6;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jzr=npp[N[m]]+5;	//Re(Aiz)の項
							//Jzi=npp[N[m]]+6;	//Im(Aiz)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+8;	//Imφの項
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiz)の項
								if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jzr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagzr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvi=1;
								}///
							}
							if(flagzr==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jzr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る(すでにあるはず）
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値//ｊω法だとBが増えない

					/////φ,実部 //P77式(4.27)にdtをかける co3がかかる項の虚実に注意
					
					I=npp[N[n]]+7;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1;	//Re(Aix)の項
							Jxi=npp[N[m]]+2;	//Im(Aix)の項
							Jyr=npp[N[m]]+3;	//Re(Aiy)の項
							Jyi=npp[N[m]]+4;	//Im(Aiy)の項
							Jzr=npp[N[m]]+5;	//Re(Aiz)の項
							Jzi=npp[N[m]]+6;	//Im(Aiz)の項
							Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+8;	//Imφの項

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)の項
								if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxr=1;
								}
								//Re(Aiy)の項
								if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyr=1;
								}
								//Re(Aiz)の項
								if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzr=1;
								}
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvi=1;
								}
							}
							if(flagxr==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxr;
							}
							if(flagyr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyr;
							}
							if(flagzr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzr;
							}
							
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvi;
							}
							///Bの計算
							//他の項と同じく、old_Aが必要ないのでBは更新されない
							//仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);

					/////φ,虚部
					
					I=npp[N[n]]+8;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)の項
							Jxi=npp[N[m]]+2;	//Im(Aix)の項
							//Jyr=npp[N[m]]+3;	//Re(Aiy)の項
							Jyi=npp[N[m]]+4;	//Im(Aiy)の項
							//Jzr=npp[N[m]]+5;	//Re(Aiz)の項
							Jzi=npp[N[m]]+6;	//Im(Aiz)の項
							Jvr=npp[N[m]]+7;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)の項
								if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxi=1;
								}
								//Aiyの項
								if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyi=1;
								}
								//Aizの項
								if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzi=1;
								}
								//Reφの項 G=α(φi-jφr)
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvr=1;
								}
							}
							if(flagxi==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxi;
							}
							if(flagyi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyi;
							}
							if(flagzi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzi;
							}
							
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvr;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}//////*/
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	cout<<"G,ROWの順番並び替え開始"<<endl;
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック //非対称行列なのでチェックしない
	//check_matrix_symmetry(pn,NUM,ROW,G);

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

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	//*///

	//行列の視覚化
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	

	//CG法実行
	double *XX=new double [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);

	//XXに格納された解を各節点に振る
	//complex<double> *Ac[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ（複素数）
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) AR[A_X][i]=XX[n];
		else if(i==-1) AI[A_X][ppn[n-1]]=XX[n];
		else if(i==-2) AR[A_Y][ppn[n-2]]=XX[n];
		else if(i==-3) AI[A_Y][ppn[n-3]]=XX[n];
		else if(i==-4) AR[A_Z][ppn[n-4]]=XX[n];
		else if(i==-5) AI[A_Z][ppn[n-5]]=XX[n];
		else if(i==-6) VR[ppn[n-6]]=XX[n];//電位φ
		else if(i==-7) VI[ppn[n-6]]=XX[n];//電位φ
	}	

	delete [] XX;
	
	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	double Am=0.0;
	double phi=0.0;
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
			phi=atan(AI[D][i]/AR[D][i]);
			A[D][i]=Am*cos(omega*t+phi);
		}
	}

	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///渦電流と強制電流を足した値を出力させる
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_Aの更新 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_Aを粒子に渡す
	cout<<"old_Aを粒子に記憶";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_Aのﾌｧｲﾙ出力
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

//////////
//jω法,虚数→実数の順に格納
void Avector3D_node_eddy2_jw2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{   
	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始（ｊω法-2）"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=(3*(node-NN)+conducter_num)*2;///未知数 複素数なので各成分に実部と虚部がある
    int pnI=3*(node-NN)+conducter_num;//未知数のうち虚数部の数
	int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
    int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=node;i++)//行列はIm(A1x),Im(A1y),Im(A1z),Im(φ1),Im(A2x),・・・,Im(φn),Re(A1x),Re(A1y)の順に格納
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;//Imx
			ppn[num+1]=-1;//Imy
			ppn[num+2]=-2;//Imz
			ppn[num+pnI]=-3;//Rex
			ppn[num+1+pnI]=-4;//Rey
			ppn[num+2+pnI]=-5;//Rez
			npp[i]=num;
			num+=3;//4つめ以降は実部として分離されているので注意
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-6;//Imφ
						ppn[num+pnI]=-7;
						num+=1;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-6;//Imφ
						ppn[num+pnI]=-7;
						num+=1;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///行列の最大幅計算  もっとも幅の大きいのはφの要素だろうと仮定している。確実ではないことに注意
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=8*mat_w;//X,Y,Z,φ成分にそれぞれRe,Imがあるので×8
    //////

	////配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
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
    //*/
    

	////行列の幅計算 各行ごとに幅を求め、メモリを削減する
	int mat_w=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=8*mat_w;//X,Y,Z,φ成分にそれぞれRe,Imがあるので×8

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=8;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=4;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=6;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=3;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////配列確保
	double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }

	//delete [] width_mat;	
	//delete [] width_node;	

    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    //*/
    
	cout<<"全体行列作成開始 pn="<<pn<<" pnI="<<pnI<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	//int J,J2,J3,flag,flag2,flag3;
	int Jxr,Jxi,Jyr,Jyi,Jzr,Jzi,flagxr,flagxi,flagyr,flagyi,flagzr,flagzi;
    for(int je=1;je<=nelm;je++)
    {   
		//if(je%10==0) cout<<je<<endl;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0);
			j0y=complex<double> (current[A_Y][je],0);
			j0z=complex<double> (current[A_Z][je],0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向、虚部

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=delta/4*j0x.imag();//Imx
				B[I]+=delta/4*j0y.imag();//Imy
				B[I+1]+=delta/4*j0z.imag();//Imz
				B[I-1+pnI]+=delta/4*j0x.real();//Rex
				B[I+pnI]+=delta/4*j0y.real();//Rey
				B[I+1+pnI]+=delta/4*j0z.real();//Rez
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxi=npp[N[m]]+1;	//Im(Aix)の項
						flagxi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flagxi=1;
							}
						}
						if(flagxi==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=Jxi;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y方向,虚部
				I=npp[N[n]]+2;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyi=npp[N[m]]+2;			//Im(Aiy)の項
						flagyi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flagyi=1;
							}
						}
						if(flagyi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyi;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Z方向,虚部
				I=npp[N[n]]+3;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzi=npp[N[m]]+3;	//Aixの項
						flagzi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flagzi=1;
							}
						}
						if(flagzi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzi;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////X方向、実部

				I=npp[N[n]]+1+pnI;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
						flagxr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flagxr=1;
							}
						}
						if(flagxr==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=Jxr;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}


				/////Y方向,実部
				I=npp[N[n]]+2+pnI;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyr=npp[N[m]]+2+pnI;			//Re(Aiy)の項
						flagyr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flagyr=1;
							}
						}
						if(flagyr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyr;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				
				/////*/

				/////Z方向,実部
				I=npp[N[n]]+3+pnI;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzr=npp[N[m]]+3+pnI;	//Aixの項
						flagzr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flagzr=1;
							}
						}
						if(flagzr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzr;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	cout<<"ベクトルポテンシャル項完了"<<endl;

	//////渦電流項計算
	//int J4,flag4;
	int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==2)
		{
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			//double co2=delta6*sig[je];//   σ/(36V)
			//double co3=sig[je]*dt*delta6;

			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=(delta/20.0)*omega*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			double co2=delta6*sig[je];//   σ/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向,虚部

					int I=npp[N[n]]+1;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
							//Jxi=npp[N[m]]+2;	//Im(Aix)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)の項
								if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jxr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagxr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvi=1;
								}///
							}
							if(flagxr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jxr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Ｙ方向,虚部

					I=npp[N[n]]+2;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)の項
							//Jyi=npp[N[m]]+4;	//Im(Aiy)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiy)の項
								if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jyr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagyr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvi=1;
								}///
							}
							if(flagyr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jyr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;

					/////Ｚ方向,虚部

					I=npp[N[n]]+3;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)の項
							//Jzi=npp[N[m]]+6;	//Im(Aiz)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiz)の項
								if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jzr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagzr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvi=1;
								}///
							}
							if(flagzr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jzr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値//ｊω法だとBが増えない


					/////φ,虚部//co3がかかる項の虚実に注意
					I=npp[N[n]]+4;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
							Jxi=npp[N[m]]+1;	//Im(Aix)の項
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)の項
							Jyi=npp[N[m]]+2;	//Im(Aiy)の項
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)の項
							Jzi=npp[N[m]]+3;	//Im(Aiz)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)の項
								if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxi=1;
								}
								//Aiyの項
								if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyi=1;
								}
								//Aizの項
								if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzi=1;
								}
								//Reφの項 G=α(φi-jφr)
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvr=1;
								}
							}
							if(flagxi==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxi;
							}
							if(flagyi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyi;
							}
							if(flagzi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzi;
							}
							
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvr;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);

					/////X方向,実部

					 I=npp[N[n]]+1+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)の項
							Jxi=npp[N[m]]+1;	//Im(Aix)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)の項
								if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jxi)
									{
										//cout<<"対角項に足されている？"<<endl;
										G[I][h]-=co*2;//起こらないはず
									}
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagxi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvr=1;
								}///
							}
							if(flagxi==0)
							{  
								
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jxi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					

					/////Y方向,実部

					I=npp[N[n]]+2+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jyr=npp[N[m]]+3;	//Re(Aiy)の項
							Jyi=npp[N[m]]+2;	//Im(Aiy)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiy)の項
								if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jyi) G[I][h]-=co*2;//起こらないはず
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagyi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvr=1;
								}///
							}
							if(flagyi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jyi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					//支配方程式のX成分から得られるBの値//ｊω法だとBが増えない

					/////Z方向,実部

					I=npp[N[n]]+3+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jzr=npp[N[m]]+5;	//Re(Aiz)の項
							Jzi=npp[N[m]]+3;	//Im(Aiz)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiz)の項
								if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jzi) G[I][h]-=co*2;//起こらないはず
									else G[I][h]-=co;//G=co(-Ai+jAr)
									flagzi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvr=1;
								}///
							}
							if(flagzi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jzi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					
					/////φ,実部 //co3がかかる項の虚実に注意
					
					I=npp[N[n]]+4+pnI;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
							Jxi=npp[N[m]]+1;	//Im(Aix)の項
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)の項
							Jyi=npp[N[m]]+2;	//Im(Aiy)の項
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)の項
							Jzi=npp[N[m]]+3;	//Im(Aiz)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)の項
								if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxr=1;
								}
								//Re(Aiy)の項
								if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyr=1;
								}
								//Re(Aiz)の項
								if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzr=1;
								}
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvi=1;
								}
							}
							if(flagxr==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxr;
							}
							if(flagyr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyr;
							}
							if(flagzr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzr;
							}
							
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvi;
							}
							///Bの計算
							//他の項と同じく、old_Aが必要ないのでBは更新されない
							//仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
					
				}
			}
		}
	}//////*/
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック //非対称行列なのでチェックしない
	//check_matrix_symmetry(pn,NUM,ROW,G);

	//対称性チェック
	/*////
	int countImIm=0;//Im行Im列
	int countImRe=0;
	int countReRe=0;
	int countReIm=0;
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						if(i<=pnI && ROW[i][j]<=pnI) countImIm++;
						if(i>pnI && ROW[i][j]>pnI) countReRe++;

						if(i<=pnI && ROW[i][j]>pnI) countImRe++;
						if(i>pnI && ROW[i][j]<=pnI) countReIm++;

						//cout<<"対称性ｴﾗｰ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
						//cout<<"対称性ｴﾗｰ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						if(ppn[i-1]>0 && ppn[J-1]==-3) cout<<"対称性ｴﾗｰ (AxI,AxR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-1 && ppn[J-1]==-4) cout<<"対称性ｴﾗｰ (AyI,AyR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-2 && ppn[J-1]==-5) cout<<"対称性ｴﾗｰ (AzI,AzR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-6 && ppn[J-1]==-7) cout<<"対称性ｴﾗｰ(φI,φR) ("<<i<<","<<J<<")="<<G[i][j]<<" ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else cout<<"対称性ｴﾗｰ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}
	if(countImIm>0) cout<<"虚数行,虚数列の項で非対称 num="<<countImIm<<endl;
	if(countImRe>0) cout<<"虚数行,実数列の項で非対称 num="<<countImRe<<endl;
	if(countReIm>0) cout<<"実数行,虚数列の項で非対称 num="<<countReIm<<endl;
	if(countReRe>0) cout<<"実数行,実数列の項で非対称 num="<<countReRe<<endl;
	////*/

	
	/*////第一象限のみの対称性を調べる
	cout<<"第一象限の対称性チェック"<<endl;
	double **G2=new double *[pnI+1];///全体行列
    for(int i=1;i<=pnI;i++) G2[i]=new double [width_mat[i]+1];
	int **ROW2=new int *[pnI+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pnI;i++) ROW2[i]=new int [width_mat[i]+1];
	int *NUM2=new int [pnI+1]; ///各行の、非ゼロ要素数

	 for(int i=1;i<=pnI;i++)//初期化
    {
        NUM2[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G2[i][j]=0;
			ROW2[i][j]=0;
		}
    }

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			if(ROW[i][j]>pnI)
			{
				NUM2[i]=NUM2[i]+1;
				ROW2[i][NUM2[i]]=ROW[i][j]-pnI;
				G2[i][NUM2[i]]=G[i][j];
			}
		}
	}

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM2[i];j++)
		{
			int J=ROW2[i][j];
			int flag=0;
			for(int k=1;k<=NUM2[J];k++)
			{
				if(ROW2[J][k]==i)
				{
					flag=1;
					if(G2[i][j]!=G2[J][k])
					{
						cout<<"対称性ｴﾗｰ ("<<i<<","<<J<<")="<<G2[i][j]<<"  ("<<J<<","<<ROW2[J][k]<<")="<<G2[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}

	for(int i=0;i<=pnI;i++) delete [] G2[i];
    delete [] G2;
	for(int i=0;i<=pnI;i++) delete [] ROW2[i];
    delete [] ROW2;
	delete [] NUM2;
	////////*/

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

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	//*///

	//行列の視覚化
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
  

	//CG法実行
	double *XX=new double [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);

	//XXに格納された解を各節点に振る
	//complex<double> *Ac[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ（複素数）
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) AI[A_X][i]=XX[n];
		else if(i==-1) AI[A_Y][ppn[n-1]]=XX[n];
		else if(i==-2) AI[A_Z][ppn[n-2]]=XX[n];
		else if(i==-3) AR[A_X][ppn[n-pnI]]=XX[n];
		else if(i==-4) AR[A_Y][ppn[n-pnI-1]]=XX[n];
		else if(i==-5) AR[A_Z][ppn[n-pnI-2]]=XX[n];
		else if(i==-6) VI[ppn[n-3]]=XX[n];//電位φ
		else if(i==-7) VR[ppn[n-3-pnI]]=XX[n];//電位φ
	}	

	delete [] XX;
	
	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
			phi=atan(AI[D][i]/AR[D][i]);
			A[D][i]=Am*cos(omega*t+phi);
		}
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*t+phi);

	}

	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///渦電流と強制電流を足した値を出力させる
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_Aの更新 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_Aを粒子に渡す
	cout<<"old_Aを粒子に記憶";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_Aのﾌｧｲﾙ出力
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] width_mat;	
	delete [] width_node;	
}

//////////
//////////
//jω法3,複素数行列を作成 //FEMをするたび、静的要素も含めて全領域でベクトルポテンシャルを計算する
void Avector3D_node_eddy2_jw3(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw)
{ 
	/////
	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始（ｊω法-3）"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//ディリクレ型境値

    ///ディリクレ型境界条件入力

	ifstream old("old_A.dat");
	if(!old) cout<<"cannot open old_A.dat"<<endl;
	old.unsetf(ifstream::dec);
	old.setf(ifstream::skipws);

	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old>>old_A[D][i];

	old.close();	
			
    for(int i=1;i<=node;i++)
    {
		if(CON->get_static_dirichlet()==OFF)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
		if(CON->get_static_dirichlet()==ON)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
    }//////////////
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである

	
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num;///未知数 複素数のまま格納するので時間差分と変わらない
    //int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
    //int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
	vector <int> ppn;
	vector <int> npp;
	npp.reserve(node+1);
    int num=0; 

    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)
		{
			//ppn[num]=i;
			//ppn[num+1]=-1;
			//ppn[num+2]=-2;
			ppn.push_back(i);
			ppn.push_back(-1);
			ppn.push_back(-2);
			npp[i]=num;
			num+=3;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }

	cout<<"num="<<num<<"pn="<<pn<<endl;
	////行列の幅計算 各行ごとに幅を求め、メモリを削減する
	int mat_w=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=4*mat_w;//
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////配列確保
	complex<double> **G=new complex<double> *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double>(0.0,0.0);
			ROW[i][j]=0;
		}
    }

	delete [] width_mat;	
	delete [] width_node;	

    complex<double> *B=new complex<double> [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//初期化
    ////
    
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];


		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=complex<double> (delta*j0x.real()/4,delta*j0x.imag()/4);
				B[I]+=complex<double> (delta*j0y.real()/4,delta*j0y.imag()/4);
				B[I+1]+=complex<double> (delta*j0z.real()/4,delta*j0z.imag()/4);
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aixの項
							    flag=1;
							}
						}
						if(flag==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aixの項
						    ROW[I][H]=J;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B[I-1]-=complex<double> (B_temp,0);
					}
				}

				/////Y方向
				I=npp[N[n]]+1+1;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiyの項
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B[I-1]-=complex<double> (B_temp,0);
					}
				}
				/////*/

				/////Z方向
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						J3=J2+1;		//Aizの項
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aizの項
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aizの項
						    ROW[I][H]=J3;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B[I-1]-=complex<double> (B_temp,0);
					}
				}
			}
		}
	}

	//////渦電流項計算
	int J4,flag4;
	//int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			//double co2=delta6*sig[je];//   σ/(36V)
			//double co3=sig[je]*dt*delta6;

			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=(delta/20.0)*omega*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			double co2=delta6*sig[je];//   σ/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J4=J+3;			//φの項
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aixの項
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag=1;
								}
							
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}///
							}
							if(flag==0)//相当するＪが存在しなかったら作る(ここではそんなはずないけど)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								G_temp=co;
								if(I==J) G[I][H]+=complex<double>(0,G_temp*2);
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Y方向

					I=npp[N[n]]+1+1;/////Y方向
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J4=J+3;			//φの項
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J2) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag2=1;
								}

								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co;
								if(I==J2) G[I][H]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J2;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J2) Sy+=2*Ay;
							//else Sy+=Ay;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sy;//支配方程式のX成分から得られるBの値

					/////Z方向
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J3) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co;
								if(I==J3) G[I][H]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
								else G[I][H]+=complex<double>(0,G_temp);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J3) Sz+=2*Az;
							//else Sz+=Az;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sz;//支配方程式のX成分から得られるBの値

					/////φ
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aixの項
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag=1;
								}
								//Aiyの項
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag2=1;
								}
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3は(1/j)がかかっている。1/j=-j
									G[I][h]+=complex<double>(0,-G_temp);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3は(1/j)がかかっている。1/j=-j
								G[I][H]+=complex<double>(0,-G_temp);
								ROW[I][H]=J4;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}///////
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

	cout<<"非零数="<<number<<endl;

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix_complex(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"行列の対称性チェック"<<endl;
	check_matrix_symmetry_complex(pn,NUM,ROW,G);

	 complex<double> *val = new complex<double> [number];
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
    /////////////////////

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*////////////
	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	//////*/
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG法実行
	complex<double> *XX=new complex<double> [pn];//行列の答え格納
    
	//if(CON->get_FEMCG()==1) ICCG3D2_complex(CON,val,ind,ptr,pn,B,number,XX);//
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	
	//XXに格納された解を各節点に振る
	//complex<double> *Ac[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ（複素数）
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	ofstream ar("old_AR.dat");
	ofstream ai("old_AI.dat");
	ofstream vr("old_VR.dat");
	ofstream vi("old_VI.dat");

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			AR[D][i]=0;
			AI[D][i]=0;
		}
		VR[i]=0;
		VI[i]=0;
	}

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[A_X][i]=XX[n].real();
			AI[A_X][i]=XX[n].imag();
		}
		else if(i==-1)
		{
			AR[A_Y][ppn[n-1]]=XX[n].real();
			AI[A_Y][ppn[n-1]]=XX[n].imag();
		}
		else if(i==-2)
		{
			AR[A_Z][ppn[n-2]]=XX[n].real();
			AI[A_Z][ppn[n-2]]=XX[n].imag();
		}
		else if(i==-3)
		{
			VR[ppn[n-3]]=XX[n].real();
			VI[ppn[n-3]]=XX[n].imag();
		}
	}	

	delete [] XX;
	
	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	///Am,φのﾌｧｲﾙ出力 φは電位ではなく位相の遅れ
	ofstream a("Am.dat");
	ofstream p("phi.dat");
	
	double Am[3];
	double phi[3];

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0.0;
			phi[D]=0.0;
		}
		
		if(NODE[i].boundary_condition==0)
		{
			for(int D=0;D<3;D++)
			{
				Am[D]=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(AI[D][i],AR[D][i]);
				//A[D][i]=Am[D];//波高を出力
				//if(CON->get_jw_Faverage()==ON) A[D][i]=Am[D]*sqrt(2.0);
				//else A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
				A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
			}
		}
		a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
		ar<<AR[A_X][i]<<" "<<AR[A_Y][i]<<" "<<AR[A_Z][i]<<endl;
		ai<<AI[A_X][i]<<" "<<AI[A_Y][i]<<" "<<AI[A_Z][i]<<endl;
		vr<<VR[i]<<endl;
		vi<<VI[i]<<endl;
	}
	a.close();
	p.close();
	ar.close();
	ai.close();
	vr.close();
	vi.close();

	//Vm,phi
	for(int i=1;i<=node;i++)
	{
		double Vm=0.0;
		double phi=0.0;
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}
			if(flagi==ON)
			{
				Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
				phi=atan2(VI[i],VR[i]);
				//phi=atan(VI[i]/VR[i]);
				V[i]=Vm*cos(omega*TIME+phi);
			}
		}
	}

	/*/////old_Aを粒子に渡す
	cout<<"old_Aを粒子に記憶";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	/*/////

	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//渦電流損[W]
	double *Je_loss_n=new double[node+1];//渦電流損[W]

	//for(int i=1;i<=node;i++) Je_loss_n[i]=0;

	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	calc_node_eddy_current_jw(CON,NODE,ELEM,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t ,TIME,sig,omega);//Jeには波の高さが入る

	//渦電流損を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//本当の体積

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//節点の渦電流損を対応する粒子へ渡す
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//

	/*//old_Aのﾌｧｲﾙ出力
	if(t==1)
	{
		cout<<"初期ステップのベクトルポテンシャル記憶"<<endl;
		ofstream g("old_A.dat");
		for(int i=1;i<=node;i++) g<<A[A_X][i]<<" "<<A[A_Y][i]<<" "<<A[A_Z][i]<<endl;
		g.close();
	}
	///*/
	
	
	//このステップで作成したNODE,ELEMをNODE_jw,ELEM_jwに記憶させる。ここで求めた波高と遅れを以降のステップで利用してBを求めるため
	NODE_jw.clear();
	ELEM_jw.clear();
	NODE_jw.resize(node+1);
	ELEM_jw.resize(nelm+1);
	for(int i=1;i<=node;i++) NODE_jw[i]=NODE[i];
	for(int i=1;i<=nelm;i++) ELEM_jw[i]=ELEM[i];

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    //delete [] npp;
    //delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	//////
}

//////////
//jω法4,行方向と列方向の格納順番を変える(没、うまくいかない）
void Avector3D_node_eddy2_jw4(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{   
	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始（ｊω法-4）"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON->get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=(3*(node-NN)+conducter_num)*2;///未知数 複素数なので各成分に実部と虚部がある
    int pnI=3*(node-NN)+conducter_num;//未知数のうち虚数部の数
	int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
    int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=node;i++)//行列はIm(A1x),Im(A1y),Im(A1z),Im(φ1),Im(A2x),・・・,Im(φn),Re(A1x),Re(A1y)の順に格納
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;//Imx
			ppn[num+1]=-1;//Imy
			ppn[num+2]=-2;//Imz
			ppn[num+pnI]=-3;//Rex
			ppn[num+1+pnI]=-4;//Rey
			ppn[num+2+pnI]=-5;//Rez
			npp[i]=num;
			num+=3;//4つめ以降は実部として分離されているので注意
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-6;//Imφ
						ppn[num+pnI]=-7;
						num+=1;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-6;//Imφ
						ppn[num+pnI]=-7;
						num+=1;
					}
			}
		}
		else npp[i]=pn+1;
    }

	/*///行列の最大幅計算  もっとも幅の大きいのはφの要素だろうと仮定している。確実ではないことに注意
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=8*mat_w;//X,Y,Z,φ成分にそれぞれRe,Imがあるので×8
    //////

	////配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
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
    //*/
    

	////行列の幅計算 各行ごとに幅を求め、メモリを削減する
	int mat_w=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=8*mat_w;//X,Y,Z,φ成分にそれぞれRe,Imがあるので×8

	double mat_ave=0.0;
	for(int i=1;i<=node;i++) mat_ave+=width_node[i];
	mat_ave=mat_ave/node;
	cout<<"mat_ave="<<mat_ave<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==2)
			{
			}

			if(flagi==ON)
			{
				width_node[i]*=8;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=4;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=6;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				for(int j=1;j<=3;j++) width_mat[j+count_mat+pnI]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////配列確保
	double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }

	//delete [] width_mat;	
	//delete [] width_node;	

    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    //*/
    
	cout<<"全体行列作成開始 pn="<<pn<<" pnI="<<pnI<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	//int J,J2,J3,flag,flag2,flag3;
	int Jxr,Jxi,Jyr,Jyi,Jzr,Jzi,flagxr,flagxi,flagyr,flagyi,flagzr,flagzi;
    for(int je=1;je<=nelm;je++)
    {   
		//if(je%10==0) cout<<je<<endl;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0);
			j0y=complex<double> (current[A_Y][je],0);
			j0z=complex<double> (current[A_Z][je],0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];
	
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向、実部

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=delta/4*j0x.real();//Rex
				B[I]+=delta/4*j0y.real();//Rey
				B[I+1]+=delta/4*j0z.real();//Rez
				B[I-1+pnI]+=delta/4*j0x.imag();//Imx
				B[I+pnI]+=delta/4*j0y.imag();//Imy
				B[I+1+pnI]+=delta/4*j0z.imag();//Imz
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
						flagxr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flagxr=1;
							}
						}
						if(flagxr==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=Jxr;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}


				/////Y方向,実部
				I=npp[N[n]]+2;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyr=npp[N[m]]+2+pnI;			//Re(Aiy)の項
						flagyr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flagyr=1;
							}
						}
						if(flagyr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyr;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				
				/////*/

				/////Z方向,実部
				I=npp[N[n]]+3;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzr=npp[N[m]]+3+pnI;	//Aixの項
						flagzr=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flagzr=1;
							}
						}
						if(flagzr==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzr;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////X方向、虚部

				 I=npp[N[n]]+1+pnI;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jxi=npp[N[m]]+1;	//Im(Aix)の項
						flagxi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flagxi=1;
							}
						}
						if(flagxi==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=Jxi;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y方向,虚部
				I=npp[N[n]]+2+pnI;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jyi=npp[N[m]]+2;			//Im(Aiy)の項
						flagyi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flagyi=1;
							}
						}
						if(flagyi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jyi;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Z方向,虚部
				I=npp[N[n]]+3+pnI;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						Jzi=npp[N[m]]+3;	//Aixの項
						flagzi=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flagzi=1;
							}
						}
						if(flagzi==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=Jzi;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}

	cout<<"ベクトルポテンシャル項完了"<<endl;

	//////渦電流項計算
	//int J4,flag4;
	int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==2)
		{
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			//double co2=delta6*sig[je];//   σ/(36V)
			//double co3=sig[je]*dt*delta6;

			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=(delta/20.0)*omega*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			double co2=delta6*sig[je];//   σ/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向,虚部

					int I=npp[N[n]]+1+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
							//Jxi=npp[N[m]]+2;	//Im(Aix)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)の項
								if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jxr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagxr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvi=1;
								}///
							}
							if(flagxr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jxr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Ｙ方向,虚部

					I=npp[N[n]]+2+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)の項
							//Jyi=npp[N[m]]+4;	//Im(Aiy)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiy)の項
								if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jyr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagyr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvi=1;
								}///
							}
							if(flagyr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jyr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;

					/////Ｚ方向,虚部

					I=npp[N[n]]+3+pnI;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)の項
							//Jzi=npp[N[m]]+6;	//Im(Aiz)の項
							//Jvr=npp[N[m]]+7;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aiz)の項
								if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jzr) G[I][h]+=co*2;//起こらないはず
									else G[I][h]+=co;//G=co(-Aui+jAur)
									flagzr=1;
								}
							
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvi=1;
								}///
							}
							if(flagzr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzr) G[I][H]+=co*2;//起こらないはず
								else G[I][H]+=co;
								ROW[I][H]=Jzr;
							}
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvi;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値//ｊω法だとBが増えない


					/////φ,虚部//co3がかかる項の虚実に注意
					I=npp[N[n]]+4+pnI;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
							Jxi=npp[N[m]]+1;	//Im(Aix)の項
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)の項
							Jyi=npp[N[m]]+2;	//Im(Aiy)の項
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)の項
							Jzi=npp[N[m]]+3;	//Im(Aiz)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項
							

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)の項
								if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxi=1;
								}
								//Aiyの項
								if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyi=1;
								}
								//Aizの項
								if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzi=1;
								}
								//Reφの項 G=α(φi-jφr)
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvr=1;
								}
							}
							if(flagxi==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxi;
							}
							if(flagyi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyi;
							}
							if(flagzi==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzi;
							}
							
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]-=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvr;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);

					/////X方向,実部

					 I=npp[N[n]]+1;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jxr=npp[N[m]]+1;	//Re(Aix)の項
							Jxi=npp[N[m]]+1;	//Im(Aix)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagxr=0;
							flagxi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aix)の項
								if(Jxi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jxi)
									{
										//cout<<"対角項に足されている？"<<endl;
										G[I][h]-=co*2;//起こらないはず
									}
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagxi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flagvr=1;
								}///
							}
							if(flagxi==0)
							{  
								
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jxi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jxi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					

					/////Y方向,実部

					I=npp[N[n]]+2;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jyr=npp[N[m]]+3;	//Re(Aiy)の項
							Jyi=npp[N[m]]+2;	//Im(Aiy)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagyr=0;
							flagyi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiy)の項
								if(Jyi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jyi) G[I][h]-=co*2;//起こらないはず
									else G[I][h]-=co;//G=co(-Aui+jAur)
									flagyi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flagvr=1;
								}///
							}
							if(flagyi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jyi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jyi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					//支配方程式のX成分から得られるBの値//ｊω法だとBが増えない

					/////Z方向,実部

					I=npp[N[n]]+3;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							//Jzr=npp[N[m]]+5;	//Re(Aiz)の項
							Jzi=npp[N[m]]+3;	//Im(Aiz)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							//Jvi=npp[N[m]]+8;	//Imφの項
							flagzr=0;
							flagzi=0;
							flagvr=0;
							flagvi=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Im(Aiz)の項
								if(Jzi==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==Jzi) G[I][h]-=co*2;//起こらないはず
									else G[I][h]-=co;//G=co(-Ai+jAr)
									flagzi=1;
								}
							
								//Reφの項
								if(Jvr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flagvr=1;
								}///
							}
							if(flagzi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==Jzi) G[I][H]-=co*2;//起こらないはず
								else G[I][H]-=co;
								ROW[I][H]=Jzi;
							}
							if(flagvr==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=Jvr;
							}

							//old_Aはｊω法の場合必要ないので、ここでBは変更されない

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					
					/////φ,実部 //co3がかかる項の虚実に注意
					
					I=npp[N[n]]+4;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							Jxr=npp[N[m]]+1+pnI;	//Re(Aix)の項
							Jxi=npp[N[m]]+1;	//Im(Aix)の項
							Jyr=npp[N[m]]+2+pnI;	//Re(Aiy)の項
							Jyi=npp[N[m]]+2;	//Im(Aiy)の項
							Jzr=npp[N[m]]+3+pnI;	//Re(Aiz)の項
							Jzi=npp[N[m]]+3;	//Im(Aiz)の項
							Jvr=npp[N[m]]+4+pnI;	//Reφの項
							Jvi=npp[N[m]]+4;	//Imφの項

							flagxr=0; flagxi=0; flagyr=0; flagyi=0; flagzr=0; flagzi=0; flagvr=0; flagvi=0;

							for(int h=1;h<=NUM[I];h++)
							{
								//Re(Aix)の項
								if(Jxr==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flagxr=1;
								}
								//Re(Aiy)の項
								if(Jyr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flagyr=1;
								}
								//Re(Aiz)の項
								if(Jzr==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flagzr=1;
								}
								//Imφの項
								if(Jvi==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flagvi=1;
								}
							}
							if(flagxr==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=Jxr;
							}
							if(flagyr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=Jyr;
							}
							if(flagzr==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=Jzr;
							}
							
							if(flagvi==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=Jvi;
							}
							///Bの計算
							//他の項と同じく、old_Aが必要ないのでBは更新されない
							//仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
					
				}
			}
		}
	}//////*/
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"対称性チェック"<<endl;
	check_matrix_symmetry(pn,NUM,ROW,G);

	//対角成分の値チェック
	cout<<"対角成分の値チェック"<<endl;
	int count_plus=0;
	int count_minus=0;
	int count_zero=0;
	for(int i=1;i<=pn;i++)
	{
		int flag=OFF;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			if(i==J)
			{
				if(G[i][j]==0)
				{
					count_zero++;
					if(ppn[i-1]==-6) cout<<"対角成分が零(φi） i="<<i<<endl;
					else if(ppn[i-1]==-7) cout<<"対角成分が零(φr） i="<<i<<endl;
					else cout<<"対角成分が零(A） i="<<i<<endl;
				}
				if(G[i][j]>0) count_plus++;
				if(G[i][j]<0) count_minus++;
				flag=ON;

			}
		}
		if(flag==OFF) count_zero++;
	}
	cout<<"正の対角項の数="<<count_plus<<endl;
	cout<<"負の対角項の数="<<count_minus<<endl;
	cout<<"零の対角項の数="<<count_zero<<endl;
	///

	//対称性チェック
	/*////
	int countImIm=0;//Im行Im列
	int countImRe=0;
	int countReRe=0;
	int countReIm=0;
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						if(i<=pnI && ROW[i][j]<=pnI) countImIm++;
						if(i>pnI && ROW[i][j]>pnI) countReRe++;

						if(i<=pnI && ROW[i][j]>pnI) countImRe++;
						if(i>pnI && ROW[i][j]<=pnI) countReIm++;

						//cout<<"対称性ｴﾗｰ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
						//cout<<"対称性ｴﾗｰ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						if(ppn[i-1]>0 && ppn[J-1]==-3) cout<<"対称性ｴﾗｰ (AxI,AxR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-1 && ppn[J-1]==-4) cout<<"対称性ｴﾗｰ (AyI,AyR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-2 && ppn[J-1]==-5) cout<<"対称性ｴﾗｰ (AzI,AzR) ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else if(ppn[i-1]==-6 && ppn[J-1]==-7) cout<<"対称性ｴﾗｰ(φI,φR) ("<<i<<","<<J<<")="<<G[i][j]<<" ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
						else cout<<"対称性ｴﾗｰ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}
	if(countImIm>0) cout<<"虚数行,虚数列の項で非対称 num="<<countImIm<<endl;
	if(countImRe>0) cout<<"虚数行,実数列の項で非対称 num="<<countImRe<<endl;
	if(countReIm>0) cout<<"実数行,虚数列の項で非対称 num="<<countReIm<<endl;
	if(countReRe>0) cout<<"実数行,実数列の項で非対称 num="<<countReRe<<endl;
	////*/

	
	/*////第一象限のみの対称性を調べる
	cout<<"第一象限の対称性チェック"<<endl;
	double **G2=new double *[pnI+1];///全体行列
    for(int i=1;i<=pnI;i++) G2[i]=new double [width_mat[i]+1];
	int **ROW2=new int *[pnI+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pnI;i++) ROW2[i]=new int [width_mat[i]+1];
	int *NUM2=new int [pnI+1]; ///各行の、非ゼロ要素数

	 for(int i=1;i<=pnI;i++)//初期化
    {
        NUM2[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G2[i][j]=0;
			ROW2[i][j]=0;
		}
    }

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			if(ROW[i][j]>pnI)
			{
				NUM2[i]=NUM2[i]+1;
				ROW2[i][NUM2[i]]=ROW[i][j]-pnI;
				G2[i][NUM2[i]]=G[i][j];
			}
		}
	}

	for(int i=1;i<=pnI;i++)
	{
		for(int j=1;j<=NUM2[i];j++)
		{
			int J=ROW2[i][j];
			int flag=0;
			for(int k=1;k<=NUM2[J];k++)
			{
				if(ROW2[J][k]==i)
				{
					flag=1;
					if(G2[i][j]!=G2[J][k])
					{
						cout<<"対称性ｴﾗｰ ("<<i<<","<<J<<")="<<G2[i][j]<<"  ("<<J<<","<<ROW2[J][k]<<")="<<G2[J][k]<<endl;
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
	}

	for(int i=1;i<=pnI;i++) delete [] G2[i];
    delete [] G2;
	for(int i=1;i<=pnI;i++) delete [] ROW2[i];
    delete [] ROW2;
	delete [] NUM2;
	////////*/

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

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	//*///

	//行列の視覚化
	ofstream m("matrix.dat");
	int mat_min=1;
	int mat_max=20000;
	for(int i=1;i<=pn;i++)
	{
		 for(int j=1;j<=NUM[i];j++)
		{
			double matX=ROW[i][j];
			double matY=-i;
			if(matX<mat_max && matY>-1*mat_max && matX>mat_min && matY<-1*mat_min) m<<matX<<" "<<matY<<endl;
		}
	}
	m.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
  

	//CG法実行
	double *XX=new double [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);

	//XXに格納された解を各節点に振る
	//complex<double> *Ac[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ（複素数）
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0) AI[A_X][i]=XX[n];
		else if(i==-1) AI[A_Y][ppn[n-1]]=XX[n];
		else if(i==-2) AI[A_Z][ppn[n-2]]=XX[n];
		else if(i==-3) AR[A_X][ppn[n-pnI]]=XX[n];
		else if(i==-4) AR[A_Y][ppn[n-pnI-1]]=XX[n];
		else if(i==-5) AR[A_Z][ppn[n-pnI-2]]=XX[n];
		else if(i==-6) VI[ppn[n-3]]=XX[n];//電位φ
		else if(i==-7) VR[ppn[n-3-pnI]]=XX[n];//電位φ
	}	

	delete [] XX;
	
	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
			phi=atan(AI[D][i]/AR[D][i]);
			A[D][i]=Am*cos(omega*TIME+phi);
		}
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*TIME+phi);

	}

	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/*///渦電流と強制電流を足した値を出力させる
	for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Je[D][i]+=current[D][i];
	check_J0Je(CON, node, nelm,NODE,ELEM,Je);
	//*/

	//old_Aの更新 
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///

	//old_Aを粒子に渡す
	cout<<"old_Aを粒子に記憶";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	
	///old_Aのﾌｧｲﾙ出力
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	/*if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] width_mat;	
	delete [] width_node;	
}

//jω法5,複素数行列を作成 //２ステップ目以降、静的要素のベクトルポテンシャル値は１ステップ目の値をディリクレ値として利用する
void Avector3D_node_eddy2_jw5(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw)
{ 
	/////
	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始（ｊω法-5 ）"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界節点数
	int NN_V=0;//ディリクレ型境界節点数のうち、渦電流を考える節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	
	complex<double> *PHAT_A[3];
	for(int D=0;D<3;D++) PHAT_A[D]=new complex<double> [node+1];//Aのディリクレ型境値
	complex<double> *PHAT_V = new complex<double> [node+1];//Vのディリクレ境界値
	

    ///ディリクレ型境界条件入力
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	//初期化
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			AR[D][i]=0;
			AI[D][i]=0;
		}
		VR[i]=0;
		VI[i]=0;
	}

	/////////1ステップ目の解読み込み
	cout<<"静的要素のディリクレ値読み込み"<<endl;
	//AR
	ifstream ar("old_AR.dat");
	if(!ar) cout<<"cannot open old_AR.dat"<<endl;
	ar.unsetf(ifstream::dec);
	ar.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ar>>AR[D][i];
	ar.close();

	//AI
	ifstream ai("old_AI.dat");
	if(!ai) cout<<"cannot open old_AI.dat"<<endl;
	ai.unsetf(ifstream::dec);
	ai.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ai>>AI[D][i];
	ai.close();

	//VR
	ifstream vr("old_VR.dat");
	if(!vr) cout<<"cannot open old_VR.dat"<<endl;
	vr.unsetf(ifstream::dec);
	vr.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vr>>VR[i];
	vr.close();

	//VI
	ifstream vi("old_VI.dat");
	if(!vi) cout<<"cannot open old_VI.dat"<<endl;
	vi.unsetf(ifstream::dec);
	vi.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vi>>VI[i];
	vi.close();
	/////
	cout<<"ディリクレ値読み込み完了"<<endl;		
	int count_r=0;
    for(int i=1;i<=node;i++)
    {
		/*///この関数にはstatic_dirichlet=OFFのとき飛んで来ない
		if(CON->get_static_dirichlet()==OFF)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
		////*/
		if(CON->get_static_dirichlet()==ON)
		{
			/*///
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			///*/
			////else
			
			{
				if(NODE[i].remesh==OFF)//初期ステップ以外はリメッシュしない節点について最初に求めたベクトルポテンシャルをディリクレ値として与える
				{
					count_r++;
					NODE[i].boundary_condition=3;//静的要素の境界条件をもとの0から3に書き換える。固定境界(1,2)も変えるが問題ないはず
					dn[i]=NN;
					PHAT_A[A_X][NN]=complex<double>(AR[A_X][i],AI[A_X][i]);
					PHAT_A[A_Y][NN]=complex<double>(AR[A_Y][i],AI[A_Y][i]);
					PHAT_A[A_Z][NN]=complex<double>(AR[A_Z][i],AI[A_Z][i]);
					PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					NN++;

					
					if(CON->get_Je_crucible()==0)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!流体が静的節点"<<endl;
					}
					if(CON->get_Je_crucible()==1)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!流体が静的節点"<<endl;
						if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) NN_V++;
					}
					if(CON->get_J0eqJe()==1) NN_V=0;
				
				}
				else
				{
					dn[i]=node+1;
				}
			}
		}
    }//////////////
    cout<<"ﾃﾞｨﾘｸﾚ節点数＝"<<NN<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN+NN_V<<" NN_V="<<NN_V<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
	cout<<"non_remesh節点数="<<count_r<<endl;

	
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			//if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==FLUID &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==FLUID) conducter_num++;
			//if(NODE[i].material==CRUCIBLE) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num+100;///未知数 複素数のまま格納するので時間差分と変わらない //なぜかここで止まるので少し余裕を取って配列を取ってみる
	cout<<"pn="<<pn<<endl;
    int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
    int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 

	cout<<"pn="<<pn<<" ppn,nppの割り振り-----";
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			npp[i]=num;
			num+=3;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						ppn[num]=-3;
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						ppn[num]=-3;
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }
	cout<<"ok"<<endl;

	////行列の幅計算 各行ごとに幅を求め、メモリを削減する
	int mat_w=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=4*mat_w;//

	cout<<"calc_matrix終了"<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		
		if(NODE[i].boundary_condition==0)
		{
			cout<<i<<endl;
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID)
				{
					flagi=ON;
					
				}
				if(NODE[i].material==CRUCIBLE)
				{
					cout<<"error! るつぼがboundary"<<endl;				
					flagi=ON;
				}
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    //////

	////配列確保
	cout<<"G,ROW,NUMの確保開始"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double>(0.0,0.0);
			ROW[i][j]=0;
		}
    }

	delete [] width_mat;	
	delete [] width_node;	

    complex<double> *B=new complex<double> [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//初期化
    ////
    
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_dir;//解行列に入るディリクレ値
	complex<double> B_N;//上にかかる形状関数由来の値
	B_dir=complex<double> (0.0,0.0);
	B_N=complex<double> (0.0,0.0);

    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		//if(ELEM[je].material==COIL)
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];


		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=complex<double> (delta*j0x.real()/4,delta*j0x.imag()/4);
				B[I]+=complex<double> (delta*j0y.real()/4,delta*j0y.imag()/4);
				B[I+1]+=complex<double> (delta*j0z.real()/4,delta*j0z.imag()/4);
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aixの項
							    flag=1;
							}
						}
						if(flag==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aixの項
						    ROW[I][H]=J;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double>(B_temp,0); 
					    B[I-1]-=PHAT_A[A_X][NN]*B_N;
					}
				}

				/////Y方向
				I=npp[N[n]]+1+1;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiyの項
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B_N=complex<double> (B_temp,0);
						B[I-1]-=PHAT_A[A_Y][NN]*B_N;
					}
				}
				/////*/

				/////Z方向
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						J3=J2+1;		//Aizの項
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aizの項
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aizの項
						    ROW[I][H]=J3;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double> (B_temp,0);
					    B[I-1]-=PHAT_A[A_Z][NN]*B_N;
					}
				}
			}
		}
	}

	//////渦電流項計算
	int J4,flag4;
	//int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			//double co2=delta6*sig[je];//   σ/(36V)
			//double co3=sig[je]*dt*delta6;

			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=(delta/20.0)*omega*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			double co2=delta6*sig[je];//   σ/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J4=J+3;			//φの項
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aixの項
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag=1;
								}
							
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}///
							}
							if(flag==0)//相当するＪが存在しなかったら作る(ここではそんなはずないけど)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								G_temp=co;
								if(I==J) G[I][H]+=complex<double>(0,G_temp*2);
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co;
							if(I==J) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//φ
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Y方向

					I=npp[N[n]]+1+1;/////Y方向
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J4=J+3;			//φの項
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J2) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag2=1;
								}

								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co;
								if(I==J2) G[I][H]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J2;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J2) Sy+=2*Ay;
							//else Sy+=Ay;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Ay
							B_temp=co;
							if(I==J2) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//φ
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sy;//支配方程式のX成分から得られるBの値

					/////Z方向
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J3) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co;
								if(I==J3) G[I][H]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
								else G[I][H]+=complex<double>(0,G_temp);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J3) Sz+=2*Az;
							//else Sz+=Az;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Az
							B_temp=co;
							if(I==J3) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//φ
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
						}
					}////////*/
					//B[I-1]+=co*Sz;//支配方程式のX成分から得られるBの値

					/////φ
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aixの項
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag=1;
								}
								//Aiyの項
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag2=1;
								}
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3は(1/j)がかかっている。1/j=-j
									G[I][h]+=complex<double>(0,-G_temp);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3は(1/j)がかかっている。1/j=-j
								G[I][H]+=complex<double>(0,-G_temp);
								ROW[I][H]=J4;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//Ay
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//Az
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//φ
							B_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
							B_N=complex<double>(0,-B_temp);
							B[I-1]-=PHAT_V[NN]*B_N;

							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}///////
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

	cout<<"非零数="<<number<<endl;

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix_complex(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"行列の対称性チェック"<<endl;
	check_matrix_symmetry_complex(pn,NUM,ROW,G);

	 complex<double> *val = new complex<double> [number];
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
    /////////////////////

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*////////////
	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	//////*/
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG法実行
	complex<double> *XX=new complex<double> [pn];//行列の答え格納
    
	//if(CON->get_FEMCG()==1) ICCG3D2_complex(CON,val,ind,ptr,pn,B,number,XX);//
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);

	
	//XXに格納された解を各節点に振る
	//complex<double> *Ac[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ（複素数）
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];

	/*////
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	/////*/
	
	cout<<"複素数解の割り振り----";
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[A_X][i]=XX[n].real();
			AI[A_X][i]=XX[n].imag();
		}
		else if(i==-1)
		{
			AR[A_Y][ppn[n-1]]=XX[n].real();
			AI[A_Y][ppn[n-1]]=XX[n].imag();
		}
		else if(i==-2)
		{
			AR[A_Z][ppn[n-2]]=XX[n].real();
			AI[A_Z][ppn[n-2]]=XX[n].imag();
		}
		else if(i==-3)
		{
			VR[ppn[n-3]]=XX[n].real();
			VI[ppn[n-3]]=XX[n].imag();
		}
	}	

	delete [] XX;

	cout<<"ok"<<endl;
	
	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	///Am,φのﾌｧｲﾙ出力 φは電位ではなく位相の遅れ
	ofstream a("Am.dat");
	ofstream p("phi.dat");
	
	double Am[3];
	double phi[3];

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0.0;
			phi[D]=0.0;
		}
		
		//if(NODE[i].boundary_condition==0)
		{
			for(int D=0;D<3;D++)
			{
				Am[D]=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(AI[D][i],AR[D][i]);
				//A[D][i]=Am[D]/sqrt(2.0);//実効値で計算
				A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
			}
		}
		a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	a.close();
	p.close();

	//Vm,phi
	for(int i=1;i<=node;i++)
	{
		double Vm=0.0;
		double phi=0.0;
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}
			if(flagi==ON)
			{
				Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
				phi=atan2(VI[i],VR[i]);
				//phi=atan(VI[i]/VR[i]);
				V[i]=Vm*cos(omega*TIME+phi);
			}
		}
	}
	

	/*/////old_Aを粒子に渡す
	cout<<"old_Aを粒子に記憶";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	/*/////

	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//渦電流損[W]
	double *Je_loss_n=new double[node+1];//渦電流損[W]

	//for(int i=1;i<=node;i++) Je_loss_n[i]=0;

	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	calc_node_eddy_current_jw(CON,NODE,ELEM,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t ,TIME,sig,omega);//Jeには波の高さが入る

	//渦電流損を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//本当の体積

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//節点の渦電流損を対応する粒子へ渡す
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//

	///old_Aのﾌｧｲﾙ出力
	if(t==1)
	{
		cout<<"初期ステップのベクトルポテンシャル記憶"<<endl;
		ofstream g("old_A.dat");
		for(int i=1;i<=node;i++) g<<A[A_X][i]<<" "<<A[A_Y][i]<<" "<<A[A_Z][i]<<endl;
		g.close();
	}
	
	
	/////このステップで作成したNODE,ELEMをNODE_jw,ELEM_jwに記憶させる。ここで求めた波高と遅れを以降のステップで利用してBを求めるため
	
	NODE_jw.clear();
	ELEM_jw.clear();
	/*////
	NODE_jw.resize(node+1);
	ELEM_jw.resize(nelm+1);
	for(int i=1;i<=node;i++) NODE_jw[i]=NODE[i];
	for(int i=1;i<=nelm;i++) ELEM_jw[i]=ELEM[i];
	//////*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT_A[D];
	delete [] PHAT_V;
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	//////
}

//jω法5,複素数行列を作成 //２ステップ目以降、静的要素のベクトルポテンシャル値は１ステップ目の値をディリクレ値として利用する//ver2では未知数に関係する配列をvector化
void Avector3D_node_eddy2_jw5_ver2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **A,int *jnb,int **nei,double **old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta)
{ 
	/////
	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始（外部領域境界化ver2 ）"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界節点数
	int NN_V=0;//ディリクレ型境界節点数のうち、渦電流を考える節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	
	complex<double> *PHAT_A[3];
	for(int D=0;D<3;D++) PHAT_A[D]=new complex<double> [node+1];//Aのディリクレ型境値
	complex<double> *PHAT_V = new complex<double> [node+1];//Vのディリクレ境界値
	

    ///ディリクレ型境界条件入力
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	//初期化
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			AR[D][i]=0;
			AI[D][i]=0;
		}
		VR[i]=0;
		VI[i]=0;
	}

	/////////1ステップ目の解読み込み
	cout<<"静的要素のディリクレ値読み込み"<<endl;
	//AR
	ifstream ar("old_AR.dat");
	if(!ar) cout<<"cannot open old_AR.dat"<<endl;
	ar.unsetf(ifstream::dec);
	ar.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ar>>AR[D][i];
	ar.close();

	//AI
	ifstream ai("old_AI.dat");
	if(!ai) cout<<"cannot open old_AI.dat"<<endl;
	ai.unsetf(ifstream::dec);
	ai.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) ai>>AI[D][i];
	ai.close();

	//VR
	ifstream vr("old_VR.dat");
	if(!vr) cout<<"cannot open old_VR.dat"<<endl;
	vr.unsetf(ifstream::dec);
	vr.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vr>>VR[i];
	vr.close();

	//VI
	ifstream vi("old_VI.dat");
	if(!vi) cout<<"cannot open old_VI.dat"<<endl;
	vi.unsetf(ifstream::dec);
	vi.setf(ifstream::skipws);
	for(int i=1;i<=node;i++) vi>>VI[i];
	vi.close();
	/////
	cout<<"ディリクレ値読み込み完了"<<endl;		
	int count_r=0;
    for(int i=1;i<=node;i++)
    {
		/*///この関数にはstatic_dirichlet=OFFのとき飛んで来ない
		if(CON->get_static_dirichlet()==OFF)
		{
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT[A_X][NN]=0;
				PHAT[A_Y][NN]=0;
				PHAT[A_Z][NN]=0;
				A[A_X][i]=0;
				A[A_Y][i]=0;
				A[A_Z][i]=0;
				NN++;
			}
			else
			{
				dn[i]=node+1;
			}
		}
		////*/
		if(CON->get_static_dirichlet()==ON)
		{
			/*///
			if(NODE[i].boundary_condition==2)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			else if(NODE[i].boundary_condition==1)
			{    
				dn[i]=NN;
				PHAT_A[A_X][NN]=complex<double>(0,0);
				PHAT_A[A_Y][NN]=complex<double>(0,0);
				PHAT_A[A_Z][NN]=complex<double>(0,0);
				PHAT_V[NN]=complex<double>(0,0);
				NN++;
			}
			///*/
			////else
			
			{
				//if(NODE[i].remesh==OFF)//初期ステップ以外はリメッシュしない節点について最初に求めたベクトルポテンシャルをディリクレ値として与える
				if(i<=node_sta)
				{
					count_r++;
					NODE[i].boundary_condition=3;//静的要素の境界条件をもとの0から3に書き換える。固定境界(1,2)も変えるが問題ないはず
					dn[i]=NN;
					PHAT_A[A_X][NN]=complex<double>(AR[A_X][i],AI[A_X][i]);
					PHAT_A[A_Y][NN]=complex<double>(AR[A_Y][i],AI[A_Y][i]);
					PHAT_A[A_Z][NN]=complex<double>(AR[A_Z][i],AI[A_Z][i]);
					PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					NN++;

					
					if(CON->get_Je_crucible()==0)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!流体が静的節点"<<endl;
					}
					if(CON->get_Je_crucible()==1)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!流体が静的節点"<<endl;
						if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) NN_V++;
					}
					if(CON->get_J0eqJe()==1) NN_V=0;
				
				}
				else
				{
					dn[i]=node+1;
				}
			}
		}
    }//////////////
    cout<<"ﾃﾞｨﾘｸﾚ節点数＝"<<NN<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN+NN_V<<" NN_V="<<NN_V<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
	cout<<"non_remesh節点数="<<count_r<<endl;

	
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0) if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
		if(CON->get_Je_crucible()==1)
		{
			//if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==FLUID &&jnb[i]!=0) conducter_num++;
			//if(NODE[i].material==CRUCIBLE &&jnb[i]!=0) conducter_num++;
			if(NODE[i].material==FLUID) conducter_num++;
			//if(NODE[i].material==CRUCIBLE) conducter_num++;
		}
	}
	if(CON->get_J0eqJe()==1) conducter_num=0;
	//////////////

	int pn=3*(node-NN)+conducter_num;///未知数 複素数のまま格納するので時間差分と変わらない //数え上げ失敗してる？//要素形成に失敗した節点の数×3だけずれてる
	cout<<"pn="<<pn<<endl;
    //int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
	//int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
	
	vector<int> ppn;

	vector<int> npp;
	npp.reserve(node+1); 
    
    int num=0; 

	cout<<"pn(参考値)="<<pn<<" ppn,nppの割り振り-----";
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)
		{
			//ppn[num]=i;
			//ppn[num+1]=-1;
			//ppn[num+2]=-2;
			ppn.push_back(i);
			ppn.push_back(-1);
			ppn.push_back(-2);
			npp[i]=num;
			num+=3;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}

			if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE)
					{
						//ppn[num]=-3;
						ppn.push_back(-3);
						num++;
					}
			}
		}
		else npp[i]=pn+1;
    }
	cout<<"ok"<<"num="<<num<<endl;

	pn=num;//numが真のpnのはず

	///////
	////行列の幅計算 各行ごとに幅を求め、メモリを削減する
	int mat_w=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅
	mat_w=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//メモリをきっちりもとめる。だけど若干の計算コストになる。
	mat_w=4*mat_w;//

	cout<<"calc_matrix終了"<<endl;
	
	int count_mat=0;
	for(int i=1;i<=node;i++)
	{
		
		if(NODE[i].boundary_condition==0)
		{
			
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID)
				{
					flagi=ON;
					
				}
				if(NODE[i].material==CRUCIBLE)
				{
					cout<<"error! るつぼがboundary"<<endl;				
					flagi=ON;
				}
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}

			if(flagi==ON)
			{
				width_node[i]*=4;
				for(int j=1;j<=4;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=4;
			}
			else if(flagi==OFF)
			{
				width_node[i]*=3;
				for(int j=1;j<=3;j++) width_mat[j+count_mat]=width_node[i];
				count_mat+=3;
			}
		}
	}	
    /////*/

	////配列確保
	cout<<"G,ROW,NUMの確保開始"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double>(0.0,0.0);
			ROW[i][j]=0;
		}
    }

	delete [] width_mat;	
	delete [] width_node;	

    complex<double> *B=new complex<double> [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//初期化
    ////
    
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_dir;//解行列に入るディリクレ値
	complex<double> B_N;//上にかかる形状関数由来の値
	B_dir=complex<double> (0.0,0.0);
	B_N=complex<double> (0.0,0.0);

    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		//if(ELEM[je].material==COIL)
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];


		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=complex<double> (delta*j0x.real()/4,delta*j0x.imag()/4);
				B[I]+=complex<double> (delta*j0y.real()/4,delta*j0y.imag()/4);
				B[I+1]+=complex<double> (delta*j0z.real()/4,delta*j0z.imag()/4);
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aixの項
							    flag=1;
							}
						}
						if(flag==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aixの項
						    ROW[I][H]=J;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double>(B_temp,0); 
					    B[I-1]-=PHAT_A[A_X][NN]*B_N;
					}
				}

				/////Y方向
				I=npp[N[n]]+1+1;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aiyの項
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					    B_N=complex<double> (B_temp,0);
						B[I-1]-=PHAT_A[A_Y][NN]*B_N;
					}
				}
				/////*/

				/////Z方向
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						J3=J2+1;		//Aizの項
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);//Aizの項
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);//Aizの項
						    ROW[I][H]=J3;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
						B_temp=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						B_N=complex<double> (B_temp,0);
					    B[I-1]-=PHAT_A[A_Z][NN]*B_N;
					}
				}
			}
		}
	}

	//////渦電流項計算
	int J4,flag4;
	//int Jvr,Jvi,flagvr,flagvi;
	
	for(int je=1;je<=nelm;je++)
    {
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			//double co=delta/20.0/dt*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			//double co2=delta6*sig[je];//   σ/(36V)
			//double co3=sig[je]*dt*delta6;

			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=(delta/20.0)*omega*sig[je];//σV/(20dt)　何度も計算することになるので係数化する	
			double co2=delta6*sig[je];//   σ/(36V)
			double co3=sig[je]*delta6/omega;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J4=J+3;			//φの項
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aixの項
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag=1;
								}
							
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}///
							}
							if(flag==0)//相当するＪが存在しなかったら作る(ここではそんなはずないけど)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								G_temp=co;
								if(I==J) G[I][H]+=complex<double>(0,G_temp*2);
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J) Sx+=2*Ax;
							//else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co;
							if(I==J) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//φ
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					//B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Y方向

					I=npp[N[n]]+1+1;/////Y方向
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J4=J+3;			//φの項
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J2) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag2=1;
								}

								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co;
								if(I==J2) G[I][H]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
								else G[I][H]+=complex<double>(0,G_temp);
								ROW[I][H]=J2;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J2) Sy+=2*Ay;
							//else Sy+=Ay;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Ay
							B_temp=co;
							if(I==J2) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//φ
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					//B[I-1]+=co*Sy;//支配方程式のX成分から得られるBの値

					/////Z方向
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co;
									if(I==J3) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									G[I][h]+=complex<double>(G_temp,0);
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co;
								if(I==J3) G[I][H]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
								else G[I][H]+=complex<double>(0,G_temp);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								G[I][H]+=complex<double>(G_temp,0);
								ROW[I][H]=J4;
							}//*/
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//if(I==J3) Sz+=2*Az;
							//else Sz+=Az;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Az
							B_temp=co;
							if(I==J3) B_N=complex<double> (0,B_temp*2);
							else B_N=complex<double> (0,B_temp); 
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//φ
							B_temp=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
							B_N=complex<double> (B_temp,0);
							B[I-1]-=PHAT_V[NN]*B_N;
						}
					}////////*/
					//B[I-1]+=co*Sz;//支配方程式のX成分から得られるBの値

					/////φ
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aixの項
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag=1;
								}
								//Aiyの項
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag2=1;
								}
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									G[I][h]+=complex<double>(G_temp,0);
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3は(1/j)がかかっている。1/j=-j
									G[I][h]+=complex<double>(0,-G_temp);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
								G[I][H]+=complex<double>(G_temp,0);
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);//co3は(1/j)がかかっている。1/j=-j
								G[I][H]+=complex<double>(0,-G_temp);
								ROW[I][H]=J4;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//Ax
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_X][NN]*B_N;
							//Ay
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Y][NN]*B_N;
							//Az
							B_temp=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							B_N=complex<double>(B_temp,0);
							B[I-1]-=PHAT_A[A_Z][NN]*B_N;
							//φ
							B_temp=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
							B_N=complex<double>(0,-B_temp);
							B[I-1]-=PHAT_V[NN]*B_N;

							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					//B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}///////
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

	cout<<"非零数="<<number<<endl;

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix_complex(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"行列の対称性チェック"<<endl;
	check_matrix_symmetry_complex(pn,NUM,ROW,G);

	 complex<double> *val = new complex<double> [number];
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
    /////////////////////

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*////////////
	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	//////*/
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG法実行
	complex<double> *XX=new complex<double> [pn];//行列の答え格納
    
	//if(CON->get_FEMCG()==1) ICCG3D2_complex(CON,val,ind,ptr,pn,B,number,XX);//
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);
	
	//XXに格納された解を各節点に振る
	//complex<double> *Ac[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ（複素数）
	//for(int D=0;D<3;D++) Ac[D]=new complex<double> [node+1];

	/*////
	double *AR[3];//ベクトルポテンシャル実部
	for(int D=0;D<3;D++) AR[D]=new double [node+1];
	double *AI[3];//ベクトルポテンシャル虚部
	for(int D=0;D<3;D++) AI[D]=new double [node+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	/////*/
	
	cout<<"複素数解の割り振り----";
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[A_X][i]=XX[n].real();
			AI[A_X][i]=XX[n].imag();
		}
		else if(i==-1)
		{
			AR[A_Y][ppn[n-1]]=XX[n].real();
			AI[A_Y][ppn[n-1]]=XX[n].imag();
		}
		else if(i==-2)
		{
			AR[A_Z][ppn[n-2]]=XX[n].real();
			AI[A_Z][ppn[n-2]]=XX[n].imag();
		}
		else if(i==-3)
		{
			VR[ppn[n-3]]=XX[n].real();
			VI[ppn[n-3]]=XX[n].imag();
		}
	}	

	delete [] XX;

	cout<<"ok"<<endl;
	
	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	///Am,φのﾌｧｲﾙ出力 φは電位ではなく位相の遅れ
	ofstream a("Am.dat");
	ofstream p("phi.dat");
	
	double Am[3];
	double phi[3];

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0.0;
			phi[D]=0.0;
		}
		
		//if(NODE[i].boundary_condition==0)
		{
			for(int D=0;D<3;D++)
			{
				Am[D]=sqrt(AR[D][i]*AR[D][i]+AI[D][i]*AI[D][i]);
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(AI[D][i],AR[D][i]);
				//A[D][i]=Am[D]/sqrt(2.0);//実効値で計算
				A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
			}
		}
		a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	a.close();
	p.close();

	//Vm,phi
	for(int i=1;i<=node;i++)
	{
		double Vm=0.0;
		double phi=0.0;
		if(NODE[i].boundary_condition==0)
		{
			int flagi=OFF;
			if(CON->get_Je_crucible()==0)
			{
				if(NODE[i].material==FLUID) flagi=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(NODE[i].material==FLUID || NODE[i].material==CRUCIBLE) flagi=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagi=OFF;
			}
			if(flagi==ON)
			{
				Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
				phi=atan2(VI[i],VR[i]);
				//phi=atan(VI[i]/VR[i]);
				V[i]=Vm*cos(omega*TIME+phi);
			}
		}
	}
	

	/*/////old_Aを粒子に渡す
	cout<<"old_Aを粒子に記憶";
	for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	cout<<" ok"<<endl;
	/*/////

	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//渦電流損[W]
	double *Je_loss_n=new double[node+1];//渦電流損[W]
	int *count_e=new int[node+1];//各説点まわりの渦電流がけいさんされた要素の数が格納。

	//for(int i=1;i<=node;i++) Je_loss_n[i]=0;

	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	calc_node_eddy_current_jw(CON,NODE,ELEM,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t ,TIME,sig,omega);//Jeには波の高さが入る

	//渦電流損を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++)
			{
				N[j]=ELEM[je].node[j];
				if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
				{
					Je_loss_n[N[j]]+=Je_loss_e[je];
					count_e[N[j]]=count_e[N[j]]+1;
				}
			}
		}
	}
	for(int i=1;i<=node;i++)
	{
		if(count_e[i]!=0)
		{
			Je_loss_n[i]/=count_e[i];
		}
	}
	
	//節点の渦電流損を対応する粒子へ渡す
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	delete [] count_e;
	//

	///old_Aのﾌｧｲﾙ出力
	if(t==1)
	{
		cout<<"初期ステップのベクトルポテンシャル記憶"<<endl;
		ofstream g("old_A.dat");
		for(int i=1;i<=node;i++) g<<A[A_X][i]<<" "<<A[A_Y][i]<<" "<<A[A_Z][i]<<endl;
		g.close();
	}
	
	
	/////このステップで作成したNODE,ELEMをNODE_jw,ELEM_jwに記憶させる。ここで求めた波高と遅れを以降のステップで利用してBを求めるため
	
	NODE_jw.clear();
	ELEM_jw.clear();
	/*////
	NODE_jw.resize(node+1);
	ELEM_jw.resize(nelm+1);
	for(int i=1;i<=node;i++) NODE_jw[i]=NODE[i];
	for(int i=1;i<=nelm;i++) ELEM_jw[i]=ELEM[i];
	//////*/

	delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT_A[D];
	delete [] PHAT_V;
	for(int D=0;D<3;D++) delete [] AR[D];
	for(int D=0;D<3;D++) delete [] AI[D];
	delete [] VR;
	delete [] VI;

    //delete [] npp;
    //delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	//////
}

///ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算関数(辺要素用)
void Avector3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP)
{
	cout<<"ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始---";

	double u0=4*PI*0.0000001;	//真空の透磁率
    double v0=1/u0;				//磁気抵抗率
	double j0x,j0y,j0z;			//電流密度[A/m^3]
	unsigned timeA=GetTickCount();

	//磁石の着磁方向を決定
	double MA=0;//=CON->get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON->get_magnet_B()*magnet_direction[A_X];
	double My=CON->get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON->get_magnet_B()*magnet_direction[A_Z];
	

	//////////////////////////////////////////////////////*/

	///周囲壁に固定境界条件を設定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].edge[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//節点pを含まない辺は境界辺
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//ディリクレ型境界辺数
    int *dn=new int [side_num+1]; //各辺がディリクレ型境界の何番目か。ディリクレ型境界上でないならside_num+1を格納
    double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"ﾃﾞｨﾘｸﾚ数="<<NN;
	/////////////*/
    
	    
    int pn=side_num-NN;				///未知数
    int *ppn=new int [pn];			//行列のn番目は辺ppn[n]
    int *npp=new int [side_num+1];	///各辺が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//未知数
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
    ////行列の最大幅計算
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]を含む要素数
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].edge[j];
			nume[edge]=nume[edge]+1;
		}
	}											///nume[i]がもとまった
	for(int i=1;i<=side_num;i++)
	{	//考え方:ある辺周りの要素群を考える。まず、その内の一つを空間においた段階で、辺の数は4つ。以後、要素を加えるごとに3つ辺が増える。よって次式となる
		int width=4+3*(nume[i]-1);//場合によってはこれより小さい値になるかもしれない。けど多くメモリを確保するぶんには問題ない
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

	
    ////配列確保
    double **G=new double *[pn+1];///全体行列
	for(int i=1;i<=pn;i++) G[i]=new double [wid[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
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
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//cout<<"行列作成開始 ";
    for(int je=1;je<=nelm;je++)
    {   
		//辺－節点ﾃｰﾌﾞﾙ作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//要素の重心座標
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
		
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
	
		
		double rp=RP[je];
	//	if(ELEM[je].material==FLUID) rp=CON->get_RP();//比透磁率
		////要素ﾏﾄﾘｸｽ作成開始
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//要素jeの辺番号
			if(SIDE[iside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				int I1=SIDE[iside].node[1];//isideを構成する2点
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jsideを構成する2点
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}
				///B[I-1]を計算する
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    ///行列の実際の最大幅を求める
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	

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

	delete [] wid;
    
	cout<<" 行列作成 幅："<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

	///////////////////////行列計算開始
	double *XX=new double [pn];//行列の答え格納
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

void Avector3D_edge_eddy(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{
	//ELEM.edge:要素を構成する辺６つ　EDGE.node:辺をつくる点２つ　
	/////
	int flageddy=OFF;
	cout<<"渦電流を考慮にいれた辺要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始"<<endl;

	double u0=PI*4e-7;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	double j0x,j0y,j0z;					//電流密度
	double Sx=0;
	double Sy=0;
	double Sz=0;
	//complex<double> j0x;
	//complex<double> j0y;
	//complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界辺数
    int *dn=new int [nedge+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnedge+1を格納
	double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
	
	///解析領域の境界に固定境界条件を設定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//節点pを含まない辺は境界辺
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"辺要素の固定境界設定完了"<<endl;

    ///ディリクレ型境界条件入力
    for(int i=1;i<=nedge;i++)
    {
        if(EDGE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else
		{
			dn[i]=nedge+1;
		}
    }//////////////
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<endl;
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える
	//辺要素でもφについては節点要素を用いる

	int conducter_num=0;//導体節点数
	int *conducter_flag=new int [node+1];//導体節点ならON
	for(int i=0;i<=node;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
	}
	if(CON->get_Je_crucible()==-1) conducter_num=0;
	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	int pn_e=nedge-NN;///辺要素の未知数
	int pn;
	if(flageddy==ON) pn=pn_e+conducter_num;//φも含めた全未知数
	if(flageddy==OFF)pn=pn_e;

	int *ppn=new int [pn];		///行列でn番目は辺番号ppn[n]
    int *npp=new int [nedge+node+1];	///各辺,導体節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
   
	for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//未知数
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
	
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node;i++)//φ
		{
			if(conducter_flag[i]==ON)//導体節点
			{
				ppn[num]=nedge+i;//
				npp[nedge+i]=num;//節点要素のnppは、辺要素のnppをすべて格納したあとの配列を用いる
				num++;
			}
			else npp[nedge+i]=pn+1;
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////行列の幅計算 辺要素
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]を含む要素数
	int *width_edge= new int [pn+1];

	for(int i=1;i<=pn;i++) width_edge[i]=0;

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]がもとまった
	for(int i=1;i<=nedge;i++)
	{
		int width=7+3*(nume[i]-2);
		if(width>mat_we) mat_we=width;
		width_edge[npp[i]+1]=width;	
	}

	
	cout<<"A-Aの行列幅計算終了"<<endl;

	////行列の幅計算 節点要素
	int mat_wn=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅

	for(int i=0;i<=pn;i++) width_mat[i]=0;
	mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_nodeが求まる
	
	//各行の節点要素由来の非零数を求める width_matに格納
	
	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	cout<<"nume_max="<<nume_max<<endl;
	
	if(flageddy==ON)
	{
		//A-φ
		int **ROW2=new int *[pn+1];
		for(int i=1;i<=pn;i++) ROW2[i]=new int [nume_max*5+1];
		for(int i=1;i<=pn;i++)//初期化
		{
			for(int j=1;j<=nume_max*5;j++)
			{
				ROW2[i][j]=0;
			}
		}
		int N2[4+1]; //要素の各節点番号格納	
		
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//要素jeの辺番号
					if(EDGE[iside].boundary_condition==0)///未知数
					{   
						int I=npp[iside]+1;///辺isideは行列のI番目
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N2[j]]==ON)
							{
								int J=npp[N2[j]+nedge]+1;
								int flag=0;			
			
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
								}
							}
						}
					}
				}
			}
		}
		cout<<"A-φの行列幅計算終了"<<endl;
		
		//φ-A
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				
				for(int i=1;i<=4;i++)
				{
					
					
					if(conducter_flag[N2[i]]==ON)
					{
						int I=npp[N2[i]+nedge]+1;
						if(I<=pn_e) cout<<"I<=pn_e I="<<I<<endl;
						for(int j=1;j<=6;j++)
						{	
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								int J=npp[iside]+1;///辺isideは行列のI番目
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
									//B[I-1]の計算

								}
							}		
						}
					}
				}
			}
		}

		cout<<"A-φ,φ-Aの行列幅計算終了"<<endl;
		for(int i=1;i<=pn;i++) delete [] ROW2[i];
		delete [] ROW2;


		//φ-φ
		for(int i=1;i<=node;i++)
		{
			if(conducter_flag[i]==ON)//導体節点
			{
				width_mat[npp[i+nedge]+1]+=width_node[i];
			}
		}
		cout<<"φ-φの行列幅計算終了"<<endl;
	}
	
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	delete [] nume;
	

	////配列確保
	//cout<<"全体行列用配列宣言"<<endl;
	double **G=new double *[pn+1];///全体行列
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数


    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=0.0;
			ROW[i][j]=0;
		}
    }

		
	delete [] width_node;
	delete [] width_edge;
	delete [] width_mat;

    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0.0;//初期化
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//int J,J2,J3,flag,flag2,flag3;
	//double G_temp=0.0;
	//double B_temp=0.0;

	//静磁場項
    for(int je=1;je<=nelm;je++)
    {
		//if(je%1000==0) cout<<je<<endl;
		//辺－節点ﾃｰﾌﾞﾙ作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//重心座標
		double Ys=0;
		double Zs=0;

		double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//積和　x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j]*0.25;
			Ys+=Y[j]*0.25;
			Zs+=Z[j]*0.25;
			
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
    
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=current[A_X][je];
			j0y=current[A_Y][je];
			j0z=current[A_Z][je];
		}
		else
		{
			j0x=0;
			j0y=0;
			j0z=0;
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////要素ﾏﾄﾘｸｽ作成開始
		//A-A
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//要素jeの辺番号
			if(EDGE[iside].boundary_condition==0)///未知数
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				//int I1=EDGE[iside].node[1];//isideを構成する2点
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jsideを構成する2点
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///未知数
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6*v;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6*v;
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}
				///B[I-1]を計算する
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B[I-1]+=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}

	//渦電流項
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//渦電流項の計算対象となる要素
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;

			double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//積和　x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//∂A/∂t
		
			double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			
			//cout<<"導体要素 "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//要素jeの辺番号
				if(EDGE[iside].boundary_condition==0)///未知数
				{   
					int I=npp[iside]+1;///辺isideは行列のI番目
					//int I1=EDGE[iside].node[1];//isideを構成する2点
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jsideを構成する2点
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///未知数
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									
									//∫∫∫zdV=V*Zsだが、∫∫∫z*zdVはV*Zs*Zsではない。(x,yも同様) 教科書p84
									//(Nk)x・(Nu)x

									///////////
									G[I][h]+=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G[I][h]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G[I][h]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G[I][h]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G[I][h]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									G[I][h]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									//(Nk)y・(Nu)y
									//////
									///////
									G[I][h]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G[I][h]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G[I][h]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G[I][h]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G[I][h]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G[I][h]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////
									////////
									//(Nk)z・(Nu)z
									G[I][h]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G[I][h]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G[I][h]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G[I][h]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G[I][h]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G[I][h]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									//////

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"起こらないはず"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y・(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z・(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								//////
							    
								ROW[I][H]=J;
							}
							///B[I-1]を計算する
							//仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int n=dn[jside];
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}///////////
					}///////
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-φ
				double co2=sig[je]*delta*delta6*delta6*delta6;//σV*(1/6V)^3
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//要素jeの辺番号
					if(EDGE[iside].boundary_condition==0)///未知数
					{   
						int I=npp[iside]+1;///辺isideは行列のI番目
						int I1=EDGE[iside].node[1];//isideを構成する2点
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G[I][h]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
										G[I][h]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
										G[I][h]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G[I][H]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
									G[I][H]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
									G[I][H]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
									ROW[I][H]=J;
								}
							}
						}
					}
				}
				//cout<<"A-φ eddy"<<endl;

				//φ-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"φ-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int J=npp[iside]+1;///辺isideは行列のI番目
							//cout<<"φ-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//isideを構成する2点
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G[I][h]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
										G[I][h]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
										G[I][h]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									G[I][H]+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
									G[I][H]+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
									G[I][H]+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
									ROW[I][H]=J;
								}
								//B[I-1]の計算
							}
							else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int n=dn[iside];
								B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}///////////
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"φ-A eddy"<<endl;

				//φ-φ
				double co3=sig[je]*dt*delta*delta6*delta6;
				for(int i=1;i<=4;i++)
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G[I][h]+=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										flag=1;
									}
								}
								if(flag==0)//相当するＪが存在しなかったら作る
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G[I][H]+=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									ROW[I][H]=J;
								}
							}
							else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int NN=dn[N[j]];
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}///////
					}
				}
				//cout<<"φ-φ eddy"<<endl;
			}
		}
	}

	///////
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"非ゼロ要素数"<<number<<endl;

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"行列の対称性チェック----";
	check_matrix_symmetry(pn,NUM,ROW,G);
	cout<<"完了"<<endl;;
	//解行列の値チェック
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
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
    /////////////////////

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/////////////
	///バンド幅を求める
	int *band = new int [pn+1];
	for(int i=0;i<=pn;i++) band[i]=0;
	int band_max=0;
	int m1=0;
	int m2=0;
	for(int i=1;i<=pn;i++)
	{
		//m1=abs(i-ROW[i][1]);//下三角側のバンド幅
		m2=abs(ROW[i][NUM[i]]-i);//上三角側のバンド幅
		//band[i]=m1+m2+1;
		band[i]=m2+1;
		if(band[i]>=band_max) band_max=band[i];
	}

	int bandmaxID=0;
	//unsigned long int band_sum=0;//大きすぎてlong intでも足りない
	long double band_sum=0;
	for(int i=1;i<=pn;i++)
	{
		band_sum+=band[i];
		if(band[i]==band_max) bandmaxID=i;
	}
	cout<<"最大バンド幅="<<band_max<<" i="<<bandmaxID<<endl;
	cout<<"合計バンド幅="<<band_sum<<endl;
	delete [] band;
	///////
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG法実行
	double *XX=new double [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) BiCGStab2_method(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==4) COCG(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==5) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);

	//XXに格納された解を各辺に振る

	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			A[i]=XX[n];
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			V[i]=XX[n];
		}
	}	

	delete [] XX;
	
	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	//calc_eddy_current(CON,NODE,ELEM,node,nelm,A,old_A,dt,V,Je,t,sig);
	for(int D=0;D<3;D++) delete [] Je[D];
	//

	/////
	//old_Aの更新 
	for(int i=1;i<=node;i++)
	{
		old_A[i]=A[i];
	}///

	//old_Aを粒子に渡す
	//cout<<"old_Aを粒子に記憶";
	//for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) PART[NODE[i].particleID].old_A[D]=old_A[D][i];	
	//cout<<" ok"<<endl;
	
	///old_Aのﾌｧｲﾙ出力
	ofstream g("old_A.dat");
	for(int i=1;i<=node;i++) g<<old_A[i]<<endl;
	g.close();
	
	//if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
    delete [] PHAT;

    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	
	//////
	/////*/
}

//辺要素による渦電流を考慮した磁気ベクトルポテンシャル計算関数,jw法
void Avector3D_edge_eddy_jw(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{
	//ELEM.edge:要素を構成する辺６つ　EDGE.node:辺をつくる点２つ　
	/////
	int flageddy=CON->get_A_phi();
	cout<<"渦電流を考慮にいれた辺要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始(jw)"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	double Sx=0;
	double Sy=0;
	double Sz=0;
	complex<double> Im=(0,1);
	
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界辺数
    int *dn=new int [nedge+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnedge+1を格納
	double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
	
	///解析領域の境界に固定境界条件を設定

	////
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//節点pを含まない辺は境界辺
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"辺要素の固定境界設定完了"<<endl;


	///ディリクレ型境界条件入力
	set_boundary_condition3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,dn,&NN,PHAT,A);
    
	/*//ディリクレ型境界条件入力
    for(int i=1;i<=nedge;i++)
    {
        if(EDGE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else
		{
			dn[i]=nedge+1;
		}
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<endl;
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える
	//辺要素でもφについては節点要素を用いる

	int conducter_num=0;//導体節点数
	int *conducter_flag=new int [node+1];//導体節点ならON
	for(int i=0;i<=node;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
	}
	if(CON->get_Je_crucible()==-1) conducter_num=0;
	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	int pn_e=nedge-NN;///辺要素の未知数 複素数のまま格納するので時間差分の２倍にはならない
	int pn=pn_e;//φも含めた全未知数
	if(flageddy==ON) pn=pn_e+conducter_num;//φも含めた全未知数
	//int pn=pn_e;

	//int *ppn=new int [pn];		///行列でn番目は辺番号ppn[n] 行列でn(ただし、nは辺数より大きい)番目は節点番号ppn[n](導体のみ格納)
    //int *npp=new int [nedge+node+1];	///各辺,導体節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    vector <int> ppn;
	vector <int> npp;
	npp.push_back(0);//0番目配列を埋める

	//npp.reserve(nedge+node+1);
	int num=0; 
   
	//for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//未知数
		{
			//ppn[num]=i;
			ppn.push_back(i);
			//npp[i]=num;
			npp.push_back(num);
			num++;
		}
		else npp.push_back(pn+1);   //npp[i]=pn+1;
    }
	cout<<"pn_e="<<pn_e<<" num="<<num<<endl;
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node;i++)//φ
		{
			if(conducter_flag[i]==ON)//導体節点
			{
				//ppn[num]=nedge+i;//
				ppn.push_back(i);//行列でn(ただし、nは辺数より大きい)番目は節点番号ppn[n]
				//npp[nedge+i]=num;//節点要素のnppは、辺要素のnppをすべて格納したあとの配列を用いる
				npp.push_back(num);
				num++;
			}
			else npp.push_back(pn+1); //npp[nedge+i]=pn+1;
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////行列の幅計算 辺要素
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]を含む要素数
	int *width_edge= new int [pn+1];

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]がもとまった
	for(int i=1;i<=nedge;i++)
	{
		int width=7+3*(nume[i]-2);
		if(width>mat_we) mat_we=width;
		if(npp[i]<pn_e) width_edge[npp[i]+1]=width;	
		//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
	}

	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	//cout<<"nume_max="<<nume_max<<endl;

	
	//cout<<"A-Aの行列幅計算終了"<<endl;

	//節点要素関係の行列幅変数の宣言
	int mat_wn=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *N_E_num=new int [node+1]; /////NODE[i]を含む要素数
	int *width_mat=new int [pn+1]; ///各行の幅



	for(int i=1;i<=node;i++) 
	{
		width_node[i]=0;
		N_E_num[i]=0;
	}

	for(int i=0;i<=pn;i++) 
	{
		width_mat[i]=0;
	}

	if(flageddy==ON)
	{
		////行列の幅計算
		
		for(int i=1;i<=nedge;i++) nume[i]=0;
		//各行の節点要素由来の非零数を求める width_matに格納

		//A-φ
		//int width_ap=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=6;j++)
				{
					int side=ELEM[je].edge[j];
					nume[side]=nume[side]+1;//
				}
			}
		}
		for(int i=1;i<=nedge;i++)
		{
			int width=0;//行の未知数幅
			if(nume[i]>0) width=3+nume[i];//ここでnume[]が0でないということは、その辺は導体要素に含まれている。要素が1つ追加されるたび、その辺に影響する導体節点は最大1つ増える
			if(npp[i]<pn_e) width_mat[npp[i]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"A-φの行列幅計算終了"<<endl;

		//φ-A
		//int width_pa=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=4;j++)
				{
					int N=ELEM[je].node[j];
					N_E_num[N]=N_E_num[N]+1;//
				}
			}
		}
		for(int i=1;i<=node;i++)
		{
			int width=0;//行の未知数幅
			if(N_E_num[i]>0) width=1+3*N_E_num[i];//ここでnume[]が0でないということは、その辺は導体要素に含まれている。要素が1つ追加されるたび、その節点に影響する導体辺は最大3つ増える
			if(npp[i+nedge]>=pn_e && npp[i+nedge]<pn) width_mat[npp[i+nedge]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"φ-Aの行列幅計算終了"<<endl;
		

		//φ-φ

		mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_nodeが求まる	

		for(int i=1;i<=node;i++)
		{
			if(conducter_flag[i]==ON)//導体節点
			{
				width_mat[npp[i+nedge]+1]=width_mat[npp[i+nedge]+1]+width_node[i];
			}
		}
		//cout<<"φ-φの行列幅計算終了"<<endl;
		/*////
		//A-φ
		int **ROW2=new int *[pn+1];
		for(int i=1;i<=pn;i++) ROW2[i]=new int [nume_max*4+1];
		for(int i=1;i<=pn;i++)//初期化
		{
			for(int j=1;j<=nume_max*4;j++)
			{
				ROW2[i][j]=0;
			}
		}
		int N2[4+1]; //要素の各節点番号格納	
		
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//要素jeの辺番号
					if(EDGE[iside].boundary_condition==0)///未知数
					{   
						int I=npp[iside]+1;///辺isideは行列のI番目
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N2[j]]==ON)
							{
								int J=npp[N2[j]+nedge]+1;
								int flag=0;			
			
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
								}
							}
						}
					}
				}
			}
		}
		cout<<"A-φの行列幅計算終了"<<endl;
		
		//φ-A
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				
				for(int i=1;i<=4;i++)
				{
					
					
					if(conducter_flag[N2[i]]==ON)
					{
						int I=npp[N2[i]+nedge]+1;
						if(I<=pn_e) cout<<"I<=pn_e I="<<I<<endl;
						for(int j=1;j<=6;j++)
						{	
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								int J=npp[iside]+1;///辺isideは行列のI番目
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
									//B[I-1]の計算

								}
							}		
						}
					}
				}
			}
		}



		cout<<"φ-Aの行列幅計算終了"<<endl;

		

		
		for(int i=1;i<=pn;i++) delete [] ROW2[i];
		delete [] ROW2;
		cout<<"1"<<endl;
		*/////

		
	}
	
	//行列幅の統合
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	
	for(int i=1;i<=pn;i++) 
	{
		if(width_mat[i]<=0)
		{
			cout<<"未知数幅が0以下の行あり i="<<i<<" "<<width_mat[i]<<endl;
			width_mat[i]=100;
		}
	}
	delete [] nume;
	
	////配列確保
	cout<<"全体行列用配列宣言"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///全体行列
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
	int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶

    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double> (0.0,0.0);
			ROW[i][j]=0;
		}
    }
	
	delete [] width_node;
	delete [] N_E_num;
	delete [] width_edge;
	delete [] width_mat;

    complex<double> *B=new complex<double> [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//初期化
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_N;//ディリクレ値として解行列に加算される項の形状関数部
	B_N=complex<double> (0.0,0.0);

	//静磁場項
    for(int je=1;je<=nelm;je++)
    {
		//if(je%1000==0) cout<<je<<endl;
		//辺－節点ﾃｰﾌﾞﾙ作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//重心座標
		double Ys=0;
		double Zs=0;

		double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//積和　x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
		Xs/=4;Ys/=4;Zs/=4;

		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////要素ﾏﾄﾘｸｽ作成開始
		//A-A(静磁場)
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//要素jeの辺番号
			if(EDGE[iside].boundary_condition==0)///未知数
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				//int I1=EDGE[iside].node[1];//isideを構成する2点
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jsideを構成する2点
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///未知数
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int n=dn[jside];
						B_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
						B[I-1]-=complex<double> (B_temp,0);
					}//////////*/
				}
				///強制電流項に関するB[I-1]を計算する
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B_temp=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x.real()*delta6/6;
					B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y.real()*delta6/6;
					B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z.real()*delta6/6;
					B[I-1]+=complex<double> (B_temp,0);
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}

	//渦電流項
	//if(flageddy==ON)
	{
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//渦電流項の計算対象となる要素
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;

			double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//積和　x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//∂A/∂t
		
			//double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=sig[je]*delta*delta6*delta6*delta6*delta6*omega;
			
			//cout<<"導体要素 "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//要素jeの辺番号
				if(EDGE[iside].boundary_condition==0)///未知数
				{   
					int I=npp[iside]+1;///辺isideは行列のI番目
					//int I1=EDGE[iside].node[1];//isideを構成する2点
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jsideを構成する2点
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///未知数
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									
									//∫∫∫zdV=V*Zsだが、∫∫∫z*zdVはV*Zs*Zsではない。(x,yも同様) 教科書p84
									//(Nk)x・(Nu)x

									///////////
									G_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
									G_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
									//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									///////////

									//(Nk)y・(Nu)y
									G_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////

									//(Nk)z・(Nu)z
									G_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									/////*/

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"起こらないはず"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y・(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z・(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								/////*/
							    
								ROW[I][H]=J;
							}
							///B[I-1]を計算する
							//仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//jw法においては、前のステップの値に基づく項が存在しないため、ここは必要ない

							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int n=dn[jside];

							//∫∫∫zdV=V*Zsだが、∫∫∫z*zdVはV*Zs*Zsではない。(x,yも同様) 教科書p84
									//(Nk)x・(Nu)x

							///////////
							B_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
							B_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
							//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							///////////

							//(Nk)y・(Nu)y
							B_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
							//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							////////

							//(Nk)z・(Nu)z
							B_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
							//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
							if(I==npp[jside]+1)
							{	cout<<"ディリクレ値が未知数のはずの行番号に存在"<<endl;
								B_N=complex<double>(0,PHAT[n]*B_temp*2);//coにはjがかかっているので虚数項へ追加//起こらない？
							}
							else B_N=complex<double>(0,PHAT[n]*B_temp);
							//B[I-1]-=PHAT[n]*B_N;//この記述は大丈夫？
							B[I-1]-=B_N;

									/////*/
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}//////////*/
					}//////*/
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-φ
				double co2=sig[je]*delta*delta6*delta6*delta6;//σV*(1/6V)^3
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//要素jeの辺番号
					if(EDGE[iside].boundary_condition==0)///未知数
					{   
						int I=npp[iside]+1;///辺isideは行列のI番目
						int I1=EDGE[iside].node[1];//isideを構成する2点
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//導体でないということはφ=0なので、ここでAの項のようなディリクレ値に関して考える必要はない
						}
					}
				}
				//cout<<"A-φ eddy"<<endl;

				//φ-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"φ-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int J=npp[iside]+1;///辺isideは行列のI番目
							//cout<<"φ-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//isideを構成する2点
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]の計算
							}
							else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int n=dn[iside];
								B_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								B_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								B_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHATは複素ディリクレに備えなおすべき
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"φ-A eddy"<<endl;

				//φ-φ
				//double co3=sig[je]*dt*delta*delta6*delta6;
				double co3=sig[je]*delta*delta6*delta6/omega;//これに1/jがかかっているため、複素数としては符号を変えて虚数項に足されることになる(1/j=-j)
				for(int i=1;i<=4;i++)
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//相当するＪが存在しなかったら作る
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int NN=dn[N[j]];
								B_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //ふつうの境界条件ではφが値を持たないので、ひとまず削除
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}//////*/
					}
				}
				//cout<<"φ-φ eddy"<<endl;
			}
		}
	}
	}

	///////
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"非ゼロ要素数"<<number<<endl;

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix_complex(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"行列の対称性チェック----";
	check_matrix_symmetry_complex(pn,NUM,ROW,G);
	//解行列の値チェック
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
	complex<double> *val = new complex<double> [number];
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
    /////////////////////

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG法実行
	complex<double> *XX=new complex<double> [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	//XXに格納された解を各辺に振る
	//cout<<"複素数解の割り振り"<<endl;

	double *AR = new double [nedge+1];
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int i=1;i<=nedge;i++)
	{
		
		AR[i]=0;
		AI[i]=0;
	}
	for(int i=1;i<=node;i++)
	{
		VR[i]=0;
		VI[i]=0;
	}

	ofstream ar("old_AR_e.dat");
	ofstream ai("old_AI_e.dat");
	ofstream vr("old_VR_e.dat");
	ofstream vi("old_VI_e.dat");
	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[i]=XX[n].real();
			AI[i]=XX[n].imag();
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			VR[i]=XX[n].real();
			VI[i]=XX[n].imag();
		}
	}	
	
	delete [] XX;

	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	///Am,φのﾌｧｲﾙ出力 φは電位ではなく位相の遅れ
	ofstream a("Am_e.dat");
	ofstream p("phi_e.dat");
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=nedge;i++)
	{
		if(EDGE[i].boundary_condition==0)
		{
			Am=sqrt(AR[i]*AR[i]+AI[i]*AI[i]);
			phi=atan2(AI[i],AR[i]);
			
			//if(CON->get_jw_Faverage()==ON) A[i]=Am/sqrt(2.0);
			//else A[i]=Am*cos(omega*t+phi);
			A[i]=Am*cos(omega*TIME+phi);
		}
		a<<Am<<endl;
		p<<phi<<endl;
		ar<<AR[i]<<endl;
		ai<<AI[i]<<endl;
	}
	for(int i=1;i<=node;i++)
	{
		vr<<VR[i]<<endl;
		vi<<VI[i]<<endl;
	}

	a.close();
	p.close();
	ar.close();
	ai.close();
	vr.close();
	vi.close();
	//cout<<"複素数を実数変換完了"<<endl;
	/*///
	for(int i=1;i<=node;i++)
	{
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*t+phi);
	}
	*/
	
	
	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//渦電流損[W]
	double *Je_loss_n=new double[node+1];//渦電流損[W]

	calc_edge_eddy_current_jw(CON,NODE,ELEM,EDGE,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t,TIME,sig,omega);
	
	
	//渦電流損を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//本当の体積

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//節点の渦電流損を対応する粒子へ渡す
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//
	
	//if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
    delete [] PHAT;

    delete [] B;

	delete [] AR;
	delete [] AI;
	delete [] VR;
	delete [] VI;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	
	//////
	/////*/
}

//辺要素による渦電流を考慮した磁気ベクトルポテンシャル計算関数,jw法,節点要素であるφだけ2次要素化
void Avector3D_edge_eddy_jw_with_parabolic_node_element(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega)
{
	//ELEM.edge:要素を構成する辺６つ　EDGE.node:辺をつくる点２つ　
	/////
	int flageddy=CON->get_A_phi();
	cout<<"渦電流を考慮にいれた1次辺2次節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始(jw)"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	double Sx=0;
	double Sy=0;
	double Sz=0;
	complex<double> Im=(0,1);
	
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界辺数
    int *dn=new int [nedge+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnedge+1を格納
	double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
	
	///解析領域の境界に固定境界条件を設定

	////
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//節点pを含まない辺は境界辺
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"辺要素の固定境界設定完了"<<endl;


	///ディリクレ型境界条件入力
	set_boundary_condition3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,dn,&NN,PHAT,A);
    
	/*//ディリクレ型境界条件入力
    for(int i=1;i<=nedge;i++)
    {
        if(EDGE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
			PHAT[NN]=0;
			A[i]=0.0;
	        NN++;
		}
		else
		{
			dn[i]=nedge+1;
		}
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<endl;
	
	//2次要素のため、導体辺を数える
	int *conduct_edge_flag=new int [nedge+1];//導体辺ならON
	for(int i=0;i<=nedge;i++) conduct_edge_flag[i]=OFF;

	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=6;j++)
			{
				int side=ELEM[je].edge[j];
				conduct_edge_flag[side]=ON;//導体要素に含まれる辺をチェックする
			}
		}
	}
	int para_node_num=0;
	
	for(int i=1;i<=nedge;i++)
	{
		EDGE[i].para_node_num=0;//初期化
		if(conduct_edge_flag[i]==ON)
		{
			para_node_num++;
			EDGE[i].para_node_num=node+para_node_num;
		}
	}

	cout<<"導体辺数="<<para_node_num<<endl;

	int node2=node+para_node_num;//2次節点も含めた全節点数

	
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える
	//辺要素でもφについては節点要素を用いる

	int conducter_num=0;//導体節点数 もともとのFEMメッシュで存在していた節点のみ数えている
	int *conducter_flag=new int [node2+1];//導体節点ならON
	for(int i=1;i<=node2;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
		else if(CON->get_Je_crucible()==-1)
		{
		}
	}
	for(int i=node+1;i<=node2;i++)//2次要素で追加される節点
	{
		conducter_flag[i]=ON;
	}


	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	cout<<"conducter_num2="<<conducter_num+para_node_num<<endl;
	int pn_e=nedge-NN;///辺要素の未知数 複素数のまま格納するので時間差分の２倍にはならない
	int pn=pn_e;//φも含めた全未知数
	if(flageddy==ON) pn=pn_e+conducter_num+para_node_num;//φも含めた全未知数
	//int pn=pn_e;

	//int *ppn=new int [pn];		///行列でn番目は辺番号ppn[n] 行列でn(ただし、nは辺数より大きい)番目は節点番号ppn[n](導体のみ格納)
    //int *npp=new int [nedge+node+1];	///各辺,導体節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    vector <int> ppn;
	vector <int> npp;
	npp.push_back(0);//0番目配列を埋める

	//npp.reserve(nedge+node+1);
	int num=0; 
   
	//for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//未知数
		{
			//ppn[num]=i;
			ppn.push_back(i);
			//npp[i]=num;
			npp.push_back(num);
			num++;
		}
		else npp.push_back(pn+1);   //npp[i]=pn+1;
    }
	cout<<"pn_e="<<pn_e<<" num="<<num<<endl;
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node2;i++)//φ
		{
			if(conducter_flag[i]==ON)//導体節点
			{
				//ppn[num]=nedge+i;//
				ppn.push_back(i);//行列でn(ただし、nは辺数より大きい)番目は節点番号ppn[n]
				//npp[nedge+i]=num;//節点要素のnppは、辺要素のnppをすべて格納したあとの配列を用いる
				npp.push_back(num);
				num++;
			}
			else npp.push_back(pn+1); //npp[nedge+i]=pn+1;
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////行列の幅計算 辺要素
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]を含む要素数
	int *width_edge= new int [pn+1];

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]がもとまった
	for(int i=1;i<=nedge;i++)
	{
		int width=7+3*(nume[i]-2);
		if(width>mat_we) mat_we=width;
		if(npp[i]<pn_e) width_edge[npp[i]+1]=width;	
		//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
	}

	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	//cout<<"nume_max="<<nume_max<<endl;

	
	//cout<<"A-Aの行列幅計算終了"<<endl;

	//節点要素関係の行列幅変数の宣言
	int mat_wn=0;
	int *width_node=new int [node2+1]; ///各節点の隣り合う節点の数＋１
	int *N_E_num=new int [node2+1]; /////NODE[i]を含む要素数
	int *width_mat=new int [pn+1]; ///各行の幅



	for(int i=1;i<=node2;i++) 
	{
		width_node[i]=0;
		N_E_num[i]=0;
	}

	for(int i=1;i<=pn;i++) 
	{
		width_mat[i]=0;
	}

	if(flageddy==ON)
	{
		////行列の幅計算
		
	
		//各行の節点要素由来の非零数を求める width_matに格納

		//A-φ
		//int width_ap=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=6;j++)
				{
					int side=ELEM[je].edge[j];
					nume[side]=nume[side]+1;//
				}
			}
		}
		for(int i=1;i<=nedge;i++)
		{
			int width=0;//行の未知数幅
			//if(nume[i]>0) width=3+nume[i];//ここでnume[]が0でないということは、その辺は導体要素に含まれている。要素が1つ追加されるたび、その辺に影響する導体節点は最大1つ増える
			if(nume[i]>0) width=6+4*nume[i];//ここでnume[]が0でないということは、その辺は導体要素に含まれている。要素が1つ追加されるたび、その辺に影響する導体節点は2次節点も含め最大4つ増える
			if(npp[i]<pn_e) width_mat[npp[i]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"A-φの行列幅計算終了"<<endl;

		//φ-A
		//int width_pa=0;
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=4;j++)
				{
					int N=ELEM[je].node[j];
					N_E_num[N]=N_E_num[N]+1;//
				}
				for(int j=1;j<=6;j++)//2次節点
				{
					int side=ELEM[je].edge[j];
					int N=EDGE[side].para_node_num;
					if(N>0) N_E_num[N]=N_E_num[N]+1;//
				}
			}
		}
		for(int i=1;i<=node2;i++)
		{
			int width=0;//行の未知数幅
			if(N_E_num[i]>0) width=1+3*N_E_num[i];//ここでnume[]が0でないということは、その辺は導体要素に含まれている。要素が1つ追加されるたび、その節点に影響する導体辺は最大3つ増える
			if(npp[i+nedge]>=pn_e) width_mat[npp[i+nedge]+1]+=width;	
			//if(npp[i]<pn) width_edge[npp[i]+1]=width;	
		}
		//cout<<"φ-Aの行列幅計算終了"<<endl;
		

		//φ-φ

		mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_nodeが求まる	

		for(int i=node+1;i<=node2;i++) width_node[i]=2+1;//2次節点は辺の中点なので、自分自身を含め隣り合うのは3点のはず？

		for(int i=1;i<=node2;i++)
		{
			if(conducter_flag[i]==ON)//導体節点
			{
				width_mat[npp[i+nedge]+1]+=width_node[i];
			}
		}
		//cout<<"φ-φの行列幅計算終了"<<endl;
		/////
	}
		
	
	//行列幅の統合
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	delete [] nume;
	
	////配列確保
	cout<<"全体行列用配列宣言"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///全体行列
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double> (0.0,0.0);
			ROW[i][j]=0;
		}
    }
	delete [] width_node;
	delete [] N_E_num;
	delete [] width_edge;
	delete [] width_mat;
	

    complex<double> *B=new complex<double> [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//初期化
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_N;//ディリクレ値として解行列に加算される項の形状関数部
	B_N=complex<double> (0.0,0.0);

	//静磁場項
    for(int je=1;je<=nelm;je++)
    {
		//if(je%1000==0) cout<<je<<endl;
		//辺－節点ﾃｰﾌﾞﾙ作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//重心座標
		double Ys=0;
		double Zs=0;

		double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//積和　x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
		Xs/=4;Ys/=4;Zs/=4;

		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////要素ﾏﾄﾘｸｽ作成開始
		//A-A(静磁場)
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//要素jeの辺番号
			if(EDGE[iside].boundary_condition==0)///未知数
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				//int I1=EDGE[iside].node[1];//isideを構成する2点
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jsideを構成する2点
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///未知数
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int n=dn[jside];
						B_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
						B[I-1]-=complex<double> (B_temp,0);
					}//////////*/
				}
				///強制電流項に関するB[I-1]を計算する
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B_temp=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x.real()*delta6/6;
					B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y.real()*delta6/6;
					B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z.real()*delta6/6;
					B[I-1]+=complex<double> (B_temp,0);
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}

	//渦電流項
	//if(flageddy==ON)
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//渦電流項の計算対象となる要素
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;

			double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//積和　x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//∂A/∂t
		
			//double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=sig[je]*delta*delta6*delta6*delta6*delta6*omega;
			
			//cout<<"導体要素 "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//要素jeの辺番号
				if(EDGE[iside].boundary_condition==0)///未知数
				{   
					int I=npp[iside]+1;///辺isideは行列のI番目
					//int I1=EDGE[iside].node[1];//isideを構成する2点
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jsideを構成する2点
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///未知数
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									
									//∫∫∫zdV=V*Zsだが、∫∫∫z*zdVはV*Zs*Zsではない。(x,yも同様) 教科書p84
									//(Nk)x・(Nu)x

									///////////
									G_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
									G_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
									//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									///////////

									//(Nk)y・(Nu)y
									G_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////

									//(Nk)z・(Nu)z
									G_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									/////*/

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"Aの時間微分項計算で新たに配列確保？"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y・(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z・(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								/////*/
							    
								ROW[I][H]=J;
							}
							///B[I-1]を計算する
							//仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//jw法においては、前のステップの値に基づく項が存在しないため、ここは必要ない

							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int n=dn[jside];

							//∫∫∫zdV=V*Zsだが、∫∫∫z*zdVはV*Zs*Zsではない。(x,yも同様) 教科書p84
									//(Nk)x・(Nu)x

							///////////
							B_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
							B_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
							//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							///////////

							//(Nk)y・(Nu)y
							B_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
							//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							////////

							//(Nk)z・(Nu)z
							B_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
							//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
							if(I==npp[jside]+1)
							{	cout<<"ディリクレ値が未知数のはずの行番号に存在"<<endl;
								B_N=complex<double>(0,PHAT[n]*B_temp*2);//coにはjがかかっているので虚数項へ追加//起こらない？
							}
							else B_N=complex<double>(0,PHAT[n]*B_temp);
							//B[I-1]-=PHAT[n]*B_N;//この記述は大丈夫？
							B[I-1]-=B_N;

									/////*/
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}//////////*/
					}//////*/
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-φ
				double co2=sig[je]*delta*delta6*delta6*delta6;//σV*(1/6V)^3
				double co4=sig[je]*delta6*delta6*delta6*delta6*delta*4;//(1/6V)^2 * 1/9V^2
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//要素jeの辺番号
					if(EDGE[iside].boundary_condition==0)///未知数
					{   
						int I=npp[iside]+1;///辺isideは行列のI番目
						int I1=EDGE[iside].node[1];//isideを構成する2点
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										//gradNiが、二次の形状関数を採用することによって大きく変化する。Ni=Ni*(2*Ni-1)
										//gradNi第一項 2Ni^2 
										//(Nk)x * ∂N∂x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*c[j];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*c[j];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*c[j];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*d[j];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*d[j];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*d[j];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*e[j];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*e[j];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*e[j];

										///////////////										
										//gradNi第二項 -Ni
										G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j]*(-1);
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j]*(-1);
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j]*(-1);


										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									//gradNiが、二次の形状関数を採用することによって大きく変化する。Ni=Ni*(2*Ni-1)
									//gradNi第一項 2Ni^2 
									//(Nk)x * ∂N∂x
									G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*c[j];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*c[j];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*c[j];
										
									///////////

									//(Nk)y * ∂N∂y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*d[j];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*d[j];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[j]*Zs + c[j]*(16*Xs*Zs*ZXs)/20 + d[j]*(16*Ys*Zs+YZs)/20 + e[j]*(4*Zs*Zs/5+ZZs/20))*co4*d[j];

									//(Nk)z * ∂N∂z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*co4*e[j];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[j]*Xs + c[j]*(4*Xs*Xs/5+XXs/20)  + d[j]*(16*Xs*Ys*XYs)/20 + e[j]*(16*Xs*Zs*ZXs)/20)*co4*e[j];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[j]*Ys + c[j]*(16*Xs*Ys*XYs)/20 + d[j]*(4*Ys*Ys/5+YYs/20) + e[j]*(16*Ys*Zs+YZs)/20)*co4*e[j];

									///////////////										
									//gradNi第二項 -Ni
									G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j]*(-1);
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j]*(-1);
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j]*(-1);
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//導体でないということはφ=0なので、ここでAの項のようなディリクレ値に関して考える必要はない
						}

						for(int j=1;j<=6;j++)//追加2次節点
						{	
							int side=ELEM[je].edge[j];
							int N=EDGE[side].para_node_num;
							if(conducter_flag[N]==ON)
							{
								int J=npp[N+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;	

								int I1=EDGE[iside].node[1];//isideを構成する2点
								int I2=EDGE[iside].node[2];
								int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
								int u2=table[j][2];
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										//gradNiが、二次の形状関数を採用することによって大きく変化する。Nij=4*NiNj
										//gradNi第一項 2Ni^2 
										//(Nk)x * ∂N∂x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

										//
										///NiとNjが逆転したの項
										////
										//(Nk)x * ∂N∂x
										G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
										///////////////										

										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									//gradNiが、二次の形状関数を採用することによって大きく変化する。Nij=4*NiNj
									//gradNi第一項 2Ni^2 
									//(Nk)x * ∂N∂x
									G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
									///////////

									//(Nk)y * ∂N∂y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

									//(Nk)z * ∂N∂z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

									//
									///NiとNjが逆転したの項
									////
									//(Nk)x * ∂N∂x
									G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
									///////////

									//(Nk)y * ∂N∂y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

									//(Nk)z * ∂N∂z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//導体でないということはφ=0なので、ここでAの項のようなディリクレ値に関して考える必要はない
						}
					}
				}
				//cout<<"A-φ eddy"<<endl;

				//φ-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"φ-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int J=npp[iside]+1;///辺isideは行列のI番目
							//cout<<"φ-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//isideを構成する2点
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										//gradNiが、二次の形状関数を採用することによって大きく変化する。Ni=Ni*(2*Ni-1)
										//gradNi第一項 2Ni^2 
										//(Nk)x * ∂N∂x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*c[i];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*c[i];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*c[i];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*d[i];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*d[i];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*d[i];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*e[i];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*e[i];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*e[i];

										///////////////										
										//gradNi第二項 -Ni
										G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i]*(-1);
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i]*(-1);
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i]*(-1);

										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									//gradNiが、二次の形状関数を採用することによって大きく変化する。Ni=Ni*(2*Ni-1)
									//gradNi第一項 2Ni^2 
									//(Nk)x * ∂N∂x
									G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*c[i];
									G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*c[i];
									G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*c[i];
										
									///////////

									//(Nk)y * ∂N∂y
									G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*d[i];
									G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*d[i];
									G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*d[i];

									//(Nk)z * ∂N∂z
									G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*e[i];
									G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*e[i];
									G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*e[i];

									///////////////										
									//gradNi第二項 -Ni
									G_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i]*(-1);
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i]*(-1);
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i]*(-1);

									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]の計算
							}
							else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int n=dn[iside];
								//gradNiが、二次の形状関数を採用することによって大きく変化する。Ni=Ni*(2*Ni-1)
								//gradNi第一項 2Ni^2 
								//(Nk)x * ∂N∂x
								B_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*c[i];
								B_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*c[i];
								B_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*c[i];
										
								///////////

								//(Nk)y * ∂N∂y
								B_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*d[i];
								B_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*d[i];
								B_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[i]*Zs + c[i]*(16*Xs*Zs*ZXs)/20 + d[i]*(16*Ys*Zs+YZs)/20 + e[i]*(4*Zs*Zs/5+ZZs/20))*co4*d[i];

								//(Nk)z * ∂N∂z
								B_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[i]+c[i]*Xs+d[i]*Ys+e[i]*Zs)*co4*e[i];
								B_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[i]*Xs + c[i]*(4*Xs*Xs/5+XXs/20)  + d[i]*(16*Xs*Ys*XYs)/20 + e[i]*(16*Xs*Zs*ZXs)/20)*co4*e[i];
								B_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[i]*Ys + c[i]*(16*Xs*Ys*XYs)/20 + d[i]*(4*Ys*Ys/5+YYs/20) + e[i]*(16*Ys*Zs+YZs)/20)*co4*e[i];

								///////////////										
								//gradNi第二項 -Ni
								B_temp+=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i]*(-1);
								B_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i]*(-1);
								B_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i]*(-1);

								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHATは複素ディリクレに備えなおすべき
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				for(int i=1;i<=6;i++)//追加2次節点
				{	
					int side=ELEM[je].edge[i];
					int N=EDGE[side].para_node_num;
					if(conducter_flag[N]==ON)
					{
						int I=npp[N+nedge]+1;
						int I1=EDGE[side].node[1];//isideを構成する2点
						int I2=EDGE[side].node[2];
						int u1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[i][2];
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int J=npp[iside]+1;///辺isideは行列のI番目
							//cout<<"φ-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//isideを構成する2点
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										//gradNiが、二次の形状関数を採用することによって大きく変化する。Nij=4*NiNj
										//gradNi第一項 2Ni^2 
										//(Nk)x * ∂N∂x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

										//
										///NiとNjが逆転したの項
										////
										//(Nk)x * ∂N∂x
										G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
										///////////////						

										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									//gradNiが、二次の形状関数を採用することによって大きく変化する。Nij=4*NiNj
										//gradNi第一項 2Ni^2 
										//(Nk)x * ∂N∂x
										G_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

										//
										///NiとNjが逆転したの項
										////
										//(Nk)x * ∂N∂x
										G_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
										G_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
										G_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
										///////////

										//(Nk)y * ∂N∂y
										G_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
										G_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
										G_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

										//(Nk)z * ∂N∂z
										G_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
										G_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
										G_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
										///////////////					

									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]の計算
							}
							else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int n=dn[iside];
								//gradNiが、二次の形状関数を採用することによって大きく変化する。Nij=4*NiNj
								//gradNi第一項 2Ni^2 
								//(Nk)x * ∂N∂x
								B_temp= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*c[u2];
								B_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*c[u2];
								B_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*c[u2];
										
								///////////

								//(Nk)y * ∂N∂y
								B_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*d[u2];
								B_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*d[u2];
								B_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*Zs + c[u1]*(16*Xs*Zs*ZXs)/20 + d[u1]*(16*Ys*Zs+YZs)/20 + e[u1]*(4*Zs*Zs/5+ZZs/20))*co4*d[u2];

								//(Nk)z * ∂N∂z
								B_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs)*co4*e[u2];
								B_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*Xs + c[u1]*(4*Xs*Xs/5+XXs/20)  + d[u1]*(16*Xs*Ys*XYs)/20 + e[u1]*(16*Xs*Zs*ZXs)/20)*co4*e[u2];
								B_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*Ys + c[u1]*(16*Xs*Ys*XYs)/20 + d[u1]*(4*Ys*Ys/5+YYs/20) + e[u1]*(16*Ys*Zs+YZs)/20)*co4*e[u2];

								//
								///NiとNjが逆転したの項
								////
								//(Nk)x * ∂N∂x
								B_temp+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*c[u1];
								B_temp+=(d[k1]*c[k2]-d[k2]*c[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*c[u1];
								B_temp+=(e[k1]*c[k2]-e[k2]*c[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*c[u1];
										
								///////////

								//(Nk)y * ∂N∂y
								B_temp+=(b[k1]*d[k2]-b[k2]*d[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*d[u1];
								B_temp+=(c[k1]*d[k2]-c[k2]*d[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*d[u1];
								B_temp+=(e[k1]*d[k2]-e[k2]*d[k1])*(b[u2]*Zs + c[u2]*(16*Xs*Zs*ZXs)/20 + d[u2]*(16*Ys*Zs+YZs)/20 + e[u2]*(4*Zs*Zs/5+ZZs/20))*co4*d[u1];

								//(Nk)z * ∂N∂z
								B_temp+=(b[k1]*e[k2]-b[k2]*e[k1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs)*co4*e[u1];
								B_temp+=(c[k1]*e[k2]-c[k2]*e[k1])*(b[u2]*Xs + c[u2]*(4*Xs*Xs/5+XXs/20)  + d[u2]*(16*Xs*Ys*XYs)/20 + e[u2]*(16*Xs*Zs*ZXs)/20)*co4*e[u1];
								B_temp+=(d[k1]*e[k2]-d[k2]*e[k1])*(b[u2]*Ys + c[u2]*(16*Xs*Ys*XYs)/20 + d[u2]*(4*Ys*Ys/5+YYs/20) + e[u2]*(16*Ys*Zs+YZs)/20)*co4*e[u1];
								///////////////				

								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHATは複素ディリクレに備えなおすべき
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"φ-A eddy"<<endl;

				//φ-φ
				/////
				//double co3=sig[je]*dt*delta*delta6*delta6;
				double co3=sig[je]*delta*delta6*delta6/omega;//これに1/jがかかっているため、複素数としては符号を変えて虚数項に足されることになる(1/j=-j)
				double co5=sig[je]*delta*delta6*delta6*delta6*delta6*16/omega;
				double co6=sig[je]*delta*delta6*delta6*delta6*4/omega;
				
				for(int i=1;i<=4;i++)//もとの節点
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{

										G_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
										G_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//相当するＪが存在しなかったら作る
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
										G_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int NN=dn[N[j]];
								B_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

								B_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
								B_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //ふつうの境界条件ではφが値を持たないので、ひとまず削除
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}///////

						for(int j=1;j<=6;j++)//二次節点
						{
							int side=ELEM[je].edge[i];
							int N=EDGE[side].para_node_num;
							if(conducter_flag[N]==ON)
							{
								int flag=0;
								int J=npp[N+nedge]+1;
								int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
								int u2=table[j][2];
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{

										G_temp = co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[i]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(c[i]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(d[i]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(e[i]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp = co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[i]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(c[i]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(d[i]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(e[i]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+=-co6*((c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//相当するＪが存在しなかったら作る
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[i]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(c[i]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(d[i]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(e[i]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

									G_temp = co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[i]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(c[i]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(d[i]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(e[i]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

									G_temp+=-co6*((c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+=-co6*((c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
							{
								//int NN=dn[N[j]];
								cout<<"φが未知数でない?"<<endl;
								B_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

								B_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
								B_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //ふつうの境界条件ではφが値を持たないので、ひとまず削除
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}
					}
				}

				for(int i=1;i<=6;i++)//追加二次節点
				{			
					int side=ELEM[je].edge[i];
					int Ni=EDGE[side].para_node_num;	
					int I=npp[Ni+nedge]+1;
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[Ni]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{

										G_temp = co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[j]*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
										G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(c[j]*(b[k2]*Xs+c[k2]*(4*Xs*Xs/5+XXs/20)+d[k2]*(16*Xs*Ys*XYs)/20 +e[k2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(d[j]*(b[k2]*Ys+c[k2]*(16*Xs*Ys*XYs)/20 +d[k2]*(4*Ys*Ys/5+YYs/20) +e[k2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(e[j]*(b[k2]*Zs+c[k2]*(16*Xs*Zs*ZXs)/20 +d[k2]*(4*Ys*Zs/5+YZs/20) +e[k2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp = co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[j]*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
										G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(c[j]*(b[k1]*Xs+c[k1]*(4*Xs*Xs/5+XXs/20)+d[k1]*(16*Xs*Ys*XYs)/20 +e[k1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(d[j]*(b[k1]*Ys+c[k1]*(16*Xs*Ys*XYs)/20 +d[k1]*(4*Ys*Ys/5+YYs/20) +e[k1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(e[j]*(b[k1]*Zs+c[k1]*(16*Xs*Zs*ZXs)/20 +d[k1]*(4*Ys*Zs/5+YZs/20) +e[k1]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+=-co6*((c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
										G_temp+=-co6*((c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));

										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//相当するＪが存在しなかったら作る
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[j]*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
									G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(c[j]*(b[k2]*Xs+c[k2]*(4*Xs*Xs/5+XXs/20)+d[k2]*(16*Xs*Ys*XYs)/20 +e[k2]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(d[j]*(b[k2]*Ys+c[k2]*(16*Xs*Ys*XYs)/20 +d[k2]*(4*Ys*Ys/5+YYs/20) +e[k2]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(e[j]*(b[k2]*Zs+c[k2]*(16*Xs*Zs*ZXs)/20 +d[k2]*(4*Ys*Zs/5+YZs/20) +e[k2]*(4*Zs*Zs/5+ZZs/20)));

									G_temp = co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[j]*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
									G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(c[j]*(b[k1]*Xs+c[k1]*(4*Xs*Xs/5+XXs/20)+d[k1]*(16*Xs*Ys*XYs)/20 +e[k1]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(d[j]*(b[k1]*Ys+c[k1]*(16*Xs*Ys*XYs)/20 +d[k1]*(4*Ys*Ys/5+YYs/20) +e[k1]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(e[j]*(b[k1]*Zs+c[k1]*(16*Xs*Zs*ZXs)/20 +d[k1]*(4*Ys*Zs/5+YZs/20) +e[k1]*(4*Zs*Zs/5+ZZs/20)));

									G_temp+=-co6*((c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
									G_temp+=-co6*((c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));

									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
							{
								//cout<<"追加2次要素節点がディリクレ？"<<endl;
								int NN=dn[N[j]];
								B_temp = co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[j]*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
								B_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(c[j]*(b[k2]*Xs+c[k2]*(4*Xs*Xs/5+XXs/20)+d[k2]*(16*Xs*Ys*XYs)/20 +e[k2]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(d[j]*(b[k2]*Ys+c[k2]*(16*Xs*Ys*XYs)/20 +d[k2]*(4*Ys*Ys/5+YYs/20) +e[k2]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(e[j]*(b[k2]*Zs+c[k2]*(16*Xs*Zs*ZXs)/20 +d[k2]*(4*Ys*Zs/5+YZs/20) +e[k2]*(4*Zs*Zs/5+ZZs/20)));

								B_temp = co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[j]*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
								B_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(c[j]*(b[k1]*Xs+c[k1]*(4*Xs*Xs/5+XXs/20)+d[k1]*(16*Xs*Ys*XYs)/20 +e[k1]*(16*Xs*Zs*ZXs)/20));
								B_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(d[j]*(b[k1]*Ys+c[k1]*(16*Xs*Ys*XYs)/20 +d[k1]*(4*Ys*Ys/5+YYs/20) +e[k1]*(16*Ys*Zs*YZs)/20));
								B_temp+= co5*(c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(e[j]*(b[k1]*Zs+c[k1]*(16*Xs*Zs*ZXs)/20 +d[k1]*(4*Ys*Zs/5+YZs/20) +e[k1]*(4*Zs*Zs/5+ZZs/20)));

								B_temp+=-co6*((c[j]*c[k1]+d[j]*d[k1]+e[j]*e[k1])*(b[k2]+c[k2]*Xs+d[k2]*Ys+e[k2]*Zs));
								B_temp+=-co6*((c[j]*c[k2]+d[j]*d[k2]+e[j]*e[k2])*(b[k1]+c[k1]*Xs+d[k1]*Ys+e[k1]*Zs));
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //ふつうの境界条件ではφが値を持たないので、ひとまず削除
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}///////

						for(int j=1;j<=6;j++)//二次節点
						{
							int side=ELEM[je].edge[i];
							int N=EDGE[side].para_node_num;
							if(conducter_flag[N]==ON)
							{
								int flag=0;
								int J=npp[N+nedge]+1;
								int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
								int u2=table[j][2];
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{

										G_temp = co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(b[k2]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+= co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(c[k2]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(d[k2]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k1]*c[u1]+d[k1]*d[u1]+e[k1]*e[u1])*(e[k2]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(b[k2]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(c[k2]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(d[k2]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k1]*c[u2]+d[k1]*d[u2]+e[k1]*e[u2])*(e[k2]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(b[k1]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(c[k1]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(d[k1]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k2]*c[u1]+d[k2]*d[u1]+e[k2]*e[u1])*(e[k1]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(b[k1]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(c[k1]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(d[k1]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
										G_temp+= co5*(c[k2]*c[u2]+d[k2]*d[u2]+e[k2]*e[u2])*(e[k1]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

										
										
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//相当するＪが存在しなかったら作る
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp = co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[i]*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(c[i]*(b[u2]*Xs+c[u2]*(4*Xs*Xs/5+XXs/20)+d[u2]*(16*Xs*Ys*XYs)/20 +e[u2]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(d[i]*(b[u2]*Ys+c[u2]*(16*Xs*Ys*XYs)/20 +d[u2]*(4*Ys*Ys/5+YYs/20) +e[u2]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(e[i]*(b[u2]*Zs+c[u2]*(16*Xs*Zs*ZXs)/20 +d[u2]*(4*Ys*Zs/5+YZs/20) +e[u2]*(4*Zs*Zs/5+ZZs/20)));

									G_temp = co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[i]*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(c[i]*(b[u1]*Xs+c[u1]*(4*Xs*Xs/5+XXs/20)+d[u1]*(16*Xs*Ys*XYs)/20 +e[u1]*(16*Xs*Zs*ZXs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(d[i]*(b[u1]*Ys+c[u1]*(16*Xs*Ys*XYs)/20 +d[u1]*(4*Ys*Ys/5+YYs/20) +e[u1]*(16*Ys*Zs*YZs)/20));
									G_temp+= co5*(c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(e[i]*(b[u1]*Zs+c[u1]*(16*Xs*Zs*ZXs)/20 +d[u1]*(4*Ys*Zs/5+YZs/20) +e[u1]*(4*Zs*Zs/5+ZZs/20)));

									G_temp+=-co6*((c[i]*c[u1]+d[i]*d[u1]+e[i]*e[u1])*(b[u2]+c[u2]*Xs+d[u2]*Ys+e[u2]*Zs));
									G_temp+=-co6*((c[i]*c[u2]+d[i]*d[u2]+e[i]*e[u2])*(b[u1]+c[u1]*Xs+d[u1]*Ys+e[u1]*Zs));
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int NN=dn[N];
								cout<<"φがディリクレ値?"<<endl;
								B_temp = co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]*(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs));
										B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(c[i]*(b[j]*Xs+c[j]*(4*Xs*Xs/5+XXs/20)+d[j]*(16*Xs*Ys*XYs)/20 +e[j]*(16*Xs*Zs*ZXs)/20));
										B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(d[i]*(b[j]*Ys+c[j]*(16*Xs*Ys*XYs)/20 +d[j]*(4*Ys*Ys/5+YYs/20) +e[j]*(16*Ys*Zs*YZs)/20));
										B_temp+= co5*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(e[i]*(b[j]*Zs+c[j]*(16*Xs*Zs*ZXs)/20 +d[j]*(4*Ys*Zs/5+YZs/20) +e[j]*(4*Zs*Zs/5+ZZs/20)));

										B_temp+=-co6*((c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*(b[i]+b[j]+(c[i]+c[j])*Xs+(d[i]+d[j])*Ys+(e[i]+e[j])*Zs));
										B_temp+= co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //ふつうの境界条件ではφが値を持たないので、ひとまず削除
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}
					}
				}
				////*/
				//cout<<"φ-φ eddy"<<endl;
			}
		}
	}
	

	///////
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"非ゼロ要素数"<<number<<endl;

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix_complex(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"行列の対称性チェック----";
	check_matrix_symmetry_complex(pn,NUM,ROW,G);
	cout<<"完了"<<endl;;
	//解行列の値チェック
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
	complex<double> *val = new complex<double> [number];
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
    /////////////////////

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG法実行
	complex<double> *XX=new complex<double> [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	//XXに格納された解を各辺に振る
	//cout<<"複素数解の割り振り"<<endl;

	double *AR = new double [nedge+1];
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];

	for(int i=1;i<=nedge;i++)
	{
		
		AR[i]=0;
		AI[i]=0;
	}
	for(int i=1;i<=node;i++)
	{
		VR[i]=0;
		VI[i]=0;
	}

	ofstream ar("old_AR_e.dat");
	ofstream ai("old_AI_e.dat");
	ofstream vr("old_VR_e.dat");
	ofstream vi("old_VI_e.dat");
	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[i]=XX[n].real();
			AI[i]=XX[n].imag();
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			VR[i]=XX[n].real();
			VI[i]=XX[n].imag();
		}
	}	
	
	delete [] XX;

	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	///Am,φのﾌｧｲﾙ出力 φは電位ではなく位相の遅れ
	ofstream a("Am_e.dat");
	ofstream p("phi_e.dat");
	double Am=0.0;
	double Vm=0.0;
	double phi=0.0;
	for(int i=1;i<=nedge;i++)
	{
		if(EDGE[i].boundary_condition==0)
		{
			Am=sqrt(AR[i]*AR[i]+AI[i]*AI[i]);
			phi=atan2(AI[i],AR[i]);
			
			//if(CON->get_jw_Faverage()==ON) A[i]=Am/sqrt(2.0);
			//else A[i]=Am*cos(omega*t+phi);
			A[i]=Am*cos(omega*TIME+phi);
		}
		a<<Am<<endl;
		p<<phi<<endl;
		ar<<AR[i]<<endl;
		ai<<AI[i]<<endl;
	}
	for(int i=1;i<=node;i++)
	{
		vr<<VR[i]<<endl;
		vi<<VI[i]<<endl;
	}

	a.close();
	p.close();
	ar.close();
	ai.close();
	vr.close();
	vi.close();
	//cout<<"複素数を実数変換完了"<<endl;
	/*///
	for(int i=1;i<=node;i++)
	{
		Vm=sqrt(VR[i]*VR[i]+VI[i]*VI[i]);
		phi=atan(VI[i]/VR[i]);
		V[i]=Vm*cos(omega*t+phi);
	}
	*/
	
	
	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//渦電流損[W]
	double *Je_loss_n=new double[node+1];//渦電流損[W]

	calc_edge_eddy_current_jw(CON,NODE,ELEM,EDGE,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t,TIME,sig,omega);
	
	
	//渦電流損を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//本当の体積

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//節点の渦電流損を対応する粒子へ渡す
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//
	
	//if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
    delete [] PHAT;

    delete [] B;

	delete [] AR;
	delete [] AI;
	delete [] VR;
	delete [] VI;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	delete [] conduct_edge_flag;
	
	//////
	/////*/
}


//辺要素による渦電流を考慮した磁気ベクトルポテンシャル計算関数,jw法 静的要素のディリクレ値読み込み
void Avector3D_edge_eddy_jw2(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double *A,int *jnb,int **nei,double *old_A,double dt,double *V,double **current,double *RP,vector<mpsparticle> &PART,double *sig,int t,double TIME,double omega, int node_sta, vector<edge3D> &static_EDGE,double *Am,double *phi)
{
	//ELEM.edge:要素を構成する辺６つ　EDGE.node:辺をつくる点２つ　
	/////
	int flageddy=OFF;//A-φ法にするかしないか

	cout<<"渦電流を考慮にいれた辺要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始(jw_static読み込み)"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	//double j0x,j0y,j0z;					//電流密度
	double Sx=0;
	double Sy=0;
	double Sz=0;
	complex<double> Im=(0,1);
	
	complex<double> j0x;
	complex<double> j0y;
	complex<double> j0z;
	//double sigma=CON->get_ele_conduc();	//電気伝導率
	unsigned timeA=GetTickCount();	//計算開始時刻

	int NN=0;//ディリクレ型境界辺数
	int NN_V=0;//ディリクレ型境界節点数
    int *dn=new int [nedge+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnedge+1を格納
	double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
	
	complex<double> *PHAT_A=new complex<double> [nedge+1];//Aのディリクレ型境値
	complex<double> *PHAT_V = new complex<double> [node+1];//Vのディリクレ境界値

	///解析領域の境界に固定境界条件を設定

	////
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iedge=ELEM[i].edge[k];
						int ia=EDGE[iedge].node[1];
						int ib=EDGE[iedge].node[2];
						if(ia!=p && ib!=p) EDGE[iedge].boundary_condition=1;//節点pを含まない辺は境界辺
						else EDGE[iedge].boundary_condition=0;
					}
				}
			}
			
		}
	}
	cout<<"辺要素の固定境界設定完了"<<endl;


	///ディリクレ型境界条件入力
	//set_boundary_condition3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,dn,&NN,PHAT,A);
    

	//静的要素のディリクレ値読み込み
	double *AR = new double [nedge+1];//ベクトルポテンシャル実部
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	//初期化
	for(int i=1;i<=nedge;i++)
	{
		//AR[i]=A[i];//何かしらの0でない境界条件(一様磁場など)が入っている場合を考え、ここで入れておく
		AR[i]=0;	
		AI[i]=0;
	}
	for(int i=1;i<=node;i++)
	{
		VR[i]=0;
		VI[i]=0;
	}

	if(CON->get_static_dirichlet()==ON)
	{
		/////////1ステップ目の解読み込み
		cout<<"静的要素のディリクレ値読み込み"<<endl;

		//AR
		ifstream ar("old_AR_e.dat");
		if(!ar) cout<<"cannot open old_AR_e.dat"<<endl;
		ar.unsetf(ifstream::dec);
		//ar.setf(ifstream::skipws);
		for(int i=1;i<=nedge;i++) ar>>AR[i];
		ar.close();

		//AI
		ifstream ai("old_AI_e.dat");
		if(!ai) cout<<"cannot open old_AI_e.dat"<<endl;
		ai.unsetf(ifstream::dec);
		//ai.setf(ifstream::skipws);
		for(int i=1;i<=nedge;i++) ai>>AI[i];
		ai.close();

		//VR
		ifstream vr("old_VR_e.dat");
		if(!vr) cout<<"cannot open old_VR.dat"<<endl;
		//vr.unsetf(ifstream::dec);
		//vr.setf(ifstream::skipws);
		for(int i=1;i<=node;i++) vr>>VR[i];
		vr.close();

		//VI
		ifstream vi("old_VI_e.dat");
		if(!vi) cout<<"cannot open old_VI.dat"<<endl;
		//vi.unsetf(ifstream::dec);
		//vi.setf(ifstream::skipws);
		for(int i=1;i<=node;i++) vi>>VI[i];
		vi.close();
		/////
		cout<<"ディリクレ値読み込み完了"<<endl;		
		int count_r=0;

		
		////

		for(int i=1;i<=nedge;i++)
		{
			/*///この関数にはstatic_dirichlet=OFFのとき飛んで来ない
			if(CON->get_static_dirichlet()==OFF)
			{
				if(NODE[i].boundary_condition==2)
				{    
					dn[i]=NN;
					PHAT[A_X][NN]=0;
					PHAT[A_Y][NN]=0;
					PHAT[A_Z][NN]=0;
					A[A_X][i]=0;
					A[A_Y][i]=0;
					A[A_Z][i]=0;
					NN++;
				}
				else if(NODE[i].boundary_condition==1)
				{    
					dn[i]=NN;
					PHAT[A_X][NN]=0;
					PHAT[A_Y][NN]=0;
					PHAT[A_Z][NN]=0;
					A[A_X][i]=0;
					A[A_Y][i]=0;
					A[A_Z][i]=0;
					NN++;
				}
				else
				{
					dn[i]=node+1;
				}
			}
			////*/

			
			{
				//if(NODE[i].remesh==OFF)//初期ステップ以外はリメッシュしない節点について最初に求めたベクトルポテンシャルをディリクレ値として与える
				if(EDGE[i].stat==ON)
				{
					count_r++;
					EDGE[i].boundary_condition=3;//静的要素の境界条件をもとの0から3に書き換える。固定境界(1,2)も変えるが問題ないはず
					dn[i]=NN;
					PHAT_A[NN]=complex<double>(AR[EDGE[i].static_num],AI[EDGE[i].static_num]);
					//PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					NN++;
				}
				else
				{
					dn[i]=nedge+1;
				}
			}
			//*/
		}//////////////
		
			for(int i=1;i<=node;i++)
			{
				if(i<=node_sta)
				{
					NODE[i].boundary_condition=3;
					PHAT_V[NN]=complex<double>(VR[i],VI[i]);
					if(CON->get_Je_crucible()==0)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!流体が静的節点"<<endl;
					}
					if(CON->get_Je_crucible()==1)
					{
						if(NODE[i].material==FLUID && jnb[i]!=0) cout<<"eroor!流体が静的節点"<<endl;
						if(NODE[i].material==CRUCIBLE &&jnb[i]!=0)
						{
							NN_V++;
						}
					}
					if(CON->get_J0eqJe()==1) NN_V=0;
				}
				else
				{
					//dn[i]=node+1;
				}
			}
		
		cout<<"ﾃﾞｨﾘｸﾚ辺数＝"<<NN<<"ﾃﾞｨﾘｸﾚ数＝"<<NN+NN_V<<" NN_V="<<NN_V<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
		//cout<<"non_remesh辺数="<<count_r<<endl;
	}

    //cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<endl;
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える
	//辺要素でもφについては節点要素を用いる

	int conducter_num=0;//導体節点数
	int *conducter_flag=new int [node+1];//導体節点ならON
	for(int i=0;i<=node;i++) conducter_flag[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(CON->get_Je_crucible()==0)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}

		if(CON->get_Je_crucible()==1)
		{
			if(NODE[i].material==FLUID && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
			if(NODE[i].material==CRUCIBLE && NODE[i].boundary_condition==0 &&jnb[i]!=0)
			{
				conducter_num++;
				conducter_flag[i]=ON;
			}
		}
	}
	if(CON->get_Je_crucible()==-1) conducter_num=0;
	//////////////
	cout<<"conducter_num="<<conducter_num<<endl;
	int pn_e=nedge-NN;///辺要素の未知数 複素数のまま格納するので時間差分の２倍にはならない
	int pn=pn_e;//φも含めた全未知数
	if(flageddy==ON) pn=pn_e+conducter_num;//φも含めた全未知数
	//int pn=pn_e;

	//int *ppn=new int [pn];		///行列でn番目は辺番号ppn[n] 行列でn(ただし、nは辺数より大きい)番目は節点番号ppn[n](導体のみ格納)
    //int *npp=new int [nedge+node+1];	///各辺,導体節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    vector <int> ppn;
	vector <int> npp;
	

	npp.reserve(nedge+node+1);
	int num=0; 
   
	//for(int i=0;i<=nedge+node;i++) npp[i]=0;

	for(int i=1;i<=nedge;i++)//A
    {
        if(EDGE[i].boundary_condition==0)//未知数
		{
			//ppn[num]=i;
			ppn.push_back(i);
			npp[i]=num;
			//npp.push_back(num);
			num++;
		}
		else npp[i]=pn+1;// npp.push_back(pn+1);   //
    }
	
	////
	if(flageddy==ON)
	{
		for(int i=1;i<=node;i++)//φ
		{
			if(NODE[i].boundary_condition==0 && conducter_flag[i]==ON)//導体節点
			{
				//ppn[num]=nedge+i;//
				ppn.push_back(i);//行列でn(ただし、nは辺数より大きい)番目は節点番号ppn[n]
				npp[nedge+i]=num;//節点要素のnppは、辺要素のnppをすべて格納したあとの配列を用いる
				//npp.push_back(num);
				num++;
			}
			else npp[nedge+i]=pn+1;// npp.push_back(pn+1); //
		}
	}
	/////
	cout<<"pn="<<pn<<" num="<<num<<endl;

	////行列の幅計算 辺要素
	//A-A
    int mat_we=0;
	int *nume=new int [nedge+1];				//EDGE[i]を含む要素数
	int *width_edge= new int [pn+1];

	for(int i=1;i<=nedge;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int side=ELEM[i].edge[j];
			nume[side]=nume[side]+1;
		}
	}											///nume[i]がもとまった
	for(int i=1;i<=nedge;i++)
	{
		if(EDGE[i].boundary_condition==0)
		{
			int width=7+3*(nume[i]-2);
			if(width>mat_we) mat_we=width;
			if(npp[i]<pn) width_edge[npp[i]+1]=width;
		}
	}

	int nume_max=0;
	for(int i=1;i<=nedge;i++) if(nume[i]>nume_max) nume_max=nume[i];
	cout<<"nume_max="<<nume_max<<endl;

	
	cout<<"A-Aの行列幅計算終了"<<endl;

	//節点要素関係の行列幅変数の宣言
	int mat_wn=0;
	int *width_node=new int [node+1]; ///各節点の隣り合う節点の数＋１
	int *width_mat=new int [pn+1]; ///各行の幅

	for(int i=0;i<=node;i++) 
	{
		width_node[i]=0;
	}

	for(int i=0;i<=pn;i++) 
	{
		width_mat[i]=0;
	}

	if(flageddy==ON)
	{
		////行列の幅計算 渦電流項
		//φ-φ
		
		mat_wn=calc_matrix_width2(CON,NODE,ELEM,node,nelm,jnb,nei,width_node);//width_nodeが求まる
	
		//各行の節点要素由来の非零数を求める width_matに格納
	
	
	
	
		//A-φ
		int **ROW2=new int *[pn+1];
		for(int i=1;i<=pn;i++) ROW2[i]=new int [nume_max*4+1];
		for(int i=1;i<=pn;i++)//初期化
		{
			for(int j=1;j<=nume_max*4;j++)
			{
				ROW2[i][j]=0;
			}
		}
		int N2[4+1]; //要素の各節点番号格納	
		
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//要素jeの辺番号
					if(EDGE[iside].boundary_condition==0)///未知数
					{   
						int I=npp[iside]+1;///辺isideは行列のI番目
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N2[j]]==ON)
							{
								int J=npp[N2[j]+nedge]+1;
								int flag=0;			
			
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
								}
							}
						}
					}
				}
			}
		}
		cout<<"A-φの行列幅計算終了"<<endl;
		
		//φ-A
		for(int je=1;je<=nelm;je++)
		{
			for(int j=1;j<=4;j++) N2[j]=ELEM[je].node[j];
			int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}
			if(flagje==ON)
			{	
				
				for(int i=1;i<=4;i++)
				{
					
					
					if(conducter_flag[N2[i]]==ON)
					{
						int I=npp[N2[i]+nedge]+1;
						if(I<=pn_e) cout<<"I<=pn_e I="<<I<<endl;
						for(int j=1;j<=6;j++)
						{	
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								int J=npp[iside]+1;///辺isideは行列のI番目
								for(int h=1;h<=width_mat[I];h++)
								{
									if(J==ROW2[I][h]) flag=1;
								}
								if(flag==0)
								{   
									width_mat[I]=width_mat[I]+1;
									int H=width_mat[I];
									ROW2[I][H]=J;
									//B[I-1]の計算

								}
							}		
						}
					}
				}
			}
		}



		cout<<"A-φ,φ-Aの行列幅計算終了"<<endl;

		//φ-φ
		for(int i=1;i<=node;i++)
		{
			if(conducter_flag[i]==ON)//導体節点
			{
				width_mat[npp[i]+1]+=width_node[i];
			}
		}
		cout<<"φ-φの行列幅計算終了"<<endl;

		for(int i=1;i<=pn;i++) delete [] ROW2[i];
		delete [] ROW2;
	}
	
	//行列幅の統合
	for(int i=1;i<=pn;i++) width_mat[i]+=width_edge[i];
	
	delete [] nume;
	
	
	////配列確保
	//cout<<"全体行列用配列宣言"<<endl;
	complex<double> **G=new complex<double> *[pn+1];///全体行列
    //for(int i=1;i<=pn;i++) G[i]=new double [width_mat[i]+1];
	for(int i=1;i<=pn;i++) G[i]=new complex<double> [width_mat[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    //for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
	for(int i=1;i<=pn;i++) ROW[i]=new int [width_mat[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数


    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=width_mat[i];j++)
		{
			G[i][j]=complex<double> (0.0,0.0);
			ROW[i][j]=0;
		}
    }
	
	delete [] width_node;
	delete [] width_edge;
	delete [] width_mat;

    complex<double> *B=new complex<double> [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=complex<double>(0.0,0.0);//初期化
    ////
    
	cout<<"pn_e="<<pn_e<<" pn="<<pn<<" nedge="<<nedge<<endl;
	cout<<"全体行列作成開始"<<endl;
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//int J,J2,J3,flag,flag2,flag3;
	double G_temp=0.0;
	double B_temp=0.0;
	complex<double> B_N;//ディリクレ値として解行列に加算される項の形状関数部
	B_N=complex<double> (0.0,0.0);

	//静磁場項
    for(int je=1;je<=nelm;je++)
    {
		//if(t>1 && je%1000==0) cout<<je<<endl;
		//辺－節点ﾃｰﾌﾞﾙ作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//重心座標
		double Ys=0;
		double Zs=0;

		double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
		double YYs=0;
		double ZZs=0;

		double XYs=0;//積和　x1*y1+x2*y2+...
		double YZs=0;
		double ZXs=0;
		
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
			//XXs+=X[j]*X[j];
			//YYs+=Y[j]*Y[j];
			//ZZs+=Z[j]*Z[j];
			
			//XYs+=X[j]*Y[j];
			//YZs+=Y[j]*Z[j];
			//ZXs+=Z[j]*X[j];
		}
		Xs/=4;Ys/=4;Zs/=4;

		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6;//delta6=1/6V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		if(ELEM[je].material==COIL)
		{
			//1ステップ目、すなわちj0が最大になっているときのみ成り立つ。あとでcurrentの定義を修正すること
			j0x=complex<double> (current[A_X][je],0.0);
			j0y=complex<double> (current[A_Y][je],0.0);
			j0z=complex<double> (current[A_Z][je],0.0);
		}
		else
		{
			j0x=complex<double> (0,0);
			j0y=complex<double> (0,0);
			j0z=complex<double> (0,0);
		}

		///比透磁率
		double v=v0/RP[je];
		double rp=u0*RP[je];

		////要素ﾏﾄﾘｸｽ作成開始
		//A-A(静磁場)
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].edge[i];//要素jeの辺番号
			if(EDGE[iside].boundary_condition==0)///未知数
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				//int I1=EDGE[iside].node[1];//isideを構成する2点
				//int I2=EDGE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].edge[j];
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
					//int J1=EDGE[jside].node[1];//jsideを構成する2点
					//int J2=EDGE[jside].node[2];
						
					if(EDGE[jside].boundary_condition==0)///未知数
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
								G[I][h]+=complex<double> (G_temp,0);
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    G_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
							G[I][H]+=complex<double> (G_temp,0);
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがﾃﾞｨﾘｸﾚ型境界辺なら
					{
						//if(EDGE[jside].boundary_condition!=3) cout<<"bcが3以外？"<<endl;
					    int n=dn[jside];
						B_temp=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*(2.0/3.0)*delta6*delta6*delta6/rp;
						B_N=complex<double>(B_temp,0); 
						B[I-1]-=PHAT_A[n]*B_N;
					}//////////*/
				}
				///強制電流項に関するB[I-1]を計算する
				if(ELEM[je].material==COIL)
				{	
					
					//B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					//B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					//B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
					B_temp=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x.real()*delta6/6;
					B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y.real()*delta6/6;
					B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z.real()*delta6/6;
					B[I-1]+=complex<double> (B_temp,0);
					
				}//////////
				//else if(ELEM[je].material==MAGNET)
				//{
				//	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				//}
			}
		} 
	}
	cout<<"静磁場項の作成完了"<<endl;

	//渦電流項
	//if(flageddy==ON)
	{
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//注目している要素で渦電流項を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)//渦電流項の計算対象となる要素
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;

			double XXs=0;//二乗和　x1*x1+x2*x2+x3*x3+x4*x4
			double YYs=0;
			double ZZs=0;

			double XYs=0;//積和　x1*y1+x2*y2+...
			double YZs=0;
			double ZXs=0;
			
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
				
				XXs+=X[j]*X[j];
				YYs+=Y[j]*Y[j];
				ZZs+=Z[j]*Z[j];
				
				XYs+=X[j]*Y[j];
				YZs+=Y[j]*Z[j];
				ZXs+=Z[j]*X[j];
			}
	    
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//delta6=1/6V
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
		
			//∂A/∂t
		
			//double co=sig[je]*delta*delta6*delta6*delta6*delta6/dt;
			//ｊω法の場合、1/dtをｊωに置き換えればよい。ｊは各成分作成時につじつまを合わせる
			double co=sig[je]*delta*delta6*delta6*delta6*delta6*omega;
			
			//cout<<"導体要素 "<<je<<endl;
			for(int i=1;i<=6;i++)
			{			
				int iside=ELEM[je].edge[i];//要素jeの辺番号
				if(EDGE[iside].boundary_condition==0)///未知数
				{   
					int I=npp[iside]+1;///辺isideは行列のI番目
					//int I1=EDGE[iside].node[1];//isideを構成する2点
					//int I2=EDGE[iside].node[2];
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
						//int J1=EDGE[jside].node[1];//jsideを構成する2点
						//int J2=EDGE[jside].node[2];
							
						if(EDGE[jside].boundary_condition==0)///未知数
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							
							int flag=0;
							
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									
									//∫∫∫zdV=V*Zsだが、∫∫∫z*zdVはV*Zs*Zsではない。(x,yも同様) 教科書p84
									//(Nk)x・(Nu)x

									///////////
									G_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
									G_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
									G_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
									G_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
									//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
									///////////

									//(Nk)y・(Nu)y
									G_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
									G_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
									G_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
									//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
									////////

									//(Nk)z・(Nu)z
									G_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
									G_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
									G_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
									G_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
									//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
									if(I==J) G[I][h]+=complex<double>(0,G_temp*2);//coにはjがかかっているので虚数項へ追加
									else G[I][h]+=complex<double>(0,G_temp);
									/////*/

									flag=1;
								}
							}
							if(flag==0)
							{  
								cout<<"起こらないはず"<<endl;
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+= (b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(4*Ys*4*Ys+YYs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(4*Zs*4*Zs+ZZs)*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
								//(Nk)y・(Nu)y
								G[I][H]+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*sig[je]*Xs*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*4*Zs+ZZs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
								//(Nk)z・(Nu)z
								G[I][H]+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*sig[je]*delta*delta6*delta6*delta6*delta6/dt;
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(16*Xs*Ys+XYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*4*Xs+XXs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								G[I][H]+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*4*Ys+YYs)*sig[je]*delta*delta6*delta6*delta6*delta6/(dt*20);
								//G[I][h]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
								/////*/
							    
								ROW[I][H]=J;
							}
							///B[I-1]を計算する
							//仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							//jw法においては、前のステップの値に基づく項が存在しないため、ここは必要ない

							//double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							//Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							//Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							//Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int n=dn[jside];

							//∫∫∫zdV=V*Zsだが、∫∫∫z*zdVはV*Zs*Zsではない。(x,yも同様) 教科書p84
									//(Nk)x・(Nu)x

							///////////
							B_temp=(b[k1]*c[k2]-b[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1])*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1])+(d[k1]*c[k2]-d[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Ys*co;
							B_temp+=((b[k1]*c[k2]-b[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(b[u1]*c[u2]-b[u2]*c[u1]))*Zs*co;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1])+(e[k1]*c[k2]-e[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Zs+YZs)*co/20;
							B_temp+=((d[k1]*c[k2]-d[k2]*c[k1])*(d[u1]*c[u2]-d[u2]*c[u1]))*(16*Ys*Ys+YYs)*co/20;
							B_temp+=((e[k1]*c[k2]-e[k2]*c[k1])*(e[u1]*c[u2]-e[u2]*c[u1]))*(16*Zs*Zs+ZZs)*co/20;
							//G_temp+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							///////////

							//(Nk)y・(Nu)y
							B_temp+= (b[k1]*d[k2]-b[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1])*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1])+(c[k1]*d[k2]-c[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Xs*co;
							B_temp+=((b[k1]*d[k2]-b[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(b[u1]*d[u2]-b[u2]*d[u1]))*Zs*co;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1])+(e[k1]*d[k2]-e[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(16*Xs*Zs+ZXs)*co/20;
							B_temp+=((c[k1]*d[k2]-c[k2]*d[k1])*(c[u1]*d[u2]-c[u2]*d[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((e[k1]*d[k2]-e[k2]*d[k1])*(e[u1]*d[u2]-e[u2]*d[u1]))*(4*Zs*Zs/5+ZZs/20)*co;
							//G_temp+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							////////

							//(Nk)z・(Nu)z
							B_temp+= (b[k1]*e[k2]-b[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1])*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1])+(c[k1]*e[k2]-c[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Xs*co;
							B_temp+=((b[k1]*e[k2]-b[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(b[u1]*e[u2]-b[u2]*e[u1]))*Ys*co;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1])+(d[k1]*e[k2]-d[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(16*Xs*Ys+XYs)*co/20;
							B_temp+=((c[k1]*e[k2]-c[k2]*e[k1])*(c[u1]*e[u2]-c[u2]*e[u1]))*(4*Xs*Xs/5+XXs/20)*co;
							B_temp+=((d[k1]*e[k2]-d[k2]*e[k1])*(d[u1]*e[u2]-d[u2]*e[u1]))*(4*Ys*Ys/5+YYs/20)*co;
							//G_temp+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
									
							if(I==npp[jside]+1)
							{	cout<<"ディリクレ値が未知数のはずの行番号に存在"<<endl;
								//B_N=complex<double>(0,PHAT[n]*B_temp*2);//coにはjがかかっているので虚数項へ追加//起こらない？
								B_N=complex<double>(0,B_temp*2);
							}
							else B_N=complex<double>(0,B_temp);
							//else B_N=complex<double>(0,PHAT[n]*B_temp);
							//B[I-1]-=PHAT[n]*B_N;//この記述は大丈夫？
							B[I-1]-=PHAT_A[n]*B_N;

									/////*/
							//B[I-1]-=co*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*((b[u1]*c[u2]-b[u2]*c[u1])+(d[u1]*c[u2]-c[u1]*d[u2])*Ys+(e[u1]*c[u2]-c[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*((b[u1]*d[u2]-b[u2]*d[u1])+(c[u1]*d[u2]-d[u1]*c[u2])*Xs+(e[u1]*d[u2]-d[u1]*e[u2])*Zs);
							//B[I-1]-=co*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*((b[u1]*e[u2]-b[u2]*e[u1])+(c[u1]*e[u2]-e[u1]*c[u2])*Xs+(d[u1]*e[u2]-e[u1]*d[u2])*Ys);
							
						}//////////*/
					}//////*/
					//B[I-1]+=co*(Sx+Sy+Sz);
				}
			}
			//cout<<"A-A eddy"<<endl;
			
			if(flageddy==ON)
			{
				//A-φ
				double co2=sig[je]*delta*delta6*delta6*delta6;//σV*(1/6V)^3
				for(int i=1;i<=6;i++)
				{			
					int iside=ELEM[je].edge[i];//要素jeの辺番号
					if(EDGE[iside].boundary_condition==0)///未知数
					{   
						int I=npp[iside]+1;///辺isideは行列のI番目
						int I1=EDGE[iside].node[1];//isideを構成する2点
						int I2=EDGE[iside].node[2];
						int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int k2=table[i][2];
					
						for(int j=1;j<=4;j++)
						{	
							if(conducter_flag[N[j]]==ON)
							{
								int J=npp[N[j]+nedge]+1;
								if(J==pn+1) cout<<"eroor J=pn+1"<<endl;
								int flag=0;			
								
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
										G[I][h]+=complex<double> (G_temp,0);

										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[j];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[j];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[j];
									G[I][H]+=complex<double> (G_temp,0);
									ROW[I][H]=J;
									
								}
							}
							//導体でないということはφ=0なので、ここでAの項のようなディリクレ値に関して考える必要はない
						}
					}
				}
				//cout<<"A-φ eddy"<<endl;

				//φ-A
				for(int i=1;i<=4;i++)
				{
					
					
					//cout<<"φ-A I="<<I<<"width="<<width_mat[I]<<endl;
					if(conducter_flag[N[i]]==ON)
					{
						int I=npp[N[i]+nedge]+1;
						for(int j=1;j<=6;j++)
						{
							//cout<<"jloop"<<endl;
							int iside=ELEM[je].edge[j];//要素jeの辺番号
							int J=npp[iside]+1;///辺isideは行列のI番目
							//cout<<"φ-A J="<<J<<endl;
							int J1=EDGE[iside].node[1];//isideを構成する2点
							int J2=EDGE[iside].node[2];
							int k1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
							int k2=table[j][2];
							int flag=0;
							if(EDGE[iside].boundary_condition==0)///未知数
							{  	
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
										G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
										G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
										G[I][h]+=complex<double> (G_temp,0);
										flag=1;
									}
								}
								if(flag==0)
								{   
									NUM[I]=NUM[I]+1;
									
									int H=NUM[I];
									G_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
									G_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
									G_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
									ROW[I][H]=J;
									G[I][H]+=complex<double> (G_temp,0);
								}
								//B[I-1]の計算
							}
							else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int n=dn[iside];
								B_temp=co2*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								B_temp+=co2*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								B_temp+=co2*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
								B[I-1]-=complex<double> (PHAT[n]*G_temp,0);//PHATは複素ディリクレに備えなおすべき
								//B[I-1]-=co2*PHAT[n]*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*c[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*d[i];
								//B[I-1]-=co2*PHAT[n]*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*e[i];
							}//////////*/
						}
					//B[I-1]+=co*(Sx+Sy+Sz);
					}
				}

				//cout<<"φ-A eddy"<<endl;

				//φ-φ
				//double co3=sig[je]*dt*delta*delta6*delta6;
				double co3=sig[je]*delta*delta6*delta6/omega;//これに1/jがかかっているため、複素数としては符号を変えて虚数項に足されることになる(1/j=-j)
				for(int i=1;i<=4;i++)
				{			
					int I=npp[N[i]+nedge]+1;
					Sx=0;
					Sy=0;
					Sz=0;
					if(conducter_flag[N[i]]==ON)
					{
						for(int j=1;j<=4;j++)
						{
							if(conducter_flag[N[j]]==ON)
							{
								int flag=0;
								int J=npp[N[j]+nedge]+1;
								for(int h=1;h<=NUM[I];h++)
								{
									if(J==ROW[I][h])//すでに同じJが格納されているなら
									{
										G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
										G[I][h]+=complex<double>(0,-G_temp);
										flag=1;
									}
								}
								if(flag==0)//相当するＪが存在しなかったら作る
								{   
									NUM[I]=NUM[I]+1;
									int H=NUM[I];
									G_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
									G[I][H]+=complex<double>(0,-G_temp);
									ROW[I][H]=J;
								}
							}
							else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
							{
								int NN=dn[N[j]];
								B_temp=co3*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j]);
								B_N=complex<double>(0,-B_temp);
								//B[I-1]-=PHAT_V[NN]*B_N;   //ふつうの境界条件ではφが値を持たないので、ひとまず削除
								//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
								//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
								//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
							}
						}//////*/
					}
				}
				//cout<<"φ-φ eddy"<<endl;
			}
		}
	}
	}

	///////
	cout<<"G,ROW作成完了"<<endl;

	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
	cout<<"非ゼロ要素数"<<number<<endl;

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix_complex(pn,NUM,ROW,G);

	//対称性チェック
	cout<<"行列の対称性チェック----";
	check_matrix_symmetry_complex(pn,NUM,ROW,G);
	cout<<"完了"<<endl;;
	//解行列の値チェック
	for(int i=0;i<pn;i++)
	{
		//if(B[i]!=0) <<"B/=0"<<endl;
	}
	
	complex<double> *val = new complex<double> [number];
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
    /////////////////////

	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_we+mat_wn<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	//CG法実行
	complex<double> *XX=new complex<double> [pn];//行列の答え格納
    
	if(CON->get_FEMCG()==0) COCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==1) ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==2) parallel_ICCOCG(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==3) cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==4) parallel_cs_MRTR(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON->get_FEMCG()==5) cs_ICMRTR(CON,val,ind,ptr,pn,B,number,XX);

	//XXに格納された解を各辺に振る
	cout<<"複素数解の割り振りおよび出力----";
	
	/*
	double *AR = new double [nedge+1];
	double *AI = new double [nedge+1];
	double *VR = new double [node+1];
	double *VI = new double [node+1];
	*/
	for(int n=0;n<pn_e;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			AR[i]=XX[n].real();
			AI[i]=XX[n].imag();
			if(EDGE[i].stat==ON) cout<<"静的辺要素が未知数扱い？"<<endl;
		}
	}	
	for(int n=pn_e;n<pn;n++)
	{
		int i=ppn[n];
		if(i>0)
		{
			VR[i]=XX[n].real();
			VI[i]=XX[n].imag();
		}
	}	
	delete [] XX;


	//格納した値から波高、位相差を求めてベクトルポテンシャルを求める
	///Am,φのﾌｧｲﾙ出力 φは電位ではなく位相の遅れ
	//ofstream a("Am_e.dat");
	//ofstream p("phi_e.dat");
	//double Am=0.0;
	double Vm=0.0;
	//double phi=0.0;
	for(int i=1;i<=nedge;i++)
	{
		//if(EDGE[i].boundary_condition==0)
		{
			Am[i]=sqrt(AR[i]*AR[i]+AI[i]*AI[i]);
			phi[i]=atan2(AI[i],AR[i]);
			
			//if(CON->get_jw_Faverage()==ON) A[i]=Am/sqrt(2.0);
			//else A[i]=Am*cos(omega*t+phi);
			A[i]=Am[i]*cos(omega*TIME+phi[i]);
		}
		//a<<Am<<endl;
		//p<<phi<<endl;
	}
	//cout<<"ok"<<endl;
	//a.close();
	//p.close();

	
	/////渦電流解決 
	double *Je[3];//渦電流
	for(int D=0;D<3;D++) Je[D]=new double [nelm+1];
	double *Je_loss_e=new double[nelm+1];//渦電流損[W]
	double *Je_loss_n=new double[node+1];//渦電流損[W]

	calc_edge_eddy_current_jw(CON,NODE,ELEM,EDGE,node,nelm,AR,AI,dt,VR,VI,Je,Je_loss_e,t,TIME,sig,omega);
	
	
	//渦電流損を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
    {
		double dis=0;
		double dis_sum=0;
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double delta=ELEM[je].volume/6;//本当の体積

			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				dis_sum+=dis;
			}
			for(int j=1;j<=4;j++)
			{
				dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
				Je_loss_n[N[j]]+=Je_loss_e[je]*dis/dis_sum;
			}
		}
	}
	
	//節点の渦電流損を対応する粒子へ渡す
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)
		{
			PART[NODE[i].particleID].heat_generation+=Je_loss_n[i];
		}
	}

	for(int D=0;D<3;D++) delete [] Je[D];
	delete [] Je_loss_e;
	delete [] Je_loss_n;
	//
	
	//if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	//{	
	//	for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	//}
	////*/

	delete [] dn;
	delete [] PHAT;
    delete [] PHAT_A;
	delete [] PHAT_V;

    delete [] B;

	delete [] AR;
	delete [] AI;
	delete [] VR;
	delete [] VI;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
	delete [] conducter_flag;
	
	//////
	/////*/
}

//電位、磁位などのポテンシャル計算関数
void potential_calculation(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei)
{
	double *V=new double [node+1];	//potential
    
    double *Ee[3];					//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Ee[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//各要素の透磁率または誘電率格納（現在では誘電率は使っていない）

	double fluid_rp=CON->get_r_perm();		//誘電率
	if(CON->get_EM_calc_type()==4) fluid_rp=CON->get_RP();	//透磁率
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=fluid_rp;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL)
		{
			if(CON->get_EM_calc_type()==1) RP[i]=1000;//本当は導体だから∞
			else if(CON->get_EM_calc_type()==4) RP[i]=1;//透磁率
		}
		//else cout<<"error:材質に対し、誘電率あるいは透磁率が不定"<<endl;
	}

	//電位解決
	VOLT3D(CON,NODE,ELEM, node, nelm,V,jnb, TIME,PART, fluid_number,nei,RP);


	delete [] V;
    for(int D=0;D<3;D++) delete [] Ee[D];
	delete [] RP;
}

///電位計算関数
void VOLT3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpsparticle> &PART,int fluid_number,int **nei,double *RP)
{
    double V1=0;
    double V2=CON->get_V();			//電位
	double u0=0.000001256;			//真空の透磁率
	double ep0=8.854e-12;			//真空の誘電率
	double le=CON->get_distancebp();//粒子間距離
	unsigned timeA=GetTickCount();	//計算開始時刻

	////////電位決定
	if(CON->get_EM_calc_type()==1)//電場計算なら
	{
		double tau=CON->get_dt()*CON->get_V_step();//時定数
		if(CON->get_V_con()==2)
		{
			V2=CON->get_V()*(1-exp(-TIME/tau));//電位を指数関数的に増加させる
		}
		else if(CON->get_V_con()==1)
		{
			V2=(CON->get_V()-CON->get_initial_V())*TIME/(CON->get_dt()*CON->get_V_step())+CON->get_initial_V();//電位を直線的に増加させる
		}
		if(V2==0) V2=1;//0だとエラーになるから1にする
		else if(V2>CON->get_V()) V2=CON->get_V();
		cout<<"電位計算開始 V="<<V2<<" ";
	}
	else if(CON->get_EM_calc_type()==4)//磁位計算
	{
		double R=1;//比率
		V1=0;
		if(CON->get_uniform_B_sw()==OFF) V2=CON->get_magnet_B()/u0*CON->get_magnet_H();
		if(CON->get_uniform_B_sw()==ON)  V2=CON->get_uniform_B()/u0*(CON->get_ZU()-CON->get_ZD());//一様磁場
		if(CON->get_V_con()==1) 
		{
			V2=(V2)*TIME/(CON->get_dt()*CON->get_V_step());
			R=TIME/(CON->get_dt()*CON->get_V_step());
			if(V2>CON->get_magnet_B()/u0*CON->get_magnet_H()) {V2=CON->get_magnet_B()/u0*CON->get_magnet_H();R=1;}
		}
		if(V2==0) V2=1;//0だとエラーになるから1にする
		
		cout<<"磁位計算開始 V2="<<V2<<"("<<R*CON->get_magnet_B()<<"T)"<<endl;
	}
	////////////////////////////

	///磁位計算の場合、解析領域の境界条件を自由境界条件にする。（ここで、MPSTOFEMの段階でそのようにしたらFINE3Dがうまくいかない）
	if(CON->get_EM_calc_type()==4)//磁位計算
	{
		if(CON->get_uniform_B_sw()==OFF)//通常はこっち
		{
			for(int i=1;i<=nelm;i++)
			{
				if(ELEM[i].material==AIR)
				{
					for(int k=1;k<=4;k++)
					{
						int kelm=ELEM[i].elm[k];
						if(kelm==0)
						{
							int j=k%4+1;//ielmとjelmの接する三角形を構成する節点番号 
							int m=4-(k-1)/2*2;
							int n=3-(k/2%2)*2;
							NODE[ELEM[i].node[j]].boundary_condition=0;//境界条件の抹消
							NODE[ELEM[i].node[m]].boundary_condition=0;
							NODE[ELEM[i].node[n]].boundary_condition=0;
						}
					}
				}
			}
		}
		if(CON->get_uniform_B_sw()==ON)//一様磁場の場合
		{
			double ZU=CON->get_ZU();		//解析領域上端
			double ZD=CON->get_ZD();		//解析領域下端
			double err=1e-14;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].r[A_Z]>ZU-err) NODE[i].boundary_condition=2;//上端
				else if(NODE[i].r[A_Z]<ZD+err) NODE[i].boundary_condition=1;//下端
				else NODE[i].boundary_condition=0;
			}
		}
	}////////*/


    int NN=0;					//ディリクレ型境界節点数
    int *dn=new int [node+1];	//各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
    double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
		if(NODE[i].boundary_condition==1)//上で境界条件を書き換えてるから、この文はelse ifにすること
		{    
	        dn[i]=NN;
	        PHAT[NN]=V1;
	        V[i]=V1;
	        NN++;
		}
		else if(NODE[i].boundary_condition==2)
		{   
	        dn[i]=NN;
	        PHAT[NN]=V2;
	        V[i]=V2;
	        NN++;
		}
		else dn[i]=node+1;
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<" ";
	    
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
	
	/*//解行列Ｂに節点の電荷を代入
	if(CON->get_charge()==1 && CON->get_FEM_calc_type()==1)
	{
		for(int k=1;k<=WATER_N;k++)//初めのWATER_N個の節点はTRANS[k]の粒子に相当
		{
			if(NODE[k].boundary_condition==0)//未知数ならBにその節点の領域があるから代入
			{
				int i=TRANS[k];//節点kに相当する粒子番号
				for(int n=1;n<=jnb[k];n++)
				{
					int jelm=nei[k][n];//節点kの隣接する要素
					int N[5];
					for(int j=1;j<=4;j++) N[j]=ELEM[jelm].node[j];
		
					///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
					ELEM[jelm].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
					B[k]+=PART[i].angle/ep0/4*ELEM[jelm].volume;
					//if(PART[i].angle!=0)cout<<"EE"<<endl;
				}
				
			}
		}
	}///////////*/
	
	double *charge=new double [nelm+1];
	for(int n=1;n<=nelm;n++) charge[n]=0;
	
    /////////全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
		double delta6=ELEM[je].volume;//体積の6倍
		
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
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
	
		/////比誘電率を定義
		double ep=RP[je];
		/*if(ELEM[je].material==FLUID)
		{
			if(CON->get_EM_calc_type()==1) ep=CON->get_r_perm();
			else if(CON->get_EM_calc_type()==4) ep=CON->get_RP();//磁位計算だから比透磁率
		}*/
	
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
								G[I][h]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
								//if(I==85839 && h==1) cout<<h<<" "<<G[I][h]<<" "<<delta<<" "<<(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])<<endl;
								flag=1;
							}
						}
						if(flag==0)
						{   
							NUM[I]=NUM[I]+1;
							int H=NUM[I];
			    
							G[I][H]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
							ROW[I][H]=J;
						}
					}
					else //N[j]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
						int n=dn[N[j]];
						B[I-1]-=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta*PHAT[n];
					}
				}
				//if(CON->get_charge()==ON) B[I-1]+=charge[je]/4*delta;
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
	if(CON->get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//else if(CON->get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		V[i]=XX[n];//電位φ
	}
	delete [] XX;
	////////////////////////////
    
    ///////電位をﾌｧｲﾙ出力
	ofstream fp("V.dat");
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].r[A_Y]<2*le && NODE[i].r[A_Y]>-2*le) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<V[i]<<endl;
    }
	fp.close(); 
    
    ////////////////////////

    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] charge;
}

//行列の幅計算関数(新)
int calc_matrix_width(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei)
{
	//係数行列のなかで、i行目に含まれる非ぜろ要素数は、計算点iからのびる辺の数(すなわち近隣計算点数)+1(自身)である。これの最大値を求める
	///配列確保
		
	int *temp_check=new int[node+1];	//一時的なﾁｪｯｸ配列
		
	///初期化
	for(int i=1;i<=node;i++) temp_check[i]=0;
	/////////////////////

	int maxwidth=0;

	//if(CON->get_FEM_calc_type()==1 || CON->get_FEM_calc_type()==4 || CON->get_iteration_count()==0)
	{
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].boundary_condition==0)
			{
				int width=0;
				vector<int> NEI2;//各節点の隣接する節点番号格納(neiは節点-要素であるが、NEI2は節点-節点)
				for(int k=1;k<=jnb[i];k++)
				{
					int jelm=nei[i][k];//節点iに隣接する要素番号
					for(int j=1;j<=4;j++)
					{
						int p=ELEM[jelm].node[j];
						if(p!=i && temp_check[p]==0 && NODE[p].boundary_condition==0)
						{	
							width++;
							NEI2.push_back(p);
							temp_check[p]=1;//もう見た
						}
					}
				}
				for(int k=0;k<NEI2.size();k++) temp_check[NEI2[k]]=0;//初期化
				if(width>maxwidth) maxwidth=width;
			}
		}
	}

	delete [] temp_check;

	return maxwidth+1;//widthは自分から伸びる辺の数に等しい。行列の幅は自分自身もいれてwidth+1
}

int calc_matrix_width2(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,int *width)
{
	//係数行列のなかで、i行目に含まれる非ぜろ要素数は、計算点iからのびる辺の数(すなわち近隣計算点数)+1(自身)である。これの最大値を求める
	///配列確保
		
	int *temp_check=new int[node+1];	//一時的なﾁｪｯｸ配列
		
	///初期化
	for(int i=1;i<=node;i++) temp_check[i]=0;
	/////////////////////

	int maxwidth=0;

	//if(CON->get_FEM_calc_type()==1 || CON->get_FEM_calc_type()==4 || CON->get_iteration_count()==0)
	{
		for(int i=1;i<=node;i++)
		{
			width[i]=0;//初期化
			
			if(NODE[i].boundary_condition==0)
			{
				int width_temp=0;
				vector<int> NEI2;//各節点の隣接する節点番号格納(neiは節点-要素であるが、NEI2は節点-節点)
				for(int k=1;k<=jnb[i];k++)
				{
					int jelm=nei[i][k];//節点iに隣接する要素番号
					for(int j=1;j<=4;j++)
					{
						int p=ELEM[jelm].node[j];
						if(p!=i && temp_check[p]==0 && NODE[p].boundary_condition==0)
						{	
							width_temp++;
							NEI2.push_back(p);
							temp_check[p]=1;//もう見た
						}
					}
				}
				for(int k=0;k<NEI2.size();k++) temp_check[NEI2[k]]=0;//初期化
				if(width_temp>maxwidth) maxwidth=width_temp;
				
				width[i]=width_temp+1;
			}
			
			
		}
	}

	delete [] temp_check;

	return maxwidth+1;//widthは自分から伸びる辺の数に等しい。行列の幅は自分自身もいれてwidth+1
}

//行列並び替え関数
void arrange_matrix(int pn,int *NUM,int **ROW,double **G)
{
	///G[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
			for(int m=1;m<j;m++)
			{
				if(ROW[i][j]<ROW[i][m])
				{
					tempG=G[i][m];
					tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////
}
//並び替え、複素数
void arrange_matrix_complex(int pn,int *NUM,int **ROW,complex<double>**G)
{
	///G[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	
    for(int i=1;i<=pn;i++)
    {
        complex<double> tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
			for(int m=1;m<j;m++)
			{
				if(ROW[i][j]<ROW[i][m])
				{
					tempG=G[i][m];
					tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////
}

///ICCG法
void ICCG3D2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	int count;						//数え上げ変数
	double accel=CON->get_CGaccl();	//加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *matrix=new double[num2];//係数行列を保存(対角成分を含む下三角行列) 非ゼロ要素だけを1次配列として保存
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
	    {
			if(ind[m]<=k)
			{
				matrix[num2]=val[m];
				ind2[num2]=ind[m];
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];			//列方向にみた、各列の要素数
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];						//ゼロ要素の値 VAL[i][k]はi列のk番目の非ゼロ要素
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];							//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	double alp,beta;
	double rLDLt_r;
	double E=1;//誤差
    double *r=new double[pn];
	double *AP = new double [pn];
	double *P = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *L=new double[num2];//不完全これスキー分解後の下三角行列L格納
	double *D1 = new double [pn];//D行列
	
	/////不完全コレスキｰ分解
	Incomplete_Cholesky_Decomposition(CON,L,D1,matrix,ptr2,ind2,pn,num2);//LとD1に値が書き込まれる
	
	delete [] matrix;

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=L[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0;
		r[n]=B[n];
	}

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/L[0]; //式（3.77） 
		else
		{
		    double sum=0;
			for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=L[m]*y[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
			y[i]=(r[i]-sum)/L[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

    cout<<"ICCG法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	count=0;
	double ep=CON->get_FEMCGep();//収束判定
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alpを求める
		rLDLt_r=0;

		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
		}
		double PAP=0;
		for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		//cout<<"alp="<<alp<<" PAP="<<PAP<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		for(int n=0;n<pn;n++) X[n]+=alp*P[n];
		
		//////////////// r=r-alp*AP
		for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		//////////////////誤差
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];
		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/L[0]; //式（3.77） 新
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//対角成分は除くからptr[i+1]-1
			    {
					sum+=L[m]*y[ind2[m]];//式（3.78）
			    }
			    int m=ptr2[i+1]-1;
				y[i]=(r[i]-sum)/L[m];
			}
		}////y[i]がもとまった。
	
		/////////LDLt_r[i]を求める
		for(int i=pn-1;i>=0;i--)
		{
		    double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;

	delete [] L;
}

//COCG法　複素数用
void COCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	complex<double> rr=(0,0);
	double E=1;//誤差
	double BB=0;
	complex<double> alp;
	complex<double> beta;
    complex<double> *r=new complex<double>[pn];
	complex<double> *AP = new complex<double> [pn];
	complex<double> *P = new complex<double> [pn];

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
	}
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;

	for(int n=0;n<pn;n++) BB= BB+norm(B[n]);
	BB=sqrt(BB);
	//cout<<"BB="<<BB<<endl;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];
	//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
	
	ofstream c("convergence.dat");

	 cout<<"COCG法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;
			//cout<<count<<"回目"<<endl;
			complex<double> rr_be;
			rr_be=rr;
			//cout<<"rr="<<rr<<endl;
			//////////////alpを求める
			for(int n=0;n<pn;n++)
			{    
				AP[n]=complex<double>(0.0,0.0);
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
			}
			//cout<<"AP[0]="<<AP[0]<<" AP[1]="<<AP[1]<<endl;//AP[0]=4+i
			complex<double> PAP=(0.0,0.0);
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
			//for(int n=0;n<pn;n++)  PAP+=conj(P[n])*conj(AP[n]);
			alp=rr/PAP;
			//cout<<"alp="<<alp<<" PAP="<<PAP<<endl;
			//////////////////////
		
			//////////////// 解更新　X(k+1)=X(k)+alp*P
			for(int n=0;n<pn;n++) X[n]+=alp*P[n];
			//cout<<"X[0]="<<X[0]<<" X[1]="<<X[1]<<endl;
			//////////////////////////////
			
			//////////////// r=r-alp*AP
			for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
			/////////////////////////////
			
			///////////////////////beta
			rr=complex<double> (0.0,0.0);
			//#pragma omp parallel for
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;
			//cout<<"rrk+1="<<rr<<endl;

			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			beta=rr/rr_be;
			//cout<<"beta="<<beta<<endl;
			///////////////////////

			//////////////////誤差
			//#pragma omp parallel for reduction(+:E)
			E=0.0;
			for(int n=0;n<pn;n++) E+=norm(r[n]);
			E=sqrt(E);
			E/=BB;
			c<<count<<" "<<E<<endl;
			//cout<<"E="<<E<<endl;
			if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
			////////////////////////
			
			///////////////////// P=r+beta*P
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];
			//cout<<"P[0]="<<P[0]<<" P[1]="<<P[1]<<endl;
		}
	}

	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] AP;
	delete [] P;
	c.close();
}

//前処理つきCOCG
void ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	//int count=0;
	complex<double> rr;
	rr=complex<double> (0,0);
	double E=1;//誤差
	double BB=0;
	complex<double> alp;
	complex<double> beta;
    complex<double> *r=new complex<double>[pn];
	complex<double> *AP = new complex<double> [pn];
	complex<double> *P = new complex<double> [pn];

	//収束過程を見たい
	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=complex<double> (0,0);
		r[n]=B[n];
		P[n]=r[n];
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	BB=sqrt(BB);
	cout<<"BB="<<BB<<endl;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];

	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	double accel_re=CON->get_CGaccl();
	complex<double> accel;//CON->get_CGaccl();//加速ファクタ
	accel=complex<double> (accel_re,0);
	double accel2=0.001;//複素シフト用加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	complex<double> *val2=new complex<double> [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;

	complex<double> one;
	one=complex<double> (1.0,0);
	complex<double> Im;
	Im=complex<double> (0.0,1.0);
	complex<double> sum;
	sum=complex<double> (0.0,0.0);

	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//加速ﾌｧｸﾀ
				//if(ind[m]==k) val2[num2]+=accel2*Im;//複素シフト
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	complex<double> **VAL=new complex<double> *[pn];//ゼロ要素の値
	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
	int **IND = new int *[pn];//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	complex<double> rLDLt_r;
	complex<double> *y=new complex<double> [pn];
	complex<double> *LDLt_r= new complex<double> [pn];
	complex<double> *D1 = new complex<double> [pn];//D行列
	
	/////不完全コレスキｰ分解
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
	        int i=ind2[m];//列番号
	        if(i==0)
			{
				val2[m]=val2[m];
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				sum=complex<double> (0,0);
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i行目のなかから列の一致するものを探している。少し手間か？
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=complex<double> (0,0);
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
				
				D1[k]=one/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
            }
	    }
	}    
	///不完全ｺﾚｽｷｰ分解完了/////////*/

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 
		else
		{
		    sum=complex<double> (0,0);
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
	    sum=complex<double> (0,0);
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

	cout<<"ICCOCG法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	//cout<<"ICCG法:未知数="<<pn<<" ---";
	//unsigned int time=GetTickCount();
	int count=0;
	rLDLt_r=complex<double> (0,0);
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//最初のrLDLt_rだけここで求める
	//while(E>CON->get_FEMCGep())
	double ep=CON->get_FEMCGep();
	while(E>ep && count<=40000)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		if(count%10000==0) ep*=10;
		//////////////alpを求める
		complex<double> PAP;
		PAP=complex<double> (0,0);
		//#pragma omp parallel for reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{
			AP[n]=complex<double>(0,0);
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		//cout<<"alp="<<alp<<endl;
		//////////////////////
		E=0;
		//#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++)
		{
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
			r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
			E+=norm(r[n]);						//更新後の誤差
		}
		E=sqrt(E);
		E/=BB;
		//cout<<"E="<<E<<endl;
		c<<count<<" "<<E<<endl;
		if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
		if(E<CON->get_FEMCGep()) cout<<"E<ε E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		
		beta=one/rLDLt_r;
		rLDLt_r=complex<double>(0,0);
		
        /////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 新
			else
			{
			   sum= complex<double> (0,0);
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//対角成分は除くからptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//式（3.78）
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]がもとまった。
	
		/////////LDLt_r[i]を求める
		for(int i=pn-1;i>=0;i--)
		{
		   sum= complex<double> (0,0);
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		//for(int n=0;n<pn;n++) rLDLt_r+=conj(r[n])*conj(LDLt_r[n]);
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}

	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	c.close();

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	
	delete [] r;
	delete [] AP;
	delete [] P;
}

//openMPによって並列化したiccocg
void parallel_ICCOCG(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	//int count=0;
	complex<double> rr;
	rr=complex<double> (0,0);
	double E=1;//誤差
	double BB=0;
	complex<double> alp;
	complex<double> beta;
    complex<double> *r=new complex<double>[pn];
	complex<double> *AP = new complex<double> [pn];
	complex<double> *P = new complex<double> [pn];
	complex<double> *sum3 = new complex<double> [pn];

	//収束過程を見たい
	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=complex<double> (0,0);
		r[n]=B[n];
		P[n]=r[n];
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	BB=sqrt(BB);
	cout<<"BB="<<BB<<endl;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];

	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	double accel_re=CON->get_CGaccl();
	complex<double> accel;//CON->get_CGaccl();//加速ファクタ
	accel=complex<double> (accel_re,0);
	double accel2=0.001;//複素シフト用加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	complex<double> *val2=new complex<double> [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;

	complex<double> one;
	one=complex<double> (1.0,0);
	complex<double> Im;
	Im=complex<double> (0.0,1.0);
	complex<double> sum;
	sum=complex<double> (0.0,0.0);

	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//加速ﾌｧｸﾀ
				//if(ind[m]==k) val2[num2]+=accel2*Im;//複素シフト
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	complex<double> **VAL=new complex<double> *[pn];//ゼロ要素の値
	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
	int **IND = new int *[pn];//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	complex<double> rLDLt_r;
	complex<double> *y=new complex<double> [pn];
	complex<double> *LDLt_r= new complex<double> [pn];
	complex<double> *D1 = new complex<double> [pn];//D行列
	
	/////不完全コレスキｰ分解
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
	        int i=ind2[m];//列番号
	        if(i==0)
			{
				val2[m]=val2[m];
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				sum=complex<double> (0,0);
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i行目のなかから列の一致するものを探している。少し手間か？
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=complex<double> (0,0);
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
				
				D1[k]=one/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
            }
	    }
	}    
	///不完全ｺﾚｽｷｰ分解完了/////////*/

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 
		else
		{
		    sum=complex<double> (0,0);
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
	    sum=complex<double> (0,0);
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

	cout<<"parallel_ICCOCG法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	//cout<<"ICCG法:未知数="<<pn<<" ---";
	//unsigned int time=GetTickCount();
	int count=0;
	rLDLt_r=complex<double> (0,0);
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//最初のrLDLt_rだけここで求める
	//while(E>CON->get_FEMCGep())
	double ep=CON->get_FEMCGep();

	//#ifdef _OPENMP
	//#pragma omp parallel
	//#endif
	while(E>ep && count<=40000)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		if(count%10000==0) ep*=10;
		//////////////alpを求める
		complex<double> PAP=0.;

		#pragma omp parallel shared(PAP)
		{
			complex< double > priv_PAP=0.;
			
			#pragma omp for
			for(int n=0;n<pn;n++)
			{
				AP[n]=complex<double>(0,0);
				for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
				priv_PAP+=P[n]*AP[n];
			}
			#pragma omp critical
			{
			  PAP += priv_PAP;
			}
		}
		#pragma omp parallel// private(sum) 
		{
			//#pragma omp for reduction(+:PAP)
			//for(int n=0;n<pn;n++) PAP+=P[n]*AP[n];
			
			#pragma omp single
			{
				alp=rLDLt_r/PAP;
				E=0;
			}
			
			#pragma omp for reduction(+:E)
			for(int n=0;n<pn;n++)
			{
				X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
				r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
				//E+=norm(r[n]);						//更新後の誤差
				E+=abs(r[n])*abs(r[n]);
			}
		//}
			//#pragma omp for reduction(+:E)
			//for(int n=0;n<pn;n++) E+=norm(r[n]);						//更新後の誤差
			#pragma omp single
			{
				E=sqrt(E);
				E/=BB;	
				//cout<<"E="<<E<<endl;
				c<<count<<" "<<E<<endl;
				if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
				if(E<CON->get_FEMCGep()) cout<<"E<ε E="<<E<<endl;
				///////////////////////beta
				beta=one/rLDLt_r;
				rLDLt_r=complex<double>(0,0);

				y[0]=r[0]/val2[0];
			}
		}	
			/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
			//#pragma omp parallel for 
			for(int i=1;i<pn;i++)
			{
				sum3[i]=complex<double>(0,0);
			    //complex<double> priv_sum1=0.; 
				//#pragma omp parallel for 
				for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//対角成分は除くからptr[i+1]-1
				{
					sum3[i]+=val2[m]*y[ind2[m]];//式（3.78）
				}
				y[i]=(r[i]-sum3[i])/val2[ptr2[i+1]-1];
				//int m2=ptr2[i+1]-1;
				//y[i]=(r[i]-sum3[i])/val2[ptr2[i+1]-1];
			}////y[i]がもとまった。
			
			//#pragma omp parallel for 

		//}
			/////////LDLt_r[i]を求める
		//	#pragma omp for 
			for(int i=pn-1;i>=0;i--)
			{
				sum= complex<double> (0,0);
				//#pragma omp parallel for   
				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
				LDLt_r[i]=y[i]-D1[i]*sum;	
			}
			/////////////////*/
		//}
		
			//for(int n=0;n<pn;n++) rLDLt_r+=conj(r[n])*conj(LDLt_r[n]);
		#pragma omp parallel shared(rLDLt_r)
		{
			complex< double > priv_rLDLt_r=0.;
			
			#pragma omp for
			for(int n=0;n<pn;n++)
			{
				priv_rLDLt_r+=r[n]*LDLt_r[n];
			}
			#pragma omp critical
			{
			  rLDLt_r += priv_rLDLt_r;
			}
		}

		beta=beta*rLDLt_r;
			/////////////////*/
			
			///////////////////// P=r+beta*P
		#pragma omp parallel for
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}

	

	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	c.close();

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	
	delete [] r;
	delete [] AP;
	delete [] P;
	delete [] sum3;
}

//////mps115-1-test1からそのまま移植した並列iccg
void parallel_ICCG3D2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	//収束過程を見たい
	ofstream c("convergence.dat");

	double accel=CON->get_CGaccl();//加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *val2=new double [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	double BB=0;
	for(int k=0;k<pn;k++)
	{
		BB+=B[k]*B[k];
	}
	BB=sqrt(BB);

	cout<<"BB="<<BB<<endl;

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//加速ﾌｧｸﾀ
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];//ゼロ要素の値
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	double alp,beta;
	double rLDLt_r;
	double E=BB;//誤差
    double *r=new double[pn];
	for(int n=0;n<pn;n++) X[n]=0;
	
	double *AP = new double [pn];
	double *P = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *D1 = new double [pn];//D行列
	
	/////不完全コレスキｰ分解
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
	        int i=ind2[m];//列番号
	        if(i==0)
			{
				val2[m]=val2[m];
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				double sum=0;
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i行目のなかから列の一致するものを探している。少し手間か？
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				D1[k]=1/val2[m];
            }
	    }
	}    
	///不完全ｺﾚｽｷｰ分解完了/////////*/

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0;
		r[n]=B[n];
	}

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 
		else
		{
		    double sum=0;
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

    cout<<"parallel_ICCG法スタート  -----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep;
	ep=CON->get_FEMCGep();//収束判定
	while(E>ep && count<=40000)
	{
		count++;
		if(count==pn) cout<<"count>pn E="<<E<<endl;
		if(count==10000) cout<<"count=10000 E="<<E<<endl;
		if(count%10000==0) ep*=10;
		//////////////alpを求める
		rLDLt_r=0;
		double PAP=0;
		//for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//sirial
		#pragma omp parallel for reduction(+:rLDLt_r) reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{	
			//printf("%d\n",omp_get_thread_num());
			rLDLt_r+=r[n]*LDLt_r[n];
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];//sirial
		alp=rLDLt_r/PAP;
		
//		cout<<"alp="<<alp<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		E=0;
		#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++) 
		{
			X[n]+=alp*P[n];	// X(k+1)=X(k)+alp*P
			r[n]-=alp*AP[n];// r=r-alp*AP
			E+=r[n]*r[n];	//誤差
		}
		E=sqrt(E)/BB;
		if(count%10==0) c<<count<<" "<<E<<endl;
		//c<<count<<" "<<E<<endl;
		if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
		//////////////////////////////
		
		//////////////// r=r-alp*AP
		//for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		/*/////////////////誤差
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		///////////////////////*/
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 新
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//対角成分は除くからptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//式（3.78）
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]がもとまった。
	
		/////////LDLt_r[i]を求める
		for(int i=pn-1;i>=0;i--)
		{
		    double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
	}
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	c.close();
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
}

//MRTR法　複素数用// "MRTR法の複素対称線型方程式への拡張(2007)"
void cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//誤差
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> v;

	complex<double> cAr_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(複素共役をバーで表記)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	v=complex<double>(1.0,0.0);//1回目の反復において、0でないことに注意

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
	}

	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
	//for(int n=0;n<pn;n++) BB+=norm(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;

	/*/////1ステップ目の計算///////////////

	for(int n=0;n<pn;n++)//Ar,Ar_c計算
	{    
		Ar[n]=complex<double>(0.0,0.0);
		cAr[n]=complex<double>(0.0,0.0);
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			Ar[n]+=val[j]*r[ind[j]];
			cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
		}
	}

	cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	
	//#pragma omp parallel for
	for(int n=0;n<pn;n++)  
	{
		cAr_r+=conj(cAr[n])*r[n];
		cAr_Ar+=conj(cAr[n])*Ar[n];
	}

	zeta=cAr_r/cAr_Ar;
	eta=complex<double>(0.0,0.0);

	//////vの計算//////////////
	for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

	//////pの計算//////////////
	for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
	zeta_old=zeta;

	///////Xの計算////////////
	for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
	//////yの計算////////////
	for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

	//////rの計算////////////
	for(int n=0;n<pn;n++) r[n]-=y[n];
	///////*/


	
	//for(int n=0;n<pn;n++) rr+=norm(r[n]);
	//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
	//E=sqrt(rr/BB);

	//cout<<"反復回数="<<1<<" E1="<<E<<endl;
	//count++;
	

	 cout<<"cs_MRTR法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			for(int n=0;n<pn;n++)//Ar,Ar_c計算
			{    
				Ar[n]=complex<double>(0.0,0.0);
				cAr[n]=complex<double>(0.0,0.0);
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
					cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
				}
			}

			cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
			cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			cy_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			cAr_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				cAr_r+=conj(cAr[n])*r[n];
				cAr_Ar+=conj(cAr[n])*Ar[n];
				cy_Ar+=y[n]*Ar[n];
				cAr_y+=conj(cAr[n])*y[n];
			}

			zeta=v*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);
			eta=-cy_Ar*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);

			//////vの計算//////////////
			v=complex<double>(0.0,0.0);
			for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

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
			//for(int n=0;n<pn;n++) rr+=norm(r[n]);
			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	c.close();
}

//対角スケーリングつきMRTR法　複素数用// "MRTR法の複素対称線型方程式への拡張(2007)"
void DS_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//誤差
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> v;

	complex<double> cAr_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(複素共役をバーで表記)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	v=complex<double>(1.0,0.0);//1回目の反復において、0でないことに注意

	//対角スケーリング
	for(int n=0;n<pn;n++)
	{    
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			val[j]=val[j]*r[ind[j]];
		}
	}


	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;

	/*/////1ステップ目の計算///////////////

	for(int n=0;n<pn;n++)//Ar,Ar_c計算
	{    
		Ar[n]=complex<double>(0.0,0.0);
		cAr[n]=complex<double>(0.0,0.0);
		for(int j=ptr[n];j<ptr[n+1];j++)
		{
			Ar[n]+=val[j]*r[ind[j]];
			cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
		}
	}

	cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	
	//#pragma omp parallel for
	for(int n=0;n<pn;n++)  
	{
		cAr_r+=conj(cAr[n])*r[n];
		cAr_Ar+=conj(cAr[n])*Ar[n];
	}

	zeta=cAr_r/cAr_Ar;
	eta=complex<double>(0.0,0.0);

	//////vの計算//////////////
	for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

	//////pの計算//////////////
	for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
	zeta_old=zeta;

	///////Xの計算////////////
	for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
	//////yの計算////////////
	for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];

	//////rの計算////////////
	for(int n=0;n<pn;n++) r[n]-=y[n];
	///////*/


	
	//for(int n=0;n<pn;n++) rr+=norm(r[n]);
	//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
	//E=sqrt(rr/BB);

	//cout<<"反復回数="<<1<<" E1="<<E<<endl;
	//count++;
	

	 cout<<"cs_MRTR法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			for(int n=0;n<pn;n++)//Ar,Ar_c計算
			{    
				Ar[n]=complex<double>(0.0,0.0);
				cAr[n]=complex<double>(0.0,0.0);
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
					cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
				}
			}

			cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
			cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			cy_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			cAr_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				cAr_r+=conj(cAr[n])*r[n];
				cAr_Ar+=conj(cAr[n])*Ar[n];
				cy_Ar+=y[n]*Ar[n];
				cAr_y+=conj(cAr[n])*y[n];
			}

			zeta=v*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);
			eta=-cy_Ar*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);

			//////vの計算//////////////
			v=complex<double>(0.0,0.0);
			for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];

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
			for(int n=0;n<pn;n++) rr+=norm(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	c.close();
}

//並列MRTR法　複素数用// "MRTR法の複素対称線型方程式への拡張(2007)"
void parallel_cs_MRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//誤差
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> v;

	complex<double> cAr_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(複素共役をバーで表記)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	v=complex<double>(1.0,0.0);//1回目の反復において、0でないことに注意

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
	}

	for(int n=0;n<pn;n++) BB+=norm(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;
	
	/*////
	#pragma omp parallel shared(PAP)
		{
			complex< double > priv_PAP=0.;
			
			#pragma omp for
			for(int n=0;n<pn;n++)
			{
				AP[n]=complex<double>(0,0);
				for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
				priv_PAP+=P[n]*AP[n];
			}
			#pragma omp critical
			{
			  PAP += priv_PAP;
			}
		}
		#pragma omp parallel// private(sum) 
		{
			//#pragma omp for reduction(+:PAP)
			//for(int n=0;n<pn;n++) PAP+=P[n]*AP[n];
			
			#pragma omp single
			{
				alp=rLDLt_r/PAP;
				E=0;
			}
			
			#pragma omp for reduction(+:E)
			for(int n=0;n<pn;n++)
			{
				X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
				r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
				E+=norm(r[n]);						//更新後の誤差
			}
		//}
			//#pragma omp for reduction(+:E)
			//for(int n=0;n<pn;n++) E+=norm(r[n]);						//更新後の誤差
			#pragma omp single
	////////*/
	cout<<"parallel_cs_MRTR法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			#pragma omp parallel 
			{
				#pragma omp for
				for(int n=0;n<pn;n++)//Ar,Ar_c計算
				{    
					Ar[n]=complex<double>(0.0,0.0);
					cAr[n]=complex<double>(0.0,0.0);
					for(int j=ptr[n];j<ptr[n+1];j++)
					{
						Ar[n]+=val[j]*r[ind[j]];
						cAr[n]+=conj(val[j]) * conj(r[ind[j]]);
					}
				}
			
				#pragma omp single
				{
					cAr_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
					cAr_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
					cy_Ar=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
					cAr_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
				}
			}
			//#pragma omp parallel for
			#pragma omp parallel shared(cAr_r,cAr_Ar,cy_Ar,cAr_y)
			{
				#pragma omp for
				for(int n=0;n<pn;n++)  
				{
					cAr_r+=conj(cAr[n])*r[n];
					cAr_Ar+=conj(cAr[n])*Ar[n];
					cy_Ar+=y[n]*Ar[n];
					cAr_y+=conj(cAr[n])*y[n];
				}
			}

			zeta=v*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);
			eta=-cy_Ar*cAr_r/(v*cAr_Ar-cy_Ar*cAr_y);

			//////vの計算//////////////
			v=complex<double>(0.0,0.0);

			#pragma omp parallel shared(v)
			{
				#pragma omp for
				for(int n=0;n<pn;n++) v+=zeta*r[n]*Ar[n];
			}

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
			for(int n=0;n<pn;n++) rr+=norm(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	c.close();
}

//ic付きMRTR法
void cs_ICMRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	//complex<double> rr=(0,0);
	double E=1;//誤差
	double BB=0;
	double rr=0;
	
	complex<double> alp;
	complex<double> beta;
	complex<double> zeta;
	complex<double> zeta_old;
	complex<double> eta;
	complex<double> nu;

	complex<double> cAr_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

	complex<double> cw_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cv_w;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cy_w;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	complex<double> cw_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

    complex<double> *r= new complex<double>[pn];
	complex<double> *Ar = new complex<double> [pn];  //   _ _
	complex<double> *cAr = new complex<double> [pn];//cAr:A r(複素共役をバーで表記)
	complex<double> *P = new complex<double> [pn];
	complex<double> *y = new complex<double> [pn];
	complex<double> *u = new complex<double> [pn];
	complex<double> *v = new complex<double> [pn];
	complex<double> *w = new complex<double> [pn];
	complex<double> *z = new complex<double> [pn];

	zeta=complex<double>(0.0,0.0);
	zeta_old=complex<double>(0.0,0.0);
	eta=complex<double>(0.0,0.0);
	nu=complex<double>(1.0,0.0);//1回目の反復において、0でないことに注意

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=complex<double> (0.0,0.0);
		r[n]=B[n];
		P[n]=r[n];
		y[n]=complex<double> (0.0,0.0);
		u[n]=complex<double> (0.0,0.0);
		v[n]=complex<double> (0.0,0.0);
		w[n]=complex<double> (0.0,0.0);
		z[n]=complex<double> (0.0,0.0);
	}

	//////前処理/////
	double accel_re=CON->get_CGaccl();
	complex<double> accel;//CON->get_CGaccl();//加速ファクタ
	accel=complex<double> (accel_re,0);
	double accel2=0.001;//複素シフト用加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	complex<double> *val2=new complex<double> [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;

	complex<double> one;
	one=complex<double> (1.0,0);
	complex<double> Im;
	Im=complex<double> (0.0,1.0);
	complex<double> sum;
	sum=complex<double> (0.0,0.0);

	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//加速ﾌｧｸﾀ
				//if(ind[m]==k) val2[num2]+=accel2*Im;//複素シフト
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	complex<double> **VAL=new complex<double> *[pn];//ゼロ要素の値
	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
	int **IND = new int *[pn];//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	complex<double> rLDLt_r;
	complex<double> *y2=new complex<double> [pn];
	complex<double> *LDLt_r= new complex<double> [pn];
	complex<double> *D1 = new complex<double> [pn];//D行列
	
	/////不完全コレスキｰ分解
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
	        int i=ind2[m];//列番号
	        if(i==0)
			{
				val2[m]=val2[m];
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				sum=complex<double> (0,0);
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i行目のなかから列の一致するものを探している。少し手間か？
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=complex<double> (0,0);
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
				
				D1[k]=one/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
            }
	    }
	}    
	///不完全ｺﾚｽｷｰ分解完了/////////*/

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y2[0]=r[0]/val2[0]; //式（3.77） 
		else
		{
		    sum=complex<double> (0,0);
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
		    y2[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
	    sum=complex<double> (0,0);
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y2[i]-D1[i]*sum;	
	}
	/////////////////*///LDLtがM-1に相当？つまり、LDLt_rはu0と等しい
	
	for(int n=0;n<pn;n++) u[n]=LDLt_r[n];
		/////////////////*/
	//////////////////////

	//for(int n=0;n<pn;n++) BB+=norm(B[n]);
	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
	//BB=sqrt(BB);
	cout<<"||r0||2="<<sqrt(BB)<<endl;


	 cout<<"cs_ICMRTR法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON->get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			count++;

			/////vの計算
			for(int n=0;n<pn;n++)
			{    
				v[n]=complex<double>(0.0,0.0);//vに相当
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					v[n]+=val[j]*u[ind[j]];
				}
			}

			/*///
			for(int n=0;n<pn;n++)//Ar,Ar_c計算
			{    
				Ar[n]=complex<double>(0.0,0.0);//vに相当
				cAr[n]=complex<double>(0.0,0.0);//conj(v)に相当
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*u[ind[j]];
					cAr[n]+=conj(val[j]) * conj(u[ind[j]]);
				}
			}
			///*/

			////wの計算

			/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
			for(int i=0;i<pn;i++)
			{
				if(i==0) y2[0]=v[0]/val2[0]; //式（3.77） 
				else
				{
					sum=complex<double> (0,0);
					/////////        
					for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//式（3.78）
					int m=ptr2[i+1]-1;
					y2[i]=(v[i]-sum)/val2[m];
				}
			}////y[i]がもとまった。
			for(int i=pn-1;i>=0;i--)
			{
				sum=complex<double> (0,0);
				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
				LDLt_r[i]=y2[i]-D1[i]*sum;	
			}
			/////////////////*///LDLtがM-1に相当？つまり、LDLt_rはu0と等しい 同様に、LDL_vはw？
	
			for(int n=0;n<pn;n++) w[n]=LDLt_r[n];

			cw_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
			cv_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			cy_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			cw_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	
			//内積の計算
			for(int n=0;n<pn;n++)  
			{
				cw_r+=w[n]*r[n];
				cv_w+=v[n]*w[n];
				cy_w+=y[n]*w[n];
				cw_y+=w[n]*y[n];
			}

			zeta=nu*cw_r/(nu*cv_w-cy_w*cw_y);
			eta=-cy_w*cw_r/(nu*cv_w-cy_w*cw_y);

			//////nuの計算//////////////
			nu=complex<double>(0.0,0.0);
			for(int n=0;n<pn;n++) nu+=zeta*r[n]*w[n];

			//////pの計算//////////////
			for(int n=0;n<pn;n++) P[n]=u[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////Xの計算////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////yの計算////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*v[n];

			//////rの計算////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//////zの計算////////////
			for(int n=0;n<pn;n++) z[n]=eta*z[n]+zeta*w[n];

			//////uの計算////////////
			for(int n=0;n<pn;n++) u[n]-=z[n];

			//誤差評価
			rr=0;
			//for(int n=0;n<pn;n++) rr+=norm(r[n]);
			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
	
			
		}
	}

	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	delete [] y2;

	delete [] r;
	delete [] Ar;
	delete [] cAr;	
	delete [] P;
	delete [] y;
	delete [] u;
	delete [] v;
	delete [] w;
	delete [] z;

	c.close();
}

//不完全ｺﾚｽｷｰ分解
void Incomplete_Cholesky_Decomposition(mpsconfig *CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2)
{
	//matrix[i]:係数行列のうち、対角を含む下三角が格納されている
	//ptr2[k]:matrix[i]に対するｱｸｾｽ配列
	//L[i],D[i]:不完全ｺﾚｽｷｰ分解格納配列
	//pn:未知数
	//num2:matrix[i]の要素数

	//cout<<"不完全ｺﾚｽｷｰ分解開始---";
	unsigned int timeC=GetTickCount();

	int UP=1;							//
	int DOWN=2;
	int flag=DOWN;						//flag==UPなら、加速係数を1より大きくしていく。DOWNなら小さくしていく
	int loopsw;							//while分のスイッチ

	double accel=1;//CON->get_CGaccl();		//加速係数
	double maxvalue=0;					//matrix中の最大値
	int Line=1;							//L行列のなかで最大値を持つ行番号
	int Column=1;						//L行列のなかで最大値を持つ列番号
	
	
	//1回目:まずは普通に不完全ｺﾚｽｷｰ分解を行う()。そこでD1[i]<0なら加速係数を増加する(UP)。D1[i]>0なら加速係数を小さくする(DOWN)
	L[0]=matrix[0]*accel;
	D1[0]=1/L[0];			//L[0]は係数行列の値と変わらないので、ゼロでもなければ負でもないと予想される
	if(D1[0]<0) cout<<"error in 不完全ｺﾚｽｷｰ分解 D1[0]<0"<<endl;
	maxvalue=D1[0]*D1[0];
	for(int k=1;k<pn;k++)	//k=0はもう済んだ。k>=1からloop
	{	
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0からk-1までの要素
		{
			int i=ind2[m];//列番号
		
			double sum=0;
			for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i行目のなかから列の一致するものを探している。少し手間か？
			L[m]=matrix[m]-sum;//i=0のときはsum=0となり自動的にL[m]=matrix[m]となる
			if(L[m]*L[m]>maxvalue)
			{
				maxvalue=L[m]*L[m];
				Line=k; Column=i;
			}
		}
		int m=ptr2[k+1]-1;				//このmは必ず対角成分に該当する
				
		double akk=matrix[m];			//係数行列の対角成分.値を保存.
		double sum=0;
		for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
		L[m]=matrix[m]*accel-sum;		//加速係数をかける
		if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkkが非常に小さいとき、このようにする。ここで、どのみちLkk<0なら収束しづらいから、Lkkが負の極小値でも正の値にしてしまってよいのでは？
		D1[k]=1/L[m];

		if(L[m]*L[m]>maxvalue)
		{
			maxvalue=L[m]*L[m];
			Line=k; Column=k;
		}
		
		if(L[m]<0)						//L[m]<0のときD1[m]<0となり、収束が極端に遅くなる.これを防ぐために加速係数を大きくする
		{
			flag=UP;
		}
	}
	

	if(flag==UP)		//加速係数を１より大きくしていく場合
	{
		loopsw=ON;
		while(loopsw==ON)
		{
			accel+=0.05;
			L[0]=matrix[0]*accel;
			D1[0]=1/L[0];			//L[0]は係数行列の値と変わらないので、ゼロでもなければ負でもないと予想される
			maxvalue=D1[0]*D1[0];
			for(int k=1;k<pn;k++)	//k=0はもう済んだ。k>=1からloop
			{	
				for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0からk-1までの要素
				{
					int i=ind2[m];//列番号
		
					double sum=0;		//accelの変更によりLとD1の値が変わったので、sumも計算しなおし。
					for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i行目のなかから列の一致するものを探している。少し手間か？
					L[m]=matrix[m]-sum;//i=0のときはsum=0となり自動的にL[m]=matrix[m]となる

					if(L[m]*L[m]>maxvalue)
					{
						maxvalue=L[m]*L[m];
						Line=k; Column=i;
					}
				}
				int m=ptr2[k+1]-1;				//このmは必ず対角成分に該当する
				
				double akk=matrix[m];			//係数行列の対角成分.値を保存.
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
				L[m]=matrix[m]*accel-sum;		//加速係数をかける
				if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkkが非常に小さいとき、このようにする。ここで、どのみちLkk<0なら収束しづらいから、Lkkが負の極小値でも正の値にしてしまってよいのでは？
				D1[k]=1/L[m];
				if(L[m]*L[m]>maxvalue)
				{
					maxvalue=L[m]*L[m];
					Line=k; Column=k;
				}
			}
			loopsw=OFF;
			for(int k=0;k<pn;k++) if(D1[k]<0) loopsw=ON;	//1つでも負のDがあればもう一度実行する
		}

	}
	else if(flag==DOWN)		//加速係数を１より小さくしていく場合
	{
		loopsw=ON;
		int flag2=OFF;
		while(loopsw==ON)
		{
			if(flag2==OFF) accel-=0.05;
			else if(flag2==ON)
			{
				accel+=0.05;
				loopsw=OFF;
			}
			
			L[0]=matrix[0]*accel;
			D1[0]=1/L[0];			//L[0]は係数行列の値と変わらないので、ゼロでもなければ負でもないと予想される
			maxvalue=D1[0]*D1[0];
			for(int k=1;k<pn;k++)	//k=0はもう済んだ。k>=1からloop
			{	
				for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0からk-1までの要素
				{
					int i=ind2[m];//列番号
		
					double sum=0;		//accelの変更によりLとD1の値が変わったので、sumも計算しなおし。
					for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i行目のなかから列の一致するものを探している。少し手間か？
					L[m]=matrix[m]-sum;//i=0のときはsum=0となり自動的にL[m]=matrix[m]となる

					if(L[m]*L[m]>maxvalue)
					{
						maxvalue=L[m]*L[m];
						Line=k; Column=i;
					}
				}
				int m=ptr2[k+1]-1;				//このmは必ず対角成分に該当する
				
				double akk=matrix[m];			//係数行列の対角成分.値を保存.
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
				L[m]=matrix[m]*accel-sum;		//加速係数をかける
				if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkkが非常に小さいとき、このようにする。ここで、どのみちLkk<0なら収束しづらいから、Lkkが負の極小値でも正の値にしてしまってよいのでは？
				D1[k]=1/L[m];
				if(L[m]*L[m]>maxvalue)
				{
					maxvalue=L[m]*L[m];
					Line=k; Column=k;
				}
			}
			//loopsw=OFF;
			for(int k=0;k<pn;k++) if(D1[k]<0) flag2=ON;	//1つでも負のDがあれば、小さくしすぎたと判断して、1つ前の加速係数を用いて再度計算する
			if(Line!=Column) flag2=ON;
		}
	}
	//if(Line!=Column) cout<<"value="<<sqrt(maxvalue)<<" Line="<<Line<<" Column="<<Column<<" D["<<Line<<"]="<<D1[Line]<<endl;

	/*/対角優位性をチェック
	for(int k=0;k<pn;k++)
	{	
		double val=0;
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++) val+=fabs(L[m]);
		val*=2;
		if(val>L[ptr2[k+1]-1]) cout<<k<<" "<<val<<" "<<D1[k]<<endl;
	}//*/

	//cout<<"加速係数="<<accel<<" time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
	cout<<"γ="<<accel<<" ";
	
}

/*/不完全ｺﾚｽｷｰ分解(旧）
void Incomplete_Cholesky_Decomposition(mpsconfig *CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2)
{
	//matrix[i]:係数行列のうち、対角を含む下三角が格納されている
	//ptr2[k]:matrix[i]に対するｱｸｾｽ配列
	//L[i],D[i]:不完全ｺﾚｽｷｰ分解格納配列
	//pn:未知数
	//num2:matrix[i]の要素数

	cout<<"不完全ｺﾚｽｷｰ分解開始---";
	unsigned int timeC=GetTickCount();
	double accel=CON->get_CGaccl();		//加速係数
	double accel2=0;					//自動計算された加速係数
	double maxvalue=0;					//matrix中の最大値
	int Line=1;							//最大値を持つ行番号
	int Column=1;						//最大値を持つ列番号
	
	//1回目:まずは普通に不完全ｺﾚｽｷｰ分解を行う。そこでD1[i]<0なら加速係数を修正し、2回目を行う
	int onemoreflag=OFF;
	L[0]=matrix[0]*accel;
	D1[0]=1/L[0];			//L[0]は係数行列の値と変わらないので、ゼロでもなければ負でもないと予想される
	if(D1[0]<0) cout<<"error in 不完全ｺﾚｽｷｰ分解 D1[0]<0"<<endl;
	maxvalue=D1[0]*D1[0];
	for(int k=1;k<pn;k++)	//k=0はもう済んだ。k>=1からloop
	{	
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0からk-1までの要素
		{
			int i=ind2[m];//列番号
		
			double sum=0;
			for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i行目のなかから列の一致するものを探している。少し手間か？
			L[m]=matrix[m]-sum;//i=0のときはsum=0となり自動的にL[m]=matrix[m]となる
			if(L[m]*L[m]>maxvalue)
			{
				maxvalue=L[m]*L[m];
				Line=k; Column=i;
			}
		}
		int m=ptr2[k+1]-1;				//このmは必ず対角成分に該当する
				
		double akk=matrix[m];			//係数行列の対角成分.値を保存.
		double sum=0;
		for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
		L[m]=matrix[m]*accel-sum;		//加速係数をかける
		if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkkが非常に小さいとき、このようにする。ここで、どのみちLkk<0なら収束しづらいから、Lkkが負の極小値でも正の値にしてしまってよいのでは？
		D1[k]=1/L[m];
		
		if(L[m]<0)						//L[m]<0のときD1[m]<0となり、収束が極端に遅くなる.これを防ぐために加速係数を大きくする
		{
			accel2=sum/akk*1.1;
			if(accel2>accel) accel=accel2;
			onemoreflag=ON;				//このスイッチがONなら、不完全ｺﾚｽｷｰ分解をもう一度行う。
		}
		if(L[m]*L[m]>maxvalue)
		{
			maxvalue=L[m]*L[m];
			Line=k; Column=k;
		}
	}
	if(Line!=Column) cout<<"value="<<maxvalue<<" Line="<<Line<<" Column="<<Column<<endl;
	if(onemoreflag==ON)
	{
		L[0]=matrix[0]*accel;
		D1[0]=1/L[0];			//L[0]は係数行列の値と変わらないので、ゼロでもなければ負でもないと予想される
		for(int k=1;k<pn;k++)	//k=0はもう済んだ。k>=1からloop
		{	
			for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0からk-1までの要素
			{
				int i=ind2[m];//列番号
		
				double sum=0;		//accelの変更によりLとD1の値が変わったので、sumも計算しなおし。
				for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i行目のなかから列の一致するものを探している。少し手間か？
				L[m]=matrix[m]-sum;//i=0のときはsum=0となり自動的にL[m]=matrix[m]となる
			}
			int m=ptr2[k+1]-1;				//このmは必ず対角成分に該当する
				
			double akk=matrix[m];			//係数行列の対角成分.値を保存.
			double sum=0;
			for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
			L[m]=matrix[m]*accel-sum;		//加速係数をかける
			if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkkが非常に小さいとき、このようにする。ここで、どのみちLkk<0なら収束しづらいから、Lkkが負の極小値でも正の値にしてしまってよいのでは？
			D1[k]=1/L[m];
		}
	}	
	cout<<"加速係数="<<accel<<" time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
	///不完全ｺﾚｽｷｰ分解完了//////////
}
///*/

//行列の対称性検査関数
void check_matrix_symmetry(int pn,int *NUM,int **ROW,double **G)
{
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(abs((G[i][j]-G[J][k])/G[i][j])>1e-8)
					{
						//cout<<"対称性ｴﾗｰ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
						cout<<"対称性ｴﾗｰ ("<<i<<","<<J<<")="<<G[i][j]<<"  ("<<J<<","<<ROW[J][k]<<")="<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0) cout<<"DDD i="<<i<<endl;
		}
	}
}

void check_matrix_symmetry_complex(int pn,int *NUM,int **ROW,complex<double> **G)
{
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						cout<<"対称性ｴﾗｰ "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;
					}
				}
			}
			if(flag==0) cout<<"DDD i="<<i<<endl;
		}
	}
	cout<<"対称性チェック完了"<<endl;
}

void diagonal_distribution(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int pn,int *NUM,int **ROW,double **G)
{
	ofstream di("dia.dat");
	ofstream xd("exdia.dat");

	double *dia=new double [pn];	//対角成分
	double *exdia=new double [pn];	//対角成分以外の和

	double sum_di=0;
	double sum_xd=0;

	for(int i=1;i<=pn;i++)
	{
		double sum=0;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			if(i==J)
			{
				dia[i]=G[i][j];
				di<<i<<" "<<dia[i]<<endl;
				sum_di+=dia[i];
			}
			else
			{
				sum+=abs(G[i][j]);
			}
		}
		exdia[i]=sum;
		xd<<i<<" "<<sum<<endl;
		sum_xd+=sum;
	}
	cout<<"係数行列のデータ作成"<<endl;
	di.close();
	xd.close();

	cout<<"r_sum="<<sum_xd<<endl;
	sum_di/=pn; sum_xd/=pn; //平均値

	//分散の計算
	double s2di=0;
	double s2xd=0;
	for(int i=1;i<=pn;i++)
	{
		s2di+=(sum_di-dia[i])*(sum_di-dia[i]);
		s2xd+=(sum_xd-exdia[i])*(sum_xd-exdia[i]);
	}

	s2di/=pn;
	s2xd/=pn;
	
	cout<<"s2D="<<s2di<<" s2r="<<s2xd<<" aveD="<<sum_di<<" aver="<<sum_xd<<endl;

	//値が大きなところだけを抜粋する //平均値の10倍以上だけ抽出
	ofstream di2("dia_p.dat");
	ofstream xd2("exdia_p.dat");

	for(int i=1;i<=pn;i++)
	{
		if(dia[i]>sum_di*5) di2<<i<<" "<<dia[i]<<endl;
		if(exdia[i]>sum_xd*5) xd2<<i<<" "<<exdia[i]<<endl;
	}

	di2.close();
	xd2.close();

	delete [] dia;
	delete [] exdia;

}

//非対称行列解法
void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :非ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解
	
	int count=0;
	double pk=1;
	double E=1;//誤差
	double alp,beta,rr,w,ita;

	double *r=new double [pn];	//残差
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
	beta=0;
	w=0;
	ita=0;

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0;
		r[n]=B[n];
	}

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
	}
	double rr0=0;
	for(int n=0;n<pn;n++) rr0+=r[n]*r[n];
	cout<<"rr0="<<rr0<<endl;
	 cout<<"BiCGstab2法スタート  -----未知数="<<pn<<"  ---";
	 unsigned int time=GetTickCount();
	double ep=CON->get_FEMCGep();//収束判定
	while(E>ep)// EP=CON->get_CGep();//収束判定(convergence test)
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
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			//for(int m=0;m<pn;m++) AP[n]+=A[n][m]*P[m];
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

		for(int n=0;n<pn;n++)
		{
			Ae[n]=0;
			//for(int m=0;m<pn;m++) Ae[n]+=A[n][m]*e[m];
			for(int m=ptr[n];m<ptr[n+1];m++) Ae[n]+=val[m]*e[ind[m]];
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
			X[n]+=alp*P[n]+Z[n];
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
		E=rr/rr0;
		//E=sqrt(rr);
		//cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
	}
	
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;

	delete [] r;
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
}

//動磁場計算関数
void calc_transitional_EM_field(mpsconfig *CON,int node,int nelm,int nedge,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,int *jnb,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int **nei,int particle_node,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,int node_sta,vector<edge3D> &static_EDGE)
{
	cout<<"動磁場解析開始"<<endl;
	double *V=new double[node+1];
	for(int i=1;i<=node;i++) V[i]=0;				//初期化　Vは渦電流の電位格納に使う
	double *current[3];								//各要素の電流密度[A/m2]格納
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];
	double *Be[3];					//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//各要素の透磁率
	double *sigma=new double [nelm+1];	//各要素の導電率

	//比透磁率決定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=CON->get_RP();
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==CRUCIBLE) RP[i]=1;
		else cout<<"error:材質に対し、透磁率が不定"<<endl;
	}

	//導電率決定
	for(int i=1;i<=nelm;i++)
	{
		if(CON->get_temperature_depend()==OFF)
		{
			if(ELEM[i].material==FLUID) sigma[i]=CON->get_ele_conduc();
			else if(ELEM[i].material==AIR) sigma[i]=0;//使うことはないはずだが、一応初期化
			else if(ELEM[i].material==COIL) sigma[i]=CON->get_ele_conduc2();
			else if(ELEM[i].material==CRUCIBLE) sigma[i]=CON->get_ele_conduc2();
			else cout<<"error:材質に対し、導電率が不定"<<endl;
		}
	}

	if(CON->get_temperature_depend()==ON)
	{
		cout<<"温度依存の導電率計算";
		double *sig_f=new double[fluid_number];
		double *sig_n=new double[node+1];
		calc_physical_property(CON,PART,fluid_number,sig_f,fluid_number,4);//抵抗率の温度依存
		for(int i=0;i<fluid_number;i++) sig_f[i]=1/sig_f[i];
		for(int in=1;in<=node;in++) 
		{
			if(NODE[in].particleID>=0)
			{
				sig_n[in]=sig_f[NODE[in].particleID];
				if(NODE[in].material!=FLUID) cout<<"流体粒子が対応していない点で温度場計算"<<endl;
			}
		}

		//導電率を要素から節点に分配する
		for(int je=1;je<=nelm;je++)
		{
			double dis=0;
			double dis_sum=0;
			int N[4+1]; //要素の各節点番号格納
			double X[4+1];
			double Y[4+1];
			double Z[4+1];
			

			if(ELEM[je].material==FLUID)
			{	
				for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
	    
				double Xs=0;//重心座標
				double Ys=0;
				double Zs=0;
				for(int j=1;j<=4;j++)
				{
					X[j]=NODE[N[j]].r[A_X];
					Y[j]=NODE[N[j]].r[A_Y];
					Z[j]=NODE[N[j]].r[A_Z];
					Xs+=X[j]*0.25;
					Ys+=Y[j]*0.25;
					Zs+=Z[j]*0.25;
				}
		
				double delta=ELEM[je].volume/6;//本当の体積

				for(int j=1;j<=4;j++)
				{
					dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
					dis_sum+=dis;
				}
				for(int j=1;j<=4;j++)
				{
					dis=sqrt((Xs-X[j])*(Xs-X[j])+(Ys-Y[j])*(Ys-Y[j])+(Zs-Z[j])*(Zs-Z[j]));
					sigma[je]+=sig_n[N[j]]*dis/dis_sum;
				}
			}
			if(ELEM[je].material==AIR) sigma[je]=0;//使うことはないはずだが、一応初期化
			if(ELEM[je].material==COIL) sigma[je]=CON->get_ele_conduc2();
			if(ELEM[je].material==CRUCIBLE) sigma[je]=CON->get_ele_conduc2();

			//チェック
			if(ELEM[je].material==FLUID)
			{
				if(sigma[je]>1e8) cout<<"導電率大きすぎ？"<<endl;
				if(sigma[je]<1e6) cout<<"導電率小さすぎ？"<<endl;
			}

		}

		
		delete [] sig_f;
		delete [] sig_n;
		cout<<" ok"<<endl;
	}

	

	cout<<"電流修正"<<endl;
	double f=CON->get_Hz();	//交流の周波数[Hz]
	double omega=2*PI*f;//角周波数[rad/sec]
	if(CON->get_J_input_way()==1)					//電流密度を他のソフトから読み込む
	{
		inport_J0_density(CON, node, nelm,NODE,ELEM,current);		
		
		/*///
		////変換
		double I0=CON->get_I0();			//電流の振幅
		for(int i=1;i<=nelm;i++)
		{
			for(int D=0;D<3;D++)
			{
				current[D][i]*=I0/900;//900はmagnet側で指定している強制電流。
			}
		}
		///*/

		/*///
		//交流変換
		//double I0=CON->get_I0();			//電流の振幅
		//double II=I0*cos(omega*TIME);	
		//double II=I0*sin(2*PI*f*(TIME+CON->get_dt()));	
		for(int i=1;i<=nelm;i++)
		{
			for(int D=0;D<3;D++)
			{
				//current[D][i]*=II/I0;
				current[D][i]*=cos(omega*TIME);
			}
		}
		////*/

		////
		////大きさ変換
		//double I0=CON->get_I0();
		double I0=CON->get_I0()*sqrt(2.0);	//電流の振幅 //実験では実効値を測定しているので、√2をかける
		for(int i=1;i<=nelm;i++)
		{
			for(int D=0;D<3;D++)
			{
				current[D][i]*=I0/900;//900はmagnet側で指定している強制電流の大きさ。
			}
		}
		///*/

		check_J0(CON, node, nelm,NODE,ELEM,current);

		//強制電流の発散をチェック ∫(gradNi・J0)dv=0
		double div=0;
		int N[4+1]; //要素の各節点番号格納
		double X[4+1];
		double Y[4+1];
		double Z[4+1];
		//double b[4+1];
		double c[4+1];
		double d[4+1];
		double e[4+1];
		double j0x, j0y, j0z;
		for(int je=1;je<=nelm;je++)
		{   
			if(ELEM[je].material==COIL)
			{
				for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		    
				for(int j=1;j<=4;j++)
				{
					X[j]=NODE[N[j]].r[A_X];
					Y[j]=NODE[N[j]].r[A_Y];
					Z[j]=NODE[N[j]].r[A_Z];
				}
				
				///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
				ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
			
				double delta6=ELEM[je].volume;//体積の6倍
			
				delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
			
				double delta=ELEM[je].volume/6;//本当の体積
			
				for(int i=1;i<=4;i++)
				{
					int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
					int m=j%4+1;
					int n=m%4+1;
			    
					c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
					d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
					e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
					if(i%2!=0)//iが奇数なら
					{
						c[i]*=-1;
						d[i]*=-1;
						e[i]*=-1;
					}
				}
				/////////
				

				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
				}
				else
				{
					j0x=0;
					j0y=0;
					j0z=0;
				}
				
				for(int n=1;n<=4;n++)
				{
					div+=c[n]*j0x+d[n]*j0y+e[n]*j0z;
				}
			}
		}
		cout<<"強制電流発散="<<div<<endl;
	}
	else
	{
		//calc_current(CON,NODE,ELEM,EDGE,node,nelm,nedge,jnb,nei,branch_num,current);//修正する必要あり(boundary_condition)
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	cout<<"電流密度の設定完了"<<endl;

	if(CON->get_FEM_elm_type()==0)//節点要素
	{
		double *A[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ
		for(int D=0;D<3;D++) A[D]=new double [node+1];
		double *old_A[3];		
		for(int D=0;D<3;D++) old_A[D]=new double [node+1]; //1step前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ

		for(int i=1;i<=node;i++) for(int D=0;D<3;D++)
		{
			A[D][i]=0; old_A[D][i]=0;
		}
				
		/////old_Aに値を格納
		if(t==1 && CON->get_restart()==OFF)				//最初のｽﾃｯﾌﾟはゼロで初期化
		{
			for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old_A[D][i]=0;	//初期化
		}
		else//それ以外はﾌｧｲﾙから読み込み
		{	
			ifstream old("old_A.dat");
			if(!old) cout<<"cannot open old_A.dat"<<endl;
			old.unsetf(ifstream::dec);
			old.setf(ifstream::skipws);

			for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old>>old_A[D][i];
					
			old.close();

			//読み込んだ値ではリメッシュ領域内部の導体の節点番号がずれているので、粒子に記憶させた値を渡す
			for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) old_A[D][i]=PART[NODE[i].particleID].old_A[D];
		}

		//ベクトルポテンシャル解決
		if(CON->get_m_A()==0) Avector3D_node_eddy2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t);
		//if(CON->get_m_A()==1) Avector3D_node_eddy2_jw(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//節点ごとに実数→虚数
		//if(CON->get_m_A()==2) Avector3D_node_eddy2_jw2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//虚数全部→実数全部
		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==OFF) Avector3D_node_eddy2_jw3(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//複素数、
		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==ON)
		{
			int dir=1;
			if(dir==0)
			{
				if(t==1) Avector3D_node_eddy2_jw3(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//複素数
				//if(t>1) Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//複素数,静的要素ディリクレ化
				if(t>1) Avector3D_node_eddy2_jw5_ver2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw,node_sta);//複素数,静的要素ディリクレ化
			}
			else
			{
				 //Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//以前の解析で計算した静的要素の値を最初からディリクレ値として利用
				 Avector3D_node_eddy2_jw5_ver2(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw,node_sta);//複素数,静的要素ディリクレ化
			}

		}
		//if(CON->get_m_A()==4) Avector3D_node_eddy2_jw4(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//行方向、列方向で入れ替え
	
		//磁束密度解決 
		Bflux3D_node(CON,NODE,ELEM,node,nelm,A,Be,t,ON);

		int N[4+1]; //要素の各節点番号格納
		int *count_e=new int [node+1];	//各説点まわりにいくつ計算対象の要素があるか。幾何的な関係はすでに求められているが、流体境界部の外と内のどちらで平均するかで変わってくるのでここでもとめる
		for(int i=1;i<=node;i++) count_e[i]=0;
		double B_sum=0;

		//磁束密度を要素から節点に分配する
		for(int je=1;je<=nelm;je++)
		{
			int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
			if(CON->get_Je_crucible()==0)
			{
				if(ELEM[je].material==FLUID) flagje=ON;
			}
			else if(CON->get_Je_crucible()==1)
			{
				if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
			}
			else if(CON->get_Je_crucible()==-1)
			{
				flagje=OFF;
			}

			if(flagje==ON)
			{	
				for(int j=1;j<=4;j++)
				{
					N[j]=ELEM[je].node[j];
					if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
					{
						B_sum=sqrt((Be[A_X][je]*Be[A_X][je]+Be[A_Y][je]*Be[A_Y][je]+Be[A_Z][je]*Be[A_Z][je])/2.0);//実効値で出力
						NODE[N[j]].B+=B_sum;
						count_e[N[j]]=count_e[N[j]]+1;
					}
				}
			}
		}
		for(int i=1;i<=node;i++)
		{
			if(count_e[i]!=0)
			{
				NODE[i].B/=count_e[i];
			}
		}
				
		//////

		//節点力法
		NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

		
		////
		if(CON->get_m_A()==1 && CON->get_jw_Faverage()==ON)//１周期あたりに働く力積を評価するため、電磁力の平均値を求める
		{
			cout<<"１周期あたりの電磁力の平均計算"<<endl;
			double TIME2;//電磁力の周期は磁場や電流の周期の半分。TIME2はTIMEから電磁力の周期の半分だけ進めた値
			//TIME2=TIME+0.000025/4;
			TIME2=TIME+1/(4*CON->get_Hz());
			cout<<"TIME2="<<TIME2<<endl;
			double *F2[DIMENTION];//半周期ずれた時間の電磁力を格納				
			for(int D=0;D<DIMENTION;D++) F2[D]=new double [fluid_number];
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F2[D][i]=0;//初期化
			
			calc_jw_field_node(CON,NODE,ELEM,dt,TIME2,PART,fluid_number,F2,t,node,nelm);

			double Fsum=0;
			double Fsum2=0;
			 for(int i=0;i<fluid_number;i++)
			{
				Fsum+=sqrt(F[A_X][i]*F[A_X][i]+F[A_Y][i]*F[A_Y][i]+F[A_Z][i]*F[A_Z][i]);
				Fsum2+=sqrt(F2[A_X][i]*F2[A_X][i]+F2[A_Y][i]*F2[A_Y][i]+F2[A_Z][i]*F2[A_Z][i]);
			 }

			if(t==1)
			{
				ofstream fout("Fsum1.dat");
				fout.close();
			}

			ofstream fout2("Fsum1.dat",ios :: app);
			fout2<<Fsum<<endl;
			fout2.close();

			if(t==1)
			{
				ofstream fouts("Fsum2.dat");
				fouts.close();
			}

			ofstream fouts2("Fsum2.dat",ios :: app);
			fouts2<<Fsum2<<endl;
			fouts2.close();

			
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=(F[D][i]+F2[D][i])/2;//電磁力の平均値
			
			

			for(int D=0;D<DIMENTION;D++) delete [] F2[D];
		}
		////*/
				

		//電磁力スムージング　fem_smnが-1なら表面、0なら使わない、1なら全部
		smoothingF3D(CON,PART,fluid_number,F,t);


		delete [] count_e;
		for(int D=0;D<3;D++) delete [] A[D];
		for(int D=0;D<3;D++)delete [] old_A[D];
	}

	if(CON->get_FEM_elm_type()==1)//辺要素
	{
		/////
		//辺要素作成
		int *branch_num=new int[node+1];	//各節点が隣接する節点数(各節点から延びる辺の数)
		int max=1000;						//節点に隣接する最大節点数
		int **nei2=new int* [node+1];		//各節点の隣接する節点番号格納(neiは節点-要素であるが、nei2は節点-節点)
		for(int i=1;i<=node;i++) nei2[i]=new int [max];
		
		////辺要素生成 (節点要素を使用する場合でも、電流密度を求めるときに辺要素がほしい)
		int KTJ=node;				//このあと動的節点(流体)を格納しないといけないから、KTJを増加
		KTJ+=CON->get_add_points();
		int KTE=12*KTJ;
		nedge=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,EDGE,branch_num,nei2,KTE,static_EDGE,t,node_sta);

		//if(ELEM[i].edge

		for(int i=1;i<=node;i++) delete [] nei2[i];
		delete [] nei2;
		delete [] branch_num;
		//////*/

		double *A=new double [nedge+1];	//辺におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ
		double *Am=new double [nedge+1];	//辺におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの瞬時値(波高)	
		double *phi=new double [nedge+1];	//辺におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの位相遅れ
		double *old_A=new double [nedge+1]; //1step前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ

		for(int i=1;i<=nedge;i++) 
		{
			A[i]=0;	//辺におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ
			Am[i]=0;//辺におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの瞬時値(波高)	
			phi=0;	//辺におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの位相遅れ
			old_A=0; //1step前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ
		}
				
		/*////old_Aに値を格納
		if(t==1 && CON->get_restart()==OFF)				//最初のｽﾃｯﾌﾟはゼロで初期化
		{
			for(int i=1;i<=node;i++) old_A[i]=0;	//初期化
		}
		else//それ以外はﾌｧｲﾙから読み込み
		{	
			//１ステップ前の辺って何？ //リメッシュしても渦電流計算部が変形しない場合は考えられるが、流体計算だと・・・
			ifstream old("old_A.dat");
			if(!old) cout<<"cannot open old_A.dat"<<endl;
			old.unsetf(ifstream::dec);
			old.setf(ifstream::skipws);

			for(int i=1;i<=nedge;i++) old>>old_A[i];
					
			old.close();

			//for(int i=1;i<=node;i++) for(int D=0;D<3;D++) if(NODE[i].particleID>=0) old_A[D][i]=PART[NODE[i].particleID].old_A[D];
		}
		*/

		//Avector3D(CON,NODE,ELEM,EDGE,node,nelm,nedge,A,jnb,branch_num,current,RP);
		//ベクトルポテンシャル解決
		if(CON->get_m_A()==0) Avector3D_edge_eddy(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//辺要素
		//else if(CON->get_m_A()==1) Avector3D_edge_eddy_jw(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//辺要素、ｊω

		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==OFF) 
		{
			if(CON->get_parabolic_node_element()==OFF) Avector3D_edge_eddy_jw(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//辺要素、ｊω
			else if(CON->get_parabolic_node_element()==ON) Avector3D_edge_eddy_jw_with_parabolic_node_element(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);
		}
		if(CON->get_m_A()==1 && CON->get_static_dirichlet()==ON)
		{
			int dir=1;
			if(dir==0)
			{
				if(t==1) Avector3D_edge_eddy_jw(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega);//辺要素、ｊω
				//if(t>1) Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//複素数,静的要素ディリクレ化
				if(t>1) Avector3D_edge_eddy_jw2(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,node_sta,static_EDGE,Am,phi);//辺要素、ｊω静的要素ディリクレ化
			}
			else
			{
				 //Avector3D_node_eddy2_jw5(CON, node, nelm,NODE,ELEM,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,NODE_jw,ELEM_jw);//以前の解析で計算した静的要素の値を最初からディリクレ値として利用
				 Avector3D_edge_eddy_jw2(CON, node, nelm,nedge,NODE,ELEM,EDGE,A,jnb,nei,old_A, dt,V,current,RP,PART,sigma,t,TIME,omega,node_sta,static_EDGE,Am,phi);//辺要素、ｊω静的要素ディリクレ化
			}

		}
		//磁束密度解決 
		Bflux3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,A,Be,t,ON);
		
		double *val=new double[nelm+1];
		for(int i=1;i<=nelm;i++)
		{
			double B=sqrt(Be[A_X][i]*Be[A_X][i]+Be[A_Y][i]*Be[A_Y][i]+Be[A_Z][i]*Be[A_Z][i]);
			val[i]=B;
		}

		/*
		if(t==1)
		{		
			int flux=0;//セグメント断面
			data_avs2flux(CON,node,nelm,NODE,ELEM,val,t,flux);//断面図,磁束
			flux=1;//スリット断面
			data_avs2flux(CON,node,nelm,NODE,ELEM,val,t,flux);//断面図,磁束
		}
		*/

		delete [] val;



		//節点力法
		NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

		/////
		if(CON->get_m_A()==1 && CON->get_jw_Faverage()==ON)//１周期あたりに働く力積を評価するため、電磁力の平均値を求める
		{
			cout<<"１周期あたりの電磁力の平均計算"<<endl;
			double TIME2;//電磁力の周期は磁場や電流の周期の半分。TIME2はTIMEから電磁力の周期の半分だけ進めた値
			TIME2=TIME+1/(4*CON->get_Hz());
			cout<<"TIME2="<<TIME2<<endl;
			double *F2[DIMENTION];//半周期ずれた時間の電磁力を格納				
			for(int D=0;D<DIMENTION;D++) F2[D]=new double [fluid_number];
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F2[D][i]=0;//初期化
			
			calc_jw_field_edge(CON,NODE,ELEM,EDGE,dt,TIME2,PART,fluid_number,F2,t,node,nelm,nedge,Am,phi);

			double Fsum=0;
			double Fsum2=0;
			 for(int i=0;i<fluid_number;i++)
			{
				Fsum+=sqrt(F[A_X][i]*F[A_X][i]+F[A_Y][i]*F[A_Y][i]+F[A_Z][i]*F[A_Z][i]);
				Fsum2+=sqrt(F2[A_X][i]*F2[A_X][i]+F2[A_Y][i]*F2[A_Y][i]+F2[A_Z][i]*F2[A_Z][i]);
			 }

			if(t==1)
			{
				ofstream fout("Fsum1.dat");
				fout.close();
			}

			ofstream fout2("Fsum1.dat",ios :: app);
			fout2<<Fsum<<endl;
			fout2.close();

			if(t==1)
			{
				ofstream fouts("Fsum2.dat");
				fouts.close();
			}

			ofstream fouts2("Fsum2.dat",ios :: app);
			fouts2<<Fsum2<<endl;
			fouts2.close();

			
			for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=(F[D][i]+F2[D][i])/2;//電磁力の平均値
			
			

			for(int D=0;D<DIMENTION;D++) delete [] F2[D];
		}
		////*/


		//電磁力スムージング　fem_smnが-1なら表面、0なら使わない、1なら全部
		smoothingF3D(CON,PART,fluid_number,F,t);

		delete [] A;
		delete [] Am;
		delete [] phi;
		delete [] old_A;

		
		
	}

	for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) PART[i].F[D]=F[D][i];//求めた電磁力をPARTに格納
	output_F_scalar_with_AVS_for_linear(CON,NODE,ELEM,t,PART,node);

	delete [] V;
	delete [] RP;
	delete [] sigma;
	

	for(int D=0;D<3;D++) delete [] current[D];
	for(int D=0;D<3;D++) delete [] Be[D];
}

//jw法用計算関数(節点要素)
void calc_jw_field_node(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm)
{
	int node11=(int) NODE.size()-1;
	//int nelm=(int) ELEM.size()-1;
	//cout<<"node="<<node<<" NODE.size()="<<node11<<endl;

	ifstream a("Am.dat");
	if(!a) cout<<"cannot open Am.dat"<<endl;
	ifstream p("phi.dat");
	if(!p) cout<<"cannot open phi.dat"<<endl;

	double Am[3];
	double phi[3];

	double *A[3];									//節点におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ
	for(int D=0;D<3;D++) A[D]=new double [node+1];
	double *Be[3];					//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	for(int i=0;i<=node;i++) for(int D=0;D<3;D++) A[D][i]=0.0;

	double omega=2*PI*CON->get_Hz();

	a.unsetf(ifstream::dec);
	a.setf(ifstream::skipws);
	
	p.unsetf(ifstream::dec);
	p.setf(ifstream::skipws);			

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++)
		{
			Am[D]=0;
			phi[D]=0;
		}
		for(int D=0;D<3;D++)
		{
			a>>Am[D];
			p>>phi[D];
			A[D][i]=Am[D]*cos(omega*TIME+phi[D]);
		}
	}
	a.close();
	p.close();


	//磁束密度解決 
	Bflux3D_node(CON,NODE,ELEM,node,nelm,A,Be,t,OFF);

	double *RP=new double [nelm+1];	//各要素の透磁率

	//透磁率決定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=CON->get_RP();
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==CRUCIBLE) RP[i]=1;
		else cout<<"error:材質に対し、透磁率が不定"<<endl;
	}

	int *jnb=new int[node+1];///各節点に隣接する要素数格納
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//各節点の周辺要素番号格納
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	//節点力法
	NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

	//電磁力スムージング　fem_smnが-1なら表面、0なら使わない、1なら全部
	//smoothingF3D(CON,PART,fluid_number,F,t);

	for(int D=0;D<3;D++) delete [] A[D];
	delete [] RP;
	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
	for(int D=0;D<3;D++) delete [] Be[D];
}

//jw法用計算関数(辺要素)
void calc_jw_field_edge(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,vector<edge3D> &EDGE,double dt,double TIME,vector<mpsparticle> &PART,int fluid_number,double **F,int t,int node, int nelm,int nedge,double *Am, double *phi)
{
	//int node11=(int) NODE.size()-1;
	//int nelm=(int) ELEM.size()-1;
	//cout<<"node="<<node<<" NODE.size()="<<node11<<endl;

	//ifstream a("Am_e.dat");
	//if(!a) cout<<"cannot open Am_e.dat"<<endl;
	//ifstream p("phi_e.dat");
	//if(!p) cout<<"cannot open phi_e.dat"<<endl;

	//double Am;
	//double phi;

	
	double *A=new double [nedge+1];	//辺におけるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ
	double *Be[3];					//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	for(int i=0;i<=node;i++) A[i]=0.0;

	double omega=2*PI*CON->get_Hz();

	//a.unsetf(ifstream::dec);
	//a.setf(ifstream::skipws);
	
	//p.unsetf(ifstream::dec);
	//p.setf(ifstream::skipws);			

	for(int i=1;i<=nedge;i++)
	{
		//Am=0;
		//phi=0;
		
		//a>>Am;
		//p>>phi;
		A[i]=Am[i]*cos(omega*TIME+phi[i]);
		
	}
	//a.close();
	//p.close();


	//磁束密度解決 
	Bflux3D_edge(CON,NODE,ELEM,EDGE,node,nelm,nedge,A,Be,t,OFF);

	double *RP=new double [nelm+1];	//各要素の透磁率

	//透磁率決定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID) RP[i]=CON->get_RP();
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==CRUCIBLE) RP[i]=1;
		else cout<<"error:材質に対し、透磁率が不定"<<endl;
	}

	int *jnb=new int[node+1];///各節点に隣接する要素数格納
    set_jnb3D(NODE,ELEM,node,nelm,jnb);
    int **nei=new int* [node+1];
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];//各節点の周辺要素番号格納
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	//節点力法
	NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number,t);

	//電磁力スムージング　fem_smnが-1なら表面、0なら使わない、1なら全部
	//smoothingF3D(CON,PART,fluid_number,F,t);

	delete [] A;
	delete [] RP;
	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
	for(int D=0;D<3;D++) delete [] Be[D];
}


///電流密度計算関数(辺要素用)
void calc_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,int *jnb,int **nei,int *branch_num,double **current)
{
	//II:電流[A]

	cout<<"電流密度計算開始"<<endl;
	
	double p=1.68e-8;//銅の電気抵抗率[Ωm]

	int side_num2=0;//コイルを構成する辺数
	for(int i=1;i<=nedge;i++)
	{
		int ia=EDGE[i].node[1];
		int ib=EDGE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL) side_num2++;
	}///side_num2がもとまった

	int *side_id=new int [side_num2+1];//コイルを構成する辺番号格納
	side_num2=0;
	for(int i=1;i<=nedge;i++)
	{
		int ia=EDGE[i].node[1];
		int ib=EDGE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL)
		{
			side_num2++;
			side_id[side_num2]=i;//辺番号iを格納
		}
	}///コイルを構成する辺番号を管理


	//////////////////////////ｺｲﾙ辺の電流に関する境界条件を設定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm!=0)
				{
					if(ELEM[kelm].material==COIL) 
					{
						///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
						int p=ELEM[i].node[j];//境界三角に属さない節点
						for(int k=1;k<=6;k++)
						{
							
							int iside=ELEM[i].edge[k];
							
							int ia=EDGE[iside].node[1];
							int ib=EDGE[iside].node[2];
							if(ib<ia)
							{
								int temp=ia;
								ia=ib;
								ib=temp;
							}///これで必ずia<ibとなった
							
							if(ia!=p && ib!=p)//辺isideはコイル境界上ということ
							{
								if(NODE[ia].boundary_condition==11 || NODE[ib].boundary_condition==11)
								{
									EDGE[iside].boundary_condition=11;//11をひとつでも含んだ辺は自然境界条件
									//cout<<(NODE[ia].r[A_X]+NODE[ib].r[A_X])/2<<" "<<(NODE[ia].r[A_Y]+NODE[ib].r[A_Y])/2<<" "<<(NODE[ia].r[A_Z]+NODE[ib].r[A_Z])/2<<endl;
								}
								else
								{
									EDGE[iside].boundary_condition=10;//T=0となる辺
								}
								//T=0の面および自然境界条件面に固定境界辺を設置する
								if(NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==22) EDGE[iside].boundary_condition=21;//辺は21→22の方向
								else if(NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==21) EDGE[iside].boundary_condition=22;//辺は22→21の方向
								
							}
						}
					}
				}
			}
		}
	}//コイル端面に自然境界、側面にT=0の固定境界を敷いた。
	////////////////*/

	///境界条件出力　うまくいかないときにみる
	//data_avs_J_boundary(node,nelm,NODE,ELEM);

	for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
		//自然境界条件を未知数あつかいする
        if(EDGE[i].boundary_condition==11) EDGE[i].boundary_condition=0;
	}
	

	double II=CON->get_J0();
	double *T=new double [nedge+1];//電流ベクトルﾎﾟﾃﾝｼｬﾙ
    int NN=0;//ディリクレ型境界辺数
    int *dn=new int [nedge+1]; //各辺がディリクレ型境界の何番目か。ディリクレ型境界上でないならside_num2+1を格納
    double *PHAT=new double [CON->get_max_DN()];//ディリクレ型境値
    
	for(int k=1;k<=nedge;k++)T[k]=0;
    ///ディリクレ型境界条件入力
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(EDGE[i].boundary_condition==10)
		{
			dn[i]=NN;//i番目の辺はNN番目のディリクレ境界辺
	        PHAT[NN]=0;
	        T[i]=0;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==21)
		{    
			int ia=EDGE[i].node[1];
			int ib=EDGE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//辺の長さ
	        dn[i]=NN;
	        PHAT[NN]=II;
	        T[i]=II;
	        NN++;
		}
		else if(EDGE[i].boundary_condition==22)
		{   
			int ia=EDGE[i].node[1];
			int ib=EDGE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//辺の長さ
	        dn[i]=NN;
	        PHAT[NN]=-II;
	        T[i]=-II;
	        NN++;
		}
		else dn[i]=side_num2+1;
    }
	cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<endl;
	/////////////*/

	
    //////int pn=side_num-NN;				///未知数
	int pn=side_num2-NN;				///未知数
    int *ppn=new int [pn];			//行列のn番目は辺ppn[n]
    int *npp=new int [nedge+1];	///各辺が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(EDGE[i].boundary_condition==0)//未知数 
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////行列の最大幅計算
    int mat_w=0;
	for(int k=1;k<=side_num2;k++)
	{
		int i=side_id[k];
		int ia=EDGE[i].node[1];
		int ib=EDGE[i].node[2];
		int width=branch_num[ia]+branch_num[ib];//行列の幅
		if(width>mat_w) mat_w=width;
	}
	//mat_w*=5;
	////////////
	
    ////配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
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
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	
    for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
		
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//要素の重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;Ys/=4;Zs/=4;
			////////////////////////////
	
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
	        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
		
			////要素ﾏﾄﾘｸｽ作成開始
			for(int i=1;i<=6;i++)
			{	
				int iside=ELEM[je].edge[i];//要素jeの辺番号
				if(EDGE[iside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
				{   
					int I=npp[iside]+1;///辺isideは行列のI番目
					int I1=EDGE[iside].node[1];//isideを構成する2点
					int I2=EDGE[iside].node[2];
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
				    for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].edge[j];
						
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
							
						if(EDGE[jside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							int flag=0;
							//if(J<=I){
							int J1=EDGE[jside].node[1];//jsideを構成する2点
							int J2=EDGE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
								    flag=1;
								}
							}
							if(flag==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
							    ROW[I][H]=J;
							}
							//}
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
						{
						    int n=dn[jside];
						    B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6*PHAT[n];
						}//////////*/
					}
				}
			}
		}
    }
    ///////////////////////*/
	

    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
   
    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
		    for(int m=1;m<j;m++)
		    {
		        if(ROW[i][j]<ROW[i][m])
				{
				    tempG=G[i][m];
				    tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////

	///対称性チェック
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			for(int k=1;k<=NUM[J];k++) if(ROW[J][k]==i) if(G[i][j]!=G[J][k])
			{
				cout<<"matrix isn't symmetric   "<<G[i][j]<<" "<<G[J][k]<<endl;
			}
		}
	}///////////*/

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
    
    cout<<"行列作成終了  "<<endl;
    
    
	//CG3D(val,ind,ptr,pn,ppn,B,T);//CG法実行
	//ICCG3D(val,ind,ptr,pn,ppn,B,T,number);//ICCG法実行
	ICCG3D2(CON,val,ind,ptr,pn,B,number,T);
    ///////////
	
	denryu_side(CON,NODE,ELEM,EDGE,node,nelm,nedge,T,current);
    
	delete [] side_id;
	delete [] T;
    ///////////////////////*/
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

////辺要素電流密度計算関数
void denryu_side(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *T,double **current)
{
	int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係

	ofstream fp("j.dat");
	for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double Xs=0;//要素の重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;
			Ys/=4;
			Zs/=4;
			////////////////////////////
	
			double delta6=ELEM[je].volume;//体積の6倍(正しい体積の値はすでにベクトルポテンシャルを求める際に計算している)
		
			delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++) current[D][je]=0;//初期化
	
			for(int i=1;i<=6;i++)
			{
				int s=ELEM[je].edge[i];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
	
				current[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*T[s];
				current[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*T[s];
				current[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*T[s];
				
			}
	
			for(int D=0;D<3;D++) current[D][je]*=delta6*delta6*2;
			
			//if(Zs>0 && Zs<0.0001) fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;
			if(Zs>0.0005 && Zs<0.001)
			fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;		
		}
		
	}	
	fp.close();
}

//辺要素作成関数
int make_edge_element(mpsconfig *CON,vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,vector<edge3D> &EDGE,int *branch_num,int **nei2,int KTE,vector<edge3D> &static_EDGE,int t,int node_sta)
{
	unsigned timeA=GetTickCount();		//計算開始時刻

	cout<<"辺要素生成開始"<<endl;
	cout<<"node="<<node<<" ele="<<nelm<<endl;
	////辺要素生成
	int *check=new int [node+1];		//各節点がﾁｪｯｸしおわったかどうか
	int *flag=new int [node+1];		//各節点がﾁｪｯｸしおわったかどうか
	int *temp_check=new int[node+1];	//一時的なﾁｪｯｸ配列
	int max=1000;						//節点に隣接する最大節点数
	int *ele_side_num=new int[nelm+1];	//各要素の第何辺までがもとまっているか
	int side_num=0;						//全辺数格納
		
	///初期化
	for(int i=1;i<=node;i++) 
	{
		check[i]=0;
		temp_check[i]=0;
		branch_num[i]=0;//初期化
		flag[i]=0;
	}
	for(int i=1;i<=nelm;i++) ele_side_num[i]=0;
	/////////////////////

	//辺番号と辺-節点情報生成
	for(int i=1;i<=node;i++)
	{
		for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];//節点iに隣接する要素番号
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p!=i && temp_check[p]==0)
				{
					branch_num[i]=branch_num[i]+1;//節点iの隣接する節点数をプラス
					temp_check[p]=1;//もう見た
					nei2[i][branch_num[i]]=p;
					if(check[p]==0)
					{	
						side_num++;
						EDGE[side_num].node[1]=i;//このｱﾙｺﾞﾘｽﾞﾑ下においては、常にi<pである.なぜならiより小さな番号の節点はすでにcheck=1だから
						EDGE[side_num].node[2]=p;

						///節点ベースの境界条件を辺ベースに拡張
						if(NODE[i].boundary_condition==NODE[p].boundary_condition) EDGE[side_num].boundary_condition=NODE[i].boundary_condition;//両端が同じ境界条件なら、その辺もその境界条件に従う。仮に両方未知数でも問題ない
						else EDGE[side_num].boundary_condition=0;//それ以外は未知数
					}
				}
			}
		}
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//初期化
	}

	/*///
	//辺番号と辺-節点情報生成//この時点で静的要素を構成する節点は決定しており、節点番号も若い方から順に振られている。また、その番号はステップが進行しても変化しないはず。静的節点のみで構成される辺から順に番号を振る
	for(int i=1;i<=node;i++)
	{		
		if(i<=node_sta) 
		{
			//flag[i]=ON;
			for(int k=1;k<=jnb[i];k++)
			{
				int jelm=nei[i][k];//節点iに隣接する要素番号
				for(int j=1;j<=4;j++)
				{
					int p=ELEM[jelm].node[j];
					if(p<=node_sta)
					{
						//flag[p]=ON;
				
						if(p!=i && temp_check[p]==0)
						{
							branch_num[i]=branch_num[i]+1;//節点iの隣接する節点数をプラス
							temp_check[p]=1;//もう見た
							nei2[i][branch_num[i]]=p;
							if(check[p]==0)
							{	
								//if(flag[i]==ON && flag[p]==ON)
								{
									side_num++;
									EDGE[side_num].node[1]=i;//このｱﾙｺﾞﾘｽﾞﾑ下においては、常にi<pである.なぜならiより小さな番号の節点はすでにcheck=1だから
									EDGE[side_num].node[2]=p;

									///節点ベースの境界条件を辺ベースに拡張
									if(NODE[i].boundary_condition==NODE[p].boundary_condition) EDGE[side_num].boundary_condition=NODE[i].boundary_condition;//両端が同じ境界条件なら、その辺もその境界条件に従う。仮に両方未知数でも問題ない
									else EDGE[side_num].boundary_condition=0;//それ以外は未知数
								}
							}
						}
					}
				}
			}
		}
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//初期化
	}

	for(int i=1;i<=node;i++) check[i]=0;

	for(int i=1;i<=node;i++)
	{		
		if(i<=node_sta) flag[i]=ON;
		for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];//節点iに隣接する要素番号
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p<=node_sta) flag[p]=ON;
				
				if(p!=i && temp_check[p]==0)
				{
					branch_num[i]=branch_num[i]+1;//節点iの隣接する節点数をプラス
					temp_check[p]=1;//もう見た
					nei2[i][branch_num[i]]=p;
					if(check[p]==0)
					{	
						if(flag[i]==OFF || flag[p]==OFF)
						{
							side_num++;
							EDGE[side_num].node[1]=i;//このｱﾙｺﾞﾘｽﾞﾑ下においては、常にi<pである.なぜならiより小さな番号の節点はすでにcheck=1だから
							EDGE[side_num].node[2]=p;

							///節点ベースの境界条件を辺ベースに拡張
							if(NODE[i].boundary_condition==NODE[p].boundary_condition) EDGE[side_num].boundary_condition=NODE[i].boundary_condition;//両端が同じ境界条件なら、その辺もその境界条件に従う。仮に両方未知数でも問題ない
							else EDGE[side_num].boundary_condition=0;//それ以外は未知数
						}
					}
				}
				flag[p]=OFF;
			}
		}
		flag[i]=OFF;
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//初期化
	}
	///*/
	

	for(int i=1;i<=side_num;i++) EDGE[i].static_num=0;

	
	
	///////////*/

	///要素-辺情報生成
	for(int i=1;i<=nelm;i++) for(int j=1;j<=6;j++) ELEM[i].edge[j]=0;
	for(int i=1;i<=side_num;i++)
	{
		int ia=EDGE[i].node[1];//辺iを構成する節点番号(ia<ib)
		int ib=EDGE[i].node[2];
		for(int k=1;k<=jnb[ia];k++)
		{
			int jelm=nei[ia][k];//節点iaの隣接する要素番号
			for(int j=1;j<=4;j++)
			{
				if(ELEM[jelm].node[j]==ib)
				{
					ele_side_num[jelm]=ele_side_num[jelm]+1;
					int a=ele_side_num[jelm];
					ELEM[jelm].edge[a]=i;
				}
			}
		}
	}////////////

	//要素のremesh情報から、辺要素が静的か動的かを判断する
	//cout<<"静的辺要素の設定開始"<<endl;
	
	for(int i=1;i<=side_num;i++) EDGE[i].stat=OFF;//初期化
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].remesh==OFF)//リメッシュしない要素の持つ辺はリメッシュされない = 静的
		{
			for(int j=1;j<=6;j++)
			{
				EDGE[ELEM[i].edge[j]].stat=ON;
			}
		}
		
	}////////////
	
	
	/*/1stepの場合、static_edgeに記録しておく
	if(t==1)
	{
		int count=0;
		static_EDGE.resize(side_num+1);
		for(int i=1;i<=side_num;i++)
		{
			EDGE[i].static_num=i;
			static_EDGE[i]=EDGE[i];
			if(static_EDGE[i].stat==ON) count++;
		}
		static_EDGE[0].static_num=count;//静的辺要素の数を格納
		//cout<<"side_num="<<side_num<<" static_EDGE.size()="<<static_EDGE.size()-1<<endl;
	}

	if(t>1)
	{
		int count=0;
		
		for(int i=1;i<=static_EDGE.size()-1;i++)
		{
			if(static_EDGE[i].stat==ON)
			{
				//if(i%50000==0) cout<<i<<endl;
				if((NODE[EDGE[i].node[1]].r!=NODE[static_EDGE[i].node[1]].r) || (NODE[EDGE[i].node[2]].r!=NODE[static_EDGE[i].node[2]].r) )
				{
					//cout<<i<<endl;
					for(int j=1;j<=side_num;j++)
					{
						if(NODE[EDGE[j].node[1]].r==NODE[static_EDGE[i].node[1]].r && NODE[EDGE[j].node[2]].r==NODE[static_EDGE[i].node[2]].r)//1ステップ目と同一の辺。静的要素部分にしか存在しないはず
						{
							EDGE[j].static_num=i;
						}
					}
				}
				else  EDGE[i].static_num=i;
				count++;
				if(count==static_EDGE[0].static_num) break;
			}
		}
	}

	*/

	delete [] check;
	delete [] flag;
	delete [] temp_check;
	delete [] ele_side_num;

	//cout<<"辺数="<<side_num<<" 最大辺数="<<KTE<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	cout<<"辺数="<<side_num<<" 最大辺数="<<KTE<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	return side_num;
}

///境界条件適用関数(３Ｄ辺要素用)
void set_boundary_condition3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A)
{
	int N=0;	//数え上げ変数

	///ディリクレ型境界条件入力
	if(CON->get_uniform_B_sw()==OFF)
	{
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{    
		        dn[i]=N;
		        PHAT[N]=0;
		        A[i]=0;
		        N++;
			}
			else if(SIDE[i].boundary_condition==2)
			{   
				/*int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				A[i]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        dn[i]=NN;
		        PHAT[NN]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        //A[i]=0;
		        NN++;*/
				dn[i]=N;
		        PHAT[N]=0;
		        A[i]=0;
		        N++;
			}
			else dn[i]=side_num+1;
		}
	}
	if(CON->get_uniform_B_sw()==ON)	//解析領域全体に一様磁場を与える場合
	{
		double B=CON->get_uniform_B();//一様磁場の大きさ[Ｔ］
		double R[3];					//辺の長さ格納(X,Y,Z方向)
		double r[3];					//辺の中点座標格納
		double err=1e-14;
		
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{   
				int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				for(int D=0;D<3;D++)
				{
					R[D]=NODE[ib].r[D]-NODE[ia].r[D];
					r[D]=(NODE[ib].r[D]+NODE[ia].r[D])*0.5;//中点
				}

			//	if(r[A_Z]<CON->get_ZU()-err && r[A_Z]>CON->get_ZD()+err)
				{

					double L=sqrt(R[A_X]*R[A_X]+R[A_Y]*R[A_Y]+R[A_Z]*R[A_Z]);//辺の長さ
				
					A[i]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);	//なぜこうなるのかはストークスの定理を利用すればわかる。
					dn[i]=N;
					PHAT[N]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);
					N++;
				}
				//else 
				//{
				////	SIDE[i].boundary_condition=0;//自然境界条件
				//	dn[i]=side_num+1;
				//}
			}
			else dn[i]=side_num+1;
		}
	}

	*NN=N;
}

//電流密度読み込み関数
void inport_J0_density(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	cout<<"current_density.txtより強制電流分布を読み込み--";
	int id;

	ifstream fp("current_density.dat");
	if(!fp) cout<<"cannot open current_density.dat"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);

	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==COIL)						//current_density.datには、コイル要素のみ出力されている仕様にしておくこと
		{												//その場合、コイル要素は静的要素でなければならない。動的なら他のソフトから読み込めない
			fp>>id;
			for(int D=0;D<3;D++)
			{
				fp>>current[D][i];		//
			}
		}
		else for(int D=0;D<3;D++) current[D][i]=0;
	}
	fp.close();
	cout<<"ok"<<endl;
}

//電流密度出力関数
void check_J0(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	int coil_num=0;							//コイル要素数
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) coil_num++;

	ofstream fout2("J0.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << coil_num <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./J0 filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./J0 filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./J0 filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./J0 filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./J0 filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./J0 filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("J0");

	fout<<"e-x e-y e-z x y z"<<endl;
	double times=1;//1e-12;
	for(int i=1;i<=nelm;i++)
    {
		if(ELEM[i].material==COIL)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			fout<<current[A_X][i]*times<<" "<<current[A_Y][i]*times<<" "<<current[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			
		}
	}
	fout.close();
}

void check_J0Je(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	int coil_num=0;							//コイル要素数
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL ||ELEM[i].material==FLUID) coil_num++;

	ofstream fout2("J0Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << coil_num <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./J0Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./J0Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./J0Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./J0Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./J0Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./J0Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("J0Je");

	fout<<"e-x e-y e-z x y z"<<endl;
	double times=1;//1e-12;
	for(int i=1;i<=nelm;i++)
    {
		if(ELEM[i].material==COIL ||ELEM[i].material==FLUID)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			fout<<current[A_X][i]*times<<" "<<current[A_Y][i]*times<<" "<<current[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			
		}
	}
	fout.close();
}


////節点要素磁束密度計算関数
void Bflux3D_node(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **B,int t, int flag)
{
	cout<<"磁束密度計算----------";

	int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];

	double times=CON->get_B_times();
	double le=CON->get_distancebp();
	int plot_type=CON->get_plot_B_type();//1:ﾍﾞｸﾄﾙ　2:スカラー
	//ofstream fp2("test.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) B[D][je]=0;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		double delta6=ELEM[je].volume;//体積の6倍

		delta6=1/delta6;

		double Xs=0;
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];Ys+=Y[j];Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		//係数作成
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int i=1;i<=4;i++)
		{
			B[A_X][je]+=d[i]*A[A_Z][N[i]]-e[i]*A[A_Y][N[i]];
			B[A_Y][je]+=e[i]*A[A_X][N[i]]-c[i]*A[A_Z][N[i]];
			B[A_Z][je]+=c[i]*A[A_Y][N[i]]-d[i]*A[A_X][N[i]];
		}

		for(int D=0;D<3;D++) B[D][je]*=delta6;

		//if(Ys>-5*le && Ys<5*le) fp2<<Xs<<" "<<Zs<<" "<<B[A_X][je]*100<<" "<<B[A_Z][je]*100<<endl;

	}
	cout<<"ok"<<endl;
//	fp2.close();

	//磁束密度出力
	if(flag==ON)
	{
		int flagB=OFF;
		if(CON->get_m_A()==0) flagB=ON;
		if(CON->get_m_A()==1)
		{
			flagB=ON;
		}
		cout<<"磁束密度出力開始----";
		if(flagB==ON)
		{
			//磁束密度出力
			ofstream fp("Bflux.dat");
			
			//double Xmin=CON->get_XL()+le; double Xmax=CON->get_XR()-le;//解析領域
			//double Zmin=CON->get_ZD()+le; double Zmax=CON->get_ZU()-le;

			double Xmin=CON->get_XL()/2+le; double Xmax=CON->get_XR()/2-le;//解析領域
			double Zmin=0.05+le; double Zmax=CON->get_ZU()-le;
			
			double Rmax=CON->get_RU()-le;
			
			double dx; 
			if(plot_type==1) dx=1*le;
			if(plot_type==2) dx=0.01;
			//double dx=0.01;
			int Nx=(int)((Xmax-Xmin)/dx);//各方向の分割数
			int Nr=(int)(2*Rmax/dx);
			int Nz;
			if(plot_type==1) Nz=(int)((Zmax-Zmin)/dx);
			if(plot_type==2) Nz=15;
			int serch=nelm;//locate関数で最初に探索する要素番号
			
			if(CON->get_region_shape()==1)		//円筒領域のときは変数を書き換えて処理
			{
				Nx=Nr;
				Xmin=-Rmax;
			}


			if(plot_type==1)		//ﾍﾞｸﾄﾙ表示
			{
				for(int n=0;n<Nx;n++)
				{
					for(int m=0;m<Nz;m++)
					{
						double xp=dx*n+Xmin;//出力する点の座標
						double yp=0;
						double zp=dx*m+Zmin;
						int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
						fp<<xp<<" "<<zp<<" "<<B[A_X][loc]*times<<" "<<B[A_Z][loc]*times<<endl;
						serch=loc;
						//if(loc==0) cout<<"EE"<<endl;
					}
				}
			}
			else if(plot_type==2)	//スカラー表示
			{
				for(int n=0;n<1;n++)
				{
					for(int m=0;m<Nz;m++)
					{
						double xp=0;//出力する点の座標
						//double xp=dx*n+Xmin;//出力する点の座標
						double yp=0;
						double zp=dx*m+0.101;
						int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
						double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc])/sqrt(2.0);//実効値出力
						fp<<xp<<" "<<zp<<" "<<BB<<endl;
						serch=loc;
					}
				}
			}
			fp.close();//*/
		

			/////microAVS用の磁束密度出力

			int count=0;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
				if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)// &&ELEM[i].remesh==OFF)
				{
					count++;
				}
			}

			ofstream fout2("Bflux.fld");
			fout2 << "# AVS field file" << endl;
			fout2 << "ndim=1" << endl;
			//fout2 << "dim1=" << fluid_number <<endl;
			fout2 << "dim1=" << count<<endl;
			fout2 << "nspace=3" << endl;
			fout2 << "veclen=3" << endl;
			fout2 << "data=float" << endl;
			fout2 << "field=irregular" << endl;
			fout2 << "label=e-x e-y e-z" << endl << endl;
			fout2 << "variable 1 file=./Bflux filetype=ascii skip=1 offset=0 stride=6" << endl;
			fout2 << "variable 2 file=./Bflux filetype=ascii skip=1 offset=1 stride=6" << endl;
			fout2 << "variable 3 file=./Bflux filetype=ascii skip=1 offset=2 stride=6" << endl;
			fout2 << "coord    1 file=./Bflux filetype=ascii skip=1 offset=3 stride=6" << endl;
			fout2 << "coord    2 file=./Bflux filetype=ascii skip=1 offset=4 stride=6" << endl;
			fout2 << "coord    3 file=./Bflux filetype=ascii skip=1 offset=5 stride=6" << endl;
			fout2.close();

			ofstream fout("Bflux");
			fout<<"e-x e-y e-z x y z"<<endl;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)// && ELEM[i].material==OFF)
				{
					fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
			fout.close();

			//スリット部
			int count2=0;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
				if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
				{
					count2++;
				}
			}

			ofstream fout3("Bfluxs.fld");
			fout3 << "# AVS field file" << endl;
			fout3 << "ndim=1" << endl;
			//fout2 << "dim1=" << fluid_number <<endl;
			fout3 << "dim1=" << count2<<endl;
			fout3 << "nspace=3" << endl;
			fout3 << "veclen=3" << endl;
			fout3 << "data=float" << endl;
			fout3 << "field=irregular" << endl;
			fout3 << "label=e-x e-y e-z" << endl << endl;
			fout3 << "variable 1 file=./Bfluxs filetype=ascii skip=1 offset=0 stride=6" << endl;
			fout3 << "variable 2 file=./Bfluxs filetype=ascii skip=1 offset=1 stride=6" << endl;
			fout3 << "variable 3 file=./Bfluxs filetype=ascii skip=1 offset=2 stride=6" << endl;
			fout3 << "coord    1 file=./Bfluxs filetype=ascii skip=1 offset=3 stride=6" << endl;
			fout3 << "coord    2 file=./Bfluxs filetype=ascii skip=1 offset=4 stride=6" << endl;
			fout3 << "coord    3 file=./Bfluxs filetype=ascii skip=1 offset=5 stride=6" << endl;
			fout3.close();

			ofstream fout4("Bfluxs");
			fout4<<"e-x e-y e-z x y z"<<endl;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
				if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
				{
					fout4<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
			fout4.close();

			int flag2=0;
			if(CON->get_EM_interval()>1) flag2=ON;
			//else if(t==1 || t%10==0) flag=ON;


			///////////////////////////////
			if(flag2==ON)
			{
				char filename[25];
				sprintf_s(filename,"Bflux%d.fld", t);
				ofstream fout2(filename);
				fout2 << "# AVS field file" << endl;
				fout2 << "ndim=1" << endl;
				fout2 << "dim1=" << count <<endl;
				//fout2 << "dim1=" << BOnum <<endl;
				fout2 << "nspace=3" << endl;
				fout2 << "veclen=3" << endl;
				fout2 << "data=float" << endl;
				fout2 << "field=irregular" << endl;
				fout2 << "label=e-x e-y e-z" << endl << endl;
				fout2 << "variable 1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
				fout2 << "variable 2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
				fout2 << "variable 3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
				fout2 << "coord    1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
				fout2 << "coord    2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
				fout2 << "coord    3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
				fout2.close();

				char filename2[25];
				sprintf_s(filename2,"Bflux%d", t);
				ofstream fout(filename2);
				fout<<"e-x e-y e-z x y z"<<endl;
				for(int i=1;i<=nelm;i++)
				{
					double r[3]={0,0,0};
					for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
					if(r[A_Y]>-0.1*le && r[A_Y]<0.1*le)
					//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
					{
						fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
					}
				}
				fout.close();
			}
		}
		cout<<"ok"<<endl;
	}

}

//辺要素用磁束密度計算
void Bflux3D_edge(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE,int node,int nelm,int nedge,double *A,double **B,int t,int flag)
{
	cout<<"磁束密度計算開始----";
	unsigned timeA=GetTickCount();
	int plot_type=CON->get_plot_B_type();	//ﾌｧｲﾙ出力形式　1=ﾍﾞｸﾄﾙ 2=スカラー
	double times=CON->get_B_times();		//ﾌｧｲﾙ出力時の倍率
	double u0=4*PI*1e-7;					//真空の透磁率
	double le=CON->get_distancebp();

	double *Xg=new double [nelm+1];			//要素の重心座標
	double *Yg=new double [nelm+1];
	double *Zg=new double [nelm+1];


	//#pragma omp parallel for
	for(int je=1;je<=nelm;je++)
    {   
		//cout<<je<<" "<<omp_get_thread_num()<<endl;//各iの計算を担当しているｽﾚｯﾄﾞ番号出力
		int N[4+1];								//要素の各節点番号格納
		double X[4+1];
		double Y[4+1];
		double Z[4+1];
		double c[4+1];
		double d[4+1];
		double e[4+1];
		int table[6+1][2+1];					//要素を構成する辺とその要素内節点番号関係

		//辺－節点ﾃｰﾌﾞﾙ作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].edge[i];
			int ia=EDGE[iside].node[1];
			int ib=EDGE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//要素の重心座標
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j]*0.25;
			Ys+=Y[j]*0.25;
			Zs+=Z[j]*0.25;
		}
		Xg[je]=Xs; Yg[je]=Ys; Zg[je]=Zs;	//重心代入
		////////////////////////////

		double delta6=ELEM[je].volume;//体積の6倍(正しい体積の値はすでにベクトルポテンシャルを求める際に計算している)
	
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int D=0;D<3;D++) B[D][je]=0;//初期化

		for(int i=1;i<=6;i++)
		{
			int s=ELEM[je].edge[i];
			int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
			int k2=table[i][2];

			B[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*A[s];
			B[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*A[s];
			B[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*A[s];
		}

		for(int D=0;D<3;D++) B[D][je]*=delta6*delta6*2;
	}

	if(flag==ON)
	{
		//磁束密度出力
		ofstream fp("Bflux.dat");
	
		//double Xmin=CON->get_XL()+le; double Xmax=CON->get_XR()-le;//解析領域
		//double Zmin=CON->get_ZD()+le; double Zmax=CON->get_ZU()-le;

		double Xmin=CON->get_XL()/2+le; double Xmax=CON->get_XR()/2-le;//解析領域
		double Zmin=0.05+le; double Zmax=CON->get_ZU()-le;
	
		double Rmax=CON->get_RU()-le;
	
		//double dx=1*le;
		//int Nx=(int)((Xmax-Xmin)/dx);//各方向の分割数
		//int Nr=(int)(2*Rmax/dx);
		//int Nz=(int)((Zmax-Zmin)/dx);
		//int serch=nelm;//locate関数で最初に探索する要素番号
		double dx; 
		if(plot_type==1) dx=1*le;
		if(plot_type==2) dx=0.01;
		//double dx=0.01;
		int Nx=(int)((Xmax-Xmin)/dx);//各方向の分割数
		int Nr=(int)(2*Rmax/dx);
		int Nz;
		if(plot_type==1) Nz=(int)((Zmax-Zmin)/dx);
		if(plot_type==2) Nz=15;
		int serch=nelm;//locate関数で最初に探索する要素番号
	
		if(CON->get_region_shape()==1)		//円筒領域のときは変数を書き換えて処理
		{
			Nx=Nr;
			Xmin=-Rmax;
		}


		if(plot_type==1)		//ﾍﾞｸﾄﾙ表示
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//出力する点の座標
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
					fp<<xp<<" "<<zp<<" "<<B[A_X][loc]*times<<" "<<B[A_Z][loc]*times<<endl;
					serch=loc;
					//if(loc==0) cout<<"EE"<<endl;
				}
			}
		}
		else if(plot_type==2)	//スカラー表示
		{
			for(int n=0;n<1;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=0;//出力する点の座標
					//double xp=dx*n+Xmin;//出力する点の座標
					double yp=0;
					double zp=dx*m+0.101;
					int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
					double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc])/sqrt(2.0);//実効値出力
					fp<<xp<<" "<<zp<<" "<<BB<<endl;
					serch=loc;
				}
			}
		}
		/*//////
		else if(plot_type==2)	//スカラー表示
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//出力する点の座標
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,serch,xp,yp,zp);
					double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc]);
					fp<<xp<<" "<<zp<<" "<<BB<<endl;
					serch=loc;
				}
			}
		}
		/////*/
		fp.close();//*/

		/////microAVS用の磁束密度出力

		int count=0;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
			//if(ELEM[i].material==CRUCIBLE)
			if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
			//if(r[A_Y]>-0.01 && r[A_Y]<0.01)
			{
				count++;
			}
		}

		ofstream fout2("Bflux.fld");
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count<<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Bflux filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Bflux filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Bflux filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Bflux filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Bflux filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Bflux filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		ofstream fout("Bflux");
		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
			//if(ELEM[i].material==CRUCIBLE)
			//if(r[A_Y]>-0.01 && r[A_Y]<0.01)
			{
				fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			}
		}
		fout.close();

		//スリット部
		int count2=0;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
			if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
			{
				count2++;
			}
		}

		ofstream fout3("Bfluxs.fld");
		fout3 << "# AVS field file" << endl;
		fout3 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout3 << "dim1=" << count2<<endl;
		fout3 << "nspace=3" << endl;
		fout3 << "veclen=3" << endl;
		fout3 << "data=float" << endl;
		fout3 << "field=irregular" << endl;
		fout3 << "label=e-x e-y e-z" << endl << endl;
		fout3 << "variable 1 file=./Bfluxs filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout3 << "variable 2 file=./Bfluxs filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout3 << "variable 3 file=./Bfluxs filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout3 << "coord    1 file=./Bfluxs filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout3 << "coord    2 file=./Bfluxs filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout3 << "coord    3 file=./Bfluxs filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout3.close();

		ofstream fout4("Bfluxs");
		fout4<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
			if(r[A_Y]>sin(PI/24)*r[A_X]-0.1*le && r[A_Y]<sin(PI/24)*r[A_X]+0.1*le)
			{
				fout4<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			}
		}
		fout4.close();

		int flag2=0;
		if(CON->get_EM_interval()>1) flag2=ON;
		else if(t==1 || t%10==0) flag2=ON;


		///////////////////////////////
		if(flag2==ON)
		{
			char filename[25];
			sprintf_s(filename,"Bflux%d.fld", t);
			ofstream fout2(filename);
			fout2 << "# AVS field file" << endl;
			fout2 << "ndim=1" << endl;
			fout2 << "dim1=" << count <<endl;
			//fout2 << "dim1=" << BOnum <<endl;
			fout2 << "nspace=3" << endl;
			fout2 << "veclen=3" << endl;
			fout2 << "data=float" << endl;
			fout2 << "field=irregular" << endl;
			fout2 << "label=e-x e-y e-z" << endl << endl;
			fout2 << "variable 1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
			fout2 << "variable 2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
			fout2 << "variable 3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
			fout2 << "coord    1 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
			fout2 << "coord    2 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
			fout2 << "coord    3 file=./Bflux"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
			fout2.close();

			char filename2[25];
			sprintf_s(filename2,"Bflux%d", t);
			ofstream fout(filename2);
			fout<<"e-x e-y e-z x y z"<<endl;
			for(int i=1;i<=nelm;i++)
			{
				double r[3]={0,0,0};
				for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
				//if(r[A_Y]<0 && ELEM[i].material==CRUCIBLE)
				if(r[A_Y]>-0.5*le && r[A_Y]<0.5*le)
				{
					fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
			fout.close();
		}
	}

	delete [] Xg;
	delete [] Yg;
	delete [] Zg;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		
}

///節点力法計算関数
void NODE_F3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpsparticle> &PART,double **F,int fluid_number,int t)
{
    cout<<"節点力法による電磁力計算--------";
    double ep0=8.854e-12;	//真空の誘電率。
    double u0=12.5e-7;		//真空の透磁率
    int N[4+1];				//要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double Fz=0;			//Z方向の合力[N]
	double Fsum=0;			//合力[N]
	unsigned timeA=GetTickCount();//計算開始時刻

	double *Fn[3];
	for(int D=0;D<3;D++) Fn[D]=new double [node+1];//NNをnodeに変更。
    
    //for(int i=1;i<=node;i++) 
	if(CON->get_EM_calc_type()==3) //磁場解析
    {
        for(int I=1;I<=node;I++)
        {
			if(NODE[I].material==FLUID)
			{
			 for(int D=0;D<3;D++) Fn[D][I]=0;//初期化
			
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//節点iが隣接する要素番号
					
					//if(ELEM[jelm].material==AIR){
					///マクスウェルの応力テンソル
					double Txx=(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Txy=2*Ee[A_X][jelm]*Ee[A_Y][jelm];
					double Txz=2*Ee[A_X][jelm]*Ee[A_Z][jelm];
					double Tyx=Txy;
					double Tyy=(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Tyz=2*Ee[A_Y][jelm]*Ee[A_Z][jelm];
					double Tzx=Txz;
					double Tzy=Tyz;
					double Tzz=(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					//////////
				
					/////係数c,d,e計算
					for(int j=1;j<=4;j++)
					{
						N[j]=ELEM[jelm].node[j];
	    				X[j]=NODE[N[j]].r[A_X];
	    				Y[j]=NODE[N[j]].r[A_Y];
	    				Z[j]=NODE[N[j]].r[A_Z];
					}
					int i=0;///節点iは要素jelmの第J番目の節点
					for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;
					int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
					int m=j%4+1;
					int n=m%4+1;
					//delta6は相殺されるのでいらない
					double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
					double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
					double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
					if( i & 1 )//iが奇数なら
					{
						c*=-1;
						d*=-1;
						e*=-1;
					}
					///////////////
		    
					double u=RP[jelm];
						
					Fn[A_X][I]+=(Txx*c+Txy*d+Txz*e)/(2*u);
					Fn[A_Y][I]+=(Tyx*c+Tyy*d+Tyz*e)/(2*u);
					Fn[A_Z][I]+=(Tzx*c+Tzy*d+Tzz*e)/(2*u);
					//}
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/(6.00000000000*u0);
				//if(Fn[A_Z][I]<0) Fn[A_Z][I]=0;
				Fz+=Fn[A_Z][I];
				Fsum+=sqrt(Fn[A_X][I]*Fn[A_X][I]+Fn[A_Y][I]*Fn[A_Y][I]+Fn[A_Z][I]*Fn[A_Z][I]);
				
			}
		}
	}
    //cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	/*//////
	if(t==1)
	{
		ofstream fout("Fz.dat");
		fout.close();
	}

	ofstream fout2("Fz.dat",ios :: app);
	fout2<<Fz<<endl;
	fout2.close();

	if(t==1)
	{
		ofstream fouts("Fsum.dat");
		fouts.close();
	}

	ofstream fouts2("Fsum.dat",ios :: app);
	fouts2<<Fsum<<endl;
	fouts2.close();
    /////*/

	/*
	int *bound=new int[node+1];
	for(int n=1;n<=node;n++)
	{
		bound[n]=OFF;
		if(NODE[n].material==FLUID)
		{
			for(int k=1;k<=jnb[n];k++)
			{
				int jelm=nei[n][k];//節点iが隣接する要素番号
				if(ELEM[jelm].material==AIR) bound[n]=ON;
			}
		}
	}
	*/
	
	///F更新
	for(int n=1;n<=node;n++)
	{
		if(NODE[n].material==FLUID)
		{
			//if(bound[n]==ON)
			{
				int i=NODE[n].particleID;//i番目の流体粒子は、n番目の節点に相当
				if(i>=0) for(int D=0;D<3;D++) F[D][i]=Fn[D][n]; //-1番目の配列にアクセスするおそれがあるので条件付け
				
			}
		}
	}
	///節点力法の力の単位は[N]なので、ここからさらに粒子≠節点の粒子に関しても力をもとめてやる必要はない

	//delete [] bound;

	/*///ファイル出力
	ofstream fp("Fn.dat");
	double le=CON->get_distancebp();
	double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=1;i<=node;i++)//流体節点のみ出力
    {
		if(NODE[i].material==FLUID)
		{
			if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
		}
	}
	//凡例出力
	fp<<0.015<<" "<<0.155" "<<1.0e-002*times<<" "<<0<<endl;

    fp.close();///////////////////

	////ファイル出力//(スリットを通る断面)
	ofstream fps("Fnslit.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=1;i<=node;i++)//流体節点のみ出力
    {
		if(NODE[i].material==FLUID)
		{
			if(NODE[i].r[A_Y]>-le*0.5+sin(PI/24)*NODE[i].r[A_X] && NODE[i].r[A_Y]<le*0.5+sin(PI/24)*NODE[i].r[A_X])
			{
				fps<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<(Fn[A_X][i]*cos(PI/24)+Fn[A_Y][i]*sin(PI/24))*times<<" "<<Fn[A_Z][i]*times<<endl;
			}
		}
	}
    fps.close();///////////////////

	////ファイル出力(Z断面)
	ofstream fpz("Fnz.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=1;i<=node;i++)//流体節点のみ出力
    {
		if(NODE[i].material==FLUID)
		{
			if(NODE[i].r[A_Z]>-le*0.5+0.14125 && NODE[i].r[A_Z]<le*0.5+0.14125) fpz<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Y][i]*times<<endl;	
		}
	}
    fpz.close();///////////////////

	//if(t=1 || t%10==0)
	if(CON->get_EM_interval()>1)
	{
		char filename[20];
		sprintf_s(filename,"Fn%d.dat", t);
		ofstream fp(filename);

		for(int i=1;i<=node;i++)//流体節点のみ出力
		{
			if(NODE[i].material==FLUID)
			{
				if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
			}
		}
		fp.close();///////////////////
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Fn%d.dat", t);
		ofstream fp(filename);

		for(int i=1;i<=node;i++)//流体節点のみ出力
		{
			if(NODE[i].material==FLUID)
			{
				if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
			}
		}
		fp.close();///////////////////
	}
	///*/

	for(int D=0;D<3;D++) delete [] Fn[D];
	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

}

void calc_eddy_current(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **old_A,double dt,double *V,double **Je,int t, double *sigma)
{
	//double sigma=CON->get_ele_conduc();//電気伝導率
	double H=CON->get_height();			//ﾌｧｲﾙ出力用高さパラメータ
	double times=CON->get_eddy_times();//4e-12;
	double le=CON->get_distancebp();

	///Je=-σ(dA/dt+gradφ) ---(★)

	cout<<"渦電流密度計算----------";

	int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	double dA[3];//ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの差分を格納（X,Y,Z)
	double gradV[3];//gradφ格納

	
	//ofstream fout("Je.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) Je[D][je]=0;

		if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		{
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6;
	
			double Xs=0;
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			//係数作成
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;	
			
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++)
			{
				dA[D]=0;
				gradV[D]=0;
			}
			for(int j=1;j<=4;j++)
			{
				for(int D=0;D<3;D++) dA[D]+=(b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*(A[D][N[j]]-old_A[D][N[j]]);
			}
			for(int D=0;D<3;D++) dA[D]*=delta6/dt;//式（★）の第1項計算完了

			for(int j=1;j<=4;j++)//gradφの計算
			{
				gradV[A_X]+=c[j]*V[N[j]];
				gradV[A_Y]+=d[j]*V[N[j]];
				gradV[A_Z]+=e[j]*V[N[j]];
			}	
			for(int D=0;D<3;D++) gradV[D]*=delta6;//式（★）の第2項計算完了

			//for(int D=0;D<3;D++) Je[D][je]=-sigma*(dA[D]+gradV[D]);
			for(int D=0;D<3;D++) Je[D][je]=-sigma[je]*(dA[D]+gradV[D]);

			//if(Zs>H && Zs<H+le) fout<<Xs<<" "<<Ys<<" "<<Je[A_X][je]*times<<" "<<Je[A_Y][je]*times<<endl;
			//sqrt(jx*jx+jy*jy)
		}
	}
	
	//出力
	int count=0;//Je.fld
	int count2=0;//Jez.fld
	double h=CON->get_height();
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID) count++;
		//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//溶融金属の中心からの断面
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
					count++;
			}
		}
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-le && r[A_Z]<0.16125+le)
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			count2++;
		}
		}
	}

	ofstream fout2("Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count/*nelm*/ <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("Je");
	fout<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//溶融金属の中心からの断面
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
				{
					fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}	
		}
		}
	}
	fout.close();

	ofstream fout3("Jez.fld");
	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout3 << "dim1=" << count2/*nelm*/ <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Jez filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Jez filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Jez filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Jez filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Jez filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Jez filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	ofstream fout4("Jez");
	fout4<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			fout4<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout4.close();

	///読み込み用ファイル作成
	ofstream g("Je.dat");
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==CRUCIBLE) g<<Je[A_X][i]<<" "<<Je[A_Y][i]<<" "<<Je[A_Z][i]<<endl;
	g.close();

	/////////////////////////////
	if(CON->get_EM_interval()>1 && t==1)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count/*nelm*/ <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count/*nelm*/ <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}



	cout<<"ok"<<endl;
	//fout.close();

}

void calc_node_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **AR,double **AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t,double TIME, double *sigma,double omega)
{
	//double sigma=CON->get_ele_conduc();//電気伝導率
	double H=CON->get_height();			//ﾌｧｲﾙ出力用高さパラメータ
	double times=CON->get_eddy_times();//4e-12;
	double le=CON->get_distancebp();

	///Je=-σ(dA/dt+gradφ) ---(★)

	cout<<"渦電流密度計算(jw法)----------";
	
	complex<double> *Jec[3];
	for(int D=0;D<3;D++) Jec[D]=new complex<double> [nelm+1];
	
    for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Jec[D][i]=complex<double>(0.0,0.0);//初期化

	int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
//	double dA[3];//ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの差分を格納（X,Y,Z)
	double temp_AR[3];
	double temp_AI[3];
	double temp_VR[3];
	double temp_VI[3];

	
	//ofstream fout("Je.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) Je[D][je]=0;

		if(ELEM[je].material==FLUID)// || ELEM[je].material==CRUCIBLE)
		{
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6;
	
			double Xs=0;
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double co=omega*sigma[je]*delta6;
			double co2=sigma[je]*delta6;
			//係数作成
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;	
			
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++)
			{
				 temp_AR[D]=0;
				 temp_AI[D]=0;
				 temp_VR[D]=0;
				 temp_VI[D]=0;
			}
			for(int j=1;j<=4;j++)
			{

				for(int D=0;D<3;D++)
				{
					temp_AR[D]+=co*((b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*(AR[D][N[j]]));
					temp_AI[D]+=co*((b[j]+c[j]*Xs+d[j]*Ys+e[j]*Zs)*(AI[D][N[j]]));
				}
				
				temp_VR[A_X]+=co2*c[j]*VR[N[j]];
				temp_VR[A_Y]+=co2*d[j]*VR[N[j]];
				temp_VR[A_Z]+=co2*e[j]*VR[N[j]];
				temp_VI[A_X]+=co2*c[j]*VI[N[j]];
				temp_VI[A_Y]+=co2*d[j]*VI[N[j]];
				temp_VI[A_Z]+=co2*e[j]*VI[N[j]];
				
				/*
				temp_VR[A_X]+=co*c[j]*VR[N[j]];
				temp_VR[A_Y]+=co*d[j]*VR[N[j]];
				temp_VR[A_Z]+=co*e[j]*VR[N[j]];
				temp_VI[A_X]+=co*c[j]*VI[N[j]];
				temp_VI[A_Y]+=co*d[j]*VI[N[j]];
				temp_VI[A_Z]+=co*e[j]*VI[N[j]];
				*/
			}
				
			//for(int D=0;D<3;D++) gradV[D]*=delta6;//式（★）の第2項計算完了
			for(int D=0;D<3;D++)
			{
				double t_R=-temp_AI[D]+temp_VR[D];
				double t_I=temp_AR[D]+temp_VI[D];
				//double t_R=-temp_AI[D]-temp_VI[D];
				//double t_I=temp_AR[D]+temp_VR[D];
				Jec[D][je]=complex<double>(-t_R,-t_I);
			}

		}
			
	}

	double Jem[3];
	double phi[3];

	for(int je=1;je<=nelm;je++)
	{
		Je_loss[je]=0.0;
		for(int D=0;D<3;D++)
		{
			Jem[D]=0.0;
			phi[D]=0.0;
		}
		
		//if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		if(ELEM[je].material==FLUID)
		{
			double Jem_sum=0.0;
			for(int D=0;D<3;D++)
			{
				double Re=Jec[D][je].real();
				double Im=Jec[D][je].imag();
				Jem[D]=sqrt(Re*Re+Im*Im);
				Jem_sum+=Jem[D]*Jem[D];
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(Im,Re);
				Je[D][je]=Jem[D]*cos(omega*TIME+phi[D]);
			}
			Jem_sum=sqrt(Jem_sum);
			//渦電流損の計算
			Je_loss[je]=Jem_sum*Jem_sum*ELEM[je].volume/(2*sigma[je]*6);
		}
		//a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		//p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	//a.close();
	//p.close();

	cout<<"ok"<<endl;

	//
	//渦電流を要素から節点に分配する(出力用)
	//int N[4+1]; //要素の各節点番号格納
	//double X[5];
	//double Y[5];
	//double Z[5];
	double Je_sum=0;

	
	int *count_e=new int [node+1];	//各説点まわりにいくつ計算対象の要素があるか。幾何的な関係はすでに求められているが、流体境界部の外と内のどちらで平均するかで変わってくるのでここでもとめる

	//渦電流密度を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++)
			{
				N[j]=ELEM[je].node[j];
				if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
				{
					Je_sum=sqrt((Je[A_X][je]*Je[A_X][je]+Je[A_Y][je]*Je[A_Y][je]+Je[A_Z][je]*Je[A_Z][je])/2.0);//実効値で出力
					NODE[N[j]].Je+=Je_sum;
					count_e[N[j]]=count_e[N[j]]+1;
				}
			}
		}
	}
	
	for(int i=1;i<=node;i++)
	{
		if(count_e[i]!=0)
		{
			NODE[i].Je/=count_e[i];
		}
	}
	/*/////////
	//出力
	cout<<"渦電流密度出力------";
	int count=0;//Je.fld
	int count2=0;//Jez.fld
	double h=CON->get_height();
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID) count++;
		//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//溶融金属の中心からの断面
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
					count++;
			}
		}
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-le && r[A_Z]<0.16125+le)
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			count2++;
		}
		}
	}

	ofstream fout2("Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count<<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("Je");
	fout<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//溶融金属の中心からの断面
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
				{
					fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}	
		}
		}
	}
	fout.close();

	ofstream fout3("Jez.fld");
	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout3 << "dim1=" << count2 <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Jez filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Jez filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Jez filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Jez filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Jez filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Jez filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	ofstream fout4("Jez");
	fout4<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		{
			fout4<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout4.close();

	///読み込み用ファイル作成
	ofstream g("Je.dat");
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==CRUCIBLE) g<<Je[A_X][i]<<" "<<Je[A_Y][i]<<" "<<Je[A_Z][i]<<endl;
	g.close();

	/////////////////////////////
	if(CON->get_EM_interval()>1 && t==1)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	cout<<"ok"<<endl;
	/////*/



	
	//fout.close();
	for(int D=0;D<DIMENTION;D++) delete [] Jec[D];
	delete [] count_e;

}

//辺要素周波数応答解析用渦電流計算関数　高橋　三次元有限要素法　p79の式(4.36)はjwがかっこの外に出ているが間違い。本来はAkの項だけにかかることに注意
void calc_edge_eddy_current_jw(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<edge3D> &EDGE, int node,int nelm,double *AR,double *AI,double dt,double *VR,double *VI,double **Je,double *Je_loss,int t,double TIME, double *sigma,double omega)
{
	//double sigma=CON->get_ele_conduc();//電気伝導率
	double H=CON->get_height();			//ﾌｧｲﾙ出力用高さパラメータ
	double times=CON->get_eddy_times();//4e-12;
	double le=CON->get_distancebp();

	///Je=-σ(dA/dt+gradφ) ---(★)

	cout<<"辺要素渦電流密度計算(jw法)----------";
	
	complex<double> *Jec[3];
	for(int D=0;D<3;D++) Jec[D]=new complex<double> [nelm+1];
	
    for(int D=0;D<3;D++) for(int i=1;i<=nelm;i++) Jec[D][i]=complex<double>(0.0,0.0);//初期化

	int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
//	double dA[3];//ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの差分を格納（X,Y,Z)
	double temp_AR[3];
	double temp_AI[3];
	double temp_VR[3];
	double temp_VI[3];

	
	//ofstream fout("Je.dat");
    for(int je=1;je<=nelm;je++)
    {   
		
		for(int D=0;D<3;D++)
		{
			Je[D][je]=0;  temp_AR[D]=0; temp_AI[D]=0; temp_VR[D]=0; temp_VI[D]=0;
		}

		
		if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].edge[i];
				int ia=EDGE[iside].node[1];
				int ib=EDGE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;//辺iの端にある節点は、要素jeのj番目の節点
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6;
	
			double Xs=0;
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
	
			double co=omega*sigma[je]*delta6*delta6;//jがかかっているが。各項の計算でつじつまを合わせる
			double co2=sigma[je]*delta6;
			//係数作成
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;	
			
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////

			for(int i=1;i<=6;i++)//div(∂A∂t)
			{	
				int iside=ELEM[je].edge[i];//要素jeの辺番号
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];

				temp_AR[A_X]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*(-AI[iside]);//渦電流のx成分実部、div(∂A)側
				temp_AI[A_X]+=co*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*AR[iside];//

				temp_AR[A_Y]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*(-AI[iside]);
				temp_AI[A_Y]+=co*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*AR[iside];

				temp_AR[A_Z]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*(-AI[iside]);
				temp_AI[A_Z]+=co*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*AR[iside];
				
			}

			for(int j=1;j<=4;j++)//gradφ
			{	
				temp_VR[A_X]+=co2*c[j]*VR[N[ j]];
				temp_VI[A_X]+=co2*c[j]*VI[N[ j]];
				
				temp_VR[A_Y]+=co2*d[j]*VR[N[ j]];
				temp_VI[A_Y]+=co2*d[j]*VI[N[ j]];

				temp_VR[A_Z]+=co2*e[j]*VR[N[ j]];
				temp_VI[A_Z]+=co2*e[j]*VI[N[ j]];
			}
				
			//for(int D=0;D<3;D++) gradV[D]*=delta6;//式（★）の第2項計算完了
			for(int D=0;D<3;D++)
			{
				double t_R=temp_AR[D]+temp_VR[D];
				double t_I=temp_AI[D]+temp_VI[D];
				//double t_R=-temp_AI[D]-temp_VI[D];
				//double t_I=temp_AR[D]+temp_VR[D];
				Jec[D][je]=complex<double>(-t_R,-t_I);
			}

		}
			
	}

	double Jem[3];
	double phi[3];

	for(int je=1;je<=nelm;je++)
	{
		Je_loss[je]=0.0;
		for(int D=0;D<3;D++)
		{
			Jem[D]=0.0;
			phi[D]=0.0;
		}
		
		if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE)
		//if(ELEM[je].material==FLUID)
		{
			double Jem_sum=0.0;
			for(int D=0;D<3;D++)
			{
				double Re=Jec[D][je].real();
				double Im=Jec[D][je].imag();
				Jem[D]=sqrt(Re*Re+Im*Im);
				Jem_sum+=Jem[D]*Jem[D];
				//phi=atan(AI[D][i]/AR[D][i]);
				phi[D]=atan2(Im,Re);
				Je[D][je]=Jem[D]*cos(omega*TIME+phi[D]);
			}
			Jem_sum=sqrt(Jem_sum);
			//渦電流損の計算
			Je_loss[je]=Jem_sum*Jem_sum*ELEM[je].volume/(2*sigma[je]*6);
		}
		//a<<Am[A_X]<<" "<<Am[A_Y]<<" "<<Am[A_Z]<<endl;
		//p<<phi[A_X]<<" "<<phi[A_Y]<<" "<<phi[A_Z]<<endl;
	}
	//a.close();
	//p.close();

	cout<<"ok"<<endl;

	//
	//渦電流を要素から節点に分配する(出力用)
	//int N[4+1]; //要素の各節点番号格納
	//double X[5];
	//double Y[5];
	//double Z[5];
	double Je_sum=0;

	
	int *count_e=new int [node+1];	//各説点まわりにいくつ計算対象の要素があるか。幾何的な関係はすでに求められているが、流体境界部の外と内のどちらで平均するかで変わってくるのでここでもとめる

	//渦電流密度を要素から節点に分配する
	for(int je=1;je<=nelm;je++)
	{
		int flagje=OFF;//注目している要素で渦電流損を計算するかしないか
		if(CON->get_Je_crucible()==0)
		{
			if(ELEM[je].material==FLUID) flagje=ON;
		}
		else if(CON->get_Je_crucible()==1)
		{
			if(ELEM[je].material==FLUID || ELEM[je].material==CRUCIBLE) flagje=ON;
		}
		else if(CON->get_Je_crucible()==-1)
		{
			flagje=OFF;
		}

		if(flagje==ON)
		{	
			for(int j=1;j<=4;j++)
			{
				N[j]=ELEM[je].node[j];
				if(NODE[N[j]].material==FLUID || NODE[N[j]].material==CRUCIBLE)
				{
					Je_sum=sqrt((Je[A_X][je]*Je[A_X][je]+Je[A_Y][je]*Je[A_Y][je]+Je[A_Z][je]*Je[A_Z][je])/2.0);//実効値で出力
					NODE[N[j]].Je+=Je_sum;
					count_e[N[j]]=count_e[N[j]]+1;
				}
			}
		}
	}
	
	for(int i=1;i<=node;i++)
	{
		if(count_e[i]!=0)
		{
			NODE[i].Je/=count_e[i];
		}
	}
	//////////
	//出力
	cout<<"渦電流密度出力------";
	int count=0;//Je.fld
	int count2=0;//Jez.fld
	double h=CON->get_height();
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID) count++;
		//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
		//if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//溶融金属の中心からの断面
		//if(r[A_Y]<=0)
		//if(r[A_Y]>=0 && r[A_X]>=0)
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//for(int j=1;j<=4;j++)
			{
				//int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
					count++;
			}
		}
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-le && r[A_Z]<0.16125+le)
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>-0.01 && r[A_Z]<0.01)
			//if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
			{
				count2++;
			}
		}
	}

	ofstream fout2("Je.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count<<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Je filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Je filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Je filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Je filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Je filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Je filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("Je");
	fout<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(r[A_Z]<=0.13125+le*0.5 &&  r[A_Z]>=0.13125-le*0.5)//溶融金属の中心からの断面
		//if(r[A_Y]>=0 && r[A_X]>=0)
		{
			if(ELEM[i].material==FLUID ||ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X] && r[A_Y]<0)
		{
			//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
			//for(int j=1;j<=4;j++)
			{
				//int jelm=ELEM[i].elm[j];
				//if(ELEM[jelm].material==AIR)
				{
					fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}	
		}
		}
	}
	fout.close();

	ofstream fout3("Jez.fld");
	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout3 << "dim1=" << count2 <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Jez filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Jez filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Jez filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Jez filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Jez filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Jez filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	ofstream fout4("Jez");
	fout4<<"e-x e-y e-z x y z"<<endl;
	
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
		//if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
		//if(ELEM[i].material==CRUCIBLE && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		//if(ELEM[i].material==CRUCIBLE && r[A_Z]>0.16125-0.5*le && r[A_Z]<0.16125+0.5*le && r[A_X]>0 && r[A_Y]<sin(PI/24)*r[A_X] && r[A_Y]>-sin(PI/24)*r[A_X])
		if(ELEM[i].material==CRUCIBLE && r[A_Z]>-0.01 && r[A_Z]<0.01)
		{
			fout4<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout4.close();

	///読み込み用ファイル作成
	ofstream g("Je.dat");
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==CRUCIBLE) g<<Je[A_X][i]<<" "<<Je[A_Y][i]<<" "<<Je[A_Z][i]<<endl;
	g.close();

	/////////////////////////////
	if(CON->get_EM_interval()>1 && t==1)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	else if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"Je%d.fld", t);
		ofstream fout2(filename);
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		//fout2 << "dim1=" << fluid_number <<endl;
		fout2 << "dim1=" << count <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./Je"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./Je"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./Je"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./Je"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./Je"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./Je"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		char filename2[20];
		sprintf_s(filename2,"Je%d", t);
		ofstream fout(filename2);
		fout<<"e-x e-y e-z x y z"<<endl;
		
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			//if(r[A_Y]<=0 &&  r[A_Y]>-0.001 && -0.02<r[A_X] && r[A_X]<0.02 && ELEM[i].material==FLUID)
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				//if(r[A_Z]<=h+le*0.5 &&  r[A_Z]>=h-le*0.5)//溶融金属の中心からの断面
				{
				fout<<Je[A_X][i]*times<<" "<<Je[A_Y][i]*times<<" "<<Je[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
				}
			}
		}
		fout.close();
	}
	cout<<"ok"<<endl;
	/////*/

	if(CON->get_model_number()==25)
	{
		//渦電流理論解
		//cout<<"全電流の理論解計算"<<endl;//高橋「三次元有限要素法」p.89 100*100*10の銅版
		double Ir=0; //全渦電流
		double sum=0;
		for(int i=1;i<=100000;i++)//Σ内部の無限和を計算
		{
			sum+=pow(-1.0,i)/(pow((2*i+1),3.0)*cosh(2*i+1)*PI*0.5);	
		}
		Ir=-0.5*CON->get_uniform_B()*(-2*PI*CON->get_Hz()*1*CON->get_ele_conduc2()*(0.1*0.5)*(0.1*0.5)*0.01*(1-32*sum/(pow(PI,3.0))));//*sin(2*PI*CON->get_Hz()*TIME)
		cout<<"全渦電流(理論解)[A]="<<Ir<<endl;

		//解析解
		Ir=0;
		double Ir2=0; //全渦電流
		for(int i=1;i<=nelm;i++)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			if(ELEM[i].material==FLUID || ELEM[i].material==CRUCIBLE)
			{
				if(-0.01<r[A_X]<0.01 && r[A_Y]>=0)
				{
					Ir+=abs(Je[A_X][i]);
					Ir2+=sqrt(Je[A_X][i]*Je[A_X][i]+Je[A_Y][i]*Je[A_Y][i]+Je[A_Z][i]*Je[A_Z][i]);
				}
			}
		}
		Ir*=0.05*0.01/2;
		Ir2*=0.05*0.01/2;
		cout<<"全渦電流(解析解,X成分のみ)[A]="<<Ir<<"　(解析解,全成分)[A]="<<Ir2<<endl;
	}


	
	//fout.close();
	for(int D=0;D<DIMENTION;D++) delete [] Jec[D];
	delete [] count_e;

}

//メッシュスムージング
void laplacian_smoothing(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int node0,int nelm0)
{
	cout<<"ラプラシアンスムージング "<<*node<<endl;
	int *flag=new int[*node+1];
	for(int i=0;i<=*node;i++) flag[i]=OFF;
	int *flag2=new int[*nelm+1];
	for(int j=0;j<=*nelm;j++) flag2[j]=OFF;
	int count=0;

	for(int i=node0;i<=*node;i++)
	{
		
		if(NODE[i].material==AIR && NODE[i].BD_node==OFF)//流体節点でもリメッシュ境界節点でもない節点だけスムージング
		{
			for(int D=0;D<3;D++) NODE[i].r[D]=0;
			for(int j=1;j<=*nelm;j++)
			{
				for(int k=1;k<=4;k++) if(ELEM[j].node[k]==i) flag2[j]=ON; //注目している節点を用いて構成される要素ならON

				if(flag2[j]==ON) for(int k=1;k<=4;k++) flag[ELEM[j].node[k]]=ON;
				flag2[j]=OFF;
			}
			flag[i]=OFF; //注目している点もflagがONになっているが、欲しいのはまわりの点なのでflagをOFF
				
			for(int l=1;l<=*node;l++)
			{
				if(flag[l]==ON)
				{
					for(int D=0;D<3;D++) NODE[i].r[D]+=NODE[l].r[D];
					count++;
				}
			}
				
			for(int D=0;D<3;D++) NODE[i].r[D]/=count;
			count=0;
			for(int l=1;l<=*node;l++) flag[l]=OFF;
		}
	}
    delete [] flag;
	delete [] flag2;
}

//メッシュスムージング
void laplacian_smoothing2(vector<point3D> &NODE,vector<element3D> &ELEM,int *node,int *nelm,mpsconfig *CON,int *jnb,int **nei,int node0)
{
	//cout<<"ラプラシアンスムージング"<<endl;
	int *flag=new int[*node+1];
	for(int i=0;i<=*node;i++) flag[i]=OFF;
	int *flag2=new int[*nelm+1];
	for(int j=0;j<=*nelm;j++) flag2[j]=OFF;
	int count=0;
	double *r[3];
	for(int D=0;D<3;D++) r[D]=new double [*node+1];

	for(int i=node0;i<=*node;i++)
	{
		for(int D=0;D<3;D++) r[D][i]=NODE[i].r[D];
		if(NODE[i].material==AIR)//流体節点でもリメッシュ境界節点でもない節点だけスムージング。node0以降はリメッシュ領域に投入した節点群のはず
		{
			for(int D=0;D<3;D++) r[D][i]=0;
			for(int j=1;j<=jnb[i];j++)
			{
				int jelm=nei[i][j];//節点iに隣接する要素番号
				for(int k=1;k<=4;k++) if(ELEM[jelm].node[k]==i) flag2[jelm]=ON; //注目している節点を用いて構成される要素ならON

				if(flag2[jelm]==ON) for(int k=1;k<=4;k++) flag[ELEM[jelm].node[k]]=ON;
				flag2[jelm]=OFF;
			}
			flag[i]=OFF; //注目している点もflagがONになっているが、欲しいのはまわりの点なのでflagをOFF
				
			for(int l=1;l<=*node;l++)
			{
				if(flag[l]==ON)
				{
					for(int D=0;D<3;D++) r[D][i]+=NODE[l].r[D];
					count++;
					
				}
			}
				
			for(int D=0;D<3;D++) r[D][i]/=count;
			count=0;
			for(int l=1;l<=*node;l++) flag[l]=OFF;
		}
	}
	for(int i=node0;i<=*node;i++) for(int D=0;D<3;D++) NODE[i].r[D]=r[D][i];
    delete [] flag;
	delete [] flag2;
	for(int D=0;D<3;D++) delete [] r[D];
}

//節点番号の並び替え関数
void node_sorting(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei)
{
	unsigned int timeA=GetTickCount();
	cout<<"節点番号並び替え"<<endl;

	////節点-節点関係
	vector<vector<int>> NEI2;//節点-節点関係 NEI2[i][j]　節点番号iに、j番目に結合している節点番号が格納される
	NEI2.resize(node+1);
	int *temp=new int[node+1];//チェック用配列

	for(int i=1;i<=node;i++)//注目している節点に繋がっている節点を調べる
	{
		//cout<<i<<endl;
		for(int k=1;k<=jnb[i];k++)
		{
			
			int jelm=nei[i][k];//節点iに隣接する要素番号
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p!=i && temp[p]==0)
				{
					NEI2[i].push_back(p);
					temp[p]=1;//もう見た
				}
			}
		}
		//for(int k=0;k<=node;k++) temp[k]=0;//初期化
		for(int k=0;k<NEI2[i].size();k++) temp[NEI2[i][k]]=0;//初期化
	}
	cout<<"節点-節点関係作成"<<endl;

	int nei2max=0;
	for(int i=1;i<=node;i++) if(NEI2[i].size()>nei2max) nei2max=(int) NEI2[i].size();
	cout<<"最大結合数="<<nei2max;

	/*//
	int nei2min=node;
	for(int i=1;i<=node;i++) if(NEI2[i].size()<nei2min) nei2min=NEI2[i].size();
	int L1ID=0;
	int L1flag=OFF;
	for(int i=1;i<=node;i++)
	{
		if(NEI2[i].size()==nei2min)
		{
			if(L1flag==OFF)
			{
				L1ID=i;
				L1flag=ON;
			}
		}
	}
	//*/
	
	////レベル設定
	int *L=new int[node+1];//節点のレベル
	for(int i=0;i<=node;i++) L[i]=0;	

	int count=node;
	for(int i=1;i<=node;i++) if(jnb[i]==0) count--;//jnb=0の節点はレベルがずっと0なので、レベルセットする節点数から除外する

	L[1]=1; //節点番号1の節点をレベル1とする
	//L[L1ID]=1;
	count--;

	int level=1;
	while(count!=0)
	{
		for(int i=1;i<=node;i++)
		{
			if(L[i]==level)
			{
				for(int k=0;k<NEI2[i].size();k++)
				{
					if(L[NEI2[i][k]]==0)//隣の節点のレベルがまだ決まってない
					{
						L[NEI2[i][k]]=L[i]+1;
						count--;
					}
				}
			}
		}
		level++;
	}

	//レベルが上の節点との結合数を求める（ソート用）
	int *nei_U=new int[node+1];//自分より上のレベルの節点との結合数
	for(int i=0;i<=node;i++) nei_U[i]=0;
	for(int i=1;i<=node;i++)
	{
		for(int k=0;k<NEI2[i].size();k++)
		{
			if(L[i]<L[NEI2[i][k]]) nei_U[i]+=1;
		}
	}
	int nei_Umax=0;
	for(int i=1;i<=node;i++) if(nei_Umax<nei_U[i]) nei_Umax=nei_U[i];

	int Lmax=0;
	for(int i=1;i<=node;i++) if(Lmax<L[i]) Lmax=L[i];

	//レベルごとの節点数を求める
	int *N=new int[Lmax+1];//レベルごとの節点数
	for(int i=0;i<=Lmax;i++) N[i]=0;
	for(int i=0;i<=Lmax;i++) for(int j=1;j<=node;j++) if(L[j]==i) N[i]+=1;

	int Nmax=0;
	for(int i=1;i<=Lmax;i++) if(Nmax<N[i]) Nmax=N[i];
	cout<<"Lmax="<<Lmax<<" Nmax="<<Nmax<<endl;

	vector<vector<int>> Ln;//Ln[i][j] i:レベル　j:レベルiに所属する節点番号を小さいものから格納
	Ln.resize(Lmax+1);
	for(int i=0;i<(int)Ln.size();i++) Ln[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++)
	{
		int k=0;
		for(int j=1;j<=node;j++)
		{
			if(L[j]==i)
			{
				k++;
				Ln[i][k]=j;		
			}
		}
	}

	/*//チェック
	for(int i=1;i<=node;i++)
	{
		int countL=0;
		int countR=0;
		if(L[i]==0 && jnb[i]>0) cout<<"レベルセットされていない節点あり i="<<i<<endl;
		for(int k=0;k<NEI2[i].size();k++)
		{
			//if(L[NEI2[i][k]]==L[i]) cout<<"同じレベルの節点が隣接 i="<<i<<endl;
			if(L[NEI2[i][k]]==L[i]) countL++;
			if(L[NEI2[i][k]]>=L[i]) countR++;
			if(L[NEI2[i][k]]>L[i]+1 ||L[NEI2[i][k]]<L[i]-1) cout<<"2以上離れたレベルの節点が隣接 i="<<i<<endl;
		}
		if(countL==NEI2[i].size() && jnb[i]>0) cout<<"同じレベルの節点のみ隣接 i="<<i<<" L[i]="<<L[i]<<endl;
		if(countR==NEI2[i].size() && jnb[i]>0 && L[i]>1) cout<<"下位レベルの節点に隣接していない i="<<i<<" L[i]="<<L[i]<<endl;

	}
	//*/

	//並び替え開始
	cout<<"並び替え開始"<<endl;
	int *flag=new int[node+1];
	int *flag2=new int[node+1];
	vector<vector<int>> Ln2;//新しい節点番号に基づいて作ったLn
	Ln2.resize(Lmax+1);
	for(int i=0;i<(int)Ln2.size();i++) Ln2[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++) for(int j=0;j<=N[i];j++) Ln2[i][j]=0;
	
	Ln2[1][1]=1;
	//Ln2[1][1]=L1ID;

	int sum=1;
	for(int w=0;w<=node;w++) flag2[w]=OFF;
	//int i=2;
	for(int i=2;i<=Lmax;i++)
	{
		sum+=N[i];
		//cout<<"レベル"<<i<<"に属する節点数"<<N[i]<<endl;
		int endj=1;
		int k=1;
		int swap=0;
		
		while(endj<=N[i])
		{
			for(int iflag=0;iflag<=node;iflag++) flag[iflag]=OFF;
	
			int num=0;
			for(int m=1;m<=N[i-1];m++) if(Ln[i-1][m]==Ln2[i-1][k]) swap=m;

			//レベルi-1のk番目の節点に結合している節点を選び出す num:選んだ節点の数
			for(int j=1;j<=N[i];j++)
			{			
				if(flag2[Ln[i][j]]==OFF)
				{
					for(int knei=0;knei<NEI2[Ln[i][j]].size();knei++)
					{
						//if(NEI2[Ln[i][j]][knei]==Ln[i-1][k])
						if(NEI2[Ln[i][j]][knei]==Ln[i-1][swap])//並び替えた後のものにkを代入して比較する
						{
							flag[Ln[i][j]]=ON;//今回のループで並び替える節点はこのフラグがONになる
							num++;
							flag2[Ln[i][j]]=ON;//このフラグがONになった節点は次のループで並び替えの対象としない
						}
					}
				}
			}
			int end_old=endj;
			endj+=num;
			//cout<<"num="<<num<<endl;			
			//cout<<"endj_old="<<end_old<<endl;
			//cout<<"N[i-1]="<<N[i-1]<<" k="<<k<<" endj="<<endj<<" N[i]="<<N[i]<<endl;

			//選んだ節点をLn[i][j]～Ln[i][endj+num-1]に割り当てる
			
			int je=end_old;

			/*//
			for(int j=1;j<=node;j++)
			{
				if(flag[j]==ON)
				{
					Ln2[j][je]=n;
				}
			}
			//*/
			
			//割り当てた節点を、結合数をもとにソートする（結合数が小さいものから番号付け) 
			//結合数が同じ場合、節点番号が小さいほうを先にする
			//しないほうがいい？
			///
			while(je<endj)
			{
				int nei_Umin=nei_Umax;
				for(int n=1;n<=node;n++) if(flag[n]==ON) if(nei_U[n]<nei_Umin) nei_Umin=nei_U[n];
				int flag_Umin=OFF;
				for(int n=1;n<=node;n++)
				{
					if(nei_U[n]==nei_Umin)
					{
						if(flag_Umin==OFF)//最小の結合数を持つ節点が複数あるときは早いもの勝ちにする
						{
							if(flag[n]==ON)
							{
								Ln2[i][je]=n;
								flag_Umin=ON;
								flag[n]=OFF;
							}
						}
					}
				}
				je++;
			}
			//*/
			k++;
		}
		//cout<<"endj="<<endj<<" N[i]="<<N[i]<<endl;
	}

	/*//
	cout<<"Ln"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln[2][i]<<" "<<nei_U[Ln[2][i]]<<endl;
	}
	cout<<"Ln2"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln2[2][i]<<" "<<nei_U[Ln2[2][i]]<<endl;
	}
	//*/

	cout<<"新節点番号の算出完了"<<endl;

	//新しい番号をつけていく
	//レベルの低いもの、ソートしたレベルごとの配列の頭のほうから順につける
	int *newID=new int[node+1];//新しい節点番号　並び替え前の節点番号がi
	for(int i=0;i<=node;i++) newID[i]=0;
	int countID=0;
	for(int i=1;i<=Lmax;i++)
	{
		for(int j=1;j<=N[i];j++)
		{
			countID++;
			newID[Ln2[i][j]]=countID;
		}
	}
	//jnb=0の節点は末尾につける
	for(int j=1;j<=N[0];j++)
	{
		countID++;
		newID[Ln[0][j]]=countID;
	}

	//節点に関係する各情報の修正
	cout<<"節点番号変更に伴う情報の修正開始"<<endl;
	modify_node_info(CON, node, nelm,NODE,ELEM,jnb,nei,newID);

	delete [] L;
	delete [] flag;
	delete [] flag2;
	delete [] temp;
	delete [] N;
	delete []newID;

	cout<<"並び替え完了－－time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//節点番号の並び替え関数 改良版をここに？
void node_sorting2(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei)
{
	unsigned int timeA=GetTickCount();
	cout<<"節点番号並び替え"<<endl;

	////節点-節点関係
	vector<vector<int>> NEI2;//節点-節点関係 NEI2[i][j]　節点番号iに、j番目に結合している節点番号が格納される
	NEI2.resize(node+1);
	int *temp=new int[node+1];//チェック用配列

	for(int i=1;i<=node;i++)//注目している節点に繋がっている節点を調べる
	{
		//cout<<i<<endl;
		for(int k=1;k<=jnb[i];k++)
		{
			
			int jelm=nei[i][k];//節点iに隣接する要素番号
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];
				if(p!=i && temp[p]==0)
				{
					NEI2[i].push_back(p);
					temp[p]=1;//もう見た
				}
			}
		}
		//for(int k=0;k<=node;k++) temp[k]=0;//初期化
		for(int k=0;k<NEI2[i].size();k++) temp[NEI2[i][k]]=0;//初期化
	}
	cout<<"節点-節点関係作成"<<endl;

	

	/*//
	int nei2min=node;
	for(int i=1;i<=node;i++) if(NEI2[i].size()<nei2min) nei2min=NEI2[i].size();
	int L1ID=0;
	int L1flag=OFF;
	for(int i=1;i<=node;i++)
	{
		if(NEI2[i].size()==nei2min)
		{
			if(L1flag==OFF)
			{
				L1ID=i;
				L1flag=ON;
			}
		}
	}
	//*/
	
	////レベル設定
	int *L=new int[node+1];//節点のレベル
	for(int i=0;i<=node;i++) L[i]=0;	

	int count=node;
	for(int i=1;i<=node;i++) if(jnb[i]==0) count--;//jnb=0の節点はレベルがずっと0なので、レベルセットする節点数から除外する

	L[1]=1; //節点番号1の節点をレベル1とする
	//L[L1ID]=1;
	count--;

	int level=1;
	while(count!=0)
	{
		for(int i=1;i<=node;i++)
		{
			if(L[i]==level)
			{
				for(int k=0;k<NEI2[i].size();k++)
				{
					if(L[NEI2[i][k]]==0)//隣の節点のレベルがまだ決まってない
					{
						L[NEI2[i][k]]=L[i]+1;
						count--;
					}
				}
			}
		}
		level++;
	}

	//レベルが上の節点との結合数を求める（ソート用）
	int *nei_U=new int[node+1];//自分より上のレベルの節点との結合数
	for(int i=0;i<=node;i++) nei_U[i]=0;
	for(int i=1;i<=node;i++)
	{
		for(int k=0;k<NEI2[i].size();k++)
		{
			if(L[i]<L[NEI2[i][k]]) nei_U[i]+=1;
		}
	}
	int nei_Umax=0;
	for(int i=1;i<=node;i++) if(nei_Umax<nei_U[i]) nei_Umax=nei_U[i];

	int Lmax=0;
	for(int i=1;i<=node;i++) if(Lmax<L[i]) Lmax=L[i];

	//レベルごとの節点数を求める
	int *N=new int[Lmax+1];//レベルごとの節点数
	for(int i=0;i<=Lmax;i++) N[i]=0;
	for(int i=0;i<=Lmax;i++) for(int j=1;j<=node;j++) if(L[j]==i) N[i]+=1;

	int Nmax=0;
	for(int i=1;i<=Lmax;i++) if(Nmax<N[i]) Nmax=N[i];
	cout<<"Lmax="<<Lmax<<" Nmax="<<Nmax<<endl;

	vector<vector<int>> Ln;//Ln[i][j] i:レベル　j:レベルiに所属する節点番号を小さいものから格納
	Ln.resize(Lmax+1);
	for(int i=0;i<(int)Ln.size();i++) Ln[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++)
	{
		int k=0;
		for(int j=1;j<=node;j++)
		{
			if(L[j]==i)
			{
				k++;
				Ln[i][k]=j;		
			}
		}
	}

	/*//チェック
	for(int i=1;i<=node;i++)
	{
		int countL=0;
		int countR=0;
		if(L[i]==0 && jnb[i]>0) cout<<"レベルセットされていない節点あり i="<<i<<endl;
		for(int k=0;k<NEI2[i].size();k++)
		{
			//if(L[NEI2[i][k]]==L[i]) cout<<"同じレベルの節点が隣接 i="<<i<<endl;
			if(L[NEI2[i][k]]==L[i]) countL++;
			if(L[NEI2[i][k]]>=L[i]) countR++;
			if(L[NEI2[i][k]]>L[i]+1 ||L[NEI2[i][k]]<L[i]-1) cout<<"2以上離れたレベルの節点が隣接 i="<<i<<endl;
		}
		if(countL==NEI2[i].size() && jnb[i]>0) cout<<"同じレベルの節点のみ隣接 i="<<i<<" L[i]="<<L[i]<<endl;
		if(countR==NEI2[i].size() && jnb[i]>0 && L[i]>1) cout<<"下位レベルの節点に隣接していない i="<<i<<" L[i]="<<L[i]<<endl;

	}
	//*/

	//並び替え開始

	cout<<"並び替え開始"<<endl;
	int *flag=new int[node+1];
	int *flag2=new int[node+1];
	vector<vector<int>> Ln2;//新しい節点番号に基づいて作ったLn
	Ln2.resize(Lmax+1);
	for(int i=0;i<(int)Ln2.size();i++) Ln2[i].resize(N[i]+1);
	for(int i=0;i<=Lmax;i++) for(int j=0;j<=N[i];j++) Ln2[i][j]=0;
	
	Ln2[1][1]=1;
	//Ln2[1][1]=L1ID;

	int sum=1;
	for(int w=0;w<=node;w++) flag2[w]=OFF;
	//int i=2;
	for(int i=2;i<=Lmax;i++)
	{
		sum+=N[i];
		cout<<"レベル"<<i<<"に属する節点数"<<N[i]<<endl;
		int endj=1;
		int k=1;
		int swap=0;
		
		while(endj<=N[i])
		{
			for(int iflag=0;iflag<=node;iflag++) flag[iflag]=OFF;
	
			int num=0;
			for(int m=1;m<=N[i-1];m++) if(Ln[i-1][m]==Ln2[i-1][k]) swap=m;

			//レベルi-1のk番目の節点に結合している節点を選び出す num:選んだ節点の数
			for(int j=1;j<=N[i];j++)
			{			
				if(flag2[Ln[i][j]]==OFF)
				{
					for(int knei=0;knei<NEI2[Ln[i][j]].size();knei++)
					{
						//if(NEI2[Ln[i][j]][knei]==Ln[i-1][k])
						if(NEI2[Ln[i][j]][knei]==Ln[i-1][swap])//並び替えた後のものにkを代入して比較する
						{
							flag[Ln[i][j]]=ON;//今回のループで並び替える節点はこのフラグがONになる
							num++;
							flag2[Ln[i][j]]=ON;//このフラグがONになった節点は次のループで並び替えの対象としない
						}
					}
				}
			}
			int end_old=endj;
			endj+=num;
			//cout<<"num="<<num<<endl;			
			//cout<<"endj_old="<<end_old<<endl;
			//cout<<"N[i-1]="<<N[i-1]<<" k="<<k<<" endj="<<endj<<" N[i]="<<N[i]<<endl;

			//選んだ節点をLn[i][j]～Ln[i][endj+num-1]に割り当てる

			int je=end_old;
			for(int n=1;n<=node;n++)
			{
				if(flag[n]==ON)
				{
					Ln2[i][je]=n;
					je++;
				}
			}
			k++;
		}
		//cout<<"endj="<<endj<<" N[i]="<<N[i]<<endl;
	}

	/*//
	cout<<"Ln"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln[2][i]<<" "<<nei_U[Ln[2][i]]<<endl;
	}
	cout<<"Ln2"<<endl;
	for(int i=1;i<=14;i++)
	{
		cout<<Ln2[2][i]<<" "<<nei_U[Ln2[2][i]]<<endl;
	}
	//*/

	cout<<"新節点番号の算出完了"<<endl;

	//新しい番号をつけていく
	//レベルの低いもの、ソートしたレベルごとの配列の頭のほうから順につける
	int *newID=new int[node+1];//新しい節点番号　並び替え前の節点番号がi
	for(int i=0;i<=node;i++) newID[i]=0;
	int countID=0;
	for(int i=1;i<=Lmax;i++)
	{
		for(int j=1;j<=N[i];j++)
		{
			countID++;
			newID[Ln2[i][j]]=countID;
		}
	}
	//jnb=0の節点は末尾につける
	for(int j=1;j<=N[0];j++)
	{
		countID++;
		newID[Ln[0][j]]=countID;
	}

	//節点に関係する各情報の修正
	cout<<"節点番号変更に伴う情報の修正開始"<<endl;
	modify_node_info(CON, node, nelm,NODE,ELEM,jnb,nei,newID);

	delete [] L;
	delete [] flag;
	delete [] flag2;
	delete [] temp;
	delete [] N;
	delete []newID;

	cout<<"並び替え完了－－time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}


//情報修正関数
void modify_node_info(mpsconfig *CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,int **nei,int *newID)
{
	vector<point3D> N_NODE;
	N_NODE.resize(node+1);

	//NODEクラスの内容を一時配列に記憶
	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) N_NODE[newID[i]].r[D]=NODE[i].r[D];
		N_NODE[newID[i]].material=NODE[i].material;
		N_NODE[newID[i]].boundary_condition=NODE[i].boundary_condition;
		N_NODE[newID[i]].particleID=NODE[i].particleID;
		N_NODE[newID[i]].remesh=NODE[i].remesh;
		N_NODE[newID[i]].BD_node=NODE[i].BD_node;
	}

	for(int i=1;i<=node;i++)
	{
		for(int D=0;D<3;D++) NODE[i].r[D]=N_NODE[i].r[D];
		NODE[i].material=N_NODE[i].material;
		NODE[i].boundary_condition=N_NODE[i].boundary_condition;
		NODE[i].particleID=N_NODE[i].particleID;
		NODE[i].remesh=N_NODE[i].remesh;
		NODE[i].BD_node=N_NODE[i].BD_node;
	}

	//要素は構成する節点番号のみ修正
	for(int i=1;i<=nelm;i++)
	{
		int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];

		ELEM[i].node[1]=newID[ia];
		ELEM[i].node[2]=newID[ib];
		ELEM[i].node[3]=newID[ic];
		ELEM[i].node[4]=newID[ip];
	}

	fill3D(NODE,ELEM,nelm);
	
	//jnb,neiは並び替えが終わったあと求めなおしている

}

//力や各物理量をコンター図で表示する関数
void output_F_scalar_with_AVS_for_linear(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int t,vector<mpsparticle> &PART,int node)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	//t=1;//いまはわざと毎ステップ上書き
	//int number=int (PART.size());
	//int node=int(NODE.size()-1);
	//int nelm=int(ELEM.size());

	//sprintf_s(filename,"pressure/pressure%d",t);//フォルダを作成して管理する場合はこちら

	cout<<"粒子の物理量出力開始"<<endl;

	////磁束密度
	//sprintf_s(filename,"PART.B%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/B%d",t);
	ofstream fout(filename);
	if(!fout)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				if(PART[p].type==FLUID)
				{
				
					double x=PART[p].r[A_X];	//rは非常に小さい値なので10^5倍しておく
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].B;
					fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout.close();


	//sprintf_s(filename,"PART.B%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/B%d.fld",t);
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
	fout2 << "label=Bflux" << endl << endl;
	fout2 << "variable 1 file=B" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    1 file=B" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    2 file=B" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    3 file=B" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2.close();

	////磁束密度断面
	n=0;
	//sprintf_s(filename,"PART.B%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Bsec%d",t);
	ofstream fouts(filename);
	if(!fouts)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				if(PART[p].type==FLUID && PART[p].r[A_Z]<0.13125+CON->get_distancebp() && PART[p].r[A_Z]>0.13125-CON->get_distancebp())
				{
				
					double x=PART[p].r[A_X];	//rは非常に小さい値なので10^5倍しておく
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].B;
					fouts << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fouts.close();


	//sprintf_s(filename,"PART.B%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Bsec%d.fld",t);
	ofstream fout2s(filename);
	if(!fout2s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout2s << "# AVS field file" << endl;
	fout2s << "ndim=1" << endl;
	fout2s << "dim1=" << n <<endl;
	fout2s << "nspace=3" << endl;
	fout2s << "veclen=1" << endl;
	fout2s << "data=float" << endl;
	fout2s << "field=irregular" << endl;
	fout2s << "label=Bflux" << endl << endl;
	fout2s << "variable 1 file=Bsec" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2s << "coord    1 file=Bsec" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2s << "coord    2 file=Bsec" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2s << "coord    3 file=Bsec" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2s.close();

	///渦電流密度
	n=0;
	//sprintf_s(filename,"PART.J%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Jesec%d",t);
	ofstream fout3(filename);
	if(!fout3)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				if(PART[p].type==FLUID)
				
				{
				
					double x=PART[p].r[A_X];	//rは非常に小さい値なので10^5倍しておく
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].Je;
					fout3 << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout3.close();

	//sprintf_s(filename,"PART.J%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Jesec%d.fld",t);
	ofstream fout4(filename);
	if(!fout4)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout4 << "# AVS field file" << endl;
	fout4 << "ndim=1" << endl;
	fout4 << "dim1=" << n <<endl;
	fout4 << "nspace=3" << endl;
	fout4 << "veclen=1" << endl;
	fout4 << "data=float" << endl;
	fout4 << "field=irregular" << endl;
	fout4 << "label=eddy " << endl << endl;
	fout4 << "variable 1 file=Je" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4 << "coord    1 file=Je" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4 << "coord    2 file=Je" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4 << "coord    3 file=Je" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4.close();

	///渦電流密度 断面
	n=0;
	//sprintf_s(filename,"PART.J%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Jesec%d",t);
	ofstream fout3s(filename);
	if(!fout3s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID && PART[p].surface==ON)
				//if(PART[p].type==FLUID)
				if(PART[p].type==FLUID && PART[p].r[A_Z]<0.13125+CON->get_distancebp() && PART[p].r[A_Z]>0.13125-CON->get_distancebp())
				{
				
					double x=PART[p].r[A_X];	//rは非常に小さい値なので10^5倍しておく
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=NODE[i].Je;
					fout3s << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout3s.close();

	//sprintf_s(filename,"PART.J%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Jesec%d.fld",t);
	ofstream fout4s(filename);
	if(!fout4s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout4s << "# AVS field file" << endl;
	fout4s << "ndim=1" << endl;
	fout4s << "dim1=" << n <<endl;
	fout4s << "nspace=3" << endl;
	fout4s << "veclen=1" << endl;
	fout4s << "data=float" << endl;
	fout4s << "field=irregular" << endl;
	fout4s << "label=eddy " << endl << endl;
	fout4s << "variable 1 file=Jesec" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4s << "coord    1 file=Jesec" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4s << "coord    2 file=Jesec" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4s << "coord    3 file=Jesec" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout4s.close();

	///電磁力
	n=0;
	//sprintf_s(filename,"PART.F%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/F%d",t);
	ofstream fout5(filename);
	if(!fout5)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				if(PART[p].type==FLUID)
				{
				
					double x=PART[p].r[A_X];	//rは非常に小さい値なので10^5倍しておく
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=sqrt(PART[p].F[A_X]*PART[p].F[A_X]+PART[p].F[A_Y]*PART[p].F[A_Y]+PART[p].F[A_Z]*PART[p].F[A_Z]);
					fout5 << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout5.close();

	//sprintf_s(filename,"PART.F%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/F%d.fld",t);
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
	fout6 << "label=lorentz" << endl << endl;
	fout6 << "variable 1 file=F" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6 << "coord    1 file=F" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6 << "coord    2 file=F" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6 << "coord    3 file=F" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6.close();

	///電磁力断面
	n=0;
	//sprintf_s(filename,"PART.F%d",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Fsec%d",t);
	ofstream fout5s(filename);
	if(!fout5s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}
	if(CON->get_dimention()==3)
	{
		for(int i=1;i<=node;i++)
		{
			int p=NODE[i].particleID;
			if(p>=0)
			{
				//if(p==1565024) cout<<i<<" "<<NODE[i].material<<endl;
				//cout<<p<<endl;
				//if(PART[p].type==FLUID)
				if(PART[p].type==FLUID && PART[p].r[A_Z]<0.13125+CON->get_distancebp() && PART[p].r[A_Z]>0.13125-CON->get_distancebp())
				{
				
					double x=PART[p].r[A_X];	//rは非常に小さい値なので10^5倍しておく
					double y=PART[p].r[A_Y];//*1.0E+05;
					double z=PART[p].r[A_Z];//*1.0E+05;
					double P=sqrt(PART[p].F[A_X]*PART[p].F[A_X]+PART[p].F[A_Y]*PART[p].F[A_Y]+PART[p].F[A_Z]*PART[p].F[A_Z]);
					fout5s << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	fout5s.close();

	//sprintf_s(filename,"PART.F%d.fld",t);//他のファイルと同じ階層に生成するならこちら
	sprintf_s(filename,"plot_scalar/Fsec%d.fld",t);
	ofstream fout6s(filename);
	if(!fout6s)
	{
		cout << "cannot open" << filename << endl;
		exit(EXIT_FAILURE);
	}

	fout6s << "# AVS field file" << endl;
	fout6s << "ndim=1" << endl;
	fout6s << "dim1=" << n <<endl;
	fout6s << "nspace=3" << endl;
	fout6s << "veclen=1" << endl;
	fout6s << "data=float" << endl;
	fout6s << "field=irregular" << endl;
	fout6s << "label=lorentz" << endl << endl;
	fout6s << "variable 1 file=Fsec" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6s << "coord    1 file=Fsec" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6s << "coord    2 file=Fsec" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6s << "coord    3 file=Fsec" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout6s.close();

	///
	/*
	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	fp<<"6 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	//fp<<"2 1 1"<<endl;	//この行の詳細はヘルプを参照
	fp<<"6 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	fp<<"speed,m/s"<<endl;
	fp<<"surface_tension,N/m^3"<<endl;
	fp<<"B,T"<<endl;
	fp<<"Je,A/m^2"<<endl;
	//fp<<"Et,V/m"<<endl;
	fp<<"Fn,N/m^3"<<endl;
	//fp<<"P,N/m^2"<<endl;
	fp<<"value1,??"<<endl;

	

	//各節点の情報値入力
	for(int i=0;i<node;i++)
	{
		int p=NODE[i].particleID;
		double speed=0;
		double F=0;
		double vn=0;
		double potential=0;
		double Pst=0;//
		
		//速度計算
		if(p>=0)
		{
			for(int D=0;D<3;D++)
			{
				speed+=PART[p].u[D]*PART[p].u[D];
				F+=PART[p].F[D]*PART[p].F[D];
				potential+=PART[p].potential[D]*PART[p].potential[D];
			}
			Pst=PART[p].dir_Pst;
		}
		speed=sqrt(speed);
		F=sqrt(F)*CON->get_density()/CON->get_particle_mass();
		potential=sqrt(F)*CON->get_density();
		////

		//double P=Pst-NODE[i].Fn;
		if(p>=0) fp<<count<<" "<<speed<<" "<<potential<<" "<<NODE[i].B<<" "<<NODE[i].Je<<" "<<NODE[i].F<<" "<<NODE[i].value1<<endl;
	}
	*/
	
	/*fp<<"3 1 1 1"<<endl;
	fp<<"potential,V"<<endl;
	fp<<"En,V/m"<<endl;
	fp<<"Fn,N/m^2"<<endl;
	
	//要素情報出力　要素情報はmicroAVSの可視化メソッドバー内の、「要素データの塗りつぶし」を押せば見れる
	for(int i=0;i<nelm;i++) fp<<i+1<<"  "<<ELEM[i].potential<<" "<<ELEM[i].En<<" "<<ELEM[i].Fn<<endl;
	*/
	//cout<<"OK"<<endl;
	//fp.close();
}

//力などの物理量をコンター図で表示する関数
void output_F_scalar_movie_with_AVS_for_linear(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int t,vector<mpsparticle> &PART,double TIME)
{
	/*
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	cout<<"Now F_scalar movie writing-----";
	int nelm=int (ELEM.size());		//要素数
	int node=int (NODE.size());		//節点数
	int STEP=CON->get_step()/(CON->get_EM_interval()*CON->get_mesh_output_interval())+1;				//出力する総ステップ数
	
	if(t==1) 
	{
		ofstream fp("F_scalar_movie.inp");
		fp<<STEP<<endl;
		fp<<"data_geom"<<endl;
		fp.close();

		ofstream fq("step_for_F_scalar_movie.dat");			//この関数が呼び出された回数をファイルに記憶し、以後、関数が呼び出されるたびにカウントを増やしていく。この値はstepで使用する。
		fq<<1<<endl;
		fq.close();
	}

	//step_nowの値をファイルより決定
	int step_now;											//この関数が呼び出された回数
	ifstream fin("step_for_F_scalar_movie.dat");
	if(!fin) cout<<"cannot open step_for_F_scalar_movie.dat"<<endl;
	fin>>step_now;
	fin.close();

	//mainファイル書き込み
	ofstream fp("F_scalar_movie.inp",ios :: app);
	fp<<"step"<<step_now<<" TIME="<<TIME<<endl;
	fp<<node<<" "<<nelm<<endl;	//節点数と要素数出力
	
	//節点番号とその座標の出力 
	for(int i=0;i<node;i++) fp<<i+1<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=0;i<nelm;i++)
	{
		fp<<i+1<<"  0 tri ";
		for(int j=0;j<3;j++)	fp<<ELEM[i].node[j]+1<<" ";
		fp<<endl;
	}

	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	fp<<"8 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	//fp<<"2 1 1"<<endl;	//この行の詳細はヘルプを参照
	fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	fp<<"speed,m/s"<<endl;
	fp<<"surface_tension,N/m^2"<<endl;
	fp<<"potential,V"<<endl;
	fp<<"En,V/m"<<endl;
	fp<<"Et,V/m"<<endl;
	fp<<"Fn,N/m^2"<<endl;
	fp<<"P,N/m^2"<<endl;
	fp<<"value1,??"<<endl;

	

	//各節点の情報値入力
	for(int i=0;i<node;i++)
	{
		int p=NODE[i].particle;
		double speed=0;
		double vn=0;
		double Pst=0;//
		if(p>=0)
		{
			for(int D=0;D<3;D++)
			{
				speed+=PART[p].u[D]*PART[p].u[D];
			}
			Pst=PART[p].dir_Pst;
		}
		speed=sqrt(speed);
		double P=Pst-NODE[i].Fn;
		fp<<i+1<<" "<<speed<<" "<<Pst<<" "<<NODE[i].potential<<" "<<NODE[i].slop1<<" "<<NODE[i].Et<<" "<<NODE[i].Fn<<" "<<P<<" "<<NODE[i].value1<<endl;
	}

	/*fp<<"3 1 1 1"<<endl;
	fp<<"potential,V"<<endl;
	fp<<"En,V/m"<<endl;
	fp<<"Fn,N/m^2"<<endl;
	
	//要素情報出力　要素情報はmicroAVSの可視化メソッドバー内の、「要素データの塗りつぶし」を押せば見れる
	for(int i=0;i<nelm;i++) fp<<i+1<<"  "<<ELEM[i].potential<<" "<<ELEM[i].En<<" "<<ELEM[i].Fn<<endl;
	/
	fp.close();

	step_now++;											//呼び出し回数を次回に向けて＋＋しておく
	ofstream fq("step_for_F_scalar_movie.dat");			//その値をファイル出力
	fq<<step_now<<endl;
	fq.close();

	cout<<"OK"<<endl;
	*/
	

	
}