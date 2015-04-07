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


#define LOOP 0
#define UNLOOP 1

void set_BEM_model(mpsconfig *CON,vector<point2D> &static_NODE,vector<element2D> &static_ELEM,int *s_node_num,int *s_elem_num,vector<mpsparticle> &PART,int fluid_number,int particle_number)
{
	//静的節点及び要素の情報生成
	set_BEM_static_model(CON,static_NODE,static_ELEM,s_node_num,s_elem_num);
}

void set_BEM_static_model(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,int *s_node_num,int *s_elem_num)
{
	//流体以外の節点・要素情報を作成

	int node_num=0;	//節点数
	int elemnum=0;	//要素数
	int ele_type=CON->get_BEM_elm_type();	//要素タイプ 0:一定 1:線形
	int B_flag=UNLOOP;		//B_flagがLOOPなら境界は繋がっている。UNLOOPならつながってない
	int divN;	//1辺の分割数
	double divL;	//分割幅

	point2D NODE01;
	element2D ELEM01;

	double Xmin=CON->get_minX();//解析領域
	double Xmax=CON->get_maxX();
	double Ymin=CON->get_minY();
	double Ymax=CON->get_maxY();

	double V1=0;//CON->get_V();
	double V2=CON->get_V();		//ポテンシャル値

	//下辺作成
	divN=50;
	divL=(Xmax-Xmin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmin+n*divL;
		NODE[node_num].r[A_Y]=Ymin;
		NODE[node_num].potential=V1;
		NODE[node_num].slop1=0;						//初期化
		NODE[node_num].boundary_condition=Diric;
		if(n==0) NODE[node_num].C=0.25;				//内角
		else NODE[node_num].C=0.5;
		NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
		node_num++;
	}
	//右辺作成
	divL=(Ymax-Ymin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmax;
		NODE[node_num].r[A_Y]=n*divL+Ymin;
		NODE[node_num].potential=0;
		NODE[node_num].slop1=0;						//初期化
		NODE[node_num].boundary_condition=Neumn;
		NODE[node_num].C=0.5;
		if(n==0)
		{
			NODE[node_num].C=0.25;				//内角
			NODE[node_num].boundary_condition=Diric;
		}
		NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
		node_num++;
	}
	//上辺作成
	divL=(Xmax-Xmin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmax-n*divL;
		NODE[node_num].r[A_Y]=Ymax;
		NODE[node_num].potential=V2;
		NODE[node_num].slop1=0;						//初期化
		NODE[node_num].boundary_condition=Diric;
		if(n==0) NODE[node_num].C=0.25;				//内角
		else NODE[node_num].C=0.5;
		NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
		node_num++;
	}
	//左辺作成
	divL=(Ymax-Ymin)/divN;
	for(int n=0;n<divN;n++)
	{
		NODE.push_back(NODE01);
		NODE[node_num].r[A_X]=Xmin;
		NODE[node_num].r[A_Y]=Ymax-n*divL;
		NODE[node_num].potential=0;
		NODE[node_num].slop1=0;						//初期化
		NODE[node_num].boundary_condition=Neumn;
		NODE[node_num].C=0.5;
		if(n==0)
		{
			NODE[node_num].C=0.25;				//内角
			NODE[node_num].boundary_condition=Diric;
			NODE[node_num].potential=V2;
		}
		NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
		node_num++;
	}
	elemnum=node_num;			//2Dでは節点数=要素数
	for(int i=0;i<node_num;i++) NODE[i].particle=-1;	//対応するものがないので－1をダミーとして格納
	B_flag=LOOP;

	//if(B_flag==LOOP)
	{
		//要素情報作成
		for(int i=0;i<elemnum;i++)
		{
			ELEM.push_back(ELEM01);
			ELEM[i].node[0]=i;
			ELEM[i].node[1]=i+1;
			if(i+1==node_num) ELEM[i].node[1]=0;//最後の要素の終点は節点0
			int n1=ELEM[i].node[0];
			int n2=ELEM[i].node[1];
			double X=NODE[n2].r[A_X]-NODE[n1].r[A_X];
			double Y=NODE[n2].r[A_Y]-NODE[n1].r[A_Y];
			double L=sqrt(X*X+Y*Y);			//要素長さ
			ELEM[i].L=L;
			ELEM[i].direct[A_X]=Y/L;
			ELEM[i].direct[A_Y]=-X/L;
			ELEM[i].r[A_X]=0.5*(NODE[n1].r[A_X]+NODE[n2].r[A_X]);
			ELEM[i].r[A_Y]=0.5*(NODE[n1].r[A_Y]+NODE[n2].r[A_Y]);
			ELEM[i].material=WALL;

			if(NODE[n1].boundary_condition==NODE[n2].boundary_condition) ELEM[i].boundary_condition=NODE[n1].boundary_condition;	//両端が同じ境界条件ならそれに習う
			else
			{
				ELEM[i].boundary_condition=Neumn;//両端で境界条件が異なる場合、それはﾃﾞｨﾘｸﾚとﾉｲﾏﾝの堺である。そのときはﾉｲﾏﾝ型とする。
				if(ele_type==CONSTANT)
				{
					int N=n1;
					if(NODE[n2].boundary_condition==Neumn) N=n2;
					ELEM[i].node[0]=N;		//一定要素の場合は、第一節点としてﾉｲﾏﾝ型の方を格納 これはELEM[i]はﾉｲﾏﾝ型で、その法線微分値は節点Nのそれに等しいことを意味する
				}
			}
		}
	}

	cout<<"固定節点数="<<node_num<<" 固定要素数="<<elemnum<<endl;
	
	*s_node_num=node_num;
	*s_elem_num=elemnum;

}

//一定要素用の、動的節点・要素情報生成関数
void set_BEM2D_dynaic_model_for_CONSTANT(mpsconfig *CON,vector<mpsparticle> &PART,vector<point2D> &NODE,vector<element2D> &ELEM,int *d_node_num,int *d_elem_num,int particle_number,int fluid_number)
{
	if(CON->get_dimention()==3) cout<<"set_BEM2D_dynaic_model_for_CONSTANTは3Dでは使用できません"<<endl;

	int node_num=0;
	int elem_num=0;
	double le=CON->get_distancebp();
	point2D NODE01;
	element2D ELEM01;
	
	int CLOSE_FLAG=ON;			//領域がきちんと閉じていたらON　そうでないならOFF

	double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];//内向き法線ベクトル

	int *BEMID=new int[particle_number];		//BEMに節点として出力するならその節点番号を入力　出力しないなら-1を格納

	//法線ベクトル生成
	for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON && PART[i].toBEM==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
	}

	for(int i=0;i<particle_number;i++) BEMID[i]=-1;		//初期化

	//節点情報入力
	for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON && PART[i].toBEM==ON)
		{
			NODE.push_back(NODE01);
			for(int D=0;D<2;D++) NODE[node_num].r[D]=PART[i].r[D];
			//NODE[node_num].boundary_condition=Diric;
			NODE[node_num].boundary_condition=BOTH;		//２媒質のときはこちら
			NODE[node_num].potential=0;					//0?
			NODE[node_num].slop1=0;					//初期化
			NODE[node_num].slop2=0; 				//初期化
			NODE[node_num].C=0.5;					//一定要素は0.5
			NODE[node_num].L=le;					//これでいい？
			NODE[node_num].particle=i;				//対応する粒子番号記憶
			BEMID[i]=node_num;							//BEMに出力するしるし
			node_num++;
		}
	}///

	for(int n=0;n<node_num;n++)
	{
		int i=NODE[n].particle;//対応する粒子
		double nx=direct[A_X][i];				//内向き単位法線ベクトルであることに注意
		double ny=direct[A_Y][i];
		int J=i;								//要素の他端を構成する粒子番号　iで初期化
		double mindis=100;
		for(int k=0;k<PART[i].N3;k++)			//近隣粒子のなかで、法線ベクトルの右側に存在する最短距離の粒子を探す
		{
			int j=PART[i].NEI3[k];
			if(BEMID[j]>=0)		//対応するBEM節点が存在するなら
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double dis=sqrt(X*X+Y*Y);
				double COS=(X*nx+Y*ny)/(dis*1.0);//粒子iの法線ベクトルと、rijのなす角度(cosθ)
				NODE[n].C=acos(COS);				//内角
				//cout<<360*NODE[n].C/(2*PI)<<endl;
				double O_product=X*ny-Y*nx;			//外積　rij×direct
				if(O_product>0)	//外積が正ということは、粒子jは法線ベクトルの右側に存在する
				{
					if(dis<mindis)
					{
						J=j;
						mindis=dis;
					}
				}
			}
		}
		if(J!=i)		//一緒に要素を構成する粒子が見つかったなら
		{
			double X=PART[J].r[A_X]-PART[i].r[A_X];
			double Y=PART[J].r[A_Y]-PART[i].r[A_Y];
			double dis=mindis;
			
			//要素作成
			ELEM.push_back(ELEM01);
			ELEM[elem_num].node[0]=n;
			ELEM[elem_num].node[1]=BEMID[J];
			for(int D=0;D<2;D++) ELEM[elem_num].r[D]=0.5*(PART[i].r[D]+PART[J].r[D]);//2点間の中点
			ELEM[elem_num].L=dis;
			ELEM[elem_num].direct[A_X]=-Y/dis;		//内側の領域なので、解析領域端で定義した法線ベクトルとは定義が反対となる
			ELEM[elem_num].direct[A_Y]=X/dis;
			ELEM[elem_num].map=0;
			ELEM[elem_num].material=FLUID;
			ELEM[elem_num].boundary_condition=NODE[n].boundary_condition;
			
			elem_num++;
		}
		else
		{
			cout<<"領域がとじていない？"<<endl;
			cout<<n<<" "<<i<<endl;
			CLOSE_FLAG=OFF;
		}
	}//*/

	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	cout<<"動的節点数="<<node_num<<" 動的要素数="<<elem_num<<endl;
	*d_node_num=node_num;
	*d_elem_num=elem_num;

	delete [] BEMID;

	if(CLOSE_FLAG==OFF)
	{
		ofstream fl("checkclose1.dat");
		ofstream fl2("checkclose2.dat");
		for(int n=0;n<node_num;n++) fl<<NODE[n].r[A_X]<<" "<<NODE[n].r[A_Y]<<" "<<n<<endl;//全節点位置
		
		for(int n=0;n<elem_num;n++) fl2<<ELEM[n].r[A_X]<<" "<<ELEM[n].r[A_Y]<<" "<<n<<endl;
		fl.close();
		fl2.close();
	}
	
}

void couple_NODE_and_ELEM(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,vector<point2D> &s_NODE,vector<element2D> &s_ELEM,vector<point2D> &dy_NODE,vector<element2D> &dy_ELEM,vector<REGION> &region)
{
	//s_NODEとdy_NODEを合体させてNODEに格納する
	int s_node_num=int (s_NODE.size());
	int dy_node_num=int (dy_NODE.size());
	int s_elem_num=int (s_ELEM.size());
	int dy_elem_num=int (dy_ELEM.size());
	int countN=0;
	int countE=0;

	point2D NODE01;
	element2D ELEM01;
	REGION temp;
	region.push_back(temp);
	region[0].start=0;					//最初の計算領域は必ず計算点０から始まる

	int *s_Nid=new int[s_node_num];		//移動先の節点番号格納。s_NODE[i]だった節点はNODE[s_Nid[i]]に格納される
	int *dy_Nid=new int[dy_node_num];


	//節点情報生成
	for(int i=0;i<s_node_num;i++)
	{
		NODE.push_back(NODE01);
		NODE[countN].r[A_X]=s_NODE[i].r[A_X];
		NODE[countN].r[A_Y]=s_NODE[i].r[A_Y];
		NODE[countN].boundary_condition=s_NODE[i].boundary_condition;
		NODE[countN].potential=s_NODE[i].potential;
		NODE[countN].slop1=s_NODE[i].slop1;
		NODE[countN].slop2=s_NODE[i].slop2;
		NODE[countN].C=s_NODE[i].C;
		NODE[countN].L=s_NODE[i].L;
		NODE[countN].particle=s_NODE[i].particle;
		s_Nid[i]=countN;			//i番目のs_NODEはcountN番目の節点
		countN++;
	}
	region.push_back(temp);//第２領域
	region[1].start=countN;
	for(int i=0;i<dy_node_num;i++)
	{
		NODE.push_back(NODE01);
		NODE[countN].r[A_X]=dy_NODE[i].r[A_X];
		NODE[countN].r[A_Y]=dy_NODE[i].r[A_Y];
		NODE[countN].boundary_condition=dy_NODE[i].boundary_condition;
		NODE[countN].potential=dy_NODE[i].potential;
		NODE[countN].slop1=dy_NODE[i].slop1;
		NODE[countN].slop2=dy_NODE[i].slop2;
		NODE[countN].C=dy_NODE[i].C;
		NODE[countN].L=dy_NODE[i].L;
		NODE[countN].particle=dy_NODE[i].particle;
		dy_Nid[i]=countN;			//i番目のdy_NODEはcountN番目の節点
		countN++;
	}
	region[0].end=countN;
	region[1].end=countN;
	

	//要素情報作成
	for(int i=0;i<s_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		ELEM[countE].node[0]=s_ELEM[i].node[0];
		ELEM[countE].node[1]=s_ELEM[i].node[1];
		ELEM[countE].r[A_X]=s_ELEM[i].r[A_X];
		ELEM[countE].r[A_Y]=s_ELEM[i].r[A_Y];
		ELEM[countE].boundary_condition=s_ELEM[i].boundary_condition;
		ELEM[countE].L=s_ELEM[i].L;
		ELEM[countE].direct[A_X]=s_ELEM[i].direct[A_X];
		ELEM[countE].direct[A_Y]=s_ELEM[i].direct[A_Y];
		ELEM[countE].map=s_ELEM[i].map;
		ELEM[countE].material=s_ELEM[i].material;
		countE++;
	}
	for(int i=0;i<countE;i++)
	{
		int n1=ELEM[i].node[0];//この段階で格納されている節点番号は、s_NODEのものであり、NODEのものではない
		int n2=ELEM[i].node[1];
		n1=s_Nid[n1];			//節点番号をNODE基準のものに書き換え
		n2=s_Nid[n2];
		ELEM[i].node[0]=n1;
		ELEM[i].node[1]=n2;
	}
	for(int i=0;i<dy_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		ELEM[countE].node[0]=dy_ELEM[i].node[0];
		ELEM[countE].node[1]=dy_ELEM[i].node[1];
		ELEM[countE].r[A_X]=dy_ELEM[i].r[A_X];
		ELEM[countE].r[A_Y]=dy_ELEM[i].r[A_Y];
		ELEM[countE].boundary_condition=dy_ELEM[i].boundary_condition;
		ELEM[countE].L=dy_ELEM[i].L;
		ELEM[countE].direct[A_X]=dy_ELEM[i].direct[A_X];
		ELEM[countE].direct[A_Y]=dy_ELEM[i].direct[A_Y];
		ELEM[countE].map=dy_ELEM[i].map;
		ELEM[countE].material=dy_ELEM[i].material;
		countE++;
	}
	for(int i=s_elem_num;i<countE;i++)
	{
		int n1=ELEM[i].node[0];//この段階で格納されている節点番号は、dy_NODEのものであり、NODEのものではない
		int n2=ELEM[i].node[1];
		n1=dy_Nid[n1];			//節点番号をNODE基準のものに書き換え
		n2=dy_Nid[n2];
		ELEM[i].node[0]=n1;
		ELEM[i].node[1]=n2;
		//cout<<n1<<" "<<n2<<endl;
	}
	
	cout<<"全節点数="<<countN<<" 全要素数="<<countE<<endl;

	//ファイル出力して座標を確認
	ofstream fp("BEM_input.dat");
	double times=0.0001;
	for(int i=0;i<countE;i++)
	{
		int n=ELEM[i].node[0];
	//	fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<endl;
		fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<"\t"<<ELEM[i].direct[A_X]*times<<"\t"<<ELEM[i].direct[A_Y]*times<<endl;
		//fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<"\t"<<NODE[n].potential<<endl;
	}//*/
	//for(int i=0;i<countN;i++) fp<<NODE[i].r[A_X]<<"\t"<<NODE[i].r[A_Y]<<"\t"<<NODE[i].potential<<endl;
	fp.close();


	delete [] s_Nid;	
	delete [] dy_Nid;
}


