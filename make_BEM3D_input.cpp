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

//直方体解析領域作成関数
void make_cube_region(mpsconfig *CON,vector<BEMpoint3D> &NODE,int *s_node_num);


void set_BEM3D_static_model(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,int *s_node_num,int *s_elem_num)
{
	//流体以外の節点・要素情報を作成

	int node_num=0;								//節点数
	int elemnum=0;								//要素数
	int BEM_elm_type=CON->get_BEM_elm_type();	//要素タイプ 0:一定 1:線形
	int B_flag=UNLOOP;							//B_flagがLOOPなら境界は繋がっている。UNLOOPならつながってない

	if(CON->get_region_shape()==0) make_cube_region(CON,NODE,s_node_num);//解析領域は直方体
	else cout<<"解析領域形状が未定"<<endl;

	int KTJ=*s_node_num;
	int KTE=12*KTJ;		//最大要素数　3次元式もとめよ
	int nelm=0;			//現在の要素数

	vector <point3D> NODE2;
	point3D NODE2_0;
	for(int i=0;i<KTJ+8+1;i++) NODE2.push_back(NODE2_0); //接点の座標(+8しているのはｽｰﾊﾟｰﾎﾞｯｸｽの頂点数)
	for(int i=0;i<KTJ;i++)
	{
		for(int D=0;D<3;D++) NODE2[i+1].r[D]=NODE[i].r[D];//座標をコピー
		NODE2[i+1].material=AIR;
	}
	
	
	//要素クラス作成
	vector <element3D> ELEM2;
	element3D ELEM2_0;
	for(int i=0;i<KTE;i++) ELEM2.push_back(ELEM2_0);

	/////////////節点座標の正規化
    double xmin=NODE2[1].r[A_X];
    double ymin=NODE2[1].r[A_Y];
    double zmin=NODE2[1].r[A_Z];
    double xmax=xmin;
    double ymax=ymin;
    double zmax=zmin;

    ///座標の最大、最小値を求める
    for(int i=2;i<=*s_node_num;i++)
    {
        if(NODE2[i].r[A_X]<xmin) xmin=NODE2[i].r[A_X];
		else if(NODE2[i].r[A_X]>xmax) xmax=NODE2[i].r[A_X];
	
		if(NODE2[i].r[A_Y]<ymin) ymin=NODE2[i].r[A_Y];
		else if(NODE2[i].r[A_Y]>ymax) ymax=NODE2[i].r[A_Y];
	
		if(NODE2[i].r[A_Z]<zmin) zmin=NODE2[i].r[A_Z];
		else if(NODE2[i].r[A_Z]>zmax) zmax=NODE2[i].r[A_Z];
    }
    ////

    double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X軸方向の寸法
    double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y軸方向の寸法
    double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z軸方向の寸法
    double rmax=rax;		///最大寸法
    if(ray>rmax) rmax=ray;
    if(raz>rmax) rmax=raz;      //ここはelseにしたらダメ

    ///座標変換
    double rrm=1.000000/rmax;///こういう書き方をすることで、数値誤差を減らせる・・？
    for(int i=1;i<=*s_node_num;i++)
    {   //   A/Bという計算をしたとき、Ａの値によって微妙に1/Bという倍率がちがってくるのではないかと考えて、下のような書き方にしている
        NODE2[i].r[A_X]=(NODE2[i].r[A_X]-xmin)*rrm;
		NODE2[i].r[A_Y]=(NODE2[i].r[A_Y]-ymin)*rrm;
		NODE2[i].r[A_Z]=(NODE2[i].r[A_Z]-zmin)*rrm;
    }
    rax*=rrm;
    ray*=rrm;
    raz*=rrm;
    /////

	int FINE_sw=OFF;
	delaun3D(CON,NODE2,ELEM2, KTJ, KTE, rax, ray, raz,s_node_num,&nelm, FINE_sw,rrm);
	/////ﾒｯｼｭ生成完了

	///メッシュ生成を確認
	double *val=new double[KTJ+1];
	for(int i=1;i<=*s_node_num;i++) val[i]=1;
	data_avs(*s_node_num,nelm,NODE2,ELEM2,KTJ,val,CON);
	delete [] val;

	//表面要素ELEMを生成
	BEMelement3D ELEM0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int jelm=ELEM2[i].elm[j];
			if(jelm==0)					//表面なら
			{
				ELEM.push_back(ELEM0);
				int ia=ELEM2[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
				int ib=ELEM2[i].node[4-(j-1)/2*2];
				int ic=ELEM2[i].node[3-(j/2%2)*2];
				
				ia-=1;		//デローニ分割のプログラム内では、節点番号は1から始まるが、BEMでは0から始まるので、ここで-1した節点番号が、該当する節点番号である。
				ib-=1;
				ic-=1;

				ELEM[elemnum].node[0]=ia;
				ELEM[elemnum].node[1]=ib;
				ELEM[elemnum].node[2]=ic;
				double iaic[3];//ia→icのﾍﾞｸﾄﾙ成分格納
				double iaib[3];//ia→ibのﾍﾞｸﾄﾙ成分格納
				for(int D=0;D<3;D++)
				{
					iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
					iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
				}
				///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのﾍﾞｸﾄﾙとなる
				double S[3];//上記のﾍﾞｸﾄﾙ成分格納
				S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
				S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
				S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
				
				double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
				ELEM[elemnum].S=SS;
				////面積Sがもとまった
				for(int D=0;D<3;D++) ELEM[elemnum].direct[D]=S[D]/SS;//外向き単位法線ﾍﾞｸﾄﾙ
				for(int D=0;D<3;D++) ELEM[elemnum].r[D]=(NODE[ia].r[D]+NODE[ib].r[D]+NODE[ic].r[D])/3;	//重心
				ELEM[elemnum].map=0;			//初期化
				elemnum++;
			}
		}
	}////要素情報生成完了

	//要素の境界条件および材質作成
	for(int i=0;i<elemnum;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];
		int BD=Diric;			//境界条件
		for(int j=0;j<3;j++) if(NODE[N[j]].boundary_condition==Neumn) BD=Neumn;//ひとつでもノイマン型の節点を含んでいれば、その面はノイマン型
		ELEM[i].boundary_condition=BD;
		ELEM[i].material=AIR;			//ここで作成される表面要素はすべてAIR   ???
	}

	*s_elem_num=elemnum;
}

//直方体解析領域作成関数
void make_cube_region(mpsconfig *CON,vector<BEMpoint3D> &NODE,int *s_node_num)
{
	int node_num=0;
	int divN[3];								//各辺の分割数
	double divL[3];								//分割幅
	double V1=0;//CON->get_V();
	double V2=CON->get_V();						//ポテンシャル値

	BEMpoint3D NODE01;

	double Xmin=CON->get_minX();				//解析領域
	double Xmax=CON->get_maxX();
	double Ymin=CON->get_minY();
	double Ymax=CON->get_maxY();
	double Zmin=CON->get_minZ();
	double Zmax=CON->get_maxZ();

	
	divN[A_X]=10;
	divN[A_Y]=10;
	divN[A_Z]=10;
	divL[A_X]=(Xmax-Xmin)/divN[A_X];
	divL[A_Y]=(Ymax-Ymin)/divN[A_Y];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];
	//底面
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin;					//解析領域の底面
			NODE[node_num].potential=V1;
			NODE[node_num].slop1=0;						//初期化
			NODE[node_num].boundary_condition=Diric;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
			node_num++;
		}
	}
	//上面
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmax;					//解析領域の上面
			NODE[node_num].potential=V2;
			NODE[node_num].slop1=0;						//初期化
			NODE[node_num].boundary_condition=Diric;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
			node_num++;
		}
	}

	
	//側面Y
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//初期化
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
			node_num++;
		}
	}
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymax;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//初期化
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
			node_num++;
		}
	}

	//側面
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//初期化
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
			node_num++;
		}
	}
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmax;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];					
			NODE[node_num].potential=0;
			NODE[node_num].slop1=0;						//初期化
			NODE[node_num].boundary_condition=Neumn;
			NODE[node_num].C=0.5;
			NODE[node_num].particle=-1;					//対応する粒子が存在しないから-1を格納
			node_num++;
		}
	}
	*s_node_num=node_num;
}

//一定要素用の、動的節点・要素情報生成関数
void set_BEM3D_dynaic_model_for_CONSTANT(mpsconfig *CON,vector<mpsparticle> &PART,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,int *d_node_num,int *d_elem_num,int particle_number,int fluid_number)
{
	int node_num=0;
	//int elem_num=0;
	double le=CON->get_distancebp();
	BEMpoint3D NODE01;
	
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
			for(int D=0;D<3;D++) NODE[node_num].r[D]=PART[i].r[D];
		//	NODE[node_num].boundary_condition=Diric;
			NODE[node_num].boundary_condition=BOTH;		//２媒質のときはこちら
			NODE[node_num].potential=0;					//0?
			NODE[node_num].slop1=0;					//初期化
			NODE[node_num].slop2=0; 				//初期化
			NODE[node_num].C=0.5;					//一定要素は0.5
			NODE[node_num].particle=i;				//対応する粒子番号記憶
			BEMID[i]=node_num;							//BEMに出力するしるし
			node_num++;
		}
	}///

	///一時的に、すこし内側にも節点を作成しておく
	for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON && PART[i].toBEM==ON)
		{
			NODE.push_back(NODE01);
			for(int D=0;D<3;D++) NODE[node_num].r[D]=PART[i].r[D]+direct[D][i]*le;//少し内側の座標を入力
			NODE[node_num].boundary_condition=NOCALC;		//あとで消すという意味で、NOCALC
			NODE[node_num].potential=0;					//0?
			NODE[node_num].slop1=0;					//初期化
			NODE[node_num].slop2=0; 				//初期化
			NODE[node_num].C=0.5;					//一定要素は0.5
			NODE[node_num].particle=-1;				//対応する粒子は存在しない
			node_num++;
		}
	}/////節点情報作成完了*/

	//delaun3D用の要素クラス作成
	int KTJ=node_num;
	int KTE=12*KTJ;		//最大要素数　3次元式もとめよ
	int nelm=0;			//現在の要素数
	vector <point3D> NODE2;
	point3D NODE2_0;
	for(int i=0;i<KTJ+8+1;i++) NODE2.push_back(NODE2_0); //接点の座標(+8しているのはｽｰﾊﾟｰﾎﾞｯｸｽの頂点数)
	for(int i=0;i<KTJ;i++)
	{
		for(int D=0;D<3;D++) NODE2[i+1].r[D]=NODE[i].r[D];//座標をコピー
		NODE2[i+1].material=FLUID;
	}

	//要素クラス作成
	vector <element3D> ELEM2;
	element3D ELEM2_0;
	for(int i=0;i<KTE;i++) ELEM2.push_back(ELEM2_0);

	/////////////節点座標の正規化
    double xmin=NODE2[1].r[A_X];
    double ymin=NODE2[1].r[A_Y];
    double zmin=NODE2[1].r[A_Z];
    double xmax=xmin;
    double ymax=ymin;
    double zmax=zmin;

    ///座標の最大、最小値を求める
    for(int i=2;i<=node_num;i++)
    {
        if(NODE2[i].r[A_X]<xmin) xmin=NODE2[i].r[A_X];
		else if(NODE2[i].r[A_X]>xmax) xmax=NODE2[i].r[A_X];
	
		if(NODE2[i].r[A_Y]<ymin) ymin=NODE2[i].r[A_Y];
		else if(NODE2[i].r[A_Y]>ymax) ymax=NODE2[i].r[A_Y];
	
		if(NODE2[i].r[A_Z]<zmin) zmin=NODE2[i].r[A_Z];
		else if(NODE2[i].r[A_Z]>zmax) zmax=NODE2[i].r[A_Z];
    }
    ////

    double rax=sqrt((xmax-xmin)*(xmax-xmin));	///X軸方向の寸法
    double ray=sqrt((ymax-ymin)*(ymax-ymin));	///Y軸方向の寸法
    double raz=sqrt((zmax-zmin)*(zmax-zmin));	///Z軸方向の寸法
    double rmax=rax;		///最大寸法
    if(ray>rmax) rmax=ray;
    if(raz>rmax) rmax=raz;      //ここはelseにしたらダメ

    ///座標変換
    double rrm=1.000000/rmax;///こういう書き方をすることで、数値誤差を減らせる・・？
    for(int i=1;i<=node_num;i++)
    {   //   A/Bという計算をしたとき、Ａの値によって微妙に1/Bという倍率がちがってくるのではないかと考えて、下のような書き方にしている
        NODE2[i].r[A_X]=(NODE2[i].r[A_X]-xmin)*rrm;
		NODE2[i].r[A_Y]=(NODE2[i].r[A_Y]-ymin)*rrm;
		NODE2[i].r[A_Z]=(NODE2[i].r[A_Z]-zmin)*rrm;
    }
    rax*=rrm;
    ray*=rrm;
    raz*=rrm;
    /////

	int FINE_sw=OFF;
	delaun3D(CON,NODE2,ELEM2, KTJ, KTE, rax, ray, raz,&node_num,&nelm, FINE_sw,rrm);
	/////ﾒｯｼｭ生成完了

	for(int i=1;i<=nelm;i++) ELEM2[i].material=FLUID;			//ここで作成される要素はすべてFLUID

	//ここで、表面が凹なところに不要なメッシュがきられている場合は、NOCALCを含まない要素のなかで、辺の長さが著しく長い要素の材質を空気にするという処理を追加すること

	///メッシュ生成を確認
	double *val=new double[KTJ+1];
	for(int i=1;i<=node_num;i++) val[i]=1;
	data_avs(node_num,nelm,NODE2,ELEM2,KTJ,val,CON);
	delete [] val;

	//表面要素ELEMを生成
	BEMelement3D ELEM0;
	int elemnum=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=4;j++)
		{
			int jelm=ELEM2[i].elm[j];
			int flag=OFF;
			if(jelm==0) flag=ON;
			//else if(ELEM2[jelm].material==AIR) flag=ON;
			if(flag==ON)					//表面なら
			{
				ELEM.push_back(ELEM0);
				int ia=ELEM2[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
				int ib=ELEM2[i].node[4-(j-1)/2*2];
				int ic=ELEM2[i].node[3-(j/2%2)*2];
				
				ia-=1;		//デローニ分割のプログラム内では、節点番号は1から始まるが、BEMでは0から始まるので、ここで-1した節点番号が、該当する節点番号である。
				ib-=1;
				ic-=1;

				ELEM[elemnum].node[0]=ia;
				ELEM[elemnum].node[1]=ib;
				ELEM[elemnum].node[2]=ic;
				double iaic[3];//ia→icのﾍﾞｸﾄﾙ成分格納
				double iaib[3];//ia→ibのﾍﾞｸﾄﾙ成分格納
				for(int D=0;D<3;D++)
				{
					iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
					iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
				}
				///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのﾍﾞｸﾄﾙとなる
				double S[3];//上記のﾍﾞｸﾄﾙ成分格納
				S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
				S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
				S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
				
				double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
				ELEM[elemnum].S=SS;
				////面積Sがもとまった
				for(int D=0;D<3;D++) ELEM[elemnum].direct[D]=S[D]/SS;//外向き単位法線ﾍﾞｸﾄﾙ
				for(int D=0;D<3;D++) ELEM[elemnum].r[D]=(NODE[ia].r[D]+NODE[ib].r[D]+NODE[ic].r[D])/3;	//重心
				ELEM[elemnum].map=0;			//初期化
				elemnum++;
				//cout<<elemnum<<endl;
			}
		}
	}////要素情報生成完了

	//NOCALC節点の消去
	int erasenum=node_num/2;//半分の節点がNOCALCだから、この数だけ消去を実行
	for(int i=0;i<erasenum;i++) NODE.pop_back();
	node_num-=erasenum;

	//要素の境界条件および材質作成
	for(int i=0;i<elemnum;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];
		int BD=NODE[N[0]].boundary_condition;			//いまのところ、境界条件はN[0]もN[1]も同じと仮定して、N[0]のを使用
		ELEM[i].boundary_condition=BD;
		ELEM[i].material=FLUID;			//ここで作成される表面要素はすべてFLUID  
	}

	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	cout<<"動的節点数="<<node_num<<" 動的要素数="<<elemnum<<endl;
	*d_node_num=node_num;
	*d_elem_num=elemnum;

	delete [] BEMID;
	
}


void couple3D_NODE_and_ELEM(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &s_NODE,vector<BEMelement3D> &s_ELEM,vector<BEMpoint3D> &dy_NODE,vector<BEMelement3D> &dy_ELEM,vector<REGION> &region)
{
	//s_NODEとdy_NODEを合体させてNODEに格納する
	int s_node_num=int (s_NODE.size());
	int dy_node_num=int (dy_NODE.size());
	int s_elem_num=int (s_ELEM.size());
	int dy_elem_num=int (dy_ELEM.size());
	int countN=0;
	int countE=0;

	BEMpoint3D NODE01;
	BEMelement3D ELEM01;
	REGION temp;
	region.push_back(temp);
	region[0].start=0;					//最初の計算領域は必ず計算点０から始まる

	int *s_Nid=new int[s_node_num];		//移動先の節点番号格納。s_NODE[i]だった節点はNODE[s_Nid[i]]に格納される
	int *dy_Nid=new int[dy_node_num];


	//節点情報生成
	for(int i=0;i<s_node_num;i++)
	{
		NODE.push_back(NODE01);
		for(int D=0;D<3;D++) NODE[countN].r[D]=s_NODE[i].r[D];
		NODE[countN].boundary_condition=s_NODE[i].boundary_condition;
		NODE[countN].potential=s_NODE[i].potential;
		NODE[countN].slop1=s_NODE[i].slop1;
		NODE[countN].slop2=s_NODE[i].slop2;
		NODE[countN].C=s_NODE[i].C;
		NODE[countN].particle=s_NODE[i].particle;
		s_Nid[i]=countN;			//i番目のs_NODEはcountN番目の節点
		countN++;
	}
	//region.push_back(temp);//第２領域
	//region[1].start=countN;
	for(int i=0;i<dy_node_num;i++)
	{
		NODE.push_back(NODE01);
		for(int D=0;D<3;D++) NODE[countN].r[D]=dy_NODE[i].r[D];
		NODE[countN].boundary_condition=dy_NODE[i].boundary_condition;
		NODE[countN].potential=dy_NODE[i].potential;
		NODE[countN].slop1=dy_NODE[i].slop1;
		NODE[countN].slop2=dy_NODE[i].slop2;
		NODE[countN].C=dy_NODE[i].C;
		NODE[countN].particle=dy_NODE[i].particle;
		dy_Nid[i]=countN;			//i番目のdy_NODEはcountN番目の節点
		countN++;
	}
	//region[0].end=countN;
	//region[1].end=countN;
	

	//要素情報作成
	for(int i=0;i<s_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		for(int j=0;j<3;j++) ELEM[countE].node[j]=s_ELEM[i].node[j];
		for(int D=0;D<3;D++) ELEM[countE].r[D]=s_ELEM[i].r[D];
		ELEM[countE].boundary_condition=s_ELEM[i].boundary_condition;
		ELEM[countE].S=s_ELEM[i].S;
		for(int D=0;D<3;D++) ELEM[countE].direct[D]=s_ELEM[i].direct[D];
		ELEM[countE].map=s_ELEM[i].map;
		ELEM[countE].material=s_ELEM[i].material;
		countE++;
	}
	for(int i=0;i<countE;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];	//この段階で格納されている節点番号は、s_NODEのものであり、NODEのものではない
		for(int j=0;j<3;j++) N[j]=s_Nid[N[j]];		//節点番号をNODE基準のものに書き換え
		for(int j=0;j<3;j++) ELEM[i].node[j]=N[j];
	}
	region.push_back(temp);//第２領域
	region[1].start=countE;
	for(int i=0;i<dy_elem_num;i++)
	{
		ELEM.push_back(ELEM01);
		for(int j=0;j<3;j++) ELEM[countE].node[j]=dy_ELEM[i].node[j];
		for(int D=0;D<3;D++) ELEM[countE].r[D]=dy_ELEM[i].r[D];
		ELEM[countE].boundary_condition=dy_ELEM[i].boundary_condition;
		ELEM[countE].S=dy_ELEM[i].S;
		for(int D=0;D<3;D++) ELEM[countE].direct[D]=dy_ELEM[i].direct[D];
		ELEM[countE].map=dy_ELEM[i].map;
		ELEM[countE].material=dy_ELEM[i].material;
		countE++;
	}
	for(int i=s_elem_num;i<countE;i++)
	{
		int N[3];
		for(int j=0;j<3;j++) N[j]=ELEM[i].node[j];	//この段階で格納されている節点番号は、dy_NODEのものであり、NODEのものではない
		for(int j=0;j<3;j++) N[j]=dy_Nid[N[j]];		//節点番号をNODE基準のものに書き換え
		for(int j=0;j<3;j++) ELEM[i].node[j]=N[j];
	}
	region[0].end=countE;
	region[1].end=countE;
	cout<<"全節点数="<<countN<<" 全要素数="<<countE<<endl;

	//ファイル出力して座標を確認
	ofstream fp("BEM_input.dat");
	double times=0.0001;
	for(int i=0;i<countE;i++)
	{
		int n=ELEM[i].node[0];
		fp<<ELEM[i].r[A_X]<<"\t"<<ELEM[i].r[A_Y]<<"\t"<<ELEM[i].r[A_Z]<<endl;
	}//*/
	//for(int i=0;i<countN;i++) fp<<NODE[i].r[A_X]<<"\t"<<NODE[i].r[A_Y]<<"\t"<<NODE[i].potential<<endl;
	fp.close();


	delete [] s_Nid;	
	delete [] dy_Nid;
}