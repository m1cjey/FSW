#include "stdafx.h"	

//メモリリーク検出
//#define _CRTDBG_MAP_ALLOC

#include"header.h"		//主要なヘッダーファイルはまとめてこのなか。
//#include"define.h"		//#define 格納
//#include"PART.h"		//class PART定義
//#include"CONFIG.h"		//CONF.cppはインクルードしなくても、CONFIG.hをインクルードするだけでよい
#include"BEMclass.h"	//FEM3D関係のclass 定義
//#include"FEM3Dclass.h"	//FEM3D関係のclass 定義
//#include<omp.h>
//#include<vector>

#include"function.h"

//////
//メモリリーク検出
//#include <crtdbg.h>
//#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)
////////*/

int _tmain(int argc, _TCHAR* argv[])
{
	//メモリリーク検出
	 _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	 _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	 _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	 ///////////////

	mpsconfig CON;//解析条件格納ｸﾗｽ
	
    int particle_number;//全粒子数
    int fluid_number;	//流体粒子数
	int out;			//fluid_number<=i<outがINWALL群で、out<=iがOUTWALL群になる。 なってないなら自動であとで並び替える
	//int count;			//カウント用変数
	int count_avs=0;
	int order_sw=OFF;	//粒子を並び替える必要があるかどうかのフラグ
    double dt=CON.get_dt();
    int STEP=CON.get_step();
    double g[DIMENTION]={0,0,0};
    if(CON.get_dimention()==2)  g[A_Y]=CON.get_g();
    if(CON.get_dimention()==3)  g[A_Z]=CON.get_g();
    /////MPSにおける定数計算
	double n0=initial_pnd(CON.get_re(),CON.get_dimention(),CON.get_model_set_way()); //初期粒子密度n0
    double N0=initial_pnd(CON.get_re2(),CON.get_dimention(),CON.get_model_set_way());//ﾗﾌﾟﾗｼｱﾝ用初期粒子数密度
    double n0_4=initial_pnd(CON.get_re4(),CON.get_dimention(),CON.get_model_set_way());//freeon用初期粒子密度
    double lamda=calclamda(&CON); //ﾗﾌﾟﾗｼｱﾝ用λ
	CON.set_Cst(calc_Cst(&CON));
    double TIME=0;		//解析時間
	double Umax=0;		//最大速度
	double mindis;		//ｸｰﾗﾝ用の最低粒子間距離
	int t_old=0;		//EM_distance関係の出力用。直前にEMしたときのステップ数を記憶
    printf("n0= %10.8f\n",n0);
    printf("N0= %10.8f\n",N0);
    printf("lamda= %10.8f\n",lamda);
	cout<<"Cst= "<<CON.get_Cst()<<", Cst*"<<CON.get_C_times()<<"="<<CON.get_Cst()*CON.get_C_times()<<endl;
    //////////////////////////////////*/

	//mindis=CON.get_distancebp();

	/*
	////////
	cout<<"particlemovieの変換"<<endl;
	ifstream fin("particle_movie_0522_termo.mgf");
	if(!fin) cout<<"cannot open particle_movie_0522_termo"<<endl;
	fin.unsetf(ifstream::dec);
	fin.setf(ifstream::skipws);

	string b;
	for(int i=0;i<2;i++) getline(fin, b);//最初の6行分進める

	ofstream fout("particle_movie_mod.mgf");
	if(!fout) cout<<"cannot open particle_movie_mod"<<endl;
	
	fout<<"# Micro AVS Geom:2.00"<<endl;
	fout<<800<<endl;	

	

	/////
	///////////////
	
	


	for(int i=0;i<800;i++)
	{
		int num=0;
		double level=0;
		double x=0;
		double y=0;
		double z=0;
		double d=0;
		double red=0;
		double green=0;
		double blue=0;
		double red_m=0;
		double green_m=0;
		double blue_m=0;
		double T=0;

		double Tmax=303;
		double Tmin=293;

		//読み込み
		if(i==0) for(int j=0;j<4;j++) getline(fin, b);//最初の4行分進める
		else for(int j=0;j<5;j++) getline(fin, b);//最初の4行分進める
		fin>>num;

		//出力
		fout<<"step"<<i+1<<endl;
		fout<<"sphere"<<endl;
		fout<<"time="<<(i+1)*0.0002<<endl;
		fout<<"color"<<endl;
		fout<<num<<endl;
		
		for(int j=0;j<num;j++)
		{
			fin>>x;
			fin>>y;
			fin>>z;
			fin>>d;
			fin>>red;
			fin>>green;
			fin>>blue;

			//温度コンターの変換
			T=sqrt(red)*100+293;
			if(T>Tmax) T=Tmax;
			level=(T-Tmin)/(Tmax-Tmin);
			red_m=level*level;
			green_m=-4*(level*(level-1));
			blue_m=1-red_m;
 
			fout<<x<<" "<<y<<" "<<z<<" ";//座標出力
				
			fout<<d<<" ";//粒子の大きさ出力
	
			fout<<red_m<<" "<<green_m<<" "<<blue_m<<endl;//色出力
		}
		
	}
	fin.close();
	fout.close();
	*/
	///////
	/////

	/////初期粒子配置書き込み　restart=ONの場合は粒子数読み込み
    if(CON.get_restart()==OFF)
	{   
		if(CON.get_model_inherit()==0)
		{
			if(CON.get_model_set_way()==0)		set_initial_placement(&CON,&particle_number); //正方格子によるモデルセット。速いが表面形状は階段状
			else if(CON.get_model_set_way()==1)	set_initial_placement_using_MD(&CON,&particle_number);//分子動力学によるモデルセット。計算時間はかかるが、表面をきれいに表現
		}
		else cout<<"前回の初期粒子モデルの引継ぎ"<<endl;
		
	}
    else if(CON.get_restart()==ON)
    {
		ifstream fin("number.dat");
		if(!fin) cout<<"cannot open number.dat"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
    
		fin>>particle_number;
		fin>>TIME;
		fin.close();
		cout<<"前回における解析を引き継ぎ"<<endl;
    }
    //////////////////////
	if(CON.get_model_inherit()==1)
	{
		ifstream fa("particle_number.dat");
		if(!fa) cout<<"cannot open number.dat"<<endl;
		fa>>particle_number;
		fa.close();
	}
    cout<<"粒子数は"<<particle_number<<endl;

	/////壁境界ポリゴン化のテスト
	int wall_poly_number=0;//壁境界のポリゴン数
	wall_poly_number=3;//読み込みなどでもとめること
	double wall_a[3];//2次元の直線方程式ax+by+c=0　3DFEMを参考にnodeやedgeのような書き方をすること
	double wall_b[3];
	double wall_c[3];
	if(CON.get_wall_poly()==1)
	{
		//壁境界の設定 //メッシュなどから読み込むことを考えると、この位置はまずいかも
		if(CON.get_model_number()==1)
		{
			if(CON.get_dimention()==2)
			{
				wall_a[0]=1; wall_b[0]=0; wall_c[0]=-0.011;//左側内壁
				wall_a[1]=1; wall_b[1]=0; wall_c[1]=0.011;//右側内壁
				wall_a[2]=0; wall_b[2]=1; wall_c[2]=-0.021;//下側内壁
			}
		}
	}
	////壁重み関数の計算



	

	/////INDEX関係
    //int *INDEX=new int[CON.get_number_of_mesh()];	//各格子に含まれる粒子数を格納
    cout<<"X_mesh="<<CON.get_X_mesh()<<" Y_mesh="<<CON.get_Y_mesh()<<" Z_mesh="<<CON.get_Z_mesh()<<endl;
    cout<<"number_of_mesh="<<CON.get_number_of_mesh()<<endl;
    ///////////////////*/

	clock_t t1=clock();		//経過時間を秒で表現するには、CLOCKS_PER_SECで割る必要がある

	/////各種のPART[i](particle[i])を作成
	mpsparticle PART0;
	vector<mpsparticle> PART;
	for(int i=0;i<particle_number;i++) PART.push_back(PART0);
    //////////////////////////////

	////BEM用のclass作成
	vector<point2D> NODE;
	vector<element2D> ELEM;
	vector<BEMpoint3D> BEMNODE3D;		//固定節点情報class
	vector<BEMelement3D> BEMELEM3D;	//固定要素情報class
	int s_node_num;					//固定節点数,動的節点数
	int s_elem_num;	

	//////////////////////////

	////FEM用のclass作成
	vector<point3D> NODE_FEM3D;
	vector<element3D> ELEM_FEM3D;
	vector<edge3D> EDGE_FEM3D;
	vector<point3D> NODE_jw;
	vector<element3D> ELEM_jw;
	int node_FEM3D=0;
	int nelm_FEM3D=0;	//全節点数,要素数
	int nedge_FEM3D=0;

	//////////////////////////


	///粒子データをファイルから読み取り
	input_particle_data(&CON,PART,particle_number,1);//最後の引数に1を渡したらinitialデータを読み取り、それ以外なら前ｽﾃｯﾌﾟデータを読み取る

	//各粒子数をカウント or 並び替え
	calc_numbers_of_particles_and_change_the_order(&CON,PART,particle_number,&fluid_number,&out,&order_sw);

	//double *Un[DIMENTION];
	//for(int D=0;D<DIMENTION;D++) Un[D]=new double [fluid_number];			//nｽﾃｯﾌﾟ時の速度(陽解析直前)を記憶. 
	double *previous_Un[DIMENTION];
	for(int D=0;D<DIMENTION;D++) previous_Un[D]=new double [fluid_number];	//n-1ｽﾃｯﾌﾟ時の速度(陽解析直前)を記憶. 仮の速度を求めるのに使用する
	double *PND2=new double [particle_number];		//圧力計算用に、各ｽﾃｯﾌﾟ開始時の粒子数密度を記憶

	if(CON.get_set_zero_speed()==ON && CON.get_restart()==ON)
	{
		for(int i=0;i<fluid_number;i++) for(int D=0;D<CON.get_dimention();D++) PART[i].u[D]=0;//restartでも速度は初期化。開発用
	}

	double distance=0.0;//EM_distance>0のときに用いる判定用距離
	for(int i=0;i<particle_number;i++)
	{
		PART[i].T=CON.get_initialT();//初期化
		PART[i].heat_gene_before1=0;
		PART[i].heat_generation=0;
		PART[i].L=CON.get_distancebp();
		PART[i].dir_Pst=0.0;
		PART[i].dir_Pem=0.0;

	}

	//最大速度を求めておく。普通はゼロだが、restartしたときとかにここで求めておかないと、1stepめのcuran数算出がおかしくなる。
	for(int i=0;i<particle_number;i++)
	{ 
		double speed=0;//粒子速度
		for(int D=0;D<DIMENTION;D++) speed+=PART[i].u[D]*PART[i].u[D];
		if(speed>Umax)  Umax=speed;		
	}
	Umax=sqrt(Umax);
	
	unsigned int time0=GetTickCount();	//解析始めの時間を記憶
	//計算スタート
	for(int t=1;t<=STEP;t++)
	{	
		unsigned int timet=GetTickCount();
			
		cout<<"陽解析 start:step="<<t<<" 流体粒子数="<<fluid_number<<" 全粒子数="<<particle_number<<endl;

		//各粒子数をカウント or 並び替え
		if(order_sw==ON) calc_numbers_of_particles_and_change_the_order(&CON,PART,particle_number,&fluid_number,&out,&order_sw);//内部でorder_swをOFFに戻している
		
		////////////////////陽解析/////////////////////////////

		//近隣粒子関係計算
		calc_neighbor_relation(&CON,PART,particle_number,n0_4,fluid_number,out);
		
		/*/////陽解析の前にreloadINDEX	
		reload_INDEX(&CON,PART,particle_number,INDEX);//格子内の粒子数更新
		int **MESH = new int *[CON.get_number_of_mesh()];
		count=0;
		for(int i=0;i<CON.get_number_of_mesh();i++)
		{       
			count+=INDEX[i];
			MESH[i]=new int [INDEX[i]];
		}
		if(count>particle_number) cout<<"INDEX error 粒子数増加?"<<endl;
		if(count<particle_number) cout<<"INDEX error 粒子数減少?"<<endl;
		reload_INDEX2(&CON,PART,particle_number,MESH);
		////////////////////////*/

		////////壁重み関数の計算
		//if(CON.get_wall_poly()==1 && CON.get_model_number()==0) calc_wallZ(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,fluid_number,out);

		
		unsigned int timeA=GetTickCount();
		/*////
		if(CON.get_freeon()==1) freeon(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,fluid_number,out);//表面判定
		else if(CON.get_freeon()==2) freeon2(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,fluid_number,out,t);//表面判定
		//else if(CON.get_freeon()==4) freeon4(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,n0,fluid_number,e0,out,t);
		else cout<<"表面判定未解決"<<endl;
		if(CON.get_surface_judge2()==ON)
		{
			surface_judge2(&CON,PART,fluid_number,particle_number);
		}
		cout<<"陽解析前の粒子依存関係計算終了　－－time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		////*/

		////最低粒子間距離をもとめる
		//double min0=CON.get_distancebp();//最低粒子間距離
		mindis=CON.get_distancebp();
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
				if(dis<mindis)
				{
					type1=PART[i].type;
					surface1=PART[i].surface;
					type2=PART[j].type;
					surface2=PART[j].surface;
					mindis=dis;
					X1=PART[i].r[A_X];
					Y1=PART[i].r[A_Y];
					Z1=PART[i].r[A_Z];
				}
			}
			/////最低粒子間距離がもとまった
			//cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
			//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
		}
		/////最低粒子間距離がもとまった
		cout<<"mindis="<<mindis<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;

		//ｸｰﾗﾝ数によるdt決定
		culan(&CON,PART,fluid_number,t,&dt,mindis,Umax,g);

		//初期粒子数密度分布を出力する
		if(t==1) output_particle_density(&CON, PART, fluid_number, n0, particle_number, t);

		//if(CON.get_modify_position()==ON) modify_position(&CON,PART, fluid_number,dt);

		///粘性項配列
		double *laplacian[DIMENTION];
		for(int D=0;D<DIMENTION;D++) laplacian[D]=new double [fluid_number];
		///表面張力関係配列
		double *potential[DIMENTION];
		for(int D=0;D<DIMENTION;D++) potential[D]= new double [fluid_number];
		//粒子の外力
		double *F[DIMENTION];					
		for(int D=0;D<DIMENTION;D++) F[D]=new double [fluid_number];
		for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=0;//初期化

		double *Un[DIMENTION];
		for(int D=0;D<DIMENTION;D++) Un[D]=new double [fluid_number];			//nｽﾃｯﾌﾟ時の速度(陽解析直前)を記憶.
		for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) Un[D][i]=PART[i].u[D];//陽解析前の速度を記憶

		//////粘性項計算
		calc_viscous_term(&CON,PART,fluid_number,particle_number,dt,N0,laplacian,lamda,t);
		
		//表面張力計算
		calc_surface_tension(&CON,PART,fluid_number,dt,particle_number,n0,potential,t);
		
		//FEM
		//cout<<"density="<<CON.get_density()<<endl;
		if(CON.get_EM_method()!=0)
		{
			if(CON.get_EM_distance()==OFF)
			{
				if(t==1 || (t-1)%CON.get_EM_interval()==0)
				{
					unsigned int timeF=GetTickCount();
					if(CON.get_EM_method()==1) FEM3D_calculation(&CON,&node_FEM3D,&nelm_FEM3D,&nedge_FEM3D,NODE_FEM3D,ELEM_FEM3D,EDGE_FEM3D,fluid_number,F, t, TIME,PART, fluid_number, particle_number,dt,NODE_jw,ELEM_jw,CON);
					else if(CON.get_EM_method()==2)
					{
						if(CON.get_dimention()==2) BEM2D(&CON,NODE,ELEM,&s_node_num,&s_elem_num,PART,fluid_number,particle_number,F,t);
						if(CON.get_dimention()==3) BEM3D(&CON,BEMNODE3D,BEMELEM3D,&s_node_num,&s_elem_num,PART,fluid_number,particle_number,F,t);
					}
					else if(CON.get_EM_method()==3) Magnetic_Moment_Method(&CON,PART,F, n0, lamda, fluid_number, particle_number);

					cout<<"FEM終了　－－time="<<(GetTickCount()-timeF)*0.001<<"[sec]"<<endl;

					for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) PART[i].F[D]=F[D][i];//求めた電磁力をPARTに格納

					//plot_F(&CON,PART,fluid_number,F,t);//電磁力出力
				
					//FEMに入らないstepでも直前のFEMで求めた電磁力を用いるため、ファイル出力
					//ポスト処理で消した粒子に対応できないので却下
					//ofstream g("F_FEM.dat");
					//for(int i=0;i<fluid_number;i++) g<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
					//g.close();
				}
				else//FEMをしないときはファイルから読み込み//ポスト処理で消した粒子に対応できないので却下。PART.Fに格納した値を使う
				{
					if(CON.get_T_field()==ON) for(int i=0;i<particle_number;i++) PART[i].heat_generation=PART[i].heat_gene_before1;//渦電流損を計算しないときでも値を覚えておく
				}
			}
			else//粒子の移動距離によってFEMを行うかどうか決める
			{
				if(t==1 || distance>CON.get_EM_distance()*CON.get_distancebp())
				{	
					distance=0.0;//リセット
					if(t!=1) cout<<"distance>"<<CON.get_EM_distance()<<"leのためFEM開始"<<endl;
					unsigned int timeF=GetTickCount();
					if(CON.get_EM_method()==1) FEM3D_calculation(&CON,&node_FEM3D,&nelm_FEM3D,&nedge_FEM3D,NODE_FEM3D,ELEM_FEM3D,EDGE_FEM3D,fluid_number,F, t, TIME,PART, fluid_number, particle_number,dt,NODE_jw,ELEM_jw,CON);
					else if(CON.get_EM_method()==2)
					{
						if(CON.get_dimention()==2) BEM2D(&CON,NODE,ELEM,&s_node_num,&s_elem_num,PART,fluid_number,particle_number,F,t);
						if(CON.get_dimention()==3) BEM3D(&CON,BEMNODE3D,BEMELEM3D,&s_node_num,&s_elem_num,PART,fluid_number,particle_number,F,t);
					}
					else if(CON.get_EM_method()==3) Magnetic_Moment_Method(&CON,PART,F, n0, lamda, fluid_number, particle_number);

					cout<<"FEM終了　－－time="<<(GetTickCount()-timeF)*0.001<<"[sec]"<<endl;

					for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) PART[i].F[D]=F[D][i];//求めた電磁力をPARTに格納

					/////////////////////////////CPU time と解析物理時間の関係性をプロット
					if(t==1)
					{
						ofstream t1("EMstep.dat");		//横軸ステップ数、縦軸dtのグラフ
						t1.close();
					}
				
					ofstream t1("EMstep.dat",ios :: app);
					t1<<t<<" "<<t-t_old<<endl;
					t1.close();
					t_old=t;

		////////////////////

					//plot_F(&CON,PART,fluid_number,F,t);//電磁力出力
				
					//FEMに入らないstepでも直前のFEMで求めた電磁力を用いるため、ファイル出力
					//ofstream g("F_FEM.dat");
					//for(int i=0;i<fluid_number;i++) g<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
					//g.close();
				}
				else//FEMをしないときはファイルから読み込み
				{
					if(CON.get_T_field()==ON) for(int i=0;i<particle_number;i++) PART[i].heat_generation=PART[i].heat_gene_before1;//渦電流損を計算しないときでも値を覚えておく
				}
			}
		}
		

		//電磁力出力
		if(CON.get_EM_method()!=0)		plot_F(&CON,PART,fluid_number,F,t);


		/////仮の速度および位置決定
		renewal_u_and_r_in_positive(&CON,PART,fluid_number,t,dt,&Umax,potential,laplacian,g,previous_Un,F);


		for(int D=0;D<DIMENTION;D++) delete [] potential[D];
		for(int D=0;D<DIMENTION;D++) delete [] laplacian[D];
		for(int D=0;D<DIMENTION;D++) delete [] F[D];

		cout<<"陽解析終了 umax="<<sqrt(Umax)<<"  limit U="<<0.2*mindis/dt<<endl;
		

		///////////////////粒子が動いたのでINDEX更新

		if(CON.get_temporary_r()==ON)//仮の位置を計算しないのなら、ここでfreeonはしなくてよい
		{
			calc_neighbor_relation(&CON,PART,particle_number,n0_4,fluid_number,out);//近隣粒子関係計算
		}
		
		/*///
		for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH[i];
		delete [] MESH;

		reload_INDEX(&CON,PART,particle_number,INDEX);//格子の粒子数更新
	
		int **MESH2 = new int *[CON.get_number_of_mesh()];
		for(int i=0;i<CON.get_number_of_mesh();i++)  MESH2[i]=new int [INDEX[i]];
	
		reload_INDEX2(&CON,PART,particle_number,MESH2);
		if(CON.get_temporary_r()==ON)//仮の位置を計算しないのなら、ここでfreeonはしなくてよい
		{
			unsigned int timeB=GetTickCount();
			if(CON.get_freeon3sw()==OFF)
			{
				if(CON.get_freeon()==1) freeon(&CON,PART,particle_number,n0_4,INDEX,MESH2,&mindis,fluid_number,out);//表面判定
				else if(CON.get_freeon()==2) freeon2(&CON,PART,particle_number,n0_4,INDEX,MESH2,&mindis,fluid_number,out,t);//表面判定
				//else if(CON.get_freeon()==4) freeon4(&CON,PART,particle_number,n0_4,INDEX,MESH2,&mindis,n0,fluid_number,e0,out,t);
			}
			if(CON.get_freeon3sw()==ON) freeon3(&CON,PART,particle_number,out);//粒子数密度のみ再計算
			cout<<"陰解析前の粒子依存関係計算終了 time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;
		}
		////*/

		//温度場
		if(CON.get_T_field()==ON)
		{
			if(CON.get_model_number()==26) 
			{
				set_Q_boundary(&CON, PART, fluid_number, particle_number,dt,t);//境界の温度そのものが実験などで既にもとまっている場合、それを対応する境界に位置する粒子に与える
			}
			if(CON.get_T_solution()==0) calc_Temperature(&CON,PART,fluid_number,particle_number,N0,lamda,dt,t);
			if(CON.get_T_solution()==1) calc_temperature_implicity(&CON,PART,fluid_number,particle_number,N0,lamda,dt,t);
		}

		
		///陰解析(内部で速度・位置修正あり)
		if(CON.get_iteration_count()==1)
		{
			if(CON.get_P_twice()==OFF) negative1(&CON,PART,fluid_number,particle_number,out,t, dt, lamda, N0,PND2,n0,Un);
			//if(CON.get_P_twice()==ON) negative1_twice(&CON,PART, fluid_number, particle_number, out, t, dt, lamda, N0, PND2, n0,Un, n0_4);
			if(CON.get_P_twice()>0) negative1_twice(&CON,PART, fluid_number, particle_number, out, t, dt, lamda, N0, PND2, n0,Un, n0_4);
		}
		//if(CON.get_iteration_count()>1)  negative2(&CON,PART,fluid_number,particle_number,t, dt, lamda, n0,INDEX,MESH2,N0,out,Un);
		else if(CON.get_iteration_count()==-1) negative3(&CON,PART, fluid_number, particle_number, t, dt, lamda, N0,NODE_FEM3D,ELEM_FEM3D, out);
		////*/
		
		//移動粒子を移動させる場合
		if(CON.get_move_prtcl()==ON) move_particle(&CON,PART,particle_number,fluid_number,dt);

		if(CON.get_modify_position()!=OFF) modify_position(&CON,PART, fluid_number,dt,particle_number);
		
		///MESH消去
		//for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH2[i];
		//delete [] MESH2;

		//陰解析後の最大速度測定////////////////
		double umax2=0;//最大速度
		int umax_ID=0;	//最大速度を有する粒子番号
		for(int i=0;i<particle_number;i++)
		{ 
			double speed=0;//粒子速度
			for(int D=0;D<DIMENTION;D++) speed+=PART[i].u[D]*PART[i].u[D];
			if(speed>umax2) 
			{
				umax2=speed;
				umax_ID=i;
			}	
		}
		//最大速度を有する粒子の周辺粒子の速度が、最大速度に比べて1m/s以上小さい場合は、最大速度を修正する
		double umax3=0;	//最大速度粒子の周辺粒子のなかでの最大速度
		for(int kk=0;kk<PART[umax_ID].N;kk++)
		{
			int j=PART[umax_ID].NEI[kk];
			double speed=0;//粒子速度
			for(int D=0;D<DIMENTION;D++) speed+=PART[j].u[D]*PART[j].u[D];
			if(speed>umax3)  umax3=speed;
		}
		umax3=sqrt(umax3);
		umax2=sqrt(umax2);
		if(umax2-umax3>1 && umax2>0)
		{	
			
			cout<<"最大速度の修正あり umax2="<<umax2<<"->"<<umax3<<endl;
			for(int D=0;D<DIMENTION;D++) PART[umax_ID].u[D]=PART[umax_ID].u[D]*umax3/umax2;
			umax2=umax3;
		}

		cout<<"陰解析終了 umax="<<umax2<<endl;
		if(umax2>Umax) Umax=umax2;

		distance+=Umax*dt;	//最大移動距離を更新
		if(CON.get_EM_method()>0 && CON.get_EM_distance()>0) cout<<"distance="<<distance<<" "<<CON.get_EM_distance()<<"le="<<CON.get_EM_distance()*CON.get_distancebp()<<endl;
			////////////////////////////////////////

		TIME+=dt;///時間更新
		
		double Pmax=0;
		for(int i=0;i<fluid_number;i++)
		{ 
			if(PART[i].P>=Pmax) Pmax=PART[i].P;
		}
		cout<<"Pmax="<<Pmax<<endl;

		if(t==1)
	{
		ofstream t1("Pmax.dat");		//横軸ステップ数、縦軸dtのグラフ
		t1.close();
	}


	ofstream t1("Pmax.dat",ios :: app);
	t1<<t<<" "<<Pmax<<endl;
	t1.close();
		
		cout<<"解析物理時間="<<TIME<<" time="<<(GetTickCount()-timet)*0.001<<"[sec]"<<endl;

		ofstream foutc("step.dat");
		foutc<<t<<endl;
		foutc.close();

		///ポスト処理：各物理量出力＆ｸｰﾗﾝ数によるdt改変&microAVS出力
		post_processing(&CON,PART,fluid_number,particle_number,dt,Umax,mindis,t,TIME,time0,count_avs);

		if(t==1 || t%CON.get_interval()==0) output_alldata_AVS(&CON,PART,fluid_number,particle_number,dt,Umax,mindis,t,TIME,time0,count_avs);
		
		//条件を満たした粒子を削除する
		//if(CON.get_check_region()==ON)	delete_particle(CON,PART,&particle_number,&fluid_number,n0_4,t);

		//restart用ﾌｧｲﾙ生成
		post_processing3(&CON,PART,fluid_number,particle_number,t,TIME);

		//領域外の粒子を検査
		order_sw=check_position(&CON,PART, fluid_number,&particle_number);//領域外粒子を検知すれば、order_sw=ONになる

		//条件を満たした粒子を削除する
		//if(CON.get_check_region()==ON)	delete_particle(CON,PART,&particle_number,&fluid_number,n0_4,t);

		//何か特別なことを測定・調査する関数.何かは自分でプログラムを作成して決める
		if(CON.get_check_something()==ON) check_something(&CON,PART, fluid_number, n0, particle_number, t);

		
		cout<<endl;

		for(int D=0;D<DIMENTION;D++) delete [] Un[D];
		//_CrtDumpMemoryLeaks();//メモリリークの検出
		
	}

	
	//delete [] INDEX;
	delete [] PND2;
	//for(int D=0;D<DIMENTION;D++) delete [] Un[D];
	for(int D=0;D<DIMENTION;D++) delete [] previous_Un[D];

	clock_t t2=clock();
	cout<<"CPU time="<<(t2-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
	MessageBeep(MB_ICONEXCLAMATION);//作業の終了を知らせるBEEP音

	//new int;
	_CrtDumpMemoryLeaks();//メモリリークの検出

	return 0;
}

//重み関数
double kernel(double r,double dis)
{
    //return (1-dis/r)*(1-dis/r);
	return r/dis-1;
	//return r*r/(dis*dis)-1;
	//return r*r*r/(dis*dis*dis)-1;
	//return (1-dis/r)*(1-dis/r);
	//return 1;
	//return log(r/dis);
}

//重み関数
double kernel2(double r,double dis,double d)
{
	//return (1-dis/r)*(1-dis/r);
	return r/dis-1;
    //return r*r*r/(dis*dis*dis)-1;
	//return r*r*r*r/(dis*dis*dis*dis)-1;
	//return pow(r,d)/pow(dis,d);
	//return log(r/dis);
}

double kernel_in_WLSM(double dis, double R)
{
	//R:影響半径
	double r=dis/R;
	double val=1-6*r*r+8*r*r*r-3*r*r*r*r;
	return val;

}


//初期粒子密度の計算
double initial_pnd(double r,int dimention,int calc_type)
{
	int size = (int)(r+1);//計算領域
	double dis;//距離
	double pnd=0;
	int count=0;
	if(dimention==2)
	{
		if(calc_type==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis<=r )
					{
						pnd+=kernel(r,dis);
						count++;
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
					if(dis!=0 && dis<=r )
					{
						pnd+=kernel(r,dis);
						count++;
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
						if(dis!=0 && dis<=r )
						{
							pnd+=kernel(r,dis);
							count++;
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
						if(dis!=0 && dis<=r )
						{
							pnd+=kernel(r,dis);
							count++;
						}
					}			
				}
			}
		}
	}
	cout<<"n0のcount="<<count<<endl;
	return pnd;  
}

///ﾗﾌﾟﾗｼｱﾝ用変数λ計算関数
double calclamda(mpsconfig *CON)
{
	int dimention=CON->get_dimention();	//解析次元
	int Ini_place=CON->get_model_set_way();	//初期粒子配置方法　0=正方 1=細密
	double R=CON->get_re2();			//ラプラシアン用影響半径
	int size = (int)(R+1);//計算領域
	int count=0;
	double dis;//距離      
	double w;
	double pnd=0;
	double lam=0;
	if(dimention==2)
	{
		if(Ini_place==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis<=R)
					{
						double length=dis*CON->get_distancebp();
						w=kernel(R,dis);
						pnd+=w;
						lam+=length*length*w;
						count++;
					}
				}				
			}
		}
		else if(Ini_place==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
					dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis<=R )
					{
						double length=dis*CON->get_distancebp();
						w=kernel(R,dis);
						pnd+=w;
						lam+=length*length*w;
						count++;
					}			
				}
			}
		}
	}
	else if(dimention==3)
	{
		if(Ini_place==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis=sqrt((double)(i*i+j*j+k*k));
						if(dis!=0 && dis<=R)
						{
							double length=dis*CON->get_distancebp();
							w=kernel(R,dis);
							pnd+=w;
							lam+=length*length*w;
							count++;
						}
					}
				}				
			}
		}
		else if(Ini_place==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-size;i<=size;i++)
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
						if(dis!=0 && dis<=R )
						{
							double length=dis*CON->get_distancebp();
							w=kernel(R,dis);
							pnd+=w;
							lam+=length*length*w;
							count++;
						}
					}			
				}
			}
		}
	}
	lam/=pnd;
	cout<<"λのcount="<<count<<endl;
	return lam;  
}

///粒子データ読み取り関数
void input_particle_data(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int t)
{
	if(t==1)//最初はinitial_input.datから読み込み
	{
		ifstream fin("initial_input.dat");
		if(!fin) cout<<"cannot open mps_input.dat"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
		for(int i=0;i<particle_number;i++)
		{       
			fin>>PART[i].id;
			for(int D=0;D<DIMENTION;D++) fin>>PART[i].r[D];
			for(int D=0;D<DIMENTION;D++) fin>>PART[i].u[D];
			fin>>PART[i].P;
			fin>>PART[i].h;
			fin>>PART[i].val;
			fin>>PART[i].type;
			fin>>PART[i].materialID;
			fin>>PART[i].surface;
			fin>>PART[i].toBEM;
			//old_Aの初期値を与える
			for(int D=0;D<DIMENTION;D++) PART[i].old_A[D]=0;

		}
		fin.close();	
	}
	///////////////////////*/
	
	if(t!=1) //2STEP以降はmps_input.datから読み込み
	{
		ifstream fin("mps_input.dat");
		if(!fin) cout<<"cannot open mps_input.dat"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
		for(int i=0;i<particle_number;i++)
		{
			fin>>PART[i].id;
			for(int D=0;D<DIMENTION;D++) fin>>PART[i].r[D];
			for(int D=0;D<DIMENTION;D++) fin>>PART[i].u[D];
			fin>>PART[i].P;
			fin>>PART[i].h;
			fin>>PART[i].val;
			fin>>PART[i].type;
			fin>>PART[i].materialID;
			fin>>PART[i].surface;
			fin>>PART[i].toBEM;
		}
		fin.close();
	}
}

//粒子数カウント関数 & 並び替え
void calc_numbers_of_particles_and_change_the_order(mpsconfig *CON,vector <mpsparticle> &PART,int particle_number,int *fluid_number,int *out,int *order_sw)
{
	//各粒子数をカウント
	int fluid_num=0;						//流体粒子数
	int inwall_num=0;
	int out_num=0;
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID) fluid_num++;
		else if(PART[i].type==INWALL) inwall_num++;
	}
	
	out_num=fluid_num+inwall_num;	//fluid_number<=i<outがINWALL群で、out<=iがOUTWALL群になる。
	*out=out_num;
	*fluid_number=fluid_num;

	//並び替え
	//mpsparticle PART_temp;			//並び替え用粒子クラス
	for(int i=0;i<fluid_num;i++)
	{
		if(PART[i].type!=FLUID)
		{
			cout<<"並び替え必要あり"<<endl;
		}
	}
	*order_sw=OFF;
}

//INDEX更新関数
void reload_INDEX(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int *INDEX)
{       
	///説明:各格子のなかに含まれる粒子数を数える  格子番号は２次元では左下で0。X方向につれ＋１で、右上で最大(X_mesh*Y_mesh) ３次元ではＺ方向にも増えていく
	
	double width=CON->get_distancebp()*CON->get_dx();		//格子幅
	for(int i=0;i<CON->get_number_of_mesh();i++) INDEX[i]=0;
	for(int i=0;i<particle_number;i++)
	{
		int X=(int)((PART[i].r[A_X]-CON->get_minX())/width);//X方向に何個目の格子か 
		int Y=(int)((PART[i].r[A_Y]-CON->get_minY())/width);//Y方向に何個目の格子か
		int Z=(int)((PART[i].r[A_Z]-CON->get_minZ())/width);//Z方向に何個目の格子か
		int number=Z*CON->get_X_mesh()*CON->get_Y_mesh()+Y*CON->get_X_mesh()+X;//粒子iを含む格子の番号
        PART[i].index=number;
		INDEX[number]++;	
	}
}

//INDEX更新関数その２ MESHに粒子番号を格納する
void reload_INDEX2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int **MESH)
{
	int *count=new int [CON->get_number_of_mesh()];
	for(int i=0;i<CON->get_number_of_mesh();i++) count[i]=0;//初期化

	for(int i=0;i<particle_number;i++)
	{
		int number=PART[i].index;	
        MESH[number][count[number]]=i;
		count[number]++;
	}
	delete [] count;
}

//法線ベクトル作成関数
void direct_f(mpsconfig *CON,vector<mpsparticle> &PART,int i,double *direct[DIMENTION])
{
	double R=CON->get_re3()*CON->get_distancebp();//法線ﾍﾞｸﾄﾙ計算に利用する影響半径

	double px=PART[i].r[A_X]+CON->get_distancebp();//x+le
	double mx=PART[i].r[A_X]-CON->get_distancebp();//x-le
	double py=PART[i].r[A_Y]+CON->get_distancebp();//y+le
	double my=PART[i].r[A_Y]-CON->get_distancebp();//y-le
	double pz=PART[i].r[A_Z]+CON->get_distancebp();//z+le
	double mz=PART[i].r[A_Z]-CON->get_distancebp();//z-le
	
	double pnd_px=pnd_for_direct(CON,PART,px,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_mx=pnd_for_direct(CON,PART,mx,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_py=pnd_for_direct(CON,PART,PART[i].r[A_X],py,PART[i].r[A_Z],R,i);
	double pnd_my=pnd_for_direct(CON,PART,PART[i].r[A_X],my,PART[i].r[A_Z],R,i);
	double pnd_pz=pnd_for_direct(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],pz,R,i);
	double pnd_mz=pnd_for_direct(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],mz,R,i);

	direct[A_X][i]=(pnd_px-pnd_mx)/(2*CON->get_distancebp());
	direct[A_Y][i]=(pnd_py-pnd_my)/(2*CON->get_distancebp());
	direct[A_Z][i]=(pnd_pz-pnd_mz)/(2*CON->get_distancebp());
	
	double a=sqrt(direct[A_X][i]*direct[A_X][i]+direct[A_Y][i]*direct[A_Y][i]+direct[A_Z][i]*direct[A_Z][i]);
	if(a!=0)
	{ 
		direct[A_X][i]/=a;
		direct[A_Y][i]/=a;
		direct[A_Z][i]/=a;
	}
}

//法線ベクトル作成関数
void direct_f2(mpsconfig *CON,vector<mpsparticle> &PART,int i,double *direct[DIMENTION])
{
	double R=CON->get_re3()*CON->get_distancebp();//法線ﾍﾞｸﾄﾙ計算に利用する影響半径

	double px=PART[i].r[A_X]+CON->get_distancebp();//x+le
	double mx=PART[i].r[A_X]-CON->get_distancebp();//x-le
	double py=PART[i].r[A_Y]+CON->get_distancebp();//y+le
	double my=PART[i].r[A_Y]-CON->get_distancebp();//y-le
	double pz=PART[i].r[A_Z]+CON->get_distancebp();//z+le
	double mz=PART[i].r[A_Z]-CON->get_distancebp();//z-le
	
	double pnd_px=pnd_for_direct2(CON,PART,px,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_mx=pnd_for_direct2(CON,PART,mx,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_py=pnd_for_direct2(CON,PART,PART[i].r[A_X],py,PART[i].r[A_Z],R,i);
	double pnd_my=pnd_for_direct2(CON,PART,PART[i].r[A_X],my,PART[i].r[A_Z],R,i);
	double pnd_pz=pnd_for_direct2(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],pz,R,i);
	double pnd_mz=pnd_for_direct2(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],mz,R,i);

	direct[A_X][i]=(pnd_px-pnd_mx)/(2*CON->get_distancebp());
	direct[A_Y][i]=(pnd_py-pnd_my)/(2*CON->get_distancebp());
	direct[A_Z][i]=(pnd_pz-pnd_mz)/(2*CON->get_distancebp());
	
	double a=sqrt(direct[A_X][i]*direct[A_X][i]+direct[A_Y][i]*direct[A_Y][i]+direct[A_Z][i]*direct[A_Z][i]);
	if(a!=0)
	{ 
		direct[A_X][i]/=a;
		direct[A_Y][i]/=a;
		direct[A_Z][i]/=a;
	}
}

//////direct用粒子数密度測定
double pnd_for_direct(mpsconfig *CON,vector<mpsparticle> &PART,double x,double y,double z,double R,int i)
{
	///粒子iの位置とはずれるので、MESHを使用したほうがよい
	///教科書では個数をかぞえているが、ここでは重み関数を用いる。
	//R=CON->get_re3()*CON->get_distancebp();
	double spnd=0;
	

	for(int k=0;k<PART[i].N3;k++)
	{       
		int j=PART[i].NEI3[k];
		double X=PART[j].r[A_X]-x;
		double Y=PART[j].r[A_Y]-y;
		double Z=PART[j].r[A_Z]-z;
		double dis=sqrt(X*X+Y*Y+Z*Z);
		//if(dis<R) spnd++;   //教科書どうり
		if(dis<R)
		{
			double w=(1-dis/R)*(1-dis/R);
			spnd+=w;
		}
	}
	return spnd;
}

//////direct用粒子数密度測定
double pnd_for_direct2(mpsconfig *CON,vector<mpsparticle> &PART,double x,double y,double z,double R,int i)
{
	///粒子iの位置とはずれるので、MESHを使用したほうがよい
	///教科書では個数をかぞえているが、ここでは重み関数を用いる。
	//R=CON->get_re3()*CON->get_distancebp();
	double spnd=0;
	

	for(int k=0;k<PART[i].N3;k++)
	{       
		int j=PART[i].NEI3[k];
		if(PART[j].type==FLUID)
		{
			double X=PART[j].r[A_X]-x;
			double Y=PART[j].r[A_Y]-y;
			double Z=PART[j].r[A_Z]-z;
			double dis=sqrt(X*X+Y*Y+Z*Z);
			//if(dis<R) spnd++;   //教科書どうり
			if(dis<R)
			{
				double w=(1-dis/R)*(1-dis/R);
				spnd+=w;
			}
		}
	}
	return spnd;
}

///ポスト処理関数
void post_processing(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double Umax,double mindis,int t,double TIME, unsigned int time0,int count_avs)
{
	double le=CON->get_distancebp();
	unsigned int timeA=GetTickCount();	//計算開始時間

	cout<<"各物理量ﾌﾟﾛｯﾄ開始: ";
	//cout<<"各物理量ﾌﾟﾛｯﾄ開始"<<endl;
	
	//_CrtDumpMemoryLeaks();
	
	//AVSに粒子データ出力/////////////////////////////////
	
	//cout<<"粒子avsデータ出力開始----";
	if(CON->get_curan()==0) if(t==1 || t%CON->get_interval()==0) particle_movie_AVS(CON,PART,fluid_number,particle_number,t,TIME);
	else if(CON->get_curan()>0)
	{
		if(t==1 || TIME>CON->get_interval()*dt*count_avs) particle_movie_AVS(CON,PART,fluid_number,particle_number,t,TIME);
		
		count_avs++;
		ofstream f("culan_AVSstep.dat");
		f<<count_avs<<endl;
		f.close();
	}
	cout<<"ok"<<endl;
	
		
	///速度をプロット
//	plot_speed(CON ,PART,particle_number,fluid_number);
	plot_speed_each(CON ,PART,particle_number,fluid_number,t);

	//////座標ﾌﾟﾛｯﾄ/////////////////////////
	//////座標ﾌﾟﾛｯﾄ/////////////////////////
	ofstream gnu1("0.dat");//解析終了後の全粒子座標をプロット
	ofstream gnu2("suf.dat");//表面粒子だけをプロット
	ofstream gnu4("suf2.dat");//流体表面粒子だけ
	if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{ 
			gnu1<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			if(PART[i].surface==ON)
			{
				gnu2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
				if(PART[i].type==FLUID) gnu4<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			}
		}
	}
	else if(CON->get_dimention()==3)
	{
		
		for(int i=0;i<particle_number;i++)
		{ 
			if(PART[i].surface==ON && PART[i].type==FLUID)//流体だけ表示したいとき
			{
				gnu2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
				if(fabs(PART[i].r[A_Y])<0.5*le) gnu4<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			}
		}
		
	}
	gnu1.close();
	gnu2.close();
	gnu4.close();

	//////////////////////////////////////////////////*/

	/////////////////////////////CPU time と解析物理時間の関係性をプロット
	if(t==1)
	{
		ofstream t1("dt.dat");		//横軸ステップ数、縦軸dtのグラフ
		t1.close();

		ofstream t2("CPU_time1.dat");//横軸CPU time、縦軸TIMEのグラフ
		t2.close();

		ofstream t3("CPU_time2.dat");//横軸CPU time、縦軸stepのグラフ
		t3.close();

	}


	ofstream t1("dt.dat",ios :: app);
	t1<<t<<" "<<dt<<endl;
	t1.close();

	ofstream t2("CPU_time1.dat",ios :: app);
	t2<<(GetTickCount()-time0)*0.001<<" "<<TIME<<endl;
	t2.close();

	ofstream t3("CPU_time2.dat",ios :: app);
	t3<<(GetTickCount()-time0)*0.001<<" "<<t<<endl;
	t3.close();
	////////////////////
	
	//////////////////////////////平均粒子密度&圧力を表示
	double ave_n0=0;
	double ave_P=0;
	int count=0;
	for(int i=0;i<fluid_number;i++) 
	{
	    if(PART[i].surface==OFF)
	    {
	        ave_n0+=PART[i].PND;
			ave_P+=PART[i].P;
			count++;
	    }
	}
	if(count!=0){ ave_n0/=count;ave_P/=count;}
	cout<<"average n0="<<ave_n0<<" average P="<<ave_P<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	//////////////////////////////////////////////////*/
}

///microAVS用粒子動画出力関数
void particle_movie_AVS(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double T)
{
	//t:タイムｽﾃｯﾌﾟ　T:総合時間
	double le=CON->get_distancebp();
	double times=1;

	if(t==1)
	{
		ofstream fout("particle_movie.mgf");

		fout<<"# Micro AVS Geom:2.00"<<endl;
		fout<<CON->get_step()/CON->get_interval()+1<<endl;//microAVSに出力する総ｽﾃｯﾌﾟ数。ファイル出力はCON->get_interval()回に1回と最初に行う。
		
		fout.close();
	}


	ofstream avs("particle_movie.mgf",ios :: app);
	avs<<"step"<<t/CON->get_interval()+1<<endl;
	avs<<"sphere"<<endl;
	avs<<"time="<<T<<endl;
	avs<<"color"<<endl;

	if(t==1)
	{
		ofstream fout2("particle_movie_section.mgf");

		fout2<<"# Micro AVS Geom:2.00"<<endl;
		fout2<<CON->get_step()/CON->get_interval()+1<<endl;//microAVSに出力する総ｽﾃｯﾌﾟ数。ファイル出力はCON->get_interval()回に1回と最初に行う。
		
		fout2.close();
	}


	ofstream avs2("particle_movie_section.mgf",ios :: app);
	avs2<<"step"<<t/CON->get_interval()+1<<endl;
	avs2<<"sphere"<<endl;
	avs2<<"time="<<T<<endl;
	avs2<<"color"<<endl;

	if(t==1)
	{
		ofstream fout3("particle_movie_trace.mgf");

		fout3<<"# Micro AVS Geom:2.00"<<endl;
		fout3<<CON->get_step()/CON->get_interval()+1<<endl;//microAVSに出力する総ｽﾃｯﾌﾟ数。ファイル出力はCON->get_interval()回に1回と最初に行う。
		
		fout3.close();
	}


	ofstream avs3("particle_movie_trace.mgf",ios :: app);
	avs3<<"step"<<t/CON->get_interval()+1<<endl;
	avs3<<"sphere"<<endl;
	avs3<<"time="<<T<<endl;
	avs3<<"color"<<endl;

	double red,green,blue;	//粒子の色を表現する3原色

/////////流体粒子のうち、壁粒子付近にあるものは表面でなくても表示させる//////////

	///freeon2関数の説明:flag1[i]の導入により高速化。ただし、1CPUなら早いが、ﾏﾙﾁCPUによる並列化は厳しい
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	double d=2;
	if(CON->get_dimention()==3) d=3;

	int out=particle_number;

	if(CON->get_T_field()==ON && CON->get_insulate()==1)
	{
		out=particle_number;	//非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
	}

	int *flag1=new int [out];		//検査フラグ。0:未検査　1:検査済み
	int *flag2=new int [out];		//流体粒子のうち、壁粒子に近いものはON。出力用
	///初期化
	//#pragmallel for
	for(int i=0;i<out;i++)
	{
		flag1[i]=0;
		flag2[i]=0;
	}
	
	if(CON->get_plot_nearbywall()==ON)
	{
		//INDEX,MESH作成
		int *INDEX=new int[CON->get_number_of_mesh()];	//各格子に含まれる粒子数を格納
		reload_INDEX(CON,PART,particle_number,INDEX);//格子内の粒子数更新
	
		int **MESH = new int *[CON->get_number_of_mesh()];
		int count=0;
		for(int i=0;i<CON->get_number_of_mesh();i++)
		{       
			count+=INDEX[i];
			MESH[i]=new int [INDEX[i]];
		}
		reload_INDEX2(CON,PART,particle_number,MESH);

		//cout<<"INDEX,MESHreload完了"<<endl;
	
		for(int i=0;i<out;i++)
		{  
			////粒子数測定
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
								if(flag1[j]==0 && j!=i)//まだ検査してないなら
								{
									double X=PART[j].r[A_X]-PART[i].r[A_X];
									double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
									double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
									double dis=sqrt(X*X+Y*Y+Z*Z);
							
									if(dis<=CON->get_re()*le && PART[j].type==INWALL)
									{       
										flag2[i]=ON;
									}
								}
							}
						}
					}
				}
			}
			flag1[i]=1;//検査終了
		}

		//cout<<"表示粒子数計測完了"<<endl;

		for(int i=0;i<CON->get_number_of_mesh();i++) delete [] MESH[i];
		delete [] MESH;
		delete [] INDEX;
	}
	

////////////////////////////////////////////////////

	
	if(CON->get_AVS()==0)	//粒子の動きだけを表示。
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		int num2=0;//表示する粒子数
		if(CON->get_dimention()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON->get_dimention()==3) 
		{
			for(int i=0;i<particle_number;i++) 
			{
				//if(PART[i].surface==ON)
				{
					if(CON->get_AVS_HALF()==ON)
					{
						if(PART[i].type==INWALL || PART[i].type==OUTWALL) if(PART[i].r[A_X]>=0) num++;//3次元の場合、内部流体は表示しない
						if(PART[i].type==FLUID)
						{
							if(PART[i].surface==ON)
							{
								num++;
							}
							if(PART[i].surface==OFF && flag2[i]==ON)
							{
								num++;
							}
							//if(PART[i].r[A_Y]<0.5*le && PART[i].r[A_Y]>-0.5*le)
							if(PART[i].r[A_X]<0.5*le && PART[i].r[A_X]>-0.5*le)
							//if(PART[i].r[A_X]<0.5*le && PART[i].r[A_X]>-0.5*le)
							{
								num2++;
							}
						}
					}
					if(CON->get_AVS_HALF()==OFF)
					{
						if(PART[i].type==INWALL || PART[i].type==OUTWALL) num++;
						if(PART[i].type==FLUID && PART[i].surface==ON) num++;

						if(PART[i].r[A_X]<0.5*le && PART[i].r[A_X]>-0.5*le)
						//if(PART[i].r[A_X]<0.5*le && PART[i].r[A_X]>-0.5*le)
						{
							num2++;
						}
					}
				}
			}
		}
		
		avs<<num<<endl;
		avs2<<num2<<endl;
		/////////////////////////////////

		////粒子出力
		if(CON->get_dimention()==2)
		{    
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
				else if(PART[i].type==INWALL && PART[i].surface==OFF) {red=1;green=0;blue=0;}
				//else if(PART[i].type==INWALL) {red=1;green=0;blue=0;}
				//else if(PART[i].type==OUTWALL) {red=0.5;green=0.5;blue=0;}
	        	else {red=0.5;green=0.5;blue=0;}//壁粒子
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
				avs<<CON->get_distancebp()/2<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}        
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				//if(PART[i].surface==ON) 
				{
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		//else if(PART[i].type==INWALL) {red=1;green=0;blue=0;}
					//else if(PART[i].type==OUTWALL) {red=0.5;green=0.5;blue=0;}
					else {red=0.5;green=0.5;blue=0;}//壁粒子
					//if(PART[i].type==FLUID && PART[i].T<CON->get_MP()) {red=1;green=0;blue=0;}

					if(CON->get_model_number()==19) //FSWの場合
					{
						if(PART[i].type==FLUID &&PART[i].r[A_X]>0) {red=1;green=0;blue=0;}
						else if(PART[i].type==FLUID &&PART[i].r[A_X]<=0) {red=0;green=0;blue=1;}

					}

			   
					if(CON->get_AVS_HALF()==ON)
					{
						
						if(PART[i].type==INWALL || PART[i].type==OUTWALL)
						{
							if(PART[i].r[A_X]>=0)
							{
							avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
							avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
							avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
							}
						}
						
						if(PART[i].type==FLUID)
						{
							if(PART[i].surface==ON)
							{
								avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
								avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
								avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
							}
							if(PART[i].surface==OFF && flag2[i]==ON)
							{
								avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
								avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
								avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
							}
							//if(PART[i].r[A_Y]<0.5*le && PART[i].r[A_Y]>-0.5*le)
							if(PART[i].r[A_X]<0.5*le && PART[i].r[A_X]>-0.5*le)
							{
								avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
								avs2<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
								avs2<<red<<" "<<green<<" "<<blue<<endl;//色出力
								
							}
						}
					}
					if(CON->get_AVS_HALF()==OFF)
					{
						if(PART[i].type==INWALL || PART[i].type==OUTWALL)
						{
							avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
							avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
							avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
						}
						if(PART[i].type==FLUID && PART[i].surface==ON)
						{
							avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
							avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
							avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
						}
					}
				}        
			}
		}
	}

	/*				
	if(CON->get_AVS()==0)	//粒子の動きだけを表示。
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		if(CON->get_dimention()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].surface==ON || PART[i].type!=FLUID) num++;//3次元の場合、内部流体は表示しない
			}
		}
		
		avs<<num<<endl;
		/////////////////////////////////

		////粒子出力
		if(CON->get_dimention()==2)
		{    
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//壁粒子
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
				avs<<CON->get_distancebp()/2<<" ";//粒子の大きさ出力

				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}        
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].surface==ON || PART[i].type!=FLUID) 
				{
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//壁粒子
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力

					avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}        
			}
		}
	}
	*/
	else if(CON->get_AVS()==1)	//粒子の圧力をコンター表示
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		if(CON->get_dimention()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON->get_dimention()==3) num=fluid_number;//3次元では流体粒子のみを表示
		
		avs<<num<<endl;
		/////////////////////////////////


		double maxP=0;//圧力の最大値と最小値をもとめる
		double minP=0;
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].P>maxP) maxP=PART[i].P;
			if(PART[i].P<minP) minP=PART[i].P;
		}
		/////maxP,minPが求まった。

		double width=maxP-minP;
		if(width<0.0001) width=0.0001;
			
		if(CON->get_dimention()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
		
	}
	else if(CON->get_AVS()==2)	//粒子の温度をコンター表示
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		if(CON->get_dimention()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON->get_dimention()==3) num=fluid_number;//3次元では流体粒子のみを表示
		
		avs<<num<<endl;
		avs3<<num<<endl;
		/////////////////////////////////

		double le=CON->get_distancebp();
		double mass=CON->get_particle_mass();
		double T;//粒子iの温度
		double width=CON->get_maxT()-CON->get_minT();//温度の幅

		if(CON->get_dimention()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				T=PART[i].T;
				/*
				if(PART[i].type==FLUID)
				{
					double hs0=mass*CON->get_Cp()*CON->get_MP();//融解開始点のエンタルピー
					double hs1=hs0+CON->get_latent_H()*mass;//融解終了点のエンタルピー
					if(PART[i].h<hs0) T=PART[i].h/mass/CON->get_Cp();
					else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON->get_MP();
					else if(hs1<PART[i].h) T=CON->get_MP()+(PART[i].h-hs1)/mass/CON->get_Cp();
				}
				else 
				{
					/////
					if(CON->get_model_number()==26)//釜の場合、壁の場所に応じて物性地をかえる
					{
						double air_mass=CON->get_particle_volume()*1.205;//空気の重さ
						double air_Cp=1006;//空気の比熱
						if(CON->get_dimention()==2)
						{
							if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08)//上壁
							{
								T=PART[i].h/air_mass/air_Cp;//固体
							}
							else T=PART[i].h/CON->get_particle_volume()*CON->get_wall_density()/CON->get_wall_Cp();//固体
						}
						else if(CON->get_dimention()==3)
						{
						}
					}
					else T=PART[i].h/(CON->get_particle_volume()*CON->get_wall_density())/CON->get_wall_Cp();//固体
					
				}
				//////*/
    ///////////////
				double level=(T-CON->get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<PART[i].L*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力

				//トレーサーを混ぜた動画用の色設定
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				//else if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08) {red=0;green=1;blue=0;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=1;green=0;blue=0;}
				else if(PART[i].type==INWALL && PART[i].surface==OFF) {red=0.5;green=0.5;blue=0;}
	        	else if(PART[i].type==INWALL && PART[i].surface==ON) {red=1;green=0;blue=0;}
				
				else {red=0.5;green=0.5;blue=0;}//壁粒子

				//if(PART[i].type==FLUID && i%50==0) {red=1;green=0;blue=0;}//トレーサー
				if(t==1)
				{
					for(int i=0;i<particle_number;i++)
					{
						PART[i].trace=OFF;
						if(PART[i].type==FLUID && i%100==0) PART[i].trace=ON;
					}
				}
				if(PART[i].trace==ON) {red=0;green=1;blue=0;}//トレーサー

				avs3<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs3<<PART[i].L*times<<" ";//粒子の大きさ出力
	
				avs3<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				double hs0=mass*CON->get_Cp()*CON->get_MP();//融解開始点のエンタルピー
				double hs1=hs0+CON->get_latent_H()*mass;//融解終了点のエンタルピー
				if(PART[i].h<hs0) T=PART[i].h/mass/CON->get_Cp();
				else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON->get_MP();
				else if(hs1<PART[i].h) T=CON->get_MP()+(PART[i].h-hs1)/mass/CON->get_Cp();
				double level=(T-CON->get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<PART[i].L*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
	}
	else if(CON->get_AVS()==3)//壁粒子は非表示
	{
		//////////////////////////////////表示する粒子数を計算
		int num=fluid_number;//表示する粒子数
		avs<<num<<endl;
		/////////////////////////////////

		////粒子出力
		if(CON->get_dimention()==2) cout<<"error in AVS() 2Dは非対応"<<endl;      
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//壁粒子
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
				avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力

				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力        
			}
		}
	}
	else if(CON->get_AVS()==4)//流体表面粒子のみ
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		if(CON->get_dimention()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON->get_dimention()==3) for(int i=0;i<particle_number;i++) if(PART[i].type==FLUID && PART[i].surface==ON) num++;//3次元の場合、内部流体は表示しない	
		avs<<num<<endl;
		/////////////////////////////////

		////粒子出力
		if(CON->get_dimention()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//壁粒子
			   
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
				avs<<CON->get_distancebp()/2<<" ";//粒子の大きさ出力

				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
			
		}
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].surface==ON) 
				{
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//壁粒子
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs<<CON->get_distancebp()/2<<" ";//粒子の大きさ出力

					avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}        
			}
		}
	}
	if(CON->get_AVS()==5)	//壁だけを表示.ただの確認用
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		if(CON->get_dimention()==2) for(int i=0;i<particle_number;i++) if(PART[i].type==INWALL || PART[i].type==OUTWALL) num++;//2次元では全粒子を表示
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL)
				{
					//if(PART[i].toBEM==MOVE)
					{
						num++;//3次元の場合、流体は表示しない
					}
				}
			}
		}
		
		avs<<num<<endl;
		/////////////////////////////////

		////粒子出力
		if(CON->get_dimention()==2)
		{    
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL)
				{
					if(PART[i].type==INWALL) red=0;green=1;blue=0;//壁粒子
					if(PART[i].type==OUTWALL) red=0.5;green=0.5;blue=0;//壁粒子
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs<<CON->get_distancebp()/2<<" ";//粒子の大きさ出力

					avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}
			}
		}        
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL)
				{
					//if(PART[i].toBEM==MOVE)
					{
	        			if(PART[i].type==INWALL) red=0;green=1;blue=0;//壁粒子
						if(PART[i].type==OUTWALL) red=0.5;green=0.5;blue=0;//壁粒子
			    
						avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
						avs<<CON->get_distancebp()<<" ";//粒子の大きさ出力

						avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
					}
				}        
			}
		}
	}
	else if(CON->get_AVS()==6 && CON->get_dimention()==3)	//特定の粒子の動きだけを表示。粒子の選択はプログラムを直接変更するしかない
	{
		int num=0;//表示する粒子数
		int num2=0;//表示する粒子数
		int *output=new int [particle_number];//ONなら出力　OFFなら出力しない
		int *output2=new int [particle_number];//ONなら出力　OFFなら出力しない
		for(int i=0;i<particle_number;i++) 
		{
				output[i]=OFF;//初期化
				output2[i]=OFF;//初期化
		}

		vector<int> ID;//AVSに出力する粒子のidをファイルからこの配列に入力する
		
		/*///
		if(t==1 && CON->get_restart()==OFF)//最初のステップ時にファイルを生成
		{
			ofstream fout2("id_for_AVS.dat");//AVSで追跡する粒子のidを出力
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].r[A_X]>-0.005 && PART[i].r[A_X]<-0.002)
				{
					if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le)
					{
						if(PART[i].r[A_Z]>0.0015 && PART[i].r[A_Z]<0.0055)
						{
							fout2<<i<<endl;
						}
					}
				}
			}
			fout2.close();
		}
		
		
		//特定粒子のidを読み込み(1step時には書き込んですぐ読み込む形になる)
		ifstream fin("id_for_AVS.dat");
		if(!fin) cout<<"cannot open mps_input.dat"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
		int value=0;
		while(!fin.eof())
		{       
			fin>>value;
			ID.push_back(value);
		}
		fin.close();
	
		//output[i]入力
		for(int k=0;k<ID.size()-1;k++)//なんか上のやりかただと、size-1しないとうまくいかない。
		{
			output[ID[k]]=ON;
			num++;
		}///////////////////////////////////*/

		for(int i=0;i<fluid_number;i++)
		{
			//if(PART[i].r[A_Y]<0 && PART[i].r[A_Y]>-0.0075)//XZ平面図
			//if(PART[i].r[A_X]<-le )//XY平面図
			//if(PART[i].r[A_X]>0)
			//if(PART[i].type==INWALL)
			{
				output[i]=ON;
				num++;
			}
			if(PART[i].r[A_Y]<0.006+0.5*le && PART[i].r[A_Y]>0.006-0.5*le)
			{
				output2[i]=ON;
				num2++;
			}
		}

		//output[i]入力
		for(int i=fluid_number;i<particle_number;i++)
		{
			
			if(CON->get_tool_angle()==0)//ツール回転なしの場合
			{
				//if(PART[i].toBEM==MOVE && PART[i].r[A_Z]<=0.006-0.2*le)//プローブのみ表示
				//if(PART[i].toBEM==MOVE && abs(PART[i].r[A_Z])<=0.003-0.2*le)//プローブのみ表示
				//if(PART[i].toBEM==MOVE && PART[i].r[A_Y]<0) 
				if(PART[i].toBEM==MOVE)//ツールのみ表示 
				//if(PART[i].toBEM==MOVE && PART[i].r[A_X]>0)
				//if(PART[i].type==FLUID)//非表示 
				{
					output[i]=ON;
					num++;
				}
			}

			if(CON->get_tool_angle()>0)//ツールをx軸まわりに回転させる。進行方向が+y
			{
				double theta=PI/180*CON->get_tool_angle();//回転する角度
				
				double z=sin(-theta)*PART[i].r[A_Y]+cos(-theta)*PART[i].r[A_Z];//回転させる前のツール位置を求める
				
				if(t==1)//初期ステップの配置で判断し、あとは最初の判定に従って表示
				{
					if(PART[i].toBEM==MOVE && z<=0.006-0.2*le)//プローブのみ表示
					//if(PART[i].toBEM==MOVE && PART[i].r[A_Y]<0) 
					//if(PART[i].toBEM==MOVE)//ツールのみ表示 
					//if(PART[i].toBEM==MOVE && PART[i].r[A_X]>0)
					{
						PART[i].color=1;
					}
					else PART[i].color=2;
				}
				if(PART[i].color==1) 
				{
					output[i]=ON;
					num++;
				}
			}
			
			//断面図用output設定
			if(PART[i].toBEM==MOVE)//ツールのみ表示 
			{
				if(PART[i].r[A_Y]<0.006+0.5*le && PART[i].r[A_Y]>0.006-0.5*le)
				{
					output2[i]=ON;
					num2++;
				}
			}
		}///
		
		avs<<num<<endl;
		avs2<<num2<<endl;
		/////////////////////////////////

		////粒子出力
		
		for(int i=0;i<particle_number;i++)
		{
			if(CON->get_model_number()==19) //FSWの場合、最初の配置で決めた粒子の配色をずっとつかいたい
			{
				//if(PART[i].type!=FLUID) PART[i].color=3; 
				if(output[i]==ON) 
				{	
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		//else if(PART[i].type==INWALL) {red=0;green=1;blue=0;}//壁粒子
					else {red=0.5;green=0.5;blue=0;}//壁粒子

					if(PART[i].type==FLUID)
					{
						//if(PART[i].materialID==1)  {red=0;green=0;blue=1;}
						//else if(PART[i].materialID==2)  {red=1;green=0;blue=0;}
						if(t==1)
						{
							if(PART[i].r[A_X]>0)  
							{
								red=1;green=0;blue=0;
								PART[i].color=1;//1なら赤
							}
							else 
							{
								red=0;green=0;blue=1;
								PART[i].color=2;//2なら青
							}
						}
						else
						{
							
							if(PART[i].color==1) {red=1;green=0;blue=0;}
							else if(PART[i].color==2) {red=0;green=0;blue=1;}

						}
					}
					//else PART[i].color=3;
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs<<le<<" ";//粒子の大きさ出力

					avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}    

				if(output2[i]==ON) 
				{
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//壁粒子

					if(PART[i].type==FLUID)
					{
						//if(PART[i].materialID==1)  {red=0;green=0;blue=1;}
						//else if(PART[i].materialID==2)  {red=1;green=0;blue=0;}
						//if(PART[i].r[A_X]>0)  {red=1;green=0;blue=0;}
						//else {red=0;green=0;blue=1;}
						if(t==1)
						{
							if(PART[i].r[A_X]>0)  
							{
								red=1;green=0;blue=0;
								PART[i].color=1;//1なら赤
							}
							else 
							{
								red=0;green=0;blue=1;
								PART[i].color=2;//2なら青
							}
						}
						else
						{
							
							if(PART[i].color==1) {red=1;green=0;blue=0;}
							else if(PART[i].color==2) {red=0;green=0;blue=1;}

						}
					}
			    
					avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs2<<le<<" ";//粒子の大きさ出力

					avs2<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}  
			}
			else 
			{
				if(output[i]==ON) 
				{	
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//壁粒子

					if(PART[i].type==FLUID)
					{
						//if(PART[i].materialID==1)  {red=0;green=0;blue=1;}
						//else if(PART[i].materialID==2)  {red=1;green=0;blue=0;}
						if(PART[i].r[A_X]>0)  {red=1;green=0;blue=0;}
						else {red=0;green=0;blue=1;}
					}
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs<<le<<" ";//粒子の大きさ出力

					avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}    

				if(output2[i]==ON) 
				{
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//壁粒子

					if(PART[i].type==FLUID)
					{
						//if(PART[i].materialID==1)  {red=0;green=0;blue=1;}
						//else if(PART[i].materialID==2)  {red=1;green=0;blue=0;}
						if(PART[i].r[A_X]>0)  {red=1;green=0;blue=0;}
						else {red=0;green=0;blue=1;}
					}
			    
					avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs2<<le<<" ";//粒子の大きさ出力

					avs2<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}  
			}
		}

		delete [] output;
		delete [] output2;
	}


	delete [] flag1;
	delete [] flag2;

	avs.close();
	avs2.close();
		////////////////////
}

//速度プロット関数
void plot_speed(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int fluid_number)
{
	
	double le=CON->get_distancebp()*0.5;
	double times=CON->get_speedtimes();
	int d=CON->get_dimention();
	int NUM;								//AVSに出力する粒子数
	int startID=0;							//最初に出力する粒子のid
	int num=0;								//数えあげ変数
	int face=CON->get_speed_face();			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON->get_speed_face_p();	//3D解析時のspeed.datの出力面の座標
	int d1,d2,d3;								//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	
	//AVS出力粒子数NUM計算
	if(CON->get_speed_plot_particle()==1) NUM=particle_number;	//全粒子出力
	else if(CON->get_speed_plot_particle()==2) NUM=fluid_number;//流体粒子のみ出力
	else if(CON->get_speed_plot_particle()==3)					//壁粒子のみ出力
	{
		NUM=particle_number; 
		startID=fluid_number;
	}
	
	ofstream vec("speed.dat");//絶対速度
	
	if(d==2)
	{
		for(int i=startID;i<NUM;i++)
		{
			vec<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
		}
	}
	else if(d==3)
	{
		//int d1,d2;				//出力に必要な次元
		if(face==0) {d1=A_Y; d2=A_Z; d3=A_X;}
		else if(face==1) {d1=A_X; d2=A_Z; d3=A_Y;}
		if(CON->get_ax_sym_modify()==OFF)
		{
			for(int i=startID;i<NUM;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];
					double z=PART[i].r[d2];
					double u=PART[i].u[d1];
					double w=PART[i].u[d2];
					vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
					if(x>xmax) xmax=x;
					if(z>ymax) ymax=z;
				}
			}
		}
		else if(CON->get_ax_sym_modify()==ON)
		{
			for(int i=startID;i<NUM;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];//出力に関与する軸 便宜上、変数名はxとなっているが、そうとは限らないことに注意
					double z=PART[i].r[d2];//出力に関与する軸
					double u=PART[i].u[d1];//出力に関与する軸
					double w=PART[i].u[d2];//出力に関与する軸

					double y=PART[i].r[d3];//出力に関与しない軸
					double v=PART[i].u[d3];//出力に関与しない軸

					double r=sqrt(x*x+y*y);//原点からの距離
			
					if(r>0)
					{
						double SIN,COS;
						if(x>0)
						{
							SIN=-y/r;
							COS=x/r;
						}
						if(x<0)
						{
							SIN=y/r;
							COS=-x/r;
						}
						double x2=COS*x-SIN*y;//回転後の座標　欲しいのはx2のみ。y2はいらない
						double u2=COS*u-SIN*v;//回転後の速度　欲しいのはu2のみ。v2はいらない
						vec<<x2<<" "<<z<<" "<<u2*times<<" "<<w*times<<endl;
						//if(SIN*x+COS*y!=0) cout<<SIN*x+COS*y<<endl;
						if(x2>xmax) xmax=x2;
						if(z>ymax) ymax=z;
					}
					else vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				}
			}
		}
	}
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON->get_legend_speed()>0) vec<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//最後に凡例出力
	vec.close();
	/////////////////////////////

	/////////重心に対する相対速度出力
	if(CON->get_relative_speed()==ON)
	{
		double U=0;								//平均速度
		double V=0;
		ofstream vec2("relative_speed.dat");

		///平均速度計算
		if(d==2) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[A_X]; V+=PART[i].u[A_Y];}
		else if(d==3) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[d1]; V+=PART[i].u[d2];}//d1,d2はすでにもとまっている
		if(fluid_number>0) {U/=fluid_number; V/=fluid_number;}
		/////////////////////

		if(d==2)
		{
			for(int i=0;i<fluid_number;i++)
			{
				double x=PART[i].r[A_X];
				double y=PART[i].r[A_Y];
				double u=PART[i].u[A_X]-U;
				double v=PART[i].u[A_Y]-V;
				vec2<<x<<" "<<y<<" "<<u*times<<" "<<v*times<<endl;
			}
		}
		else if(d==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];
					double z=PART[i].r[d2];
					double u=PART[i].u[d1]-U;
					double w=PART[i].u[d2]-V;
					vec2<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				}
			}
		}
		vec2.close();
	}/////////////////////


	if(CON->get_flat_speed_plot()==ON && d==3) //水平方向の速度出力
	{
		ofstream flat("flat_speed.dat");		
		double face_p2=CON->get_flat_speed_p();
		//for(int i=startID;i<NUM;i++)
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Z]>face_p2-le && PART[i].r[A_Z]<face_p2+le) flat<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
		}
		flat.close();
	}

	if(CON->get_speed_AVS()==ON && d==3)
	{
		num=0;
		for(int i=startID;i<NUM;i++)  num++;//if(PART[i].r[face]<face_p)
		
		ofstream fout2("speed_dist.fld");
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		fout2 << "dim1=" << num <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./speedvec filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./speedvec filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./speedvec filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./speedvec filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./speedvec filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./speedvec filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		ofstream fout("speedvec");
		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=startID;i<NUM;i++)
		{
			//if(PART[i].r[face]<face_p)
			{
				fout<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<" "<<PART[i].u[A_Z]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		fout.close();
	}
}

void plot_speed_each(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int fluid_number,int t)
{
	
	double le=CON->get_distancebp()*0.5;
	double times=CON->get_speedtimes();
	int d=CON->get_dimention();
	int NUM;								//AVSに出力する粒子数
	int startID=0;							//最初に出力する粒子のid
	int num=0;								//数えあげ変数
	int face=CON->get_speed_face();			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON->get_speed_face_p();	//3D解析時のspeed.datの出力面の座標
	int d1,d2,d3;								//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	
	//AVS出力粒子数NUM計算
	if(CON->get_speed_plot_particle()==1) NUM=particle_number;	//全粒子出力
	else if(CON->get_speed_plot_particle()==2) NUM=fluid_number;//流体粒子のみ出力
	else if(CON->get_speed_plot_particle()==3)					//壁粒子のみ出力
	{
		NUM=particle_number; 
		startID=fluid_number;
	}
	
	ofstream vec("speed.dat");//絶対速度
	
	if(d==2)
	{
		for(int i=startID;i<NUM;i++)
		{
			vec<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
		}
	}
	else if(d==3)
	{

		//int d1,d2;				//出力に必要な次元
		if(face==0) {d1=A_Y; d2=A_Z; d3=A_X;}
		else if(face==1) {d1=A_X; d2=A_Z; d3=A_Y;}
		if(CON->get_ax_sym_modify()==OFF)
		{
			for(int i=startID;i<NUM;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];
					double z=PART[i].r[d2];
					double u=PART[i].u[d1];
					double w=PART[i].u[d2];
					vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
					if(x>xmax) xmax=x;
					if(z>ymax) ymax=z;
				}
			}
		}
		else if(CON->get_ax_sym_modify()==ON)
		{
			for(int i=startID;i<NUM;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];//出力に関与する軸 便宜上、変数名はxとなっているが、そうとは限らないことに注意
					double z=PART[i].r[d2];//出力に関与する軸
					double u=PART[i].u[d1];//出力に関与する軸
					double w=PART[i].u[d2];//出力に関与する軸

					double y=PART[i].r[d3];//出力に関与しない軸
					double v=PART[i].u[d3];//出力に関与しない軸

					double r=sqrt(x*x+y*y);//原点からの距離
			
					if(r>0)
					{
						double SIN,COS;
						if(x>0)
						{
							SIN=-y/r;
							COS=x/r;
						}
						if(x<0)
						{
							SIN=y/r;
							COS=-x/r;
						}
						double x2=COS*x-SIN*y;//回転後の座標　欲しいのはx2のみ。y2はいらない
						double u2=COS*u-SIN*v;//回転後の速度　欲しいのはu2のみ。v2はいらない
						vec<<x2<<" "<<z<<" "<<u2*times<<" "<<w*times<<endl;
						//if(SIN*x+COS*y!=0) cout<<SIN*x+COS*y<<endl;
						if(x2>xmax) xmax=x2;
						if(z>ymax) ymax=z;
					}
					else vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				}
			}
		}
	}
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON->get_legend_speed()>0) vec<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//最後に凡例出力
	vec.close();

	xmax=-100;						//出力粒子の最大横座標
	ymax=-100;						//出力粒子の最大縦座標
	
	if(t==1 || t%10==0)
	{
		char filename[20];
		sprintf_s(filename,"speed%d.dat", t);
		ofstream vec2(filename);
		if(d==2)
		{
			for(int i=startID;i<NUM;i++)
			{
				vec2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
				if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
				if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
			}
		}
		else if(d==3)
		{
			//int d1,d2;				//出力に必要な次元
			if(face==0) {d1=A_Y; d2=A_Z; d3=A_X;}
			else if(face==1) {d1=A_X; d2=A_Z; d3=A_Y;}
			else if(face==2)	{d1=A_X; d2=A_Y; d3=A_Z;}
			if(CON->get_ax_sym_modify()==OFF)
			{
				for(int i=startID;i<NUM;i++)
				{
					if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
					{
						double x=PART[i].r[d1];
						double z=PART[i].r[d2];
						double u=PART[i].u[d1];
						double w=PART[i].u[d2];
						vec2<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
						if(x>xmax) xmax=x;
						if(z>ymax) ymax=z;
					}
				}
			}
			else if(CON->get_ax_sym_modify()==ON)
			{
				for(int i=startID;i<NUM;i++)
				{
					if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
					{
						double x=PART[i].r[d1];//出力に関与する軸 便宜上、変数名はxとなっているが、そうとは限らないことに注意
						double z=PART[i].r[d2];//出力に関与する軸
						double u=PART[i].u[d1];//出力に関与する軸
						double w=PART[i].u[d2];//出力に関与する軸

						double y=PART[i].r[d3];//出力に関与しない軸
						double v=PART[i].u[d3];//出力に関与しない軸

						double r=sqrt(x*x+y*y);//原点からの距離
			
						if(r>0)
						{
							double SIN,COS;
							if(x>0)
							{
								SIN=-y/r;
								COS=x/r;
							}
							if(x<0)
							{
								SIN=y/r;
								COS=-x/r;
							}
							double x2=COS*x-SIN*y;//回転後の座標　欲しいのはx2のみ。y2はいらない
							double u2=COS*u-SIN*v;//回転後の速度　欲しいのはu2のみ。v2はいらない
							vec2<<x2<<" "<<z<<" "<<u2*times<<" "<<w*times<<endl;
							//if(SIN*x+COS*y!=0) cout<<SIN*x+COS*y<<endl;
							if(x2>xmax) xmax=x2;
							if(z>ymax) ymax=z;
						}
						else vec2<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
					}
				}
			}
		}

		//凡例出力
		xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
		ymax+=4*le;

		if(CON->get_legend_speed()>0) vec2<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//最後に凡例出力
	
		vec2.close();///////////////////


		if(CON->get_model_number()==19&&CON->get_process_type()==2)
		{
			face_p=0.0;//ツール中央を通る速度を別途表示
			sprintf_s(filename,"speedn%d.dat", t);
			ofstream vec3(filename);
			if(d==2)
			{
				for(int i=startID;i<NUM;i++)
				{
					vec3<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
					if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
					if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
				}
			}
			else if(d==3)
			{
				//int d1,d2;				//出力に必要な次元
				if(face==0) {d1=A_Y; d2=A_Z; d3=A_X;}
				else if(face==1) {d1=A_X; d2=A_Z; d3=A_Y;}
				if(CON->get_ax_sym_modify()==OFF)
				{
					for(int i=startID;i<NUM;i++)
					{
						if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
						{
							double x=PART[i].r[d1];
							double z=PART[i].r[d2];
							double u=PART[i].u[d1];
							double w=PART[i].u[d2];
							vec3<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
							if(x>xmax) xmax=x;
							if(z>ymax) ymax=z;
						}
					}
				}
				else if(CON->get_ax_sym_modify()==ON)
				{
					for(int i=startID;i<NUM;i++)
					{
						if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
						{
							double x=PART[i].r[d1];//出力に関与する軸 便宜上、変数名はxとなっているが、そうとは限らないことに注意
							double z=PART[i].r[d2];//出力に関与する軸
							double u=PART[i].u[d1];//出力に関与する軸
							double w=PART[i].u[d2];//出力に関与する軸

							double y=PART[i].r[d3];//出力に関与しない軸
							double v=PART[i].u[d3];//出力に関与しない軸

							double r=sqrt(x*x+y*y);//原点からの距離
			
							if(r>0)
							{
								double SIN,COS;
								if(x>0)
								{
									SIN=-y/r;
									COS=x/r;
								}
								if(x<0)
								{
									SIN=y/r;
									COS=-x/r;
								}
								double x2=COS*x-SIN*y;//回転後の座標　欲しいのはx2のみ。y2はいらない
								double u2=COS*u-SIN*v;//回転後の速度　欲しいのはu2のみ。v2はいらない
								vec3<<x2<<" "<<z<<" "<<u2*times<<" "<<w*times<<endl;
								//if(SIN*x+COS*y!=0) cout<<SIN*x+COS*y<<endl;
								if(x2>xmax) xmax=x2;
								if(z>ymax) ymax=z;
							}
							else vec3<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
						}
					}
				}
			}

			//凡例出力
			xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
			ymax+=4*le;

			if(CON->get_legend_speed()>0) vec3<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//最後に凡例出力
	
			vec3.close();///////////////////
		}
	}
	/////////////////////////////
}

void plot_F(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *F[DIMENTION],int t)
{
	int d=CON->get_dimention();
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le;
	
	double Fz=0;
	double Fsum=0;
	 for(int i=0;i<fluid_number;i++)
    {
		Fz+=F[A_Z][i];
		Fsum+=sqrt(PART[i].F[A_X]*PART[i].F[A_X]+PART[i].F[A_Y]*PART[i].F[A_Y]+PART[i].F[A_Z]*PART[i].F[A_Z]);
	 }

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


	////ファイル出力
	ofstream fp("Fn.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();
	double times=CON->get_times()*CON->get_density()/CON->get_particle_mass();
	double cross_section=CON->get_speed_face_p();

    for(int i=0;i<fluid_number;i++)//流体節点のみ出力
    {
		//if(PART[i].r[A_Y]>-le*0.5&& PART[i].r[A_Y]<le*0.5) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F[A_X][i]*times<<" "<<F[A_Z][i]*times<<endl;
		if(PART[i].r[A_Y]>cross_section-le*0.5&& PART[i].r[A_Y]<cross_section+le*0.5) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].F[A_X]*times<<" "<<PART[i].F[A_Z]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
	}

	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	//凡例出力
	if(CON->get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力
	

    fp.close();///////////////////

	////ファイル出力//(スリットを通る断面)
	xmax=-100;						//出力粒子の最大横座標
	ymax=-100;						//出力粒子の最大縦座標
	ofstream fps("Fnslit.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=0;i<fluid_number;i++)//流体節点のみ出力
    {
		if(PART[i].r[A_Y]>-le*0.5+sin(PI/24)*PART[i].r[A_X] && PART[i].r[A_Y]<le*0.5+sin(PI/24)*PART[i].r[A_X])
		{
			fps<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<(PART[i].F[A_X]*cos(PI/24)+PART[i].F[A_Y]*sin(PI/24))*times<<" "<<PART[i].F[A_Z]*times<<endl;
		}
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		
	}

	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	//凡例出力
	if(CON->get_legend_F()>0) fps<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力
    fps.close();
	///////////////////

	////ファイル出力(Z断面)
	xmax=-100;						//出力粒子の最大横座標
	ymax=-100;						//出力粒子の最大縦座標
	ofstream fpz("Fnz.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=0;i<fluid_number;i++)//流体節点のみ出力
    {
		if(PART[i].r[A_Z]>-le*0.5+0.14125 && PART[i].r[A_Z]<le*0.5+0.14125) fpz<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].F[A_X]*times<<" "<<PART[i].F[A_Y]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
	}

	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	//凡例出力
	if(CON->get_legend_F()>0) fpz<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力
    
    fpz.close();///////////////////

	//if(t=1 || t%10==0)
	/////////////一定ステップごとに電磁力出力
	xmax=-100;						//出力粒子の最大横座標
	ymax=-100;						//出力粒子の最大縦座標
	if((t-1)%CON->get_EM_interval()==0)
	{
		char filename[20];
		sprintf_s(filename,"Fn%d.dat", t);
		ofstream fp(filename);

		for(int i=0;i<fluid_number;i++)//流体節点のみ出力
		{
			if(PART[i].r[A_Y]>-le*0.5&& PART[i].r[A_Y]<le*0.5) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].F[A_X]*times<<" "<<PART[i].F[A_Z]*times<<endl;	
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}

		//凡例出力
		xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
		ymax+=4*le;

		if(CON->get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力
	
		fp.close();///////////////////
	}
	
}

//粒子情報を出力する関数
void output_alldata_AVS(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double Umax,double mindis,int t,double TIME, unsigned int time0,int count_avs)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」

	if(t==1)
	{
		ofstream fout("alldata_part.dat");

		fout<<"# Micro AVS Geom:2.00"<<endl;
		fout<<CON->get_step()/CON->get_interval()+1<<endl;//microAVSに出力する総ｽﾃｯﾌﾟ数。ファイル出力はCON->get_interval()回に1回と最初に行う。
		
		fout.close();
	}


	ofstream fp("alldata_part.dat",ios :: app);
	fp<<"step"<<t/CON->get_interval()+1<<endl;
	fp<<"sphere"<<endl;
	fp<<"time="<<TIME<<endl;
	fp<<"color"<<endl;
	
	 
	//節点番号とその座標の出力 
	//for(int i=0;i<particle_number;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
	/*/要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=0;i<nelm;i++)
	{
		fp<<i+1<<"  0 tri ";
		for(int j=0;j<3;j++)	fp<<ELEM[i].node[j]+1<<" ";
		fp<<endl;
	}
	////*/


	////*/

	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	//fp<<"6 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	//fp<<"2 1 1"<<endl;	//この行の詳細はヘルプを参照
	//fp<<"6 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	//fp<<"6 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	////fp<<"first_id"<<endl;
	//fp<<"speed,m/s"<<endl;
	//fp<<"P,N/m^2"<<endl;
	//fp<<"h,W/m^3"<<endl;
	//fp<<"surface"<<endl;
	//fp<<"color"<<endl;//粒子の表示色

	if(t==1) for(int i=0;i<particle_number;i++) PART[i].firstID=i;//最初の状態での粒子番号を記憶

	//各節点の情報値入力
	for(int i=0;i<particle_number;i++)
	{
		int p=i;
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
		
		if(PART[i].type==FLUID) fp<<" "<<PART[i].firstID<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].u[A_X]<<" "<<PART[i].u[A_Y]<<" "<<PART[i].u[A_Z]<<" "<<PART[i].P<<" "<<" "<<PART[i].T<<" "<<PART[i].surface<<" "<<PART[i].type<<" "<<PART[i].color<<endl;
		else fp<<" "<<PART[i].firstID<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].u[A_X]<<" "<<PART[i].u[A_Y]<<" "<<PART[i].u[A_Z]<<" "<<PART[i].P<<" "<<" "<<PART[i].T<<" "<<PART[i].surface<<" "<<PART[i].type<<" "<<3<<endl;
	}

	/*fp<<"3 1 1 1"<<endl;
	fp<<"potential,V"<<endl;
	fp<<"En,V/m"<<endl;
	fp<<"Fn,N/m^2"<<endl;
	
	//要素情報出力　要素情報はmicroAVSの可視化メソッドバー内の、「要素データの塗りつぶし」を押せば見れる
	for(int i=0;i<nelm;i++) fp<<i+1<<"  "<<ELEM[i].potential<<" "<<ELEM[i].En<<" "<<ELEM[i].Fn<<endl;
	*/
	cout<<"OK"<<endl;
	fp.close();
}

//ｸｰﾗﾝ数
void culan(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double *dt,double mindis,double Umax,double *g)
{
	double C=CON->get_curan();		//ｸｰﾗﾝ条件数
	double le=mindis;				//最短粒子間距離
	//double le=CON->get_distancebp()/sqrt(2.0);	//平均粒子間距離

	///////////クーラン数//////////////////
	if(C>0)
	{    
		double newdt=*dt;	//新しいdt
		
		if(Umax!=0) newdt=C*le/Umax;
		if(newdt>CON->get_dt()) *dt=CON->get_dt();
		else *dt=newdt;

		//if(*dt<=CON->get_dt()/10) *dt=CON->get_dt()/10;

		if(*dt!=CON->get_dt()) cout<<"ｸｰﾗﾝ数制限 dt="<<*dt<<endl;
		else cout<<"設定値　dt="<<*dt<<endl;
	}

	///拡散数の正確な定義を調べて書きなおせ
	if(CON->get_vis()!=0 && CON->get_vis_calc_type()==POSITIVE)
	{
		if(*dt>0.25*le*le/CON->get_vis()) cout<<"拡散数違反"<<endl;
	}
	/////////////////////////////
	
}
/*//ｸｰﾗﾝ数
void culan(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double *dt,double mindis,double Umax,double *g)
{
	double C=CON->get_curan();		//ｸｰﾗﾝ条件数
	double le=mindis;				//最短粒子間距離

	///////////クーラン数//////////////////
	if(C>0)
	{      
		double newdt=*dt;	//新しいdt
		
		if(Umax!=0) newdt=C*le/Umax;
		if(newdt>CON->get_dt()) *dt=CON->get_dt();
		else *dt=newdt;
		
		if(*dt!=CON->get_dt()) cout<<"ｸｰﾗﾝ数制限 dt="<<*dt<<endl;
	}

	///拡散数の正確な定義を調べて書きなおせ
	if(CON->get_vis()!=0 && CON->get_vis_calc_type()==POSITIVE)
	{
		if(*dt>0.25*le*le/CON->get_vis()) cout<<"拡散数違反"<<endl;
	}
	/////////////////////////////
	
}
///*/

///ラプラシアン計算関数
void u_laplacian_f(mpsconfig *CON,vector<mpsparticle> &PART,double *laplacian[DIMENTION],double n0,double lamda,int fluid_number,double dt)
{
	cout<<"速度拡散項計算-------";

	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;						//ﾗﾌﾟﾗｼｱﾝ用影響半径
	int d=CON->get_dimention();						//次元
	unsigned int timeA=GetTickCount();				//計算開始時刻

	///拡散数に違反していないか確認
	double limit=0.25;
	if(d==3) limit=0.125;
	if(dt>limit*le*le/CON->get_vis()) cout<<"拡散数違反 dt<"<<limit*le*le/CON->get_vis()<<" ";
	/////////////*/
   
	for(int i=0;i<fluid_number;i++)
	{
		double lam=0;	//各粒子の正確なλ
		double W=0;	//各粒子の正確な粒子数密度
		for(int D=0;D<DIMENTION;D++) laplacian[D][i]=0.0;//初期化
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k]; 
			if(CON->get_wall_adheision()==0)//ﾌﾘｰｽﾘｯﾌﾟ
			{
				if(PART[j].type==FLUID)
				{ 
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);//粒子間の距離
					
					double w=kernel(R,dis);//重み関数
					W+=w;
					lam+=dis*dis*w;
					if(CON->get_laplacian()==0 || CON->get_laplacian()==1)
					{
						for(int D=0;D<DIMENTION;D++)
						{
							double u=PART[j].u[D]-PART[i].u[D];
							laplacian[D][i]+=u*w;
						}
					}
					else if(CON->get_laplacian()==2)
					{ 
						laplacian[A_X][i]+=(PART[j].u[A_X]-PART[i].u[A_X])*w/(dis*dis);
						laplacian[A_Y][i]+=(PART[j].u[A_Y]-PART[i].u[A_Y])*w/(dis*dis);
						laplacian[A_Z][i]+=(PART[j].u[A_Z]-PART[i].u[A_Z])*w/(dis*dis);
					}
				}
			}
			else if(CON->get_wall_adheision()==1 || CON->get_wall_adheision()==2)//ﾉﾝｽﾘｯﾌﾟ
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);//粒子間の距離
				
				double w=kernel(R,dis);//重み関数
				W+=w;
				lam+=dis*dis*w;
				if(CON->get_laplacian()==0 || CON->get_laplacian()==1)
				{
					for(int D=0;D<DIMENTION;D++)
					{
						double u=PART[j].u[D]-PART[i].u[D];
						if(CON->get_wall_adheision()==2)
						{   //本来の定義どうりのﾉﾝｽﾘｯﾌﾟ
							if(PART[j].type==INWALL || PART[j].type==OUTWALL) u*=2;
						}
						laplacian[D][i]+=u*w;
					}
				}
				else if(CON->get_laplacian()==2)
				{ 
					for(int D=0;D<DIMENTION;D++)
					{
						double u=PART[j].u[D]-PART[i].u[D];
						if(CON->get_wall_adheision()==2)
						{   //本来の定義どうりのﾉﾝｽﾘｯﾌﾟ
							if(PART[j].type==INWALL || PART[j].type==OUTWALL) u*=2;
						}
						laplacian[D][i]+=u*w/(dis*dis);
					}
				}
			}
		}
		if(W!=0) lam/=W;
		else if(W==0) lam=lamda;
		for(int D=0;D<DIMENTION;D++)
		{        
			if(CON->get_laplacian()==0)	 laplacian[D][i]=laplacian[D][i]*2*d/(n0*lamda);
			if(CON->get_laplacian()==1 &&W!=0) laplacian[D][i]=laplacian[D][i]*2*d/(W*lam);
			if(CON->get_laplacian()==2 &&W!=0) laplacian[D][i]=laplacian[D][i]*2*d/W;
		}
	}
   ///////*/

	cout<<"ok  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
   
}

///粘性項陰解析関数
void visterm_negative(mpsconfig *CON,vector<mpsparticle> &PART,double *laplacian[DIMENTION],double n0,double lamda,int fluid_number,int particle_number,double dt,int t)
{
	//現在は3Ｄのみ対応 //20130505 修正。2Dにも適用可能のはず
	//解くべき式は	ui(t+1)-ui(t)=Δt*v*(2d/λn0)Σ(uj-ui)w
	//				⇔Δt*v*(2d/λn0)Σ(uj-ui)w-ui(t+1)=-ui(t)
	//				⇔Σ(uj-ui)w-λn0/(2vdΔt)ui(t+1)=-λn0/(2vdΔt)ui(t)
	//				⇔vΣ(uj-ui)w-λn0/(2dΔt)ui(t+1)=-λn0/(2dΔt)ui(t)

	cout<<"速度拡散項計算---";

	double R=CON->get_re2()*CON->get_distancebp();	//ﾗﾌﾟﾗｼｱﾝ用影響半径
	int d=CON->get_dimention();						//次元
	unsigned int timeA=GetTickCount();				//計算開始時刻
	int count=0;
	//int pn=fluid_number*d;							//未知数:粒子数×速度D(次元)成分 現在は３Ｄのみ対応
	int pn=fluid_number;
	double co=lamda*n0/(2*d*dt);///CON->get_vis();	//計算によく現れる係数
	double co2=2*d*dt/(lamda*n0)*CON->get_vis();
	double vis0=CON->get_vis();

	double *vis=new double [fluid_number];			//各粒子の動粘性項計算
	double *B   = new double[pn];					//解行列

	double *B_vec[DIMENTION];
    for(int D=0;D<DIMENTION;D++) B_vec[D]=new double [fluid_number];//各方向の解行列

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			laplacian[D][i]=0;		//初期化
			B_vec[D][i]=0;
		}
	}
	///各粒子の動粘性を計算

	if(CON->get_model_number()==19)
	{
		calc_vis_value(CON,PART,fluid_number,vis,dt,t,particle_number);//FSWモデルの場合
		for(int i=0;i<particle_number;i++)	PART[i].vis=0;
		for(int i=0;i<fluid_number;i++)	PART[i].vis=vis[i];
		if(t==0||t%CON->get_interval()==0)	output_viscousity_avs(CON,PART,t,particle_number,fluid_number);
	}
	else for(int i=0;i<fluid_number;i++) vis[i]=vis0;

	if(CON->get_temperature_depend()==ON) calc_physical_property(CON,PART,fluid_number,vis,particle_number,5);//動粘性の温度依存
	//////////////////////////////////////////////////////////////////

	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) laplacian[D][i]=0;//初期化

	//解行列B[]作成
	for(int i=0;i<fluid_number;i++)//未知数の順番はPART[i].u[A_X],PART[i+1].u[A_X]・・・PART[i].u[A_Y],PART[i+1].u[A_Y],
	{
		for(int D=0;D<d;D++)
		{
			B_vec[D][i]=-co*PART[i].u[D]*PART[i].PND2/n0;//n0ではなくwi使用
		}
	}
	/////*/

	int number=0;			//係数行列の非ゼロ要素数
	for(int i=0;i<fluid_number;i++)
	{
		int num=1;//粒子i周辺の流体粒子数 初期値が1なのは自分自身をｶｳﾝﾄしているから
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k];
			if(j<fluid_number) num++;
		}
		//
		number+=num;
	}///numberがもとまった
	
    double *val = new double [number];
	int *ind = new int [number];//非ゼロ要素の列番号格納配列
	int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
	
	/////////////////////val,ind ,ptrに値を格納
	int index=0;
	for(int n=0;n<pn;n++)
	{
		//int D=n/fluid_number;//n番目の未知数はu[D]に関する未知数　流体が存在しないモデルでこの関数を回さないよう注意
	    
		ptr[n]=index;
		int i=n;//ｎ番目の未知数に該当する粒子番号
		
		int KK=index;
		double AA=0;
	    ind[index]=n;
	    index++;

	    for(int k=0;k<PART[i].N2;k++)
	    {
	        int j=PART[i].NEI2[k]; 
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double w=kernel(R,dis);

			if(j<fluid_number)//粒子jが流体なら
			{
				double v=2*vis[i]*vis[j]/(vis[i]+vis[j]);
				val[index]=w*v;
				ind[index]=j;
				index++;
				AA+=w*v;
			}
			else
			{
				for(int D=0;D<d;D++) B_vec[D][n]-=w*PART[j].u[D]*vis[i];//粒子jが流体でないなら
				AA+=w*vis[i];
			}
	    }
		//val[KK]=-AA-co;
		val[KK]=-AA-co*PART[i].PND2/n0;//n0ではなくwi使用
	}
	ptr[pn]=number;//この最後の行がないと、任意のnでfor(int j=ptr[n];j<ptr[n+1];j++) のようなことができない
	////////////////////*/	 


	///ICCG法でとく場合、係数行列をきちんと並び変えしなくてはならない
	if(CON->get_vis_solver()==1)
	{
		//ここで作られる行列は対角成分がすべて負なので、これを正にするために、係数行列と解行列に-1をかける
		for(int n=0;n<pn;n++)
		{
			for(int j=ptr[n];j<ptr[n+1];j++) val[j]*=-1;
			for(int D=0;D<d;D++) B_vec[D][n]*=-1;
		}

		//#pragma omp parallel for
		for(int n=0;n<pn;n++)
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
	}//////////////////////////////////*/

	double *r=new double[pn];
	double *X=new double[pn];		//行列の答え格納
	double *AP = new double [pn];
	double *P = new double [pn];

	for(int D=0;D<d;D++)
	{
		count=0;
		for(int n=0;n<pn;n++) B[n]=B_vec[D][n];
		/////////////////////////初期値//////////////////
		if(CON->get_initial_u()==OFF)
		{
			for(int n=0;n<pn;n++) 
			{
				 X[n]=0;
				 r[n]=B[n];
				 P[n]=r[n];
			}
		}
		else if(CON->get_initial_u()==ON) //初期値を与える場合　そんなに効果ない感じ
		{
			for(int n=0;n<pn;n++)
			{
				X[n]=PART[n].u[D];
			}
			for(int n=0;n<pn;n++)
			{
				double AX=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AX+=val[j]*X[ind[j]];
				r[n]=B[n]-AX;
				P[n]=r[n];
			}
		}
		//////////////////////////////////////////////

		
		//if(CON->get_vis_solver()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,1e-8); //CG法により行列を解く
		if(CON->get_vis_solver()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,1e-8); //CG法により行列を解く
		else if(CON->get_vis_solver()==1) iccg(CON,val,ind,ptr,pn,B,number,X,r,P,1e-8,&count);//ICCGによる行列計算開始
		else if(CON->get_vis_solver()==2) MRTR(CON,r, pn,X,&count,1e-8,val,ind,ptr);

		for(int n=0;n<pn;n++)
		{	
			//PART[i].u[D]=XX[n];		//こっちだと速度をここで決めてしまう
			laplacian[D][n]=(X[n]-PART[n].u[D])/(dt*vis0);//こっちだとrenewal関数ないで速度更新する形になる
		}
		cout<<"反復回数:"<<count<<"  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		if(count==1) for(int i=0;i<fluid_number;i++) laplacian[D][i]=0;//反復回数が１というのは明らかにエラーなので、その値は使用しない
	}

	delete [] val;
	delete [] ind;
	delete [] ptr;
    delete [] B;

	delete [] r;
	delete [] X;
	delete [] AP;
	delete [] P;

	delete [] vis;
	for(int D=0;D<DIMENTION;D++) delete [] B_vec[D];
}

//CG法
void CG_method(mpsconfig *CON,double *r,double *P,double *AP,double *val,int *ind,int *ptr,int pn,double *X,int *countN,double EP)
{
	cout<<"CG法スタート------";
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double rr=0;
	double E=1;//誤差
	double alp,beta;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];
	
	if(CON->get_omp_P()==OFF)//通常版
	{
		while(E>EP)// EP=CON->get_CGep();//収束判定(convergence test)
		{
			count++;
			//////////////alpを求める
			for(int n=0;n<pn;n++)
			{      
				AP[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
			}
			double PAP=0;
			for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
			alp=rr/PAP;
		//	cout<<"alp="<<alp<<" rr="<<rr<<" PAP="<<PAP<<endl;
			//////////////////////
		
			//////////////// 解更新　X(k+1)=X(k)+alp*P
			for(int n=0;n<pn;n++) X[n]+=alp*P[n];
			//////////////////////////////
			
			//////////////// r=r-alp*AP
			for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
			/////////////////////////////
			
			///////////////////////beta
			beta=1.0/rr;
			rr=0;
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			beta=beta*rr;
			///////////////////////

			//////////////////誤差
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
			////////////////////////
			
			///////////////////// P=r+beta*P
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];
		}
	}
	else if(CON->get_omp_P()==ON)//openMPを使用する場合
	{
		while(E>EP)
		{
			count++;
			//////////////alpを求める
			double PAP=0;
			#pragma omp parallel for reduction(+:PAP)
			for(int n=0;n<pn;n++)
			{
				AP[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
				PAP+=P[n]*AP[n];
			}
			alp=rr/PAP;
			//////////////////////

			//////////////
			E=0;//誤差
			beta=1.0/rr;
			rr=0;
			#pragma omp parallel for reduction(+:rr)
			for(int n=0;n<pn;n++) 
			{
				X[n]+=alp*P[n];// 解更新　X(k+1)=X(k)+alp*P
				r[n]-=alp*AP[n];// r=r-alp*AP
				rr+=r[n]*r[n];
			}
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
		
			beta=beta*rr;///beta
			
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];/// P=r+beta*P
		}
	}
	*countN=count;//反復回数を渡す
}

///ICCG法
void iccg(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X,double *r,double *P,double EP,int *count2)
{
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	double accel=0.87;//CON->get_CGaccl();//加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *val2=new double [num2];
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
	double E=1;//誤差
	double *AP = new double [pn];
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

	cout<<"ICCG法:未知数="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=EP;//収束判定
	rLDLt_r=0;
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//最初のrLDLt_rだけここで求める
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alpを求める
		double PAP=0;
		#pragma omp parallel for reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		//cout<<"alp="<<alp<<endl;
		//////////////////////
		E=0;
		#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++)
		{
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
			r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
			E+=r[n]*r[n];						//更新後の誤差
		}
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
		
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
	//cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
		
	delete [] AP;

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
	*count2=count;//反復回数を格納して返す
}

///粘性項計算関数
void calc_viscous_term(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double N0,double *laplacian[DIMENTION],double lamda,int t)
{
	int d=CON->get_dimention();

	if(CON->get_nensei()>0)
	{
		if(CON->get_vis_calc_type()==POSITIVE) u_laplacian_f(CON,PART,laplacian,N0,lamda,fluid_number,dt);//陽解法
		if(CON->get_vis_calc_type()==NEGATIVE) //陰解法
		{
			int flag=ON;
			if(fluid_number==0) flag=OFF;			//流体が存在しなければ計算しない
			if(t==1 && CON->get_restart()==OFF)//全粒子速度が0の場合、ＣＧ法でalpha=∞となりエラーになる。それを避けたい
			{
				flag=OFF;
				for(int i=0;i<particle_number;i++) for(int D=0;D<d;D++) if(PART[i].u[D]!=0) flag=ON;
			}
			else if(t==1 && CON->get_restart()==ON && CON->get_set_zero_speed()==ON) flag=OFF;
			if(flag==ON)  visterm_negative(CON,PART,laplacian,N0,lamda,fluid_number,particle_number,dt,t);
			if(flag==OFF) for(int D=0;D<CON->get_dimention();D++) for(int i=0;i<fluid_number;i++) laplacian[D][i]=0;//初期化
		}
	}
	else for(int D=0;D<d;D++) for(int i=0;i<fluid_number;i++) laplacian[D][i]=0;//計算しない場合も初期化だけしておく		
}







//ガウスの消去法 解は最終的にBのなかへ
void gauss(double *matrix,double *B,int N)
{
	for(int k=0;k<N;k++)
	{
		double akk=matrix[k*N+k];
		
		for(int i=0;i<N;i++)
		{
			if(i!=k)
			{
				double A=matrix[i*N+k]/akk;
				//for(int j=0;j<N;j++)
				for(int j=k;j<N;j++)
				{					
					matrix[i*N+j]-=A*matrix[k*N+j];
				}
				B[i]-=A*B[k];
			}
		}
	}
	for(int k=0;k<N;k++) B[k]/=matrix[k*N+k];

}

/////仮の速度および位置決定
void renewal_u_and_r_in_positive(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double dt,double *Umax,double **potential,double **laplacian,double *g,double **previous_Un,double **F)
{

	double U=0;						//最大速度
	double vis=CON->get_vis();
	double *vis_T=new double [fluid_number];
	double mass=CON->get_particle_mass();	//粒子の質量
	double *mass_T=new double [fluid_number];	//粒子の質量
	int d=CON->get_dimention();
	int sw=CON->get_temporary_r();	//ＯＮなら仮の位置を計算する
	double T0=CON->get_initialT();//基準温度
	double alp=CON->get_CTE();//線膨張係数
	

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			mass_T[i]=mass;//温度依存性の考慮の有無に関わらず、ひとまずconfigの値に初期化する。calc_physical_propertyに行っても変更されない場合があるため。
			vis_T[i]=vis;
		}
	}

	double *old_U[DIMENTION];
	for(int D=0;D<DIMENTION;D++) old_U[D]=new double [fluid_number];//変更前の速度を記憶しておく

	if(t==1) for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENTION;D++) previous_Un[D][i]=0;//t=1のときは初期化   

	if(CON->get_temperature_depend()==ON)
	{
		calc_physical_property(CON,PART,fluid_number,mass_T,fluid_number,1);//密度の温度依存　この時点でmass_Tに入るのは密度(kg/m3)
		for(int i=0;i<fluid_number;i++) mass_T[i]=mass*mass_T[i]/CON->get_density();
		//for(int i=0;i<fluid_number;i++) mass_T[i]=mass;

		//calc_physical_property(CON,PART,fluid_number,vis_T,fluid_number,5);//動粘性の温度依存 粘性項の陰解析にて、
	}
		
			

	//potential[D][i]を場合によってはゼロに初期化する
	if(CON->get_dir_for_P()==1 || CON->get_dir_for_P()==3) //表面粒子の表面張力は圧力値として計算されているので,ここでは考慮しないよう初期化する
	{
		for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) potential[D][i]=0;
	}
	//////////////////////*/

	/////////////速度更新
	for(int i=0;i<fluid_number;i++)
	{        
		double speed=0;//粒子速度
		for(int D=0;D<d;D++)
		{   
			old_U[D][i]=PART[i].u[D];
			
			//if(CON->get_temperature_depend()==OFF) PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]+F[D][i]/mass);
			//if(CON->get_temperature_depend()==ON) PART[i].u[D]+=dt*(vis_T[i]*laplacian[D][i]+potential[D][i]+g[D]+F[D][i]/mass_T[i]);
			if(CON->get_temperature_depend()==OFF) 
			{
				if(CON->get_buoyant()==OFF) PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]+PART[i].F[D]/mass);
				if(CON->get_buoyant()==ON)//浮力のブシネスク近似。
				{
					double T0=CON->get_initialT();//基準温度
					double alp=CON->get_CTE();//線膨張係数
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1+(mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0))+PART[i].F[D]/mass);//これが正しいはず
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass_T[i])+PART[i].F[D]/mass);//おかしい？massTで割る根拠がない。
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass)+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
				}

					
					
			}
			if(CON->get_temperature_depend()==ON) 
			{
				if(CON->get_buoyant()==OFF)//ブシネスク近似を用いない。密度変化があると連続の式が速度発散0にはならなくなるからこちらを使うならいろいろ直す必要あり
				{
					PART[i].u[D]+=dt*(vis_T[i]*laplacian[D][i]+potential[D][i]+g[D]+PART[i].F[D]/mass_T[i]);
				}
				if(CON->get_buoyant()==ON)//浮力のブシネスク近似。
				{
					double T0=CON->get_initialT();//基準温度
					double alp=CON->get_CTE();//線膨張係数
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1+(mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0))+PART[i].F[D]/mass);//これが正しいはず
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass_T[i])+PART[i].F[D]/mass);//おかしい？massTで割る根拠がない。
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass)+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
				}
			}
			
			//PART[i].u[D]=previous_Un[D][i]+dt*(vis*laplacian[D][i]+potential[D][i]+g[D]);//蛙とび法
			speed+=PART[i].u[D]*PART[i].u[D];
		}
		if(speed>U) U=speed;	
	}
	*Umax=U;


	//表面だったら速度の更新をやめる
	//if(CON->get_fix_surface()==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].surface==ON) PART[i].u[D]=0;

	//位置更新
	if(sw==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*0.5*(PART[i].u[D]+old_U[D][i]);//台形則
	
	//previous_Un(1step前の速度情報)の値を更新
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++)previous_Un[D][i]=old_U[D][i];
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++)previous_Un[D][i]=PART[i].u[D];
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++)previous_Un[D][i]=0;
	
	
	if(CON->get_buoyant()==ON)
	{
		ofstream vec("gf.dat");//絶対速度
		double le=CON->get_distancebp()*0.5;
		double times=0.1;
		int d=CON->get_dimention();
		int NUM=0;								//AVSに出力する粒子数
		int startID=0;							//最初に出力する粒子のid
		int num=0;								//数えあげ変数
		int face=CON->get_speed_face();			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ 2=XY
		double face_p=CON->get_speed_face_p();	//3D解析時のspeed.datの出力面の座標
		int d1,d2,d3;								//3D解析時の出力に必要な次元
		double xmax=-100;						//出力粒子の最大横座標
		double ymax=-100;						//出力粒子の最大縦座標
	
		//AVS出力粒子数NUM計算
			
		NUM=fluid_number;//流体粒子のみ出力
			
	
			
	
		if(d==2)
		{
			for(int i=startID;i<NUM;i++)
			{
				vec<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<0<<" "<<-g[A_Y]*alp*(PART[i].T-T0)<<endl;//浮力項を出力
				if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
				if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
			}
		}
		else if(d==3)
		{
			//int d1,d2;				//出力に必要な次元
			if(face==0) {d1=A_Y; d2=A_Z; d3=A_X;}
			else if(face==1) {d1=A_X; d2=A_Z; d3=A_Y;}
						
			for(int i=startID;i<NUM;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];
					double z=PART[i].r[d2];
					double u=PART[i].u[d1];
					double w=PART[i].u[d2];
					vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
					if(x>xmax) xmax=x;
					if(z>ymax) ymax=z;
				}
			}
						
		}
		xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
		ymax+=4*le;
		vec<<xmax<<" "<<ymax<<" "<<0<<" "<<0.1*g[A_Y]*times<<endl;//最後に凡例出力
		vec.close();
	}
	for(int D=0;D<DIMENTION;D++) delete [] old_U[D];
	delete [] vis_T;
	delete [] mass_T;
	
}

///速度発散計算関数
double divergence(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0)
{
    double W=0;										//粒子数密度
    double R=CON->get_distancebp()*CON->get_re();	//影響半径
    double div=0;									//発散の値

	for(int k=0;k<PART[i].N;k++)
    {    
        int j=PART[i].NEI[k]; 
        double X=PART[j].r[A_X]-PART[i].r[A_X];
		double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
		double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
		double dis=sqrt(X*X+Y*Y+Z*Z);
			       
		double w=kernel(R,dis);
		
		div+=(PART[j].u[A_X]-PART[i].u[A_X])*X*w/(dis*dis);
		div+=(PART[j].u[A_Y]-PART[i].u[A_Y])*Y*w/(dis*dis);
		div+=(PART[j].u[A_Z]-PART[i].u[A_Z])*Z*w/(dis*dis);
		W+=w;
    }
    if(W!=0)
	{
		div*=CON->get_dimention()/W;
	}
    return div;
}

///速度発散計算関数(WLSM法)
double divergence2(mpsconfig *CON,vector<mpsparticle> &PART,int i,int surface_sw)
{
	//surface_sw: ONなら表面を無視する、OFFなら表面も考慮にいれた発散を計算する。ON,OFFの定義を混乱しないように
    double div=0;//発散の値
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N=0;					//係数行列の元
	int order=CON->get_divU_order();				//近似曲面のｵｰﾀﾞｰ。 1=線形 2=二次
    
	//係数行列の大きさの決定
	if(d==2)
	{
		if(order==1) N=2;
		else if(order==2) N=5;
		else if(order==3) N=9;
	}
	else if(d==3)
	{
		if(order==1) N=4;
		else if(order==2) N=10;
	}
	////////////////////////////////

	double *matrix=new double [N*N];	//N×Nの係数行列
	double *matrix2=new double [N*N];	//matrixのバックアップ
	double *B1=new double [N];			//Nの解行列
	double *B2=new double [N];			//Nの解行列
	double *B3=new double [N];			//Nの解行列

	for(int n=0;n<N*N;n++) matrix[n]=0;	//初期化
	for(int n=0;n<N;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}

	if(d==2 && order==1)				//二次元
	{
		if(PART[i].N>=2)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// leで割るのは正規化のため //　Lは粒子の基準長さ(可変解像度関係)
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
				double U=(PART[j].u[A_X]-PART[i].u[A_X])/PART[i].L;
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y])/PART[i].L;
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				if(dis>1) w=1/(dis*dis*dis*dis);
				//if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);

				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//表面粒子の速度を無視する
					
				matrix[0]+=X*X*w;			//ΣXjwj
				matrix[1]+=X*Y*w;		//ΣXjYjwj
				matrix[3]+=Y*Y*w;			//ΣYjwj
				
				B1[0]+=U*X*w;//ΣujXjwj
				B1[1]+=U*Y*w;//ΣujYjwj
				B2[0]+=V*X*w;//ΣvjXjwj
				B2[1]+=V*Y*w;//ΣvjYjwj
			}
			
			matrix[2]=matrix[1];		//ΣXjYjwj

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//値を保存

			gauss(matrix,B1,N);//ガウスの消去法で解く

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//ガウスの消去法で解く

			double a_u,b_u;
			double a_v,b_v;
			a_u=B1[0]; a_v=B2[0];//X方向微分
			b_u=B1[1]; b_v=B2[1];//Y方向微分
			
			div=(B1[0]+B2[1]);
			if(div+div!=2*div) cout<<"速度の発散が非実数 i="<<i<<endl;

			/*/誤差を計算
			double Q1=0; double Q2=0;
			int count_for_Q=0;
			
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					double X=(PART[j].r[A_X]-PART[i].r[A_X]);
					double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
					double dis=sqrt(X*X+Y*Y);
					
					double w=1;
					if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
					double dQ1=(a_u*X+b_u*Y+PART[i].u[A_X]-PART[j].u[A_X]);
					double dQ2=(a_v*X+b_v*Y+PART[i].u[A_Y]-PART[j].u[A_Y]);
					Q1+=dQ1*dQ1*w;
					Q2+=dQ2*dQ2*w;
					count_for_Q++;
				}
			}
			
			if(count_for_Q>0) {Q1/=count_for_Q; Q2/=count_for_Q;}
			//if(fabs(PART[i].u[A_X])>1e-6) Q1/=fabs(PART[i].u[A_X]);
			//if(fabs(PART[i].u[A_Y])>1e-6) Q2/=fabs(PART[i].u[A_Y]);
			cout<<i<<" "<<Q1<<" "<<Q2<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<endl;
			///*/	
		}
		else div=0;	//周辺粒子が少なすぎる場合はゼロにすればよい。どうせ圧力勾配も計算されないから。
	}
	if(d==2 && order==2)//二次元2次式
	{
		if(PART[i].N>=4)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// leで割るのは打ち切り誤差防止
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
				double U=(PART[j].u[A_X]-PART[i].u[A_X])/PART[i].L;
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y])/PART[i].L;
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				if(dis>1) w=1/(dis*dis*dis*dis);
				//if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//表面粒子の速度を無視する
					
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

				
				B1[0]+=U*X*w;//ΣujXjwj
				B1[1]+=U*Y*w;//ΣujYjwj
				B1[2]+=U*X*X*w;//ΣujXj^2wj
				B1[3]+=U*X*Y*w;//ΣujXjYjwj
				B1[4]+=U*Y*Y*w;//ΣujYj^2wj

				B2[0]+=V*X*w;//ΣvjXjwj
				B2[1]+=V*Y*w;//ΣvjYjwj
				B2[2]+=V*X*X*w;//ΣvjXj^2wj
				B2[3]+=V*X*Y*w;//ΣvjXjYjwj
				B2[4]+=V*Y*Y*w;//ΣvjYj^2wj
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

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//値を保存
			gauss(matrix,B1,N);//ガウスの消去法で解く
			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//ガウスの消去法で解く

			double a_u,b_u,c_u,d_u,e_u;
			double a_v,b_v,c_v,d_v,e_v;
			a_u=B1[0]; a_v=B2[0];//X方向微分
			b_u=B1[1]; b_v=B2[1];//Y方向微分
			c_u=B1[2]; c_v=B2[2];
			d_u=B1[3]; d_v=B2[3];
			e_u=B1[4]; e_v=B2[4];
			
			div=(B1[0]+B2[1]);
			if(div+div!=2*div) cout<<"速度の発散が非実数 i="<<i<<endl;
			//cout<<a_u<<" "<<b_u<<" "<<a_v<<" "<<b_v<<endl;

			/*/誤差を計算
			double Q1=0; double Q2=0;
			int count_for_Q=0;
			
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type!=OUTWALL)
				{
					double X=(PART[j].r[A_X]-PART[i].r[A_X]);
					double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
					double dis=sqrt(X*X+Y*Y);
					
					double w=1;
					if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
					double dQ1=(a_u*X+b_u*Y+c_u*X*X+d_u*X*Y+e_u*Y*Y+PART[i].u[A_X]-PART[j].u[A_X]);
					double dQ2=(a_v*X+b_v*Y+c_v*X*X+d_v*X*Y+e_v*Y*Y+PART[i].u[A_Y]-PART[j].u[A_Y]);
					Q1+=dQ1*dQ1*w;
					Q2+=dQ2*dQ2*w;
					count_for_Q++;
				}
			}
			
			if(count_for_Q!=0) {Q1/=count_for_Q; Q2/=count_for_Q;}
			//if(fabs(PART[i].u[A_X])>1e-6) Q1/=fabs(PART[i].u[A_X]);
			//if(fabs(PART[i].u[A_Y])>1e-6) Q2/=fabs(PART[i].u[A_Y]);
			cout<<i<<" "<<Q1<<" "<<Q2<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<endl;
			///*/
		}
		else if(PART[i].N>=2)
		{
			div=0;
			cout<<"周辺粒子数が2<=N<4なのでdiv=0にした。1次近似などにすべき？"<<endl;
		}
		else div=0;	//周辺粒子が少なすぎる場合はゼロにすればよい。どうせ圧力勾配も計算されないから。
	}
	if(d==2 && order==3)//2次元3次式
	{
		if(PART[i].N>=6)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは打ち切り誤差防止
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double U=(PART[j].u[A_X]-PART[i].u[A_X]);
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				//if(dis>1) w=r*r*r*r/(dis*dis*dis*dis);
				if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//表面粒子の速度を無視する
					
				matrix[0]+=X*X*w;			//ΣXjwj
				matrix[1]+=X*Y*w;		//ΣXjYjwj
				matrix[2]+=X*X*X*w;			//ΣXj^3wj
				matrix[3]+=X*X*Y*w;			//ΣXj^2Yjwj
				matrix[4]+=X*Y*Y*w;			//ΣXjYj^2wj
				matrix[5]+=X*X*X*X*w;
				matrix[6]+=X*X*X*Y*w;
				matrix[7]+=X*X*Y*Y*w;
				matrix[8]+=X*Y*Y*Y*w;

				matrix[10]+=Y*Y*w;			
				matrix[11]+=X*X*Y*w;		
				matrix[13]+=Y*Y*Y*w;
				matrix[17]+=Y*Y*Y*Y*w;
	
				matrix[23]+=X*X*X*X*X*w;
				matrix[24]+=X*X*X*X*Y*w;
				matrix[25]+=X*X*X*Y*Y*w;
				matrix[26]+=X*X*Y*Y*Y*w;
				matrix[35]+=X*Y*Y*Y*Y*w;

				matrix[44]+=Y*Y*Y*Y*Y*w;
				matrix[50]+=X*X*X*X*X*X*w;
				matrix[51]+=X*X*X*X*X*Y*w;
				matrix[52]+=X*X*X*X*Y*Y*w;
				matrix[53]+=X*X*X*Y*Y*Y*w;

				matrix[62]+=X*X*Y*Y*Y*Y*w;
				matrix[71]+=X*Y*Y*Y*Y*Y*w;
				matrix[80]+=Y*Y*Y*Y*Y*Y*w;

				
				B1[0]+=U*X*w;//ΣujXjwj
				B1[1]+=U*Y*w;//ΣujYjwj
				B1[2]+=U*X*X*w;//ΣujXj^2wj
				B1[3]+=U*X*Y*w;//ΣujXjYjwj
				B1[4]+=U*Y*Y*w;//ΣujYj^2wj
				B1[5]+=U*X*X*X*w;
				B1[6]+=U*X*X*Y*w;
				B1[7]+=U*X*Y*Y*w;
				B1[8]+=U*Y*Y*Y*w;

				B2[0]+=V*X*w;//ΣvjXjwj
				B2[1]+=V*Y*w;//ΣvjYjwj
				B2[2]+=V*X*X*w;//ΣvjXj^2wj
				B2[3]+=V*X*Y*w;//ΣvjXjYjwj
				B2[4]+=V*Y*Y*w;//ΣvjYj^2wj
				B2[5]+=V*X*X*X*w;
				B2[6]+=V*X*X*Y*w;
				B2[7]+=V*X*Y*Y*w;
				B2[8]+=V*Y*Y*Y*w;
			}
			
			matrix[9]=matrix[1];		//ΣXjYjwj
			matrix[11]=matrix[3]; matrix[19]=matrix[3]; matrix[27]=matrix[3];
			matrix[12]=matrix[4]; matrix[28]=matrix[4]; matrix[36]=matrix[4];
			matrix[14]=matrix[6]; matrix[46]=matrix[6]; matrix[54]=matrix[6]; matrix[21]=matrix[6];	matrix[29]=matrix[6];//X*X*X*Y
			matrix[15]=matrix[7]; matrix[55]=matrix[7]; matrix[63]=matrix[7]; matrix[22]=matrix[7];	matrix[30]=matrix[7]; matrix[38]=matrix[7];//X*X*Y*Y
			matrix[16]=matrix[8]; matrix[64]=matrix[8]; matrix[72]=matrix[8]; matrix[31]=matrix[8];	matrix[39]=matrix[8];//X*Y*Y*Y
			matrix[18]=matrix[2];	//X*X*X
			matrix[20]=matrix[5];	//X*X*X*X
			matrix[47]=matrix[23];	//X*X*X*X*X
			matrix[32]=matrix[24];	matrix[48]=matrix[24]; matrix[56]=matrix[24];	//X*X*X*X*Y		
			matrix[33]=matrix[25];	matrix[41]=matrix[25]; matrix[49]=matrix[25]; matrix[57]=matrix[25]; matrix[63]=matrix[25];//X*X*X*Y*Y
			matrix[34]=matrix[26];  matrix[34]=matrix[42]; matrix[34]=matrix[58]; matrix[34]=matrix[64]; matrix[34]=matrix[74];//X*X*Y*Y*Y
			matrix[43]=matrix[35];  matrix[67]=matrix[35]; matrix[75]=matrix[35];//X*Y*Y*Y*Y
			matrix[37]=matrix[13];
			matrix[40]=matrix[17];	matrix[73]=matrix[17];//Y*Y*Y*Y
			matrix[76]=matrix[44];
			matrix[45]=matrix[5];	//X*X*X*X
			matrix[59]=matrix[51];	//X*X*X*X*X*Y
			matrix[60]=matrix[52];	matrix[68]=matrix[52];	//X*X*X*X*Y*Y	
			matrix[61]=matrix[53];  matrix[69]=matrix[53]; matrix[77]=matrix[53]; 
			matrix[70]=matrix[62];  matrix[78]=matrix[62];
			matrix[79]=matrix[71];

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//値を保存

			gauss(matrix,B1,N);//ガウスの消去法で解く

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//ガウスの消去法で解く

			double a_u,b_u,c_u,d_u,e_u,f_u,g_u,h_u,i_u;
			double a_v,b_v,c_v,d_v,e_v,f_v,g_v,h_v,i_v;
			a_u=B1[0]; a_v=B2[0];//X方向微分
			b_u=B1[1]; b_v=B2[1];//Y方向微分
			c_u=B1[2]; c_v=B2[2];//X方向2階微分
			d_u=B1[3]; d_v=B2[3];//XY方向微分
			e_u=B1[4]; e_v=B2[4];//Y方向2階微分
			f_u=B1[5]; f_v=B2[5];
			g_u=B1[6]; g_v=B2[6];
			h_u=B1[7]; h_v=B2[7];
			i_u=B1[8]; i_v=B2[8];
		
			div=(B1[0]+B2[1]);
			if(div+div!=2*div) cout<<"erroe"<<endl;
		}
		else if(PART[i].N>=2)
		{
			div=0;
			cout<<"周辺粒子数が2<=N<6なのでdiv=0にした。1次近似などにすべき？"<<endl;
		}
		else div=0;	//周辺粒子が少なすぎる場合はゼロにすればよい。どうせ圧力勾配も計算されないから。
	}
	else if(d==3 && order==1)//3次元1次式
	{
		if(PART[i].N>5)
		{
			//div=cacl_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//最後の引数は未知数
			div=cacl_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//最後の引数は未知数
		}
		else 
		{
			div=0;
			//cout<<"警告 近隣粒子数が5以下("<<PART[i].N<<")です"<<endl;
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
		
		int nei=0;
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<=le) nei++;
		}

		//if(PART[i].N>8)
		if(nei>=6)			//le以下の粒子が６以上あれば２次近似
		{
			//div=cacl_WLSM_divu_D3_order2(CON,PART,matrix,B1,B2,B3,i,9);
			div=cacl_WLSM_divu_D3_order2_2(CON,PART,matrix,B1,B2,B3,i,10);
		}
		else if(PART[i].N>5)
		{
			//div=cacl_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//この場合、最後の引数はN=3を渡すことに注意
			div=cacl_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//最後の引数は未知数
		}
		else 
		{
			div=0;
			//cout<<"警告 近隣粒子数が9以下("<<PART[i].N<<")です"<<endl;
		}
	}

	delete [] matrix;
	delete [] matrix2;
	delete [] B1;
	delete [] B2;
	delete [] B3;

    return div;
}

///速度発散計算関数(WLSM法) 近似曲面が自身の関数値を通るという仮定をしない場合
double divergence3(mpsconfig *CON,vector<mpsparticle> &PART,int i,int surface_sw)
{
	//surface_sw: ONなら表面を無視する、OFFなら表面も考慮にいれた発散を計算する。ON,OFFの定義を混乱しないように
    double div=0;//発散の値
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N0=0;							//係数行列の元
	int order=CON->get_divU_order();				//近似曲面のｵｰﾀﾞｰ。 1=線形 2=二次
    
	//係数行列の大きさの決定
	if(d==2)
	{
		if(order==1) N0=3;
		else if(order==2) N0=6;
		else if(order==3) N0=10;
	}
	else if(d==3)
	{
		if(order==1) N0=4;
		else if(order==2) N0=10;
		else if(order==3) N0=20;
	}
	////////////////////////////////

	double *matrix=new double [N0*N0];	//N×Nの係数行列
	double *matrix2=new double [N0*N0];	//matrixのバックアップ
	double **MAT= new double*[N0];		//N×Nの係数行列(配列は2次元)
	for(int n=0;n<N0;n++) MAT[n]=new double[N0];
	double *base=new double[N0];			//基底ベクトル格納
	double *B1=new double [N0];			//Nの解行列
	double *B2=new double [N0];			//Nの解行列
	double *B3=new double [N0];			//Nの解行列

	for(int n=0;n<N0*N0;n++) matrix[n]=0;	//初期化
	for(int n=0;n<N0;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}
	for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//初期化

	int N=N0;
	int jnb=PART[i].N;	//周辺粒子数

	if(surface_sw==ON)	//表面粒子を無視するなら(表面粒子の速度は速度発散ゼロを満たしていないから)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			if(PART[j].type==FLUID && PART[j].surface==ON) jnb--; 
		}
	}

	if(d==2 )				//二次元
	{
		if(jnb<15) {order=2; N=6;}		//周辺粒子数が3次近似するのに十分でないなら、2次近似する(保険もかけて少し大きめにとってある)
		if(jnb<6) {order=1; N=3;}		//周辺粒子数が２時近似するのに十分でないなら、１次近似する
		if(CON->get_divU_order()==1) {order=1; N=3;}
		if(jnb>=2)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// leで割るのは正規化のため
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
				double U=(PART[j].u[A_X])/PART[i].L;
				double V=(PART[j].u[A_Y])/PART[i].L;
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				if(dis>1) w=1/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//表面粒子の速度を無視する
				
				if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//基底ベクトル
				else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=1;}
				else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y; base[9]=1;}

				//行列作成
				for(int n=0;n<N;n++)
				{
					for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
					B1[n]+=base[n]*w*U;
					B2[n]+=base[n]*w*V;
				}
			}
			
			MAT[N-1][N-1]+=1;		//一番右下の配列に自分自身の寄与を加える
			B1[N-1]+=PART[i].u[A_X]/PART[i].L;	//粒子iの寄与。
			B2[N-1]+=PART[i].u[A_Y]/PART[i].L;	//粒子iの寄与。
			//これ以外はX=0になるので加算する必要はない*/

			/*double *eigen=new double[N];
			calc_eigen_by_jacobi(MAT, N,eigen);	//固有値を求める
			double min=1000; double max=0;
			for(int n=0;n<N;n++) 
			{
				if(eigen[n]<min) min=eigen[n];
				if(eigen[n]>max) max=eigen[n];
			}
			delete [] eigen;*/

			/*double maxD=MAT[0][0]; int maxid=0;
			double minD=MAT[0][0];	int minid=0;
			for(int n=0;n<N;n++)
			{
				if(maxD<MAT[n][n]) {maxD=MAT[n][n]; maxid=n;}
				if(minD>MAT[n][n]) {minD=MAT[n][n]; minid=n;}
			}
			//maxD=0; minD=0;
			for(int n=0;n<N;n++) if(n!=maxid) {maxD+=MAT[maxid][n]; minD-=MAT[minid][n];}
			
			//cout<<maxD/minD<<endl;/*/
					
			for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//値をmatrixに転送

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//値を保存

			gauss(matrix,B1,N);//ガウスの消去法で解く

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//ガウスの消去法で解く

			div=(B1[0]+B2[1]);

			//cout<<div<<" "<<N<<" "<<jnb<<" "<<max/min<<endl;

			if(div+div!=2*div) cout<<"速度の発散が非実数 i="<<i<<endl;
			else
			{
				/*double Q1=0; double Q2=0;		//Uの誤差
				int num=1;
				B1[N-1]*=PART[i].L; 
				B2[N-1]*=PART[i].L; 
				for(int k=0;k<PART[i].N;k++)
				{
					int j=PART[i].NEI[k];
					double X=(PART[j].r[A_X]-PART[i].r[A_X]);
					double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
					double U=(PART[j].u[A_X]);
					double V=(PART[j].u[A_Y]);
					double dis=sqrt(X*X+Y*Y);
					
					double w=1;
					if(dis>1) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
					if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//表面粒子の速度を無視する
					if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//基底ベクトル
					else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=1;}
					else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y; base[9]=1;}
					
					double uu=0; double vv=0;			//近時曲面上の速度
					for(int n=0;n<N;n++) uu+=base[n]*B1[n];
					for(int n=0;n<N;n++) vv+=base[n]*B2[n];
					Q1+=(uu-U)*(uu-U)*w;
					Q2+=(vv-V)*(vv-V)*w;
					if(w!=0) num++;
				}
				Q1+=(B1[N-1]-PART[i].u[A_X])*(B1[N-1]-PART[i].u[A_X]);
				Q2+=(B2[N-1]-PART[i].u[A_Y])*(B2[N-1]-PART[i].u[A_Y]);
				Q1/=num*PART[i].u[A_X]; Q2/=num*PART[i].u[A_Y];

				if(Q1>0.01 || Q2>0.01) cout<<i<<" "<<Q1<<" "<<Q2<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<endl;*/
			}
		}
		else div=0;	//周辺粒子が少なすぎる場合はゼロにすればよい。どうせ圧力勾配も計算されないから。
	}
	else if(d==3)
	{
		///係数行列は1次近似なら
		///   ΣΔx2    ΣΔxΔy  ΣΔxΔz ΣΔx a = ΣΔxfj  
		///  ΣΔxΔy    ΣΔy2   ΣΔyΔz ΣΔy b = ΣΔyfj 
		///  ΣΔxΔz   ΣΔyΔz  ΣΔz2  ΣΔz  c = ΣΔzfj 
		///  ΣΔx      ΣΔy     ΣΔz     Σ1  d = Σfj

		//2次近似なら
		//P=aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+kΔzΔx+lとおくと、
		///係数行列は
		///   ΣΔx2      ΣΔxΔy    ΣΔxΔz    ΣΔx3       ΣΔxΔy2    ΣΔxΔz2    ΣΔx2Δy     ΣΔxΔyΔz  ΣΔx2Δz    ΣΔx    a = ΣΔxf  
		///   ΣΔxΔy    ΣΔy2      ΣΔyΔz    ΣΔx2Δy    ΣΔy3       ΣΔyΔz2    ΣΔxΔy2     ΣΔy2Δz    ΣΔxΔyΔz  ΣΔy    b = ΣΔyf
		///   ΣΔxΔz    ΣΔyΔz    ΣΔz2      ΣΔx2Δz    ΣΔy2Δz    ΣΔz3       ΣΔxΔyΔz   ΣΔyΔz2    ΣΔxΔz2    ΣΔz    c = ΣΔzf
		///   ΣΔx3      ΣΔx2Δy   ΣΔx2Δz   ΣΔx4       ΣΔx2Δy2   ΣΔx2Δz2   ΣΔx3Δy     ΣΔx2ΔyΔz ΣΔx3Δz    ΣΔx2   d = ΣΔx2f
		///   ΣΔxΔy2   ΣΔy3      ΣΔy2Δz   ΣΔx2Δy2   ΣΔy4       ΣΔy2Δz2   ΣΔxΔy3     ΣΔy3Δz    ΣΔxΔy2Δz ΣΔy2   e = ΣΔy2f
		///   ΣΔxΔz2   ΣΔyΔz2   ΣΔz3      ΣΔx2Δz2   ΣΔy2Δz2   ΣΔz4       ΣΔxΔyΔz2  ΣΔyΔz3    ΣΔxΔz3	 ΣΔz2   f = ΣΔz2f
		///   ΣΔx2Δy   ΣΔxΔy2   ΣΔxΔyΔz ΣΔx3Δy    ΣΔxΔy3    ΣΔxΔyΔz2 ΣΔx2Δy2    ΣΔxΔy2Δz ΣΔx2ΔyΔz ΣΔxΔy g = ΣΔxΔyf
		///   ΣΔxΔyΔz ΣΔy2Δz   ΣΔyΔz2   ΣΔx2ΔyΔz ΣΔy3Δz    ΣΔyΔz3    ΣΔxΔy2Δz  ΣΔy2Δz2   ΣΔxΔyΔz2 ΣΔyΔz h = ΣΔyΔzf
		///   ΣΔx2Δz   ΣΔxΔyΔz ΣΔxΔz2   ΣΔx3Δz    ΣΔxΔy2Δz  ΣΔxΔz3   ΣΔx2Δy Δz ΣΔxΔyΔz2 ΣΔx2Δz2   ΣΔzΔx k = ΣΔxΔzf
		///   ΣΔx       ΣΔy ΣΔz ΣΔx2      ΣΔy2       ΣΔz        ΣΔxΔy     ΣΔyΔz      ΣΔxΔz     ΣΔzΔx     Σ1      l = Σfj

		//3次近似なら
		//P=aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+kΔzΔx+lΔx3+mΔy3+nΔz3+oΔx2Δy+pΔxΔy2+qΔy2Δz+rΔyΔz2+sΔx2Δz+tΔz2Δx+uΔxΔyΔz+vとおく
		
		if(order==3) 
		{
			if(jnb<20) {order=2; N=10;}//粒子数が少ないときは2次近似にする
		}
		if(order==2) 
		{
			if(jnb<10) {order=1; N=4;}//それでも粒子数が少ないときは１次近似にする
		}
		
		
		if(jnb>=3)	//周辺粒子が２つ以下だったら計算しない
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは打ち切り誤差防止
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double U=(PART[j].u[A_X]);
				double V=(PART[j].u[A_Y]);
				double W=(PART[j].u[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
					
				double w=1;
				if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//表面粒子の速度を無視する

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
			
			MAT[N-1][N-1]+=1;		//一番右下の配列に自分自身の寄与を加える
			B1[N-1]+=PART[i].u[A_X];	//粒子iの寄与。
			B2[N-1]+=PART[i].u[A_Y];	//粒子iの寄与。
			B3[N-1]+=PART[i].u[A_Z];	//粒子iの寄与。
			//これ以外はX=0になるので加算する必要はない
					
			for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//値をmatrixに転送

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//値を保存

			gauss(matrix,B1,N);//ガウスの消去法で解く

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//ガウスの消去法で解く

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B3,N);//ガウスの消去法で解く

			
			double dudx=B1[0];	//1次近似でも2次近似でも、未知数の順番的にこうなるようにしてある
			double dvdy=B2[1];
			double dwdz=B3[2];
			//cout<<dudx<<" "<<dvdy<<" "<<dwdz<<endl;
				
			div=(dudx+dvdy+dwdz);
			
		}
		else div=0;
	}

	delete [] matrix;
	delete [] matrix2;
	delete [] B1;
	delete [] B2;
	delete [] B3;
	for(int n=0;n<N0;n++) delete [] MAT[n];
    delete [] MAT;
	delete [] base;

    return div;
}

///速度発散計算関数(入部らの方法)　テンソル積と重み平均を用いる
double divergence4(mpsconfig *CON,vector<mpsparticle> &PART,int i)
{
	//2次元しかできてない
    double div=0;//発散の値
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N=0;					//係数行列の元
	int order=1;				//近似曲面のｵｰﾀﾞｰ。 1=線形 2=二次
    
	//係数行列の大きさの決定
	if(order==1) N=d;
	//else //まだできてない
	////////////////////////////////

	double *matrix=new double [N*N];	//N×Nの係数行列
	double *B1=new double [N];			//Nの解行列
	double *B2=new double [N];			//Nの解行列
	double *B3=new double [N];			//Nの解行列

	for(int n=0;n<N*N;n++) matrix[n]=0;	//初期化
	for(int n=0;n<N;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}

	if(d==2 && order==1)				//二次元
	{
		if(PART[i].N>=2)
		{
			double W=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは正規化のため
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double U=(PART[j].u[A_X]-PART[i].u[A_X]);
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
				double dis=sqrt(X*X+Y*Y);
				X/=dis;	Y/=dis;				//単位ベクトルにする
				double w=1;
				//if(dis>1) w=1/(dis*dis*dis*dis);
				if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
					
				matrix[0]+=X*X*w;
				matrix[1]+=X*Y*w;
				matrix[3]+=Y*Y*w;
				
				B1[A_X]+=U*w*X/dis;
				B1[A_Y]+=U*w*Y/dis;

				B2[A_X]+=V*w*X/dis;
				B2[A_Y]+=V*w*Y/dis;

				W+=w;
			}
			
			matrix[2]=matrix[1];

			for(int n=0;n<N*N;n++) matrix[n]/=W;
			for(int n=0;n<N;n++)
			{
				B1[n]/=W;
				B2[n]/=W;
			}

			//matrixの逆行列を求める 求まった逆行列はmatrixの中を上書きして格納される
			calc_inverse_matrix(CON,PART, N, matrix);

			double dudx=matrix[0]*B1[0]+matrix[1]*B1[1];
			double dudy=matrix[2]*B1[0]+matrix[3]*B1[1];
			double dvdx=matrix[0]*B2[0]+matrix[1]*B2[1];
			double dvdy=matrix[2]*B2[0]+matrix[3]*B2[1];
			
			div=(dudx+dvdy);
			if(div+div!=2*div) cout<<"速度の発散が非実数 i="<<i<<endl;

		}
		else div=0;	//周辺粒子が少なすぎる場合はゼロにすればよい。どうせ圧力勾配も計算されないから。
	}
	else cout<<"まだできてない"<<endl;
	

	delete [] matrix;
	delete [] B1;
	delete [] B2;
	delete [] B3;

    return div;
}

//5次の連立方程式の解1,2を返す関数
void return_X_for5N(double *matrix,int N,double *B1,double *B2,double *dudx,double *dudy)
{
	double a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55;
	double b1,b2,b3,b4,b5;
	double c1,c2,c3,c4,c5;

	a11=matrix[0];a12=matrix[1];a13=matrix[2];a14=matrix[3];a15=matrix[4];
	a21=matrix[5];a22=matrix[6];a23=matrix[7];a24=matrix[8];a25=matrix[9];
	a31=matrix[10];a32=matrix[11];a33=matrix[12];a34=matrix[13];a35=matrix[14];
	a41=matrix[15];a42=matrix[16];a43=matrix[17];a44=matrix[18];a45=matrix[19];
	a51=matrix[20];a52=matrix[21];a53=matrix[22];a54=matrix[23];a55=matrix[24];

	b1=B1[0];b2=B1[1];b3=B1[2];b4=B1[3];b5=B1[4];
	c1=B2[0];c2=B2[1];c3=B2[2];c4=B2[3];c5=B2[4];
	
	double determinant=(a11*a22*a33*a44*a55-a11*a22*a33*a45*a54-a11*a22*a34*a43*a55+a11*a22*a34*a45*a53+a11*a22*a35*a43*a54-a11*a22*a35*a44*a53-a11*a23*a32*a44*a55+a11*a23*a32*a45*a54+a11*a23*a34*a42*a55-a11*a23*a34*a45*a52-a11*a23*a35*a42*a54+a11*a23*a35*a44*a52+a11*a24*a32*a43*a55-a11*a24*a32*a45*a53-a11*a24*a33*a42*a55+a11*a24*a33*a45*a52+a11*a24*a35*a42*a53-a11*a24*a35*a43*a52-a11*a25*a32*a43*a54+a11*a25*a32*a44*a53+a11*a25*a33*a42*a54-a11*a25*a33*a44*a52-a11*a25*a34*a42*a53+a11*a25*a34*a43*a52-a12*a21*a33*a44*a55+a12*a21*a33*a45*a54+a12*a21*a34*a43*a55-a12*a21*a34*a45*a53-a12*a21*a35*a43*a54+a12*a21*a35*a44*a53+a12*a23*a31*a44*a55-a12*a23*a31*a45*a54-a12*a23*a34*a41*a55+a12*a23*a34*a45*a51+a12*a23*a35*a41*a54-a12*a23*a35*a44*a51-a12*a24*a31*a43*a55+a12*a24*a31*a45*a53+a12*a24*a33*a41*a55-a12*a24*a33*a45*a51-a12*a24*a35*a41*a53+a12*a24*a35*a43*a51+a12*a25*a31*a43*a54-a12*a25*a31*a44*a53-a12*a25*a33*a41*a54+a12*a25*a33*a44*a51+a12*a25*a34*a41*a53-a12*a25*a34*a43*a51+a13*a21*a32*a44*a55-a13*a21*a32*a45*a54-a13*a21*a34*a42*a55+a13*a21*a34*a45*a52+a13*a21*a35*a42*a54-a13*a21*a35*a44*a52-a13*a22*a31*a44*a55+a13*a22*a31*a45*a54+a13*a22*a34*a41*a55-a13*a22*a34*a45*a51-a13*a22*a35*a41*a54+a13*a22*a35*a44*a51+a13*a24*a31*a42*a55-a13*a24*a31*a45*a52-a13*a24*a32*a41*a55+a13*a24*a32*a45*a51+a13*a24*a35*a41*a52-a13*a24*a35*a42*a51-a13*a25*a31*a42*a54+a13*a25*a31*a44*a52+a13*a25*a32*a41*a54-a13*a25*a32*a44*a51-a13*a25*a34*a41*a52+a13*a25*a34*a42*a51-a14*a21*a32*a43*a55+a14*a21*a32*a45*a53+a14*a21*a33*a42*a55-a14*a21*a33*a45*a52-a14*a21*a35*a42*a53+a14*a21*a35*a43*a52+a14*a22*a31*a43*a55-a14*a22*a31*a45*a53-a14*a22*a33*a41*a55+a14*a22*a33*a45*a51+a14*a22*a35*a41*a53-a14*a22*a35*a43*a51-a14*a23*a31*a42*a55+a14*a23*a31*a45*a52+a14*a23*a32*a41*a55-a14*a23*a32*a45*a51-a14*a23*a35*a41*a52+a14*a23*a35*a42*a51+a14*a25*a31*a42*a53-a14*a25*a31*a43*a52-a14*a25*a32*a41*a53+a14*a25*a32*a43*a51+a14*a25*a33*a41*a52-a14*a25*a33*a42*a51+a15*a21*a32*a43*a54-a15*a21*a32*a44*a53-a15*a21*a33*a42*a54+a15*a21*a33*a44*a52+a15*a21*a34*a42*a53-a15*a21*a34*a43*a52-a15*a22*a31*a43*a54+a15*a22*a31*a44*a53+a15*a22*a33*a41*a54-a15*a22*a33*a44*a51-a15*a22*a34*a41*a53+a15*a22*a34*a43*a51+a15*a23*a31*a42*a54-a15*a23*a31*a44*a52-a15*a23*a32*a41*a54+a15*a23*a32*a44*a51+a15*a23*a34*a41*a52-a15*a23*a34*a42*a51-a15*a24*a31*a42*a53+a15*a24*a31*a43*a52+a15*a24*a32*a41*a53-a15*a24*a32*a43*a51-a15*a24*a33*a41*a52+a15*a24*a33*a42*a51);
	
	*dudx=(b1*a22*a33*a44*a55-b1*a22*a33*a45*a54-b1*a22*a34*a43*a55+b1*a22*a34*a45*a53+b1*a22*a35*a43*a54-b1*a22*a35*a44*a53-b1*a23*a32*a44*a55+b1*a23*a32*a45*a54+b1*a23*a34*a42*a55-b1*a23*a34*a45*a52-b1*a23*a35*a42*a54+b1*a23*a35*a44*a52+b1*a24*a32*a43*a55-b1*a24*a32*a45*a53-b1*a24*a33*a42*a55+b1*a24*a33*a45*a52+b1*a24*a35*a42*a53-b1*a24*a35*a43*a52-b1*a25*a32*a43*a54+b1*a25*a32*a44*a53+b1*a25*a33*a42*a54-b1*a25*a33*a44*a52-b1*a25*a34*a42*a53+b1*a25*a34*a43*a52-a12*b2*a33*a44*a55+a12*b2*a33*a45*a54+a12*b2*a34*a43*a55-a12*b2*a34*a45*a53-a12*b2*a35*a43*a54+a12*b2*a35*a44*a53+a12*a23*b3*a44*a55-a12*a23*b3*a45*a54-a12*a23*a34*b4*a55+a12*a23*a34*a45*b5+a12*a23*a35*b4*a54-a12*a23*a35*a44*b5-a12*a24*b3*a43*a55+a12*a24*b3*a45*a53+a12*a24*a33*b4*a55-a12*a24*a33*a45*b5-a12*a24*a35*b4*a53+a12*a24*a35*a43*b5+a12*a25*b3*a43*a54-a12*a25*b3*a44*a53-a12*a25*a33*b4*a54+a12*a25*a33*a44*b5+a12*a25*a34*b4*a53-a12*a25*a34*a43*b5+a13*b2*a32*a44*a55-a13*b2*a32*a45*a54-a13*b2*a34*a42*a55+a13*b2*a34*a45*a52+a13*b2*a35*a42*a54-a13*b2*a35*a44*a52-a13*a22*b3*a44*a55+a13*a22*b3*a45*a54+a13*a22*a34*b4*a55-a13*a22*a34*a45*b5-a13*a22*a35*b4*a54+a13*a22*a35*a44*b5+a13*a24*b3*a42*a55-a13*a24*b3*a45*a52-a13*a24*a32*b4*a55+a13*a24*a32*a45*b5+a13*a24*a35*b4*a52-a13*a24*a35*a42*b5-a13*a25*b3*a42*a54+a13*a25*b3*a44*a52+a13*a25*a32*b4*a54-a13*a25*a32*a44*b5-a13*a25*a34*b4*a52+a13*a25*a34*a42*b5-a14*b2*a32*a43*a55+a14*b2*a32*a45*a53+a14*b2*a33*a42*a55-a14*b2*a33*a45*a52-a14*b2*a35*a42*a53+a14*b2*a35*a43*a52+a14*a22*b3*a43*a55-a14*a22*b3*a45*a53-a14*a22*a33*b4*a55+a14*a22*a33*a45*b5+a14*a22*a35*b4*a53-a14*a22*a35*a43*b5-a14*a23*b3*a42*a55+a14*a23*b3*a45*a52+a14*a23*a32*b4*a55-a14*a23*a32*a45*b5-a14*a23*a35*b4*a52+a14*a23*a35*a42*b5+a14*a25*b3*a42*a53-a14*a25*b3*a43*a52-a14*a25*a32*b4*a53+a14*a25*a32*a43*b5+a14*a25*a33*b4*a52-a14*a25*a33*a42*b5+a15*b2*a32*a43*a54-a15*b2*a32*a44*a53-a15*b2*a33*a42*a54+a15*b2*a33*a44*a52+a15*b2*a34*a42*a53-a15*b2*a34*a43*a52-a15*a22*b3*a43*a54+a15*a22*b3*a44*a53+a15*a22*a33*b4*a54-a15*a22*a33*a44*b5-a15*a22*a34*b4*a53+a15*a22*a34*a43*b5+a15*a23*b3*a42*a54-a15*a23*b3*a44*a52-a15*a23*a32*b4*a54+a15*a23*a32*a44*b5+a15*a23*a34*b4*a52-a15*a23*a34*a42*b5-a15*a24*b3*a42*a53+a15*a24*b3*a43*a52+a15*a24*a32*b4*a53-a15*a24*a32*a43*b5-a15*a24*a33*b4*a52+a15*a24*a33*a42*b5)/determinant;
	*dudy=(a11*c2*a33*a44*a55-a11*c2*a33*a45*a54-a11*c2*a34*a43*a55+a11*c2*a34*a45*a53+a11*c2*a35*a43*a54-a11*c2*a35*a44*a53-a11*a23*c3*a44*a55+a11*a23*c3*a45*a54+a11*a23*a34*c4*a55-a11*a23*a34*a45*c5-a11*a23*a35*c4*a54+a11*a23*a35*a44*c5+a11*a24*c3*a43*a55-a11*a24*c3*a45*a53-a11*a24*a33*c4*a55+a11*a24*a33*a45*c5+a11*a24*a35*c4*a53-a11*a24*a35*a43*c5-a11*a25*c3*a43*a54+a11*a25*c3*a44*a53+a11*a25*a33*c4*a54-a11*a25*a33*a44*c5-a11*a25*a34*c4*a53+a11*a25*a34*a43*c5-c1*a21*a33*a44*a55+c1*a21*a33*a45*a54+c1*a21*a34*a43*a55-c1*a21*a34*a45*a53-c1*a21*a35*a43*a54+c1*a21*a35*a44*a53+c1*a23*a31*a44*a55-c1*a23*a31*a45*a54-c1*a23*a34*a41*a55+c1*a23*a34*a45*a51+c1*a23*a35*a41*a54-c1*a23*a35*a44*a51-c1*a24*a31*a43*a55+c1*a24*a31*a45*a53+c1*a24*a33*a41*a55-c1*a24*a33*a45*a51-c1*a24*a35*a41*a53+c1*a24*a35*a43*a51+c1*a25*a31*a43*a54-c1*a25*a31*a44*a53-c1*a25*a33*a41*a54+c1*a25*a33*a44*a51+c1*a25*a34*a41*a53-c1*a25*a34*a43*a51+a13*a21*c3*a44*a55-a13*a21*c3*a45*a54-a13*a21*a34*c4*a55+a13*a21*a34*a45*c5+a13*a21*a35*c4*a54-a13*a21*a35*a44*c5-a13*c2*a31*a44*a55+a13*c2*a31*a45*a54+a13*c2*a34*a41*a55-a13*c2*a34*a45*a51-a13*c2*a35*a41*a54+a13*c2*a35*a44*a51+a13*a24*a31*c4*a55-a13*a24*a31*a45*c5-a13*a24*c3*a41*a55+a13*a24*c3*a45*a51+a13*a24*a35*a41*c5-a13*a24*a35*c4*a51-a13*a25*a31*c4*a54+a13*a25*a31*a44*c5+a13*a25*c3*a41*a54-a13*a25*c3*a44*a51-a13*a25*a34*a41*c5+a13*a25*a34*c4*a51-a14*a21*c3*a43*a55+a14*a21*c3*a45*a53+a14*a21*a33*c4*a55-a14*a21*a33*a45*c5-a14*a21*a35*c4*a53+a14*a21*a35*a43*c5+a14*c2*a31*a43*a55-a14*c2*a31*a45*a53-a14*c2*a33*a41*a55+a14*c2*a33*a45*a51+a14*c2*a35*a41*a53-a14*c2*a35*a43*a51-a14*a23*a31*c4*a55+a14*a23*a31*a45*c5+a14*a23*c3*a41*a55-a14*a23*c3*a45*a51-a14*a23*a35*a41*c5+a14*a23*a35*c4*a51+a14*a25*a31*c4*a53-a14*a25*a31*a43*c5-a14*a25*c3*a41*a53+a14*a25*c3*a43*a51+a14*a25*a33*a41*c5-a14*a25*a33*c4*a51+a15*a21*c3*a43*a54-a15*a21*c3*a44*a53-a15*a21*a33*c4*a54+a15*a21*a33*a44*c5+a15*a21*a34*c4*a53-a15*a21*a34*a43*c5-a15*c2*a31*a43*a54+a15*c2*a31*a44*a53+a15*c2*a33*a41*a54-a15*c2*a33*a44*a51-a15*c2*a34*a41*a53+a15*c2*a34*a43*a51+a15*a23*a31*c4*a54-a15*a23*a31*a44*c5-a15*a23*c3*a41*a54+a15*a23*c3*a44*a51+a15*a23*a34*a41*c5-a15*a23*a34*c4*a51-a15*a24*a31*c4*a53+a15*a24*a31*a43*c5+a15*a24*c3*a41*a53-a15*a24*c3*a43*a51-a15*a24*a33*a41*c5+a15*a24*a33*c4*a51)/determinant;
	
}

//divergence2における、3次元1次近似を行う関数
double cacl_WLSM_divu_D3_order1(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///係数行列は
	///   ΣΔx2    ΣΔxΔy  ΣΔxΔz  a = ΣΔxΔf  
	///  ΣΔxΔy    ΣΔy2   ΣΔyΔz  b = ΣΔyΔf 
	///  ΣΔxΔz   ΣΔyΔz   ΣΔz2   c = ΣΔzΔf 

	double le=CON->get_distancebp();
	double matrix_val[9];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは打ち切り誤差防止
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
					
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
		weight[k]=w;
					
		matrix[0]+=X*X*w;			//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
					
		matrix[4]+=Y*Y*w;			//ΣYj^2wj
		matrix[5]+=Y*Z*w;		//ΣYjZjwj

		matrix[8]+=Z*Z*w;			//ΣZj^2wj
			
		B1[0]+=U*X*w;//ΣfjXjwj
		B1[1]+=U*Y*w;//ΣfjYjwj
		B1[2]+=U*Z*w;//ΣfjZjwj

		B2[0]+=V*X*w;//ΣfjXjwj
		B2[1]+=V*Y*w;//ΣfjYjwj
		B2[2]+=V*Z*w;//ΣfjZjwj

		B3[0]+=W*X*w;//ΣfjXjwj
		B3[1]+=W*Y*w;//ΣfjYjwj
		B3[2]+=W*Z*w;//ΣfjZjwj
	}
			
	matrix[3]=matrix[1];		//ΣXjYjwj
	matrix[6]=matrix[2];		//ΣXjZjwj
	matrix[7]=matrix[5];		//ΣYjZjwj

	for(int L=0;L<9;L++) matrix_val[L]=matrix[L];//行列の値を保存

	/*double dudx=0;//こっちのほうが若干早い。けど誤差評価したいならガウス
	double dvdy=0;
	double dwdz=0;
	double determinant=(matrix[0]*matrix[4]*matrix[8]-matrix[0]*matrix[5]*matrix[7]-matrix[1]*matrix[3]*matrix[8]+matrix[1]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*matrix[7]-matrix[2]*matrix[4]*matrix[6]);//行列式
			
	dudx=(B1[0]*matrix[4]*matrix[8]-B1[0]*matrix[5]*matrix[7]-matrix[1]*B1[1]*matrix[8]+matrix[1]*matrix[5]*B1[2]+matrix[2]*B1[1]*matrix[7]-matrix[2]*matrix[4]*B1[2])/determinant;
	dvdy=(matrix[0]*B2[1]*matrix[8]-matrix[0]*matrix[5]*B2[2]-B2[0]*matrix[3]*matrix[8]+B2[0]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*B2[2]-matrix[2]*B2[1]*matrix[6])/determinant;
	dwdz=(matrix[0]*matrix[4]*B3[2]-matrix[0]*B3[1]*matrix[7]-matrix[1]*matrix[3]*B3[2]+matrix[1]*B3[1]*matrix[6]+B3[0]*matrix[3]*matrix[7]-B3[0]*matrix[4]*matrix[6])/determinant;
	*/	

	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<3;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else B1[0]=0; //flagがOFFならどのみちdudxはゼロ。
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else B2[1]=0;	//flagがOFFならどのみちdvdyはゼロ。
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else B3[2]=0;	//flagがOFFならどのみちdwdzはゼロ。

	//計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の誤差
	double W=0;//重みの総和
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	/////*/

	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];

	double div=(dudx+dvdy+dwdz);

	delete [] weight;

	return div;
}

//divergence2における、3次元1次近似を行う関数ver.2
double cacl_WLSM_divu_D3_order1_2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///係数行列は
	///   ΣΔx2    ΣΔxΔy  ΣΔxΔz ΣΔx a = ΣΔxfj  
	///  ΣΔxΔy    ΣΔy2   ΣΔyΔz ΣΔy b = ΣΔyfj 
	///  ΣΔxΔz   ΣΔyΔz  ΣΔz2  ΣΔz  c = ΣΔzfj 
	///  ΣΔx      ΣΔy     ΣΔz     Σ1  d = Σfj

	double le=CON->get_distancebp();
	double matrix_val[16];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは打ち切り誤差防止
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		//double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		//double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		//double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double U=(PART[j].u[A_X]);
		double V=(PART[j].u[A_Y]);
		double W=(PART[j].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
					
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
		weight[k]=w;
					
		matrix[0]+=X*X*w;			//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*w;
					
		matrix[5]+=Y*Y*w;			//ΣYj^2wj
		matrix[6]+=Y*Z*w;		//ΣYjZjwj
		matrix[7]+=Y*w;

		matrix[10]+=Z*Z*w;			//ΣZj^2wj
		matrix[11]+=Z*w;

		matrix[15]+=w;
			
		B1[0]+=U*X*w;//ΣfjXjwj
		B1[1]+=U*Y*w;//ΣfjYjwj
		B1[2]+=U*Z*w;//ΣfjZjwj
		B1[3]+=U*w;//Σfjwj

		B2[0]+=V*X*w;//ΣfjXjwj
		B2[1]+=V*Y*w;//ΣfjYjwj
		B2[2]+=V*Z*w;//ΣfjZjwj
		B2[3]+=V*w;//Σfjwj

		B3[0]+=W*X*w;//ΣfjXjwj
		B3[1]+=W*Y*w;//ΣfjYjwj
		B3[2]+=W*Z*w;//ΣfjZjwj
		B3[3]+=W*w;//Σfjwj
	}
			
	matrix[4]=matrix[1];	
	matrix[8]=matrix[2];		
	matrix[9]=matrix[6];
	matrix[12]=matrix[3];
	matrix[13]=matrix[7];
	matrix[14]=matrix[11];

	matrix[15]+=1;//自分自身
	B1[3]+=PART[i].u[A_X];
	B2[3]+=PART[i].u[A_Y];
	B3[3]+=PART[i].u[A_Z];

	for(int L=0;L<16;L++) matrix_val[L]=matrix[L];//行列の値を保存

	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<4;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B2[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B3[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else for(int k=0;k<4;k++) B1[k]=0; //flagがOFFならどのみちdudxはゼロ。
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else for(int k=0;k<4;k++) B2[k]=0;	//flagがOFFならどのみちdvdyはゼロ。
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else for(int k=0;k<4;k++) B3[k]=0;	//flagがOFFならどのみちdwdzはゼロ。

	//計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の誤差
	double err[3];
	double W=1;//重みの総和
	err[A_X]=B1[3]-PART[i].u[A_X];//自身の誤差
	err[A_Y]=B2[3]-PART[i].u[A_Y];
	err[A_Z]=B3[3]-PART[i].u[A_Z];
	Q[A_X]+=err[A_X]*err[A_X];
	Q[A_Y]+=err[A_Y]*err[A_Y];
	Q[A_Z]+=err[A_Z]*err[A_Z];
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]);
		double V=(PART[j].u[A_Y]);
		double W=(PART[j].u[A_Z]);
		double w=weight[k];
		W+=w;
		
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]-U;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]-V;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]-W;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	/////*/

	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];

	double div=(dudx+dvdy+dwdz);

	delete [] weight;

	return div;
}

//divergence2における、3次元2次近似を行う関数
double cacl_WLSM_divu_D3_order2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	double le=CON->get_distancebp();
	double matrix_val[81];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
						
		double w=1;
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


		B1[0]+=U*X*w;		//a
		B1[1]+=U*Y*w;		//b
		B1[2]+=U*Z*w;		//c
		B1[3]+=U*X*X*w;		//d
		B1[4]+=U*Y*Y*w;		//e
		B1[5]+=U*Z*Z*w;		//f
		B1[6]+=U*X*Y*w;		//g
		B1[7]+=U*Y*Z*w;		//h
		B1[8]+=U*X*Z*w;		//i

		B2[0]+=V*X*w;		//a
		B2[1]+=V*Y*w;		//b
		B2[2]+=V*Z*w;		//c
		B2[3]+=V*X*X*w;		//d
		B2[4]+=V*Y*Y*w;		//e
		B2[5]+=V*Z*Z*w;		//f
		B2[6]+=V*X*Y*w;		//g
		B2[7]+=V*Y*Z*w;		//h
		B2[8]+=V*X*Z*w;		//i

		B3[0]+=W*X*w;		//a
		B3[1]+=W*Y*w;		//b
		B3[2]+=W*Z*w;		//c
		B3[3]+=W*X*X*w;		//d
		B3[4]+=W*Y*Y*w;		//e
		B3[5]+=W*Z*Z*w;		//f
		B3[6]+=W*X*Y*w;		//g
		B3[7]+=W*Y*Z*w;		//h
		B3[8]+=W*X*Z*w;		//i
		
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

	for(int L=0;L<81;L++) matrix_val[L]=matrix[L];//行列の値を保存
			
	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;
	int Yflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Zflag=OFF;
	for(int k=0;k<9;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flagがOFFならどのみちdudxはゼロ。

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flagがOFFならどのみちdvdyはゼロ。

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flagがOFFならどのみちdwdzはゼロ。

	////計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の標準偏差
	double W=0;		//重みの総和
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]*X*X+B1[4]*Y*Y+B1[5]*Z*Z+B1[6]*X*Y+B1[7]*Y*Z+B1[8]*X*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]*X*X+B2[4]*Y*Y+B2[5]*Z*Z+B2[6]*X*Y+B2[7]*Y*Z+B2[8]*X*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]*X*X+B3[4]*Y*Y+B3[5]*Z*Z+B3[6]*X*Y+B3[7]*Y*Z+B3[8]*X*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	//////


	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];
			
	double div=dudx+dvdy+dwdz;

	delete [] weight;

	return div;
}

//divergence2における、3次元2次近似を行う関数ver.2 未知数がひとつ多い
double cacl_WLSM_divu_D3_order2_2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	//P=aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+iΔzΔx+Pとおくと、
	///係数行列は
	///   ΣΔx2      ΣΔxΔy    ΣΔxΔz    ΣΔx3       ΣΔxΔy2    ΣΔxΔz2    ΣΔx2Δy     ΣΔxΔyΔz  ΣΔx2Δz    ΣΔx    a = ΣΔxf  
	///   ΣΔxΔy    ΣΔy2      ΣΔyΔz    ΣΔx2Δy    ΣΔy3       ΣΔyΔz2    ΣΔxΔy2     ΣΔy2Δz    ΣΔxΔyΔz  ΣΔy    b = ΣΔyf
	///   ΣΔxΔz    ΣΔyΔz    ΣΔz2      ΣΔx2Δz    ΣΔy2Δz    ΣΔz3       ΣΔxΔyΔz   ΣΔyΔz2    ΣΔxΔz2    ΣΔz    c = ΣΔzf
	///   ΣΔx3      ΣΔx2Δy   ΣΔx2Δz   ΣΔx4       ΣΔx2Δy2   ΣΔx2Δz2   ΣΔx3Δy     ΣΔx2ΔyΔz ΣΔx3Δz    ΣΔx2   d = ΣΔx2f
	///   ΣΔxΔy2   ΣΔy3      ΣΔy2Δz   ΣΔx2Δy2   ΣΔy4       ΣΔy2Δz2   ΣΔxΔy3     ΣΔy3Δz    ΣΔxΔy2Δz ΣΔy2   e = ΣΔy2f
	///   ΣΔxΔz2   ΣΔyΔz2   ΣΔz3      ΣΔx2Δz2   ΣΔy2Δz2   ΣΔz4       ΣΔxΔyΔz2  ΣΔyΔz3    ΣΔxΔz3	 ΣΔz2   f = ΣΔz2f
	///   ΣΔx2Δy   ΣΔxΔy2   ΣΔxΔyΔz ΣΔx3Δy    ΣΔxΔy3    ΣΔxΔyΔz2 ΣΔx2Δy2    ΣΔxΔy2Δz ΣΔx2ΔyΔz ΣΔxΔy g = ΣΔxΔyf
	///   ΣΔxΔyΔz ΣΔy2Δz   ΣΔyΔz2   ΣΔx2ΔyΔz ΣΔy3Δz    ΣΔyΔz3    ΣΔxΔy2Δz  ΣΔy2Δz2   ΣΔxΔyΔz2 ΣΔyΔz h = ΣΔyΔzf
	///   ΣΔx2Δz   ΣΔxΔyΔz ΣΔxΔz2   ΣΔx3Δz    ΣΔxΔy2Δz  ΣΔxΔz3   ΣΔx2Δy Δz ΣΔxΔyΔz2 ΣΔx2Δz2   ΣΔzΔx i = ΣΔxΔzf
	///   ΣΔx       ΣΔy ΣΔz ΣΔx2      ΣΔy2       ΣΔz        ΣΔxΔy     ΣΔyΔz      ΣΔxΔz     ΣΔzΔx     Σ1      P = Σfj

	double le=CON->get_distancebp();
	double matrix_val[100];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
						
		double w=1;
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
		matrix[9]+=X*w;		//ΣXj^2Zjwj
	
		matrix[11]+=Y*Y*w;		
		matrix[12]+=Y*Z*w;		
		matrix[13]+=X*X*Y*w;		
		matrix[14]+=Y*Y*Y*w;		
		matrix[15]+=Y*Z*Z*w;
		matrix[16]+=X*Y*Y*w;
		matrix[17]+=Y*Y*Z*w;
		matrix[19]+=Y*w;
					
		matrix[22]+=Z*Z*w;			
		matrix[23]+=X*X*Z*w;
		matrix[24]+=Y*Y*Z*w;
		matrix[25]+=Y*Y*Y*w;
		matrix[29]+=Z*w;
					
		matrix[33]+=X*X*X*X*w;
		matrix[34]+=X*X*Y*Y*w;
		matrix[35]+=X*X*Z*Z*w;	
		matrix[36]+=X*X*X*Y*w;	
		matrix[37]+=X*X*Y*Z*w;	
		matrix[38]+=X*X*X*Z*w;	
					
		matrix[44]+=Y*Y*Y*Y*w;
		matrix[45]+=Y*Y*Z*Z*w;
		matrix[46]+=X*Y*Y*Y*w;
		matrix[47]+=Y*Y*Y*Z*w;
		matrix[48]+=X*Y*Y*Z*w;

		matrix[55]+=Z*Z*Z*Z*w;	//6行目
		matrix[56]+=X*Y*Z*Z*w;
		matrix[57]+=Y*Z*Z*Z*w;
		matrix[58]+=X*Z*Z*Z*w;

		matrix[99]+=w;
		//7～9行目はすべて既存の要素から転用が可能

		B1[0]+=U*X*w;		//a
		B1[1]+=U*Y*w;		//b
		B1[2]+=U*Z*w;		//c
		B1[3]+=U*X*X*w;		//d
		B1[4]+=U*Y*Y*w;		//e
		B1[5]+=U*Z*Z*w;		//f
		B1[6]+=U*X*Y*w;		//g
		B1[7]+=U*Y*Z*w;		//h
		B1[8]+=U*X*Z*w;		//i
		B1[9]+=U*w;

		B2[0]+=V*X*w;		//a
		B2[1]+=V*Y*w;		//b
		B2[2]+=V*Z*w;		//c
		B2[3]+=V*X*X*w;		//d
		B2[4]+=V*Y*Y*w;		//e
		B2[5]+=V*Z*Z*w;		//f
		B2[6]+=V*X*Y*w;		//g
		B2[7]+=V*Y*Z*w;		//h
		B2[8]+=V*X*Z*w;		//i
		B2[9]+=V*w;	

		B3[0]+=W*X*w;		//a
		B3[1]+=W*Y*w;		//b
		B3[2]+=W*Z*w;		//c
		B3[3]+=W*X*X*w;		//d
		B3[4]+=W*Y*Y*w;		//e
		B3[5]+=W*Z*Z*w;		//f
		B3[6]+=W*X*Y*w;		//g
		B3[7]+=W*Y*Z*w;		//h
		B3[8]+=W*X*Z*w;		//i
		B3[9]+=W*w;
	}
	matrix[10]=matrix[1];
	matrix[18]=matrix[7];

	matrix[20]=matrix[2];
	matrix[21]=matrix[12];
	matrix[24]=matrix[16];
	matrix[26]=matrix[7];
	matrix[27]=matrix[15];
	matrix[28]=matrix[5];

	for(int k=0;k<=2;k++) matrix[30+k]=matrix[3+10*k];//30～32要素
	matrix[39]=matrix[0];

	for(int k=0;k<=3;k++) matrix[40+k]=matrix[4+10*k];//40～43要素
	matrix[49]=matrix[11];

	for(int k=0;k<=4;k++) matrix[50+k]=matrix[5+10*k];//50～54要素
	matrix[59]=matrix[22];

	for(int k=0;k<=5;k++) matrix[60+k]=matrix[6+10*k];//60～65要素
	matrix[66]=matrix[34];
	matrix[67]=matrix[48];
	matrix[68]=matrix[37];
	matrix[69]=matrix[1];

	for(int k=0;k<=6;k++) matrix[70+k]=matrix[7+10*k];//70～76要素
	matrix[77]=matrix[54];
	matrix[78]=matrix[56];
	matrix[79]=matrix[12];
	
	for(int k=0;k<=7;k++) matrix[80+k]=matrix[8+10*k];//80～87要素
	matrix[88]=matrix[35];
	matrix[89]=matrix[20];

	for(int k=0;k<=8;k++) matrix[90+k]=matrix[9+10*k];//90～98要素

	matrix[99]+=1;//自身
	B1[9]+=PART[i].u[A_X];
	B2[9]+=PART[i].u[A_Y];
	B3[9]+=PART[i].u[A_Z];

	for(int L=0;L<100;L++) matrix_val[L]=matrix[L];//行列の値を保存
			
	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;
	int Yflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Zflag=OFF;
	for(int k=0;k<10;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flagがOFFならどのみちdudxはゼロ。

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flagがOFFならどのみちdvdyはゼロ。

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flagがOFFならどのみちdwdzはゼロ。

	////計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の標準偏差
	double W=0;		//重みの総和
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]*X*X+B1[4]*Y*Y+B1[5]*Z*Z+B1[6]*X*Y+B1[7]*Y*Z+B1[8]*X*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]*X*X+B2[4]*Y*Y+B2[5]*Z*Z+B2[6]*X*Y+B2[7]*Y*Z+B2[8]*X*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]*X*X+B3[4]*Y*Y+B3[5]*Z*Z+B3[6]*X*Y+B3[7]*Y*Z+B3[8]*X*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	//////


	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];
			
	double div=dudx+dvdy+dwdz;

	delete [] weight;

	return div;
}

///quickMPS用ポスト処理関数
void post_processing3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double TIME)
{
	
	///restart用ﾌｧｲﾙ出力
	if(t==CON->get_step() || t%CON->get_autosave()==0)
	{
		////////restart用に粒子数と粒子データを記録
		ofstream hoge1("number.dat");
		hoge1<<particle_number<<endl;
		hoge1<<TIME<<endl;
		hoge1.close();
		
		FILE *hoge6;
		hoge6=fopen("initial_input.dat","w");//mps_input.datとは区別しないと、restartが失敗したとき困る
		for(int i=0;i<particle_number;i++)///流体解析粒子から先に記述
		{
			fprintf( hoge6, "%d\t",i);
			fprintf( hoge6, "%5.15f\t",PART[i].r[A_X]);
			fprintf( hoge6, "%5.15f\t",PART[i].r[A_Y]);
			fprintf( hoge6, "%5.15f\t",PART[i].r[A_Z]);
			fprintf( hoge6, "%5.15f\t",PART[i].u[A_X]);  //速度x成分
			fprintf( hoge6, "%5.15f\t",PART[i].u[A_Y]);  //速度y成分
			fprintf( hoge6, "%5.15f\t",PART[i].u[A_Z]);  //速度z成分
			fprintf( hoge6, "%5.15f\t",PART[i].P); //圧力
			fprintf( hoge6, "%5.15f\t",PART[i].h); //エンタルピー
			fprintf( hoge6, "%5.15f\t",PART[i].val);
			fprintf( hoge6, "%d\n",PART[i].type);
			fprintf( hoge6, "%d\n",PART[i].materialID);
			fprintf( hoge6, "%d\n",PART[i].surface);
			fprintf( hoge6, "%d\n",PART[i].toBEM);
		}
		fclose(hoge6);
		//////////////////////////////////////////////
	}
}

//粒子が解析領域の外にでていないかチェック
int check_position(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int *particle_number)
{
	int sw=OFF;	//本関数で返す値。OFFなら粒子数に変化なし。ONなら変化したという印
	int num=0;	//消滅する粒子数
	double le=CON->get_distancebp();
	int *flag=new int [fluid_number];	//ONなら領域内 OFFなら領域外

	double Xmax=CON->get_maxX(); double Xmin=CON->get_minX();
	double Ymax=CON->get_maxY(); double Ymin=CON->get_minY();
	double Zmax=CON->get_maxZ(); double Zmin=CON->get_minZ();

	double dx=CON->get_dx()*le;	//格子幅

	Xmax-=dx; Ymax-=dx; Zmax-=2*dx;		//保険をかけて1格子分内側に境界をとる。これより外側なら粒子を消す
	Xmin+=dx; Ymin+=dx; Zmin+=2*dx;

	vector<mpsparticle>::iterator p,p0;//反復子
	p0=PART.begin();

	for(int i=0;i<fluid_number;i++) 
	{
		flag[i]=ON;
		if(PART[i].r[A_X]<Xmin || PART[i].r[A_X]>Xmax) flag[i]=OFF;
		else if(PART[i].r[A_Y]<Ymin || PART[i].r[A_Y]>Ymax) flag[i]=OFF;
		else if(PART[i].r[A_Z]<Zmin || PART[i].r[A_Z]>Zmax) flag[i]=OFF;
	}//flag[i]が求まった

	int min_nei=3;//5;//CON->get_min_nei();
	if(CON->get_dimention()==3) min_nei=0;//5;
	for(int i=0;i<fluid_number;i++) if(PART[i].N<1) flag[i]=OFF;//周辺粒子数が少ない粒子も削除?

	//速度に上限値を設け、上回っていたら削除
	double limit_U=CON->get_max_speed();			//速度がこれを上回る粒子は削除する
	for(int i=0;i<fluid_number;i++)
	{
		double speed=0;
		for(int D=0;D<3;D++) speed+=PART[i].u[D]*PART[i].u[D];
		speed=sqrt(speed);
		if(speed>limit_U) flag[i]=OFF;
	}////*/

	//近接しすぎている場合、番号の若いほうを削除
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double r[3];
			for(int D=0;D<3;D++) r[D]=PART[j].r[D]-PART[i].r[D];
			double dis=sqrt(r[A_X]*r[A_X]+r[A_Y]*r[A_Y]+r[A_Z]*r[A_Z]);
			double L=PART[i].L;
			if(L>PART[j].L) L=PART[j].L;
			if(dis<L*0.2)
			{
				if(PART[j].type!=FLUID) flag[i]=OFF;//粒子jが壁粒子なら、流体側を消す。
				else	//流体同士の場合
				{
					if(PART[i].surface==ON)
					{
						if(PART[j].surface==ON) flag[i]=OFF;
						else flag[j]=OFF;	//流体iが表面で、流体jが内部なら、jを消す
					}
					else flag[i]=OFF;
				}
			}
		}
	}////////*/


	//for(int i=0;i<fluid_number;i++) if(PART[i].val!=0) flag[i]=ON;

	for(int i=0;i<fluid_number;i++) if(flag[i]==OFF) num++;
	/*///
	//消失する粒子の電荷を周辺粒子へ分配 (valが電荷密度でしかない場合は、そのまま消失させればよい)
	if(CON->get_conductive_appro()==ON && CON->get_EM_calc_type()==1)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(flag[i]==OFF)
			{
				int jnb=0;									//粒子iの周辺粒子のうち、flag=ONな粒子数
				for(int k=0;k<PART[i].N;k++)  if(flag[PART[i].NEI[k]]==ON) jnb++;
				if(jnb>0)
				{
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(flag[j]==ON) PART[j].val+=PART[i].val/jnb;			//導体近似する場合は、消失する粒子の電荷を近隣粒子に分配しておく。
					}
					PART[i].val=0;
				}
			}
		}
	}
	/////*/

	
	if(num>0)//領域外粒子を検知したなら
	{
		sw=ON;
		int *erase_id=new int[num];
		int count=0;
		double val=0;
		for(int i=0;i<fluid_number;i++)
		{
			if(flag[i]==OFF)
			{
				erase_id[count]=i;//消すべきidを記憶
				val+=PART[i].val;
				count++;
			}
		}
		if(val!=0) cout<<"val="<<val<<"が消失"<<endl;
		
		for(int i=0;i<num;i++)
		{
			p=PART.begin();
			p+=erase_id[i];
			
			PART.erase(p);
			for(int j=i+1;j<num;j++)
			{
				erase_id[j]=erase_id[j]-1;//1つ値をさげる
			}
		}
		delete [] erase_id;

		//idがずれたから戻す
		for(int i=0;i<PART.size();i++) if(PART[i].id!=i) PART[i].id=i; 
	}

	if(sw==ON)
	{
		cout<<"領域外粒子を探知 "<<num<<"個の粒子を消去--現在の粒子数="<<PART.size()<<endl;
		*particle_number=*particle_number-num;
	}
	
	delete [] flag;

	return sw;
}

//粒子が近づきすぎるのを防ぐ関数
void modify_position(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double dt,int particle_number)
{
	double le=CON->get_distancebp();

	double *as_i=new double[fluid_number];		//粒子iとして、この関数に参加するか、しないか
	double *as_j=new double[particle_number];		//粒子iとして、この関数に参加するか、しないか
	double *new_r[DIMENTION];
    for(int D=0;D<DIMENTION;D++) new_r[D]=new double [fluid_number];//新しい位置ベクトル

	if(CON->get_modify_position()==1)
	{
		for(int i=0;i<fluid_number;i++) as_i[i]=ON;		//すべての粒子で計算
		for(int i=0;i<particle_number;i++) as_j[i]=ON;
	}
	else if(CON->get_modify_position()==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON) 
			{
				as_i[i]=ON;		//表面粒子で計算
				as_j[i]=ON;		//表面粒子で計算
			}
			else as_i[i]=OFF;
		}
		for(int i=fluid_number;i<particle_number;i++) as_j[i]=ON;//壁はすべて考慮
	}

	

	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) new_r[D][i]=PART[i].r[D];	//初期化
	for(int i=0;i<fluid_number;i++)
	{
		if(as_i[i]==ON)
		{
			double mindis=le;
			mindis=100;			//距離の短いものを遠ざけるだけでなく、遠いものをleにしたいときはこっち
			int J=i;			//最近接粒子
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(as_j[j]==ON)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					if(dis<mindis)
					{
						mindis=dis;
						J=j;
					}
				}
			}
			if(J>i)
			{
				double L=le-mindis;//開くべき距離
				double dL[DIMENTION];
				for(int D=0;D<DIMENTION;D++) dL[D]=PART[J].r[D]-PART[i].r[D];
				if(J!=i && PART[J].type==FLUID)//leより近接している流体粒子があったなら
				{
					for(int D=0;D<DIMENTION;D++)
					{
						double dU=0.5*L/dt;	//変化すべき速度
					//	PART[J].r[D]+=dL[D]/mindis*dU*dt;
						new_r[D][J]+=dL[D]/mindis*dU*dt;

					//	PART[i].r[D]-=dL[D]/mindis*dU*dt;
						new_r[D][i]-=dL[D]/mindis*dU*dt;
					}				
				}
				else if(J!=i && PART[J].type!=FLUID)//leより近接している壁粒子があったなら
				{
					for(int D=0;D<DIMENTION;D++)
					{
						double dU=L/dt;	//変化すべき速度
					//	PART[i].r[D]-=dL[D]/mindis*dU*dt;
						new_r[D][i]-=dL[D]/mindis*dU*dt;
					}				
				}
			}
		}
	}
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) PART[i].r[D]=new_r[D][i];	//代入

	delete [] as_i;
	delete [] as_j;
	for(int D=0;D<DIMENTION;D++) delete [] new_r[D];
}

//各粒子の動粘性係数計算関数
void calc_vis_value(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *vis,double dt,int t,int particle_number)
{
	double R[3];
	double le=CON->get_distancebp();
	double r=CON->get_re()*le;
	int d=CON->get_dimention();
	double RR=8.314;						//ガス定数8.314[J/mol/K]
	double TT=CON->get_roomT();				//室温[K]
	double *Q=new double[fluid_number];		//活量 [J/mol] A6061:145000, A1050:156888,
	double *alpha=new double[fluid_number];	//[1/Pa]		A6061:0.045*1e-6, A1050:0.037*1e-6,
	double *A=new double[fluid_number];		//[1/sec]	A6061:/8.8632*1e6? それともexp(19.3)?, A1050:exp(26.69),
	double *N=new double[fluid_number];		//指数			A6061:3.55, A1050:3.84,
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1)

		{
			Q[i]=158300;			//A1100のパラメータ
			alpha[i]=0.045*1e-6;
			A[i]=exp(24.67);
			N[i]=5.66;
		}
		else if(PART[i].materialID==2)
		{
			Q[i]=166900;			//A5056のパラメータ
			alpha[i]=0.015*1e-6;
			A[i]=exp(23.05);
			N[i]=4.82;
		}
	}

	double V=get_volume(CON);				//粒子の体積
	double *density=new double [fluid_number];//各粒子の密度格納
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		else if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	double *mass=new double [fluid_number];	//各粒子の質量格納
	for(int i=0;i<fluid_number;i++) mass[i]=density[i]*V;
		
	double co=0.9;								//仕事が熱になる(ロス)割合(0～1)
	double *Cp=new double [fluid_number];		//比熱
	double *MP=new double [fluid_number];		//融点
	double *latent_H=new double [fluid_number];		//融点

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1)
		{
			Cp[i]=CON->get_Cp();
			MP[i]=CON->get_MP();
			latent_H[i]=CON->get_latent_H();
		}
		else if(PART[i].materialID==2)
		{
			Cp[i]=CON->get_Cp2();
			MP[i]=CON->get_MP2();
			latent_H[i]=CON->get_latent_H2();
		}
	}

	double *sigma=new double [fluid_number];//各粒子の相当応力格納
	double *ep=new double [fluid_number];	//各粒子の相当ひずみ速度格納
	double *T=new double [fluid_number];	//各粒子の温度格納[K]
	double *heat=new double [fluid_number];	//各粒子の発熱量格納[J]

	//unsigned int timeA=GetTickCount();
	//#pragma omp parallel for
	for(int i=0;i<fluid_number;i++)
	{
		double hs0=mass[i]*Cp[i]*MP[i];//融解開始点のエンタルピー
		double hs1=hs0+latent_H[i]*mass[i];    //融解終了点のエンタルピー	

		///温度Ｔの計算
		if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];///固体
		else if(hs0<=PART[i].h && PART[i].h<=hs1) T[i]=MP[i];//融点
		else if(hs1<PART[i].h) T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];//液体
		////////////*/
		
		double W=0;
		double nensei=1e8;
		double U[3][3];//∂u[i]/∂xj格納
		for(int n=0;n<3;n++) for(int m=0;m<3;m++) U[n][m]=0;//初期化

		for(int k=0;k<PART[i].N;k++)
		{       
			int j=PART[i].NEI[k];
			
			for(int D=0;D<3;D++) R[D]=PART[j].r[D]-PART[i].r[D];
			double dis=sqrt(R[A_X]*R[A_X]+R[A_Y]*R[A_Y]+R[A_Z]*R[A_Z]);
		
			double w=kernel(r,dis);
			W+=w;
		
			for(int n=0;n<3;n++) for(int m=0;m<3;m++) U[n][m]+=(PART[j].u[n]-PART[i].u[n])*R[m]*w/(dis*dis);
	    }
		if(W!=0) for(int n=0;n<3;n++) for(int m=0;m<3;m++) U[n][m]*=d/W;
		
		double E[3][3];//εij=0.5*(∂ui/∂xj+∂uj/∂xi)格納
		for(int n=0;n<3;n++) for(int m=0;m<3;m++) E[n][m]=0.5*(U[n][m]+U[m][n]);

		////ひずみ計算
		double ep1=0;
		for(int n=0;n<3;n++) for(int m=0;m<3;m++) ep1+=E[n][m]*E[n][m];
		ep1*=2.0/3;
		ep[i]=sqrt(ep1);
		/////
		
		double Z=ep[i]*exp(Q[i]/(RR*T[i]));		//Zener-Hollomon parameter

		double G=Z/A[i];
		double H=pow(G,1.0/N[i])+sqrt(pow(G,2.0/N[i])+1);
		sigma[i]=log(H)/alpha[i];
		
		if(ep[i]>0) nensei=sigma[i]/(3*ep[i]);
		
		vis[i]=nensei/density[i];					//動粘性係数

		//温度場を考慮するなら、粘性による発熱を計算する。
		//heat[i]=V*sigma[i]*ep[i]*dt*co;//塑性仕事増分dW=σ(dε)のco%が熱になると仮定	
		//PART[i].h+=heat[i];
		//塑性仕事増分dW=σ(dε)のco%が熱になると仮定
		PART[i].heat_generation+=sigma[i]*ep[i]*co;		//heat_generationに発熱量[W]を格納


		/*/熱になったぶん、粒子の運動エネルギーをさげる
		double U=0;//粒子の速度
		for(int D=0;D<d;D++) U+=PART[i].u[D]*PART[i].u[D];
		U=sqrt(U);
		double B=sqrt(2*heat/(V*density));//減速される速度
			
		if(U>1e-10 && U>B)
		{cout<<"速度減速"<<endl;
			for(int D=0;D<d;D++) PART[i].u[D]-=PART[i].u[D]/U*B;
		}
		///*/

	}////*/
	//cout<<(GetTickCount()-timeA)*0.001<<endl;

	//動粘性値出力
	if(t==1 || t%10==0)
	{
		ofstream fp("dy_vis.dat");
		ofstream fq("dy_vis2.dat");
		ofstream fr("heat.dat");
		ofstream fs("ep_speed.dat");
		ofstream ft("sigma.dat");
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Z]>0.004-0.5*le && PART[i].r[A_Z]<0.004+0.5*le)
			{
				fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<vis[i]<<endl;
				fr<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<heat[i]/V<<endl;//単位体積あたりの発熱量を出力[J/m3]
				fs<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<ep[i]<<endl;//単位体積あたりの発熱量を出力[J/m3]
				ft<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<sigma[i]<<endl;//単位体積あたりの発熱量を出力[J/m3]
			}
			if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fq<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<vis[i]<<endl;
		}
		fp.close();
		fq.close();
		fr.close();
		fs.close();
		ft.close();
	}
	/////////////////*/

	delete [] sigma;
	delete [] ep;
	delete [] T;
	delete [] heat;
	delete [] density;
	delete [] mass;
	delete [] Cp;
	delete [] MP;
	delete [] latent_H;

	delete [] Q;
	delete [] alpha;
	delete [] A;
	delete [] N;
}

//温度AVSファイル出力関数
void output_viscousity_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	double cross_section=CON->get_speed_face_p();
	//t=1;//いまはわざと毎ステップ上書き

	//sprintf_s(filename,"pressure/pressure%d",t);//フォルダを作成して管理する場合はこちら
	sprintf_s(filename,"vis_XZ%d",t);//他のファイルと同じ階層に生成するならこちら
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
			if(PART[i].type==FLUID)
			{
				if(PART[i].r[A_Y]<cross_section+0.5*le && PART[i].r[A_Y]>cross_section-0.5*le)	
				//if(PART[i].r[A_Y]<0.006+0.5*le && PART[i].r[A_Y]>0.006-0.5*le)	
				{
					double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
					double y=PART[i].r[A_Y]*1.0E+05;
					double z=PART[i].r[A_Z]*1.0E+05;
					double P=PART[i].vis;
					fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
					n++;
				}
			}
		}
	}
	else if(CON->get_dimention()==2)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type!=FLUID)
			{
			double x=PART[i].r[A_X]*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
			double y=PART[i].r[A_Y]*1.0E+05;
			double z=PART[i].r[A_Z]*1.0E+05;

			//double x=PART[i].r[A_X];//*1.0E+05;	//rは非常に小さい値なので10^5倍しておく
			//double y=PART[i].r[A_Y];//*1.0E+05;
			//double z=PART[i].r[A_Z];//*1.0E+05;
			double P=PART[i].vis;
			//double P=PART[i].heat_gene_before1;
			fout << P << "\t" << x << "\t" << y << "\t" << z << endl;
			n++;
		}
		}
	}
	fout.close();
	//sprintf_s(filename,"pressure/pressure%d.fld",t);//フォルダを作成して管理する場合はこちら
	sprintf_s(filename,"vis_XZ%d.fld",t);//他のファイルと同じ階層に生成するならこちら
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
	fout2 << "label=viscousity" << endl << endl;
	//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//フォルダを作成して管理する場合はこちら
	fout2 << "variable 1 file=vis_XZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    1 file=vis_XZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    2 file=vis_XZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2 << "coord    3 file=vis_XZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//他のファイルと同じ階層に生成するならこちら
	fout2.close();
}



//各粒子の物性値を温度に応じて変化させる//産総研のデータベース　http://riodb.ibase.aist.go.jp/TPDB/DBGVsupport/detail/aluminum.html
void calc_physical_property(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *val,int particle_number,int sw)
{
	double MP=CON->get_MP();//融点
	
	for(int i=0;i<fluid_number;i++)
	{
		double T=PART[i].T;

		if(CON->get_material()==Al)
		{
			if(sw==1)//密度
			{
				if(PART[i].T<=MP) val[i]=-0.207*PART[i].T+2757; //産総研のデータから線形近似
				if(PART[i].T>MP) val[i]=2377.23-0.311*(PART[i].T-933.47);//産総研の近似式
			}
			if(sw==2)//粘性
			{
				if(PART[i].T<=MP)
				{
					//cout<<"粘性:固体?"<<endl;
					//val[i]=1e10;//固体
					val[i]=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));
				}
					if(PART[i].T>MP) val[i]=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));//産総研の近似式
			}
			if(sw==3)//表面張力係数
			{
				if(PART[i].T<=MP) val[i]=0;//固体の場合、表面張力自体が存在しない
				if(PART[i].T>MP) val[i]=-0.0001*T+0.988;//アルミニウム技術便覧のデータから線形近似

			}
			if(sw==4)//抵抗率
			{
				if(PART[i].T<=MP) val[i]=-T*T*3e-14+T*9E-11-2e-9; //産総研のデータから2次多項式近似
				if(PART[i].T>MP) val[i]=-T*T*2e-14+T*2E-10+8e-8;//産総研のデータから2次多項式近似
			}
			if(sw==5)//動粘性
			{
				if(PART[i].T<=MP)
				{
					//cout<<"動粘性:固体?"<<endl;
					double val_d=-0.207*PART[i].T+2757; //産総研のデータから線形近似
					double val_n=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));
					//double val_n=1e10;
					val[i]=val_n/val_d;
				}
				if(PART[i].T>MP)
				{
					double val_d=2377.23-0.311*(PART[i].T-933.47);//産総研の近似式
					double val_n=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));
					val[i]=val_n/val_d;
				}
			}
		}
		else if(CON->get_material()==H2O)
		{
			if(sw==1)//密度
			{
				if(T<=MP) {}//個体なら何もしない
				//else if(T>MP) val[i]=(0.99986775 + (T-273.15)*6.7866875E-5 -pow(T-273.15,2.0)*9.09099173E-6 + pow(T-273.15,3.0)*1.02598151E-7 -pow(T-273.15,4.0)*1.3502904E-9 + pow(T-273.15,5.0)*1.32674392E-11 - pow(T-273.15,6.0)*6.461418E-14)*999.975; //産総研のデータから線形近似 DENS=(0.99986775+6.7866875E-5*(T-273.15)-9.09099173E-6*(T-273.15)^2+1.02598151E-7*(T-273.15)^3-1.3502904E-9*(T-273.15)^4+1.32674392E-11*(T-273.15)^5-6.461418E-14*(T-273.15)^6)*999.975
				else if(T>MP) val[i] = -T*T*T*T* 1.67669860531849E-6 + T*T*T*2.15183081454928E-3 - T*T*1.03632161524970 + T*2.21554244733677E2 - 1.67188884200996E4;
				//if(i==0) cout<<CON->get_density()<<" "<<val[i]<<endl;
			}
			if(sw==2)//粘性
			{
			}
			if(sw==3)//表面張力係数
			{
			}
			if(sw==4)//抵抗率
			{
			}
			if(sw==5)//動粘性
			{
				if(T<=MP) {cout<<"H2O固体"<<endl;}//個体なら何もしない
				//else if(T>MP) val[i]=0.0005 - T*5E-6 + pow(T,2.0)*2E-8 - pow(T,3.0)*5E-11 +pow(T,4.0)*3E-14;
				else if(T>MP) val[i]= T*T*T*T*3.40023563631489E-14 - T*T*T*4.64957159378912E-11 + T*T*2.38796292580695E-8 - T*5.46534377544010E-6 + 4.71249468774775E-4;
				//if(i==0) cout<<"T="<<T<<" vis="<<CON->get_vis()<<" "<<val[i]<<endl;
			}
		}
	}

}

///移動粒子移動関数
void move_particle(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double dt)
{
	cout<<"移動壁粒子を移動"<<endl;
	int direction=CON->get_move_u_dirct();	//移動粒子を移動させる方向 現在は±X方向=±1,±Y方向=±2,±Z方向=±3
	double speed=CON->get_move_speed();		//移動粒子の移動速度[m/s]
	int D;									//移動方向の本プログラムにおける対応する次元 A_X=0;A_Y=1;A_Z=2;

	if(direction>0) D=direction-1;
	else if(direction<0)
	{
		D=-direction-1;
		speed*=-1;				//速度を反転
	}

	for(int i=fluid_number;i<particle_number;i++)
	{
		if(PART[i].toBEM==MOVE)
		{
			PART[i].r[D]+=speed*dt;//粒子を移動
		}
	}
}

///特殊ファイル出力関数
void output_special_graph(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double TIME,int t)
{
	if(t==1 && CON->get_restart()==OFF && CON->get_model_number()==20)
	{
		ofstream fout("length.dat");//File作成
		fout.close();
	}


	if(t%20==0 || t==1)
	{
		if(CON->get_model_number()==20)
		{
			ofstream ip("length.dat",ios :: app);
		
			double maxh=-100;
			double minh=100;
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].type==BOFLUID)
				{
					if(PART[i].r[A_Z]>maxh) maxh=PART[i].r[A_Z];
					if(PART[i].r[A_Z]<minh) minh=PART[i].r[A_Z];
				}
			}
			ip<<TIME<<" "<<maxh-minh<<endl;

			ip.close();
		}
	}
}

//粒子体積計算関数
double get_volume(mpsconfig *CON)
{
	double V=0;//体積
	double le=CON->get_distancebp();
	if(CON->get_model_set_way()==0)	//正方格子のとき
	{
		if(CON->get_dimention()==2){V=le*le;}
		else	{V=le*le*le;}
	}
	else if(CON->get_model_set_way()==1)	//細密格子のとき
	{
		if(CON->get_dimention()==2){V=sqrt(3.0)/2*le*le;}
		else V=le*le*le/sqrt(2.0);
	}	
	else cout<<"モデルの積み方が不定です 体積を計算できません"<<endl;
	return V;
}

//粒子数密度のずれ測定関数
void output_particle_density(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double n0,int particle_number,int t)
{
	///初期粒子数密度分布を出力する

	ofstream fp("initial_n0.dat");
	double le=CON->get_distancebp();

	if(CON->get_dimention()==2) {for(int i=0;i<particle_number;i++) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].PND<<endl;}
	else if(CON->get_dimention()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].PND<<endl;
		}
	}
	fp.close();
}

void delete_particle(mpsconfig &CON,vector<mpsparticle> &PART,int *particle_number,int *fluid_number,double n0_4,int t)
{
	double le=CON.get_distancebp();
	double begin_num=*fluid_number;
	int count=0;

	double Xmax=CON.get_maxX(); double Xmin=CON.get_minX();
	double Ymax=CON.get_maxY(); double Ymin=CON.get_minY();
	double Zmax=CON.get_maxZ(); double Zmin=CON.get_minZ();

	double ZM=0.0;
	for(int i=0;i<*particle_number;i++)
	{
		if(PART[i].r[A_Z]>ZM) ZM=PART[i].r[A_Z];
	}

	int i=0;
	while(1)
	{
		int flag=0;
		double X=PART[i].r[A_X]; double Y=PART[i].r[A_Y]; double Z=PART[i].r[A_Z];

		if(PART[i].type==FLUID)
		{
			//領域外に出た場合 (電極の中に埋もれても削除)
			if(X<Xmin || X>Xmax) flag=1;
			if(Y<Ymin || Y>Ymax) flag=2;
			if(Z<Zmin || Z>Zmax) flag=3;//電極の中に埋もれても削除
			if(CON.get_model_number()==14)	if(PART[i].r[A_Z]<-le*2.0) flag=4;//(静電無化の場合)電極の中に埋もれても削除
			if(CON.get_model_number()==19)	if(PART[i].r[A_Z]>ZM*1.1) flag=8;//(fswの場合)ツールより上にあったら削除
			if(CON.get_model_number()==20 || CON.get_model_number()==21 || CON.get_model_number()==22)//(CC溶解の場合)るつぼにめり込んだら削除
			{
				if(Z>=0)
				{
					if(X*X+Y*Y>(0.03+le)*(0.03+le))
					{
						flag=4;
					}
				}
				else if(Z<0)//内壁円柱部
				{
					if(X*X+Y*Y+Z*Z>(0.03+le)*(0.03+le))
					{
						flag=4;
					}
				}

			}
			if(CON.get_model_number()==1)//(CC溶解の場合)るつぼにめり込んだら削除
			{
				if(X>0.0155 || X<-0.0155 ) flag=4;
				if(Y>0.0255 || Y<-0.0255 ) flag=4;
			}
			//孤立飛散粒子の削除
			///*if(CON.get_process_iso()==1)*/	if(PART[i].fly==ISOLATION) flag=5;
			//粒子位置が不定(-1.#IND)のとき
			if(2*PART[i].r[A_X]!=PART[i].r[A_X]+PART[i].r[A_X]) flag=6;
			else if(2*PART[i].r[A_Y]!=PART[i].r[A_Y]+PART[i].r[A_Y]) flag=6;
			else if(2*PART[i].r[A_Z]!=PART[i].r[A_Z]+PART[i].r[A_Z]) flag=6;

			//ﾚｲﾘｰ分裂専用 初期状態の飛び出た粒子を消す
			//if(CON.get_model_number()==20 && CON.get_current_step()==1 && PART[i].PND4<35)	flag=6;
			
			//速度が過大
			double speed=0;//粒子速度
			for(int D=0;D<DIMENTION;D++) speed+=PART[i].u[D]*PART[i].u[D];
			speed=sqrt(speed);
			if(speed>CON.get_max_speed()) flag=7;
			
		}

		if(flag>0)
		{
			//座標出力ﾃﾝﾌﾟﾚ    <<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z];
			if(flag==1)	cout<<"Ｘ領域外";
			if(flag==2)	cout<<"Ｙ領域外";
			if(flag==3)	cout<<"Ｚ領域外";
			if(flag==4)	cout<<"壁に侵入";
			if(flag==5)	cout<<"孤立飛散";
			if(flag==6)	cout<<"座標不定";
			if(flag==7)	cout<<"最大速度が過大";
			if(flag==8)	cout<<"ツールより上に存在";
			//cout<<" i="<<i<<"/"<<(int)PART.size()<<" ｸﾞﾙｰﾌﾟ"<<PART[i].group<<endl;
			//PART[i]を削除
			vector<mpsparticle>::iterator it=PART.begin();//イテレータ初期化
			it+=i;				//iを指定
			it=PART.erase(it);	//削除
			//cout<<" ｸﾞﾙｰﾌﾟ"<<PART[i].group<<" -> 削除"<<endl;	//←なぜかここで止まる。なぜ？ コメントアウトすると通る。
			*fluid_number-=1;
			*particle_number=(int)PART.size();

			count++;

			//※削除した場合は配列が詰められるのでi++は行わない
		}
		else i++;

		if(i>=(int)PART.size())	break;
	}

	if(count>0)	cout<<"粒子"<<count<<"個を削除\t"<<"流体粒子："<<begin_num<<" -> "<<begin_num-count<<endl;


	//流体粒子数の推移をグラフに出力
	if(t==1)//1ステップ目はファイルをリセット
	{
		ofstream freset("fluid_number.dat");
		freset.close();
	}
	ofstream fout("fluid_number.dat",ios::app);
	fout<<t<<"\t"<<begin_num-count<<endl;
	fout.close();
	
	//流体最大高さを出力
	if(t==1)//1ステップ目はファイルをリセット
	{
		ofstream freset2("fluid_h.dat");
		freset2.close();
	}
	ofstream fout2("fluid_h.dat",ios::app);
	double Zmax2=0;
	for(int i=0;i<*fluid_number;i++)
	{
		if(PART[i].r[A_Z]>Zmax2) Zmax2=PART[i].r[A_Z];
	}
	fout2<<t<<"\t"<<Zmax2<<endl;
	fout2.close();
}

//何かを測定・調査する関数（何かは自分でプログラムをいじって決める）
void check_something(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double n0,int particle_number,int t)
{
	double le=CON->get_distancebp();
	if(CON->get_model_number()==19 && t%10==0)//FSW
	{
		int overN=0;//基準面より上にきた流体粒子数
		int underN=0;//基準面より下にきたツール粒子数
		double standard_H=6e-3;//基準面のZ座標
		for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Z]>standard_H+le*0.5) overN++;
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].toBEM==MOVE)
			{
				if(PART[i].r[A_Z]<standard_H+le*0.5) underN++;
			}
		}
		cout<<"overN/underN="<<overN<<"/"<<underN<<endl;
	}
}

///壁重み関数の計算関数
void calc_wallZ(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,double *mindis,int fluid_number,int out)
{
	cout<<"壁重み関数の計算開始"<<endl;
    ///壁を均一に配置した場合の壁重み関数を計算する。//model_number=0に均一な壁粒子配置を与えているので、0盤目の粒子の座標を変えながら調べる
	double le=CON->get_distancebp();
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	
	//cout<<PART[0].r[A_Y]<<endl;
	//密度を求め、表面判定を行う。表面ならP=0にする
	//omp_set_num_threads(8);//ｽﾚｯﾄﾞ数指定
	//#pragma omp parallel for
	int sol_wz=4000;//riwをどれだけ変えるか
	for(int wz=1;wz<sol_wz;wz++)
	{
		for(int i=0;i<1;i++)//OUTWALL以外の粒子。//OUTWALLの粒子数密度などはいらない
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
			if(PART[i].type==FLUID) PART[i].r[A_Y]=le*0.001*wz;
			for( int j=0;j<out;j++)    
			{	
				     
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
			
			if(i==0)
			{
				//流体粒子数の推移をグラフに出力
				if(wz==1)//1ステップ目はファイルをリセット
				{
					ofstream fout("wallZ_pnd1.dat");
					fout.close();
				}
				if(wz==1)//1ステップ目はファイルをリセット
				{
					ofstream fout2("wallZ_pnd2.dat");
					fout2.close();
				}
				if(wz==1)//1ステップ目はファイルをリセット
				{
					ofstream fout4("wallZ_pnd4.dat");
					fout4.close();
				}

				ofstream fout("wallZ_pnd1.dat",ios::app);
				if(pnd>0) fout<<PART[0].r[A_Y]<<" "<<pnd<<endl;
				fout.close();

				ofstream fout2("wallZ_pnd2.dat",ios::app);
				if(pnd2>0) fout2<<PART[0].r[A_Y]<<" "<<pnd2<<endl;
				fout2.close();

				ofstream fout4("wallZ_pnd4.dat",ios::app);
				if(pnd4>0) fout4<<PART[0].r[A_Y]<<" "<<pnd4<<endl;
				fout4.close();
			}

			
			//PART[i].PND=pnd;
			//PART[i].PND2=pnd2;
			//PART[i].N=N;
			//PART[i].N2=N2;
			//PART[i].N3=N3;
			//if(PART[i].N3>800) cout<<"error in freeon over. Change more than "<<PART[i].N3<<" coods:"<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			////////////////////
	    
		}
	}
	cout<<"壁重みの計算完了"<<endl;

}