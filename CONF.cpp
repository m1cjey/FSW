#include "stdafx.h"	
#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。
#include"define.h"	//#define 格納
#include"CONFIG.h"	//class CON定義

mpsconfig::mpsconfig()
{
	///////解析条件
	dt=5.0e-4;//0.0001;//0.000025*4;//////0.000025*10;//0.0005;//0.000025/2;//0.0000125    //jw法の場合、コイル周期の整数倍にすること         
	step=13760;
	dimention=3; 
	
	maxX=0.02;//0.1;//0.01;
	minX=-0.02;//-0.1;//-0.01;
	maxY=0.03;////0.1;//0.01;             //indexの関係上、Z方向には余裕をもつこと。
	minY=-0.02;//-0.1;//-0.01;
	maxZ=0.02;////0.1;//0.3;//0.01;
	minZ=-0.01;//-0.1;//-0.01;

	material=FSWA1;//Fe,Al,H2Oは流体。FSWA1は固体のA1100

    ///////流体の物性値
	if(material==Fe)//Fe(1536℃)
	{
		density=7015;//8000;//1000;           //water:997.04  ｴﾀﾉｰﾙ:798
		nensei=0.0055;//0.008;//0.001;			//water:0.001 ｴﾀﾉｰﾙ:0.01084 
		sigma=1.872;//1.85;//0.07196;//0.07196;			//water:0.07196 ｴﾀﾉｰﾙ:0.02361 SUS404:1.85
		Cp=640/10;     			//water:4.2[kJ/(kg・K)] 鋼:800J SUS404:645J/(kgK)
		k=28;       			//water:0.6[W/(m・K)]
		latent_H=0;        		//water:334000J/kg  鋼209.3J
		MP=273;		//融点[K]
		CTE=0.002;//2.1e-4;//2.1e-4;	//線膨張係数 SUS410:10.4e-6  水:2.1e-4
		ele_conduc=1/(1.386e-6);
	}

	if(material==Al)//Al(660℃)
	{
		density=2385;//8000;//1000;           //water:997.04  ｴﾀﾉｰﾙ:798
		nensei=0.0013;//0.13;//0.0013;//0.008;//0.001;			//water:0.001 ｴﾀﾉｰﾙ:0.01084 
		sigma=0.865;//0.914;//1.85;//0.07196;//0.07196;			//water:0.07196 ｴﾀﾉｰﾙ:0.02361 SUS404:1.85
		Cp=900;//1080;//640/10;     			//water:4.2[kJ/(kg・K)] 鋼:800J SUS404:645J/(kgK)
		k=238;//94;//28;       			//water:0.6[W/(m・K)]
		latent_H=396567.45;//0;        		//water:334000J/kg  鋼209.3J
		MP=933;		//融点[K]
		CTE=0.002;//2.1e-4;//2.1e-4;	//線膨張係数 SUS410:10.4e-6  水:2.1e-4
		ele_conduc=1/(0.2477e-6);//3.74e7;//1/(0.2477e-6);//3.74e7;
		//ele_conduc=1/(1.386e-6);
	}

	if(material==H2O)//水
	{
		density=997.05;           //water:997.04  ｴﾀﾉｰﾙ:798
		nensei=0.001;			//water:0.001 ｴﾀﾉｰﾙ:0.01084 
		sigma=0.07196;			//water:0.07196 ｴﾀﾉｰﾙ:0.02361 SUS404:1.85
		Cp=4200;//1080;//640/10;     			//water:4.2[kJ/(kg・K)] 鋼:800J SUS404:645J/(kgK)
		k=0.6;//94;//28;       			//water:0.6[W/(m・K)]
		latent_H=334000;//396567.45;//0;        		//water:334000J/kg  鋼209.3J
		MP=273;		//融点[K]
		CTE=2.1e-4;	//線膨張係数 SUS410:10.4e-6  水:2.1e-4
		ele_conduc=1/(0.2477e-6);//3.74e7;//1/(0.2477e-6);//3.74e7;
		//ele_conduc=1/(1.386e-6);
	}

	if(material==FSWA1)//A1100
	{
		density=2710;           //water:997.04  ｴﾀﾉｰﾙ:798
		nensei=100000;			//water:0.001 ｴﾀﾉｰﾙ:0.01084 
		sigma=0.07196;			//water:0.07196 ｴﾀﾉｰﾙ:0.02361 SUS404:1.85
		Cp=900;//1080;//640/10;     			//water:4.2[kJ/(kg・K)] 鋼:800J SUS404:645J/(kgK)
		k=234;//94;//28;       			//water:0.6[W/(m・K)]
		latent_H=396567.45;//0;        		//water:334000J/kg  鋼209.3J//model作成時に潜熱込みで初期状態のエンタルピーを考えている(=流体とみなす)ので、固体の場合0にしておく。相変態やるときはどうする？
		MP=933;	//融点[K]
		CTE=2.1e-4;	//線膨張係数 SUS410:10.4e-6  水:2.1e-4
		ele_conduc=1/(0.2477e-6);//3.74e7;//1/(0.2477e-6);//3.74e7;
		//ele_conduc=1/(1.386e-6);
	}

	///////流体2の物性値//FSWの場合、materialID=1がA1100,materialID=2がA5056
	density2=2710;//2640;       //water:997.04  ｴﾀﾉｰﾙ:798 磁性流体：1400
	nensei2=0.0055;//1000;		//water:0.00893 ｴﾀﾉｰﾙ:0.01084 磁性流体：0.03
	sigma2=1.872;//0.03;		//water:0.07196 ｴﾀﾉｰﾙ:0.02361 SUS404:1.85
	Cp2=900;     		//water:4.2[kJ/(kg・K)] 鋼:800J SUS404:645J/(kgK)
	k2=234;//117;       		//熱伝導率[W/(m・K)]　Al:240 water:0.6
	latent_H2=0;        	//water:334000J/kg  鋼209.3J
	MP2=933;			//融点[K]
	CTE2=0.002;//2.1e-4;//2.1e-4;	//線膨張係数 SUS410:10.4e-6  水:2.1e-4
	
	////////粒子配置用
	fluidwidth=25;   //40    mm単位
	fluid_h=0;	//流体円中部の高さ　CC用
	distancebp=0.0005;//0.001;//0.0015;//0.002;//0.001;//0.0005
	wlength=2;
	height=-0.004;//-0.004;//0.06;//0.102;//0.18;//0.005;    
	tool_angle=0;//FSWにおいて、ツールを傾ける角度(弧度)
	tool_type=0;//FSWにおける、ツール形状0:デフォルト(円柱)　1:円錐 2:円柱裏表
	process_type=2;//FSWにおける過程選択		//0:plunge 1:traverse 2:plung→traverse
	dwelling_time=0.0;//1.0;
	change_step=8000;
	airwall=OFF;//IH釜において、流体上部に壁を配置するかどうか

	///////粒子法用パラメータ
	re=2.1;//2.5;//2.1
	re2=2.1;//2.5;//2.1
	re3=3.1;//3.1
	re4=3.0;//3.1
	beta = 0.90;//0.85;
	dx=4;					//index用格子幅。最大粒子半径を含む整数にすること
    times=2e-8;//50;
	
	////////表面張力関係                                      
	surface_tension=0;      //0=OFF 1=形状 2=粒子間ポテンシャル
	smooth=OFF;				//ｽﾑｰｼﾞﾝｸﾞ　1=ON 0=OFF
	SFT=0;//0;					//1=ON 0=OFF 表面張力(surface tension)の温度(T)依存性スイッチ
	smn=0;					//表面張力のｽﾑｰｼﾞﾝｸﾞの回数
	smth_sumP=0;			//curvのｽﾑｰｼﾞﾝｸﾞ回数 0ならOFF
    wall_direct=1;			//for ver2~5. BDWALLの法線ﾍﾞｸﾄﾙ計算方法　0=無視 1=通常 2=水平方向 3=垂直方向
	non_slip=ON;			//for ver.1 粘着条件 ONならBDWALL上に仮想的に液体が存在する設定。その際はBDWALLの配置に留意
	suf_modify=ON;			//曲率計算で得られた曲面に一致するように粒子位置を修正するか、しないか
	interpolate_curv=OFF;	//曲率が計算されなかった粒子の曲率を、周辺粒子の値から補間するか、しないか
	Cst=0;					//ポテンシャル係数 //いじらない
	C_type=2;				//ポテンシャル係数の決定方法　0:１粒子との総和　1:近藤らの方法（論文）　2:二分割した影響半径の総当り
	C_times=1;				//ポテンシャル係数の倍数
	wall_C=1;//0.5;				//壁との親和力係数

	////////圧力計算
	iteration_count=1;		//圧力計算を行う反復回数 通常は1に設定	0でFEM	
	solution=1;             //0=CG method 1=ICCG 2=ICCG2はﾒﾓﾘ節約
	B_of_P=1;				//圧力計算の解行列 0=教科書 1=速度発散 2=0+1 3=速度発散2 4=3+PND 5=(ni-nk)+(nk-n0)
	w_div=10;				//発散の重み B_of_Pが2のとき、速度発散をどれくらい重視するか
	divU_method=3;			//速度発散の計算方法 1=MPS,2=WLSM,3=自身を通らないWLSM 4=入部ら
	divU_order=3;			//WSLMでdivUを計算するときのオーダー MAX=3
	HL_sw=OFF;				//ﾗﾌﾟﾗｼｱﾝに対し高次の離散化を行う 1=ON 0=OFF
	dir_for_P=1;			//圧力計算の際、表面に非ゼロのディリクレ値を代入するか、しないか 0=OFF 1=表面張力 2=電磁力 3=両方
	initialP=ON;            //圧力計算時に初期値として現圧力を代入するかどうか　1=ON 0=OFF
	P_twice=2;			//粒子数密度一定で位置更新後、速度発散で再度圧力を計算するかどうか 0でoff 整数回ごとに1回粒子数密度の計算をする
	CGep=1e-3;				//圧力計算時における収束判定
	omp_P=OFF;				//圧力計算時のCG法でopenMPを使用するか 1=ON 0=OFF
	negativeP=ON;			//負圧　1=許可 0=不許可
	set_P=OFF;              //0=OFF 1=静水圧 2=関数
	Pgrad=5;				//圧力勾配計算法 1=教科書 2=表面のみ法線 3=2+minP=0のとき表面反発 4:WLSM 5=WLSM(自身を通らない )6=入部
	minP=OFF;				//最少圧力ｽｲｯﾁ 1=ON 0=OFF
	ave_Pdirect=OFF;		//圧力勾配用法線ﾍﾞｸﾄﾙの平均化 1=ON 0=OFF
	Pgrad_order=1;			//圧力勾配ver.4における精度 1=線形 2=2次
	artP_sw=OFF;			//人工圧力計算フラグ 0=OFF 1=artP等方 2=Monaghan
	artP=0;				//人工圧力値 通常は0に設定
	Pgrad_times=20;			//圧力勾配をﾌﾟﾛｯﾄする際の、表面張力に対する倍率 通常は1に設定
	P_AVS=10;				//microAVS用の圧力ファイルを出力するstep間隔。0ならOFF
	Pgrad_limit=1;//0.1;		//圧力勾配による位置修正がleの何％以上にならないよう限界値を設ける 限界値を設定したくない場合は大きな値を代入
	gridless_P=ON;			//圧力をグリッドレス法でとくかどうか//0:mps 1:グリッドレス
	wall_poly=0;			//壁境界をポリゴンで表すかどうか0:壁粒子1:ポリゴン
	interpolate_surface_u=ON;	//圧力勾配による速度修正後の表面粒子の速度を内部粒子の速度で補完するかしないか。
	
	////////温度場関係
	T_field=1;              //温度場解析　　1=ON 0=OFF
	temperature_depend=OFF;		//物性値が温度に依存して変化するかどうか アルミニウム、水のみ対応
	insulate=1;             //壁との断熱状態 0=断熱　1=非断熱 2=局所的に断熱、どこが断熱なのかはソースを確認のこと
	T_laplacian=0;          //温度のﾗﾌﾟﾗｼｱﾝ。　0=教科書　1=λ[i] 2=発散・勾配
	wall_density=2700;//8940;//7850;		//壁の密度[kg/m^3]
	wall_Cp=880;//385;//460;			//壁の比熱[J/(kg・K)]
	wall_k=236;//24;//394;//24;				//壁の熱伝導率[W/(m・K)]
	
	roomT=293;				//室温(20℃) int型でよい
	initialT=293;//MP+10;	//293;		//初期温度
	air_cool=ON;			//空気との熱伝達を考慮するかしないか 1=ON 0=OFF
	T_expansion=OFF;		//熱膨張計算 1=ON 0=OFF
	T_CG=0;					//陰的のとき、温度場解析のソルバー　0:CG 1:iccg
	T_solution=0;			//0:陽的 1:陰的
	T_CGep=1.0e-3;			//収束判定
	buoyant=OFF;			//浮力(密度のブジネスク近似)  1=ON 0=OFF
	TplotZ=0.004;			//3D解析において、XY平面の温度を出力するときのＺ座標
	T_AVS=200;				//microAVS用の温度ファイルを出力するstep間隔。0ならOFF
	output_temperature_face=1;	//0=YZ平面 1=XZ平面	2=XY平面

	////////電磁力計算
	EM_method=0;			//電磁場の解法 0=OFF 1=FEM 2=BEM 3=磁気ﾓｰﾒﾝﾄ法
	EM_calc_type=1;        //0=ﾃﾞﾛｰﾆのみ 1=電場 2=静磁場 3=動磁場 4=磁位
	EM_interval=5;//5;        //電磁場計算を何ｽﾃｯﾌﾟに一回行うか。通常は1に設定
	EM_distance=0.1;			//電磁力計算を行う頻度を、粒子の総移動距離で判定する 0:OFF 格納された値*粒子間距離より粒子の移動距離が大きければfem
	region_shape=0;			//解析領域形状　0=立方体 1=円筒
	XR=0.1;//0.2;//1;//0.2;//0.2;//0.2;//2.00255;//0.2;//0.01;				 //解析領域
	XL=-0.1;//-0.2;//-1;//0.2;//-0.2;//-1;//-0.2;//-2.00255;//-0.2;//-0.01;//2.00255
	YU=0.1;//0.2;//1;//0.2;//0.2;//1;//0.2;//2.00255;//0.2;//0.01;
	YD=-0.1;//-0.2;//-1;//-0.2;//-0.2;//-1;//-0.2;//-2.00255;//-0.2;//-0.01;	
	ZU=0.1;//0.3;//1;//0.2;//0.3;//1;//0.300;//0.30;//2.0;//0.30;//0.3069;//0.3;//0.01;2.0					 //液滴 0.01 コイル:0.15 るつぼ:-0.0003 1mm空域：0.30725 空域なし：0.3
	ZD=-0.1;//-0.1;//-1;//-0.2;//-0.1;//-1;//-0.1;//-0.1;//-2.0;//-0.1;//-0.0931;//-0.1;//-0.01;//-0.003*2;  //液滴 -0.01 コイル:-0.15 るつぼ:-0.0002  1mm空域：-0.09275 空域なし：-0.1
	RU=0.1;//0.1;				//解析領域が円筒形となるときのその半径
	FEM_smn=1;				//電磁力ｽﾑｰｼﾞﾝｸﾞ回数　0ならOFF ﾏｲﾅｽなら表面のみ ﾌﾟﾗｽなら内部も。

	///////メッシュ
	mesher=1;				//メッシュ作成方法 0:吉川さんのプログラム 1:tetgen
	mesh_input=0;			//メッシュをどうやってつくるか 0:自分 1:Magnet //mesher=1のときmagnetのメッシュを読み込む
	remesh_sw=ON;			//ON:remesh領域を考慮し、そこだけremesh OFF:FULL領域を毎回分割
	boxalpha=2.0;			//ｽｰﾊﾟｰﾎﾞｯｸｽの大きさ 1〜10の実数にすること。通常は2に設定
	fine=1;					//節点の自動追加を行うかどうか 0:OFF 1:辺ﾍﾞｰｽ 2:重心ﾍﾞｰｽ
	co_fine=3;//3.5;			//節点自動追加における係数.辺ﾍﾞｰｽの場合は最少辺と最大辺の比率のしきい値
	add_points=100000;//40000;       //自動節点追加数
	poly_memory=50000;		//poly3D()で確保するメモリの数
	CDT_sw=ON;				//制約付きデローニ分割にするかどうか
	defer_f=ON;				//最初のデローニで、不適切な流体節点の追加を後回しにするかどうか
	mesh_smth=0;			//メッシュスムージング回数　0:OFF
	material_judge=0;		//要素の最大辺長さに基づく材質判定を追加で行うかどうか　0:OFF 正の値n:le*nより長い辺があればそれは空気
	del_length=2.0;			//TetGen用 長い要素を削除する閾値 leの何倍以上の辺を持つ要素を消すか
	BD_fluid=OFF;				//リメッシュ境界節点と流体節点で構成される要素を流体とするかどうか　るつぼと流体の間の空気をつぶせるが、粒子と完全な対応にならない
	output_wall_way=0;		//壁粒子を出力してメッシュを作成するかどうか　0=流体のみ 1=流体と壁すべて　2=流体とINWALLのみ
	layer_thin=0;			//空気層の厚さ distancebp*layer_thinが実際の厚さ
	layer_num=0;			//空気層の層数

	////////BEM,FEM関係
	BEM_elm_type=CONSTANT;	//要素ﾀｲﾌﾟ 0:一定要素 1:線形要素
	FEM_elm_type=0;			//要素ﾀｲﾌﾟ 0:節点要素 1:辺要素
	max_DN=25000;            //ﾃﾞｨﾘｸﾚ型境界条件をとる最大節点(辺)数
	FEMCG=5;				//FEMにおける行列解法 1:ICCG 2:並列iccg(parallel-iccg) 3:MRTR 4:対角MRTR 5:ICMRTR
	CGaccl=1.02;//1.02;	//1.3			//ICCG法における加速ファクタ　1のときﾌｧｸﾀＯＦＦ
	FEMCGep=1e-2;//1e-5;    		//ICCGの収束判定
	FEMtimes=1;				//電磁力をﾌﾟﾛｯﾄする際の、表面張力に対する倍率 通常は1に設定
	legend_F=0.002;//2e-7;			//Fn.datの凡例に出力する力 [N] //
	node_sort=0;			//節点番号並び替え　1:テキストどおり　2:結合数に基づくソートをしない、あとは1と同じ
	thinout_fluid=0;		//流体粒子を間引いて節点に渡すかどうか 0:OFF n(自然数):表皮深さ以外のところは、n個に1個受け渡す //球にしかまだ対応していない
	
	///////電界計算
	V=4000;					//電界計算用電圧　もうひとつは０とおいている。
	V_step=1000;//1500;
	r_perm=80;				//比誘電率
	V_con=0;				//電圧条件　0:ﾊﾟﾙｽ　1:ﾘﾆｱ　2:時定数
	initial_V=2000;			//電圧条件ﾘﾆｱのときの初期電圧
	E_times=0;//1e-11;			//電界出力倍率 0なら出力しない
	eleforce=4;				//静電力計算方法 1:節点力法 2:積分面 3:表面力 4:divT
	charge=0;				//電荷考慮 0=OFF 1=電荷密度 2=クーロン力

	///////磁場計算
	J_input_way=1;			//電流密度入手方法 0:自分 1:ソフト
	J0=100;//180000000;		//強制電流密度[A/m2]
	I0=900;//900;//1000;//200;			//
	RP=1.0;//1.28;					//比透磁率
	//ele_conduc=0.7215e6;//3.77e7;//0.7215e6;//3.77e7;//0.7215e6;//1e7;			//電気伝導率（鉄）
	ele_conduc2=5.81e7;//5.814e7;			//電気伝導率（銅）
	Hz=40000;//38500;//40000;//10000			//周波数
	div_Hz=4;				//１周期の分割数(解析精度) 4の倍数がいい 40?
	jheat=0;				//渦電流による発熱を考慮するか　0=OFF 1=ON
	m_force=2;				//電磁力計算方式 0=節点力法 1=体積力 2=kelvin 3=積分面 4=divT 5=VAG 6=積分つき節点力法 7=MC
	NLBHI=0;				//体積力において、要素Ｂから要素Ｈを求める際に非線形性を考慮するか、しないか(non linier B H inverter)
	NLMH=OFF;					//Ｍの算出に非線形性を考慮するか、しないか
	magnet_H=0.01;			//永久磁石の高さ
	magnet_B=0.145;			//永久磁石の強さ[T]
	uniform_B_sw=OFF;		//解析領域中に一様磁場を発生させるか否か 0=OFF 1=ON
	uniform_B=0.1;			//一様磁場の大きさ[T]
	B_times=1;//1;			//ﾌｧｲﾙ出力する際の、磁束密度の倍率
	plot_B_type=1;			//磁束密度出力タイプ 1=ﾍﾞｸﾄﾙ 2=スカラー 0=OFF
	Je_crucible=0;			//渦電流の計算対象　0:FLUID 1:FLUID,CRUCIBLE -1:渦電流が発生しない磁場解析
	J0eqJe=0;				//Je.datに書かれた渦電流を強制電流として与える(未完成)
	m_A=1;					//0:step by step 1:jω法
	jw_interval=5;			//何周期分jw法で求めた波形を使いまわすか　//jw_FaverageがOFFのときのみ利用
	jw_Faverage=ON;			//jw法で求めたベクトルポテンシャルの波形から1周期あたりの電磁力の平均を求めるかどうか　0=OFF,1=ON //OFFは現在使用不可
	static_dirichlet=OFF;	//静的要素(Magnet作成部など)ではベクトルポテンシャルの値を一度求めたらそれ以降その値を使うかどうか //jw法のみ
	eddy_times=1e-10;
	A_phi=ON;				//動磁場解析において、A-φ法を使うかどうか 20130523 辺要素においてA法は対応していない。このスイッチ実装までA法だと思っていた方法は、渦電流連続の式を解いていなかったのでおかしい(φ-Aの項)
	parabolic_node_element=OFF;	//A-φ法において、辺1次節点2次要素を使うかどうか				

	////////各種スイッチ　どのような計算を考慮するか、しないか
	g=-9.8;//-9.8;//-9.8;//-9.8;					//-9.8//
	restart=0;				//1=ON 0=OFF
	autosave=50;			//ｵｰﾄｾｰﾌﾞ間隔。無効にしたいときは大きな数字を代入しておく
	curan=0;				//ｸｰﾗﾝ数条件　0ならＯＦＦ
	modify_position=ON;	//粒子間距離がleより小さい場合これを修正するか、しないか
	vis_calc_type=1;		//粘性項計算手法　0=POSITIVE=陽解法 1=NEGATIVE=陰解法
	wall_adheision=2;		//0=ﾌﾘｰｽﾘｯﾌﾟ 1=ﾉﾝｽﾘｯﾌﾟ  2=2*ﾉﾝｽﾘｯﾌﾟ
	laplacian=2;            //0=教科書　1=λ[i] 2=発散・勾配
	vis_solver=1;			//粘性項を陰解析で解く際の行列ｿﾙﾊﾞｰ 0:CG 1:ICCG
	initial_u=OFF;			//粘性項を陰的に解く際に、初期値として現在速度を入力するか、しないか
	temporary_r=OFF;		//陽解析後に仮の位置を計算するかしないか。 1=ON 0=OFF 通常はONに設定
	fix_center=0;			//1=ON 0=OFF
	freeon=1;				//粒子依存関係関数　1:並列化可能 2:並列不可 4:GPU
	freeon3sw=1;			//freeon3を計算するかしないか 1=ON 0=OFF
	surface_judge2=ON;		//surface_judge2()を使用するかしないか  1=ON 0=OFF
	move_prtcl=ON;			//移動粒子を考慮するかしないか 1=ON 0=OFF
	move_u_dirct=-3;//-3;//2;			//移動粒子を移動する方向　現在は±X方向=±1,±Y方向=±2,±Z方向=±3
	move_speed=1.0e-3;//1.0e-3//1.5e-3;//1.5e-3;//1e-3;//8.333*1e-3;//12.5*1e-3;	//移動粒子の移動速度[m/s]
	move_speed2=8.333e-3;//500mm/min吉川さんの論文引用
	check_something=OFF;		//check_something()を実行するかしないか 1=ON 0=OFF
	check_region=ON;			//条件を満たした粒子を消去するかしないか
	max_speed=20;//5;			//chec_regionがONになっているとき、粒子の速度がこの値以上なら削除する
	AVS_HALF=OFF;				//壁粒子を半分だけ出力する
	adaptive_sw=OFF;		//解像度を可変にするか、しないか 1=ON 0=OFF
	threshold=1;//1e-3;			//アダプティブにする際の、圧力誤差の閾値
	fix_surface=0;			//表面流体を固定するかどうか
	output_forward=OFF;
	output_backward=ON;
	output_another_face=ON;

	model_number=19;
	model_set_way=1;		//modelをセットする方法　0=正方格子 1=MD
	model_inherit=0;		//直前の解析で使用していた粒子モデルを使いまわすかどうか
	ea_model=0;				//静電噴霧において使用するモデル 0:円柱電極 1:製品モデル(cadから製作したstlファイルを読み込み)

	////////速度ﾌﾟﾛｯﾄ変数
	speed_plot_particle=2;	//速度をﾌﾟﾛｯﾄする粒子の種類 1=すべて 2=fluid 3=壁
	speedtimes=1e-3;//2e-2;//1e-2;//1e-3;		//速度ﾌﾟﾛｯﾄ時の、座標に対する速度の倍率
	speed_face=1;			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ	2=XY	//粘性、温度分布の出力面を兼ね備えている
	speed_face2=2;			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ	2=XY	//粘性、温度分布の出力面を兼ね備えている
	speed_face_p=0.0;//-1.0e-3;//0.006;//0.0;		//3D解析時のspeed.datの出力面の座標
	ax_sym_modify=OFF;		//3D時のspeed.datに関して、軸対称による出力修正を行うか否か　1=ON 0=OF
	flat_speed_plot=OFF;	//OFF//水平方向の速度(XY面)をﾌﾟﾛｯﾄするかしないか1=ON 0=OFF		//speed_eachにXY平面の出力を付け加えたため不要かも
	flat_speed_p=0.0;		//0.004//flat_speed.datの出力面の座標	////speed_eachにXY平面の出力を付け加えたため不要かも
	relative_speed=OFF;		//重心に対する相対速度を出力するかしないか 1=ON 0=OFF
	speed_AVS=ON;			//microAVSによる3D速度分布出力するかしないか 1=ON 0=OFF
	legend_speed=0.1;		//speed.datの凡例に出力する速度[m/s]
	set_zero_speed=ON;		//restart時に速度をゼロセットするかしないか  1=ON 0=OFF
	
	///////GPU関係
	M_form=ELL;				//CG法における係数行列の格納方法 CSR_scl,CSR_vec,ELLの3つをｻﾎﾟｰﾄ
	MAX_thread=512;			//ひとつのSMあたりの最大ｽﾚｯﾄﾞ数　ふつうは512

	///////ﾌｧｲﾙ出力変数 
	interval=200;//20;			//２以上の整数にすること
	AVS=6;                  //0:普通　1:圧力　2:温度 3:壁非表示 4：表面のみ 5:壁(in,outの区別つき) 6:特定
	maxT=343;
	minT=293;
	P_interval=10000;
	ST_interval=10000;
	PND_interval=10000;
	mesh_output_interval=1;
	plot_nearbywall=OFF;	//壁に近い流体粒子を出力するかどうか。三次元CC解析でるつぼを割った中身の流体を見たいときに用いる
}