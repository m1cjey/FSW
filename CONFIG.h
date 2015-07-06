////新しい粒子を追加したときいじらないといけないのは①u_laplacian_fのﾉﾝｽﾘｯﾌﾟ②P_gradient③P_gradient2④calc_Temperatureのk[i]
////⑤ freeon
#ifndef config
#define config

#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。

class mpsconfig
{       
	///////解析条件

	double dt;//時間刻み幅 
	int    step;
	int    dimention;//次元
	double maxX;//解析領域
	double minX;
	double maxY;
	double minY;
	double maxZ;
	double minZ;

	int material;
	
	///////流体1の物性値
	
	double density;//密度
	double nensei;//粘性係数
	double sigma;//表面張力係数
	double Cp;//定圧比熱
	double k;//熱伝導率
	double latent_H;//潜熱(latent_heat)
	double MP;//融点(melting point)
	double CTE;//熱膨張係数[1/K] 

	int temperature_depend;	//物性値が温度依存の変化をするかどうか

	///////流体2の物性値
	
	double density2;//密度
	double nensei2;//粘性係数
	double sigma2;//表面張力係数
	double Cp2;//定圧比熱
	double k2;//熱伝導率
	double latent_H2;//潜熱(latent_heat)
	double MP2;//融点(melting point)
	double CTE2;//熱膨張係数[1/K] 
	
	////////粒子配置用
	
	int fluidwidth;//流体の代表長さ
	int fluid_h;//流体の代表長さ
	double distancebp;//初期粒子間距離
	double wlength;//左右の壁の距離。流体の何倍か
	double height;//流体の重心高さ
	double tool_angle;//FSWにおいて、ツールを傾ける角度
	int tool_type;
	int process_type;
	double dwelling_time;
	int airwall;
	int change_step;
	
	///////粒子法用パラメータ
	
	double re;		//一般的な粒子半径re
	double re2;		//ラブラシアン用のre
	double re3;		//表面張力用のre
	double re4;		//freeon
	double beta;		//β
	int    dx;		//index用格子幅。最大粒子半径を含む整数にすること
	double times;		//ﾍﾞｸﾄﾙﾌﾟﾛｯﾄ用の倍率
	
	////////表面張力関係
	
	int    surface_tension; //表面張力スイッチ　0=OFF 1=教科書 2=粒子間引力
	int    smooth;			//スムージングスイッチ　1=ON 0=OFF
	int    SFT;				//1=ON 0=OFF 表面張力(surface tension)の温度(T)依存性スイッチ
	int    smn;				//ｽﾑｰｼﾞﾝｸﾞの回数
	int    smth_sumP;       //for ver3. curvのｽﾑｰｼﾞﾝｸﾞ 0=OFF 1=ON
    int    wall_direct;		//for ver2. BDWALLの法線ﾍﾞｸﾄﾙ計算方法　0=無視 1=通常 2=水平方向 3=垂直方向
	int    non_slip;		//for ver.5 粘着条件 ONならBDWALL上に仮想的に液体が存在する設定。その際はBDWALLの配置に留意
	int    suf_modify;		//曲率計算で得られた曲面に一致するように粒子位置を修正するか、しないか
	int    interpolate_curv;//曲率が計算されなかった粒子の曲率を、周辺粒子の値から補間
	double Cst;				//ポテンシャル係数
	int    C_type;			//ポテンシャル係数の決定方法
	double C_times;			//ポテンシャル係数の倍数
	double wall_C;			//壁との親和力係数
	
	////////圧力関係

	int iteration_count;//圧力計算を行う回数 通常は1に設定
	int    solution;	//連立方程式の解放　0=CG method 1=ICCG method
	int    B_of_P;  	//圧力計算の解行列 0=教科書 1=速度発散 2=0+1
	double w_div;		//発散の重み 
	int    divU_method;	//速度発散の計算方法 1=MPS,2=WLSM,3=自身を通らないWLSM 4=入部ら
	int    divU_order;	//WSLMでdivUを計算するときのオーダー MAX=3
	int    HL_sw;		//ﾗﾌﾟﾗｼｱﾝに対し高次の離散化を行う 1=ON 0=OFF
	int	   dir_for_P;	//圧力計算の際、表面に非ゼロのディリクレ値を代入するか、しないか
	int    initialP;	//圧力計算時に初期値として現圧力を代入するかどうか　1=ON 0=OFF
	int    P_twice;
	double CGep;		//圧力計算時における収束判定
	int	   omp_P;		//圧力計算時のCG法でopenMPを使用するか 1=ON 0=OFF
	int    negativeP;	//負圧　1=許可 0=不許可
	int    set_P;		//圧力強制決定スイッチ 1=ON 0=OFF
	int    Pgrad;		//圧力勾配計算法 0=教科書 1=表面のみ法線
	int    minP;		////最少圧力ｽｲｯﾁ 1=ON 0=OFF
	int    ave_Pdirect;	//圧力勾配用法線ﾍﾞｸﾄﾙの平均化 1=ON 0=OFF
	int    Pgrad_order;	//圧力勾配ver.4における精度 1=線形 2=2次
	int	   artP_sw;		//人工圧力計算フラグ 0=OFF 1=artP等方 2=Monaghan
	double artP;		//人工圧力値[Pa]
	double Pgrad_times;	//圧力勾配をﾌﾟﾛｯﾄする際の、表面張力に対する倍率 通常は1に設定
	int    P_AVS;		////microAVS用の圧力ファイルを出力するstep間隔。0ならOFF
	double Pgrad_limit;
	int    gridless_P;	//グリッドレス法による圧力計算を行うか(ON)、MPSによる計算を行うか(OFF)
	int	   wall_poly;
	int    interpolate_surface_u; //圧力勾配による速度修正後の表面粒子の速度を内部粒子の速度で補完するかしないか。

	////////温度場関係
	
	int    T_field;		//温度場解析　　1=ON 0=OFF
	int    insulate;	//壁との断熱状態 0=断熱　1=非断熱
	int    T_laplacian;	//温度のﾗﾌﾟﾗｼｱﾝ。　0=教科書　1=λ[i] 2=発散・勾配
	double wall_density;		//壁の密度[kg/m^3]
	double wall_Cp;			//壁の比熱[J/(kg・K)]
	double wall_k;				//壁の熱伝導率[W/(m・K)]
	double    initialT;		//初期温度
	double    roomT;		//室温 int型でよい
	int    air_cool;	//空気との熱伝達を考慮するかしないか 1=ON 0=OFF
	int    T_expansion;	//熱膨張計算　1=ON 0=OFF
	int    T_CG;
	int    T_solution;
	double T_CGep;
	int    buoyant;		//浮力(密度のﾌﾞｼﾞﾈｽｸ近似)  1=ON 0=OFF
	double TplotZ;		//3D解析において、XY平面の温度を出力するときのＺ座標
	int    T_AVS;		//microAVS用の温度ファイルを出力するstep間隔。0ならOFF
	int output_temperature_face;
	int output_equivalent_strain_rate_face;
	int output_forward;
	int output_backward;
	int output_another_face;

	////////電磁場の解法
	int		EM_method;			//電磁場の解法 0=OFF 1=FEM 2=BEM 3=磁気ﾓｰﾒﾝﾄ法
	int		EM_calc_type;		//0=OFF 1=電場 2=磁場 3=磁場(渦電流)
	int		EM_interval;    //電磁場計算を何ｽﾃｯﾌﾟに一回行うか。通常は1に設定
	double  EM_distance;	//電磁力計算を行う頻度を、粒子の総移動距離で判定する 0:OFF この値*leを移動距離が超えたらFEMを実行
	int    region_shape;	//解析領域形状　0=立方体 1=円筒
	double XL;				//解析領域左端
	double XR;              //解析領域右端
	double YU;              //解析領域上端
	double YD;              //解析領域下端
	double ZU;				//解析領域Z方向上端
	double ZD;				//解析領域Z方向下端
	double RU;				//解析領域が円筒形となるときのその半径

	////////メッシュ
	int    mesher;			//メッシュ作成方法 0:吉川さんのデローニ分割プログラム
	int    mesh_input;		//meshの作成方法 0:自分 1:Magnet
	int    remesh_sw;		//ON:remesh領域を考慮し、そこだけremesh OFF:FULL領域を毎回分割
	double boxalpha;		//ｽｰﾊﾟｰﾎﾞｯｸｽの大きさ 1～10の実数にすること。通常は2に設定
	int    fine;			//節点の自動追加を行うかどうか 0:OFF 1:ON
	double co_fine;			//節点自動追加における係数.大きいほど物体から離れたら粗くなる
	int    add_points;		//自動節点追加数
	int    poly_memory;		//poly3D()で確保するメモリの数
	int	   CDT_sw;
	int	   defer_f;
	int	   mesh_smth;
	double material_judge;
	double del_length;
	int	   BD_fluid;
	int    output_wall_way;	//壁粒子を出力してメッシュを作成するかどうか　0=流体のみ 1=流体と壁すべて　2=流体とINWALLのみ
	double layer_thin;			//空気層の厚さ distancebp*layer_thinが実際の厚さ
	int    layer_num;			//空気層の層数

	////////FEM関係
	int	   BEM_elm_type;		//要素ﾀｲﾌﾟ 0:節点要素 1:辺要素
	int	   FEM_elm_type;	//要素ﾀｲﾌﾟ 0:節点要素 1:辺要素
	int    FEM_smn;			//電磁力ｽﾑｰｼﾞﾝｸﾞ回数　0ならOFF
	int    max_DN;			//ﾃﾞｨﾘｸﾚ型境界条件をとる最大節点数
	int	   FEMCG;			//FEMにおける行列解法 0:CG 1:ICCG
	double CGaccl;			//CG,ICCG法における加速ファクタ　CGaccelerator
	double FEMCGep;			//ICCGの収束判定
	double FEMtimes;		//電磁力をﾌﾟﾛｯﾄする際の、表面張力に対する倍率 通常は1に設定
	double legend_F;		//F.datの凡例に出力する力[N]
	int	   node_sort;
	int	   thinout_fluid;
	
	///////電界計算
	double V;		//電圧
	int    V_step;		//電圧印加開始ステップ
	double r_perm;		//比誘電率 (relative permittivity)
	int    V_con;		//電圧条件　0:ﾊﾟﾙｽ　1:ﾘﾆｱ　2:時定数
	double initial_V;	//電圧条件ﾘﾆｱのときの初期電圧
	double E_times;		//ﾌｧｲﾙ出力する際の、電界の倍率
	int    eleforce;	//静電力計算方法 1:節点力法 2:積分面
	int	   charge;		//電荷を考慮するかしないか。0=OFF 1=ON
	
	///////磁場計算
	
	int    J_input_way;	//電流密度入手方法 0:自分 1:ソフト
	double J0;//強制電流密度[A/m2]
	double I0;//強制電流値[A]
	double RP;//比透磁率(Relative Permeability)
	double ele_conduc;//電気伝導率
	double ele_conduc2;//電気伝導率
	double  Hz;///交流の周波数
	int    div_Hz;//１周期の分割数(解析精度)
	int    jheat;		//渦電流による発熱を考慮するか　0=OFF 1=ON
	int	   m_force;		//電磁力計算方式 0=ﾏｸｽｳｪﾙの応力 1=体積力
	int    NLBHI;		//体積力において、要素Ｂから要素Ｈを求める際に非線形性を考慮するか、しないか(non linier B H inverter)
	int    NLMH;		//Ｍの算出に非線形性を考慮するか、しないか
	double magnet_H;	//永久磁石の大きさ
	double magnet_B;	//永久磁石の強さ[T]
	int    uniform_B_sw;	//解析領域中に一様磁場を発生させるか否か 0=OFF 1=ON
	double uniform_B;	//一様磁場の大きさ[T]
	double B_times;		//ﾌｧｲﾙ出力する際の、磁束密度の倍率
	int    plot_B_type;	//磁束密度出力タイプ 1=ﾍﾞｸﾄﾙ 2=スカラー
	int	   Je_crucible; //渦電流の計算対象　0=溶融金属 1=溶融金属とるつぼ（コイルも）
	int	   J0eqJe;
	int    m_A;
	int    jw_interval;
	int    jw_Faverage;
	int    static_dirichlet;
	double eddy_times;
	int    A_phi;
	int    parabolic_node_element;
	
	////////各種スイッチ　どのような計算を考慮するか、しないか
	
	double g;		//重力加速度　通常は-9.8とすること
	int    restart;		//restart用スイッチ
	int    autosave;    //ｵｰﾄｾｰﾌﾞ間隔。無効にしたいときは大きな数字を代入しておく
	double curan;		//クーラン数条件 0ならOFF
	int    modify_position;//粒子間距離がleより小さい場合これを修正するか、しないか
	int    vis_calc_type;//粘性項計算手法　0=陽解法 1=陰解法
	int    wall_adheision;	//壁の粘性状態　1=ﾉﾝｽﾘｯﾌﾟ 0=ﾌﾘｰｽﾘｯﾌﾟ
	int    laplacian;	//ﾗﾌﾟﾗｼｱﾝ。　0=教科書　1=λ[i] 2=発散・勾配
	int	   vis_solver;  //粘性項を陰解析で解く際の行列ｿﾙﾊﾞｰ 0:CG 1:ICCG
	int    initial_u;	//粘性項を陰的に解く際に、初期値として現在速度を入力するか、しないか
	int    temporary_r; //陽解析後に仮の位置を計算するかしないか。 1=ON 0=OFF 通常はONに設定
	
	int    fix_center;	//1=ON 0=OFF
	int    freeon;		//粒子依存関係関数　1:並列化可能 2:並列不可
	int    freeon3sw;	//freeon3を計算するかしないか 1=ON 0=OFF
	int    surface_judge2;//surface_judge2()を使用するかしないか  1=ON 0=OFF
	int    move_prtcl;	//移動粒子を考慮するかしないか 1=ON 0=OFF
	int    move_u_dirct;//移動粒子を移動する方向　現在は±X方向=±1,±Y方向=±2,±Z方向=±3
	double move_speed;	//移動粒子の移動速度[m/s]
	double move_speed2; //移動粒子の移動速度[m/s]
	int	   check_region;//条件を満たした粒子を削除するかしないか
	double max_speed;	//chec_regionがONになっているとき、粒子の速度がこの値以上なら削除する
	int    check_something;//check_something()を実行するかしないか 1=ON 0=OFF
	int	   AVS_HALF;
	int    adaptive_sw;	//解像度を可変にするか、しないか 1=ON 0=OFF
	double threshold;	//アダプティブにする際の、圧力誤差の閾値
	int    fix_surface;
	int output_viscousity_face;
	
	int    model_number;
	int    model_set_way;		//modelをセットする方法　0=正方格子 1=MD
	int    model_inherit;		//前の解析で使った粒子モデルを引き継ぐ
	int	   ea_model;

	///////速度ﾌﾟﾛｯﾄ変数
	int    speed_plot_particle;	//速度をﾌﾟﾛｯﾄする粒子の種類 1=すべて 2=fluid
	double speedtimes;			//速度ﾌﾟﾛｯﾄ時の、座標に対する速度の倍率
	int    speed_face;			//speed.datの出力面 0=YZ平面 1=XZ 2=XY
	double speed_face_p;		//3D解析時のspeed.datの出力面の座標
	int    ax_sym_modify;		//3D時のspeed.datに関して、軸対称による出力修正を行うか否か　1=ON 0=OF
	int    flat_speed_plot;		//水平方向の速度をﾌﾟﾛｯﾄするかしないか
	double flat_speed_p;		//flat_speed.datの出力面の座標
	int	   relative_speed;		//重心に対する相対速度を出力するかしないか 1=ON 0=OFF
	int    speed_AVS;			//microAVSによる3D速度分布出力するかしないか 1=ON 0=OFF
	double legend_speed;		//speed.datの凡例に出力する速度[m/s]
	int    set_zero_speed;//restart時に速度をゼロセットするかしないか  1=ON 0=OFF

	///////GPU関係
	int    M_form;				//CG法における係数行列の格納方法 CSR_scl,CSR_vec,ELLの3つをｻﾎﾟｰﾄ
	int    MAX_thread;			//ひとつのSMあたりの最大ｽﾚｯﾄﾞ数　ふつうは512
	
	///////ﾌｧｲﾙ出力変数
	
	int    interval;//AVS用ステップ間隔
	int    AVS;	//0:普通　1:圧力　2:温度
	double maxT;    //AVS(2)における最大温度
	double minT;    //AVS(2)における最小温度
	int P_interval;
	int ST_interval;
	int PND_interval;
	int mesh_output_interval;
	int plot_nearbywall;
	
	public:
	mpsconfig();
	double get_dt()		{return dt;}
	int get_step()	{return step;}
	int    get_dimention()	{return dimention;}
	double get_maxX()	{return maxX;}
	double get_minX()	{return minX;}
	double get_maxY()	{return maxY;}
	double get_minY()	{return minY;}
	double get_maxZ()	{return maxZ;}
	double get_minZ()	{return minZ;}

	int get_material() {return material;}
	
	double get_density()	{return density;}
	double get_nensei()	{return nensei;}
	double get_sigma()      {return sigma;}
	double get_vis()	{return nensei/density;}//動粘性係数
	double get_Cp() 	{return Cp;}
	double get_k()		{return k;}
	double get_latent_H()	{return latent_H;}
	double get_MP()		{return MP;}
	double get_CTE()	{return CTE;}

	int get_temperature_depend(){return temperature_depend;}
	double get_density2()	{return density2;}
	double get_nensei2()	{return nensei2;}
	double get_sigma2()      {return sigma2;}
	double get_vis2()	{return nensei2/density2;}//動粘性係数
	double get_Cp2() 	{return Cp2;}
	double get_k2()		{return k2;}
	double get_latent_H2()	{return latent_H2;}
	double get_MP2()		{return MP2;}
	double get_CTE2()	{return CTE2;}
	
	int get_fluidwidth() {return fluidwidth;}
	int get_fluid_h() {return fluid_h;}
	double get_distancebp() {return distancebp;}
	double get_wlength()    {return wlength;}
	double get_height()	{return height;}
	double get_tool_angle() {return tool_angle;}
	int    get_tool_type() {return tool_type;}
	int	get_process_type()	{return process_type;}
	double get_dwelling_time(){return dwelling_time;}
	int    get_airwall() {return airwall;}
	int get_change_step(){return change_step;}
	
	double get_re()		{return re;}
	double get_re2()	{return re2;}
	double get_re3()	{return re3;}
	double get_re4()	{return re4;}	
	double get_beta()	{return beta;}
	int    get_dx()		{return dx;}
	double get_times()	{return times;}
	int get_X_mesh()	{return (int)((maxX-minX)/(distancebp*dx)+0.00000000000001);}//X軸方向の格子数 丸め誤差を防ぐために0.001を足している
	int get_Y_mesh()	{return (int)((maxY-minY)/(distancebp*dx)+0.00000000000001);}//Y軸方向の格子数 丸め誤差を防ぐために0.001を足している
	int get_Z_mesh()	{return (int)((maxZ-minZ)/(distancebp*dx)+0.00000000000001);}//Z軸方向の格子数 丸め誤差を防ぐために0.001を足している
	int get_number_of_mesh(){return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);}//格子数：X_mesh*Y_mesh*Z_mesh
	
	int    get_surface_tension()	{return surface_tension;}
	int    get_smooth()	{return smooth;}
	int    get_SFT()	{return SFT;}
	int    get_smn()	{return smn;}
	int    get_smth_sumP()  {return smth_sumP;}
	int    get_wall_direct(){return wall_direct;}
	int    get_non_slip()	{return non_slip;}
	int    get_suf_modify(){return suf_modify;}
	int    get_interpolate_curv() {return interpolate_curv;}
	double get_Cst() {return Cst;}				//ポテンシャル係数
	void   set_Cst(double calc_Cst) {Cst=calc_Cst;}	//ポテンシャル係数
	int    get_C_type() {return C_type;}			//ポテンシャル係数の決定方法
	double get_C_times() {return C_times;}			//ポテンシャル係数の倍数
	double get_wall_C() {return wall_C;}			//壁との親和力係数

	int    get_iteration_count(){return iteration_count;}
	int    get_Pgrad()	{return Pgrad;}
	int    get_ave_Pdirect(){return ave_Pdirect;}
	int    get_solution()	{return solution;}
	int    get_B_of_P()	{return B_of_P;}
	double get_w_div()	{return w_div;}
	int    get_divU_method() {return divU_method;}//速度発散の計算方法 1=MPS,2=WLSM,3=自身を通らないWLSM 4=入部ら
	int    get_divU_order()	{return divU_order;}//WSLMでdivUを計算するときのオーダー MAX=3
	int    get_HL_sw()	{return HL_sw;}
	int    get_dir_for_P(){return dir_for_P;}
	int    get_initialP()	{return initialP;}
	int    get_P_twice(){return P_twice;}
	double get_CGep()	{return CGep;}
	int    get_omp_P() {return omp_P;}
	int    get_negativeP(){return negativeP;}
	int    get_minP()	{return minP;}
	int    get_set_P()	{return set_P;}
	int    get_artP_sw(){return artP_sw;}
	double get_artP()	{return artP;}
	double get_Pgrad_times(){return Pgrad_times;}
	int get_Pgrad_order(){return Pgrad_order;}
	int    get_P_AVS()	{return P_AVS;}
	double get_Pgrad_limit(){return Pgrad_limit;}
	int    get_gridless_P()	{return gridless_P;}
	int    get_wall_poly()	{return wall_poly;}
	int    get_interpolate_surface_u(){return interpolate_surface_u;}
	
	int    get_T_field()	{return T_field;}
	int    get_insulate()	{return insulate;}
	int    get_T_laplacian(){return T_laplacian;}
	double get_wall_density(){return wall_density;}
	double get_wall_Cp()    {return wall_Cp;}
	double get_wall_k()		{return wall_k;}
	double get_initialT()	{return initialT;}
	double get_roomT()	{return roomT;}
	int    get_air_cool() {return air_cool;}
	int    get_T_expansion(){return T_expansion;}
	int    get_T_solution(){return T_solution;}
    int    get_T_CG(){return T_CG;}
	double get_T_CGep(){return T_CGep;}
	int    get_buoyant()	{return buoyant;}
	double get_TplotZ()		{return TplotZ;}
	int    get_T_AVS()		{return T_AVS;}
	int get_output_temperature_face()	{return		output_temperature_face;}

	int    get_EM_method()	{return EM_method;}
	int    get_EM_calc_type() {return EM_calc_type;}
	int    get_EM_interval(){return EM_interval;}
	double get_EM_distance(){return EM_distance;}
	int    get_region_shape() {return region_shape;}
	double get_XL()		{return XL;}
	double get_XR()		{return XR;}
	double get_YU()		{return YU;}
	double get_YD()		{return YD;}
	double get_ZU()		{return ZU;}
	double get_ZD()		{return ZD;}
	double get_RU()		{return RU;}

	int    get_mesher() {return mesher;}
	int    get_mesh_input() {return mesh_input;}
	int    get_remesh_sw()	{return remesh_sw;}
	double get_boxalpha(){return boxalpha;}
	int    get_fine()	{return fine;}
	double get_co_fine(){return co_fine;}
	int    get_add_points()	{return add_points;}
	int    get_poly_memory(){return poly_memory;}
	int	   get_CDT_sw(){return CDT_sw;}
	int    get_defer_f(){return defer_f;}
	int    get_mesh_smth(){return mesh_smth;}
	double get_material_judge(){return material_judge;}
	double get_del_length(){return del_length;}
	int	   get_BD_fluid(){return BD_fluid;}
	int    get_output_wall_way()	{return output_wall_way;}
	double get_layer_thin()		{return layer_thin;}	//空気層の厚さ distancebp*layer_thinが実際の厚さ
	int    get_layer_num()				{return layer_num;}//空気層の層数

	int	   get_BEM_elm_type(){return BEM_elm_type;}
	int    get_FEM_elm_type() {return FEM_elm_type;}
	int    get_FEM_smn()	{return FEM_smn;}
	int    get_max_DN()	{return max_DN;}
	int    get_FEMCG()	{return FEMCG;}
	double    get_CGaccl()	{return CGaccl;}
	double get_FEMCGep() {return FEMCGep;}
	double get_FEMtimes(){return FEMtimes;}
	double get_legend_F(){return legend_F;}
	int	   get_node_sort(){return node_sort;}
	int	   get_thinout_fluid(){return thinout_fluid;}
	
	double get_V()		{return V;}
	int    get_V_step()	{return V_step;}
	double get_r_perm()	{return r_perm;}
	int    get_V_con()	{return V_con;}
	double get_initial_V()	{return initial_V;}
	double get_E_times(){return E_times;}
	int	   get_eleforce(){return eleforce;}
	int    get_charge()	{return charge;}
	
	int    get_J_input_way() {return J_input_way;}
	double get_J0()		{return J0;}
	double get_I0()		{return I0;}
	double get_RP()		{return RP;}
	double get_ele_conduc() {return ele_conduc;}
	double get_ele_conduc2() {return ele_conduc2;}
	double    get_Hz()		{return Hz;}
	int    get_div_Hz()	{return div_Hz;}
	int    get_jheat()	{return jheat;}
	int    get_m_force(){return m_force;}
	int	   get_NLBHI()	{return NLBHI;}
	int    get_NLMH()	{return NLMH;}
	double get_magnet_H(){return magnet_H;}
	double get_magnet_B(){return magnet_B;}
	int    get_uniform_B_sw(){return uniform_B_sw;}
	double get_uniform_B(){return uniform_B;}
	double get_B_times(){return B_times;}
	int    get_plot_B_type(){return plot_B_type;}
	int	   get_Je_crucible(){return Je_crucible;}
	int    get_J0eqJe(){return J0eqJe;}
	int    get_m_A(){return m_A;}
	int    get_jw_interval(){return jw_interval;}
	int    get_jw_Faverage(){return jw_Faverage;}
	int    get_static_dirichlet(){return static_dirichlet;}
	double get_eddy_times()  {return eddy_times;}
	int		get_A_phi() {return A_phi;}
	int    get_parabolic_node_element() {return parabolic_node_element;}

	double get_g()		{return g;}
	int    get_restart()    {return restart;}
	int    get_autosave()	{return autosave;}
	double get_curan()    {return curan;}
	int    get_modify_position() {return modify_position;}
	int    get_vis_calc_type(){return vis_calc_type;}
	int    get_wall_adheision() {return wall_adheision;}
	int    get_laplacian()	{return laplacian;}
	int    get_vis_solver() {return vis_solver;}
	int    get_initial_u()  {return initial_u;}
	int    get_temporary_r(){return temporary_r;}
	
	int    get_fix_center() {return fix_center;}
	int    get_freeon()	{return freeon;}
	int    get_freeon3sw(){return freeon3sw;}
	int    get_surface_judge2(){return surface_judge2;}
	int    get_move_prtcl(){return move_prtcl;}
	int    get_move_u_dirct(){return move_u_dirct;}
	double get_move_speed() {return move_speed;}
	double get_move_speed2(){return move_speed2;}
	int    get_check_something() {return check_something;}
	int    get_check_region() {return check_region;}
	double get_max_speed() {return max_speed;}
	int    get_set_zero_speed(){return set_zero_speed;}
	int    get_AVS_HALF(){return AVS_HALF;}
	int    get_adaptive_sw(){return adaptive_sw;}	//解像度を可変にするか、しないか 1=ON 0=OFF
	double get_threshold(){return threshold;}		//アダプティブにする際の、圧力誤差の閾値
	int get_fix_surface(){return fix_surface;}
	int get_output_viscousity_face(){return output_viscousity_face;}
	int get_output_equivalent_strain_rate_face(){return output_equivalent_strain_rate_face;}
	int get_output_forward(){return output_forward;}
	int get_output_backward(){return output_backward;}
	int get_output_another_face(){return output_another_face;}

	int    get_model_number(){return model_number;}
	int    get_model_set_way(){return model_set_way;}
	int    get_model_inherit(){return model_inherit;}
	int   get_ea_model(){return ea_model;}

	int    get_speed_plot_particle(){return speed_plot_particle;}
	double get_speedtimes(){return speedtimes;}
	int    get_speed_face() {return speed_face;}
	double get_speed_face_p(){return speed_face_p;}
	int    get_ax_sym_modify(){return ax_sym_modify;}
	int    get_flat_speed_plot(){return flat_speed_plot;}
	double get_flat_speed_p() {return flat_speed_p;}
	int    get_relative_speed(){return relative_speed;}
	int    get_speed_AVS(){return speed_AVS;}
	double get_legend_speed(){return legend_speed;}

	int    get_M_form()		{return M_form;}
	int    get_MAX_thread()	{return MAX_thread;}
	
	
	int    get_interval()	{return interval;}
	int    get_AVS()	{return AVS;}
	double get_maxT()	{return maxT;}
	double get_minT()	{return minT;}
	int		get_P_interval()	{return P_interval;}
	int		get_ST_interval()	{return ST_interval;}
	int		get_PND_interval()	{return PND_interval;}
	int    get_mesh_output_interval()	{return mesh_output_interval;}
	int    get_plot_nearbywall() {return plot_nearbywall;}

	double get_particle_mass()///////質量				//粒子の体積。初期配置の仕方に応じて計算方法を変える
	{
		double mass;
		if(model_set_way==0)	//正方格子のとき
		{
			if(dimention==2){mass=density*distancebp*distancebp;}
			else			{mass=density*distancebp*distancebp*distancebp;}
		}
		else if(model_set_way==1)	//細密格子のとき
		{
			if(dimention==2){mass=density*sqrt(3.0)/4*distancebp*distancebp;}
			else
			{
				//mass=density*sqrt(2.0)/12*distancebp*distancebp*distancebp;
				mass=density*distancebp*distancebp*distancebp/sqrt(2.0);
			}
		}
		else cout<<"モデルの積み方が不定です 質量を計算できません"<<endl;
		return mass;
	}

	double get_particle_volume(){ return (get_particle_mass()/get_density());}
	
};


#endif
