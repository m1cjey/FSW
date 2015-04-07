////TetGen呼び出し関数  TetGenの入り口

#include "header.h"

#define FULL 3
#define REMESH 4
#define FULL_INPORT 5


//TetGen入口  モデルにより分岐
void tetgen_function::call_TetGen(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	cout<<"TetGenによるメッシュ生成開始"<<endl;
	clock_t t1=clock();	//クロック数取得


	//TetGen_test(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//test
	if(CON.get_mesh_input()==0)
	{
		//////モデル毎に分岐/////////////////////////////
		if(CON.get_model_number()==14)
		{
			if(CON.get_ea_model()==0) TetGen_nanoe(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//静電霧化用
			else if(CON.get_ea_model()==1) TetGen_nanoe2(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//静電霧化用,cad読み込み
		}
		if(CON.get_model_number()==15)	TetGen_rayleigh(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//レイリー分裂用
		if(CON.get_model_number()==20 || CON.get_model_number()==21)	TetGen_CC(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//コールドクルーシブル用
		if(CON.get_model_number()==25)	TetGen_eddytest(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	
		if(CON.get_model_number()==30)	TetGen_nanoe2(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//静電霧化用,cad読み込み
		/////////////////////////////////////////////////
	}
	if(CON.get_mesh_input()==1)
	{
		if(CON.get_model_number()==20 || CON.get_model_number()==21) TetGen_inport(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
	}

	clock_t t2=clock();	//クロック数取得
	cout<<"メッシュ生成完了  CPU time="<<(t2-t1)/CLOCKS_PER_SEC<<endl;
}


//////これより下に、作りたいモデルのメッシュ生成プログラムを記述してください/////////////////////////////////////////////

//静電霧化用  メッシュ生成
void tetgen_function::TetGen_nanoe(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenioクラス宣言と初期化
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////境界節点と境界面のデータの作成////////////////
	/*-----------------------------------------------------------
	Set******Boundary関数にて，部品ごとのNODEとFACEに境界節点と境界面を格納し，
	AddBoundaryData関数にて，NODEとFACEのデータをNODEallとFACEallに格納していく．
	-----------------------------------------------------------*/

	////////////水滴領域////////////
	NODE.clear();
	FACE.clear();
	//SetFluidBoundary(CON, PART, TET, NODE, FACE);	//流体境界の設定
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//外側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//内側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetTRANS(NODE, TRANS);	//節点番号と粒子番号の対応付け
	//AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////空気領域////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////水滴表面付近追加節点////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	
	////////////平板電極////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////円柱電極////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////電極土台////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	


	//.nodeファイル作成
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"polyファイル読み込み"<<endl;
	in.load_poly("all_boundary");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	cout<<"自動分割開始"<<endl;
	tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1以下では切れない
	out.save_nodes("output");
	out.save_elements("output");

	//材質修正
	cout<<"材質修正"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集

	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//静電霧化用  メッシュ生成// cadから読み込み
void tetgen_function::TetGen_nanoe2(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenioクラス宣言と初期化
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////境界節点と境界面のデータの作成////////////////
	/*-----------------------------------------------------------
	Set******Boundary関数にて，部品ごとのNODEとFACEに境界節点と境界面を格納し，
	AddBoundaryData関数にて，NODEとFACEのデータをNODEallとFACEallに格納していく．
	-----------------------------------------------------------*/

	/*///////////水滴領域////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//流体境界の設定
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//外側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//内側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	SetTRANS(NODE, TRANS);	//節点番号と粒子番号の対応付け
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	///////////*/
	
	////////////空気領域////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////
	/*
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	////////*/

	////////////対向電極////////////
	NODE.clear();
	FACE.clear();
	//SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	SetCounterElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);
	

	////////////放電電極////////////
	NODE.clear();
	FACE.clear();
	SetDischargeElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	/*///////////電極土台////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	//////*/
	


	//.nodeファイル作成
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"polyファイル読み込み"<<endl;
	in.load_poly("all_boundary");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	cout<<"自動分割開始"<<endl;
	tetrahedralize("pqa1.67e-7AYn", &in, &out);	//1.1以下では切れない
	//tetrahedralize("pq1.2An", &in, &out);	//1.1以下では切れない
	out.save_nodes("output");
	out.save_elements("output");

	//材質修正
	cout<<"材質修正"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集

	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}


//ﾚｲﾘｰ分裂用  メッシュ生成
void tetgen_function::TetGen_rayleigh(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenioクラス宣言と初期化
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////水滴領域////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);
	SetTRANS(NODE, TRANS);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////空気領域////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////空気領域追加節点////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	//*/


	//.nodeファイル作成
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"polyファイル読み込み"<<endl;
	in.load_poly("all_boundary");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	cout<<"自動分割開始"<<endl;
	tetrahedralize("pq1.3a1.67e-7AYYn", &in, &out);	//1.1未満では切れない デフォルトはrqa1.1AYYn
	out.save_nodes("output");
	out.save_elements("output");
	//out.save_faces("output");
	//out.save_neighbors("output");
	//*/

	cout<<"材質修正"<<endl;
	//ModifyAttribute(CON, TET, NODEall, ELEMall, in, out);
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集

	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//CC用  メッシュ生成
void tetgen_function::TetGen_CC(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenioクラス宣言と初期化
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////境界節点と境界面のデータの作成////////////////
	/*-----------------------------------------------------------
	Set******Boundary関数にて，部品ごとのNODEとFACEに境界節点と境界面を格納し，
	AddBoundaryData関数にて，NODEとFACEのデータをNODEallとFACEallに格納していく．
	-----------------------------------------------------------*/

	////////////水滴領域////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//流体境界の設定
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//外側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//内側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	SetTRANS(NODE, TRANS);	//節点番号と粒子番号の対応付け
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////空気領域////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////コイル////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	/*///////////るつぼ////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////水滴表面付近追加節点////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////平板電極////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////円柱電極////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////電極土台////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.nodeファイル作成
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"polyファイル読み込み"<<endl;
	in.load_poly("all_boundary");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	cout<<"自動分割開始"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1以下では切れない
	tetrahedralize("pq1.2a1.67e-6AYYn", &in, &out);	//1.1以下では切れない
	out.save_nodes("output");
	out.save_elements("output");

	//材質修正
	cout<<"材質修正"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集
			
	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//  渦電流問題用メッシュ生成
void tetgen_function::TetGen_eddytest(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenioクラス宣言と初期化
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////境界節点と境界面のデータの作成////////////////
	/*-----------------------------------------------------------
	Set******Boundary関数にて，部品ごとのNODEとFACEに境界節点と境界面を格納し，
	AddBoundaryData関数にて，NODEとFACEのデータをNODEallとFACEallに格納していく．
	-----------------------------------------------------------*/

	////////////水滴領域//////////// この問題では水滴はいらない
	NODE.clear();
	FACE.clear();
	//SetFluidBoundary(CON, PART, TET, NODE, FACE);	//流体境界の設定
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//外側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//内側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetTRANS(NODE, TRANS);	//節点番号と粒子番号の対応付け
	//AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	/*///////////空気領域////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	//////*/

	////////////空気領域2////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary2(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	/*///////////コイル////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	////////////導電体部////////////
	NODE.clear();
	FACE.clear();
	SetConductBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CONDUCT);
	////////*/


	/*///////////るつぼ////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////水滴表面付近追加節点////////////
	NODE.clear();
	FACE.clear();
	//SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	//AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////平板電極////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////円柱電極////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////電極土台////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.nodeファイル作成
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"polyファイル読み込み"<<endl;
	in.load_poly("all_boundary");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	cout<<"自動分割開始"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1以下では切れない
	//tetrahedralize("pq1.2a1.67e-6AYYn", &in, &out);	//1.1以下では切れない
	tetrahedralize("pq1.5Aa1e-3YYn", &in, &out);	//1.1以下では切れない
	out.save_nodes("output");
	out.save_elements("output");

	//材質修正
	cout<<"材質修正"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集
			
	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//外部作成データから読み込んでメッシュ作成
void tetgen_function::TetGen_inport(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	cout<<"magnetデータを読み込みtetgenメッシュ作成"<<endl;
	//tetgenioクラス宣言と初期化
	tetgenio in, out,addin;
	in.initialize();
	out.initialize();

	//tetgenio addin;
	addin.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////境界節点と境界面のデータの作成////////////////
	/*-----------------------------------------------------------
	Set******Boundary関数にて，部品ごとのNODEとFACEに境界節点と境界面を格納し，
	AddBoundaryData関数にて，NODEとFACEのデータをNODEallとFACEallに格納していく．
	-----------------------------------------------------------*/

	/*///////////水滴領域////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//流体境界の設定
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//外側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//内側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	SetTRANS(NODE, TRANS);	//節点番号と粒子番号の対応付け
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////空気領域////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////コイル////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	/*///////////るつぼ////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////水滴表面付近追加節点////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////平板電極////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////円柱電極////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////電極土台////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.nodeファイル作成
	//MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	//MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"staticファイル読み込み"<<endl;
	//in.load_node("static");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	in.load_tetmesh("static");	//.node,.eleを読込
	//addin.load_node("fluid");
	cout<<"読込メッシュの再構成開始"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1以下では切れない
	tetrahedralize("r", &in, &out);	//読み込んだstaticからtetgen用に再構成
	out.save_nodes("output");
	out.save_elements("output");
	out.save_faces("output");


	////////////水滴領域////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//流体境界の設定
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//外側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//内側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	SetTRANS(NODE, TRANS);	//節点番号と粒子番号の対応付け
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	//in.load_node("fluid");	//.node,.eleを読込
	cout<<"流体節点の追加開始"<<endl;
	tetrahedralize("i", &in, &out);
	

	//材質修正
	//cout<<"材質修正"<<endl;
	//ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集
	
	
			
	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//CC用  メッシュ生成
void tetgen_function::TetGen_test(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenioクラス宣言と初期化
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////境界節点と境界面のデータの作成////////////////
	/*-----------------------------------------------------------
	Set******Boundary関数にて，部品ごとのNODEとFACEに境界節点と境界面を格納し，
	AddBoundaryData関数にて，NODEとFACEのデータをNODEallとFACEallに格納していく．
	-----------------------------------------------------------*/

	////////////水滴領域////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//流体境界の設定
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//外側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//内側ダミー節点を使った水滴メッシュ生成（バグが取れずに未完成）
	SetTRANS(NODE, TRANS);	//節点番号と粒子番号の対応付け
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////空気領域////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////コイル////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	/*///////////るつぼ////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////水滴表面付近追加節点////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////平板電極////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////円柱電極////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////電極土台////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.nodeファイル作成
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"polyファイル読み込み"<<endl;
	in.load_poly("all_boundary");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	cout<<"自動分割開始"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1以下では切れない
	tetrahedralize("pq1.2a1.67e-6AYYn", &in, &out);	//1.1以下では切れない
	out.save_nodes("output");
	out.save_elements("output");

	//材質修正
	cout<<"材質修正"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集
			
	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}