/*------------------------------------------------------------------------------------------------------------------------------

【自作TetGen用関数】
TetGenでメッシュ生成を行うために用いる関数は、全てtetgen_functionクラスに属しています。
まず、tetgen_functionクラスの適当なインスタンス(オブジェクト)を生成してcall_TetGenに入ってください。
以後、このクラスの関数内であれば、インスタンスの生成は不要です。

------------------------------------------------------------------------------------------------------------------------------*/

#ifndef TETFUNCH
#define TETFUNCH


class tetgen_function
{
public:

	///////////tetcall.cpp/////////////////////////////////////////////////////////////////////////////////

	void call_TetGen(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//TetGenメッシュ生成の入口
	void TetGen_nanoe(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//静電霧化用  メッシュ生成
	void TetGen_rayleigh(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//静電霧化用  メッシュ生成
	void TetGen_CC(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//コールドクルーシブル用  メッシュ生成
	void TetGen_eddytest(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//渦電流解析用テスト問題
	void TetGen_test(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//コールドクルーシブル用  メッシュ生成
	void TetGen_inport(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//コールドクルーシブル用  メッシュ生成
	void TetGen_nanoe2(mpsconfig&, vector<mpsparticle>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//静電霧化用  メッシュ生成


	///////////tetfunc.cpp/////////////////////////////////////////////////////////////////////////////////

	void GetPointList(vector<tetgen_node>&, tetgenio&, tetgenio&);
	//void GetPointList_Fluid(vector<tetgen_node>&, tetgenio&, tetgenio&);
	void GetTetrahedronList(vector<tetgen_element>&, tetgenio&, tetgenio&);
	void GetTetrahedronList_full(vector<tetgen_element>&, tetgenio&, tetgenio&);
	void GetFacetList(vector<tetgen_facet>&, tetgenio&, tetgenio&, int);
	void GetFacetListforCAD(vector<tetgen_facet>&, tetgenio&, tetgenio&, int);
	void GetFacetList_from_neigh(mpsconfig&, vector<tetgen_element>&, vector<tetgen_facet>&);
	void SelectFaceNode(mpsconfig&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void DelDummyNode(mpsconfig&, vector<tetgen_node>&, vector<tetgen_facet>&, int);

	void DelThinTetrahedron(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void DelThinTetrahedron_SharpElem(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&, vector<tetgen_facet>&, tetgenio&, tetgenio&);
	void DelTetrahedron_OutsideDummy(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void DelTetrahedron_InsideDummy(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void SetRelation_NodeNode(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void SetRelation_NodeElem(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&);
	void SetRelation_ElemElem(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&);

	void MakeNodeFile(mpsconfig&, vector<tetgen_node>&, char*);
	void MakeNodeFile_NonAttributeAndBoundary(mpsconfig&, vector<tetgen_node>&, char*);
	void MakeElemFile(mpsconfig&, vector<tetgen_element>&, char*);
	void MakeFaceFile(mpsconfig&, vector<tetgen_facet>&, char*);
	void MakePolyFile(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&, char*);
	void MakeSmeshFile(mpsconfig&, vector<tetgen_facet>&, char*);

	void SetTRANS(vector<tetgen_node>&, vector<int>&);

	void SetFluidBoundary					(mpsconfig&, vector<mpsparticle>&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetFluidBoundary_OutsideDummyNode	(mpsconfig&, vector<mpsparticle>&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetFluidBoundary_InsideDummyNode	(mpsconfig&, vector<mpsparticle>&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetAirBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetAirBoundary2(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetAirFineBoundary(mpsconfig&, vector<mpsparticle>&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetCoilBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetConductBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetCrucibleBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetPlateElectrodeBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetCounterElectrodeBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetDischargeElectrodeBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetColumnElectrodeBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetBaseBoundary(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);

	void UniteBoundaryData(mpsconfig&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<int>&);
	void AddBoundaryData(mpsconfig&, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_node>&, vector<tetgen_facet>&, int);

	void DecisionAttribute(vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void ModifyAttribute(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void ModifyAttribute_tetgenio(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void FineElement(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);

	double Distance(tetgen_node&, tetgen_node&);
	void CalcBarycentricElement(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&);

};

#endif