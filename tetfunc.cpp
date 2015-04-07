#include "header.h"


//.nodeデータ取得関数
void tetgen_function::GetPointList(vector<tetgen_node> &NODE, tetgenio &in, tetgenio &out)
{
	NODE.clear();

	tetgen_node temp;
	for(int i=0;i<out.numberofpoints;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<3;n++)	temp.r[n]=out.pointlist[i*3+n];
		temp.attribute=(int)out.pointattributelist[i];
		temp.boundary=out.pointmarkerlist[i];

		NODE.push_back(temp);
	}
}


//.eleデータ取得関数(簡易版)
void tetgen_function::GetTetrahedronList(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];
		temp.attribute=0;	//恐らく、未定義で代入すると変な数値が入ってバグるのでここでは0にしとく

		ELEM.push_back(temp);
	}
}


//.eleデータ取得関数(材質・要素要素関係含む)
void tetgen_function::GetTetrahedronList_full(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;											//節点番号
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];	//構成節点
		for(int n=0;n<4;n++)	temp.nei_elem[n]=out.neighborlist[i*4+n];	//要素-要素関係
		temp.attribute=(int)out.tetrahedronattributelist[i];				//材質
		//temp.volume=out.tetrahedronvolumelist[i];							//体積(取得できない)

		ELEM.push_back(temp);
	}
}


//.faceデータ取得関数
void tetgen_function::GetFacetList(vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out, int boundary)
{
	//※boundarymarkerは引数として与えている。PLCモードで作ったメッシュではないため境界が出力されない。

	FACE.clear();

	tetgen_facet temp;
	for(int i=0;i<out.numberoftrifaces;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<3;n++)	temp.node[n]=out.trifacelist[i*3+n];
		//temp.boundary=out.trifacemarkerlist[i];
		temp.boundary=boundary;

		FACE.push_back(temp);
	}//*/

	//for(int i=0;i<3;i++)
	//for(int n=0;n<3;n++)	cout<<out.trifacelist[i*3+n]<<endl;
}

//.faceデータ取得関数 cad(stl読み込み)用。1から始まるnode番号を
void tetgen_function::GetFacetListforCAD(vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out, int boundary)
{
	//※boundarymarkerは引数として与えている。PLCモードで作ったメッシュではないため境界が出力されない。

	FACE.clear();

	tetgen_facet temp;
	for(int i=0;i<out.numberoftrifaces;i++)
	{
		temp.id=i+out.firstnumber-1;
		for(int n=0;n<3;n++)	temp.node[n]=out.trifacelist[i*3+n]-1;
		//temp.boundary=out.trifacemarkerlist[i];
		temp.boundary=boundary;

		FACE.push_back(temp);
	}//*/

	//for(int i=0;i<3;i++)
	//for(int n=0;n<3;n++)	cout<<out.trifacelist[i*3+n]<<endl;
}


//.nodeファイル作成関数
void tetgen_function::MakeNodeFile(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	//cout<<filename<<" 出力"<<endl;

	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	}

	fout.close();
}


//.nodeファイル作成関数
void tetgen_function::MakeNodeFile_NonAttributeAndBoundary(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	cout<<filename<<" 出力"<<endl;

	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}

	fout.close();
}


//.eleファイル作成関数
void tetgen_function::MakeElemFile(mpsconfig &CON, vector<tetgen_element> &ELEM, char *filename)
{
	ofstream fout(filename);

	fout<<(int)ELEM.size()<<"\t"<<"4"<<"\t"<<"1"<<endl;
	for(int i=0;i<(int)ELEM.size();i++)
	{
		fout<<ELEM[i].id<<"\t";
		for(int n=0;n<4;n++)	fout<<ELEM[i].node[n]<<"\t";
		fout<<ELEM[i].attribute<<endl;
	}

	fout.close();
}


//.faceファイル作成関数
void tetgen_function::MakeFaceFile(mpsconfig &CON, vector<tetgen_facet> &FACE, char *filename)
{
	ofstream fout(filename);

	fout<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fout<<FACE[i].id;
		for(int n=0;n<3;n++)	fout<<" "<<FACE[i].node[n];
		fout<<" "<<FACE[i].boundary;
		fout<<endl;
	}

	fout.clear();
}


//.polyファイル作成関数
void tetgen_function::MakePolyFile(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, char *filename)
{
	//cout<<filename<<" 出力"<<endl;

	ofstream fout(filename);

	//node list (ここでは出力しない)
	fout<<"#node"<<endl;
	fout<<"0"<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	//fout<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	//for(int i=0;i<(int)NODE.size();i++)
	//{
		//fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	//}

	//faset list
	fout<<"#faset"<<endl;
	fout<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fout<<"1"<<" "<<"0"<<" "<<FACE[i].boundary<<endl;
		//fout<<"1"<<" "<<FACE[i].boundary<<endl;
		//fout<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<" "<<FACE[i].boundary<<endl;
		fout<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<endl;
	}

	//hole list
	fout<<"#hole"<<endl;
	fout<<"0"<<endl;
	fout<<endl;

	////////////////////////////region attributeの決定 (配列に格納してから出力する)/////////////////////////////
	//材質の指定は，境界内にある一点の座標を決め，そこの材質を指定することで，同じ境界内にある要素が全てその材質になる．
	//水は分裂を伴い，材質の指定が困難であるため，ここでは行わない．
	//後の材質の修正において，未定義となっている要素を水要素とする．

	vector<region_attribute_list> REGION;
	region_attribute_list temp;
	temp.id=0;
	temp.region_number=0;
	temp.region_attribute=0;
	
	//静電霧化
	if(CON.get_model_number()==14)
	{
		if(CON.get_ea_model()==0)//円柱電極モデル
		{
			//空気
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=(CON.get_ZU()+TET.height_plate+TET.thickness_plate)/2;	//解析領域の上端と平板電極の上面の中間点
			temp.region_number=AIR;
			temp.region_attribute=AIR;
			REGION.push_back(temp);
		
			//平板
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=TET.height_plate+TET.thickness_plate/2.0;	//平板の厚み方向の中間点
			temp.region_number=ELECTRODE2;
			temp.region_attribute=ELECTRODE2;
			REGION.push_back(temp);

			//円柱
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=-TET.length_column/2.0;		//円柱長さ方向の中間点
			temp.region_number=ELECTRODE1;
			temp.region_attribute=ELECTRODE1;
			REGION.push_back(temp);

			//土台
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=-TET.length_column-TET.thickness_base/2.0;		//土台の厚み方向の中間点
			temp.region_number=ELECTRODE1;
			temp.region_attribute=ELECTRODE1;
			REGION.push_back(temp);
		}
		if(CON.get_ea_model()==1)//製品電極(cadベースモデル)
		{
			//空気
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域の上端と平板電極の上面の中間点
			temp.region_number=AIR;
			temp.region_attribute=AIR;
			REGION.push_back(temp);

			//対向電極
			temp.id+=1;
			temp.r[A_X]=0.0016;
			temp.r[A_Y]=0;
			temp.r[A_Z]=0.0029;	//平板の厚み方向の中間点//原点からのz方向距離3mm、直径3mm~4mm
			temp.region_number=ELECTRODE2;
			temp.region_attribute=ELECTRODE2;
			REGION.push_back(temp);

			//放電電極
			temp.id+=1;
			temp.r[A_X]=0.00;
			temp.r[A_Y]=0;
			temp.r[A_Z]=-0.001;	//平板の厚み方向の中間点//原点からのz方向距離3mm、直径3mm~4mm
			temp.region_number=ELECTRODE1;
			temp.region_attribute=ELECTRODE1;
			REGION.push_back(temp);
		}
	}

	//レイリー分裂
	if(CON.get_model_number()==15)
	{
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);
	}

	//CC
	if(CON.get_model_number()==20 ||CON.get_model_number()==21)
	{
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		/*//るつぼ
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0.053;
		temp.r[A_Z]=0.025;
		temp.region_number=CRUCIBLE;
		temp.region_attribute=CRUCIBLE;
		REGION.push_back(temp);
		////*/

		//コイル
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0.053;
		temp.r[A_Z]=0.025;
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);
		///*/
	}

	//渦電流用簡易モデル
	if(CON.get_model_number()==25)
	{
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=TET.thickness_conduct*0.501;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//導体片
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0.0;
		temp.r[A_Z]=0.0;
		temp.region_number=CONDUCT;
		temp.region_attribute=CONDUCT;
		REGION.push_back(temp);
		/////

		/*////コイル
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0.053;
		temp.r[A_Z]=0.025;
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);
		///*/
	}
	////////////////////////////////////////////////////////////////////*/

	//region attribute list
	fout<<"#region attribute"<<endl;
	fout<<(int)REGION.size()<<endl;
	for(int i=0;i<(int)REGION.size();i++)
	{
		fout<<REGION[i].id<<" "<<REGION[i].r[A_X]<<" "<<REGION[i].r[A_Y]<<" "<<REGION[i].r[A_Z]<<" "<<REGION[i].region_number<<" "<<REGION[i].region_attribute<<endl;
	}

	fout.close();
}


//.smeshファイル作成関数
void tetgen_function::MakeSmeshFile(mpsconfig &CON, vector<tetgen_facet> &FACE, char *filename)
{
	double le=CON.get_distancebp();

	ofstream fsmesh(filename);

	//node list (ここでは出力しない)
	fsmesh<<"#node"<<endl;
	fsmesh<<"0"<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	//fsmesh<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	//for(int i=0;i<(int)NODE.size();i++)
	//{
		//fsmesh<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	//}

	//faset list
	fsmesh<<"#faset"<<endl;
	fsmesh<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fsmesh<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<" "<<FACE[i].boundary<<endl;
	}

	//hole list
	fsmesh<<"#hole"<<endl;
	fsmesh<<"0"<<endl;
	fsmesh<<endl;

	//region attribute list
	fsmesh<<"#region attribute"<<endl;
	fsmesh<<"2"<<endl;
	fsmesh<<"1"<<" "<<"0"<<" "<<"0"<<" "<<le*2<<" "<<"1"<<" "<<"1"<<endl;
	fsmesh<<"2"<<" "<<"0"<<" "<<"0"<<" "<<-le*2<<" "<<"2"<<" "<<"2"<<endl;
	
	fsmesh.close();
}


//流体境界面作成
void tetgen_function::SetFluidBoundary(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEw, vector<tetgen_facet> &FACEw)
{
	cout<<"流体境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMw;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;
	temp.attribute=WATER;
	//temp.boundary=WATER;

	double le=CON.get_distancebp();	//粒子間距離
	double rc=TET.radius_column;	//電極半径
	int type;
	int part_no;


	//粒子データの取り込み
	
	//静電無化
	if(CON.get_model_number()==14)
	{
		for(int i=0;i<(int)PART.size();i++)
		{
			part_no=i;
			for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];
			type=PART[i].type;
			if(type==FRFLUID)	temp.boundary=FRFLUID;
			if(type==BOFLUID)	temp.boundary=BOFLUID;
			
			if(temp.r[A_Z]>le*0.5)
			{
				trans.push_back(part_no);
				NODEw.push_back(temp);
				temp.id+=1;
			}
		}
	}
	//レイリー分裂
	if(CON.get_model_number()==15)
	{
		for(int i=0;i<(int)PART.size();i++)
		{
			part_no=i;
			for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];
			type=PART[i].type;
			if(type==FRFLUID)	temp.boundary=FRFLUID;
			if(type==BOFLUID)	temp.boundary=BOFLUID;
			
			trans.push_back(part_no);
			NODEw.push_back(temp);
			temp.id+=1;
		}
	}

	//CC
	if(CON.get_model_number()==20 || CON.get_model_number()==21)
	{
		for(int i=0;i<(int)PART.size();i++)
		{
			part_no=i;
			for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];
			type=PART[i].type;
			//if(type==FRFLUID)	temp.boundary=FRFLUID;
			//if(type==BOFLUID)	temp.boundary=BOFLUID;
			if(type==FLUID && PART[i].surface==OFF)	temp.boundary=FRFLUID;
			if(type==FLUID && PART[i].surface==ON)	temp.boundary=BOFLUID;

			if(PART[i].type==FLUID)//壁粒子はメッシュ作成に含めない
			{
				trans.push_back(part_no);
				NODEw.push_back(temp);
				temp.id+=1;
			}
		}
	}

	//nodeファイル作成
	MakeNodeFile(CON, NODEw, "NODEw.node");

	//.nodeファイル読み取り
	in.load_node("NODEw");

	//まずは流体節点のみで分割
	tetrahedralize("", &in, &out);
		// i 要素内へ節点を追加(今の所無くても機能している)
		// f .faceファイルに境界ではない面も含める
		// e .edgeファイルの出力(ONにするとなぜか止まってしまう)
		// n .neighファイルの出力

	//出力
	out.save_nodes("fluid_whole");
	out.save_elements("fluid_whole");

	//////////////////ここまでで流体節点のみを使って、すべての要素が繋がった凸なメッシュができた(fluid_wholeで確認可能)*/


	///////////////不要な要素の削除

	//.nodeの取得
	GetPointList(NODEw, in, out);
	//.eleの取得
	GetTetrahedronList(ELEMw, in, out);

	//長い要素の除去
	DelThinTetrahedron(CON, TET, NODEw, ELEMw, in, out);

	//節点-要素関係
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//要素-要素関係
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//要素-要素関係より流体表面取得
	GetFacetList_from_neigh(CON, ELEMw, FACEw);

	//飛び出た尖った要素の削除
	DelThinTetrahedron_SharpElem(CON, NODEw, ELEMw, FACEw, in, out);

	/////飛び出た要素を削除したので，もう一度，節点と要素の隣接関係を求め直す
	//節点-要素関係
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//要素-要素関係
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//要素-要素関係より流体表面取得
	GetFacetList_from_neigh(CON, ELEMw, FACEw);

	//粒子番号をNODEw側に代入
	for(int i=0;i<NODEw.size();i++)	NODEw[i].part_no=trans[i];
	
	//表面を構成する節点を選択し，配列番号を詰める 流体内部も粒子の節点を使う場合はコメントアウト
	//SelectFaceNode(CON, NODEw, FACEw);

	/////////////////要素確認用ファイル///////////////////////////////////
	out.save_nodes("boundary_fluid");	//fluid.2.nodeと同じファイル
	MakeElemFile(CON, ELEMw, "boundary_fluid.ele");
	MakeFaceFile(CON, FACEw, "boundary_fluid.face");
	////////////////ここまでで水滴のメッシュが切れた//////////////////////*/
}


//流体境界面作成  外側法線方向ダミー節点法
void tetgen_function::SetFluidBoundary_OutsideDummyNode(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEw, vector<tetgen_facet> &FACEw)
{
	cout<<"流体境界作成(外側法線ダミー節点法)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMw;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;
	temp.attribute=WATER;
	temp.boundary=WATER;
	
	double le=CON.get_distancebp();	//粒子間距離
	double rc=TET.radius_column;	//電極半径
	//double dummy;
	//int part_no;
	//int type;

	int num_dummy=0;	//ダミー節点数
	double maxZ=0;


	//流体節点
	for(int i=0;i<(int)PART.size();i++)
	{
		if(PART[i].r[A_Z]>le*0.5)
		//if((PART[i].type==BOFLUID && PART[i].r[A_Z]>le*0.5) || (PART[i].r[A_Z]>le*0.5 && PART[i].r[A_Z]<le*1.5))
		{
			temp.boundary=PART[i].type;
			for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];
			trans.push_back(i);
			NODEw.push_back(temp);
			temp.id+=1;

			//Z座標の最大値を更新
			if(maxZ<temp.r[A_Z])	maxZ=temp.r[A_Z];
		}
	}
	//cout <<minZ<<endl;
	//Z座標の最小値が負だったら
	//if(minZ<0)	minZ=0.1*le;	//0.1*leとしておく
	//cout <<minZ<<endl;

	//法線方向ダミー節点
	temp.attribute=AIR;
	temp.boundary=AIR;

	//ダミー節点////////////////////////////////////////////
	//法線外側
	for(int i=0;i<(int)PART.size();i++)
	{
		for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];

		if(PART[i].type==BOFLUID)
		{
			for(int d=0;d<3;d++)	temp.r[d]-=PART[i].n[d]*2.0*le;
			if(temp.r[A_Z]>le*0.5)
			{				
				NODEw.push_back(temp);
				temp.id+=1;
				trans.push_back(-1);//粒子に対応しない節点は-1を格納
				num_dummy+=1;
			}
		}
	}//*/

	//電極の上端面（流体の下）
	for(double r=le;r<rc+3*le;r+=le)//広めに置く
	{
		int nr=int(2.0*PI*r/le);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			//double theta=nt*d_theta;
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=0;
			NODEw.push_back(temp);
			temp.id+=1;
			trans.push_back(-1);//粒子に対応しない節点は-1を格納
			num_dummy+=1;
		}
	}//*/

	/*//スーパーボックス////////////////////////////////////////////
	for(int z=0;z<2;z++)
	{
		if(z==0)	temp.r[A_Z]= -5*le;
		if(z==1)	temp.r[A_Z]= maxZ+5*le;

		for(int y=0;y<2;y++)
		{
			if(y==0)	temp.r[A_Y]= -0.001;
			if(y==1)	temp.r[A_Y]=  0.001;

			for(int x=0;x<2;x++)
			{
				if(x==0)	temp.r[A_X]= -0.001;
				if(x==1)	temp.r[A_X]=  0.001;

				NODEw.push_back(temp);
				temp.id+=1;
				trans.push_back(-1);//粒子に対応しない節点は-1を格納
				num_dummy+=1;
			}
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEw, "fluid.1.node");

	//.nodeファイル読み取り
	in.load_node("fluid.1");

	//分割
	tetrahedralize("", &in, &out);
		// i 要素内へ節点を追加(今の所無くても機能している)
		// f .faceファイルに境界ではない面も含める
		// e .edgeファイルの出力(ONにするとなぜか止まってしまう)
		// n .neighファイルの出力

	//出力
	out.save_nodes("fluid.2");
	out.save_elements("fluid.2");

	//////////////////ここまでで流体表面とその外側のダミー節点を使った凸な形状のメッシュができる(fluid.2で確認可能)*/


	///////////////不要な要素の削除

	//.nodeの取得
	GetPointList(NODEw, in, out);
	//.eleの取得
	GetTetrahedronList(ELEMw, in, out);

	//ダミー要素の除去
	DelTetrahedron_OutsideDummy(CON, NODEw, ELEMw, in, out);

	//節点-要素関係
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//要素-要素関係
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//要素-要素関係より流体表面取得
	GetFacetList_from_neigh(CON, ELEMw, FACEw);

	//ダミー節点の削除
	DelDummyNode(CON, NODEw, FACEw, num_dummy);
	//表面を構成する節点を選択し，配列番号を詰める 流体内部も粒子の節点を使う場合はコメントアウト
	//SelectFaceNode(CON, NODEw, FACEw);

	//テスト
	//MakeNodeFile(CON, NODEw, "NODEw.node");
	//MakeFaceFile(CON, FACEw, "FACEw.face");


	/////////////////要素確認用ファイル///////////////////////////////////
	//out.save_nodes("boundary_fluid");	//fluid.2.nodeと同じファイル
	MakeNodeFile_NonAttributeAndBoundary(CON, NODEw, "boundary_fluid.node");
	MakeElemFile(CON, ELEMw, "boundary_fluid.ele");
	MakeFaceFile(CON, FACEw, "boundary_fluid.face");
	////////////////ここまでで水滴のメッシュが切れた//////////////////////*/
}


//流体境界面作成  内側法線方向ダミー節点法(開発途中)
void tetgen_function::SetFluidBoundary_InsideDummyNode(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEw, vector<tetgen_facet> &FACEw)
{
	cout<<"流体境界作成(内側流体法線ダミー節点法)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMw;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;
	temp.attribute=WATER;
	temp.boundary=WATER;

	double le=CON.get_distancebp();	//粒子間距離
	double rc=TET.radius_column;	//電極半径
	//int type;
	//int part_no;
	double minZ=le;


	//粒子データの取り込み
	for(int i=0;i<(int)PART.size();i++)
	{		
		if(PART[i].type==BOFLUID && PART[i].r[A_Z]>le*0.5)
		{
			temp.boundary=PART[i].type;
			for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];
			trans.push_back(i);
			NODEw.push_back(temp);
			temp.id+=1;

			//Z座標の最小値を更新
			if(minZ>temp.r[A_Z])	minZ=temp.r[A_Z];
		}
	}
	//cout <<minZ<<endl;
	//Z座標の最小値が負だったら
	if(minZ<0)	minZ=0.1*le;	//0.1*leとしておく
	//cout <<minZ<<endl;

	//内側にダミー節点を置く
	for(int i=0;i<(int)PART.size();i++)
	{
		temp.boundary=UNDEFINED;
		for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];

		if(PART[i].type==BOFLUID/* || (temp.r[A_Z]>le*0.5 && temp.r[A_Z]<le*1.5)*/)
		{
			for(int d=0;d<3;d++)	temp.r[d]+=PART[i].n[d]*0.5*le;
			if(temp.r[A_Z]>le*0.5)
			{				
				NODEw.push_back(temp);
				temp.id+=1;
				trans.push_back(-1);//粒子に対応しない節点は-1を格納
			}
		}
		if(PART[i].type==INWALL/* && sqrt(temp.r[A_X]*temp.r[A_X]+temp.r[A_Y]*temp.r[A_Y])<rc*/)//電極の真上
		{
			temp.r[A_Z]=minZ-0.1*le;
			//temp.r[A_Z]=0.5*le;
			NODEw.push_back(temp);
			temp.id+=1;
			trans.push_back(-1);//粒子に対応しない節点は-1を格納
		}//*/
	}

	/*//電極の真上(強制的に並べる)
	temp.boundary=UNDEFINED;
	for(double r=le;r<rc+0.1*le;r+=le)
	{
		int nr=int(2.0*PI*r/le);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			//double theta=nt*d_theta;
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=0.5*le;
			NODEw.push_back(temp);
			temp.id+=1;
			trans.push_back(-1);//粒子に対応しない節点は-1を格納
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEw, "fluid.1.node");

	//.nodeファイル読み取り
	in.load_node("fluid.1");

	//まずは流体節点のみで分割
	tetrahedralize("", &in, &out);
		// i 要素内へ節点を追加(今の所無くても機能している)
		// f .faceファイルに境界ではない面も含める
		// e .edgeファイルの出力(ONにするとなぜか止まってしまう)
		// n .neighファイルの出力

	//出力
	out.save_nodes("fluid.2");
	out.save_elements("fluid.2");

	//////////////////ここまでで流体節点のみを使って、すべての要素が繋がった凸なメッシュができた(fluid.2で確認可能)*/


	///////////////不要な要素の削除

	//.nodeの取得
	GetPointList(NODEw, in, out);
	//.eleの取得
	GetTetrahedronList(ELEMw, in, out);

	//ダミー要素の除去
	DelTetrahedron_InsideDummy(CON, NODEw, ELEMw, in, out);

	//節点-要素関係
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//要素-要素関係
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//要素-要素関係より流体表面取得
	GetFacetList_from_neigh(CON, ELEMw, FACEw);


	//ダミー節点の削除
	//DelDummyNode(CON, NODEw, FACEw, num_dummy);


	//粒子番号をNODEw側に代入
	for(int i=0;i<NODEw.size();i++)	NODEw[i].part_no=trans[i];
	
	//表面を構成する節点を選択し，配列番号を詰める 流体内部も粒子の節点を使う場合はコメントアウト
	SelectFaceNode(CON, NODEw, FACEw);

	//テスト
	//MakeNodeFile(CON, NODEw, "NODEw.node");
	//MakeFaceFile(CON, FACEw, "FACEw.face");


	/////////////////要素確認用ファイル///////////////////////////////////
	out.save_nodes("boundary_fluid");	//fluid.2.nodeと同じファイル
	MakeElemFile(CON, ELEMw, "boundary_fluid.ele");
	MakeFaceFile(CON, FACEw, "boundary_fluid.face");
	////////////////ここまでで水滴のメッシュが切れた//////////////////////*/
}


//長い要素の除去
void tetgen_function::DelThinTetrahedron(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"不要な要素の削除  ";

	double le=CON.get_distancebp();
	double delL;
	int flag=0;

	//del_length.datがあれば読み込み
	ifstream del("del_length.dat");
	if(!del)
	{
		cout<<"tetgen_configより削除判定辺長さを決定  ";
		delL=TET.del_length;
	}
	else
	{
		cout<<"del_length.datより削除判定辺長さを決定  ";
		del>>delL;
		flag=1;		//ファイルから読み込んだらフラグON
	}
	del.close();

	if(flag==1)//ファイルの数字を戻しておく
	{
		ofstream del2("del_length.dat");
		del2<<TET.del_length<<endl;
		del2.close();
	}

	cout<<"del_length="<<delL<<endl;


	cout<<"除去前要素数: "<<(int)ELEM.size()<<endl;

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag1=UNDEFINED;
		int flag2=UNDEFINED;
		int flag3=UNDEFINED;	//1で削除
		int del=OFF;
		int count=0;

		/*//4点が表面節点で構成されていればフラグ1ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag1=ON;
		else			flag1=OFF;//*/

		/*//4点が内部節点で構成されていればフラグ2ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==FRFLUID)	count+=1;
		}
		if(count==4)	flag2=ON;
		else			flag2=OFF;//*/

		//一つでも長い辺があればフラグON
		for(int n1=0;n1<4;n1++)
		{
			for(int n2=n1+1;n2<4;n2++)
			{
				double dis=Distance(NODE[ELEM[i].node[n1]+out.firstnumber], NODE[ELEM[i].node[n2]+out.firstnumber]);
				if(dis>le*delL)	flag3=ON;
			}
		}//*/

		if(flag3==ON)	del=ON;
		//if(flag2==ON)	del=OFF;
		//if(flag3==2)	del=ON;

		//削除
		if(del==ON)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//イテレータ初期化
			it+=i;				//iを指定
			it=ELEM.erase(it);	//削除してイテレータを返す
			del_count++;
		}
		else i++;
	}

	//要素番号の振り直し
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"除去後要素数: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"の要素を削除 ----------OK"<<endl;
}

//飛び出て尖った要素の削除
void tetgen_function::DelThinTetrahedron_SharpElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out)
{
	//要素-要素関係が求まっていることを前提とする
	cout<<"尖った要素の要素の削除  ";

	double le=CON.get_distancebp();
	int flag=0;

	//boundaryをFRFLUID(=0)に初期化
	for(int i=0;i<(int)NODE.size();i++)	NODE[i].boundary=FRFLUID;

	//FACEに含まれる節点のみBOFLUIDにする
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)	NODE[FACE[i].node[n]].boundary=BOFLUID;
	}


	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag1=UNDEFINED;	//1で削除
		int del=OFF;
		int count1=0;
		int count2=0;

		//4面のうち2面以上接する要素がないものは削除
		//4点全てがBOFLUIDであるものは削除
		for(int n=0;n<4;n++)
		{
			if(ELEM[i].nei_elem[n]==-1)	count1+=1;
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count2+=1;
		}
		if(count1>=2 || count2==4)	flag1=ON;
		else						flag1=OFF;//*/

		if(flag1==ON)	del=ON;

		//削除
		if(del==ON)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//イテレータ初期化
			it+=i;				//iを指定
			it=ELEM.erase(it);	//削除してイテレータを返す
			del_count++;
		}
		else i++;
	}

	//要素番号の振り直し
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"除去後要素数: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"の要素を削除 ----------OK"<<endl;

}

//不要な要素の除去(外側ダミー節点法用)
void tetgen_function::DelTetrahedron_OutsideDummy(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"不要な要素の削除"<<endl;
	cout<<"除去前要素数: "<<(int)ELEM.size()<<endl;

	double le=CON.get_distancebp();

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag=0;	//1で削除
		
		//1つでもダミー(空気)節点があればフラグON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==AIR)
			{
				flag=1;
				break;
			}
		}//*/

		/*//4点が表面節点で構成されていればフラグON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag=1;//*/


		if(flag==1)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//イテレータ初期化
			it+=i;				//iを指定
			it=ELEM.erase(it);	//削除してイテレータを返す
			del_count++;
		}
		else i++;
	}

	//要素番号の振り直し
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"除去後要素数: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"の要素を削除 ----------OK"<<endl;
}


//不要な要素の除去(内側ダミー節点法用)
void tetgen_function::DelTetrahedron_InsideDummy(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"不要な要素の削除"<<endl;
	cout<<"除去前要素数: "<<(int)ELEM.size()<<endl;

	double le=CON.get_distancebp();

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag=0;	//1で削除
		
		//4点が表面節点で構成されていればフラグON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag=1;//*/

		if(flag==1)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//イテレータ初期化
			it+=i;				//iを指定
			it=ELEM.erase(it);	//削除してイテレータを返す
			del_count++;
		}
		else i++;
	}

	//要素番号の振り直し
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"除去後要素数: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"の要素を削除 ----------OK"<<endl;
}


//節点-要素関係(tetgenioのedgeリストから取得) ※直前にtetgenioからデータを取得しておくこと
void tetgen_function::SetRelation_NodeNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"節点-節点関係";

	//一応初期化
	for(int i=0;i<(int)NODE.size();i++)
	{
		NODE[i].nei_node.clear();
	}

	for(int i=0;i<out.numberofedges;i++)
	{
		int node1=out.edgelist[i*2+0];
		int node2=out.edgelist[i*2+1];
		
		NODE[node1].nei_node.push_back(node2);
		NODE[node2].nei_node.push_back(node1);
	}

	cout<<"----------OK"<<endl;
}


//節点-要素関係
void tetgen_function::SetRelation_NodeElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"節点-要素関係";

	//一応初期化
	for(int i=0;i<(int)NODE.size();i++)
	{
		NODE[i].nei_node.clear();
		NODE[i].nei_elem.clear();
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			NODE[ELEM[i].node[n]].nei_elem.push_back(i);	//要素を追加
		}
	}

	/*//出力 および最大数・最小数の出力
	//int max=0;
	//int min=(int)NODE[0].nei_elem.size();

	ofstream fout("neigh_node-elem.dat");
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id;
		for(int n=0;n<(int)NODE[i].nei_elem.size();n++)
		{
			fout<<" "<<NODE[i].nei_elem[n];
		}
		fout<<endl;

		//最大最小更新
		//if(max<(int)NODE[i].nei_elem.size())	max=(int)NODE[i].nei_elem.size();
		//if(min>(int)NODE[i].nei_elem.size())	min=(int)NODE[i].nei_elem.size();
	}
	fout.clear();

	//cout<<"最大数: "<<max<<endl;
	//cout<<"最小数: "<<min<<endl;
	//*/

	cout<<"----------OK"<<endl;
}


//要素-要素関係
void tetgen_function::SetRelation_ElemElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"要素-要素関係";
	
	vector<int> nei_all;
	
	//初期化
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)	ELEM[i].nei_elem[n]=-2;	//未定義を-2としておく
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		//4つの節点の節点-要素関係にある要素を格納する
		nei_all.clear();
		
		for(int n=0;n<4;n++)
		{
			for(int j=0;j<(int)NODE[ELEM[i].node[n]].nei_elem.size();j++)
			{
				int elem=NODE[ELEM[i].node[n]].nei_elem[j];
				if(elem!=i)	nei_all.push_back(elem);	//要素iは除く
			}
		}//nei_allに格納完了(同じ要素番号が含まれている可能性あり)

		//確認用
		//if(i==10000)	for(int j=0;j<(int)nei_all.size();j++)	cout<<nei_all[j]<<endl;

		//面を探索
		for(int ni=0;ni<4;ni++)
		{
			//まずは要素iの面を指定
			int face[3];	//3つ番号で面を指定
			int c=0;//数え上げ変数

			for(int f=0;f<4;f++)
			{
				if(ni!=f)
				{
					face[c]=ELEM[i].node[f];
					c++;
				}
			}//face[3]にn番目の面が格納

			//面を探索
			int correct_nei=-1;
			for(int j=0;j<(int)nei_all.size();j++)
			{
				int count=0;//このカウントが3になれば確定

				for(int nj=0;nj<4;nj++)
				{
					int node_j=ELEM[nei_all[j]].node[nj];
					for(int f=0;f<3;f++)
					{
						if(node_j==face[f])	count++;
					}	
				}
				if(count==3)
				{
					correct_nei=nei_all[j];
					break;
				}
			}//もし見つからなかったらcorrecr_neiには-1が入っている

			ELEM[i].nei_elem[ni]=correct_nei;
		}
	}//*/

	/*//出力
	ofstream fout("neigh.dat");
	for(int i=0;i<(int)ELEM.size();i++)
	{
		fout<<ELEM[i].id;
		for(int n=0;n<4;n++)
		{
			fout<<" "<<ELEM[i].nei_elem[n];		
		}
		fout<<endl;
	}
	fout.clear();
	//*/

	cout<<"----------OK"<<endl;
}


//要素除去後の流体表面定義
void tetgen_function::GetFacetList_from_neigh(mpsconfig &CON, vector<tetgen_element> &ELEM, vector<tetgen_facet> &FACE)
{
	cout<<"要素除去後の流体表面定義";

	//初期化
	FACE.clear();
	
	int id=0;	//id用(表面の数)
	
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			if(ELEM[i].nei_elem[n]==-1)//対面が存在しない→表面
			{
				tetgen_facet temp;	
				int c=0;

				for(int f=0;f<4;f++)
				{
					if(n!=f)
					{
						temp.node[c]=ELEM[i].node[f];
						c++;
					}
				}

				temp.id=id;
				temp.boundary=WATER;
				FACE.push_back(temp);
				id++;
			}
		}
	}

	cout<<"----------OK"<<endl;
}


//表面節点以外を削除，表面節点に番号を振りなおす
void tetgen_function::SelectFaceNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE)
{
	//boundaryをフラグに使わせてもらう。
	//boundary==-1は表面を構成していない節点→削除
	//boundary==-2は表面を構成している節点→新しい節点番号を与える
	
	//とりあえず初期化
	for(int i=0;i<(int)NODE.size();i++)	NODE[i].boundary=-1;
	
	//FACEにある節点はflagを-2に
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			NODE[FACE[i].node[n]].boundary=-2;
		}
	}

	//新しい節点番号の決定
	int id=0;
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].boundary==-2)
		{
			NODE[i].boundary=id;
			id++;
		}
	}
	//表面節点にはboundaryに新しい節点番号が入る。内部節点には-1が入る

	//FACEの構成節点番号の変換
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			FACE[i].node[n]=NODE[FACE[i].node[n]].boundary;	//ここで-1や-2となるものはないはず。あれば上の処理が間違っている
		}
	}

	//内部節点の削除 NODEの節点番号の変換 boundaryを元に戻す
	int k=0;
	while(k<(int)NODE.size())
	{
		if(NODE[k].boundary==-1)
		{
			vector<tetgen_node>::iterator it=NODE.begin();	//イテレータ初期化
			it+=k;				//k番目を指定
			it=NODE.erase(it);	//削除してイテレータを返す
		}
		else
		{
			NODE[k].id=NODE[k].boundary;
			NODE[k].boundary=WATER;
			k++;
		}
	}

}


//ダミー節点を削除
void tetgen_function::DelDummyNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int num_dummy)
{
	//NODEの中には前半に流体表面節点，後半にダミー節点が固まって格納されているので，後半のダミー節点の部分のみを消せばよい
	//ダミー要素を消してから表面データを取得しているので，表面を構成する節点番号の変更は不要

	//ダミー変数の数だけpopbackで末尾の要素から消す
	for(int i=0;i<num_dummy;i++)
	{
		NODE.pop_back();
	}
}


//空気境界面作成
void tetgen_function::SetAirBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	cout<<"空気境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp_n;
	temp_n.id=0;
	temp_n.attribute=AIR;
	temp_n.boundary=AIR;

	//分割数決定  1割だけオフセットしている
	double dL=TET.fine_air;
	int lx = int((CON.get_XL()-0.1*dL)/dL);
	int ux = int((CON.get_XR()+0.1*dL)/dL);
	int ly = int((CON.get_YD()-0.1*dL)/dL);
	int uy = int((CON.get_YU()+0.1*dL)/dL);
	int lz = int((CON.get_ZD()-0.1*dL)/dL);
	int uz = int((CON.get_ZU()+0.1*dL)/dL);

	//接点データ作成
	//静電無化
	if(CON.get_model_number()==14)
	{
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//解析領域の端に来たとき節点を置く
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}//*/

	//レイリー分裂
	if(CON.get_model_number()==15)
	{
		temp_n.boundary=ELECTRODE1;

		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//解析領域の端に来たとき節点を置く
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}

	//CC
	if(CON.get_model_number()==20 || CON.get_model_number()==21)
	{
		//temp_n.boundary=ELECTRODE1;
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//解析領域の端に来たとき節点を置く
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}//*/

	//導体片
	if(CON.get_model_number()==25)
	{
		//temp_n.boundary=ELECTRODE1;
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					//if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//解析領域の端に来たとき節点を置く
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}

	//nodeファイル作成
	MakeNodeFile(CON, NODEa, "NODEa1.node");

	in.load_node("NODEa1");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_air1");
	out.save_elements("boundary_air1");
	out.save_faces("boundary_air1");
	//*/

	//境界面データ取得
	GetFacetList(FACEa, in, out, AIR);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEa1, "FACEa1.face");

}

//空気境界面作成
void tetgen_function::SetAirBoundary2(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	cout<<"空気境界(細分化用)作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp_n;
	temp_n.id=0;
	temp_n.attribute=AIR;
	temp_n.boundary=AIR;

	
	//導体片
	if(CON.get_model_number()==25)
	{
		
		//分割数決定
		double	le = CON.get_distancebp();
		double	rc = TET.radius_column;
		double	h =  TET.thickness_conduct;						//導体片厚さ
		double	dh = TET.fine_conduct_z;							//厚み方向メッシュ粗さ
		int		nh = int((h/2+0.1*dh)/dh);	//厚み方向分割数
		double	L =  TET.length_conduct;
		double	dL = TET.fine_conduct_xy;						//xy方向メッシュ粗さ
		int		nL = int((L/2+0.1*dL)/dL);	//xy方向分割数(片側)

		int lay=3;

		//その他の面
		for(int z=-nh-lay;z<=nh+lay;z++)
		{
			for(int y=-nL-lay;y<=nL+lay;y++)
			{
				for(int x=-nL-lay;x<=nL+lay;x++)
				{
					if(x<-nL || x>nL || y<-nL || y>nL || z>nh || z<-nh)	//平板電極領域の端より外側に空気節点を置く
					{
						//if(sqrt(temp.r[A_X]*temp.r[A_X]+temp.r[A_Y]*temp.r[A_Y])>rc+le || z>0)
						{
							temp_n.r[A_X]=x*dL;//+dh*(x-nL);
							temp_n.r[A_Y]=y*dL;//+dh*(y-nL);
							temp_n.r[A_Z]=z*dh;
							NODEa.push_back(temp_n);
							temp_n.id+=1;
						}
					}
				}
			}
		}//*/

		
		//分割数決定  1割だけオフセットしている
		double dL2=TET.fine_air;
		int lx = int((CON.get_XL()-0.1*dL2)/dL2);
		int ux = int((CON.get_XR()+0.1*dL2)/dL2);
		int ly = int((CON.get_YD()-0.1*dL2)/dL2);
		int uy = int((CON.get_YU()+0.1*dL2)/dL2);
		int lz = int((CON.get_ZD()-0.1*dL2)/dL2);
		int uz = int((CON.get_ZU()+0.1*dL2)/dL2);


	
		//temp_n.boundary=ELECTRODE1;
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x*dL2<(-nL-lay)*dL || x*dL2>(nL+lay)*dL || y*dL2<(-nL-lay)*dL || y*dL2>(nL+lay)*dL || z*dL2>(nh+lay)*dh || z*dL2<(-nh-lay)*dh)	//平板電極領域の端より外側に空気節点を置く
					{
						temp_n.r[A_X]=x*dL2;
						temp_n.r[A_Y]=y*dL2;
						temp_n.r[A_Z]=z*dL2;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}


	cout<<"細分化節点位置の計上完了"<<endl;
	//nodeファイル作成
	MakeNodeFile(CON, NODEa, "NODEa1.node");

	in.load_node("NODEa1");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_air2");
	out.save_elements("boundary_air2");
	out.save_faces("boundary_air2");
	//*/

	//境界面データ取得
	GetFacetList(FACEa, in, out, AIR);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEa1, "FACEa1.face");

}


//水滴表面付近追加節点
void tetgen_function::SetAirFineBoundary(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	//水滴の表面の法線方向に何層かのメッシュ層を作成する

	cout<<"水滴表面付近追加節点"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=AIR;
	temp.boundary=0;

	//分割数決定
	double le=CON.get_distancebp();
	double r=TET.radius_column+5*le;		//円筒半径
	double uz=TET.height_plate-20*le;		//円筒最大高さ
	double lz=-5*le;		//円筒最小高さ -5*le
	double dL=le;					//円筒長さ方向メッシュ粗さ


	/*//円筒部分
	for(double z=lz;z<uz+0.1*dL;z+=dL)
	{
		int nr=int(2.0*PI*r/le);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);
			temp.r[A_Z]=z;
			NODEa.push_back(temp);
			temp.id+=1;
		}
	}//*/
	
	double thick=TET.thick_layer;

	//水滴法線外側方向
	for(int i=0;i<(int)PART.size();i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==ON)//PART[i].type==BOFLUID 
		{
			for(int n=1;n<=TET.num_layer_out;n++)
			{
				temp.r[A_X]=PART[i].r[A_X]-PART[i].n[A_X]*thick*le*n;
				temp.r[A_Y]=PART[i].r[A_Y]-PART[i].n[A_Y]*thick*le*n;
				temp.r[A_Z]=PART[i].r[A_Z]-PART[i].n[A_Z]*thick*le*n;
				NODEa.push_back(temp);
				temp.id+=1;
			}
		}
	}
	//*/
		
	//水滴法線内側方向
	for(int i=0;i<(int)PART.size();i++)
	{
		if(PART[i].type==FLUID && PART[i].surface==ON)//if(PART[i].type==BOFLUID)
		{
			for(int n=1;n<=TET.num_layer_in;n++)
			{
				temp.r[A_X]=PART[i].r[A_X]+PART[i].n[A_X]*thick*le*n;
				temp.r[A_Y]=PART[i].r[A_Y]+PART[i].n[A_Y]*thick*le*n;
				temp.r[A_Z]=PART[i].r[A_Z]+PART[i].n[A_Z]*thick*le*n;
				NODEa.push_back(temp);
				temp.id+=1;
			}
		}
	}
	//*/


	//nodeファイル作成
	MakeNodeFile(CON, NODEa, "NODEa2.node");

	in.load_node("NODEa2");
	tetrahedralize("", &in, &out);
	//out.save_nodes("boundary_air2");
	//out.save_elements("boundary_air2");
	//out.save_faces("boundary_air2");
	//*/

	//境界面データ取得
	//GetFacetList(FACEa, in, out, AIR);

	//上面と下面を削除する
	int i=0;
	while(i<(int)FACEa.size())
	{
		double dis1=Distance(NODEa[FACEa[i].node[0]], NODEa[FACEa[i].node[1]]);
		double dis2=Distance(NODEa[FACEa[i].node[1]], NODEa[FACEa[i].node[2]]);
		double dis3=Distance(NODEa[FACEa[i].node[2]], NODEa[FACEa[i].node[0]]);

		/*int flag=0;//3になったら削除
		for(int n=0;n<3;n++)
		{
			if(NODEa[FACEa[i].node[n]].r[A_Z]>uz-2*le)	flag+=1;
			if(NODEa[FACEa[i].node[n]].r[A_Z]<lz+2*le)	flag+=1;
		}*/
		//if(flag>=1)
		if(dis1>2.0*le || dis2>2.0*le || dis3>2.0*le)
		{
			vector<tetgen_facet>::iterator it=FACEa.begin();	//イテレータ初期化
			it+=i;				//iを指定
			it=FACEa.erase(it);	//削除してイテレータを返す
		}
		else i++;
	}
	//id再振り分け
	for(int i=0;i<(int)FACEa.size();i++)
	{
		FACEa[i].id=i;
	}

	//.faceファイル作成
	//MakeFaceFile(CON, FACEa, "FACEa2.face");
}

//コイル境界面作成
void tetgen_function::SetCoilBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEd, vector<tetgen_facet> &FACEd)
{
	cout<<"コイル境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=COIL;
	temp.boundary=COIL;

	//ringcoil.nodeからの読み込み 
	//オプションをつけずにstlを読み込めば、境界のみの.nodeがつくられる(tetgen ○○.stl)

	//あらかじめstlから作成した.node,.faceをNODE,FACEに与える

	//node
	ifstream fp("coil_ring.node");
	if(!fp) cout<<"can not open coil.node"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//最初の3行分進める
	
	int node_num=0;
	int stuff=0;
	fp>>node_num;
	for(int i=0;i<3;i++) fp>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<node_num;i++)
	{
		temp.id=i;
		fp>>stuff;
		fp>>temp.r[A_X];
		fp>>temp.r[A_Y];
		fp>>temp.r[A_Z];
		//fp>>temp.r[A_Y];
		//fp>>temp.r[A_Z];
		NODEd.push_back(temp);
		
	}
	fp.close();
	
	//MakeNodeFile(CON, NODEd, "NODEd.node");

	//nodeファイル作成
	MakeNodeFile(CON, NODEd, "NODEd.node");

	in.load_node("NODEd");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_coil");
	out.save_elements("boundary_coil");
	out.save_faces("boundary_coil");
	///*/

	/*/face
	ifstream fp2("coil_ring.face");
	if(!fp2) cout<<"can not open coil.face"<<endl;
	fp2.unsetf(ifstream::dec);
	fp2.setf(ifstream::skipws);

	tetgen_facet tempf;

	int face_num=0;
	stuff=0;
	fp2>>face_num;
	fp2>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<face_num;i++)
	{
		tempf.id=i+out.firstnumber;
		fp2>>stuff;
		for(int n=0;n<3;n++) 
		{
			fp2>>tempf.node[n];
			tempf.node[n]-=1;//他の方法で作成したNODEなどにそろえるため、開始番号を0にする必要がある
		}
		tempf.boundary=COIL;
		FACEd.push_back(tempf);
	}
	fp2.close();
	////*/

	//境界面データ取得
	GetFacetList(FACEd, in, out, COIL);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEd, "FACEd.face");
}

//導電片境界作成
void tetgen_function::SetConductBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEb, vector<tetgen_facet> &FACEb)
{
	cout<<"導体境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=CONDUCT;//とりあえず銅の物性となるcrucibleに材質を指定
	temp.boundary=CONDUCT;

	//分割数決定
	double	le = CON.get_distancebp();
	double	rc = TET.radius_column;
	double	h =  TET.thickness_conduct;						//導体片厚さ
	double	dh = TET.fine_conduct_z;							//厚み方向メッシュ粗さ
	int		nh = int((h/2+0.1*dh)/dh);	//厚み方向分割数
	double	L =  TET.length_conduct;
	double	dL = TET.fine_conduct_xy;						//xy方向メッシュ粗さ
	int		nL = int((L/2+0.1*dL)/dL);	//xy方向分割数(片側)

	//その他の面
	for(int z=-nh;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				//if(x==-nL || x==nL || y==-nL || y==nL || z==nh || z==-nh)	//平板電極領域の端に来たとき節点を置く
				{
					//if(sqrt(temp.r[A_X]*temp.r[A_X]+temp.r[A_Y]*temp.r[A_Y])>rc+le || z>0)
					{
						temp.r[A_X]=x*dL;
						temp.r[A_Y]=y*dL;
						temp.r[A_Z]=z*dh;
						NODEb.push_back(temp);
						temp.id+=1;
					}
				}
			}
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEb, "NODEb.node");

	in.load_node("NODEb");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_conduct");
	out.save_elements("boundary_conduct");
	out.save_faces("boundary_conduct");
	//*/

	//境界面データ取得
	GetFacetList(FACEb, in, out, CONDUCT);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}

//るつぼ境界面作成
void tetgen_function::SetCrucibleBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEd, vector<tetgen_facet> &FACEd)
{
	cout<<"るつぼ境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=CRUCIBLE;
	temp.boundary=CRUCIBLE;

	//ringcoil.nodeからの読み込み 
	//オプションをつけずにstlを読み込めば、境界のみの.nodeがつくられる(tetgen ○○.stl)

	//あらかじめstlから作成した.node,.faceをNODE,FACEに与える

	//node
	ifstream fp("crucible.node");
	if(!fp) cout<<"can not open OutputMeshAndSolution-mass.txt"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//最初の3行分進める
	
	int node_num=0;
	int stuff=0;
	fp>>node_num;
	for(int i=0;i<3;i++) fp>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<node_num;i++)
	{
		temp.id=i;
		fp>>stuff;

		fp>>temp.r[A_X];
		fp>>temp.r[A_Y];
		fp>>temp.r[A_Z];
		

		/*
		fp>>temp.r[A_X];
		fp>>temp.r[A_Z];//もとのstlがy,zが入れ替わっている(シンフォニア提供CAD)のでこちらで修正
		fp>>temp.r[A_Y];
		*/
		temp.r[A_X]*=1e-3;
		temp.r[A_Y]*=1e-3;
		temp.r[A_Z]*=1e-3;
		/*/
		fp>>temp.r[A_X];
		fp>>temp.r[A_Y];
		fp>>temp.r[A_Z];
		*/
		NODEd.push_back(temp);
		
	}
	fp.close();
	
	//MakeNodeFile(CON, NODEd, "NODEd.node");

	/*/nodeファイル作成
	MakeNodeFile(CON, NODEd, "NODEd.node");

	in.load_node("NODEd");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_cc");
	out.save_elements("boundary_cc");
	out.save_faces("boundary_cc");
	///*/

	/////face
	ifstream fp2("crucible.face");
	if(!fp2) cout<<"can not open OutputMeshAndSolution-mass.txt"<<endl;
	fp2.unsetf(ifstream::dec);
	fp2.setf(ifstream::skipws);

	tetgen_facet tempf;

	int face_num=0;
	stuff=0;
	fp2>>face_num;
	fp2>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<face_num;i++)
	{
		tempf.id=i+out.firstnumber;
		fp2>>stuff;
		for(int n=0;n<3;n++) 
		{
			fp2>>tempf.node[n];
			tempf.node[n]-=1;//他の方法で作成したNODEなどにそろえるため、開始番号を0にする必要がある
		}
		tempf.boundary=CRUCIBLE;
		FACEd.push_back(tempf);
	}
	fp2.close();
	////////*/

	//境界面データ取得
	GetFacetList(FACEd, in, out, CRUCIBLE);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEd, "FACEd.face");
}


//平板電極境界面作成
void tetgen_function::SetPlateElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"平板電極境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE2;
	temp.boundary=ELECTRODE2;

	//分割数決定
	double	h = TET.height_plate;						//平板高さ
	double	dh = TET.fine_plate_t;						//厚み方向メッシュ粗さ
	int		nh = int((TET.thickness_plate+0.1*dh)/dh);	//厚み方向分割数
	double	dL = TET.fine_plate_L;						//xy方向メッシュ粗さ
	int		nL = int((TET.length_plate/2+0.1*dL)/dL);	//xy方向分割数(片側)

	//接点データ作成
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==0 || z==nh)	//平板電極領域の端に来たとき節点を置く
				{
					temp.r[A_X]=x*dL;
					temp.r[A_Y]=y*dL;
					temp.r[A_Z]=z*dh+h;
					//temp.attribute=ELECTRODE;
					//temp.boundary=ELECTRODE;
					NODEp.push_back(temp);
					temp.id+=1;
				}
			}
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEp, "NODEp.node");

	in.load_node("NODEp");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_plate");
	out.save_elements("boundary_plate");
	out.save_faces("boundary_plate");
	//*/

	//境界面データ取得
	GetFacetList(FACEp, in, out, ELECTRODE2);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//円柱電極境界面作成
void tetgen_function::SetColumnElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"円柱電極境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//分割数決定
	double le=CON.get_distancebp();
	double rc=TET.radius_column;	//円柱半径
	double L=TET.length_column;		//円柱長さ
	double dL=TET.fine_column_L;	//円柱長さ方向メッシュ粗さ
	//double z0=0;					//上面の位置

	
	//上面
	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=0;
	NODEc.push_back(temp);
	temp.id+=1;

	for(double r=le;r<rc+0.1*le;r+=le)
	{
		int nr=int(2.0*PI*r/le);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			//double theta=nt*d_theta;
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=0;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}


	//側面
	int flag=0;
	double dz=le;
	double z=0;
	rc*=1.005;	//わずかに太くする(メッシュが繋がるのを防ぐため)
	int nr=int(2.0*PI*rc/(2*le));//leが2倍されていることに注意
	double d_theta=360.0/(double)nr;
	
	while(1)
	{
		//z方向への移動
		if(flag==0)			dz*=1.05;
		else if(flag==1)	dz=dL;
		z-=dz;

		if(dz>dL)	flag=1;
		if(z<-L+le)
		{
			z=-L+le;
			//flag=2;
			break;
		}

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=rc*sin(theta*PI/180.0);
			temp.r[A_Y]=rc*cos(theta*PI/180.0);
			temp.r[A_Z]=z;
			NODEc.push_back(temp);
			temp.id+=1;
		}

		//if(flag==2)	break;
	}

	//下面
	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L+le;
	NODEc.push_back(temp);
	temp.id+=1;

	//for(double r=le;r<rc+0.1*le;r+=le)
	for(double r=rc;r>le-0.1*le;r-=2*le)//leが2倍されていることに注意
	{
		int nr=int(2.0*PI*r/(2*le));//leが2倍されていることに注意
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L+le;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}//*/


	//nodeファイル作成
	MakeNodeFile(CON, NODEc, "NODEc.node");

	in.load_node("NODEc");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_column");
	out.save_elements("boundary_column");
	out.save_faces("boundary_column");
	//*/

	//境界面データ取得
	GetFacetList(FACEc, in, out, ELECTRODE1);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEc, "FACEc.face");
}

//対向電極境界面作成
void tetgen_function::SetCounterElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"対向電極境界作成(cad読み込み)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE2;
	temp.boundary=ELECTRODE2;

	
	in.load_stl("ea2");//モデルの座標系で出力するとなぜか読み込んでくれない
	tetrahedralize("p", &in, &out);
	out.save_nodes("boundary_counter");
	out.save_elements("boundary_counter");
	out.save_faces("boundary_counter");
	//*/
	cout<<"stlからのメッシュ作成完了"<<endl;
	/////
	//NODEpに、作成されたboundaryをコピーする　//すごく無駄な気がする。もう少し賢くやるにはNODEに各部品ごとの情報を出力してboundarydataをほかの関数でまとめる仕組みから変える必要がある
	//stlで読み込むと節点番号は1からスタートするため、従来の0から始まる手打ち座標に合わせるには1ずらす必要がある。このとき、麺情報も同様にずらさないといけない
	//in.initialize();
	//out.initialize();
	////

	//上で作ったメッシュは、座標系が基モデルとは異なるstlでつくったものなので、モデルの座標系で出力したものを参照し、あうようにずらす
	ifstream fp2("ea2_dis.dat");
	if(!fp2) cout<<"ea2.dat"<<endl;
	fp2.unsetf(ifstream::dec);
	fp2.setf(ifstream::skipws);
	int stuff=0;
	string b;
	//for(int i=0;i<3;i++) getline(fp2, b);//最初の3行分進める

	double x=0,y=0,z=0 ,x1=0,y1=0,z1=0, x2=0,y2=0,z2=0;
	
		fp2>>x1;
		fp2>>y1;
		fp2>>z1;
		fp2>>x2;
		fp2>>y2;
		fp2>>z2;
		
	cout<<x1<<" "<<y1<<" "<<z1<<endl;
	x=x1-x2; y=y1-y2; z=z1-z2;//真の値と現在の座標の差
	fp2.close();
	

	//node
	ifstream fp("boundary_counter.node");
	if(!fp) cout<<"boundary_counter"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//最初の3行分進める
	
	int node_num=0;
	stuff=0;
	fp>>node_num;
	for(int i=0;i<3;i++) fp>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<node_num;i++)//0からスタート
	{
		//temp.id=i;
		
		fp>>stuff;
		fp>>temp.r[A_X];
		fp>>temp.r[A_Y];
		fp>>temp.r[A_Z];
		
		/////
		temp.r[A_X]+=x;
		temp.r[A_Y]+=y;
		temp.r[A_Z]+=z;
		////
		//fp>>temp.r[A_Y];
		//fp>>temp.r[A_Z];
		NODEp.push_back(temp);
		temp.id+=1;
		//temp.id+=1;
		
	}
	fp.close();
	///

	//nodeファイル作成
	//MakeNodeFile(CON, NODEp, "NODEp.node");

	//in.load_node("NODEp");
	//in.load_stl("ea_electro2");
	//tetrahedralize("", &in, &out);
	///out.save_nodes("boundary_counter2");
	//out.save_elements("boundary_counter2");
	//out.save_faces("boundary_counter2");

	//境界面データ取得
	GetFacetListforCAD(FACEp, in, out, ELECTRODE2);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}

//放電電極境界面作成
void tetgen_function::SetDischargeElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"放電電極境界作成(cad読み込み)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	
	in.load_stl("ea1");//モデルの座標系で出力するとなぜか読み込んでくれない
	tetrahedralize("p", &in, &out);
	out.save_nodes("boundary_discharge");
	out.save_elements("boundary_discharge");
	out.save_faces("boundary_discharge");
	//*/
	//cout<<"stlからのメッシュ作成完了"<<endl;
	/////
	//NODEpに、作成されたboundaryをコピーする　//すごく無駄な気がする。もう少し賢くやるにはNODEに各部品ごとの情報を出力してboundarydataをほかの関数でまとめる仕組みから変える必要がある
	//stlで読み込むと節点番号は1からスタートするため、従来の0から始まる手打ち座標に合わせるには1ずらす必要がある。このとき、麺情報も同様にずらさないといけない
	//in.initialize();
	//out.initialize();
	////

	//上で作ったメッシュは、座標系が基モデルとは異なるstlでつくったものなので、モデルの座標系で出力したものを参照し、あうようにずらす
	ifstream fp2("ea1_dis.dat");
	if(!fp2) cout<<"ea1.dat"<<endl;
	fp2.unsetf(ifstream::dec);
	fp2.setf(ifstream::skipws);
	int stuff=0;
	string b;
	//for(int i=0;i<3;i++) getline(fp2, b);//最初の3行分進める

	double x=0,y=0,z=0 ,x1=0,y1=0,z1=0, x2=0,y2=0,z2=0;
	
		fp2>>x1;//本当の位置
		fp2>>y1;
		fp2>>z1;
		fp2>>x2;//座標系を合わせなかった位置(stlメッシュ作成位置)
		fp2>>y2;
		fp2>>z2;
		
	cout<<x1<<" "<<y1<<" "<<z1<<endl;
	x=x1-x2; y=y1-y2; z=z1-z2;//真の値と現在の座標の差
	fp2.close();
	

	//node
	ifstream fp("boundary_discharge.node");
	if(!fp) cout<<"boundary_discharge"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//最初の3行分進める
	
	int node_num=0;
	stuff=0;
	fp>>node_num;
	for(int i=0;i<3;i++) fp>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<node_num;i++)//0からスタート
	{
		//temp.id=i;
		
		fp>>stuff;
		fp>>temp.r[A_X];
		fp>>temp.r[A_Y];
		fp>>temp.r[A_Z];
		
		/////
		temp.r[A_X]+=x;
		temp.r[A_Y]+=y;
		temp.r[A_Z]+=z;
		////
		//fp>>temp.r[A_Y];
		//fp>>temp.r[A_Z];
		NODEp.push_back(temp);
		temp.id+=1;
		//temp.id+=1;
		
	}
	fp.close();
	///

	//nodeファイル作成
	//MakeNodeFile(CON, NODEp, "NODEp.node");

	//in.load_node("NODEp");
	//in.load_stl("ea_electro2");
	//tetrahedralize("", &in, &out);
	///out.save_nodes("boundary_counter2");
	//out.save_elements("boundary_counter2");
	//out.save_faces("boundary_counter2");

	//境界面データ取得
	GetFacetListforCAD(FACEp, in, out, ELECTRODE1);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//土台境界面作成
void tetgen_function::SetBaseBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEb, vector<tetgen_facet> &FACEb)
{
	cout<<"土台境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//分割数決定
	double	le = CON.get_distancebp();
	double	rc = TET.radius_column;
	double	h = TET.length_column;						//平板高さ
	double	dh = TET.fine_base;						//厚み方向メッシュ粗さ
	int		nh = int((TET.thickness_base+0.1*dh)/dh);	//厚み方向分割数
	double	L = TET.length_base;
	double	dL = TET.fine_base;						//xy方向メッシュ粗さ
	int		nL = int((TET.length_base/2+0.1*dL)/dL);	//xy方向分割数(片側)


	//上面接続部分
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L;
	NODEb.push_back(temp);
	temp.id+=1;

	rc*=1.005;//円柱に合わせてわずかに太くする
	for(double r=rc;r>le-0.1*le;r-=2*le)//leが2倍されていることに注意
	{
		int nr=int(2.0*PI*r/(2*le));//leが2倍されていることに注意
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L;
			NODEb.push_back(temp);
			temp.id+=1;
		}
	}//*/

	//上面表面
	double r=rc+2*le;
	double s=2*le;
	while(r<sqrt(2.0)*L/2+dL)
	{
		int nr=int(2.0*PI*r/s);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L;
			if(fabs(temp.r[A_X])<L/2 && fabs(temp.r[A_Y])<L/2)
			{
				NODEb.push_back(temp);
				temp.id+=1;
			}
		}
		r*=1.15;
		s*=1.15;
	}//*/

	//その他の面
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==nh)	//平板電極領域の端に来たとき節点を置く
				{
					//if(sqrt(temp.r[A_X]*temp.r[A_X]+temp.r[A_Y]*temp.r[A_Y])>rc+le || z>0)
					{
						temp.r[A_X]=x*dL;
						temp.r[A_Y]=y*dL;
						temp.r[A_Z]=-h-z*dL;
						NODEb.push_back(temp);
						temp.id+=1;
					}
				}
			}
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEb, "NODEb.node");

	in.load_node("NODEb");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_base");
	out.save_elements("boundary_base");
	out.save_faces("boundary_base");
	//*/

	//境界面データ取得
	GetFacetList(FACEb, in, out, ELECTRODE1);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//境界節点・境界面データの結合
void tetgen_function::UniteBoundaryData(mpsconfig &CON, 
					   vector<tetgen_node> &NODE, vector<tetgen_node> &NODEa1, vector<tetgen_node> &NODEa2, vector<tetgen_node> &NODEp, vector<tetgen_node> &NODEc, vector<tetgen_node> &NODEb, vector<tetgen_node> &NODEw, 
					   vector<tetgen_facet> &FACE, vector<tetgen_facet> &FACEa1, vector<tetgen_facet> &FACEa2, vector<tetgen_facet> &FACEp, vector<tetgen_facet> &FACEc, vector<tetgen_facet> &FACEb, vector<tetgen_facet> &FACEw, 
					   vector<int> &TRANS)
{
	cout<<"境界節点・境界面データの結合"<<endl;

	NODE.clear();
	FACE.clear();
	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=0;	//節点番号のオフセット量
	int offset_f=0;	//表面番号のオフセット量


	//水滴境界
	for(int i=0;i<(int)NODEw.size();i++)//節点
	{
		TRANS.push_back(NODEw[i].part_no);	//TRANSに粒子番号を格納(boundaryに粒子番号を格納している)

		temp_n=NODEw[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=WATER;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEw.size();i++)//表面
	{
		temp_f=FACEw[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=WATER;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEw.size();
	offset_f+=(int)FACEw.size();

	//空気境界
	for(int i=0;i<(int)NODEa1.size();i++)//節点
	{
		temp_n=NODEa1[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa1.size();i++)//表面
	{
		temp_f=FACEa1[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEa1.size();
	offset_f+=(int)FACEa1.size();

	//空気境界 高解像度面
	for(int i=0;i<(int)NODEa2.size();i++)//節点
	{
		temp_n=NODEa2[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa2.size();i++)//表面
	{
		temp_f=FACEa2[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEa2.size();
	offset_f+=(int)FACEa2.size();

	//平板電極境界
	for(int i=0;i<(int)NODEp.size();i++)//節点
	{
		temp_n=NODEp[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE2;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEp.size();i++)//表面
	{
		temp_f=FACEp[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE2;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEp.size();
	offset_f+=(int)FACEp.size();

	//円柱電極境界
	for(int i=0;i<(int)NODEc.size();i++)//節点
	{
		temp_n=NODEc[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEc.size();i++)//表面
	{
		temp_f=FACEc[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEc.size();
	offset_f+=(int)FACEc.size();

	//土台境界
	for(int i=0;i<(int)NODEb.size();i++)//節点
	{
		temp_n=NODEb[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEb.size();i++)//表面
	{
		temp_f=FACEb[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	//offset_n+=(int)NODEb.size();
	//offset_f+=(int)FACEb.size();

	//cout<<"----------OK"<<endl;
}

//境界節点・境界面データの追加
void tetgen_function::AddBoundaryData(mpsconfig &CON, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int attribute)
{
	//NODE,FACEに格納されている各部品のデータを，NODEall,FACEallに格納していく．

	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=(int)NODEall.size();	//節点番号のオフセット量
	int offset_f=(int)FACEall.size();	//表面番号のオフセット量


	//節点の追加
	for(int i=0;i<(int)NODE.size();i++)//節点
	{
		temp_n=NODE[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=attribute;
		NODEall.push_back(temp_n);
	}

	//面の追加
	for(int i=0;i<(int)FACE.size();i++)//表面
	{
		temp_f=FACE[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=attribute;
		FACEall.push_back(temp_f);
	}
}

//TRANS[]の格納
void tetgen_function::SetTRANS(vector<tetgen_node> &NODE, vector<int> &TRANS)
{
	//TRANS[i]には、節点番号iに対応する粒子番号を格納する。
	//FEM3D.cpp では、節点番号が1から始まるので、TRANS[0]には宣言後に-1を入れる。（ここでは既に入っている）

	for(int i=0;i<(int)NODE.size();i++)
	{
		TRANS.push_back(NODE[i].part_no);	//TRANSに粒子番号を格納
	}
}


//材質の決定（まだ未完成）
void tetgen_function::DecisionAttribute(vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double M[4];	//節点の材質を格納
	
	for(int i=0;i<(int)ELEM.size();i++)
	{	
		if(ELEM[i].attribute==AIR)
		{
			for(int n=0;n<4;n++)	M[n]=NODE[ELEM[i].node[n]].attribute;

			if(M[0]==AIR || M[1]==AIR || M[2]==AIR || M[3]==AIR)	//1つでも空気節点があれば空気要素
			{
				ELEM[i].attribute=AIR;
				out.tetrahedronattributelist[i]=AIR;
			}
			else
			{
				ELEM[i].attribute=WATER;	//残りは水
				out.tetrahedronattributelist[i]=WATER;
			}
		}
	}
}


//材質の修正  この関数を使うときは直前にtetgenioからデータを取得しておくこと
void tetgen_function::ModifyAttribute(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	//tetgenioからデータを取得済みと仮定する
	//節点-節点関係が得られていると仮定する

	/*double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;*/

	//要素の材質の修正  この時点では水要素となる領域は材質番号が0(未定義)となっている
	for(int i=0;i<(int)ELEM.size();i++)
	{
		//未定義(0)の材質を水にする
		if(ELEM[i].attribute==0)
		{
			ELEM[i].attribute=WATER;
		}
	}

	//節点の材質を要素の材質と合わせる

	//まず全部空気にする
	for(int i=0;i<(int)ELEM.size();i++)
	{
		NODE[i].attribute=AIR;
	}
	//水節点の決定
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}
	//電極節点の決定
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}

	//節点-節点関係より節点の修正
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].attribute==AIR)
		{
			int num_air=0;
			int num_ele=0;
		}
	}
}


//材質の修正  tetgenioを直接編集
void tetgen_function::ModifyAttribute_tetgenio(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;

	//節点自動挿入された節点の材質をとりあえず空気とする
	for(int i=0;i<out.numberofpoints;i++)
	{
		if(out.pointattributelist[i]==0)
		{
			out.pointattributelist[i]=AIR;
		}
	}//*/

	//要素の材質の修正  この時点では水要素となる領域は材質番号がデフォルトで0(未定義)となっている
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		//未定義の材質を水にする
		if(out.tetrahedronattributelist[i]==0)
		{
			out.tetrahedronattributelist[i]=WATER;//0だったものがWATERになる
		}
	}
	
	//静電無化の場合　流体-円柱間と円柱-土台間に空気の隙間ができているのでそれを埋める
	if(CON.get_model_number()==14)
	{
		//for(int times=0;times<30;times++)
		{
			for(int i=0;i<out.numberoftetrahedra;i++)
			{
				//空気要素を判定して変更する
				if(out.tetrahedronattributelist[i]==AIR)
				{
					int flag=0;

					//各節点の材質による判定
					int Mnode[4];	//節点の材質を格納
					int node_A=0;	//1要素における空気節点数
					int node_W=0;	//1要素における水節点数
					int node_E1=0;	//1要素における円柱電極節点数
					int node_E2=0;	//1要素における平板電極節点数
					for(int n=0;n<4;n++)	Mnode[n]=(int)out.pointattributelist[out.tetrahedronlist[4*i+n]];
					
					//隣り合う要素による判定
					int Melem[4];	//隣り合う要素の材質を格納
					int elem_A=0;	//隣り合う空気要素数
					int elem_W=0;	//隣り合う水要素数
					int elem_E1=0;	//隣り合う電極要素数
					for(int n=0;n<4;n++)	Melem[n]=(int)out.tetrahedronattributelist[out.neighborlist[4*i+n]];
					
					for(int n=0;n<4;n++)	//1値の要素の中で水と電極の節点数を数える
					{
						if(Mnode[n]==AIR)			node_A+=1;
						if(Mnode[n]==WATER)			node_W+=1;
						if(Mnode[n]==ELECTRODE1)	node_E1+=1;
						if(Mnode[n]==ELECTRODE2)	node_E2+=1;
						if(Melem[n]==AIR)			elem_A+=1;
						if(Melem[n]==WATER)			elem_W+=1;
						if(Melem[n]==ELECTRODE1)	elem_E1+=1;
					}

					if(node_W==4)	flag=1;							//全てが水節点のときは水
					if(node_A==0 && node_W>0 && node_E1>0)	flag=1;	//空気節点が0で水節点と電極節点が1以上ときは水
					if(elem_A<=1)	flag=1;							//隣接する空気要素が1以下のときは水
					if(elem_W>0 && elem_E1>0)	flag=1;				//隣接する水要素と電極要素が1以上のときは水

					//if(node_E2==4)	flag=2;	//全てが電極節点のときは電極  これをすると円柱と土台の境にわずかな斜めの面ができてしまう

					//if((M[0]==AIR || M[0]==WATER) && M[1]!=AIR && M[2]!=AIR && M[3]!=AIR)	flag=1;//-------------①
					//if(M[0]==M[1] && M[1]==M[2] && M[2]==M[3] && M[0]==WATER)		flag=1;
					//*/

					//重心座標による判定(節点材質による判定でこぼれたものを拾う)
					double r[3]={0,0,0};
					for(int d=0;d<3;d++)
					{
						for(int n=0;n<4;n++)	r[d]+=out.pointlist[3*(out.tetrahedronlist[4*i+n])+d];
						r[d]/=4;
					}
					//if(sqrt(r[A_X]*r[A_X]+r[A_Y]*r[A_Y])<rc-le && r[A_Z]<le && r[A_Z]>0)	flag=1;	//電極の縁に針みたいな流体要素ができることがあるのでleだけ短くしておく  判定漏れ次の材質判定で拾ってもらうことに期待
					if(sqrt(r[A_X]*r[A_X]+r[A_Y]*r[A_Y])<rc && r[A_Z]>-L && r[A_Z]<-L+le)	flag=2;	//土台と円柱の境目は座標による判定でも比較的簡単
					//*/

					if(flag==1)	//フラグが1なら水
					{
						out.tetrahedronattributelist[i]=WATER;
					}
					if(flag==2)	//フラグが2なら円柱電極
					{
						out.tetrahedronattributelist[i]=ELECTRODE1;
					}
				}
			}
		}
	}

	//自動追加された節点ののattributeとboundary_markerの修正
	//この時点では流体内部や電極内部や電極内部の節点は空気節点になっているのでそれぞれの材質に修正する
	//※電極の一番上の節点は電極の節点とするため、先に水の節点決めてから電極節点を決める
	
	//水要素を構成する節点は水節点に
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		if(out.tetrahedronattributelist[i]==WATER)
		{
			for(int n=0;n<4;n++)
			{
				out.pointattributelist[out.tetrahedronlist[i*4+n]]=WATER;
				out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=WATER;
			}
		}
		else if(out.tetrahedronattributelist[i]==CONDUCT)
		{
			for(int n=0;n<4;n++)
			{
				out.pointattributelist[out.tetrahedronlist[i*4+n]]=CONDUCT;
				out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=CONDUCT;
			}
		}
		else if(out.tetrahedronattributelist[i]==AIR)
		{
			for(int n=0;n<4;n++)
			{
				out.pointattributelist[out.tetrahedronlist[i*4+n]]=AIR;
				out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=AIR;
			}
		}
		else if(out.tetrahedronattributelist[i]==COIL)
		{
			for(int n=0;n<4;n++)
			{
				out.pointattributelist[out.tetrahedronlist[i*4+n]]=COIL;
				out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=COIL;
			}
		}
	}

	
	
	//静電霧化　平板電極および円柱電極(boundary_markerについては平板と円柱で分ける)
	if(CON.get_model_number()==14)
	{
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			//円柱
			if(out.tetrahedronattributelist[i]==ELECTRODE1)
			{
				for(int n=0;n<4;n++)
				{
					out.pointattributelist[out.tetrahedronlist[i*4+n]]=ELECTRODE;
					out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=ELECTRODE1;
				}
				//材質はELECTRODEに戻しておく
				out.tetrahedronattributelist[i]=ELECTRODE;
			}
			//平板
			if(out.tetrahedronattributelist[i]==ELECTRODE2)
			{
				for(int n=0;n<4;n++)
				{
					out.pointattributelist[out.tetrahedronlist[i*4+n]]=ELECTRODE;
					out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=ELECTRODE2;
				}
				//材質はELECTRODEに戻しておく
				out.tetrahedronattributelist[i]=ELECTRODE;
			}
		}
	}

}


//要素の細分化
void tetgen_function::FineElement(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	tetgenio add;
	add.initialize();

	vector<tetgen_node> NODEadd;
	tetgen_node temp;
	temp.id=0;
	temp.attribute=0;
	temp.boundary=0;

	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		if(out.tetrahedronattributelist[i]==AIR)
		{
			double r[3]={0,0,0};
			//for(int d=0;d<3;d++)	r[d]=out.vpointlist[i*3+d];
			for(int d=0;d<3;d++)
			{
				for(int n=0;n<4;n++)	r[d]+=out.pointlist[3*(out.tetrahedronlist[4*i+n])+d];
				r[d]/=4;
			}

			if(fabs(r[A_X])<0.0005 && fabs(r[A_Y])<0.0005 && fabs(r[A_Z])<0.0005)
			{
				temp.r[A_X]=r[A_X];
				temp.r[A_Y]=r[A_Y];	
				temp.r[A_Z]=r[A_Z];
				NODEadd.push_back(temp);
				temp.id+=1;
			}
		}
	}

	MakeNodeFile(CON, NODEadd, "output-a.node");
}//*/

//節点間距離計算関数
double tetgen_function::Distance(tetgen_node &point1, tetgen_node &point2)
{
	double dis=0;

	for(int d=0;d<3;d++)	dis+=(point2.r[d]-point1.r[d])*(point2.r[d]-point1.r[d]);

	return sqrt(dis);
}


//要素重心座標計算関数
void tetgen_function::CalcBarycentricElement(mpsconfig&, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	double r[3]={0,0,0};

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int d=0;d<3;d++)	for(int n=0;n<4;n++)	r[d]+=NODE[ELEM[i].node[n]].r[d];
		for(int d=0;d<3;d++)	ELEM[i].g[d]=r[d]/4;
	}
}