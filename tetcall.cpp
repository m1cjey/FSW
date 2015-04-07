////TetGen�Ăяo���֐�  TetGen�̓����

#include "header.h"

#define FULL 3
#define REMESH 4
#define FULL_INPORT 5


//TetGen����  ���f���ɂ�蕪��
void tetgen_function::call_TetGen(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	cout<<"TetGen�ɂ�郁�b�V�������J�n"<<endl;
	clock_t t1=clock();	//�N���b�N���擾


	//TetGen_test(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//test
	if(CON.get_mesh_input()==0)
	{
		//////���f�����ɕ���/////////////////////////////
		if(CON.get_model_number()==14)
		{
			if(CON.get_ea_model()==0) TetGen_nanoe(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//�Ód�����p
			else if(CON.get_ea_model()==1) TetGen_nanoe2(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//�Ód�����p,cad�ǂݍ���
		}
		if(CON.get_model_number()==15)	TetGen_rayleigh(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//���C���[�����p
		if(CON.get_model_number()==20 || CON.get_model_number()==21)	TetGen_CC(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//�R�[���h�N���[�V�u���p
		if(CON.get_model_number()==25)	TetGen_eddytest(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	
		if(CON.get_model_number()==30)	TetGen_nanoe2(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);	//�Ód�����p,cad�ǂݍ���
		/////////////////////////////////////////////////
	}
	if(CON.get_mesh_input()==1)
	{
		if(CON.get_model_number()==20 || CON.get_model_number()==21) TetGen_inport(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
	}

	clock_t t2=clock();	//�N���b�N���擾
	cout<<"���b�V����������  CPU time="<<(t2-t1)/CLOCKS_PER_SEC<<endl;
}


//////�����艺�ɁA��肽�����f���̃��b�V�������v���O�������L�q���Ă�������/////////////////////////////////////////////

//�Ód�����p  ���b�V������
void tetgen_function::TetGen_nanoe(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////���E�ߓ_�Ƌ��E�ʂ̃f�[�^�̍쐬////////////////
	/*-----------------------------------------------------------
	Set******Boundary�֐��ɂāC���i���Ƃ�NODE��FACE�ɋ��E�ߓ_�Ƌ��E�ʂ��i�[���C
	AddBoundaryData�֐��ɂāCNODE��FACE�̃f�[�^��NODEall��FACEall�Ɋi�[���Ă����D
	-----------------------------------------------------------*/

	////////////���H�̈�////////////
	NODE.clear();
	FACE.clear();
	//SetFluidBoundary(CON, PART, TET, NODE, FACE);	//���̋��E�̐ݒ�
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�O���_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�����_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetTRANS(NODE, TRANS);	//�ߓ_�ԍ��Ɨ��q�ԍ��̑Ή��t��
	//AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////��C�̈�////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////���H�\�ʕt�ߒǉ��ߓ_////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	
	////////////���d��////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////�~���d��////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////�d�ɓy��////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	


	//.node�t�@�C���쐬
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"poly�t�@�C���ǂݍ���"<<endl;
	in.load_poly("all_boundary");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	cout<<"���������J�n"<<endl;
	tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	out.save_nodes("output");
	out.save_elements("output");

	//�ގ��C��
	cout<<"�ގ��C��"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW

	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//�Ód�����p  ���b�V������// cad����ǂݍ���
void tetgen_function::TetGen_nanoe2(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////���E�ߓ_�Ƌ��E�ʂ̃f�[�^�̍쐬////////////////
	/*-----------------------------------------------------------
	Set******Boundary�֐��ɂāC���i���Ƃ�NODE��FACE�ɋ��E�ߓ_�Ƌ��E�ʂ��i�[���C
	AddBoundaryData�֐��ɂāCNODE��FACE�̃f�[�^��NODEall��FACEall�Ɋi�[���Ă����D
	-----------------------------------------------------------*/

	/*///////////���H�̈�////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//���̋��E�̐ݒ�
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�O���_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�����_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	SetTRANS(NODE, TRANS);	//�ߓ_�ԍ��Ɨ��q�ԍ��̑Ή��t��
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	///////////*/
	
	////////////��C�̈�////////////
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

	////////////�Ό��d��////////////
	NODE.clear();
	FACE.clear();
	//SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	SetCounterElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);
	

	////////////���d�d��////////////
	NODE.clear();
	FACE.clear();
	SetDischargeElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	/*///////////�d�ɓy��////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	//////*/
	


	//.node�t�@�C���쐬
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"poly�t�@�C���ǂݍ���"<<endl;
	in.load_poly("all_boundary");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	cout<<"���������J�n"<<endl;
	tetrahedralize("pqa1.67e-7AYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	//tetrahedralize("pq1.2An", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	out.save_nodes("output");
	out.save_elements("output");

	//�ގ��C��
	cout<<"�ގ��C��"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW

	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}


//ڲذ�����p  ���b�V������
void tetgen_function::TetGen_rayleigh(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////���H�̈�////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);
	SetTRANS(NODE, TRANS);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////��C�̈�////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////��C�̈�ǉ��ߓ_////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	//*/


	//.node�t�@�C���쐬
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"poly�t�@�C���ǂݍ���"<<endl;
	in.load_poly("all_boundary");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	cout<<"���������J�n"<<endl;
	tetrahedralize("pq1.3a1.67e-7AYYn", &in, &out);	//1.1�����ł͐؂�Ȃ� �f�t�H���g��rqa1.1AYYn
	out.save_nodes("output");
	out.save_elements("output");
	//out.save_faces("output");
	//out.save_neighbors("output");
	//*/

	cout<<"�ގ��C��"<<endl;
	//ModifyAttribute(CON, TET, NODEall, ELEMall, in, out);
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW

	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//CC�p  ���b�V������
void tetgen_function::TetGen_CC(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////���E�ߓ_�Ƌ��E�ʂ̃f�[�^�̍쐬////////////////
	/*-----------------------------------------------------------
	Set******Boundary�֐��ɂāC���i���Ƃ�NODE��FACE�ɋ��E�ߓ_�Ƌ��E�ʂ��i�[���C
	AddBoundaryData�֐��ɂāCNODE��FACE�̃f�[�^��NODEall��FACEall�Ɋi�[���Ă����D
	-----------------------------------------------------------*/

	////////////���H�̈�////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//���̋��E�̐ݒ�
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�O���_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�����_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	SetTRANS(NODE, TRANS);	//�ߓ_�ԍ��Ɨ��q�ԍ��̑Ή��t��
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////��C�̈�////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////�R�C��////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	/*///////////���////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////���H�\�ʕt�ߒǉ��ߓ_////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////���d��////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////�~���d��////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////�d�ɓy��////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.node�t�@�C���쐬
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"poly�t�@�C���ǂݍ���"<<endl;
	in.load_poly("all_boundary");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	cout<<"���������J�n"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	tetrahedralize("pq1.2a1.67e-6AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	out.save_nodes("output");
	out.save_elements("output");

	//�ގ��C��
	cout<<"�ގ��C��"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW
			
	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//  �Q�d�����p���b�V������
void tetgen_function::TetGen_eddytest(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////���E�ߓ_�Ƌ��E�ʂ̃f�[�^�̍쐬////////////////
	/*-----------------------------------------------------------
	Set******Boundary�֐��ɂāC���i���Ƃ�NODE��FACE�ɋ��E�ߓ_�Ƌ��E�ʂ��i�[���C
	AddBoundaryData�֐��ɂāCNODE��FACE�̃f�[�^��NODEall��FACEall�Ɋi�[���Ă����D
	-----------------------------------------------------------*/

	////////////���H�̈�//////////// ���̖��ł͐��H�͂���Ȃ�
	NODE.clear();
	FACE.clear();
	//SetFluidBoundary(CON, PART, TET, NODE, FACE);	//���̋��E�̐ݒ�
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�O���_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�����_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetTRANS(NODE, TRANS);	//�ߓ_�ԍ��Ɨ��q�ԍ��̑Ή��t��
	//AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	/*///////////��C�̈�////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	//////*/

	////////////��C�̈�2////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary2(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	/*///////////�R�C��////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	////////////���d�̕�////////////
	NODE.clear();
	FACE.clear();
	SetConductBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CONDUCT);
	////////*/


	/*///////////���////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////���H�\�ʕt�ߒǉ��ߓ_////////////
	NODE.clear();
	FACE.clear();
	//SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	//AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////���d��////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////�~���d��////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////�d�ɓy��////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.node�t�@�C���쐬
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"poly�t�@�C���ǂݍ���"<<endl;
	in.load_poly("all_boundary");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	cout<<"���������J�n"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	//tetrahedralize("pq1.2a1.67e-6AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	tetrahedralize("pq1.5Aa1e-3YYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	out.save_nodes("output");
	out.save_elements("output");

	//�ގ��C��
	cout<<"�ގ��C��"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW
			
	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//�O���쐬�f�[�^����ǂݍ���Ń��b�V���쐬
void tetgen_function::TetGen_inport(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	cout<<"magnet�f�[�^��ǂݍ���tetgen���b�V���쐬"<<endl;
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out,addin;
	in.initialize();
	out.initialize();

	//tetgenio addin;
	addin.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////���E�ߓ_�Ƌ��E�ʂ̃f�[�^�̍쐬////////////////
	/*-----------------------------------------------------------
	Set******Boundary�֐��ɂāC���i���Ƃ�NODE��FACE�ɋ��E�ߓ_�Ƌ��E�ʂ��i�[���C
	AddBoundaryData�֐��ɂāCNODE��FACE�̃f�[�^��NODEall��FACEall�Ɋi�[���Ă����D
	-----------------------------------------------------------*/

	/*///////////���H�̈�////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//���̋��E�̐ݒ�
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�O���_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�����_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	SetTRANS(NODE, TRANS);	//�ߓ_�ԍ��Ɨ��q�ԍ��̑Ή��t��
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////��C�̈�////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////�R�C��////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	/*///////////���////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////���H�\�ʕt�ߒǉ��ߓ_////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////���d��////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////�~���d��////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////�d�ɓy��////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.node�t�@�C���쐬
	//MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	//MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"static�t�@�C���ǂݍ���"<<endl;
	//in.load_node("static");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	in.load_tetmesh("static");	//.node,.ele��Ǎ�
	//addin.load_node("fluid");
	cout<<"�Ǎ����b�V���̍č\���J�n"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	tetrahedralize("r", &in, &out);	//�ǂݍ���static����tetgen�p�ɍč\��
	out.save_nodes("output");
	out.save_elements("output");
	out.save_faces("output");


	////////////���H�̈�////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//���̋��E�̐ݒ�
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�O���_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�����_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	SetTRANS(NODE, TRANS);	//�ߓ_�ԍ��Ɨ��q�ԍ��̑Ή��t��
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	//in.load_node("fluid");	//.node,.ele��Ǎ�
	cout<<"���̐ߓ_�̒ǉ��J�n"<<endl;
	tetrahedralize("i", &in, &out);
	

	//�ގ��C��
	//cout<<"�ގ��C��"<<endl;
	//ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW
	
	
			
	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}

//CC�p  ���b�V������
void tetgen_function::TetGen_test(mpsconfig &CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenio�N���X�錾�Ə�����
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen�ݒ�N���X
	tetgen_config TET;

	//�S�̔z�񏉊���
	NODEall.clear();
	FACEall.clear();
	//�e���i�p�z��i�g���񂵁j
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////////���E�ߓ_�Ƌ��E�ʂ̃f�[�^�̍쐬////////////////
	/*-----------------------------------------------------------
	Set******Boundary�֐��ɂāC���i���Ƃ�NODE��FACE�ɋ��E�ߓ_�Ƌ��E�ʂ��i�[���C
	AddBoundaryData�֐��ɂāCNODE��FACE�̃f�[�^��NODEall��FACEall�Ɋi�[���Ă����D
	-----------------------------------------------------------*/

	////////////���H�̈�////////////
	NODE.clear();
	FACE.clear();
	SetFluidBoundary(CON, PART, TET, NODE, FACE);	//���̋��E�̐ݒ�
	//SetFluidBoundary_OutsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�O���_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	//SetFluidBoundary_InsideDummyNode(CON, PART, TET, NODEw, FACEw);	//�����_�~�[�ߓ_���g�������H���b�V�������i�o�O����ꂸ�ɖ������j
	SetTRANS(NODE, TRANS);	//�ߓ_�ԍ��Ɨ��q�ԍ��̑Ή��t��
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, WATER);
	
	////////////��C�̈�////////////
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

	////////////�R�C��////////////
	NODE.clear();
	FACE.clear();
	SetCoilBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////*/

	/*///////////���////////////
	NODE.clear();
	FACE.clear();
	SetCrucibleBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, CRUCIBLE);
	///////*/

	/*///////////���H�\�ʕt�ߒǉ��ߓ_////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	/////////////

	

	////////////���d��////////////
	NODE.clear();
	FACE.clear();
	SetPlateElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE2);

	////////////�~���d��////////////
	NODE.clear();
	FACE.clear();
	SetColumnElectrodeBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);

	////////////�d�ɓy��////////////
	NODE.clear();
	FACE.clear();
	SetBaseBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, ELECTRODE1);
	
	//////*/

	//.node�t�@�C���쐬
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.poly�t�@�C���̏o��
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"poly�t�@�C���ǂݍ���"<<endl;
	in.load_poly("all_boundary");	//.poly��ǂݍ��񂾂玩���I��.node���ǂݍ��ނ炵��
	cout<<"���������J�n"<<endl;
	//tetrahedralize("pq1.2a1.67e-7AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	tetrahedralize("pq1.2a1.67e-6AYYn", &in, &out);	//1.1�ȉ��ł͐؂�Ȃ�
	out.save_nodes("output");
	out.save_elements("output");

	//�ގ��C��
	cout<<"�ގ��C��"<<endl;
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenio�𒼐ڕҏW
			
	//�o��
	cout<<"TetView�p�t�@�C���o��"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//�ߓ_�E�v�f�f�[�^�擾
	cout<<"TetGen���ߓ_�E�v�f�f�[�^�擾"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);

}