#include "header.h"


//.node�f�[�^�擾�֐�
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


//.ele�f�[�^�擾�֐�(�ȈՔ�)
void tetgen_function::GetTetrahedronList(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];
		temp.attribute=0;	//���炭�A����`�ő������ƕςȐ��l�������ăo�O��̂ł����ł�0�ɂ��Ƃ�

		ELEM.push_back(temp);
	}
}


//.ele�f�[�^�擾�֐�(�ގ��E�v�f�v�f�֌W�܂�)
void tetgen_function::GetTetrahedronList_full(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;											//�ߓ_�ԍ�
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];	//�\���ߓ_
		for(int n=0;n<4;n++)	temp.nei_elem[n]=out.neighborlist[i*4+n];	//�v�f-�v�f�֌W
		temp.attribute=(int)out.tetrahedronattributelist[i];				//�ގ�
		//temp.volume=out.tetrahedronvolumelist[i];							//�̐�(�擾�ł��Ȃ�)

		ELEM.push_back(temp);
	}
}


//.face�f�[�^�擾�֐�
void tetgen_function::GetFacetList(vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out, int boundary)
{
	//��boundarymarker�͈����Ƃ��ė^���Ă���BPLC���[�h�ō�������b�V���ł͂Ȃ����ߋ��E���o�͂���Ȃ��B

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

//.face�f�[�^�擾�֐� cad(stl�ǂݍ���)�p�B1����n�܂�node�ԍ���
void tetgen_function::GetFacetListforCAD(vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out, int boundary)
{
	//��boundarymarker�͈����Ƃ��ė^���Ă���BPLC���[�h�ō�������b�V���ł͂Ȃ����ߋ��E���o�͂���Ȃ��B

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


//.node�t�@�C���쐬�֐�
void tetgen_function::MakeNodeFile(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	//cout<<filename<<" �o��"<<endl;

	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	}

	fout.close();
}


//.node�t�@�C���쐬�֐�
void tetgen_function::MakeNodeFile_NonAttributeAndBoundary(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	cout<<filename<<" �o��"<<endl;

	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}

	fout.close();
}


//.ele�t�@�C���쐬�֐�
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


//.face�t�@�C���쐬�֐�
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


//.poly�t�@�C���쐬�֐�
void tetgen_function::MakePolyFile(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, char *filename)
{
	//cout<<filename<<" �o��"<<endl;

	ofstream fout(filename);

	//node list (�����ł͏o�͂��Ȃ�)
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

	////////////////////////////region attribute�̌��� (�z��Ɋi�[���Ă���o�͂���)/////////////////////////////
	//�ގ��̎w��́C���E���ɂ����_�̍��W�����߁C�����̍ގ����w�肷�邱�ƂŁC�������E���ɂ���v�f���S�Ă��̍ގ��ɂȂ�D
	//���͕���𔺂��C�ގ��̎w�肪����ł��邽�߁C�����ł͍s��Ȃ��D
	//��̍ގ��̏C���ɂ����āC����`�ƂȂ��Ă���v�f�𐅗v�f�Ƃ���D

	vector<region_attribute_list> REGION;
	region_attribute_list temp;
	temp.id=0;
	temp.region_number=0;
	temp.region_attribute=0;
	
	//�Ód����
	if(CON.get_model_number()==14)
	{
		if(CON.get_ea_model()==0)//�~���d�Ƀ��f��
		{
			//��C
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=(CON.get_ZU()+TET.height_plate+TET.thickness_plate)/2;	//��͗̈�̏�[�ƕ��d�ɂ̏�ʂ̒��ԓ_
			temp.region_number=AIR;
			temp.region_attribute=AIR;
			REGION.push_back(temp);
		
			//����
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=TET.height_plate+TET.thickness_plate/2.0;	//���̌��ݕ����̒��ԓ_
			temp.region_number=ELECTRODE2;
			temp.region_attribute=ELECTRODE2;
			REGION.push_back(temp);

			//�~��
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=-TET.length_column/2.0;		//�~�����������̒��ԓ_
			temp.region_number=ELECTRODE1;
			temp.region_attribute=ELECTRODE1;
			REGION.push_back(temp);

			//�y��
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=-TET.length_column-TET.thickness_base/2.0;		//�y��̌��ݕ����̒��ԓ_
			temp.region_number=ELECTRODE1;
			temp.region_attribute=ELECTRODE1;
			REGION.push_back(temp);
		}
		if(CON.get_ea_model()==1)//���i�d��(cad�x�[�X���f��)
		{
			//��C
			temp.id+=1;
			temp.r[A_X]=0;
			temp.r[A_Y]=0;
			temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�̏�[�ƕ��d�ɂ̏�ʂ̒��ԓ_
			temp.region_number=AIR;
			temp.region_attribute=AIR;
			REGION.push_back(temp);

			//�Ό��d��
			temp.id+=1;
			temp.r[A_X]=0.0016;
			temp.r[A_Y]=0;
			temp.r[A_Z]=0.0029;	//���̌��ݕ����̒��ԓ_//���_�����z��������3mm�A���a3mm~4mm
			temp.region_number=ELECTRODE2;
			temp.region_attribute=ELECTRODE2;
			REGION.push_back(temp);

			//���d�d��
			temp.id+=1;
			temp.r[A_X]=0.00;
			temp.r[A_Y]=0;
			temp.r[A_Z]=-0.001;	//���̌��ݕ����̒��ԓ_//���_�����z��������3mm�A���a3mm~4mm
			temp.region_number=ELECTRODE1;
			temp.region_attribute=ELECTRODE1;
			REGION.push_back(temp);
		}
	}

	//���C���[����
	if(CON.get_model_number()==15)
	{
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);
	}

	//CC
	if(CON.get_model_number()==20 ||CON.get_model_number()==21)
	{
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		/*//���
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0.053;
		temp.r[A_Z]=0.025;
		temp.region_number=CRUCIBLE;
		temp.region_attribute=CRUCIBLE;
		REGION.push_back(temp);
		////*/

		//�R�C��
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0.053;
		temp.r[A_Z]=0.025;
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);
		///*/
	}

	//�Q�d���p�ȈՃ��f��
	if(CON.get_model_number()==25)
	{
		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//��C
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=TET.thickness_conduct*0.501;	//��͗̈�����9���̂Ƃ���
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//���̕�
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0.0;
		temp.r[A_Z]=0.0;
		temp.region_number=CONDUCT;
		temp.region_attribute=CONDUCT;
		REGION.push_back(temp);
		/////

		/*////�R�C��
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


//.smesh�t�@�C���쐬�֐�
void tetgen_function::MakeSmeshFile(mpsconfig &CON, vector<tetgen_facet> &FACE, char *filename)
{
	double le=CON.get_distancebp();

	ofstream fsmesh(filename);

	//node list (�����ł͏o�͂��Ȃ�)
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


//���̋��E�ʍ쐬
void tetgen_function::SetFluidBoundary(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEw, vector<tetgen_facet> &FACEw)
{
	cout<<"���̋��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMw;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;
	temp.attribute=WATER;
	//temp.boundary=WATER;

	double le=CON.get_distancebp();	//���q�ԋ���
	double rc=TET.radius_column;	//�d�ɔ��a
	int type;
	int part_no;


	//���q�f�[�^�̎�荞��
	
	//�Ód����
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
	//���C���[����
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

			if(PART[i].type==FLUID)//�Ǘ��q�̓��b�V���쐬�Ɋ܂߂Ȃ�
			{
				trans.push_back(part_no);
				NODEw.push_back(temp);
				temp.id+=1;
			}
		}
	}

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEw, "NODEw.node");

	//.node�t�@�C���ǂݎ��
	in.load_node("NODEw");

	//�܂��͗��̐ߓ_�݂̂ŕ���
	tetrahedralize("", &in, &out);
		// i �v�f���֐ߓ_��ǉ�(���̏������Ă��@�\���Ă���)
		// f .face�t�@�C���ɋ��E�ł͂Ȃ��ʂ��܂߂�
		// e .edge�t�@�C���̏o��(ON�ɂ���ƂȂ����~�܂��Ă��܂�)
		// n .neigh�t�@�C���̏o��

	//�o��
	out.save_nodes("fluid_whole");
	out.save_elements("fluid_whole");

	//////////////////�����܂łŗ��̐ߓ_�݂̂��g���āA���ׂĂ̗v�f���q�������ʂȃ��b�V�����ł���(fluid_whole�Ŋm�F�\)*/


	///////////////�s�v�ȗv�f�̍폜

	//.node�̎擾
	GetPointList(NODEw, in, out);
	//.ele�̎擾
	GetTetrahedronList(ELEMw, in, out);

	//�����v�f�̏���
	DelThinTetrahedron(CON, TET, NODEw, ELEMw, in, out);

	//�ߓ_-�v�f�֌W
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W��藬�̕\�ʎ擾
	GetFacetList_from_neigh(CON, ELEMw, FACEw);

	//��яo��������v�f�̍폜
	DelThinTetrahedron_SharpElem(CON, NODEw, ELEMw, FACEw, in, out);

	/////��яo���v�f���폜�����̂ŁC������x�C�ߓ_�Ɨv�f�̗אڊ֌W�����ߒ���
	//�ߓ_-�v�f�֌W
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W��藬�̕\�ʎ擾
	GetFacetList_from_neigh(CON, ELEMw, FACEw);

	//���q�ԍ���NODEw���ɑ��
	for(int i=0;i<NODEw.size();i++)	NODEw[i].part_no=trans[i];
	
	//�\�ʂ��\������ߓ_��I�����C�z��ԍ����l�߂� ���̓��������q�̐ߓ_���g���ꍇ�̓R�����g�A�E�g
	//SelectFaceNode(CON, NODEw, FACEw);

	/////////////////�v�f�m�F�p�t�@�C��///////////////////////////////////
	out.save_nodes("boundary_fluid");	//fluid.2.node�Ɠ����t�@�C��
	MakeElemFile(CON, ELEMw, "boundary_fluid.ele");
	MakeFaceFile(CON, FACEw, "boundary_fluid.face");
	////////////////�����܂łŐ��H�̃��b�V�����؂ꂽ//////////////////////*/
}


//���̋��E�ʍ쐬  �O���@�������_�~�[�ߓ_�@
void tetgen_function::SetFluidBoundary_OutsideDummyNode(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEw, vector<tetgen_facet> &FACEw)
{
	cout<<"���̋��E�쐬(�O���@���_�~�[�ߓ_�@)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMw;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;
	temp.attribute=WATER;
	temp.boundary=WATER;
	
	double le=CON.get_distancebp();	//���q�ԋ���
	double rc=TET.radius_column;	//�d�ɔ��a
	//double dummy;
	//int part_no;
	//int type;

	int num_dummy=0;	//�_�~�[�ߓ_��
	double maxZ=0;


	//���̐ߓ_
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

			//Z���W�̍ő�l���X�V
			if(maxZ<temp.r[A_Z])	maxZ=temp.r[A_Z];
		}
	}
	//cout <<minZ<<endl;
	//Z���W�̍ŏ��l������������
	//if(minZ<0)	minZ=0.1*le;	//0.1*le�Ƃ��Ă���
	//cout <<minZ<<endl;

	//�@�������_�~�[�ߓ_
	temp.attribute=AIR;
	temp.boundary=AIR;

	//�_�~�[�ߓ_////////////////////////////////////////////
	//�@���O��
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
				trans.push_back(-1);//���q�ɑΉ����Ȃ��ߓ_��-1���i�[
				num_dummy+=1;
			}
		}
	}//*/

	//�d�ɂ̏�[�ʁi���̂̉��j
	for(double r=le;r<rc+3*le;r+=le)//�L�߂ɒu��
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
			trans.push_back(-1);//���q�ɑΉ����Ȃ��ߓ_��-1���i�[
			num_dummy+=1;
		}
	}//*/

	/*//�X�[�p�[�{�b�N�X////////////////////////////////////////////
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
				trans.push_back(-1);//���q�ɑΉ����Ȃ��ߓ_��-1���i�[
				num_dummy+=1;
			}
		}
	}//*/

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEw, "fluid.1.node");

	//.node�t�@�C���ǂݎ��
	in.load_node("fluid.1");

	//����
	tetrahedralize("", &in, &out);
		// i �v�f���֐ߓ_��ǉ�(���̏������Ă��@�\���Ă���)
		// f .face�t�@�C���ɋ��E�ł͂Ȃ��ʂ��܂߂�
		// e .edge�t�@�C���̏o��(ON�ɂ���ƂȂ����~�܂��Ă��܂�)
		// n .neigh�t�@�C���̏o��

	//�o��
	out.save_nodes("fluid.2");
	out.save_elements("fluid.2");

	//////////////////�����܂łŗ��̕\�ʂƂ��̊O���̃_�~�[�ߓ_���g�����ʂȌ`��̃��b�V�����ł���(fluid.2�Ŋm�F�\)*/


	///////////////�s�v�ȗv�f�̍폜

	//.node�̎擾
	GetPointList(NODEw, in, out);
	//.ele�̎擾
	GetTetrahedronList(ELEMw, in, out);

	//�_�~�[�v�f�̏���
	DelTetrahedron_OutsideDummy(CON, NODEw, ELEMw, in, out);

	//�ߓ_-�v�f�֌W
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W��藬�̕\�ʎ擾
	GetFacetList_from_neigh(CON, ELEMw, FACEw);

	//�_�~�[�ߓ_�̍폜
	DelDummyNode(CON, NODEw, FACEw, num_dummy);
	//�\�ʂ��\������ߓ_��I�����C�z��ԍ����l�߂� ���̓��������q�̐ߓ_���g���ꍇ�̓R�����g�A�E�g
	//SelectFaceNode(CON, NODEw, FACEw);

	//�e�X�g
	//MakeNodeFile(CON, NODEw, "NODEw.node");
	//MakeFaceFile(CON, FACEw, "FACEw.face");


	/////////////////�v�f�m�F�p�t�@�C��///////////////////////////////////
	//out.save_nodes("boundary_fluid");	//fluid.2.node�Ɠ����t�@�C��
	MakeNodeFile_NonAttributeAndBoundary(CON, NODEw, "boundary_fluid.node");
	MakeElemFile(CON, ELEMw, "boundary_fluid.ele");
	MakeFaceFile(CON, FACEw, "boundary_fluid.face");
	////////////////�����܂łŐ��H�̃��b�V�����؂ꂽ//////////////////////*/
}


//���̋��E�ʍ쐬  �����@�������_�~�[�ߓ_�@(�J���r��)
void tetgen_function::SetFluidBoundary_InsideDummyNode(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEw, vector<tetgen_facet> &FACEw)
{
	cout<<"���̋��E�쐬(�������̖@���_�~�[�ߓ_�@)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMw;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;
	temp.attribute=WATER;
	temp.boundary=WATER;

	double le=CON.get_distancebp();	//���q�ԋ���
	double rc=TET.radius_column;	//�d�ɔ��a
	//int type;
	//int part_no;
	double minZ=le;


	//���q�f�[�^�̎�荞��
	for(int i=0;i<(int)PART.size();i++)
	{		
		if(PART[i].type==BOFLUID && PART[i].r[A_Z]>le*0.5)
		{
			temp.boundary=PART[i].type;
			for(int d=0;d<3;d++)	temp.r[d]=PART[i].r[d];
			trans.push_back(i);
			NODEw.push_back(temp);
			temp.id+=1;

			//Z���W�̍ŏ��l���X�V
			if(minZ>temp.r[A_Z])	minZ=temp.r[A_Z];
		}
	}
	//cout <<minZ<<endl;
	//Z���W�̍ŏ��l������������
	if(minZ<0)	minZ=0.1*le;	//0.1*le�Ƃ��Ă���
	//cout <<minZ<<endl;

	//�����Ƀ_�~�[�ߓ_��u��
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
				trans.push_back(-1);//���q�ɑΉ����Ȃ��ߓ_��-1���i�[
			}
		}
		if(PART[i].type==INWALL/* && sqrt(temp.r[A_X]*temp.r[A_X]+temp.r[A_Y]*temp.r[A_Y])<rc*/)//�d�ɂ̐^��
		{
			temp.r[A_Z]=minZ-0.1*le;
			//temp.r[A_Z]=0.5*le;
			NODEw.push_back(temp);
			temp.id+=1;
			trans.push_back(-1);//���q�ɑΉ����Ȃ��ߓ_��-1���i�[
		}//*/
	}

	/*//�d�ɂ̐^��(�����I�ɕ��ׂ�)
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
			trans.push_back(-1);//���q�ɑΉ����Ȃ��ߓ_��-1���i�[
		}
	}//*/

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEw, "fluid.1.node");

	//.node�t�@�C���ǂݎ��
	in.load_node("fluid.1");

	//�܂��͗��̐ߓ_�݂̂ŕ���
	tetrahedralize("", &in, &out);
		// i �v�f���֐ߓ_��ǉ�(���̏������Ă��@�\���Ă���)
		// f .face�t�@�C���ɋ��E�ł͂Ȃ��ʂ��܂߂�
		// e .edge�t�@�C���̏o��(ON�ɂ���ƂȂ����~�܂��Ă��܂�)
		// n .neigh�t�@�C���̏o��

	//�o��
	out.save_nodes("fluid.2");
	out.save_elements("fluid.2");

	//////////////////�����܂łŗ��̐ߓ_�݂̂��g���āA���ׂĂ̗v�f���q�������ʂȃ��b�V�����ł���(fluid.2�Ŋm�F�\)*/


	///////////////�s�v�ȗv�f�̍폜

	//.node�̎擾
	GetPointList(NODEw, in, out);
	//.ele�̎擾
	GetTetrahedronList(ELEMw, in, out);

	//�_�~�[�v�f�̏���
	DelTetrahedron_InsideDummy(CON, NODEw, ELEMw, in, out);

	//�ߓ_-�v�f�֌W
	SetRelation_NodeElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W
	SetRelation_ElemElem(CON, NODEw, ELEMw);
	//�v�f-�v�f�֌W��藬�̕\�ʎ擾
	GetFacetList_from_neigh(CON, ELEMw, FACEw);


	//�_�~�[�ߓ_�̍폜
	//DelDummyNode(CON, NODEw, FACEw, num_dummy);


	//���q�ԍ���NODEw���ɑ��
	for(int i=0;i<NODEw.size();i++)	NODEw[i].part_no=trans[i];
	
	//�\�ʂ��\������ߓ_��I�����C�z��ԍ����l�߂� ���̓��������q�̐ߓ_���g���ꍇ�̓R�����g�A�E�g
	SelectFaceNode(CON, NODEw, FACEw);

	//�e�X�g
	//MakeNodeFile(CON, NODEw, "NODEw.node");
	//MakeFaceFile(CON, FACEw, "FACEw.face");


	/////////////////�v�f�m�F�p�t�@�C��///////////////////////////////////
	out.save_nodes("boundary_fluid");	//fluid.2.node�Ɠ����t�@�C��
	MakeElemFile(CON, ELEMw, "boundary_fluid.ele");
	MakeFaceFile(CON, FACEw, "boundary_fluid.face");
	////////////////�����܂łŐ��H�̃��b�V�����؂ꂽ//////////////////////*/
}


//�����v�f�̏���
void tetgen_function::DelThinTetrahedron(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"�s�v�ȗv�f�̍폜  ";

	double le=CON.get_distancebp();
	double delL;
	int flag=0;

	//del_length.dat������Γǂݍ���
	ifstream del("del_length.dat");
	if(!del)
	{
		cout<<"tetgen_config���폜����Ӓ���������  ";
		delL=TET.del_length;
	}
	else
	{
		cout<<"del_length.dat���폜����Ӓ���������  ";
		del>>delL;
		flag=1;		//�t�@�C������ǂݍ��񂾂�t���OON
	}
	del.close();

	if(flag==1)//�t�@�C���̐�����߂��Ă���
	{
		ofstream del2("del_length.dat");
		del2<<TET.del_length<<endl;
		del2.close();
	}

	cout<<"del_length="<<delL<<endl;


	cout<<"�����O�v�f��: "<<(int)ELEM.size()<<endl;

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag1=UNDEFINED;
		int flag2=UNDEFINED;
		int flag3=UNDEFINED;	//1�ō폜
		int del=OFF;
		int count=0;

		/*//4�_���\�ʐߓ_�ō\������Ă���΃t���O1ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag1=ON;
		else			flag1=OFF;//*/

		/*//4�_�������ߓ_�ō\������Ă���΃t���O2ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==FRFLUID)	count+=1;
		}
		if(count==4)	flag2=ON;
		else			flag2=OFF;//*/

		//��ł������ӂ�����΃t���OON
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

		//�폜
		if(del==ON)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//�C�e���[�^������
			it+=i;				//i���w��
			it=ELEM.erase(it);	//�폜���ăC�e���[�^��Ԃ�
			del_count++;
		}
		else i++;
	}

	//�v�f�ԍ��̐U�蒼��
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"������v�f��: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"�̗v�f���폜 ----------OK"<<endl;
}

//��яo�Đ�����v�f�̍폜
void tetgen_function::DelThinTetrahedron_SharpElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out)
{
	//�v�f-�v�f�֌W�����܂��Ă��邱�Ƃ�O��Ƃ���
	cout<<"������v�f�̗v�f�̍폜  ";

	double le=CON.get_distancebp();
	int flag=0;

	//boundary��FRFLUID(=0)�ɏ�����
	for(int i=0;i<(int)NODE.size();i++)	NODE[i].boundary=FRFLUID;

	//FACE�Ɋ܂܂��ߓ_�̂�BOFLUID�ɂ���
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)	NODE[FACE[i].node[n]].boundary=BOFLUID;
	}


	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag1=UNDEFINED;	//1�ō폜
		int del=OFF;
		int count1=0;
		int count2=0;

		//4�ʂ̂���2�ʈȏ�ڂ���v�f���Ȃ����͍̂폜
		//4�_�S�Ă�BOFLUID�ł�����͍̂폜
		for(int n=0;n<4;n++)
		{
			if(ELEM[i].nei_elem[n]==-1)	count1+=1;
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count2+=1;
		}
		if(count1>=2 || count2==4)	flag1=ON;
		else						flag1=OFF;//*/

		if(flag1==ON)	del=ON;

		//�폜
		if(del==ON)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//�C�e���[�^������
			it+=i;				//i���w��
			it=ELEM.erase(it);	//�폜���ăC�e���[�^��Ԃ�
			del_count++;
		}
		else i++;
	}

	//�v�f�ԍ��̐U�蒼��
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"������v�f��: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"�̗v�f���폜 ----------OK"<<endl;

}

//�s�v�ȗv�f�̏���(�O���_�~�[�ߓ_�@�p)
void tetgen_function::DelTetrahedron_OutsideDummy(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"�s�v�ȗv�f�̍폜"<<endl;
	cout<<"�����O�v�f��: "<<(int)ELEM.size()<<endl;

	double le=CON.get_distancebp();

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag=0;	//1�ō폜
		
		//1�ł��_�~�[(��C)�ߓ_������΃t���OON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==AIR)
			{
				flag=1;
				break;
			}
		}//*/

		/*//4�_���\�ʐߓ_�ō\������Ă���΃t���OON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag=1;//*/


		if(flag==1)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//�C�e���[�^������
			it+=i;				//i���w��
			it=ELEM.erase(it);	//�폜���ăC�e���[�^��Ԃ�
			del_count++;
		}
		else i++;
	}

	//�v�f�ԍ��̐U�蒼��
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"������v�f��: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"�̗v�f���폜 ----------OK"<<endl;
}


//�s�v�ȗv�f�̏���(�����_�~�[�ߓ_�@�p)
void tetgen_function::DelTetrahedron_InsideDummy(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"�s�v�ȗv�f�̍폜"<<endl;
	cout<<"�����O�v�f��: "<<(int)ELEM.size()<<endl;

	double le=CON.get_distancebp();

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag=0;	//1�ō폜
		
		//4�_���\�ʐߓ_�ō\������Ă���΃t���OON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag=1;//*/

		if(flag==1)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//�C�e���[�^������
			it+=i;				//i���w��
			it=ELEM.erase(it);	//�폜���ăC�e���[�^��Ԃ�
			del_count++;
		}
		else i++;
	}

	//�v�f�ԍ��̐U�蒼��
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"������v�f��: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"�̗v�f���폜 ----------OK"<<endl;
}


//�ߓ_-�v�f�֌W(tetgenio��edge���X�g����擾) �����O��tetgenio����f�[�^���擾���Ă�������
void tetgen_function::SetRelation_NodeNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"�ߓ_-�ߓ_�֌W";

	//�ꉞ������
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


//�ߓ_-�v�f�֌W
void tetgen_function::SetRelation_NodeElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"�ߓ_-�v�f�֌W";

	//�ꉞ������
	for(int i=0;i<(int)NODE.size();i++)
	{
		NODE[i].nei_node.clear();
		NODE[i].nei_elem.clear();
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			NODE[ELEM[i].node[n]].nei_elem.push_back(i);	//�v�f��ǉ�
		}
	}

	/*//�o�� ����эő吔�E�ŏ����̏o��
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

		//�ő�ŏ��X�V
		//if(max<(int)NODE[i].nei_elem.size())	max=(int)NODE[i].nei_elem.size();
		//if(min>(int)NODE[i].nei_elem.size())	min=(int)NODE[i].nei_elem.size();
	}
	fout.clear();

	//cout<<"�ő吔: "<<max<<endl;
	//cout<<"�ŏ���: "<<min<<endl;
	//*/

	cout<<"----------OK"<<endl;
}


//�v�f-�v�f�֌W
void tetgen_function::SetRelation_ElemElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"�v�f-�v�f�֌W";
	
	vector<int> nei_all;
	
	//������
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)	ELEM[i].nei_elem[n]=-2;	//����`��-2�Ƃ��Ă���
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		//4�̐ߓ_�̐ߓ_-�v�f�֌W�ɂ���v�f���i�[����
		nei_all.clear();
		
		for(int n=0;n<4;n++)
		{
			for(int j=0;j<(int)NODE[ELEM[i].node[n]].nei_elem.size();j++)
			{
				int elem=NODE[ELEM[i].node[n]].nei_elem[j];
				if(elem!=i)	nei_all.push_back(elem);	//�v�fi�͏���
			}
		}//nei_all�Ɋi�[����(�����v�f�ԍ����܂܂�Ă���\������)

		//�m�F�p
		//if(i==10000)	for(int j=0;j<(int)nei_all.size();j++)	cout<<nei_all[j]<<endl;

		//�ʂ�T��
		for(int ni=0;ni<4;ni++)
		{
			//�܂��͗v�fi�̖ʂ��w��
			int face[3];	//3�ԍ��Ŗʂ��w��
			int c=0;//�����グ�ϐ�

			for(int f=0;f<4;f++)
			{
				if(ni!=f)
				{
					face[c]=ELEM[i].node[f];
					c++;
				}
			}//face[3]��n�Ԗڂ̖ʂ��i�[

			//�ʂ�T��
			int correct_nei=-1;
			for(int j=0;j<(int)nei_all.size();j++)
			{
				int count=0;//���̃J�E���g��3�ɂȂ�Ίm��

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
			}//����������Ȃ�������correcr_nei�ɂ�-1�������Ă���

			ELEM[i].nei_elem[ni]=correct_nei;
		}
	}//*/

	/*//�o��
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


//�v�f������̗��̕\�ʒ�`
void tetgen_function::GetFacetList_from_neigh(mpsconfig &CON, vector<tetgen_element> &ELEM, vector<tetgen_facet> &FACE)
{
	cout<<"�v�f������̗��̕\�ʒ�`";

	//������
	FACE.clear();
	
	int id=0;	//id�p(�\�ʂ̐�)
	
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			if(ELEM[i].nei_elem[n]==-1)//�Ζʂ����݂��Ȃ����\��
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


//�\�ʐߓ_�ȊO���폜�C�\�ʐߓ_�ɔԍ���U��Ȃ���
void tetgen_function::SelectFaceNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE)
{
	//boundary���t���O�Ɏg�킹�Ă��炤�B
	//boundary==-1�͕\�ʂ��\�����Ă��Ȃ��ߓ_���폜
	//boundary==-2�͕\�ʂ��\�����Ă���ߓ_���V�����ߓ_�ԍ���^����
	
	//�Ƃ肠����������
	for(int i=0;i<(int)NODE.size();i++)	NODE[i].boundary=-1;
	
	//FACE�ɂ���ߓ_��flag��-2��
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			NODE[FACE[i].node[n]].boundary=-2;
		}
	}

	//�V�����ߓ_�ԍ��̌���
	int id=0;
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].boundary==-2)
		{
			NODE[i].boundary=id;
			id++;
		}
	}
	//�\�ʐߓ_�ɂ�boundary�ɐV�����ߓ_�ԍ�������B�����ߓ_�ɂ�-1������

	//FACE�̍\���ߓ_�ԍ��̕ϊ�
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			FACE[i].node[n]=NODE[FACE[i].node[n]].boundary;	//������-1��-2�ƂȂ���̂͂Ȃ��͂��B����Ώ�̏������Ԉ���Ă���
		}
	}

	//�����ߓ_�̍폜 NODE�̐ߓ_�ԍ��̕ϊ� boundary�����ɖ߂�
	int k=0;
	while(k<(int)NODE.size())
	{
		if(NODE[k].boundary==-1)
		{
			vector<tetgen_node>::iterator it=NODE.begin();	//�C�e���[�^������
			it+=k;				//k�Ԗڂ��w��
			it=NODE.erase(it);	//�폜���ăC�e���[�^��Ԃ�
		}
		else
		{
			NODE[k].id=NODE[k].boundary;
			NODE[k].boundary=WATER;
			k++;
		}
	}

}


//�_�~�[�ߓ_���폜
void tetgen_function::DelDummyNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int num_dummy)
{
	//NODE�̒��ɂ͑O���ɗ��̕\�ʐߓ_�C�㔼�Ƀ_�~�[�ߓ_���ł܂��Ċi�[����Ă���̂ŁC�㔼�̃_�~�[�ߓ_�̕����݂̂������΂悢
	//�_�~�[�v�f�������Ă���\�ʃf�[�^���擾���Ă���̂ŁC�\�ʂ��\������ߓ_�ԍ��̕ύX�͕s�v

	//�_�~�[�ϐ��̐�����popback�Ŗ����̗v�f�������
	for(int i=0;i<num_dummy;i++)
	{
		NODE.pop_back();
	}
}


//��C���E�ʍ쐬
void tetgen_function::SetAirBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	cout<<"��C���E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp_n;
	temp_n.id=0;
	temp_n.attribute=AIR;
	temp_n.boundary=AIR;

	//����������  1�������I�t�Z�b�g���Ă���
	double dL=TET.fine_air;
	int lx = int((CON.get_XL()-0.1*dL)/dL);
	int ux = int((CON.get_XR()+0.1*dL)/dL);
	int ly = int((CON.get_YD()-0.1*dL)/dL);
	int uy = int((CON.get_YU()+0.1*dL)/dL);
	int lz = int((CON.get_ZD()-0.1*dL)/dL);
	int uz = int((CON.get_ZU()+0.1*dL)/dL);

	//�ړ_�f�[�^�쐬
	//�Ód����
	if(CON.get_model_number()==14)
	{
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//��͗̈�̒[�ɗ����Ƃ��ߓ_��u��
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

	//���C���[����
	if(CON.get_model_number()==15)
	{
		temp_n.boundary=ELECTRODE1;

		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//��͗̈�̒[�ɗ����Ƃ��ߓ_��u��
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
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//��͗̈�̒[�ɗ����Ƃ��ߓ_��u��
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

	//���̕�
	if(CON.get_model_number()==25)
	{
		//temp_n.boundary=ELECTRODE1;
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					//if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//��͗̈�̒[�ɗ����Ƃ��ߓ_��u��
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

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEa, "NODEa1.node");

	in.load_node("NODEa1");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_air1");
	out.save_elements("boundary_air1");
	out.save_faces("boundary_air1");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEa, in, out, AIR);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEa1, "FACEa1.face");

}

//��C���E�ʍ쐬
void tetgen_function::SetAirBoundary2(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	cout<<"��C���E(�ו����p)�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp_n;
	temp_n.id=0;
	temp_n.attribute=AIR;
	temp_n.boundary=AIR;

	
	//���̕�
	if(CON.get_model_number()==25)
	{
		
		//����������
		double	le = CON.get_distancebp();
		double	rc = TET.radius_column;
		double	h =  TET.thickness_conduct;						//���̕Ќ���
		double	dh = TET.fine_conduct_z;							//���ݕ������b�V���e��
		int		nh = int((h/2+0.1*dh)/dh);	//���ݕ���������
		double	L =  TET.length_conduct;
		double	dL = TET.fine_conduct_xy;						//xy�������b�V���e��
		int		nL = int((L/2+0.1*dL)/dL);	//xy����������(�Б�)

		int lay=3;

		//���̑��̖�
		for(int z=-nh-lay;z<=nh+lay;z++)
		{
			for(int y=-nL-lay;y<=nL+lay;y++)
			{
				for(int x=-nL-lay;x<=nL+lay;x++)
				{
					if(x<-nL || x>nL || y<-nL || y>nL || z>nh || z<-nh)	//���d�ɗ̈�̒[���O���ɋ�C�ߓ_��u��
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

		
		//����������  1�������I�t�Z�b�g���Ă���
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
					if(x*dL2<(-nL-lay)*dL || x*dL2>(nL+lay)*dL || y*dL2<(-nL-lay)*dL || y*dL2>(nL+lay)*dL || z*dL2>(nh+lay)*dh || z*dL2<(-nh-lay)*dh)	//���d�ɗ̈�̒[���O���ɋ�C�ߓ_��u��
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


	cout<<"�ו����ߓ_�ʒu�̌v�㊮��"<<endl;
	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEa, "NODEa1.node");

	in.load_node("NODEa1");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_air2");
	out.save_elements("boundary_air2");
	out.save_faces("boundary_air2");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEa, in, out, AIR);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEa1, "FACEa1.face");

}


//���H�\�ʕt�ߒǉ��ߓ_
void tetgen_function::SetAirFineBoundary(mpsconfig &CON, vector<mpsparticle> &PART, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	//���H�̕\�ʂ̖@�������ɉ��w���̃��b�V���w���쐬����

	cout<<"���H�\�ʕt�ߒǉ��ߓ_"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=AIR;
	temp.boundary=0;

	//����������
	double le=CON.get_distancebp();
	double r=TET.radius_column+5*le;		//�~�����a
	double uz=TET.height_plate-20*le;		//�~���ő卂��
	double lz=-5*le;		//�~���ŏ����� -5*le
	double dL=le;					//�~�������������b�V���e��


	/*//�~������
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

	//���H�@���O������
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
		
	//���H�@����������
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


	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEa, "NODEa2.node");

	in.load_node("NODEa2");
	tetrahedralize("", &in, &out);
	//out.save_nodes("boundary_air2");
	//out.save_elements("boundary_air2");
	//out.save_faces("boundary_air2");
	//*/

	//���E�ʃf�[�^�擾
	//GetFacetList(FACEa, in, out, AIR);

	//��ʂƉ��ʂ��폜����
	int i=0;
	while(i<(int)FACEa.size())
	{
		double dis1=Distance(NODEa[FACEa[i].node[0]], NODEa[FACEa[i].node[1]]);
		double dis2=Distance(NODEa[FACEa[i].node[1]], NODEa[FACEa[i].node[2]]);
		double dis3=Distance(NODEa[FACEa[i].node[2]], NODEa[FACEa[i].node[0]]);

		/*int flag=0;//3�ɂȂ�����폜
		for(int n=0;n<3;n++)
		{
			if(NODEa[FACEa[i].node[n]].r[A_Z]>uz-2*le)	flag+=1;
			if(NODEa[FACEa[i].node[n]].r[A_Z]<lz+2*le)	flag+=1;
		}*/
		//if(flag>=1)
		if(dis1>2.0*le || dis2>2.0*le || dis3>2.0*le)
		{
			vector<tetgen_facet>::iterator it=FACEa.begin();	//�C�e���[�^������
			it+=i;				//i���w��
			it=FACEa.erase(it);	//�폜���ăC�e���[�^��Ԃ�
		}
		else i++;
	}
	//id�ĐU�蕪��
	for(int i=0;i<(int)FACEa.size();i++)
	{
		FACEa[i].id=i;
	}

	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEa, "FACEa2.face");
}

//�R�C�����E�ʍ쐬
void tetgen_function::SetCoilBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEd, vector<tetgen_facet> &FACEd)
{
	cout<<"�R�C�����E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=COIL;
	temp.boundary=COIL;

	//ringcoil.node����̓ǂݍ��� 
	//�I�v�V������������stl��ǂݍ��߂΁A���E�݂̂�.node��������(tetgen ����.stl)

	//���炩����stl����쐬����.node,.face��NODE,FACE�ɗ^����

	//node
	ifstream fp("coil_ring.node");
	if(!fp) cout<<"can not open coil.node"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//�ŏ���3�s���i�߂�
	
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

	//node�t�@�C���쐬
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
			tempf.node[n]-=1;//���̕��@�ō쐬����NODE�Ȃǂɂ��낦�邽�߁A�J�n�ԍ���0�ɂ���K�v������
		}
		tempf.boundary=COIL;
		FACEd.push_back(tempf);
	}
	fp2.close();
	////*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEd, in, out, COIL);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEd, "FACEd.face");
}

//���d�Ћ��E�쐬
void tetgen_function::SetConductBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEb, vector<tetgen_facet> &FACEb)
{
	cout<<"���̋��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=CONDUCT;//�Ƃ肠�������̕����ƂȂ�crucible�ɍގ����w��
	temp.boundary=CONDUCT;

	//����������
	double	le = CON.get_distancebp();
	double	rc = TET.radius_column;
	double	h =  TET.thickness_conduct;						//���̕Ќ���
	double	dh = TET.fine_conduct_z;							//���ݕ������b�V���e��
	int		nh = int((h/2+0.1*dh)/dh);	//���ݕ���������
	double	L =  TET.length_conduct;
	double	dL = TET.fine_conduct_xy;						//xy�������b�V���e��
	int		nL = int((L/2+0.1*dL)/dL);	//xy����������(�Б�)

	//���̑��̖�
	for(int z=-nh;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				//if(x==-nL || x==nL || y==-nL || y==nL || z==nh || z==-nh)	//���d�ɗ̈�̒[�ɗ����Ƃ��ߓ_��u��
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

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEb, "NODEb.node");

	in.load_node("NODEb");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_conduct");
	out.save_elements("boundary_conduct");
	out.save_faces("boundary_conduct");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEb, in, out, CONDUCT);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}

//��ڋ��E�ʍ쐬
void tetgen_function::SetCrucibleBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEd, vector<tetgen_facet> &FACEd)
{
	cout<<"��ڋ��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=CRUCIBLE;
	temp.boundary=CRUCIBLE;

	//ringcoil.node����̓ǂݍ��� 
	//�I�v�V������������stl��ǂݍ��߂΁A���E�݂̂�.node��������(tetgen ����.stl)

	//���炩����stl����쐬����.node,.face��NODE,FACE�ɗ^����

	//node
	ifstream fp("crucible.node");
	if(!fp) cout<<"can not open OutputMeshAndSolution-mass.txt"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//�ŏ���3�s���i�߂�
	
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
		fp>>temp.r[A_Z];//���Ƃ�stl��y,z������ւ���Ă���(�V���t�H�j�A��CAD)�̂ł�����ŏC��
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

	/*/node�t�@�C���쐬
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
			tempf.node[n]-=1;//���̕��@�ō쐬����NODE�Ȃǂɂ��낦�邽�߁A�J�n�ԍ���0�ɂ���K�v������
		}
		tempf.boundary=CRUCIBLE;
		FACEd.push_back(tempf);
	}
	fp2.close();
	////////*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEd, in, out, CRUCIBLE);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEd, "FACEd.face");
}


//���d�ɋ��E�ʍ쐬
void tetgen_function::SetPlateElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"���d�ɋ��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE2;
	temp.boundary=ELECTRODE2;

	//����������
	double	h = TET.height_plate;						//������
	double	dh = TET.fine_plate_t;						//���ݕ������b�V���e��
	int		nh = int((TET.thickness_plate+0.1*dh)/dh);	//���ݕ���������
	double	dL = TET.fine_plate_L;						//xy�������b�V���e��
	int		nL = int((TET.length_plate/2+0.1*dL)/dL);	//xy����������(�Б�)

	//�ړ_�f�[�^�쐬
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==0 || z==nh)	//���d�ɗ̈�̒[�ɗ����Ƃ��ߓ_��u��
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

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEp, "NODEp.node");

	in.load_node("NODEp");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_plate");
	out.save_elements("boundary_plate");
	out.save_faces("boundary_plate");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEp, in, out, ELECTRODE2);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//�~���d�ɋ��E�ʍ쐬
void tetgen_function::SetColumnElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"�~���d�ɋ��E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//����������
	double le=CON.get_distancebp();
	double rc=TET.radius_column;	//�~�����a
	double L=TET.length_column;		//�~������
	double dL=TET.fine_column_L;	//�~�������������b�V���e��
	//double z0=0;					//��ʂ̈ʒu

	
	//���
	//���S
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


	//����
	int flag=0;
	double dz=le;
	double z=0;
	rc*=1.005;	//�킸���ɑ�������(���b�V�����q����̂�h������)
	int nr=int(2.0*PI*rc/(2*le));//le��2�{����Ă��邱�Ƃɒ���
	double d_theta=360.0/(double)nr;
	
	while(1)
	{
		//z�����ւ̈ړ�
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

	//����
	//���S
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L+le;
	NODEc.push_back(temp);
	temp.id+=1;

	//for(double r=le;r<rc+0.1*le;r+=le)
	for(double r=rc;r>le-0.1*le;r-=2*le)//le��2�{����Ă��邱�Ƃɒ���
	{
		int nr=int(2.0*PI*r/(2*le));//le��2�{����Ă��邱�Ƃɒ���
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


	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEc, "NODEc.node");

	in.load_node("NODEc");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_column");
	out.save_elements("boundary_column");
	out.save_faces("boundary_column");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEc, in, out, ELECTRODE1);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEc, "FACEc.face");
}

//�Ό��d�ɋ��E�ʍ쐬
void tetgen_function::SetCounterElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"�Ό��d�ɋ��E�쐬(cad�ǂݍ���)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE2;
	temp.boundary=ELECTRODE2;

	
	in.load_stl("ea2");//���f���̍��W�n�ŏo�͂���ƂȂ����ǂݍ���ł���Ȃ�
	tetrahedralize("p", &in, &out);
	out.save_nodes("boundary_counter");
	out.save_elements("boundary_counter");
	out.save_faces("boundary_counter");
	//*/
	cout<<"stl����̃��b�V���쐬����"<<endl;
	/////
	//NODEp�ɁA�쐬���ꂽboundary���R�s�[����@//���������ʂȋC������B���������������ɂ�NODE�Ɋe���i���Ƃ̏����o�͂���boundarydata���ق��̊֐��ł܂Ƃ߂�d�g�݂���ς���K�v������
	//stl�œǂݍ��ނƐߓ_�ԍ���1����X�^�[�g���邽�߁A�]����0����n�܂��ł����W�ɍ��킹��ɂ�1���炷�K�v������B���̂Ƃ��A�ˏ������l�ɂ��炳�Ȃ��Ƃ����Ȃ�
	//in.initialize();
	//out.initialize();
	////

	//��ō�������b�V���́A���W�n����f���Ƃ͈قȂ�stl�ł��������̂Ȃ̂ŁA���f���̍��W�n�ŏo�͂������̂��Q�Ƃ��A�����悤�ɂ��炷
	ifstream fp2("ea2_dis.dat");
	if(!fp2) cout<<"ea2.dat"<<endl;
	fp2.unsetf(ifstream::dec);
	fp2.setf(ifstream::skipws);
	int stuff=0;
	string b;
	//for(int i=0;i<3;i++) getline(fp2, b);//�ŏ���3�s���i�߂�

	double x=0,y=0,z=0 ,x1=0,y1=0,z1=0, x2=0,y2=0,z2=0;
	
		fp2>>x1;
		fp2>>y1;
		fp2>>z1;
		fp2>>x2;
		fp2>>y2;
		fp2>>z2;
		
	cout<<x1<<" "<<y1<<" "<<z1<<endl;
	x=x1-x2; y=y1-y2; z=z1-z2;//�^�̒l�ƌ��݂̍��W�̍�
	fp2.close();
	

	//node
	ifstream fp("boundary_counter.node");
	if(!fp) cout<<"boundary_counter"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//�ŏ���3�s���i�߂�
	
	int node_num=0;
	stuff=0;
	fp>>node_num;
	for(int i=0;i<3;i++) fp>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<node_num;i++)//0����X�^�[�g
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

	//node�t�@�C���쐬
	//MakeNodeFile(CON, NODEp, "NODEp.node");

	//in.load_node("NODEp");
	//in.load_stl("ea_electro2");
	//tetrahedralize("", &in, &out);
	///out.save_nodes("boundary_counter2");
	//out.save_elements("boundary_counter2");
	//out.save_faces("boundary_counter2");

	//���E�ʃf�[�^�擾
	GetFacetListforCAD(FACEp, in, out, ELECTRODE2);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}

//���d�d�ɋ��E�ʍ쐬
void tetgen_function::SetDischargeElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"���d�d�ɋ��E�쐬(cad�ǂݍ���)"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	
	in.load_stl("ea1");//���f���̍��W�n�ŏo�͂���ƂȂ����ǂݍ���ł���Ȃ�
	tetrahedralize("p", &in, &out);
	out.save_nodes("boundary_discharge");
	out.save_elements("boundary_discharge");
	out.save_faces("boundary_discharge");
	//*/
	//cout<<"stl����̃��b�V���쐬����"<<endl;
	/////
	//NODEp�ɁA�쐬���ꂽboundary���R�s�[����@//���������ʂȋC������B���������������ɂ�NODE�Ɋe���i���Ƃ̏����o�͂���boundarydata���ق��̊֐��ł܂Ƃ߂�d�g�݂���ς���K�v������
	//stl�œǂݍ��ނƐߓ_�ԍ���1����X�^�[�g���邽�߁A�]����0����n�܂��ł����W�ɍ��킹��ɂ�1���炷�K�v������B���̂Ƃ��A�ˏ������l�ɂ��炳�Ȃ��Ƃ����Ȃ�
	//in.initialize();
	//out.initialize();
	////

	//��ō�������b�V���́A���W�n����f���Ƃ͈قȂ�stl�ł��������̂Ȃ̂ŁA���f���̍��W�n�ŏo�͂������̂��Q�Ƃ��A�����悤�ɂ��炷
	ifstream fp2("ea1_dis.dat");
	if(!fp2) cout<<"ea1.dat"<<endl;
	fp2.unsetf(ifstream::dec);
	fp2.setf(ifstream::skipws);
	int stuff=0;
	string b;
	//for(int i=0;i<3;i++) getline(fp2, b);//�ŏ���3�s���i�߂�

	double x=0,y=0,z=0 ,x1=0,y1=0,z1=0, x2=0,y2=0,z2=0;
	
		fp2>>x1;//�{���̈ʒu
		fp2>>y1;
		fp2>>z1;
		fp2>>x2;//���W�n�����킹�Ȃ������ʒu(stl���b�V���쐬�ʒu)
		fp2>>y2;
		fp2>>z2;
		
	cout<<x1<<" "<<y1<<" "<<z1<<endl;
	x=x1-x2; y=y1-y2; z=z1-z2;//�^�̒l�ƌ��݂̍��W�̍�
	fp2.close();
	

	//node
	ifstream fp("boundary_discharge.node");
	if(!fp) cout<<"boundary_discharge"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);
	//string b;
	//for(int i=0;i<3;i++) getline(fp, b);//�ŏ���3�s���i�߂�
	
	int node_num=0;
	stuff=0;
	fp>>node_num;
	for(int i=0;i<3;i++) fp>>stuff;
	if(stuff!=0) cout<<"error"<<endl;

	for(int i=0;i<node_num;i++)//0����X�^�[�g
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

	//node�t�@�C���쐬
	//MakeNodeFile(CON, NODEp, "NODEp.node");

	//in.load_node("NODEp");
	//in.load_stl("ea_electro2");
	//tetrahedralize("", &in, &out);
	///out.save_nodes("boundary_counter2");
	//out.save_elements("boundary_counter2");
	//out.save_faces("boundary_counter2");

	//���E�ʃf�[�^�擾
	GetFacetListforCAD(FACEp, in, out, ELECTRODE1);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//�y�䋫�E�ʍ쐬
void tetgen_function::SetBaseBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEb, vector<tetgen_facet> &FACEb)
{
	cout<<"�y�䋫�E�쐬"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//����������
	double	le = CON.get_distancebp();
	double	rc = TET.radius_column;
	double	h = TET.length_column;						//������
	double	dh = TET.fine_base;						//���ݕ������b�V���e��
	int		nh = int((TET.thickness_base+0.1*dh)/dh);	//���ݕ���������
	double	L = TET.length_base;
	double	dL = TET.fine_base;						//xy�������b�V���e��
	int		nL = int((TET.length_base/2+0.1*dL)/dL);	//xy����������(�Б�)


	//��ʐڑ�����
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L;
	NODEb.push_back(temp);
	temp.id+=1;

	rc*=1.005;//�~���ɍ��킹�Ă킸���ɑ�������
	for(double r=rc;r>le-0.1*le;r-=2*le)//le��2�{����Ă��邱�Ƃɒ���
	{
		int nr=int(2.0*PI*r/(2*le));//le��2�{����Ă��邱�Ƃɒ���
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

	//��ʕ\��
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

	//���̑��̖�
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==nh)	//���d�ɗ̈�̒[�ɗ����Ƃ��ߓ_��u��
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

	//node�t�@�C���쐬
	MakeNodeFile(CON, NODEb, "NODEb.node");

	in.load_node("NODEb");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_base");
	out.save_elements("boundary_base");
	out.save_faces("boundary_base");
	//*/

	//���E�ʃf�[�^�擾
	GetFacetList(FACEb, in, out, ELECTRODE1);
	//.face�t�@�C���쐬
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//���E�ߓ_�E���E�ʃf�[�^�̌���
void tetgen_function::UniteBoundaryData(mpsconfig &CON, 
					   vector<tetgen_node> &NODE, vector<tetgen_node> &NODEa1, vector<tetgen_node> &NODEa2, vector<tetgen_node> &NODEp, vector<tetgen_node> &NODEc, vector<tetgen_node> &NODEb, vector<tetgen_node> &NODEw, 
					   vector<tetgen_facet> &FACE, vector<tetgen_facet> &FACEa1, vector<tetgen_facet> &FACEa2, vector<tetgen_facet> &FACEp, vector<tetgen_facet> &FACEc, vector<tetgen_facet> &FACEb, vector<tetgen_facet> &FACEw, 
					   vector<int> &TRANS)
{
	cout<<"���E�ߓ_�E���E�ʃf�[�^�̌���"<<endl;

	NODE.clear();
	FACE.clear();
	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=0;	//�ߓ_�ԍ��̃I�t�Z�b�g��
	int offset_f=0;	//�\�ʔԍ��̃I�t�Z�b�g��


	//���H���E
	for(int i=0;i<(int)NODEw.size();i++)//�ߓ_
	{
		TRANS.push_back(NODEw[i].part_no);	//TRANS�ɗ��q�ԍ����i�[(boundary�ɗ��q�ԍ����i�[���Ă���)

		temp_n=NODEw[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=WATER;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEw.size();i++)//�\��
	{
		temp_f=FACEw[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=WATER;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEw.size();
	offset_f+=(int)FACEw.size();

	//��C���E
	for(int i=0;i<(int)NODEa1.size();i++)//�ߓ_
	{
		temp_n=NODEa1[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa1.size();i++)//�\��
	{
		temp_f=FACEa1[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEa1.size();
	offset_f+=(int)FACEa1.size();

	//��C���E ���𑜓x��
	for(int i=0;i<(int)NODEa2.size();i++)//�ߓ_
	{
		temp_n=NODEa2[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa2.size();i++)//�\��
	{
		temp_f=FACEa2[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEa2.size();
	offset_f+=(int)FACEa2.size();

	//���d�ɋ��E
	for(int i=0;i<(int)NODEp.size();i++)//�ߓ_
	{
		temp_n=NODEp[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE2;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEp.size();i++)//�\��
	{
		temp_f=FACEp[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE2;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEp.size();
	offset_f+=(int)FACEp.size();

	//�~���d�ɋ��E
	for(int i=0;i<(int)NODEc.size();i++)//�ߓ_
	{
		temp_n=NODEc[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEc.size();i++)//�\��
	{
		temp_f=FACEc[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	offset_n+=(int)NODEc.size();
	offset_f+=(int)FACEc.size();

	//�y�䋫�E
	for(int i=0;i<(int)NODEb.size();i++)//�ߓ_
	{
		temp_n=NODEb[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEb.size();i++)//�\��
	{
		temp_f=FACEb[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//�I�t�Z�b�g�ʍX�V
	//offset_n+=(int)NODEb.size();
	//offset_f+=(int)FACEb.size();

	//cout<<"----------OK"<<endl;
}

//���E�ߓ_�E���E�ʃf�[�^�̒ǉ�
void tetgen_function::AddBoundaryData(mpsconfig &CON, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int attribute)
{
	//NODE,FACE�Ɋi�[����Ă���e���i�̃f�[�^���CNODEall,FACEall�Ɋi�[���Ă����D

	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=(int)NODEall.size();	//�ߓ_�ԍ��̃I�t�Z�b�g��
	int offset_f=(int)FACEall.size();	//�\�ʔԍ��̃I�t�Z�b�g��


	//�ߓ_�̒ǉ�
	for(int i=0;i<(int)NODE.size();i++)//�ߓ_
	{
		temp_n=NODE[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=attribute;
		NODEall.push_back(temp_n);
	}

	//�ʂ̒ǉ�
	for(int i=0;i<(int)FACE.size();i++)//�\��
	{
		temp_f=FACE[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=attribute;
		FACEall.push_back(temp_f);
	}
}

//TRANS[]�̊i�[
void tetgen_function::SetTRANS(vector<tetgen_node> &NODE, vector<int> &TRANS)
{
	//TRANS[i]�ɂ́A�ߓ_�ԍ�i�ɑΉ����闱�q�ԍ����i�[����B
	//FEM3D.cpp �ł́A�ߓ_�ԍ���1����n�܂�̂ŁATRANS[0]�ɂ͐錾���-1������B�i�����ł͊��ɓ����Ă���j

	for(int i=0;i<(int)NODE.size();i++)
	{
		TRANS.push_back(NODE[i].part_no);	//TRANS�ɗ��q�ԍ����i�[
	}
}


//�ގ��̌���i�܂��������j
void tetgen_function::DecisionAttribute(vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double M[4];	//�ߓ_�̍ގ����i�[
	
	for(int i=0;i<(int)ELEM.size();i++)
	{	
		if(ELEM[i].attribute==AIR)
		{
			for(int n=0;n<4;n++)	M[n]=NODE[ELEM[i].node[n]].attribute;

			if(M[0]==AIR || M[1]==AIR || M[2]==AIR || M[3]==AIR)	//1�ł���C�ߓ_������΋�C�v�f
			{
				ELEM[i].attribute=AIR;
				out.tetrahedronattributelist[i]=AIR;
			}
			else
			{
				ELEM[i].attribute=WATER;	//�c��͐�
				out.tetrahedronattributelist[i]=WATER;
			}
		}
	}
}


//�ގ��̏C��  ���̊֐����g���Ƃ��͒��O��tetgenio����f�[�^���擾���Ă�������
void tetgen_function::ModifyAttribute(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	//tetgenio����f�[�^���擾�ς݂Ɖ��肷��
	//�ߓ_-�ߓ_�֌W�������Ă���Ɖ��肷��

	/*double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;*/

	//�v�f�̍ގ��̏C��  ���̎��_�ł͐��v�f�ƂȂ�̈�͍ގ��ԍ���0(����`)�ƂȂ��Ă���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		//����`(0)�̍ގ��𐅂ɂ���
		if(ELEM[i].attribute==0)
		{
			ELEM[i].attribute=WATER;
		}
	}

	//�ߓ_�̍ގ���v�f�̍ގ��ƍ��킹��

	//�܂��S����C�ɂ���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		NODE[i].attribute=AIR;
	}
	//���ߓ_�̌���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}
	//�d�ɐߓ_�̌���
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}

	//�ߓ_-�ߓ_�֌W���ߓ_�̏C��
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].attribute==AIR)
		{
			int num_air=0;
			int num_ele=0;
		}
	}
}


//�ގ��̏C��  tetgenio�𒼐ڕҏW
void tetgen_function::ModifyAttribute_tetgenio(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;

	//�ߓ_�����}�����ꂽ�ߓ_�̍ގ����Ƃ肠������C�Ƃ���
	for(int i=0;i<out.numberofpoints;i++)
	{
		if(out.pointattributelist[i]==0)
		{
			out.pointattributelist[i]=AIR;
		}
	}//*/

	//�v�f�̍ގ��̏C��  ���̎��_�ł͐��v�f�ƂȂ�̈�͍ގ��ԍ����f�t�H���g��0(����`)�ƂȂ��Ă���
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		//����`�̍ގ��𐅂ɂ���
		if(out.tetrahedronattributelist[i]==0)
		{
			out.tetrahedronattributelist[i]=WATER;//0���������̂�WATER�ɂȂ�
		}
	}
	
	//�Ód�����̏ꍇ�@����-�~���ԂƉ~��-�y��Ԃɋ�C�̌��Ԃ��ł��Ă���̂ł���𖄂߂�
	if(CON.get_model_number()==14)
	{
		//for(int times=0;times<30;times++)
		{
			for(int i=0;i<out.numberoftetrahedra;i++)
			{
				//��C�v�f�𔻒肵�ĕύX����
				if(out.tetrahedronattributelist[i]==AIR)
				{
					int flag=0;

					//�e�ߓ_�̍ގ��ɂ�锻��
					int Mnode[4];	//�ߓ_�̍ގ����i�[
					int node_A=0;	//1�v�f�ɂ������C�ߓ_��
					int node_W=0;	//1�v�f�ɂ����鐅�ߓ_��
					int node_E1=0;	//1�v�f�ɂ�����~���d�ɐߓ_��
					int node_E2=0;	//1�v�f�ɂ����镽�d�ɐߓ_��
					for(int n=0;n<4;n++)	Mnode[n]=(int)out.pointattributelist[out.tetrahedronlist[4*i+n]];
					
					//�ׂ荇���v�f�ɂ�锻��
					int Melem[4];	//�ׂ荇���v�f�̍ގ����i�[
					int elem_A=0;	//�ׂ荇����C�v�f��
					int elem_W=0;	//�ׂ荇�����v�f��
					int elem_E1=0;	//�ׂ荇���d�ɗv�f��
					for(int n=0;n<4;n++)	Melem[n]=(int)out.tetrahedronattributelist[out.neighborlist[4*i+n]];
					
					for(int n=0;n<4;n++)	//1�l�̗v�f�̒��Ő��Ɠd�ɂ̐ߓ_���𐔂���
					{
						if(Mnode[n]==AIR)			node_A+=1;
						if(Mnode[n]==WATER)			node_W+=1;
						if(Mnode[n]==ELECTRODE1)	node_E1+=1;
						if(Mnode[n]==ELECTRODE2)	node_E2+=1;
						if(Melem[n]==AIR)			elem_A+=1;
						if(Melem[n]==WATER)			elem_W+=1;
						if(Melem[n]==ELECTRODE1)	elem_E1+=1;
					}

					if(node_W==4)	flag=1;							//�S�Ă����ߓ_�̂Ƃ��͐�
					if(node_A==0 && node_W>0 && node_E1>0)	flag=1;	//��C�ߓ_��0�Ő��ߓ_�Ɠd�ɐߓ_��1�ȏ�Ƃ��͐�
					if(elem_A<=1)	flag=1;							//�אڂ����C�v�f��1�ȉ��̂Ƃ��͐�
					if(elem_W>0 && elem_E1>0)	flag=1;				//�אڂ��鐅�v�f�Ɠd�ɗv�f��1�ȏ�̂Ƃ��͐�

					//if(node_E2==4)	flag=2;	//�S�Ă��d�ɐߓ_�̂Ƃ��͓d��  ���������Ɖ~���Ɠy��̋��ɂ킸���Ȏ΂߂̖ʂ��ł��Ă��܂�

					//if((M[0]==AIR || M[0]==WATER) && M[1]!=AIR && M[2]!=AIR && M[3]!=AIR)	flag=1;//-------------�@
					//if(M[0]==M[1] && M[1]==M[2] && M[2]==M[3] && M[0]==WATER)		flag=1;
					//*/

					//�d�S���W�ɂ�锻��(�ߓ_�ގ��ɂ�锻��ł��ڂꂽ���̂��E��)
					double r[3]={0,0,0};
					for(int d=0;d<3;d++)
					{
						for(int n=0;n<4;n++)	r[d]+=out.pointlist[3*(out.tetrahedronlist[4*i+n])+d];
						r[d]/=4;
					}
					//if(sqrt(r[A_X]*r[A_X]+r[A_Y]*r[A_Y])<rc-le && r[A_Z]<le && r[A_Z]>0)	flag=1;	//�d�ɂ̉��ɐj�݂����ȗ��̗v�f���ł��邱�Ƃ�����̂�le�����Z�����Ă���  ����R�ꎟ�̍ގ�����ŏE���Ă��炤���ƂɊ���
					if(sqrt(r[A_X]*r[A_X]+r[A_Y]*r[A_Y])<rc && r[A_Z]>-L && r[A_Z]<-L+le)	flag=2;	//�y��Ɖ~���̋��ڂ͍��W�ɂ�锻��ł���r�I�ȒP
					//*/

					if(flag==1)	//�t���O��1�Ȃ琅
					{
						out.tetrahedronattributelist[i]=WATER;
					}
					if(flag==2)	//�t���O��2�Ȃ�~���d��
					{
						out.tetrahedronattributelist[i]=ELECTRODE1;
					}
				}
			}
		}
	}

	//�����ǉ����ꂽ�ߓ_�̂�attribute��boundary_marker�̏C��
	//���̎��_�ł͗��̓�����d�ɓ�����d�ɓ����̐ߓ_�͋�C�ߓ_�ɂȂ��Ă���̂ł��ꂼ��̍ގ��ɏC������
	//���d�ɂ̈�ԏ�̐ߓ_�͓d�ɂ̐ߓ_�Ƃ��邽�߁A��ɐ��̐ߓ_���߂Ă���d�ɐߓ_�����߂�
	
	//���v�f���\������ߓ_�͐��ߓ_��
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

	
	
	//�Ód�����@���d�ɂ���щ~���d��(boundary_marker�ɂ��Ă͕��Ɖ~���ŕ�����)
	if(CON.get_model_number()==14)
	{
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			//�~��
			if(out.tetrahedronattributelist[i]==ELECTRODE1)
			{
				for(int n=0;n<4;n++)
				{
					out.pointattributelist[out.tetrahedronlist[i*4+n]]=ELECTRODE;
					out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=ELECTRODE1;
				}
				//�ގ���ELECTRODE�ɖ߂��Ă���
				out.tetrahedronattributelist[i]=ELECTRODE;
			}
			//����
			if(out.tetrahedronattributelist[i]==ELECTRODE2)
			{
				for(int n=0;n<4;n++)
				{
					out.pointattributelist[out.tetrahedronlist[i*4+n]]=ELECTRODE;
					out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=ELECTRODE2;
				}
				//�ގ���ELECTRODE�ɖ߂��Ă���
				out.tetrahedronattributelist[i]=ELECTRODE;
			}
		}
	}

}


//�v�f�̍ו���
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

//�ߓ_�ԋ����v�Z�֐�
double tetgen_function::Distance(tetgen_node &point1, tetgen_node &point2)
{
	double dis=0;

	for(int d=0;d<3;d++)	dis+=(point2.r[d]-point1.r[d])*(point2.r[d]-point1.r[d]);

	return sqrt(dis);
}


//�v�f�d�S���W�v�Z�֐�
void tetgen_function::CalcBarycentricElement(mpsconfig&, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	double r[3]={0,0,0};

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int d=0;d<3;d++)	for(int n=0;n<4;n++)	r[d]+=NODE[ELEM[i].node[n]].r[d];
		for(int d=0;d<3;d++)	ELEM[i].g[d]=r[d]/4;
	}
}