#include "stdafx.h"	

//���������[�N���o
//#define _CRTDBG_MAP_ALLOC

#include"header.h"		//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
//#include"define.h"		//#define �i�[
//#include"PART.h"		//class PART��`
//#include"CONFIG.h"		//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
#include"BEMclass.h"	//FEM3D�֌W��class ��`
//#include"FEM3Dclass.h"	//FEM3D�֌W��class ��`
//#include<omp.h>
//#include<vector>

#include"function.h"

//////
//���������[�N���o
//#include <crtdbg.h>
//#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)
////////*/

int _tmain(int argc, _TCHAR* argv[])
{
	//���������[�N���o
	 _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	 _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	 _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	 ///////////////

	mpsconfig CON;//��͏����i�[�׽
	
    int particle_number;//�S���q��
    int fluid_number;	//���̗��q��
	int out;			//fluid_number<=i<out��INWALL�Q�ŁAout<=i��OUTWALL�Q�ɂȂ�B �Ȃ��ĂȂ��Ȃ玩���ł��Ƃŕ��ёւ���
	//int count;			//�J�E���g�p�ϐ�
	int count_avs=0;
	int order_sw=OFF;	//���q����ёւ���K�v�����邩�ǂ����̃t���O
    double dt=CON.get_dt();
    int STEP=CON.get_step();
    double g[DIMENTION]={0,0,0};
    if(CON.get_dimention()==2)  g[A_Y]=CON.get_g();
    if(CON.get_dimention()==3)  g[A_Z]=CON.get_g();
    /////MPS�ɂ�����萔�v�Z
	double n0=initial_pnd(CON.get_re(),CON.get_dimention(),CON.get_model_set_way()); //�������q���xn0
    double N0=initial_pnd(CON.get_re2(),CON.get_dimention(),CON.get_model_set_way());//���׼�ݗp�������q�����x
    double n0_4=initial_pnd(CON.get_re4(),CON.get_dimention(),CON.get_model_set_way());//freeon�p�������q���x
    double lamda=calclamda(&CON); //���׼�ݗp��
	CON.set_Cst(calc_Cst(&CON));
    double TIME=0;		//��͎���
	double Umax=0;		//�ő呬�x
	double mindis;		//���ݗp�̍Œᗱ�q�ԋ���
	int t_old=0;		//EM_distance�֌W�̏o�͗p�B���O��EM�����Ƃ��̃X�e�b�v�����L��
    printf("n0= %10.8f\n",n0);
    printf("N0= %10.8f\n",N0);
    printf("lamda= %10.8f\n",lamda);
	cout<<"Cst= "<<CON.get_Cst()<<", Cst*"<<CON.get_C_times()<<"="<<CON.get_Cst()*CON.get_C_times()<<endl;
    //////////////////////////////////*/

	//mindis=CON.get_distancebp();

	/*
	////////
	cout<<"particlemovie�̕ϊ�"<<endl;
	ifstream fin("particle_movie_0522_termo.mgf");
	if(!fin) cout<<"cannot open particle_movie_0522_termo"<<endl;
	fin.unsetf(ifstream::dec);
	fin.setf(ifstream::skipws);

	string b;
	for(int i=0;i<2;i++) getline(fin, b);//�ŏ���6�s���i�߂�

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

		//�ǂݍ���
		if(i==0) for(int j=0;j<4;j++) getline(fin, b);//�ŏ���4�s���i�߂�
		else for(int j=0;j<5;j++) getline(fin, b);//�ŏ���4�s���i�߂�
		fin>>num;

		//�o��
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

			//���x�R���^�[�̕ϊ�
			T=sqrt(red)*100+293;
			if(T>Tmax) T=Tmax;
			level=(T-Tmin)/(Tmax-Tmin);
			red_m=level*level;
			green_m=-4*(level*(level-1));
			blue_m=1-red_m;
 
			fout<<x<<" "<<y<<" "<<z<<" ";//���W�o��
				
			fout<<d<<" ";//���q�̑傫���o��
	
			fout<<red_m<<" "<<green_m<<" "<<blue_m<<endl;//�F�o��
		}
		
	}
	fin.close();
	fout.close();
	*/
	///////
	/////

	/////�������q�z�u�������݁@restart=ON�̏ꍇ�͗��q���ǂݍ���
    if(CON.get_restart()==OFF)
	{   
		if(CON.get_model_inherit()==0)
		{
			if(CON.get_model_set_way()==0)		set_initial_placement(&CON,&particle_number); //�����i�q�ɂ�郂�f���Z�b�g�B�������\�ʌ`��͊K�i��
			else if(CON.get_model_set_way()==1)	set_initial_placement_using_MD(&CON,&particle_number);//���q���͊w�ɂ�郂�f���Z�b�g�B�v�Z���Ԃ͂����邪�A�\�ʂ����ꂢ�ɕ\��
		}
		else cout<<"�O��̏������q���f���̈��p��"<<endl;
		
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
		cout<<"�O��ɂ������͂������p��"<<endl;
    }
    //////////////////////
	if(CON.get_model_inherit()==1)
	{
		ifstream fa("particle_number.dat");
		if(!fa) cout<<"cannot open number.dat"<<endl;
		fa>>particle_number;
		fa.close();
	}
    cout<<"���q����"<<particle_number<<endl;

	/////�ǋ��E�|���S�����̃e�X�g
	int wall_poly_number=0;//�ǋ��E�̃|���S����
	wall_poly_number=3;//�ǂݍ��݂Ȃǂł��Ƃ߂邱��
	double wall_a[3];//2�����̒���������ax+by+c=0�@3DFEM���Q�l��node��edge�̂悤�ȏ����������邱��
	double wall_b[3];
	double wall_c[3];
	if(CON.get_wall_poly()==1)
	{
		//�ǋ��E�̐ݒ� //���b�V���Ȃǂ���ǂݍ��ނ��Ƃ��l����ƁA���̈ʒu�͂܂�������
		if(CON.get_model_number()==1)
		{
			if(CON.get_dimention()==2)
			{
				wall_a[0]=1; wall_b[0]=0; wall_c[0]=-0.011;//��������
				wall_a[1]=1; wall_b[1]=0; wall_c[1]=0.011;//�E������
				wall_a[2]=0; wall_b[2]=1; wall_c[2]=-0.021;//��������
			}
		}
	}
	////�Ǐd�݊֐��̌v�Z



	

	/////INDEX�֌W
    //int *INDEX=new int[CON.get_number_of_mesh()];	//�e�i�q�Ɋ܂܂�闱�q�����i�[
    cout<<"X_mesh="<<CON.get_X_mesh()<<" Y_mesh="<<CON.get_Y_mesh()<<" Z_mesh="<<CON.get_Z_mesh()<<endl;
    cout<<"number_of_mesh="<<CON.get_number_of_mesh()<<endl;
    ///////////////////*/

	clock_t t1=clock();		//�o�ߎ��Ԃ�b�ŕ\������ɂ́ACLOCKS_PER_SEC�Ŋ���K�v������

	/////�e���PART[i](particle[i])���쐬
	mpsparticle PART0;
	vector<mpsparticle> PART;
	for(int i=0;i<particle_number;i++) PART.push_back(PART0);
    //////////////////////////////

	////BEM�p��class�쐬
	vector<point2D> NODE;
	vector<element2D> ELEM;
	vector<BEMpoint3D> BEMNODE3D;		//�Œ�ߓ_���class
	vector<BEMelement3D> BEMELEM3D;	//�Œ�v�f���class
	int s_node_num;					//�Œ�ߓ_��,���I�ߓ_��
	int s_elem_num;	

	//////////////////////////

	////FEM�p��class�쐬
	vector<point3D> NODE_FEM3D;
	vector<element3D> ELEM_FEM3D;
	vector<edge3D> EDGE_FEM3D;
	vector<point3D> NODE_jw;
	vector<element3D> ELEM_jw;
	int node_FEM3D=0;
	int nelm_FEM3D=0;	//�S�ߓ_��,�v�f��
	int nedge_FEM3D=0;

	//////////////////////////


	///���q�f�[�^���t�@�C������ǂݎ��
	input_particle_data(&CON,PART,particle_number,1);//�Ō�̈�����1��n������initial�f�[�^��ǂݎ��A����ȊO�Ȃ�O�ï�߃f�[�^��ǂݎ��

	//�e���q�����J�E���g or ���ёւ�
	calc_numbers_of_particles_and_change_the_order(&CON,PART,particle_number,&fluid_number,&out,&order_sw);

	//double *Un[DIMENTION];
	//for(int D=0;D<DIMENTION;D++) Un[D]=new double [fluid_number];			//n�ï�ߎ��̑��x(�z��͒��O)���L��. 
	double *previous_Un[DIMENTION];
	for(int D=0;D<DIMENTION;D++) previous_Un[D]=new double [fluid_number];	//n-1�ï�ߎ��̑��x(�z��͒��O)���L��. ���̑��x�����߂�̂Ɏg�p����
	double *PND2=new double [particle_number];		//���͌v�Z�p�ɁA�e�ï�ߊJ�n���̗��q�����x���L��

	if(CON.get_set_zero_speed()==ON && CON.get_restart()==ON)
	{
		for(int i=0;i<fluid_number;i++) for(int D=0;D<CON.get_dimention();D++) PART[i].u[D]=0;//restart�ł����x�͏������B�J���p
	}

	double distance=0.0;//EM_distance>0�̂Ƃ��ɗp���锻��p����
	for(int i=0;i<particle_number;i++)
	{
		PART[i].T=CON.get_initialT();//������
		PART[i].heat_gene_before1=0;
		PART[i].heat_generation=0;
		PART[i].L=CON.get_distancebp();
		PART[i].dir_Pst=0.0;
		PART[i].dir_Pem=0.0;

	}

	//�ő呬�x�����߂Ă����B���ʂ̓[�������Arestart�����Ƃ��Ƃ��ɂ����ŋ��߂Ă����Ȃ��ƁA1step�߂�curan���Z�o�����������Ȃ�B
	for(int i=0;i<particle_number;i++)
	{ 
		double speed=0;//���q���x
		for(int D=0;D<DIMENTION;D++) speed+=PART[i].u[D]*PART[i].u[D];
		if(speed>Umax)  Umax=speed;		
	}
	Umax=sqrt(Umax);
	
	unsigned int time0=GetTickCount();	//��͎n�߂̎��Ԃ��L��
	//�v�Z�X�^�[�g
	for(int t=1;t<=STEP;t++)
	{	
		unsigned int timet=GetTickCount();
			
		cout<<"�z��� start:step="<<t<<" ���̗��q��="<<fluid_number<<" �S���q��="<<particle_number<<endl;

		//�e���q�����J�E���g or ���ёւ�
		if(order_sw==ON) calc_numbers_of_particles_and_change_the_order(&CON,PART,particle_number,&fluid_number,&out,&order_sw);//������order_sw��OFF�ɖ߂��Ă���
		
		////////////////////�z���/////////////////////////////

		//�ߗח��q�֌W�v�Z
		calc_neighbor_relation(&CON,PART,particle_number,n0_4,fluid_number,out);
		
		/*/////�z��͂̑O��reloadINDEX	
		reload_INDEX(&CON,PART,particle_number,INDEX);//�i�q���̗��q���X�V
		int **MESH = new int *[CON.get_number_of_mesh()];
		count=0;
		for(int i=0;i<CON.get_number_of_mesh();i++)
		{       
			count+=INDEX[i];
			MESH[i]=new int [INDEX[i]];
		}
		if(count>particle_number) cout<<"INDEX error ���q������?"<<endl;
		if(count<particle_number) cout<<"INDEX error ���q������?"<<endl;
		reload_INDEX2(&CON,PART,particle_number,MESH);
		////////////////////////*/

		////////�Ǐd�݊֐��̌v�Z
		//if(CON.get_wall_poly()==1 && CON.get_model_number()==0) calc_wallZ(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,fluid_number,out);

		
		unsigned int timeA=GetTickCount();
		/*////
		if(CON.get_freeon()==1) freeon(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,fluid_number,out);//�\�ʔ���
		else if(CON.get_freeon()==2) freeon2(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,fluid_number,out,t);//�\�ʔ���
		//else if(CON.get_freeon()==4) freeon4(&CON,PART,particle_number,n0_4,INDEX,MESH,&mindis,n0,fluid_number,e0,out,t);
		else cout<<"�\�ʔ��薢����"<<endl;
		if(CON.get_surface_judge2()==ON)
		{
			surface_judge2(&CON,PART,fluid_number,particle_number);
		}
		cout<<"�z��͑O�̗��q�ˑ��֌W�v�Z�I���@�|�|time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		////*/

		////�Œᗱ�q�ԋ��������Ƃ߂�
		//double min0=CON.get_distancebp();//�Œᗱ�q�ԋ���
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
			/////�Œᗱ�q�ԋ��������Ƃ܂���
			//cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
			//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
		}
		/////�Œᗱ�q�ԋ��������Ƃ܂���
		cout<<"mindis="<<mindis<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;

		//���ݐ��ɂ��dt����
		culan(&CON,PART,fluid_number,t,&dt,mindis,Umax,g);

		//�������q�����x���z���o�͂���
		if(t==1) output_particle_density(&CON, PART, fluid_number, n0, particle_number, t);

		//if(CON.get_modify_position()==ON) modify_position(&CON,PART, fluid_number,dt);

		///�S�����z��
		double *laplacian[DIMENTION];
		for(int D=0;D<DIMENTION;D++) laplacian[D]=new double [fluid_number];
		///�\�ʒ��͊֌W�z��
		double *potential[DIMENTION];
		for(int D=0;D<DIMENTION;D++) potential[D]= new double [fluid_number];
		//���q�̊O��
		double *F[DIMENTION];					
		for(int D=0;D<DIMENTION;D++) F[D]=new double [fluid_number];
		for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) F[D][i]=0;//������

		double *Un[DIMENTION];
		for(int D=0;D<DIMENTION;D++) Un[D]=new double [fluid_number];			//n�ï�ߎ��̑��x(�z��͒��O)���L��.
		for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) Un[D][i]=PART[i].u[D];//�z��͑O�̑��x���L��

		//////�S�����v�Z
		calc_viscous_term(&CON,PART,fluid_number,particle_number,dt,N0,laplacian,lamda,t);
		
		//�\�ʒ��͌v�Z
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

					cout<<"FEM�I���@�|�|time="<<(GetTickCount()-timeF)*0.001<<"[sec]"<<endl;

					for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) PART[i].F[D]=F[D][i];//���߂��d���͂�PART�Ɋi�[

					//plot_F(&CON,PART,fluid_number,F,t);//�d���͏o��
				
					//FEM�ɓ���Ȃ�step�ł����O��FEM�ŋ��߂��d���͂�p���邽�߁A�t�@�C���o��
					//�|�X�g�����ŏ��������q�ɑΉ��ł��Ȃ��̂ŋp��
					//ofstream g("F_FEM.dat");
					//for(int i=0;i<fluid_number;i++) g<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
					//g.close();
				}
				else//FEM�����Ȃ��Ƃ��̓t�@�C������ǂݍ���//�|�X�g�����ŏ��������q�ɑΉ��ł��Ȃ��̂ŋp���BPART.F�Ɋi�[�����l���g��
				{
					if(CON.get_T_field()==ON) for(int i=0;i<particle_number;i++) PART[i].heat_generation=PART[i].heat_gene_before1;//�Q�d�������v�Z���Ȃ��Ƃ��ł��l���o���Ă���
				}
			}
			else//���q�̈ړ������ɂ����FEM���s�����ǂ������߂�
			{
				if(t==1 || distance>CON.get_EM_distance()*CON.get_distancebp())
				{	
					distance=0.0;//���Z�b�g
					if(t!=1) cout<<"distance>"<<CON.get_EM_distance()<<"le�̂���FEM�J�n"<<endl;
					unsigned int timeF=GetTickCount();
					if(CON.get_EM_method()==1) FEM3D_calculation(&CON,&node_FEM3D,&nelm_FEM3D,&nedge_FEM3D,NODE_FEM3D,ELEM_FEM3D,EDGE_FEM3D,fluid_number,F, t, TIME,PART, fluid_number, particle_number,dt,NODE_jw,ELEM_jw,CON);
					else if(CON.get_EM_method()==2)
					{
						if(CON.get_dimention()==2) BEM2D(&CON,NODE,ELEM,&s_node_num,&s_elem_num,PART,fluid_number,particle_number,F,t);
						if(CON.get_dimention()==3) BEM3D(&CON,BEMNODE3D,BEMELEM3D,&s_node_num,&s_elem_num,PART,fluid_number,particle_number,F,t);
					}
					else if(CON.get_EM_method()==3) Magnetic_Moment_Method(&CON,PART,F, n0, lamda, fluid_number, particle_number);

					cout<<"FEM�I���@�|�|time="<<(GetTickCount()-timeF)*0.001<<"[sec]"<<endl;

					for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) PART[i].F[D]=F[D][i];//���߂��d���͂�PART�Ɋi�[

					/////////////////////////////CPU time �Ɖ�͕������Ԃ̊֌W�����v���b�g
					if(t==1)
					{
						ofstream t1("EMstep.dat");		//�����X�e�b�v���A�c��dt�̃O���t
						t1.close();
					}
				
					ofstream t1("EMstep.dat",ios :: app);
					t1<<t<<" "<<t-t_old<<endl;
					t1.close();
					t_old=t;

		////////////////////

					//plot_F(&CON,PART,fluid_number,F,t);//�d���͏o��
				
					//FEM�ɓ���Ȃ�step�ł����O��FEM�ŋ��߂��d���͂�p���邽�߁A�t�@�C���o��
					//ofstream g("F_FEM.dat");
					//for(int i=0;i<fluid_number;i++) g<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
					//g.close();
				}
				else//FEM�����Ȃ��Ƃ��̓t�@�C������ǂݍ���
				{
					if(CON.get_T_field()==ON) for(int i=0;i<particle_number;i++) PART[i].heat_generation=PART[i].heat_gene_before1;//�Q�d�������v�Z���Ȃ��Ƃ��ł��l���o���Ă���
				}
			}
		}
		

		//�d���͏o��
		if(CON.get_EM_method()!=0)		plot_F(&CON,PART,fluid_number,F,t);


		/////���̑��x����шʒu����
		renewal_u_and_r_in_positive(&CON,PART,fluid_number,t,dt,&Umax,potential,laplacian,g,previous_Un,F);


		for(int D=0;D<DIMENTION;D++) delete [] potential[D];
		for(int D=0;D<DIMENTION;D++) delete [] laplacian[D];
		for(int D=0;D<DIMENTION;D++) delete [] F[D];

		cout<<"�z��͏I�� umax="<<sqrt(Umax)<<"  limit U="<<0.2*mindis/dt<<endl;
		

		///////////////////���q���������̂�INDEX�X�V

		if(CON.get_temporary_r()==ON)//���̈ʒu���v�Z���Ȃ��̂Ȃ�A������freeon�͂��Ȃ��Ă悢
		{
			calc_neighbor_relation(&CON,PART,particle_number,n0_4,fluid_number,out);//�ߗח��q�֌W�v�Z
		}
		
		/*///
		for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH[i];
		delete [] MESH;

		reload_INDEX(&CON,PART,particle_number,INDEX);//�i�q�̗��q���X�V
	
		int **MESH2 = new int *[CON.get_number_of_mesh()];
		for(int i=0;i<CON.get_number_of_mesh();i++)  MESH2[i]=new int [INDEX[i]];
	
		reload_INDEX2(&CON,PART,particle_number,MESH2);
		if(CON.get_temporary_r()==ON)//���̈ʒu���v�Z���Ȃ��̂Ȃ�A������freeon�͂��Ȃ��Ă悢
		{
			unsigned int timeB=GetTickCount();
			if(CON.get_freeon3sw()==OFF)
			{
				if(CON.get_freeon()==1) freeon(&CON,PART,particle_number,n0_4,INDEX,MESH2,&mindis,fluid_number,out);//�\�ʔ���
				else if(CON.get_freeon()==2) freeon2(&CON,PART,particle_number,n0_4,INDEX,MESH2,&mindis,fluid_number,out,t);//�\�ʔ���
				//else if(CON.get_freeon()==4) freeon4(&CON,PART,particle_number,n0_4,INDEX,MESH2,&mindis,n0,fluid_number,e0,out,t);
			}
			if(CON.get_freeon3sw()==ON) freeon3(&CON,PART,particle_number,out);//���q�����x�̂ݍČv�Z
			cout<<"�A��͑O�̗��q�ˑ��֌W�v�Z�I�� time="<<(GetTickCount()-timeB)*0.001<<"[sec]"<<endl;
		}
		////*/

		//���x��
		if(CON.get_T_field()==ON)
		{
			if(CON.get_model_number()==26) 
			{
				set_Q_boundary(&CON, PART, fluid_number, particle_number,dt,t);//���E�̉��x���̂��̂������ȂǂŊ��ɂ��Ƃ܂��Ă���ꍇ�A�����Ή����鋫�E�Ɉʒu���闱�q�ɗ^����
			}
			if(CON.get_T_solution()==0) calc_Temperature(&CON,PART,fluid_number,particle_number,N0,lamda,dt,t);
			if(CON.get_T_solution()==1) calc_temperature_implicity(&CON,PART,fluid_number,particle_number,N0,lamda,dt,t);
		}

		
		///�A���(�����ő��x�E�ʒu�C������)
		if(CON.get_iteration_count()==1)
		{
			if(CON.get_P_twice()==OFF) negative1(&CON,PART,fluid_number,particle_number,out,t, dt, lamda, N0,PND2,n0,Un);
			//if(CON.get_P_twice()==ON) negative1_twice(&CON,PART, fluid_number, particle_number, out, t, dt, lamda, N0, PND2, n0,Un, n0_4);
			if(CON.get_P_twice()>0) negative1_twice(&CON,PART, fluid_number, particle_number, out, t, dt, lamda, N0, PND2, n0,Un, n0_4);
		}
		//if(CON.get_iteration_count()>1)  negative2(&CON,PART,fluid_number,particle_number,t, dt, lamda, n0,INDEX,MESH2,N0,out,Un);
		else if(CON.get_iteration_count()==-1) negative3(&CON,PART, fluid_number, particle_number, t, dt, lamda, N0,NODE_FEM3D,ELEM_FEM3D, out);
		////*/
		
		//�ړ����q���ړ�������ꍇ
		if(CON.get_move_prtcl()==ON) move_particle(&CON,PART,particle_number,fluid_number,dt);

		if(CON.get_modify_position()!=OFF) modify_position(&CON,PART, fluid_number,dt,particle_number);
		
		///MESH����
		//for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH2[i];
		//delete [] MESH2;

		//�A��͌�̍ő呬�x����////////////////
		double umax2=0;//�ő呬�x
		int umax_ID=0;	//�ő呬�x��L���闱�q�ԍ�
		for(int i=0;i<particle_number;i++)
		{ 
			double speed=0;//���q���x
			for(int D=0;D<DIMENTION;D++) speed+=PART[i].u[D]*PART[i].u[D];
			if(speed>umax2) 
			{
				umax2=speed;
				umax_ID=i;
			}	
		}
		//�ő呬�x��L���闱�q�̎��ӗ��q�̑��x���A�ő呬�x�ɔ�ׂ�1m/s�ȏ㏬�����ꍇ�́A�ő呬�x���C������
		double umax3=0;	//�ő呬�x���q�̎��ӗ��q�̂Ȃ��ł̍ő呬�x
		for(int kk=0;kk<PART[umax_ID].N;kk++)
		{
			int j=PART[umax_ID].NEI[kk];
			double speed=0;//���q���x
			for(int D=0;D<DIMENTION;D++) speed+=PART[j].u[D]*PART[j].u[D];
			if(speed>umax3)  umax3=speed;
		}
		umax3=sqrt(umax3);
		umax2=sqrt(umax2);
		if(umax2-umax3>1 && umax2>0)
		{	
			
			cout<<"�ő呬�x�̏C������ umax2="<<umax2<<"->"<<umax3<<endl;
			for(int D=0;D<DIMENTION;D++) PART[umax_ID].u[D]=PART[umax_ID].u[D]*umax3/umax2;
			umax2=umax3;
		}

		cout<<"�A��͏I�� umax="<<umax2<<endl;
		if(umax2>Umax) Umax=umax2;

		distance+=Umax*dt;	//�ő�ړ��������X�V
		if(CON.get_EM_method()>0 && CON.get_EM_distance()>0) cout<<"distance="<<distance<<" "<<CON.get_EM_distance()<<"le="<<CON.get_EM_distance()*CON.get_distancebp()<<endl;
			////////////////////////////////////////

		TIME+=dt;///���ԍX�V
		
		double Pmax=0;
		for(int i=0;i<fluid_number;i++)
		{ 
			if(PART[i].P>=Pmax) Pmax=PART[i].P;
		}
		cout<<"Pmax="<<Pmax<<endl;

		if(t==1)
	{
		ofstream t1("Pmax.dat");		//�����X�e�b�v���A�c��dt�̃O���t
		t1.close();
	}


	ofstream t1("Pmax.dat",ios :: app);
	t1<<t<<" "<<Pmax<<endl;
	t1.close();
		
		cout<<"��͕�������="<<TIME<<" time="<<(GetTickCount()-timet)*0.001<<"[sec]"<<endl;

		ofstream foutc("step.dat");
		foutc<<t<<endl;
		foutc.close();

		///�|�X�g�����F�e�����ʏo�́����ݐ��ɂ��dt����&microAVS�o��
		post_processing(&CON,PART,fluid_number,particle_number,dt,Umax,mindis,t,TIME,time0,count_avs);

		if(t==1 || t%CON.get_interval()==0) output_alldata_AVS(&CON,PART,fluid_number,particle_number,dt,Umax,mindis,t,TIME,time0,count_avs);
		
		//�����𖞂��������q���폜����
		//if(CON.get_check_region()==ON)	delete_particle(CON,PART,&particle_number,&fluid_number,n0_4,t);

		//restart�p̧�ِ���
		post_processing3(&CON,PART,fluid_number,particle_number,t,TIME);

		//�̈�O�̗��q������
		order_sw=check_position(&CON,PART, fluid_number,&particle_number);//�̈�O���q�����m����΁Aorder_sw=ON�ɂȂ�

		//�����𖞂��������q���폜����
		//if(CON.get_check_region()==ON)	delete_particle(CON,PART,&particle_number,&fluid_number,n0_4,t);

		//�������ʂȂ��Ƃ𑪒�E��������֐�.�����͎����Ńv���O�������쐬���Č��߂�
		if(CON.get_check_something()==ON) check_something(&CON,PART, fluid_number, n0, particle_number, t);

		
		cout<<endl;

		for(int D=0;D<DIMENTION;D++) delete [] Un[D];
		//_CrtDumpMemoryLeaks();//���������[�N�̌��o
		
	}

	
	//delete [] INDEX;
	delete [] PND2;
	//for(int D=0;D<DIMENTION;D++) delete [] Un[D];
	for(int D=0;D<DIMENTION;D++) delete [] previous_Un[D];

	clock_t t2=clock();
	cout<<"CPU time="<<(t2-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
	MessageBeep(MB_ICONEXCLAMATION);//��Ƃ̏I����m�点��BEEP��

	//new int;
	_CrtDumpMemoryLeaks();//���������[�N�̌��o

	return 0;
}

//�d�݊֐�
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

//�d�݊֐�
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
	//R:�e�����a
	double r=dis/R;
	double val=1-6*r*r+8*r*r*r-3*r*r*r*r;
	return val;

}


//�������q���x�̌v�Z
double initial_pnd(double r,int dimention,int calc_type)
{
	int size = (int)(r+1);//�v�Z�̈�
	double dis;//����
	double pnd=0;
	int count=0;
	if(dimention==2)
	{
		if(calc_type==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
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
		if(calc_type==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
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
		if(calc_type==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
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
		if(calc_type==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
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
						if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//k����Ȃ�ii��jj�����炷
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
	cout<<"n0��count="<<count<<endl;
	return pnd;  
}

///���׼�ݗp�ϐ��Ɍv�Z�֐�
double calclamda(mpsconfig *CON)
{
	int dimention=CON->get_dimention();	//��͎���
	int Ini_place=CON->get_model_set_way();	//�������q�z�u���@�@0=���� 1=�ז�
	double R=CON->get_re2();			//���v���V�A���p�e�����a
	int size = (int)(R+1);//�v�Z�̈�
	int count=0;
	double dis;//����      
	double w;
	double pnd=0;
	double lam=0;
	if(dimention==2)
	{
		if(Ini_place==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
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
		else if(Ini_place==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
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
		if(Ini_place==0)				//�����z�u�Ƃ��Đ����i�q���Ƃ����ꍇ
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
		else if(Ini_place==1)				//�����z�u�Ƃ��čז��Z�@���Ƃ����ꍇ
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
						if(j%2!=0) ii+=0.5;//j����Ȃ�ii��0.5�i�q�������炷
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//k����Ȃ�ii��jj�����炷
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
	cout<<"�ɂ�count="<<count<<endl;
	return lam;  
}

///���q�f�[�^�ǂݎ��֐�
void input_particle_data(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int t)
{
	if(t==1)//�ŏ���initial_input.dat����ǂݍ���
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
			//old_A�̏����l��^����
			for(int D=0;D<DIMENTION;D++) PART[i].old_A[D]=0;

		}
		fin.close();	
	}
	///////////////////////*/
	
	if(t!=1) //2STEP�ȍ~��mps_input.dat����ǂݍ���
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

//���q���J�E���g�֐� & ���ёւ�
void calc_numbers_of_particles_and_change_the_order(mpsconfig *CON,vector <mpsparticle> &PART,int particle_number,int *fluid_number,int *out,int *order_sw)
{
	//�e���q�����J�E���g
	int fluid_num=0;						//���̗��q��
	int inwall_num=0;
	int out_num=0;
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==FLUID) fluid_num++;
		else if(PART[i].type==INWALL) inwall_num++;
	}
	
	out_num=fluid_num+inwall_num;	//fluid_number<=i<out��INWALL�Q�ŁAout<=i��OUTWALL�Q�ɂȂ�B
	*out=out_num;
	*fluid_number=fluid_num;

	//���ёւ�
	//mpsparticle PART_temp;			//���ёւ��p���q�N���X
	for(int i=0;i<fluid_num;i++)
	{
		if(PART[i].type!=FLUID)
		{
			cout<<"���ёւ��K�v����"<<endl;
		}
	}
	*order_sw=OFF;
}

//INDEX�X�V�֐�
void reload_INDEX(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int *INDEX)
{       
	///����:�e�i�q�̂Ȃ��Ɋ܂܂�闱�q���𐔂���  �i�q�ԍ��͂Q�����ł͍�����0�BX�����ɂ�{�P�ŁA�E��ōő�(X_mesh*Y_mesh) �R�����ł͂y�����ɂ������Ă���
	
	double width=CON->get_distancebp()*CON->get_dx();		//�i�q��
	for(int i=0;i<CON->get_number_of_mesh();i++) INDEX[i]=0;
	for(int i=0;i<particle_number;i++)
	{
		int X=(int)((PART[i].r[A_X]-CON->get_minX())/width);//X�����ɉ��ڂ̊i�q�� 
		int Y=(int)((PART[i].r[A_Y]-CON->get_minY())/width);//Y�����ɉ��ڂ̊i�q��
		int Z=(int)((PART[i].r[A_Z]-CON->get_minZ())/width);//Z�����ɉ��ڂ̊i�q��
		int number=Z*CON->get_X_mesh()*CON->get_Y_mesh()+Y*CON->get_X_mesh()+X;//���qi���܂ފi�q�̔ԍ�
        PART[i].index=number;
		INDEX[number]++;	
	}
}

//INDEX�X�V�֐����̂Q MESH�ɗ��q�ԍ����i�[����
void reload_INDEX2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int **MESH)
{
	int *count=new int [CON->get_number_of_mesh()];
	for(int i=0;i<CON->get_number_of_mesh();i++) count[i]=0;//������

	for(int i=0;i<particle_number;i++)
	{
		int number=PART[i].index;	
        MESH[number][count[number]]=i;
		count[number]++;
	}
	delete [] count;
}

//�@���x�N�g���쐬�֐�
void direct_f(mpsconfig *CON,vector<mpsparticle> &PART,int i,double *direct[DIMENTION])
{
	double R=CON->get_re3()*CON->get_distancebp();//�@���޸�ٌv�Z�ɗ��p����e�����a

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

//�@���x�N�g���쐬�֐�
void direct_f2(mpsconfig *CON,vector<mpsparticle> &PART,int i,double *direct[DIMENTION])
{
	double R=CON->get_re3()*CON->get_distancebp();//�@���޸�ٌv�Z�ɗ��p����e�����a

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

//////direct�p���q�����x����
double pnd_for_direct(mpsconfig *CON,vector<mpsparticle> &PART,double x,double y,double z,double R,int i)
{
	///���qi�̈ʒu�Ƃ͂����̂ŁAMESH���g�p�����ق����悢
	///���ȏ��ł͌����������Ă��邪�A�����ł͏d�݊֐���p����B
	//R=CON->get_re3()*CON->get_distancebp();
	double spnd=0;
	

	for(int k=0;k<PART[i].N3;k++)
	{       
		int j=PART[i].NEI3[k];
		double X=PART[j].r[A_X]-x;
		double Y=PART[j].r[A_Y]-y;
		double Z=PART[j].r[A_Z]-z;
		double dis=sqrt(X*X+Y*Y+Z*Z);
		//if(dis<R) spnd++;   //���ȏ��ǂ���
		if(dis<R)
		{
			double w=(1-dis/R)*(1-dis/R);
			spnd+=w;
		}
	}
	return spnd;
}

//////direct�p���q�����x����
double pnd_for_direct2(mpsconfig *CON,vector<mpsparticle> &PART,double x,double y,double z,double R,int i)
{
	///���qi�̈ʒu�Ƃ͂����̂ŁAMESH���g�p�����ق����悢
	///���ȏ��ł͌����������Ă��邪�A�����ł͏d�݊֐���p����B
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
			//if(dis<R) spnd++;   //���ȏ��ǂ���
			if(dis<R)
			{
				double w=(1-dis/R)*(1-dis/R);
				spnd+=w;
			}
		}
	}
	return spnd;
}

///�|�X�g�����֐�
void post_processing(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double Umax,double mindis,int t,double TIME, unsigned int time0,int count_avs)
{
	double le=CON->get_distancebp();
	unsigned int timeA=GetTickCount();	//�v�Z�J�n����

	cout<<"�e��������ۯĊJ�n: ";
	//cout<<"�e��������ۯĊJ�n"<<endl;
	
	//_CrtDumpMemoryLeaks();
	
	//AVS�ɗ��q�f�[�^�o��/////////////////////////////////
	
	//cout<<"���qavs�f�[�^�o�͊J�n----";
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
	
		
	///���x���v���b�g
//	plot_speed(CON ,PART,particle_number,fluid_number);
	plot_speed_each(CON ,PART,particle_number,fluid_number,t);

	//////���W��ۯ�/////////////////////////
	//////���W��ۯ�/////////////////////////
	ofstream gnu1("0.dat");//��͏I����̑S���q���W���v���b�g
	ofstream gnu2("suf.dat");//�\�ʗ��q�������v���b�g
	ofstream gnu4("suf2.dat");//���̕\�ʗ��q����
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
			if(PART[i].surface==ON && PART[i].type==FLUID)//���̂����\���������Ƃ�
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

	/////////////////////////////CPU time �Ɖ�͕������Ԃ̊֌W�����v���b�g
	if(t==1)
	{
		ofstream t1("dt.dat");		//�����X�e�b�v���A�c��dt�̃O���t
		t1.close();

		ofstream t2("CPU_time1.dat");//����CPU time�A�c��TIME�̃O���t
		t2.close();

		ofstream t3("CPU_time2.dat");//����CPU time�A�c��step�̃O���t
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
	
	//////////////////////////////���ϗ��q���x&���͂�\��
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

///microAVS�p���q����o�͊֐�
void particle_movie_AVS(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double T)
{
	//t:�^�C���ï�߁@T:��������
	double le=CON->get_distancebp();
	double times=1;

	if(t==1)
	{
		ofstream fout("particle_movie.mgf");

		fout<<"# Micro AVS Geom:2.00"<<endl;
		fout<<CON->get_step()/CON->get_interval()+1<<endl;//microAVS�ɏo�͂��鑍�ï�ߐ��B�t�@�C���o�͂�CON->get_interval()���1��ƍŏ��ɍs���B
		
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
		fout2<<CON->get_step()/CON->get_interval()+1<<endl;//microAVS�ɏo�͂��鑍�ï�ߐ��B�t�@�C���o�͂�CON->get_interval()���1��ƍŏ��ɍs���B
		
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
		fout3<<CON->get_step()/CON->get_interval()+1<<endl;//microAVS�ɏo�͂��鑍�ï�ߐ��B�t�@�C���o�͂�CON->get_interval()���1��ƍŏ��ɍs���B
		
		fout3.close();
	}


	ofstream avs3("particle_movie_trace.mgf",ios :: app);
	avs3<<"step"<<t/CON->get_interval()+1<<endl;
	avs3<<"sphere"<<endl;
	avs3<<"time="<<T<<endl;
	avs3<<"color"<<endl;

	double red,green,blue;	//���q�̐F��\������3���F

/////////���̗��q�̂����A�Ǘ��q�t�߂ɂ�����͕̂\�ʂłȂ��Ă��\��������//////////

	///freeon2�֐��̐���:flag1[i]�̓����ɂ�荂�����B�������A1CPU�Ȃ瑁�����A���CPU�ɂ����񉻂͌�����
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	double d=2;
	if(CON->get_dimention()==3) d=3;

	int out=particle_number;

	if(CON->get_T_field()==ON && CON->get_insulate()==1)
	{
		out=particle_number;	//��f�M�̉��x����v�Z����ۂ́AOUTWALL�����ӗ��q���i�[���Ă����K�v������B�����out=particle_number�Ƃ��邱�Ƃł���ɑΏ�
	}

	int *flag1=new int [out];		//�����t���O�B0:�������@1:�����ς�
	int *flag2=new int [out];		//���̗��q�̂����A�Ǘ��q�ɋ߂����̂�ON�B�o�͗p
	///������
	//#pragmallel for
	for(int i=0;i<out;i++)
	{
		flag1[i]=0;
		flag2[i]=0;
	}
	
	if(CON->get_plot_nearbywall()==ON)
	{
		//INDEX,MESH�쐬
		int *INDEX=new int[CON->get_number_of_mesh()];	//�e�i�q�Ɋ܂܂�闱�q�����i�[
		reload_INDEX(CON,PART,particle_number,INDEX);//�i�q���̗��q���X�V
	
		int **MESH = new int *[CON->get_number_of_mesh()];
		int count=0;
		for(int i=0;i<CON->get_number_of_mesh();i++)
		{       
			count+=INDEX[i];
			MESH[i]=new int [INDEX[i]];
		}
		reload_INDEX2(CON,PART,particle_number,MESH);

		//cout<<"INDEX,MESHreload����"<<endl;
	
		for(int i=0;i<out;i++)
		{  
			////���q������
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
								if(flag1[j]==0 && j!=i)//�܂��������ĂȂ��Ȃ�
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
			flag1[i]=1;//�����I��
		}

		//cout<<"�\�����q���v������"<<endl;

		for(int i=0;i<CON->get_number_of_mesh();i++) delete [] MESH[i];
		delete [] MESH;
		delete [] INDEX;
	}
	

////////////////////////////////////////////////////

	
	if(CON->get_AVS()==0)	//���q�̓���������\���B
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		int num2=0;//�\�����闱�q��
		if(CON->get_dimention()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON->get_dimention()==3) 
		{
			for(int i=0;i<particle_number;i++) 
			{
				//if(PART[i].surface==ON)
				{
					if(CON->get_AVS_HALF()==ON)
					{
						if(PART[i].type==INWALL || PART[i].type==OUTWALL) if(PART[i].r[A_X]>=0) num++;//3�����̏ꍇ�A�������͕̂\�����Ȃ�
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

		////���q�o��
		if(CON->get_dimention()==2)
		{    
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
				else if(PART[i].type==INWALL && PART[i].surface==OFF) {red=1;green=0;blue=0;}
				//else if(PART[i].type==INWALL) {red=1;green=0;blue=0;}
				//else if(PART[i].type==OUTWALL) {red=0.5;green=0.5;blue=0;}
	        	else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
				avs<<CON->get_distancebp()/2<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
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
					else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
					//if(PART[i].type==FLUID && PART[i].T<CON->get_MP()) {red=1;green=0;blue=0;}

					if(CON->get_model_number()==19) //FSW�̏ꍇ
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
							avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
							avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
							avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
							}
						}
						
						if(PART[i].type==FLUID)
						{
							if(PART[i].surface==ON)
							{
								avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
								avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
								avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
							}
							if(PART[i].surface==OFF && flag2[i]==ON)
							{
								avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
								avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
								avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
							}
							//if(PART[i].r[A_Y]<0.5*le && PART[i].r[A_Y]>-0.5*le)
							if(PART[i].r[A_X]<0.5*le && PART[i].r[A_X]>-0.5*le)
							{
								avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
								avs2<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
								avs2<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
								
							}
						}
					}
					if(CON->get_AVS_HALF()==OFF)
					{
						if(PART[i].type==INWALL || PART[i].type==OUTWALL)
						{
							avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
							avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
							avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
						}
						if(PART[i].type==FLUID && PART[i].surface==ON)
						{
							avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
							avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
							avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
						}
					}
				}        
			}
		}
	}

	/*				
	if(CON->get_AVS()==0)	//���q�̓���������\���B
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		if(CON->get_dimention()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].surface==ON || PART[i].type!=FLUID) num++;//3�����̏ꍇ�A�������͕̂\�����Ȃ�
			}
		}
		
		avs<<num<<endl;
		/////////////////////////////////

		////���q�o��
		if(CON->get_dimention()==2)
		{    
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
				avs<<CON->get_distancebp()/2<<" ";//���q�̑傫���o��

				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
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
	        		else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��

					avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}        
			}
		}
	}
	*/
	else if(CON->get_AVS()==1)	//���q�̈��͂��R���^�[�\��
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		if(CON->get_dimention()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON->get_dimention()==3) num=fluid_number;//3�����ł͗��̗��q�݂̂�\��
		
		avs<<num<<endl;
		/////////////////////////////////


		double maxP=0;//���͂̍ő�l�ƍŏ��l�����Ƃ߂�
		double minP=0;
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].P>maxP) maxP=PART[i].P;
			if(PART[i].P<minP) minP=PART[i].P;
		}
		/////maxP,minP�����܂����B

		double width=maxP-minP;
		if(width<0.0001) width=0.0001;
			
		if(CON->get_dimention()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
		
	}
	else if(CON->get_AVS()==2)	//���q�̉��x���R���^�[�\��
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		if(CON->get_dimention()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON->get_dimention()==3) num=fluid_number;//3�����ł͗��̗��q�݂̂�\��
		
		avs<<num<<endl;
		avs3<<num<<endl;
		/////////////////////////////////

		double le=CON->get_distancebp();
		double mass=CON->get_particle_mass();
		double T;//���qi�̉��x
		double width=CON->get_maxT()-CON->get_minT();//���x�̕�

		if(CON->get_dimention()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				T=PART[i].T;
				/*
				if(PART[i].type==FLUID)
				{
					double hs0=mass*CON->get_Cp()*CON->get_MP();//�Z���J�n�_�̃G���^���s�[
					double hs1=hs0+CON->get_latent_H()*mass;//�Z���I���_�̃G���^���s�[
					if(PART[i].h<hs0) T=PART[i].h/mass/CON->get_Cp();
					else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON->get_MP();
					else if(hs1<PART[i].h) T=CON->get_MP()+(PART[i].h-hs1)/mass/CON->get_Cp();
				}
				else 
				{
					/////
					if(CON->get_model_number()==26)//���̏ꍇ�A�ǂ̏ꏊ�ɉ����ĕ����n��������
					{
						double air_mass=CON->get_particle_volume()*1.205;//��C�̏d��
						double air_Cp=1006;//��C�̔�M
						if(CON->get_dimention()==2)
						{
							if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08)//���
							{
								T=PART[i].h/air_mass/air_Cp;//�ő�
							}
							else T=PART[i].h/CON->get_particle_volume()*CON->get_wall_density()/CON->get_wall_Cp();//�ő�
						}
						else if(CON->get_dimention()==3)
						{
						}
					}
					else T=PART[i].h/(CON->get_particle_volume()*CON->get_wall_density())/CON->get_wall_Cp();//�ő�
					
				}
				//////*/
    ///////////////
				double level=(T-CON->get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<PART[i].L*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��

				//�g���[�T�[������������p�̐F�ݒ�
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				//else if(PART[i].r[A_Y]>0 && abs(PART[i].r[A_X])<=0.08) {red=0;green=1;blue=0;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=1;green=0;blue=0;}
				else if(PART[i].type==INWALL && PART[i].surface==OFF) {red=0.5;green=0.5;blue=0;}
	        	else if(PART[i].type==INWALL && PART[i].surface==ON) {red=1;green=0;blue=0;}
				
				else {red=0.5;green=0.5;blue=0;}//�Ǘ��q

				//if(PART[i].type==FLUID && i%50==0) {red=1;green=0;blue=0;}//�g���[�T�[
				if(t==1)
				{
					for(int i=0;i<particle_number;i++)
					{
						PART[i].trace=OFF;
						if(PART[i].type==FLUID && i%100==0) PART[i].trace=ON;
					}
				}
				if(PART[i].trace==ON) {red=0;green=1;blue=0;}//�g���[�T�[

				avs3<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs3<<PART[i].L*times<<" ";//���q�̑傫���o��
	
				avs3<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				double hs0=mass*CON->get_Cp()*CON->get_MP();//�Z���J�n�_�̃G���^���s�[
				double hs1=hs0+CON->get_latent_H()*mass;//�Z���I���_�̃G���^���s�[
				if(PART[i].h<hs0) T=PART[i].h/mass/CON->get_Cp();
				else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON->get_MP();
				else if(hs1<PART[i].h) T=CON->get_MP()+(PART[i].h-hs1)/mass/CON->get_Cp();
				double level=(T-CON->get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//�^�񒆂�1�̕�����
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
				
				avs<<PART[i].L*times<<" ";//���q�̑傫���o��
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
			}
		}
	}
	else if(CON->get_AVS()==3)//�Ǘ��q�͔�\��
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=fluid_number;//�\�����闱�q��
		avs<<num<<endl;
		/////////////////////////////////

		////���q�o��
		if(CON->get_dimention()==2) cout<<"error in AVS() 2D�͔�Ή�"<<endl;      
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
				avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��

				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��        
			}
		}
	}
	else if(CON->get_AVS()==4)//���̕\�ʗ��q�̂�
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		if(CON->get_dimention()==2) num=particle_number;//2�����ł͑S���q��\��
		else if(CON->get_dimention()==3) for(int i=0;i<particle_number;i++) if(PART[i].type==FLUID && PART[i].surface==ON) num++;//3�����̏ꍇ�A�������͕̂\�����Ȃ�	
		avs<<num<<endl;
		/////////////////////////////////

		////���q�o��
		if(CON->get_dimention()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			   
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
				avs<<CON->get_distancebp()/2<<" ";//���q�̑傫���o��

				avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
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
	        		else {red=0.5;green=0.5;blue=0;}//�Ǘ��q
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs<<CON->get_distancebp()/2<<" ";//���q�̑傫���o��

					avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}        
			}
		}
	}
	if(CON->get_AVS()==5)	//�ǂ�����\��.�����̊m�F�p
	{
		//////////////////////////////////�\�����闱�q�����v�Z
		int num=0;//�\�����闱�q��
		if(CON->get_dimention()==2) for(int i=0;i<particle_number;i++) if(PART[i].type==INWALL || PART[i].type==OUTWALL) num++;//2�����ł͑S���q��\��
		else if(CON->get_dimention()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL)
				{
					//if(PART[i].toBEM==MOVE)
					{
						num++;//3�����̏ꍇ�A���͕̂\�����Ȃ�
					}
				}
			}
		}
		
		avs<<num<<endl;
		/////////////////////////////////

		////���q�o��
		if(CON->get_dimention()==2)
		{    
			for(int i=0;i<particle_number;i++)
			{
				if(PART[i].type==INWALL || PART[i].type==OUTWALL)
				{
					if(PART[i].type==INWALL) red=0;green=1;blue=0;//�Ǘ��q
					if(PART[i].type==OUTWALL) red=0.5;green=0.5;blue=0;//�Ǘ��q
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs<<CON->get_distancebp()/2<<" ";//���q�̑傫���o��

					avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
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
	        			if(PART[i].type==INWALL) red=0;green=1;blue=0;//�Ǘ��q
						if(PART[i].type==OUTWALL) red=0.5;green=0.5;blue=0;//�Ǘ��q
			    
						avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
						avs<<CON->get_distancebp()<<" ";//���q�̑傫���o��

						avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
					}
				}        
			}
		}
	}
	else if(CON->get_AVS()==6 && CON->get_dimention()==3)	//����̗��q�̓���������\���B���q�̑I���̓v���O�����𒼐ڕύX���邵���Ȃ�
	{
		int num=0;//�\�����闱�q��
		int num2=0;//�\�����闱�q��
		int *output=new int [particle_number];//ON�Ȃ�o�́@OFF�Ȃ�o�͂��Ȃ�
		int *output2=new int [particle_number];//ON�Ȃ�o�́@OFF�Ȃ�o�͂��Ȃ�
		for(int i=0;i<particle_number;i++) 
		{
				output[i]=OFF;//������
				output2[i]=OFF;//������
		}

		vector<int> ID;//AVS�ɏo�͂��闱�q��id���t�@�C�����炱�̔z��ɓ��͂���
		
		/*///
		if(t==1 && CON->get_restart()==OFF)//�ŏ��̃X�e�b�v���Ƀt�@�C���𐶐�
		{
			ofstream fout2("id_for_AVS.dat");//AVS�ŒǐՂ��闱�q��id���o��
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
		
		
		//���藱�q��id��ǂݍ���(1step���ɂ͏�������ł����ǂݍ��ތ`�ɂȂ�)
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
	
		//output[i]����
		for(int k=0;k<ID.size()-1;k++)//�Ȃ񂩏�̂�肩�����ƁAsize-1���Ȃ��Ƃ��܂������Ȃ��B
		{
			output[ID[k]]=ON;
			num++;
		}///////////////////////////////////*/

		for(int i=0;i<fluid_number;i++)
		{
			//if(PART[i].r[A_Y]<0 && PART[i].r[A_Y]>-0.0075)//XZ���ʐ}
			//if(PART[i].r[A_X]<-le )//XY���ʐ}
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

		//output[i]����
		for(int i=fluid_number;i<particle_number;i++)
		{
			
			if(CON->get_tool_angle()==0)//�c�[����]�Ȃ��̏ꍇ
			{
				//if(PART[i].toBEM==MOVE && PART[i].r[A_Z]<=0.006-0.2*le)//�v���[�u�̂ݕ\��
				//if(PART[i].toBEM==MOVE && abs(PART[i].r[A_Z])<=0.003-0.2*le)//�v���[�u�̂ݕ\��
				//if(PART[i].toBEM==MOVE && PART[i].r[A_Y]<0) 
				if(PART[i].toBEM==MOVE)//�c�[���̂ݕ\�� 
				//if(PART[i].toBEM==MOVE && PART[i].r[A_X]>0)
				//if(PART[i].type==FLUID)//��\�� 
				{
					output[i]=ON;
					num++;
				}
			}

			if(CON->get_tool_angle()>0)//�c�[����x���܂��ɉ�]������B�i�s������+y
			{
				double theta=PI/180*CON->get_tool_angle();//��]����p�x
				
				double z=sin(-theta)*PART[i].r[A_Y]+cos(-theta)*PART[i].r[A_Z];//��]������O�̃c�[���ʒu�����߂�
				
				if(t==1)//�����X�e�b�v�̔z�u�Ŕ��f���A���Ƃ͍ŏ��̔���ɏ]���ĕ\��
				{
					if(PART[i].toBEM==MOVE && z<=0.006-0.2*le)//�v���[�u�̂ݕ\��
					//if(PART[i].toBEM==MOVE && PART[i].r[A_Y]<0) 
					//if(PART[i].toBEM==MOVE)//�c�[���̂ݕ\�� 
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
			
			//�f�ʐ}�poutput�ݒ�
			if(PART[i].toBEM==MOVE)//�c�[���̂ݕ\�� 
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

		////���q�o��
		
		for(int i=0;i<particle_number;i++)
		{
			if(CON->get_model_number()==19) //FSW�̏ꍇ�A�ŏ��̔z�u�Ō��߂����q�̔z�F�������Ƃ�������
			{
				//if(PART[i].type!=FLUID) PART[i].color=3; 
				if(output[i]==ON) 
				{	
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		//else if(PART[i].type==INWALL) {red=0;green=1;blue=0;}//�Ǘ��q
					else {red=0.5;green=0.5;blue=0;}//�Ǘ��q

					if(PART[i].type==FLUID)
					{
						//if(PART[i].materialID==1)  {red=0;green=0;blue=1;}
						//else if(PART[i].materialID==2)  {red=1;green=0;blue=0;}
						if(t==1)
						{
							if(PART[i].r[A_X]>0)  
							{
								red=1;green=0;blue=0;
								PART[i].color=1;//1�Ȃ��
							}
							else 
							{
								red=0;green=0;blue=1;
								PART[i].color=2;//2�Ȃ��
							}
						}
						else
						{
							
							if(PART[i].color==1) {red=1;green=0;blue=0;}
							else if(PART[i].color==2) {red=0;green=0;blue=1;}

						}
					}
					//else PART[i].color=3;
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs<<le<<" ";//���q�̑傫���o��

					avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}    

				if(output2[i]==ON) 
				{
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//�Ǘ��q

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
								PART[i].color=1;//1�Ȃ��
							}
							else 
							{
								red=0;green=0;blue=1;
								PART[i].color=2;//2�Ȃ��
							}
						}
						else
						{
							
							if(PART[i].color==1) {red=1;green=0;blue=0;}
							else if(PART[i].color==2) {red=0;green=0;blue=1;}

						}
					}
			    
					avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs2<<le<<" ";//���q�̑傫���o��

					avs2<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}  
			}
			else 
			{
				if(output[i]==ON) 
				{	
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//�Ǘ��q

					if(PART[i].type==FLUID)
					{
						//if(PART[i].materialID==1)  {red=0;green=0;blue=1;}
						//else if(PART[i].materialID==2)  {red=1;green=0;blue=0;}
						if(PART[i].r[A_X]>0)  {red=1;green=0;blue=0;}
						else {red=0;green=0;blue=1;}
					}
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs<<le<<" ";//���q�̑傫���o��

					avs<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
				}    

				if(output2[i]==ON) 
				{
					if(PART[i].type==FLUID && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if(PART[i].type==FLUID && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//�Ǘ��q

					if(PART[i].type==FLUID)
					{
						//if(PART[i].materialID==1)  {red=0;green=0;blue=1;}
						//else if(PART[i].materialID==2)  {red=1;green=0;blue=0;}
						if(PART[i].r[A_X]>0)  {red=1;green=0;blue=0;}
						else {red=0;green=0;blue=1;}
					}
			    
					avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//���W�o��
			
					avs2<<le<<" ";//���q�̑傫���o��

					avs2<<red<<" "<<green<<" "<<blue<<endl;//�F�o��
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

//���x�v���b�g�֐�
void plot_speed(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int fluid_number)
{
	
	double le=CON->get_distancebp()*0.5;
	double times=CON->get_speedtimes();
	int d=CON->get_dimention();
	int NUM;								//AVS�ɏo�͂��闱�q��
	int startID=0;							//�ŏ��ɏo�͂��闱�q��id
	int num=0;								//���������ϐ�
	int face=CON->get_speed_face();			//3D��͎���speed.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON->get_speed_face_p();	//3D��͎���speed.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;								//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	
	//AVS�o�͗��q��NUM�v�Z
	if(CON->get_speed_plot_particle()==1) NUM=particle_number;	//�S���q�o��
	else if(CON->get_speed_plot_particle()==2) NUM=fluid_number;//���̗��q�̂ݏo��
	else if(CON->get_speed_plot_particle()==3)					//�Ǘ��q�̂ݏo��
	{
		NUM=particle_number; 
		startID=fluid_number;
	}
	
	ofstream vec("speed.dat");//��Α��x
	
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
		//int d1,d2;				//�o�͂ɕK�v�Ȏ���
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
					double x=PART[i].r[d1];//�o�͂Ɋ֗^���鎲 �֋X��A�ϐ�����x�ƂȂ��Ă��邪�A�����Ƃ͌���Ȃ����Ƃɒ���
					double z=PART[i].r[d2];//�o�͂Ɋ֗^���鎲
					double u=PART[i].u[d1];//�o�͂Ɋ֗^���鎲
					double w=PART[i].u[d2];//�o�͂Ɋ֗^���鎲

					double y=PART[i].r[d3];//�o�͂Ɋ֗^���Ȃ���
					double v=PART[i].u[d3];//�o�͂Ɋ֗^���Ȃ���

					double r=sqrt(x*x+y*y);//���_����̋���
			
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
						double x2=COS*x-SIN*y;//��]��̍��W�@�~�����̂�x2�̂݁By2�͂���Ȃ�
						double u2=COS*u-SIN*v;//��]��̑��x�@�~�����̂�u2�̂݁Bv2�͂���Ȃ�
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
	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	if(CON->get_legend_speed()>0) vec<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	vec.close();
	/////////////////////////////

	/////////�d�S�ɑ΂��鑊�Α��x�o��
	if(CON->get_relative_speed()==ON)
	{
		double U=0;								//���ϑ��x
		double V=0;
		ofstream vec2("relative_speed.dat");

		///���ϑ��x�v�Z
		if(d==2) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[A_X]; V+=PART[i].u[A_Y];}
		else if(d==3) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[d1]; V+=PART[i].u[d2];}//d1,d2�͂��łɂ��Ƃ܂��Ă���
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


	if(CON->get_flat_speed_plot()==ON && d==3) //���������̑��x�o��
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
	int NUM;								//AVS�ɏo�͂��闱�q��
	int startID=0;							//�ŏ��ɏo�͂��闱�q��id
	int num=0;								//���������ϐ�
	int face=CON->get_speed_face();			//3D��͎���speed.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double face_p=CON->get_speed_face_p();	//3D��͎���speed.dat�̏o�͖ʂ̍��W
	int d1,d2,d3;								//3D��͎��̏o�͂ɕK�v�Ȏ���
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	
	//AVS�o�͗��q��NUM�v�Z
	if(CON->get_speed_plot_particle()==1) NUM=particle_number;	//�S���q�o��
	else if(CON->get_speed_plot_particle()==2) NUM=fluid_number;//���̗��q�̂ݏo��
	else if(CON->get_speed_plot_particle()==3)					//�Ǘ��q�̂ݏo��
	{
		NUM=particle_number; 
		startID=fluid_number;
	}
	
	ofstream vec("speed.dat");//��Α��x
	
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

		//int d1,d2;				//�o�͂ɕK�v�Ȏ���
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
					double x=PART[i].r[d1];//�o�͂Ɋ֗^���鎲 �֋X��A�ϐ�����x�ƂȂ��Ă��邪�A�����Ƃ͌���Ȃ����Ƃɒ���
					double z=PART[i].r[d2];//�o�͂Ɋ֗^���鎲
					double u=PART[i].u[d1];//�o�͂Ɋ֗^���鎲
					double w=PART[i].u[d2];//�o�͂Ɋ֗^���鎲

					double y=PART[i].r[d3];//�o�͂Ɋ֗^���Ȃ���
					double v=PART[i].u[d3];//�o�͂Ɋ֗^���Ȃ���

					double r=sqrt(x*x+y*y);//���_����̋���
			
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
						double x2=COS*x-SIN*y;//��]��̍��W�@�~�����̂�x2�̂݁By2�͂���Ȃ�
						double u2=COS*u-SIN*v;//��]��̑��x�@�~�����̂�u2�̂݁Bv2�͂���Ȃ�
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
	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	if(CON->get_legend_speed()>0) vec<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	vec.close();

	xmax=-100;						//�o�͗��q�̍ő剡���W
	ymax=-100;						//�o�͗��q�̍ő�c���W
	
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
			//int d1,d2;				//�o�͂ɕK�v�Ȏ���
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
						double x=PART[i].r[d1];//�o�͂Ɋ֗^���鎲 �֋X��A�ϐ�����x�ƂȂ��Ă��邪�A�����Ƃ͌���Ȃ����Ƃɒ���
						double z=PART[i].r[d2];//�o�͂Ɋ֗^���鎲
						double u=PART[i].u[d1];//�o�͂Ɋ֗^���鎲
						double w=PART[i].u[d2];//�o�͂Ɋ֗^���鎲

						double y=PART[i].r[d3];//�o�͂Ɋ֗^���Ȃ���
						double v=PART[i].u[d3];//�o�͂Ɋ֗^���Ȃ���

						double r=sqrt(x*x+y*y);//���_����̋���
			
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
							double x2=COS*x-SIN*y;//��]��̍��W�@�~�����̂�x2�̂݁By2�͂���Ȃ�
							double u2=COS*u-SIN*v;//��]��̑��x�@�~�����̂�u2�̂݁Bv2�͂���Ȃ�
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

		//�}��o��
		xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
		ymax+=4*le;

		if(CON->get_legend_speed()>0) vec2<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	
		vec2.close();///////////////////


		if(CON->get_model_number()==19&&CON->get_process_type()==2)
		{
			face_p=0.0;//�c�[��������ʂ鑬�x��ʓr�\��
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
				//int d1,d2;				//�o�͂ɕK�v�Ȏ���
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
							double x=PART[i].r[d1];//�o�͂Ɋ֗^���鎲 �֋X��A�ϐ�����x�ƂȂ��Ă��邪�A�����Ƃ͌���Ȃ����Ƃɒ���
							double z=PART[i].r[d2];//�o�͂Ɋ֗^���鎲
							double u=PART[i].u[d1];//�o�͂Ɋ֗^���鎲
							double w=PART[i].u[d2];//�o�͂Ɋ֗^���鎲

							double y=PART[i].r[d3];//�o�͂Ɋ֗^���Ȃ���
							double v=PART[i].u[d3];//�o�͂Ɋ֗^���Ȃ���

							double r=sqrt(x*x+y*y);//���_����̋���
			
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
								double x2=COS*x-SIN*y;//��]��̍��W�@�~�����̂�x2�̂݁By2�͂���Ȃ�
								double u2=COS*u-SIN*v;//��]��̑��x�@�~�����̂�u2�̂݁Bv2�͂���Ȃ�
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

			//�}��o��
			xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
			ymax+=4*le;

			if(CON->get_legend_speed()>0) vec3<<xmax<<" "<<ymax<<" "<<CON->get_legend_speed()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	
			vec3.close();///////////////////
		}
	}
	/////////////////////////////
}

void plot_F(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *F[DIMENTION],int t)
{
	int d=CON->get_dimention();
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
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


	////�t�@�C���o��
	ofstream fp("Fn.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();
	double times=CON->get_times()*CON->get_density()/CON->get_particle_mass();
	double cross_section=CON->get_speed_face_p();

    for(int i=0;i<fluid_number;i++)//���̐ߓ_�̂ݏo��
    {
		//if(PART[i].r[A_Y]>-le*0.5&& PART[i].r[A_Y]<le*0.5) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F[A_X][i]*times<<" "<<F[A_Z][i]*times<<endl;
		if(PART[i].r[A_Y]>cross_section-le*0.5&& PART[i].r[A_Y]<cross_section+le*0.5) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].F[A_X]*times<<" "<<PART[i].F[A_Z]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
	}

	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	//�}��o��
	if(CON->get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	

    fp.close();///////////////////

	////�t�@�C���o��//(�X���b�g��ʂ�f��)
	xmax=-100;						//�o�͗��q�̍ő剡���W
	ymax=-100;						//�o�͗��q�̍ő�c���W
	ofstream fps("Fnslit.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=0;i<fluid_number;i++)//���̐ߓ_�̂ݏo��
    {
		if(PART[i].r[A_Y]>-le*0.5+sin(PI/24)*PART[i].r[A_X] && PART[i].r[A_Y]<le*0.5+sin(PI/24)*PART[i].r[A_X])
		{
			fps<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<(PART[i].F[A_X]*cos(PI/24)+PART[i].F[A_Y]*sin(PI/24))*times<<" "<<PART[i].F[A_Z]*times<<endl;
		}
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		
	}

	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	//�}��o��
	if(CON->get_legend_F()>0) fps<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
    fps.close();
	///////////////////

	////�t�@�C���o��(Z�f��)
	xmax=-100;						//�o�͗��q�̍ő剡���W
	ymax=-100;						//�o�͗��q�̍ő�c���W
	ofstream fpz("Fnz.dat");
	//double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le*CON->get_FEMtimes();

    for(int i=0;i<fluid_number;i++)//���̐ߓ_�̂ݏo��
    {
		if(PART[i].r[A_Z]>-le*0.5+0.14125 && PART[i].r[A_Z]<le*0.5+0.14125) fpz<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].F[A_X]*times<<" "<<PART[i].F[A_Y]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
	}

	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	//�}��o��
	if(CON->get_legend_F()>0) fpz<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
    
    fpz.close();///////////////////

	//if(t=1 || t%10==0)
	/////////////���X�e�b�v���Ƃɓd���͏o��
	xmax=-100;						//�o�͗��q�̍ő剡���W
	ymax=-100;						//�o�͗��q�̍ő�c���W
	if((t-1)%CON->get_EM_interval()==0)
	{
		char filename[20];
		sprintf_s(filename,"Fn%d.dat", t);
		ofstream fp(filename);

		for(int i=0;i<fluid_number;i++)//���̐ߓ_�̂ݏo��
		{
			if(PART[i].r[A_Y]>-le*0.5&& PART[i].r[A_Y]<le*0.5) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].F[A_X]*times<<" "<<PART[i].F[A_Z]*times<<endl;	
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}

		//�}��o��
		xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
		ymax+=4*le;

		if(CON->get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��
	
		fp.close();///////////////////
	}
	
}

//���q�����o�͂���֐�
void output_alldata_AVS(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double Umax,double mindis,int t,double TIME, unsigned int time0,int count_avs)
{
	//�Q�l�ɂ��Ă��鏑����microAVS�̃w���v�ł��Ȃ��̃f�[�^�́H���u��\���i�q�^�f�[�^�i�A�X�L�[�j�̏����v

	if(t==1)
	{
		ofstream fout("alldata_part.dat");

		fout<<"# Micro AVS Geom:2.00"<<endl;
		fout<<CON->get_step()/CON->get_interval()+1<<endl;//microAVS�ɏo�͂��鑍�ï�ߐ��B�t�@�C���o�͂�CON->get_interval()���1��ƍŏ��ɍs���B
		
		fout.close();
	}


	ofstream fp("alldata_part.dat",ios :: app);
	fp<<"step"<<t/CON->get_interval()+1<<endl;
	fp<<"sphere"<<endl;
	fp<<"time="<<TIME<<endl;
	fp<<"color"<<endl;
	
	 
	//�ߓ_�ԍ��Ƃ��̍��W�̏o�� 
	//for(int i=0;i<particle_number;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
	/*/�v�f�ԍ��Ɨv�f�`��̎�ށA�����ėv�f���\������ߓ_�ԍ��o��
	for(int i=0;i<nelm;i++)
	{
		fp<<i+1<<"  0 tri ";
		for(int j=0;j<3;j++)	fp<<ELEM[i].node[j]+1<<" ";
		fp<<endl;
	}
	////*/


	////*/

	//fp<<"2 3"<<endl;//�ߓ_�̏��ʂ�2�ŁA�v�f�̏��ʂ�3�Ƃ������ƁB
	//fp<<"6 0"<<endl;//�ߓ_�̏��ʂ�8�ŁA�v�f�̏��ʂ�0�Ƃ������ƁB
	//fp<<"2 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	//fp<<"6 1 1 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	//fp<<"6 1 1 1 1 1 1"<<endl;	//���̍s�̏ڍׂ̓w���v���Q��
	////fp<<"first_id"<<endl;
	//fp<<"speed,m/s"<<endl;
	//fp<<"P,N/m^2"<<endl;
	//fp<<"h,W/m^3"<<endl;
	//fp<<"surface"<<endl;
	//fp<<"color"<<endl;//���q�̕\���F

	if(t==1) for(int i=0;i<particle_number;i++) PART[i].firstID=i;//�ŏ��̏�Ԃł̗��q�ԍ����L��

	//�e�ߓ_�̏��l����
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
	
	//�v�f���o�́@�v�f����microAVS�̉������\�b�h�o�[���́A�u�v�f�f�[�^�̓h��Ԃ��v�������Ό����
	for(int i=0;i<nelm;i++) fp<<i+1<<"  "<<ELEM[i].potential<<" "<<ELEM[i].En<<" "<<ELEM[i].Fn<<endl;
	*/
	cout<<"OK"<<endl;
	fp.close();
}

//���ݐ�
void culan(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double *dt,double mindis,double Umax,double *g)
{
	double C=CON->get_curan();		//���ݏ�����
	double le=mindis;				//�ŒZ���q�ԋ���
	//double le=CON->get_distancebp()/sqrt(2.0);	//���ϗ��q�ԋ���

	///////////�N�[������//////////////////
	if(C>0)
	{    
		double newdt=*dt;	//�V����dt
		
		if(Umax!=0) newdt=C*le/Umax;
		if(newdt>CON->get_dt()) *dt=CON->get_dt();
		else *dt=newdt;

		//if(*dt<=CON->get_dt()/10) *dt=CON->get_dt()/10;

		if(*dt!=CON->get_dt()) cout<<"���ݐ����� dt="<<*dt<<endl;
		else cout<<"�ݒ�l�@dt="<<*dt<<endl;
	}

	///�g�U���̐��m�Ȓ�`�𒲂ׂď����Ȃ���
	if(CON->get_vis()!=0 && CON->get_vis_calc_type()==POSITIVE)
	{
		if(*dt>0.25*le*le/CON->get_vis()) cout<<"�g�U���ᔽ"<<endl;
	}
	/////////////////////////////
	
}
/*//���ݐ�
void culan(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double *dt,double mindis,double Umax,double *g)
{
	double C=CON->get_curan();		//���ݏ�����
	double le=mindis;				//�ŒZ���q�ԋ���

	///////////�N�[������//////////////////
	if(C>0)
	{      
		double newdt=*dt;	//�V����dt
		
		if(Umax!=0) newdt=C*le/Umax;
		if(newdt>CON->get_dt()) *dt=CON->get_dt();
		else *dt=newdt;
		
		if(*dt!=CON->get_dt()) cout<<"���ݐ����� dt="<<*dt<<endl;
	}

	///�g�U���̐��m�Ȓ�`�𒲂ׂď����Ȃ���
	if(CON->get_vis()!=0 && CON->get_vis_calc_type()==POSITIVE)
	{
		if(*dt>0.25*le*le/CON->get_vis()) cout<<"�g�U���ᔽ"<<endl;
	}
	/////////////////////////////
	
}
///*/

///���v���V�A���v�Z�֐�
void u_laplacian_f(mpsconfig *CON,vector<mpsparticle> &PART,double *laplacian[DIMENTION],double n0,double lamda,int fluid_number,double dt)
{
	cout<<"���x�g�U���v�Z-------";

	double le=CON->get_distancebp();
	double R=CON->get_re2()*le;						//���׼�ݗp�e�����a
	int d=CON->get_dimention();						//����
	unsigned int timeA=GetTickCount();				//�v�Z�J�n����

	///�g�U���Ɉᔽ���Ă��Ȃ����m�F
	double limit=0.25;
	if(d==3) limit=0.125;
	if(dt>limit*le*le/CON->get_vis()) cout<<"�g�U���ᔽ dt<"<<limit*le*le/CON->get_vis()<<" ";
	/////////////*/
   
	for(int i=0;i<fluid_number;i++)
	{
		double lam=0;	//�e���q�̐��m�ȃ�
		double W=0;	//�e���q�̐��m�ȗ��q�����x
		for(int D=0;D<DIMENTION;D++) laplacian[D][i]=0.0;//������
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k]; 
			if(CON->get_wall_adheision()==0)//�ذ�د��
			{
				if(PART[j].type==FLUID)
				{ 
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);//���q�Ԃ̋���
					
					double w=kernel(R,dis);//�d�݊֐�
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
			else if(CON->get_wall_adheision()==1 || CON->get_wall_adheision()==2)//�ݽد��
			{
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);//���q�Ԃ̋���
				
				double w=kernel(R,dis);//�d�݊֐�
				W+=w;
				lam+=dis*dis*w;
				if(CON->get_laplacian()==0 || CON->get_laplacian()==1)
				{
					for(int D=0;D<DIMENTION;D++)
					{
						double u=PART[j].u[D]-PART[i].u[D];
						if(CON->get_wall_adheision()==2)
						{   //�{���̒�`�ǂ�����ݽد��
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
						{   //�{���̒�`�ǂ�����ݽد��
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

///�S�����A��͊֐�
void visterm_negative(mpsconfig *CON,vector<mpsparticle> &PART,double *laplacian[DIMENTION],double n0,double lamda,int fluid_number,int particle_number,double dt,int t)
{
	//���݂�3�c�̂ݑΉ� //20130505 �C���B2D�ɂ��K�p�\�̂͂�
	//�����ׂ�����	ui(t+1)-ui(t)=��t*v*(2d/��n0)��(uj-ui)w
	//				�̃�t*v*(2d/��n0)��(uj-ui)w-ui(t+1)=-ui(t)
	//				�̃�(uj-ui)w-��n0/(2vd��t)ui(t+1)=-��n0/(2vd��t)ui(t)
	//				��v��(uj-ui)w-��n0/(2d��t)ui(t+1)=-��n0/(2d��t)ui(t)

	cout<<"���x�g�U���v�Z---";

	double R=CON->get_re2()*CON->get_distancebp();	//���׼�ݗp�e�����a
	int d=CON->get_dimention();						//����
	unsigned int timeA=GetTickCount();				//�v�Z�J�n����
	int count=0;
	//int pn=fluid_number*d;							//���m��:���q���~���xD(����)���� ���݂͂R�c�̂ݑΉ�
	int pn=fluid_number;
	double co=lamda*n0/(2*d*dt);///CON->get_vis();	//�v�Z�ɂ悭�����W��
	double co2=2*d*dt/(lamda*n0)*CON->get_vis();
	double vis0=CON->get_vis();

	double *vis=new double [fluid_number];			//�e���q�̓��S�����v�Z
	double *B   = new double[pn];					//���s��

	double *B_vec[DIMENTION];
    for(int D=0;D<DIMENTION;D++) B_vec[D]=new double [fluid_number];//�e�����̉��s��

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			laplacian[D][i]=0;		//������
			B_vec[D][i]=0;
		}
	}
	///�e���q�̓��S�����v�Z

	if(CON->get_model_number()==19)
	{
		calc_vis_value(CON,PART,fluid_number,vis,dt,t,particle_number);//FSW���f���̏ꍇ
		for(int i=0;i<particle_number;i++)	PART[i].vis=0;
		for(int i=0;i<fluid_number;i++)	PART[i].vis=vis[i];
		if(t==0||t%CON->get_interval()==0)	output_viscousity_avs(CON,PART,t,particle_number,fluid_number);
	}
	else for(int i=0;i<fluid_number;i++) vis[i]=vis0;

	if(CON->get_temperature_depend()==ON) calc_physical_property(CON,PART,fluid_number,vis,particle_number,5);//���S���̉��x�ˑ�
	//////////////////////////////////////////////////////////////////

	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) laplacian[D][i]=0;//������

	//���s��B[]�쐬
	for(int i=0;i<fluid_number;i++)//���m���̏��Ԃ�PART[i].u[A_X],PART[i+1].u[A_X]�E�E�EPART[i].u[A_Y],PART[i+1].u[A_Y],
	{
		for(int D=0;D<d;D++)
		{
			B_vec[D][i]=-co*PART[i].u[D]*PART[i].PND2/n0;//n0�ł͂Ȃ�wi�g�p
		}
	}
	/////*/

	int number=0;			//�W���s��̔�[���v�f��
	for(int i=0;i<fluid_number;i++)
	{
		int num=1;//���qi���ӂ̗��̗��q�� �����l��1�Ȃ͎̂������g���Ă��Ă��邩��
		for(int k=0;k<PART[i].N2;k++)
		{
			int j=PART[i].NEI2[k];
			if(j<fluid_number) num++;
		}
		//
		number+=num;
	}///number�����Ƃ܂���
	
    double *val = new double [number];
	int *ind = new int [number];//��[���v�f�̗�ԍ��i�[�z��
	int *ptr = new int [pn+1];//�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	
	/////////////////////val,ind ,ptr�ɒl���i�[
	int index=0;
	for(int n=0;n<pn;n++)
	{
		//int D=n/fluid_number;//n�Ԗڂ̖��m����u[D]�Ɋւ��関�m���@���̂����݂��Ȃ����f���ł��̊֐����񂳂Ȃ��悤����
	    
		ptr[n]=index;
		int i=n;//���Ԗڂ̖��m���ɊY�����闱�q�ԍ�
		
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

			if(j<fluid_number)//���qj�����̂Ȃ�
			{
				double v=2*vis[i]*vis[j]/(vis[i]+vis[j]);
				val[index]=w*v;
				ind[index]=j;
				index++;
				AA+=w*v;
			}
			else
			{
				for(int D=0;D<d;D++) B_vec[D][n]-=w*PART[j].u[D]*vis[i];//���qj�����̂łȂ��Ȃ�
				AA+=w*vis[i];
			}
	    }
		//val[KK]=-AA-co;
		val[KK]=-AA-co*PART[i].PND2/n0;//n0�ł͂Ȃ�wi�g�p
	}
	ptr[pn]=number;//���̍Ō�̍s���Ȃ��ƁA�C�ӂ�n��for(int j=ptr[n];j<ptr[n+1];j++) �̂悤�Ȃ��Ƃ��ł��Ȃ�
	////////////////////*/	 


	///ICCG�@�łƂ��ꍇ�A�W���s���������ƕ��ѕς����Ȃ��Ă͂Ȃ�Ȃ�
	if(CON->get_vis_solver()==1)
	{
		//�����ō����s��͑Ίp���������ׂĕ��Ȃ̂ŁA����𐳂ɂ��邽�߂ɁA�W���s��Ɖ��s���-1��������
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
	double *X=new double[pn];		//�s��̓����i�[
	double *AP = new double [pn];
	double *P = new double [pn];

	for(int D=0;D<d;D++)
	{
		count=0;
		for(int n=0;n<pn;n++) B[n]=B_vec[D][n];
		/////////////////////////�����l//////////////////
		if(CON->get_initial_u()==OFF)
		{
			for(int n=0;n<pn;n++) 
			{
				 X[n]=0;
				 r[n]=B[n];
				 P[n]=r[n];
			}
		}
		else if(CON->get_initial_u()==ON) //�����l��^����ꍇ�@����ȂɌ��ʂȂ�����
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

		
		//if(CON->get_vis_solver()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,1e-8); //CG�@�ɂ��s�������
		if(CON->get_vis_solver()==0) CG_method(CON,r,P,AP,val,ind,ptr,pn,X,&count,1e-8); //CG�@�ɂ��s�������
		else if(CON->get_vis_solver()==1) iccg(CON,val,ind,ptr,pn,B,number,X,r,P,1e-8,&count);//ICCG�ɂ��s��v�Z�J�n
		else if(CON->get_vis_solver()==2) MRTR(CON,r, pn,X,&count,1e-8,val,ind,ptr);

		for(int n=0;n<pn;n++)
		{	
			//PART[i].u[D]=XX[n];		//���������Ƒ��x�������Ō��߂Ă��܂�
			laplacian[D][n]=(X[n]-PART[n].u[D])/(dt*vis0);//����������renewal�֐��Ȃ��ő��x�X�V����`�ɂȂ�
		}
		cout<<"������:"<<count<<"  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		if(count==1) for(int i=0;i<fluid_number;i++) laplacian[D][i]=0;//�����񐔂��P�Ƃ����͖̂��炩�ɃG���[�Ȃ̂ŁA���̒l�͎g�p���Ȃ�
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

//CG�@
void CG_method(mpsconfig *CON,double *r,double *P,double *AP,double *val,int *ind,int *ptr,int pn,double *X,int *countN,double EP)
{
	cout<<"CG�@�X�^�[�g------";
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double rr=0;
	double E=1;//�덷
	double alp,beta;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];
	
	if(CON->get_omp_P()==OFF)//�ʏ��
	{
		while(E>EP)// EP=CON->get_CGep();//��������(convergence test)
		{
			count++;
			//////////////alp�����߂�
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
		
			//////////////// ���X�V�@X(k+1)=X(k)+alp*P
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

			//////////////////�덷
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
			////////////////////////
			
			///////////////////// P=r+beta*P
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];
		}
	}
	else if(CON->get_omp_P()==ON)//openMP���g�p����ꍇ
	{
		while(E>EP)
		{
			count++;
			//////////////alp�����߂�
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
			E=0;//�덷
			beta=1.0/rr;
			rr=0;
			#pragma omp parallel for reduction(+:rr)
			for(int n=0;n<pn;n++) 
			{
				X[n]+=alp*P[n];// ���X�V�@X(k+1)=X(k)+alp*P
				r[n]-=alp*AP[n];// r=r-alp*AP
				rr+=r[n]*r[n];
			}
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
		
			beta=beta*rr;///beta
			
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];/// P=r+beta*P
		}
	}
	*countN=count;//�����񐔂�n��
}

///ICCG�@
void iccg(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X,double *r,double *P,double EP,int *count2)
{
	//val :�[���v�f�̒l
	//ind:��[���v�f�̗�ԍ��i�[�z��
	//ptr:�e�s�̗v�f��val�̉��Ԗڂ���͂��܂�̂����i�[
	//X[n]:��

	double accel=0.87;//CON->get_CGaccl();//�����t�@�N�^
	
	int num2=0;//�Ίp�������܂ށA���O�p�s�񂾂����l���ɂ��ꂽ��[���v�f��
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *val2=new double [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//����̧��
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//��������Ă����Ȃ��ƁA�Ō��(int m=ptr2[k];m<ptr2[k+1];m++)�݂����Ȃ��Ƃ��ł��Ȃ�

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];//�[���v�f�̒l
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];//��[���v�f�̍s�ԍ��i�[�z��
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg�@
	double alp,beta;
	double rLDLt_r;
	double E=1;//�덷
	double *AP = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *D1 = new double [pn];//D�s��
	
	/////�s���S�R���X�L�����
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
	        int i=ind2[m];//��ԍ�
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
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i�s�ڂ̂Ȃ������̈�v������̂�T���Ă���B������Ԃ��H
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
	///�s���S�ڽ����������/////////*/

	///�����ɂ����z��ɒl����
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k�s�ڂ̔�O�v�f
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j 
		else
		{
		    double sum=0;
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//���i3.78�j
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]�����Ƃ܂����B
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

	cout<<"ICCG�@:���m��="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=EP;//��������
	rLDLt_r=0;
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//�ŏ���rLDLt_r���������ŋ��߂�
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alp�����߂�
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
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P �X�V��̏ꏊ
			r[n]-=alp*AP[n];// r=r-alp*AP       �X�V��̎c��
			E+=r[n]*r[n];						//�X�V��̌덷
		}
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]�����Ƃ߁ALDLt_r[i]�����Ƃ߂�B
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //���i3.77�j �V
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//�Ίp�����͏�������ptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//���i3.78�j
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]�����Ƃ܂����B
	
		/////////LDLt_r[i]�����߂�
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
	//cout<<"������="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
		
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
	*count2=count;//�����񐔂��i�[���ĕԂ�
}

///�S�����v�Z�֐�
void calc_viscous_term(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double N0,double *laplacian[DIMENTION],double lamda,int t)
{
	int d=CON->get_dimention();

	if(CON->get_nensei()>0)
	{
		if(CON->get_vis_calc_type()==POSITIVE) u_laplacian_f(CON,PART,laplacian,N0,lamda,fluid_number,dt);//�z��@
		if(CON->get_vis_calc_type()==NEGATIVE) //�A��@
		{
			int flag=ON;
			if(fluid_number==0) flag=OFF;			//���̂����݂��Ȃ���Όv�Z���Ȃ�
			if(t==1 && CON->get_restart()==OFF)//�S���q���x��0�̏ꍇ�A�b�f�@��alpha=���ƂȂ�G���[�ɂȂ�B������������
			{
				flag=OFF;
				for(int i=0;i<particle_number;i++) for(int D=0;D<d;D++) if(PART[i].u[D]!=0) flag=ON;
			}
			else if(t==1 && CON->get_restart()==ON && CON->get_set_zero_speed()==ON) flag=OFF;
			if(flag==ON)  visterm_negative(CON,PART,laplacian,N0,lamda,fluid_number,particle_number,dt,t);
			if(flag==OFF) for(int D=0;D<CON->get_dimention();D++) for(int i=0;i<fluid_number;i++) laplacian[D][i]=0;//������
		}
	}
	else for(int D=0;D<d;D++) for(int i=0;i<fluid_number;i++) laplacian[D][i]=0;//�v�Z���Ȃ��ꍇ���������������Ă���		
}







//�K�E�X�̏����@ ���͍ŏI�I��B�̂Ȃ���
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

/////���̑��x����шʒu����
void renewal_u_and_r_in_positive(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double dt,double *Umax,double **potential,double **laplacian,double *g,double **previous_Un,double **F)
{

	double U=0;						//�ő呬�x
	double vis=CON->get_vis();
	double *vis_T=new double [fluid_number];
	double mass=CON->get_particle_mass();	//���q�̎���
	double *mass_T=new double [fluid_number];	//���q�̎���
	int d=CON->get_dimention();
	int sw=CON->get_temporary_r();	//�n�m�Ȃ牼�̈ʒu���v�Z����
	double T0=CON->get_initialT();//����x
	double alp=CON->get_CTE();//���c���W��
	

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			mass_T[i]=mass;//���x�ˑ����̍l���̗L���Ɋւ�炸�A�ЂƂ܂�config�̒l�ɏ���������Bcalc_physical_property�ɍs���Ă��ύX����Ȃ��ꍇ�����邽�߁B
			vis_T[i]=vis;
		}
	}

	double *old_U[DIMENTION];
	for(int D=0;D<DIMENTION;D++) old_U[D]=new double [fluid_number];//�ύX�O�̑��x���L�����Ă���

	if(t==1) for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENTION;D++) previous_Un[D][i]=0;//t=1�̂Ƃ��͏�����   

	if(CON->get_temperature_depend()==ON)
	{
		calc_physical_property(CON,PART,fluid_number,mass_T,fluid_number,1);//���x�̉��x�ˑ��@���̎��_��mass_T�ɓ���͖̂��x(kg/m3)
		for(int i=0;i<fluid_number;i++) mass_T[i]=mass*mass_T[i]/CON->get_density();
		//for(int i=0;i<fluid_number;i++) mass_T[i]=mass;

		//calc_physical_property(CON,PART,fluid_number,vis_T,fluid_number,5);//���S���̉��x�ˑ� �S�����̉A��͂ɂāA
	}
		
			

	//potential[D][i]���ꍇ�ɂ���Ă̓[���ɏ���������
	if(CON->get_dir_for_P()==1 || CON->get_dir_for_P()==3) //�\�ʗ��q�̕\�ʒ��͈͂��͒l�Ƃ��Čv�Z����Ă���̂�,�����ł͍l�����Ȃ��悤����������
	{
		for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) potential[D][i]=0;
	}
	//////////////////////*/

	/////////////���x�X�V
	for(int i=0;i<fluid_number;i++)
	{        
		double speed=0;//���q���x
		for(int D=0;D<d;D++)
		{   
			old_U[D][i]=PART[i].u[D];
			
			//if(CON->get_temperature_depend()==OFF) PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]+F[D][i]/mass);
			//if(CON->get_temperature_depend()==ON) PART[i].u[D]+=dt*(vis_T[i]*laplacian[D][i]+potential[D][i]+g[D]+F[D][i]/mass_T[i]);
			if(CON->get_temperature_depend()==OFF) 
			{
				if(CON->get_buoyant()==OFF) PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]+PART[i].F[D]/mass);
				if(CON->get_buoyant()==ON)//���͂̃u�V�l�X�N�ߎ��B
				{
					double T0=CON->get_initialT();//����x
					double alp=CON->get_CTE();//���c���W��
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1+(mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0))+PART[i].F[D]/mass);//���ꂪ�������͂�
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass_T[i])+PART[i].F[D]/mass);//���������HmassT�Ŋ��鍪�����Ȃ��B
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass)+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
				}

					
					
			}
			if(CON->get_temperature_depend()==ON) 
			{
				if(CON->get_buoyant()==OFF)//�u�V�l�X�N�ߎ���p���Ȃ��B���x�ω�������ƘA���̎������x���U0�ɂ͂Ȃ�Ȃ��Ȃ邩�炱������g���Ȃ炢�낢�뒼���K�v����
				{
					PART[i].u[D]+=dt*(vis_T[i]*laplacian[D][i]+potential[D][i]+g[D]+PART[i].F[D]/mass_T[i]);
				}
				if(CON->get_buoyant()==ON)//���͂̃u�V�l�X�N�ߎ��B
				{
					double T0=CON->get_initialT();//����x
					double alp=CON->get_CTE();//���c���W��
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1+(mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0))+PART[i].F[D]/mass);//���ꂪ�������͂�
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass_T[i])+PART[i].F[D]/mass);//���������HmassT�Ŋ��鍪�����Ȃ��B
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*(1-alp*(PART[i].T-T0)/mass)+PART[i].F[D]/mass);
					//PART[i].u[D]+=dt*(vis*laplacian[D][i]+potential[D][i]+g[D]*((mass_T[i]-mass)/mass_T[i])+PART[i].F[D]/mass);
				}
			}
			
			//PART[i].u[D]=previous_Un[D][i]+dt*(vis*laplacian[D][i]+potential[D][i]+g[D]);//�^�Ƃі@
			speed+=PART[i].u[D]*PART[i].u[D];
		}
		if(speed>U) U=speed;	
	}
	*Umax=U;


	//�\�ʂ������瑬�x�̍X�V����߂�
	//if(CON->get_fix_surface()==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) if(PART[i].surface==ON) PART[i].u[D]=0;

	//�ʒu�X�V
	if(sw==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*0.5*(PART[i].u[D]+old_U[D][i]);//��`��
	
	//previous_Un(1step�O�̑��x���)�̒l���X�V
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++)previous_Un[D][i]=old_U[D][i];
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++)previous_Un[D][i]=PART[i].u[D];
	//for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++)previous_Un[D][i]=0;
	
	
	if(CON->get_buoyant()==ON)
	{
		ofstream vec("gf.dat");//��Α��x
		double le=CON->get_distancebp()*0.5;
		double times=0.1;
		int d=CON->get_dimention();
		int NUM=0;								//AVS�ɏo�͂��闱�q��
		int startID=0;							//�ŏ��ɏo�͂��闱�q��id
		int num=0;								//���������ϐ�
		int face=CON->get_speed_face();			//3D��͎���speed.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
		double face_p=CON->get_speed_face_p();	//3D��͎���speed.dat�̏o�͖ʂ̍��W
		int d1,d2,d3;								//3D��͎��̏o�͂ɕK�v�Ȏ���
		double xmax=-100;						//�o�͗��q�̍ő剡���W
		double ymax=-100;						//�o�͗��q�̍ő�c���W
	
		//AVS�o�͗��q��NUM�v�Z
			
		NUM=fluid_number;//���̗��q�̂ݏo��
			
	
			
	
		if(d==2)
		{
			for(int i=startID;i<NUM;i++)
			{
				vec<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<0<<" "<<-g[A_Y]*alp*(PART[i].T-T0)<<endl;//���͍����o��
				if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
				if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];
			}
		}
		else if(d==3)
		{
			//int d1,d2;				//�o�͂ɕK�v�Ȏ���
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
		xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
		ymax+=4*le;
		vec<<xmax<<" "<<ymax<<" "<<0<<" "<<0.1*g[A_Y]*times<<endl;//�Ō�ɖ}��o��
		vec.close();
	}
	for(int D=0;D<DIMENTION;D++) delete [] old_U[D];
	delete [] vis_T;
	delete [] mass_T;
	
}

///���x���U�v�Z�֐�
double divergence(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0)
{
    double W=0;										//���q�����x
    double R=CON->get_distancebp()*CON->get_re();	//�e�����a
    double div=0;									//���U�̒l

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

///���x���U�v�Z�֐�(WLSM�@)
double divergence2(mpsconfig *CON,vector<mpsparticle> &PART,int i,int surface_sw)
{
	//surface_sw: ON�Ȃ�\�ʂ𖳎�����AOFF�Ȃ�\�ʂ��l���ɂ��ꂽ���U���v�Z����BON,OFF�̒�`���������Ȃ��悤��
    double div=0;//���U�̒l
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N=0;					//�W���s��̌�
	int order=CON->get_divU_order();				//�ߎ��Ȗʂ̵��ް�B 1=���` 2=��
    
	//�W���s��̑傫���̌���
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

	double *matrix=new double [N*N];	//N�~N�̌W���s��
	double *matrix2=new double [N*N];	//matrix�̃o�b�N�A�b�v
	double *B1=new double [N];			//N�̉��s��
	double *B2=new double [N];			//N�̉��s��
	double *B3=new double [N];			//N�̉��s��

	for(int n=0;n<N*N;n++) matrix[n]=0;	//������
	for(int n=0;n<N;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}

	if(d==2 && order==1)				//�񎟌�
	{
		if(PART[i].N>=2)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// le�Ŋ���̂͐��K���̂��� //�@L�͗��q�̊����(�ω𑜓x�֌W)
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
				double U=(PART[j].u[A_X]-PART[i].u[A_X])/PART[i].L;
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y])/PART[i].L;
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				if(dis>1) w=1/(dis*dis*dis*dis);
				//if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);

				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//�\�ʗ��q�̑��x�𖳎�����
					
				matrix[0]+=X*X*w;			//��Xjwj
				matrix[1]+=X*Y*w;		//��XjYjwj
				matrix[3]+=Y*Y*w;			//��Yjwj
				
				B1[0]+=U*X*w;//��ujXjwj
				B1[1]+=U*Y*w;//��ujYjwj
				B2[0]+=V*X*w;//��vjXjwj
				B2[1]+=V*Y*w;//��vjYjwj
			}
			
			matrix[2]=matrix[1];		//��XjYjwj

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//�l��ۑ�

			gauss(matrix,B1,N);//�K�E�X�̏����@�ŉ���

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//�K�E�X�̏����@�ŉ���

			double a_u,b_u;
			double a_v,b_v;
			a_u=B1[0]; a_v=B2[0];//X��������
			b_u=B1[1]; b_v=B2[1];//Y��������
			
			div=(B1[0]+B2[1]);
			if(div+div!=2*div) cout<<"���x�̔��U������� i="<<i<<endl;

			/*/�덷���v�Z
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
		else div=0;	//���ӗ��q�����Ȃ�����ꍇ�̓[���ɂ���΂悢�B�ǂ������͌��z���v�Z����Ȃ�����B
	}
	if(d==2 && order==2)//�񎟌�2����
	{
		if(PART[i].N>=4)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// le�Ŋ���̂͑ł��؂�덷�h�~
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
				double U=(PART[j].u[A_X]-PART[i].u[A_X])/PART[i].L;
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y])/PART[i].L;
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				if(dis>1) w=1/(dis*dis*dis*dis);
				//if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//�\�ʗ��q�̑��x�𖳎�����
					
				matrix[0]+=X*X*w;			//��Xjwj
				matrix[1]+=X*Y*w;		//��XjYjwj
				matrix[2]+=X*X*X*w;			//��Xj^3wj
				matrix[3]+=X*X*Y*w;			//��Xj^2Yjwj
				matrix[4]+=X*Y*Y*w;			//��XjYj^2wj

				matrix[6]+=Y*Y*w;			//��Yj^2wj
				matrix[9]+=Y*Y*Y*w;			//��Yj^3wj

				matrix[12]+=X*X*X*X*w;			//��Xj^4wj
				matrix[13]+=X*X*X*Y*w;			//��Xj^3Yjwj
				matrix[14]+=X*X*Y*Y*w;			//��Xj^2Yj^2wj
	
				matrix[19]+=X*Y*Y*Y*w;			//��XjYj^3wj

				matrix[24]+=Y*Y*Y*Y*w;			//��Yj^4wj

				
				B1[0]+=U*X*w;//��ujXjwj
				B1[1]+=U*Y*w;//��ujYjwj
				B1[2]+=U*X*X*w;//��ujXj^2wj
				B1[3]+=U*X*Y*w;//��ujXjYjwj
				B1[4]+=U*Y*Y*w;//��ujYj^2wj

				B2[0]+=V*X*w;//��vjXjwj
				B2[1]+=V*Y*w;//��vjYjwj
				B2[2]+=V*X*X*w;//��vjXj^2wj
				B2[3]+=V*X*Y*w;//��vjXjYjwj
				B2[4]+=V*Y*Y*w;//��vjYj^2wj
			}
			
			matrix[5]=matrix[1];		//��XjYjwj
			matrix[7]=matrix[3];		//��Xj^2Yjwj
			matrix[8]=matrix[4];		//��XjYj^2wj
			matrix[10]=matrix[2];		//��Xj^3Yjwj
			matrix[11]=matrix[3];		//��Xj^2Yjwj
			matrix[15]=matrix[3];		//��Xj^2Yjwj
			matrix[16]=matrix[4];		//��XjYj^2wj
			matrix[17]=matrix[13];		//��Xj^3Yjwj
			matrix[18]=matrix[14];		//��Xj^2Yj^2wj
			matrix[20]=matrix[4];		//��XjYj^2wj
			matrix[21]=matrix[9];		//��Yj^3wj
			matrix[22]=matrix[14];		//��Xj^2Yj^2wj
			matrix[23]=matrix[19];		//��XjYj^3wj

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//�l��ۑ�
			gauss(matrix,B1,N);//�K�E�X�̏����@�ŉ���
			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//�K�E�X�̏����@�ŉ���

			double a_u,b_u,c_u,d_u,e_u;
			double a_v,b_v,c_v,d_v,e_v;
			a_u=B1[0]; a_v=B2[0];//X��������
			b_u=B1[1]; b_v=B2[1];//Y��������
			c_u=B1[2]; c_v=B2[2];
			d_u=B1[3]; d_v=B2[3];
			e_u=B1[4]; e_v=B2[4];
			
			div=(B1[0]+B2[1]);
			if(div+div!=2*div) cout<<"���x�̔��U������� i="<<i<<endl;
			//cout<<a_u<<" "<<b_u<<" "<<a_v<<" "<<b_v<<endl;

			/*/�덷���v�Z
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
			cout<<"���ӗ��q����2<=N<4�Ȃ̂�div=0�ɂ����B1���ߎ��Ȃǂɂ��ׂ��H"<<endl;
		}
		else div=0;	//���ӗ��q�����Ȃ�����ꍇ�̓[���ɂ���΂悢�B�ǂ������͌��z���v�Z����Ȃ�����B
	}
	if(d==2 && order==3)//2����3����
	{
		if(PART[i].N>=6)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͑ł��؂�덷�h�~
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double U=(PART[j].u[A_X]-PART[i].u[A_X]);
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				//if(dis>1) w=r*r*r*r/(dis*dis*dis*dis);
				if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//�\�ʗ��q�̑��x�𖳎�����
					
				matrix[0]+=X*X*w;			//��Xjwj
				matrix[1]+=X*Y*w;		//��XjYjwj
				matrix[2]+=X*X*X*w;			//��Xj^3wj
				matrix[3]+=X*X*Y*w;			//��Xj^2Yjwj
				matrix[4]+=X*Y*Y*w;			//��XjYj^2wj
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

				
				B1[0]+=U*X*w;//��ujXjwj
				B1[1]+=U*Y*w;//��ujYjwj
				B1[2]+=U*X*X*w;//��ujXj^2wj
				B1[3]+=U*X*Y*w;//��ujXjYjwj
				B1[4]+=U*Y*Y*w;//��ujYj^2wj
				B1[5]+=U*X*X*X*w;
				B1[6]+=U*X*X*Y*w;
				B1[7]+=U*X*Y*Y*w;
				B1[8]+=U*Y*Y*Y*w;

				B2[0]+=V*X*w;//��vjXjwj
				B2[1]+=V*Y*w;//��vjYjwj
				B2[2]+=V*X*X*w;//��vjXj^2wj
				B2[3]+=V*X*Y*w;//��vjXjYjwj
				B2[4]+=V*Y*Y*w;//��vjYj^2wj
				B2[5]+=V*X*X*X*w;
				B2[6]+=V*X*X*Y*w;
				B2[7]+=V*X*Y*Y*w;
				B2[8]+=V*Y*Y*Y*w;
			}
			
			matrix[9]=matrix[1];		//��XjYjwj
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

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//�l��ۑ�

			gauss(matrix,B1,N);//�K�E�X�̏����@�ŉ���

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//�K�E�X�̏����@�ŉ���

			double a_u,b_u,c_u,d_u,e_u,f_u,g_u,h_u,i_u;
			double a_v,b_v,c_v,d_v,e_v,f_v,g_v,h_v,i_v;
			a_u=B1[0]; a_v=B2[0];//X��������
			b_u=B1[1]; b_v=B2[1];//Y��������
			c_u=B1[2]; c_v=B2[2];//X����2�K����
			d_u=B1[3]; d_v=B2[3];//XY��������
			e_u=B1[4]; e_v=B2[4];//Y����2�K����
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
			cout<<"���ӗ��q����2<=N<6�Ȃ̂�div=0�ɂ����B1���ߎ��Ȃǂɂ��ׂ��H"<<endl;
		}
		else div=0;	//���ӗ��q�����Ȃ�����ꍇ�̓[���ɂ���΂悢�B�ǂ������͌��z���v�Z����Ȃ�����B
	}
	else if(d==3 && order==1)//3����1����
	{
		if(PART[i].N>5)
		{
			//div=cacl_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//�Ō�̈����͖��m��
			div=cacl_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//�Ō�̈����͖��m��
		}
		else 
		{
			div=0;
			//cout<<"�x�� �ߗח��q����5�ȉ�("<<PART[i].N<<")�ł�"<<endl;
		}
	}
	else if(d==3 && order==2)//3����2����
	{
		//P=Pi+a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+i��z��x�Ƃ����ƁA
		///�W���s���
		///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z     a = ����x��P  
		///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z   b = ����y��P
		///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2     c = ����z��P
		///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z     d = ����x2��P
		///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z  e = ����y2��P
		///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3     f = ����z2��P
		///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z  g = ����x��y��P
		///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2  h = ����y��z��P
		///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2    g = ����x��z��P
		
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
		if(nei>=6)			//le�ȉ��̗��q���U�ȏ゠��΂Q���ߎ�
		{
			//div=cacl_WLSM_divu_D3_order2(CON,PART,matrix,B1,B2,B3,i,9);
			div=cacl_WLSM_divu_D3_order2_2(CON,PART,matrix,B1,B2,B3,i,10);
		}
		else if(PART[i].N>5)
		{
			//div=cacl_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//���̏ꍇ�A�Ō�̈�����N=3��n�����Ƃɒ���
			div=cacl_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//�Ō�̈����͖��m��
		}
		else 
		{
			div=0;
			//cout<<"�x�� �ߗח��q����9�ȉ�("<<PART[i].N<<")�ł�"<<endl;
		}
	}

	delete [] matrix;
	delete [] matrix2;
	delete [] B1;
	delete [] B2;
	delete [] B3;

    return div;
}

///���x���U�v�Z�֐�(WLSM�@) �ߎ��Ȗʂ����g�̊֐��l��ʂ�Ƃ�����������Ȃ��ꍇ
double divergence3(mpsconfig *CON,vector<mpsparticle> &PART,int i,int surface_sw)
{
	//surface_sw: ON�Ȃ�\�ʂ𖳎�����AOFF�Ȃ�\�ʂ��l���ɂ��ꂽ���U���v�Z����BON,OFF�̒�`���������Ȃ��悤��
    double div=0;//���U�̒l
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N0=0;							//�W���s��̌�
	int order=CON->get_divU_order();				//�ߎ��Ȗʂ̵��ް�B 1=���` 2=��
    
	//�W���s��̑傫���̌���
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

	double *matrix=new double [N0*N0];	//N�~N�̌W���s��
	double *matrix2=new double [N0*N0];	//matrix�̃o�b�N�A�b�v
	double **MAT= new double*[N0];		//N�~N�̌W���s��(�z���2����)
	for(int n=0;n<N0;n++) MAT[n]=new double[N0];
	double *base=new double[N0];			//���x�N�g���i�[
	double *B1=new double [N0];			//N�̉��s��
	double *B2=new double [N0];			//N�̉��s��
	double *B3=new double [N0];			//N�̉��s��

	for(int n=0;n<N0*N0;n++) matrix[n]=0;	//������
	for(int n=0;n<N0;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}
	for(int n=0;n<N0;n++) for(int m=0;m<N0;m++)  MAT[n][m]=0;//������

	int N=N0;
	int jnb=PART[i].N;	//���ӗ��q��

	if(surface_sw==ON)	//�\�ʗ��q�𖳎�����Ȃ�(�\�ʗ��q�̑��x�͑��x���U�[���𖞂����Ă��Ȃ�����)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			if(PART[j].type==FLUID && PART[j].surface==ON) jnb--; 
		}
	}

	if(d==2 )				//�񎟌�
	{
		if(jnb<15) {order=2; N=6;}		//���ӗ��q����3���ߎ�����̂ɏ\���łȂ��Ȃ�A2���ߎ�����(�ی��������ď����傫�߂ɂƂ��Ă���)
		if(jnb<6) {order=1; N=3;}		//���ӗ��q�����Q���ߎ�����̂ɏ\���łȂ��Ȃ�A�P���ߎ�����
		if(CON->get_divU_order()==1) {order=1; N=3;}
		if(jnb>=2)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X])/PART[i].L;// le�Ŋ���̂͐��K���̂���
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/PART[i].L;
				double U=(PART[j].u[A_X])/PART[i].L;
				double V=(PART[j].u[A_Y])/PART[i].L;
				double dis=sqrt(X*X+Y*Y);
					
				double w=1;
				if(dis>1) w=1/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//�\�ʗ��q�̑��x�𖳎�����
				
				if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//���x�N�g��
				else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=1;}
				else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y; base[9]=1;}

				//�s��쐬
				for(int n=0;n<N;n++)
				{
					for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
					B1[n]+=base[n]*w*U;
					B2[n]+=base[n]*w*V;
				}
			}
			
			MAT[N-1][N-1]+=1;		//��ԉE���̔z��Ɏ������g�̊�^��������
			B1[N-1]+=PART[i].u[A_X]/PART[i].L;	//���qi�̊�^�B
			B2[N-1]+=PART[i].u[A_Y]/PART[i].L;	//���qi�̊�^�B
			//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�*/

			/*double *eigen=new double[N];
			calc_eigen_by_jacobi(MAT, N,eigen);	//�ŗL�l�����߂�
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
					
			for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//�l��matrix�ɓ]��

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//�l��ۑ�

			gauss(matrix,B1,N);//�K�E�X�̏����@�ŉ���

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//�K�E�X�̏����@�ŉ���

			div=(B1[0]+B2[1]);

			//cout<<div<<" "<<N<<" "<<jnb<<" "<<max/min<<endl;

			if(div+div!=2*div) cout<<"���x�̔��U������� i="<<i<<endl;
			else
			{
				/*double Q1=0; double Q2=0;		//U�̌덷
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
					if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//�\�ʗ��q�̑��x�𖳎�����
					if(order==1) {base[0]=X; base[1]=Y; base[2]=1;}	//���x�N�g��
					else if(order==2) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=1;}
					else if(order==3) {base[0]=X; base[1]=Y; base[2]=X*X; base[3]=X*Y; base[4]=Y*Y; base[5]=X*X*X; base[6]=X*X*Y; base[7]=X*Y*Y; base[8]=Y*Y*Y; base[9]=1;}
					
					double uu=0; double vv=0;			//�ߎ��Ȗʏ�̑��x
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
		else div=0;	//���ӗ��q�����Ȃ�����ꍇ�̓[���ɂ���΂悢�B�ǂ������͌��z���v�Z����Ȃ�����B
	}
	else if(d==3)
	{
		///�W���s���1���ߎ��Ȃ�
		///   ����x2    ����x��y  ����x��z ����x a = ����xfj  
		///  ����x��y    ����y2   ����y��z ����y b = ����yfj 
		///  ����x��z   ����y��z  ����z2  ����z  c = ����zfj 
		///  ����x      ����y     ����z     ��1  d = ��fj

		//2���ߎ��Ȃ�
		//P=a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+k��z��x+l�Ƃ����ƁA
		///�W���s���
		///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z    ����x    a = ����xf  
		///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z  ����y    b = ����yf
		///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2    ����z    c = ����zf
		///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z    ����x2   d = ����x2f
		///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z ����y2   e = ����y2f
		///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3	 ����z2   f = ����z2f
		///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z ����x��y g = ����x��yf
		///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2 ����y��z h = ����y��zf
		///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2   ����z��x k = ����x��zf
		///   ����x       ����y ����z ����x2      ����y2       ����z        ����x��y     ����y��z      ����x��z     ����z��x     ��1      l = ��fj

		//3���ߎ��Ȃ�
		//P=a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+k��z��x+l��x3+m��y3+n��z3+o��x2��y+p��x��y2+q��y2��z+r��y��z2+s��x2��z+t��z2��x+u��x��y��z+v�Ƃ���
		
		if(order==3) 
		{
			if(jnb<20) {order=2; N=10;}//���q�������Ȃ��Ƃ���2���ߎ��ɂ���
		}
		if(order==2) 
		{
			if(jnb<10) {order=1; N=4;}//����ł����q�������Ȃ��Ƃ��͂P���ߎ��ɂ���
		}
		
		
		if(jnb>=3)	//���ӗ��q���Q�ȉ���������v�Z���Ȃ�
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͑ł��؂�덷�h�~
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
				double U=(PART[j].u[A_X]);
				double V=(PART[j].u[A_Y]);
				double W=(PART[j].u[A_Z]);
				double dis=sqrt(X*X+Y*Y+Z*Z);
					
				double w=1;
				if(dis>PART[i].L) w=PART[i].L*PART[i].L*PART[i].L*PART[i].L/(dis*dis*dis*dis);
				if(surface_sw==ON) if(PART[j].type==FLUID && PART[j].surface==ON) w=0;	//�\�ʗ��q�̑��x�𖳎�����

				if(order==1) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=1;}	//���x�N�g��
				else if(order==2) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=1;}
				else if(order==3) {base[0]=X; base[1]=Y; base[2]=Z; base[3]=X*X; base[4]=Y*Y; base[5]=Z*Z; base[6]=X*Y; base[7]=Z*Y; base[8]=X*Z; base[9]=X*X*X; base[10]=Y*Y*Y; base[11]=Z*Z*Z; base[12]=X*X*Y; base[13]=X*Y*Y; base[14]=Y*Y*Z; base[15]=Y*Z*Z; base[16]=X*X*Z; base[17]=X*Z*Z; base[18]=X*Y*Z; base[19]=1;}
				
				//�s��쐬
				for(int n=0;n<N;n++)
				{
					for(int m=0;m<N;m++) MAT[n][m]+=base[n]*base[m]*w;
					B1[n]+=base[n]*w*U;
					B2[n]+=base[n]*w*V;
					B3[n]+=base[n]*w*W;
				}
			}
			
			MAT[N-1][N-1]+=1;		//��ԉE���̔z��Ɏ������g�̊�^��������
			B1[N-1]+=PART[i].u[A_X];	//���qi�̊�^�B
			B2[N-1]+=PART[i].u[A_Y];	//���qi�̊�^�B
			B3[N-1]+=PART[i].u[A_Z];	//���qi�̊�^�B
			//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�
					
			for(int n=0;n<N*N;n++) matrix[n]=MAT[n/N][n%N];	//�l��matrix�ɓ]��

			for(int n=0;n<N*N;n++) matrix2[n]=matrix[n];//�l��ۑ�

			gauss(matrix,B1,N);//�K�E�X�̏����@�ŉ���

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B2,N);//�K�E�X�̏����@�ŉ���

			for(int n=0;n<N*N;n++) matrix[n]=matrix2[n];
			gauss(matrix,B3,N);//�K�E�X�̏����@�ŉ���

			
			double dudx=B1[0];	//1���ߎ��ł�2���ߎ��ł��A���m���̏��ԓI�ɂ����Ȃ�悤�ɂ��Ă���
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

///���x���U�v�Z�֐�(������̕��@)�@�e���\���ςƏd�ݕ��ς�p����
double divergence4(mpsconfig *CON,vector<mpsparticle> &PART,int i)
{
	//2���������ł��ĂȂ�
    double div=0;//���U�̒l
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimention();
	int N=0;					//�W���s��̌�
	int order=1;				//�ߎ��Ȗʂ̵��ް�B 1=���` 2=��
    
	//�W���s��̑傫���̌���
	if(order==1) N=d;
	//else //�܂��ł��ĂȂ�
	////////////////////////////////

	double *matrix=new double [N*N];	//N�~N�̌W���s��
	double *B1=new double [N];			//N�̉��s��
	double *B2=new double [N];			//N�̉��s��
	double *B3=new double [N];			//N�̉��s��

	for(int n=0;n<N*N;n++) matrix[n]=0;	//������
	for(int n=0;n<N;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}

	if(d==2 && order==1)				//�񎟌�
	{
		if(PART[i].N>=2)
		{
			double W=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
			
				double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͐��K���̂���
				double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
				double U=(PART[j].u[A_X]-PART[i].u[A_X]);
				double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
				double dis=sqrt(X*X+Y*Y);
				X/=dis;	Y/=dis;				//�P�ʃx�N�g���ɂ���
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

			//matrix�̋t�s������߂� ���܂����t�s���matrix�̒����㏑�����Ċi�[�����
			calc_inverse_matrix(CON,PART, N, matrix);

			double dudx=matrix[0]*B1[0]+matrix[1]*B1[1];
			double dudy=matrix[2]*B1[0]+matrix[3]*B1[1];
			double dvdx=matrix[0]*B2[0]+matrix[1]*B2[1];
			double dvdy=matrix[2]*B2[0]+matrix[3]*B2[1];
			
			div=(dudx+dvdy);
			if(div+div!=2*div) cout<<"���x�̔��U������� i="<<i<<endl;

		}
		else div=0;	//���ӗ��q�����Ȃ�����ꍇ�̓[���ɂ���΂悢�B�ǂ������͌��z���v�Z����Ȃ�����B
	}
	else cout<<"�܂��ł��ĂȂ�"<<endl;
	

	delete [] matrix;
	delete [] B1;
	delete [] B2;
	delete [] B3;

    return div;
}

//5���̘A���������̉�1,2��Ԃ��֐�
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

//divergence2�ɂ�����A3����1���ߎ����s���֐�
double cacl_WLSM_divu_D3_order1(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///�W���s���
	///   ����x2    ����x��y  ����x��z  a = ����x��f  
	///  ����x��y    ����y2   ����y��z  b = ����y��f 
	///  ����x��z   ����y��z   ����z2   c = ����z��f 

	double le=CON->get_distancebp();
	double matrix_val[9];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͑ł��؂�덷�h�~
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
					
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
		weight[k]=w;
					
		matrix[0]+=X*X*w;			//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
					
		matrix[4]+=Y*Y*w;			//��Yj^2wj
		matrix[5]+=Y*Z*w;		//��YjZjwj

		matrix[8]+=Z*Z*w;			//��Zj^2wj
			
		B1[0]+=U*X*w;//��fjXjwj
		B1[1]+=U*Y*w;//��fjYjwj
		B1[2]+=U*Z*w;//��fjZjwj

		B2[0]+=V*X*w;//��fjXjwj
		B2[1]+=V*Y*w;//��fjYjwj
		B2[2]+=V*Z*w;//��fjZjwj

		B3[0]+=W*X*w;//��fjXjwj
		B3[1]+=W*Y*w;//��fjYjwj
		B3[2]+=W*Z*w;//��fjZjwj
	}
			
	matrix[3]=matrix[1];		//��XjYjwj
	matrix[6]=matrix[2];		//��XjZjwj
	matrix[7]=matrix[5];		//��YjZjwj

	for(int L=0;L<9;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�

	/*double dudx=0;//�������̂ق����኱�����B���ǌ덷�]���������Ȃ�K�E�X
	double dvdy=0;
	double dwdz=0;
	double determinant=(matrix[0]*matrix[4]*matrix[8]-matrix[0]*matrix[5]*matrix[7]-matrix[1]*matrix[3]*matrix[8]+matrix[1]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*matrix[7]-matrix[2]*matrix[4]*matrix[6]);//�s��
			
	dudx=(B1[0]*matrix[4]*matrix[8]-B1[0]*matrix[5]*matrix[7]-matrix[1]*B1[1]*matrix[8]+matrix[1]*matrix[5]*B1[2]+matrix[2]*B1[1]*matrix[7]-matrix[2]*matrix[4]*B1[2])/determinant;
	dvdy=(matrix[0]*B2[1]*matrix[8]-matrix[0]*matrix[5]*B2[2]-B2[0]*matrix[3]*matrix[8]+B2[0]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*B2[2]-matrix[2]*B2[1]*matrix[6])/determinant;
	dwdz=(matrix[0]*matrix[4]*B3[2]-matrix[0]*B3[1]*matrix[7]-matrix[1]*matrix[3]*B3[2]+matrix[1]*B3[1]*matrix[6]+B3[0]*matrix[3]*matrix[7]-B3[0]*matrix[4]*matrix[6])/determinant;
	*/	

	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<3;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else B1[0]=0; //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else B2[1]=0;	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else B3[2]=0;	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	//�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̌덷
	double W=0;//�d�݂̑��a
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

//divergence2�ɂ�����A3����1���ߎ����s���֐�ver.2
double cacl_WLSM_divu_D3_order1_2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///�W���s���
	///   ����x2    ����x��y  ����x��z ����x a = ����xfj  
	///  ����x��y    ����y2   ����y��z ����y b = ����yfj 
	///  ����x��z   ����y��z  ����z2  ����z  c = ����zfj 
	///  ����x      ����y     ����z     ��1  d = ��fj

	double le=CON->get_distancebp();
	double matrix_val[16];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// le�Ŋ���̂͑ł��؂�덷�h�~
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
					
		matrix[0]+=X*X*w;			//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
		matrix[3]+=X*w;
					
		matrix[5]+=Y*Y*w;			//��Yj^2wj
		matrix[6]+=Y*Z*w;		//��YjZjwj
		matrix[7]+=Y*w;

		matrix[10]+=Z*Z*w;			//��Zj^2wj
		matrix[11]+=Z*w;

		matrix[15]+=w;
			
		B1[0]+=U*X*w;//��fjXjwj
		B1[1]+=U*Y*w;//��fjYjwj
		B1[2]+=U*Z*w;//��fjZjwj
		B1[3]+=U*w;//��fjwj

		B2[0]+=V*X*w;//��fjXjwj
		B2[1]+=V*Y*w;//��fjYjwj
		B2[2]+=V*Z*w;//��fjZjwj
		B2[3]+=V*w;//��fjwj

		B3[0]+=W*X*w;//��fjXjwj
		B3[1]+=W*Y*w;//��fjYjwj
		B3[2]+=W*Z*w;//��fjZjwj
		B3[3]+=W*w;//��fjwj
	}
			
	matrix[4]=matrix[1];	
	matrix[8]=matrix[2];		
	matrix[9]=matrix[6];
	matrix[12]=matrix[3];
	matrix[13]=matrix[7];
	matrix[14]=matrix[11];

	matrix[15]+=1;//�������g
	B1[3]+=PART[i].u[A_X];
	B2[3]+=PART[i].u[A_Y];
	B3[3]+=PART[i].u[A_Z];

	for(int L=0;L<16;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�

	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<4;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B2[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B3[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else for(int k=0;k<4;k++) B1[k]=0; //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrix��ύX�O�ɖ߂�
	}
	else for(int k=0;k<4;k++) B2[k]=0;	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else for(int k=0;k<4;k++) B3[k]=0;	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	//�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̌덷
	double err[3];
	double W=1;//�d�݂̑��a
	err[A_X]=B1[3]-PART[i].u[A_X];//���g�̌덷
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

//divergence2�ɂ�����A3����2���ߎ����s���֐�
double cacl_WLSM_divu_D3_order2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	double le=CON->get_distancebp();
	double matrix_val[81];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

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
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
		weight[k]=w;

		matrix[0]+=X*X*w;		//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
		matrix[3]+=X*X*X*w;		//��Xj^3wj
		matrix[4]+=X*Y*Y*w;		//��XjYj^2wj
		matrix[5]+=X*Z*Z*w;		//��XjZj^2wj
		matrix[6]+=X*X*Y*w;		//��Xj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//��XjYjZjwj
		matrix[8]+=X*X*Z*w;		//��Xj^2Zjwj
	
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

		matrix[50]+=Z*Z*Z*Z*w;	//6�s��
		matrix[51]+=X*Y*Z*Z*w;
		matrix[52]+=Y*Z*Z*Z*w;
		matrix[53]+=X*Z*Z*Z*w;

		//7�`9�s�ڂ͂��ׂĊ����̗v�f����]�p���\


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
	matrix[9]=matrix[1];		//��XjYjwj
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

	for(int L=0;L<81;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�
			
	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;
	int Yflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Zflag=OFF;
	for(int k=0;k<9;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	////�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̕W���΍�
	double W=0;		//�d�݂̑��a
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

//divergence2�ɂ�����A3����2���ߎ����s���֐�ver.2 ���m�����ЂƂ���
double cacl_WLSM_divu_D3_order2_2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	//P=a��x+b��y+c��z+d��x2+e��y2+f��z2+g��x��y+h��y��z+i��z��x+P�Ƃ����ƁA
	///�W���s���
	///   ����x2      ����x��y    ����x��z    ����x3       ����x��y2    ����x��z2    ����x2��y     ����x��y��z  ����x2��z    ����x    a = ����xf  
	///   ����x��y    ����y2      ����y��z    ����x2��y    ����y3       ����y��z2    ����x��y2     ����y2��z    ����x��y��z  ����y    b = ����yf
	///   ����x��z    ����y��z    ����z2      ����x2��z    ����y2��z    ����z3       ����x��y��z   ����y��z2    ����x��z2    ����z    c = ����zf
	///   ����x3      ����x2��y   ����x2��z   ����x4       ����x2��y2   ����x2��z2   ����x3��y     ����x2��y��z ����x3��z    ����x2   d = ����x2f
	///   ����x��y2   ����y3      ����y2��z   ����x2��y2   ����y4       ����y2��z2   ����x��y3     ����y3��z    ����x��y2��z ����y2   e = ����y2f
	///   ����x��z2   ����y��z2   ����z3      ����x2��z2   ����y2��z2   ����z4       ����x��y��z2  ����y��z3    ����x��z3	 ����z2   f = ����z2f
	///   ����x2��y   ����x��y2   ����x��y��z ����x3��y    ����x��y3    ����x��y��z2 ����x2��y2    ����x��y2��z ����x2��y��z ����x��y g = ����x��yf
	///   ����x��y��z ����y2��z   ����y��z2   ����x2��y��z ����y3��z    ����y��z3    ����x��y2��z  ����y2��z2   ����x��y��z2 ����y��z h = ����y��zf
	///   ����x2��z   ����x��y��z ����x��z2   ����x3��z    ����x��y2��z  ����x��z3   ����x2��y ��z ����x��y��z2 ����x2��z2   ����z��x i = ����x��zf
	///   ����x       ����y ����z ����x2      ����y2       ����z        ����x��y     ����y��z      ����x��z     ����z��x     ��1      P = ��fj

	double le=CON->get_distancebp();
	double matrix_val[100];			//matrix�̒l�͈�xgauss()�Ŏg�p����ƒl���ς��̂ŁA�ύX�O��matrix��ۑ�����
	double *weight=new double [PART[i].N];	//���ӗ��q�̏d�݂��i�[����B���ƂŌ덷�]���̂Ƃ��Čv�Z���s�v�ɂȂ�B

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
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//���̎�����dis=R�̂Ƃ��d�݂��[���ɋ߂��Ȃ�
		weight[k]=w;

		matrix[0]+=X*X*w;		//��Xj^2wj
		matrix[1]+=X*Y*w;		//��XjYjwj
		matrix[2]+=X*Z*w;		//��XjZjwj
		matrix[3]+=X*X*X*w;		//��Xj^3wj
		matrix[4]+=X*Y*Y*w;		//��XjYj^2wj
		matrix[5]+=X*Z*Z*w;		//��XjZj^2wj
		matrix[6]+=X*X*Y*w;		//��Xj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//��XjYjZjwj
		matrix[8]+=X*X*Z*w;		//��Xj^2Zjwj
		matrix[9]+=X*w;		//��Xj^2Zjwj
	
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

		matrix[55]+=Z*Z*Z*Z*w;	//6�s��
		matrix[56]+=X*Y*Z*Z*w;
		matrix[57]+=Y*Z*Z*Z*w;
		matrix[58]+=X*Z*Z*Z*w;

		matrix[99]+=w;
		//7�`9�s�ڂ͂��ׂĊ����̗v�f����]�p���\

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

	for(int k=0;k<=2;k++) matrix[30+k]=matrix[3+10*k];//30�`32�v�f
	matrix[39]=matrix[0];

	for(int k=0;k<=3;k++) matrix[40+k]=matrix[4+10*k];//40�`43�v�f
	matrix[49]=matrix[11];

	for(int k=0;k<=4;k++) matrix[50+k]=matrix[5+10*k];//50�`54�v�f
	matrix[59]=matrix[22];

	for(int k=0;k<=5;k++) matrix[60+k]=matrix[6+10*k];//60�`65�v�f
	matrix[66]=matrix[34];
	matrix[67]=matrix[48];
	matrix[68]=matrix[37];
	matrix[69]=matrix[1];

	for(int k=0;k<=6;k++) matrix[70+k]=matrix[7+10*k];//70�`76�v�f
	matrix[77]=matrix[54];
	matrix[78]=matrix[56];
	matrix[79]=matrix[12];
	
	for(int k=0;k<=7;k++) matrix[80+k]=matrix[8+10*k];//80�`87�v�f
	matrix[88]=matrix[35];
	matrix[89]=matrix[20];

	for(int k=0;k<=8;k++) matrix[90+k]=matrix[9+10*k];//90�`98�v�f

	matrix[99]+=1;//���g
	B1[9]+=PART[i].u[A_X];
	B2[9]+=PART[i].u[A_Y];
	B3[9]+=PART[i].u[A_Z];

	for(int L=0;L<100;L++) matrix_val[L]=matrix[L];//�s��̒l��ۑ�
			
	//�s����K�E�X�̏����@�ŉ����@����B�Ɋi�[�����
	int Xflag=OFF;
	int Yflag=OFF;//�K�E�X�̏����@�����邩���Ȃ����B���s�񂪃[���Ȃ�K�E�X�̏����@�����炾�߁B
	int Zflag=OFF;
	for(int k=0;k<10;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B2[k]!=0) Yflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
		if(B3[k]!=0) Zflag=ON;//B1[k]1�����ׂă[���Ȃ�Xflag��OFF�̂܂܁B
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flag��OFF�Ȃ�ǂ݂̂�dudx�̓[���B

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//�s��̒l�߂�
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dvdy�̓[���B

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flag��OFF�Ȃ�ǂ݂̂�dwdz�̓[���B

	////�v�Z�I��

	//�덷�𒲍�
	double Q[3]={0,0,0};//�e�����̕W���΍�
	double W=0;		//�d�݂̑��a
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

///quickMPS�p�|�X�g�����֐�
void post_processing3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double TIME)
{
	
	///restart�p̧�ُo��
	if(t==CON->get_step() || t%CON->get_autosave()==0)
	{
		////////restart�p�ɗ��q���Ɨ��q�f�[�^���L�^
		ofstream hoge1("number.dat");
		hoge1<<particle_number<<endl;
		hoge1<<TIME<<endl;
		hoge1.close();
		
		FILE *hoge6;
		hoge6=fopen("initial_input.dat","w");//mps_input.dat�Ƃ͋�ʂ��Ȃ��ƁArestart�����s�����Ƃ�����
		for(int i=0;i<particle_number;i++)///���̉�͗��q�����ɋL�q
		{
			fprintf( hoge6, "%d\t",i);
			fprintf( hoge6, "%5.15f\t",PART[i].r[A_X]);
			fprintf( hoge6, "%5.15f\t",PART[i].r[A_Y]);
			fprintf( hoge6, "%5.15f\t",PART[i].r[A_Z]);
			fprintf( hoge6, "%5.15f\t",PART[i].u[A_X]);  //���xx����
			fprintf( hoge6, "%5.15f\t",PART[i].u[A_Y]);  //���xy����
			fprintf( hoge6, "%5.15f\t",PART[i].u[A_Z]);  //���xz����
			fprintf( hoge6, "%5.15f\t",PART[i].P); //����
			fprintf( hoge6, "%5.15f\t",PART[i].h); //�G���^���s�[
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

//���q����͗̈�̊O�ɂłĂ��Ȃ����`�F�b�N
int check_position(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int *particle_number)
{
	int sw=OFF;	//�{�֐��ŕԂ��l�BOFF�Ȃ痱�q���ɕω��Ȃ��BON�Ȃ�ω������Ƃ�����
	int num=0;	//���ł��闱�q��
	double le=CON->get_distancebp();
	int *flag=new int [fluid_number];	//ON�Ȃ�̈�� OFF�Ȃ�̈�O

	double Xmax=CON->get_maxX(); double Xmin=CON->get_minX();
	double Ymax=CON->get_maxY(); double Ymin=CON->get_minY();
	double Zmax=CON->get_maxZ(); double Zmin=CON->get_minZ();

	double dx=CON->get_dx()*le;	//�i�q��

	Xmax-=dx; Ymax-=dx; Zmax-=2*dx;		//�ی���������1�i�q�������ɋ��E���Ƃ�B������O���Ȃ痱�q������
	Xmin+=dx; Ymin+=dx; Zmin+=2*dx;

	vector<mpsparticle>::iterator p,p0;//�����q
	p0=PART.begin();

	for(int i=0;i<fluid_number;i++) 
	{
		flag[i]=ON;
		if(PART[i].r[A_X]<Xmin || PART[i].r[A_X]>Xmax) flag[i]=OFF;
		else if(PART[i].r[A_Y]<Ymin || PART[i].r[A_Y]>Ymax) flag[i]=OFF;
		else if(PART[i].r[A_Z]<Zmin || PART[i].r[A_Z]>Zmax) flag[i]=OFF;
	}//flag[i]�����܂���

	int min_nei=3;//5;//CON->get_min_nei();
	if(CON->get_dimention()==3) min_nei=0;//5;
	for(int i=0;i<fluid_number;i++) if(PART[i].N<1) flag[i]=OFF;//���ӗ��q�������Ȃ����q���폜?

	//���x�ɏ���l��݂��A�����Ă�����폜
	double limit_U=CON->get_max_speed();			//���x����������闱�q�͍폜����
	for(int i=0;i<fluid_number;i++)
	{
		double speed=0;
		for(int D=0;D<3;D++) speed+=PART[i].u[D]*PART[i].u[D];
		speed=sqrt(speed);
		if(speed>limit_U) flag[i]=OFF;
	}////*/

	//�ߐڂ������Ă���ꍇ�A�ԍ��̎Ⴂ�ق����폜
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
				if(PART[j].type!=FLUID) flag[i]=OFF;//���qj���Ǘ��q�Ȃ�A���̑��������B
				else	//���̓��m�̏ꍇ
				{
					if(PART[i].surface==ON)
					{
						if(PART[j].surface==ON) flag[i]=OFF;
						else flag[j]=OFF;	//����i���\�ʂŁA����j�������Ȃ�Aj������
					}
					else flag[i]=OFF;
				}
			}
		}
	}////////*/


	//for(int i=0;i<fluid_number;i++) if(PART[i].val!=0) flag[i]=ON;

	for(int i=0;i<fluid_number;i++) if(flag[i]==OFF) num++;
	/*///
	//�������闱�q�̓d�ׂ����ӗ��q�֕��z (val���d�ז��x�ł����Ȃ��ꍇ�́A���̂܂܏���������΂悢)
	if(CON->get_conductive_appro()==ON && CON->get_EM_calc_type()==1)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(flag[i]==OFF)
			{
				int jnb=0;									//���qi�̎��ӗ��q�̂����Aflag=ON�ȗ��q��
				for(int k=0;k<PART[i].N;k++)  if(flag[PART[i].NEI[k]]==ON) jnb++;
				if(jnb>0)
				{
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(flag[j]==ON) PART[j].val+=PART[i].val/jnb;			//���̋ߎ�����ꍇ�́A�������闱�q�̓d�ׂ��ߗח��q�ɕ��z���Ă����B
					}
					PART[i].val=0;
				}
			}
		}
	}
	/////*/

	
	if(num>0)//�̈�O���q�����m�����Ȃ�
	{
		sw=ON;
		int *erase_id=new int[num];
		int count=0;
		double val=0;
		for(int i=0;i<fluid_number;i++)
		{
			if(flag[i]==OFF)
			{
				erase_id[count]=i;//�����ׂ�id���L��
				val+=PART[i].val;
				count++;
			}
		}
		if(val!=0) cout<<"val="<<val<<"������"<<endl;
		
		for(int i=0;i<num;i++)
		{
			p=PART.begin();
			p+=erase_id[i];
			
			PART.erase(p);
			for(int j=i+1;j<num;j++)
			{
				erase_id[j]=erase_id[j]-1;//1�l��������
			}
		}
		delete [] erase_id;

		//id�����ꂽ����߂�
		for(int i=0;i<PART.size();i++) if(PART[i].id!=i) PART[i].id=i; 
	}

	if(sw==ON)
	{
		cout<<"�̈�O���q��T�m "<<num<<"�̗��q������--���݂̗��q��="<<PART.size()<<endl;
		*particle_number=*particle_number-num;
	}
	
	delete [] flag;

	return sw;
}

//���q���߂Â�������̂�h���֐�
void modify_position(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double dt,int particle_number)
{
	double le=CON->get_distancebp();

	double *as_i=new double[fluid_number];		//���qi�Ƃ��āA���̊֐��ɎQ�����邩�A���Ȃ���
	double *as_j=new double[particle_number];		//���qi�Ƃ��āA���̊֐��ɎQ�����邩�A���Ȃ���
	double *new_r[DIMENTION];
    for(int D=0;D<DIMENTION;D++) new_r[D]=new double [fluid_number];//�V�����ʒu�x�N�g��

	if(CON->get_modify_position()==1)
	{
		for(int i=0;i<fluid_number;i++) as_i[i]=ON;		//���ׂĂ̗��q�Ōv�Z
		for(int i=0;i<particle_number;i++) as_j[i]=ON;
	}
	else if(CON->get_modify_position()==2)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON) 
			{
				as_i[i]=ON;		//�\�ʗ��q�Ōv�Z
				as_j[i]=ON;		//�\�ʗ��q�Ōv�Z
			}
			else as_i[i]=OFF;
		}
		for(int i=fluid_number;i<particle_number;i++) as_j[i]=ON;//�ǂ͂��ׂčl��
	}

	

	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) new_r[D][i]=PART[i].r[D];	//������
	for(int i=0;i<fluid_number;i++)
	{
		if(as_i[i]==ON)
		{
			double mindis=le;
			mindis=100;			//�����̒Z�����̂��������邾���łȂ��A�������̂�le�ɂ������Ƃ��͂�����
			int J=i;			//�ŋߐڗ��q
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
				double L=le-mindis;//�J���ׂ�����
				double dL[DIMENTION];
				for(int D=0;D<DIMENTION;D++) dL[D]=PART[J].r[D]-PART[i].r[D];
				if(J!=i && PART[J].type==FLUID)//le���ߐڂ��Ă��闬�̗��q���������Ȃ�
				{
					for(int D=0;D<DIMENTION;D++)
					{
						double dU=0.5*L/dt;	//�ω����ׂ����x
					//	PART[J].r[D]+=dL[D]/mindis*dU*dt;
						new_r[D][J]+=dL[D]/mindis*dU*dt;

					//	PART[i].r[D]-=dL[D]/mindis*dU*dt;
						new_r[D][i]-=dL[D]/mindis*dU*dt;
					}				
				}
				else if(J!=i && PART[J].type!=FLUID)//le���ߐڂ��Ă���Ǘ��q���������Ȃ�
				{
					for(int D=0;D<DIMENTION;D++)
					{
						double dU=L/dt;	//�ω����ׂ����x
					//	PART[i].r[D]-=dL[D]/mindis*dU*dt;
						new_r[D][i]-=dL[D]/mindis*dU*dt;
					}				
				}
			}
		}
	}
	for(int i=0;i<fluid_number;i++) for(int D=0;D<DIMENTION;D++) PART[i].r[D]=new_r[D][i];	//���

	delete [] as_i;
	delete [] as_j;
	for(int D=0;D<DIMENTION;D++) delete [] new_r[D];
}

//�e���q�̓��S���W���v�Z�֐�
void calc_vis_value(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *vis,double dt,int t,int particle_number)
{
	double R[3];
	double le=CON->get_distancebp();
	double r=CON->get_re()*le;
	int d=CON->get_dimention();
	double RR=8.314;						//�K�X�萔8.314[J/mol/K]
	double TT=CON->get_roomT();				//����[K]
	double *Q=new double[fluid_number];		//���� [J/mol] A6061:145000, A1050:156888,
	double *alpha=new double[fluid_number];	//[1/Pa]		A6061:0.045*1e-6, A1050:0.037*1e-6,
	double *A=new double[fluid_number];		//[1/sec]	A6061:/8.8632*1e6? ����Ƃ�exp(19.3)?, A1050:exp(26.69),
	double *N=new double[fluid_number];		//�w��			A6061:3.55, A1050:3.84,
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1)

		{
			Q[i]=158300;			//A1100�̃p�����[�^
			alpha[i]=0.045*1e-6;
			A[i]=exp(24.67);
			N[i]=5.66;
		}
		else if(PART[i].materialID==2)
		{
			Q[i]=166900;			//A5056�̃p�����[�^
			alpha[i]=0.015*1e-6;
			A[i]=exp(23.05);
			N[i]=4.82;
		}
	}

	double V=get_volume(CON);				//���q�̑̐�
	double *density=new double [fluid_number];//�e���q�̖��x�i�[
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].materialID==1) density[i]=CON->get_density();
		else if(PART[i].materialID==2) density[i]=CON->get_density2();
	}
	double *mass=new double [fluid_number];	//�e���q�̎��ʊi�[
	for(int i=0;i<fluid_number;i++) mass[i]=density[i]*V;
		
	double co=0.9;								//�d�����M�ɂȂ�(���X)����(0�`1)
	double *Cp=new double [fluid_number];		//��M
	double *MP=new double [fluid_number];		//�Z�_
	double *latent_H=new double [fluid_number];		//�Z�_

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

	double *sigma=new double [fluid_number];//�e���q�̑������͊i�[
	double *ep=new double [fluid_number];	//�e���q�̑����Ђ��ݑ��x�i�[
	double *T=new double [fluid_number];	//�e���q�̉��x�i�[[K]
	double *heat=new double [fluid_number];	//�e���q�̔��M�ʊi�[[J]

	//unsigned int timeA=GetTickCount();
	//#pragma omp parallel for
	for(int i=0;i<fluid_number;i++)
	{
		double hs0=mass[i]*Cp[i]*MP[i];//�Z���J�n�_�̃G���^���s�[
		double hs1=hs0+latent_H[i]*mass[i];    //�Z���I���_�̃G���^���s�[	

		///���x�s�̌v�Z
		if(PART[i].h<hs0) T[i]=PART[i].h/mass[i]/Cp[i];///�ő�
		else if(hs0<=PART[i].h && PART[i].h<=hs1) T[i]=MP[i];//�Z�_
		else if(hs1<PART[i].h) T[i]=MP[i]+(PART[i].h-hs1)/mass[i]/Cp[i];//�t��
		////////////*/
		
		double W=0;
		double nensei=1e8;
		double U[3][3];//��u[i]/��xj�i�[
		for(int n=0;n<3;n++) for(int m=0;m<3;m++) U[n][m]=0;//������

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
		
		double E[3][3];//��ij=0.5*(��ui/��xj+��uj/��xi)�i�[
		for(int n=0;n<3;n++) for(int m=0;m<3;m++) E[n][m]=0.5*(U[n][m]+U[m][n]);

		////�Ђ��݌v�Z
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
		
		vis[i]=nensei/density[i];					//���S���W��

		//���x����l������Ȃ�A�S���ɂ�锭�M���v�Z����B
		//heat[i]=V*sigma[i]*ep[i]*dt*co;//�Y���d������dW=��(d��)��co%���M�ɂȂ�Ɖ���	
		//PART[i].h+=heat[i];
		//�Y���d������dW=��(d��)��co%���M�ɂȂ�Ɖ���
		PART[i].heat_generation+=sigma[i]*ep[i]*co;		//heat_generation�ɔ��M��[W]���i�[


		/*/�M�ɂȂ����Ԃ�A���q�̉^���G�l���M�[��������
		double U=0;//���q�̑��x
		for(int D=0;D<d;D++) U+=PART[i].u[D]*PART[i].u[D];
		U=sqrt(U);
		double B=sqrt(2*heat/(V*density));//��������鑬�x
			
		if(U>1e-10 && U>B)
		{cout<<"���x����"<<endl;
			for(int D=0;D<d;D++) PART[i].u[D]-=PART[i].u[D]/U*B;
		}
		///*/

	}////*/
	//cout<<(GetTickCount()-timeA)*0.001<<endl;

	//���S���l�o��
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
				fr<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<heat[i]/V<<endl;//�P�ʑ̐ς�����̔��M�ʂ��o��[J/m3]
				fs<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<ep[i]<<endl;//�P�ʑ̐ς�����̔��M�ʂ��o��[J/m3]
				ft<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<sigma[i]<<endl;//�P�ʑ̐ς�����̔��M�ʂ��o��[J/m3]
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

//���xAVS�t�@�C���o�͊֐�
void output_viscousity_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number)
{
	char filename[30];
	int n=0;
	double le=CON->get_distancebp();
	double cross_section=CON->get_speed_face_p();
	//t=1;//���܂͂킴�Ɩ��X�e�b�v�㏑��

	//sprintf_s(filename,"pressure/pressure%d",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf_s(filename,"vis_XZ%d",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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
					double x=PART[i].r[A_X]*1.0E+05;	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
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
			double x=PART[i].r[A_X]*1.0E+05;	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
			double y=PART[i].r[A_Y]*1.0E+05;
			double z=PART[i].r[A_Z]*1.0E+05;

			//double x=PART[i].r[A_X];//*1.0E+05;	//r�͔��ɏ������l�Ȃ̂�10^5�{���Ă���
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
	//sprintf_s(filename,"pressure/pressure%d.fld",t);//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	sprintf_s(filename,"vis_XZ%d.fld",t);//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
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
	//fout2 << "variable 1 file=./pressure" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    1 file=./pressure" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    2 file=./pressure" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	//fout2 << "coord    3 file=./pressure" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//�t�H���_���쐬���ĊǗ�����ꍇ�͂�����
	fout2 << "variable 1 file=vis_XZ" << t << " " << "filetype=ascii offset=0 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    1 file=vis_XZ" << t << " " << "filetype=ascii offset=1 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    2 file=vis_XZ" << t << " " << "filetype=ascii offset=2 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2 << "coord    3 file=vis_XZ" << t << " " << "filetype=ascii offset=3 stride=4" << endl;//���̃t�@�C���Ɠ����K�w�ɐ�������Ȃ炱����
	fout2.close();
}



//�e���q�̕����l�����x�ɉ����ĕω�������//�Y�����̃f�[�^�x�[�X�@http://riodb.ibase.aist.go.jp/TPDB/DBGVsupport/detail/aluminum.html
void calc_physical_property(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *val,int particle_number,int sw)
{
	double MP=CON->get_MP();//�Z�_
	
	for(int i=0;i<fluid_number;i++)
	{
		double T=PART[i].T;

		if(CON->get_material()==Al)
		{
			if(sw==1)//���x
			{
				if(PART[i].T<=MP) val[i]=-0.207*PART[i].T+2757; //�Y�����̃f�[�^������`�ߎ�
				if(PART[i].T>MP) val[i]=2377.23-0.311*(PART[i].T-933.47);//�Y�����̋ߎ���
			}
			if(sw==2)//�S��
			{
				if(PART[i].T<=MP)
				{
					//cout<<"�S��:�ő�?"<<endl;
					//val[i]=1e10;//�ő�
					val[i]=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));
				}
					if(PART[i].T>MP) val[i]=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));//�Y�����̋ߎ���
			}
			if(sw==3)//�\�ʒ��͌W��
			{
				if(PART[i].T<=MP) val[i]=0;//�ő̂̏ꍇ�A�\�ʒ��͎��̂����݂��Ȃ�
				if(PART[i].T>MP) val[i]=-0.0001*T+0.988;//�A���~�j�E���Z�p�֗��̃f�[�^������`�ߎ�

			}
			if(sw==4)//��R��
			{
				if(PART[i].T<=MP) val[i]=-T*T*3e-14+T*9E-11-2e-9; //�Y�����̃f�[�^����2���������ߎ�
				if(PART[i].T>MP) val[i]=-T*T*2e-14+T*2E-10+8e-8;//�Y�����̃f�[�^����2���������ߎ�
			}
			if(sw==5)//���S��
			{
				if(PART[i].T<=MP)
				{
					//cout<<"���S��:�ő�?"<<endl;
					double val_d=-0.207*PART[i].T+2757; //�Y�����̃f�[�^������`�ߎ�
					double val_n=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));
					//double val_n=1e10;
					val[i]=val_n/val_d;
				}
				if(PART[i].T>MP)
				{
					double val_d=2377.23-0.311*(PART[i].T-933.47);//�Y�����̋ߎ���
					double val_n=0.001*(exp(log(10.0)*(-0.7324+803.49/PART[i].T)));
					val[i]=val_n/val_d;
				}
			}
		}
		else if(CON->get_material()==H2O)
		{
			if(sw==1)//���x
			{
				if(T<=MP) {}//�̂Ȃ牽�����Ȃ�
				//else if(T>MP) val[i]=(0.99986775 + (T-273.15)*6.7866875E-5 -pow(T-273.15,2.0)*9.09099173E-6 + pow(T-273.15,3.0)*1.02598151E-7 -pow(T-273.15,4.0)*1.3502904E-9 + pow(T-273.15,5.0)*1.32674392E-11 - pow(T-273.15,6.0)*6.461418E-14)*999.975; //�Y�����̃f�[�^������`�ߎ� DENS=(0.99986775+6.7866875E-5*(T-273.15)-9.09099173E-6*(T-273.15)^2+1.02598151E-7*(T-273.15)^3-1.3502904E-9*(T-273.15)^4+1.32674392E-11*(T-273.15)^5-6.461418E-14*(T-273.15)^6)*999.975
				else if(T>MP) val[i] = -T*T*T*T* 1.67669860531849E-6 + T*T*T*2.15183081454928E-3 - T*T*1.03632161524970 + T*2.21554244733677E2 - 1.67188884200996E4;
				//if(i==0) cout<<CON->get_density()<<" "<<val[i]<<endl;
			}
			if(sw==2)//�S��
			{
			}
			if(sw==3)//�\�ʒ��͌W��
			{
			}
			if(sw==4)//��R��
			{
			}
			if(sw==5)//���S��
			{
				if(T<=MP) {cout<<"H2O�ő�"<<endl;}//�̂Ȃ牽�����Ȃ�
				//else if(T>MP) val[i]=0.0005 - T*5E-6 + pow(T,2.0)*2E-8 - pow(T,3.0)*5E-11 +pow(T,4.0)*3E-14;
				else if(T>MP) val[i]= T*T*T*T*3.40023563631489E-14 - T*T*T*4.64957159378912E-11 + T*T*2.38796292580695E-8 - T*5.46534377544010E-6 + 4.71249468774775E-4;
				//if(i==0) cout<<"T="<<T<<" vis="<<CON->get_vis()<<" "<<val[i]<<endl;
			}
		}
	}

}

///�ړ����q�ړ��֐�
void move_particle(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double dt)
{
	cout<<"�ړ��Ǘ��q���ړ�"<<endl;
	int direction=CON->get_move_u_dirct();	//�ړ����q���ړ���������� ���݂́}X����=�}1,�}Y����=�}2,�}Z����=�}3
	double speed=CON->get_move_speed();		//�ړ����q�̈ړ����x[m/s]
	int D;									//�ړ������̖{�v���O�����ɂ�����Ή����鎟�� A_X=0;A_Y=1;A_Z=2;

	if(direction>0) D=direction-1;
	else if(direction<0)
	{
		D=-direction-1;
		speed*=-1;				//���x�𔽓]
	}

	for(int i=fluid_number;i<particle_number;i++)
	{
		if(PART[i].toBEM==MOVE)
		{
			PART[i].r[D]+=speed*dt;//���q���ړ�
		}
	}
}

///����t�@�C���o�͊֐�
void output_special_graph(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double TIME,int t)
{
	if(t==1 && CON->get_restart()==OFF && CON->get_model_number()==20)
	{
		ofstream fout("length.dat");//File�쐬
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

//���q�̐όv�Z�֐�
double get_volume(mpsconfig *CON)
{
	double V=0;//�̐�
	double le=CON->get_distancebp();
	if(CON->get_model_set_way()==0)	//�����i�q�̂Ƃ�
	{
		if(CON->get_dimention()==2){V=le*le;}
		else	{V=le*le*le;}
	}
	else if(CON->get_model_set_way()==1)	//�ז��i�q�̂Ƃ�
	{
		if(CON->get_dimention()==2){V=sqrt(3.0)/2*le*le;}
		else V=le*le*le/sqrt(2.0);
	}	
	else cout<<"���f���̐ςݕ����s��ł� �̐ς��v�Z�ł��܂���"<<endl;
	return V;
}

//���q�����x�̂��ꑪ��֐�
void output_particle_density(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double n0,int particle_number,int t)
{
	///�������q�����x���z���o�͂���

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
			//�̈�O�ɏo���ꍇ (�d�ɂ̒��ɖ�����Ă��폜)
			if(X<Xmin || X>Xmax) flag=1;
			if(Y<Ymin || Y>Ymax) flag=2;
			if(Z<Zmin || Z>Zmax) flag=3;//�d�ɂ̒��ɖ�����Ă��폜
			if(CON.get_model_number()==14)	if(PART[i].r[A_Z]<-le*2.0) flag=4;//(�Ód�����̏ꍇ)�d�ɂ̒��ɖ�����Ă��폜
			if(CON.get_model_number()==19)	if(PART[i].r[A_Z]>ZM*1.1) flag=8;//(fsw�̏ꍇ)�c�[������ɂ�������폜
			if(CON.get_model_number()==20 || CON.get_model_number()==21 || CON.get_model_number()==22)//(CC�n���̏ꍇ)��ڂɂ߂荞�񂾂�폜
			{
				if(Z>=0)
				{
					if(X*X+Y*Y>(0.03+le)*(0.03+le))
					{
						flag=4;
					}
				}
				else if(Z<0)//���ǉ~����
				{
					if(X*X+Y*Y+Z*Z>(0.03+le)*(0.03+le))
					{
						flag=4;
					}
				}

			}
			if(CON.get_model_number()==1)//(CC�n���̏ꍇ)��ڂɂ߂荞�񂾂�폜
			{
				if(X>0.0155 || X<-0.0155 ) flag=4;
				if(Y>0.0255 || Y<-0.0255 ) flag=4;
			}
			//�Ǘ���U���q�̍폜
			///*if(CON.get_process_iso()==1)*/	if(PART[i].fly==ISOLATION) flag=5;
			//���q�ʒu���s��(-1.#IND)�̂Ƃ�
			if(2*PART[i].r[A_X]!=PART[i].r[A_X]+PART[i].r[A_X]) flag=6;
			else if(2*PART[i].r[A_Y]!=PART[i].r[A_Y]+PART[i].r[A_Y]) flag=6;
			else if(2*PART[i].r[A_Z]!=PART[i].r[A_Z]+PART[i].r[A_Z]) flag=6;

			//ڲذ������p ������Ԃ̔�яo�����q������
			//if(CON.get_model_number()==20 && CON.get_current_step()==1 && PART[i].PND4<35)	flag=6;
			
			//���x���ߑ�
			double speed=0;//���q���x
			for(int D=0;D<DIMENTION;D++) speed+=PART[i].u[D]*PART[i].u[D];
			speed=sqrt(speed);
			if(speed>CON.get_max_speed()) flag=7;
			
		}

		if(flag>0)
		{
			//���W�o�������    <<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z];
			if(flag==1)	cout<<"�w�̈�O";
			if(flag==2)	cout<<"�x�̈�O";
			if(flag==3)	cout<<"�y�̈�O";
			if(flag==4)	cout<<"�ǂɐN��";
			if(flag==5)	cout<<"�Ǘ���U";
			if(flag==6)	cout<<"���W�s��";
			if(flag==7)	cout<<"�ő呬�x���ߑ�";
			if(flag==8)	cout<<"�c�[������ɑ���";
			//cout<<" i="<<i<<"/"<<(int)PART.size()<<" ��ٰ��"<<PART[i].group<<endl;
			//PART[i]���폜
			vector<mpsparticle>::iterator it=PART.begin();//�C�e���[�^������
			it+=i;				//i���w��
			it=PART.erase(it);	//�폜
			//cout<<" ��ٰ��"<<PART[i].group<<" -> �폜"<<endl;	//���Ȃ��������Ŏ~�܂�B�Ȃ��H �R�����g�A�E�g����ƒʂ�B
			*fluid_number-=1;
			*particle_number=(int)PART.size();

			count++;

			//���폜�����ꍇ�͔z�񂪋l�߂���̂�i++�͍s��Ȃ�
		}
		else i++;

		if(i>=(int)PART.size())	break;
	}

	if(count>0)	cout<<"���q"<<count<<"���폜\t"<<"���̗��q�F"<<begin_num<<" -> "<<begin_num-count<<endl;


	//���̗��q���̐��ڂ��O���t�ɏo��
	if(t==1)//1�X�e�b�v�ڂ̓t�@�C�������Z�b�g
	{
		ofstream freset("fluid_number.dat");
		freset.close();
	}
	ofstream fout("fluid_number.dat",ios::app);
	fout<<t<<"\t"<<begin_num-count<<endl;
	fout.close();
	
	//���̍ő卂�����o��
	if(t==1)//1�X�e�b�v�ڂ̓t�@�C�������Z�b�g
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

//�����𑪒�E��������֐��i�����͎����Ńv���O�������������Č��߂�j
void check_something(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double n0,int particle_number,int t)
{
	double le=CON->get_distancebp();
	if(CON->get_model_number()==19 && t%10==0)//FSW
	{
		int overN=0;//��ʂ���ɂ������̗��q��
		int underN=0;//��ʂ�艺�ɂ����c�[�����q��
		double standard_H=6e-3;//��ʂ�Z���W
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

///�Ǐd�݊֐��̌v�Z�֐�
void calc_wallZ(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,double *mindis,int fluid_number,int out)
{
	cout<<"�Ǐd�݊֐��̌v�Z�J�n"<<endl;
    ///�ǂ��ψ�ɔz�u�����ꍇ�̕Ǐd�݊֐����v�Z����B//model_number=0�ɋψ�ȕǗ��q�z�u��^���Ă���̂ŁA0�Ֆڂ̗��q�̍��W��ς��Ȃ��璲�ׂ�
	double le=CON->get_distancebp();
	int SIZE=CON->get_X_mesh()*CON->get_Y_mesh();
	
	//cout<<PART[0].r[A_Y]<<endl;
	//���x�����߁A�\�ʔ�����s���B�\�ʂȂ�P=0�ɂ���
	//omp_set_num_threads(8);//�گ�ސ��w��
	//#pragma omp parallel for
	int sol_wz=4000;//riw���ǂꂾ���ς��邩
	for(int wz=1;wz<sol_wz;wz++)
	{
		for(int i=0;i<1;i++)//OUTWALL�ȊO�̗��q�B//OUTWALL�̗��q�����x�Ȃǂ͂���Ȃ�
		{    
			//printf("%d %d\n",i,omp_get_thread_num());//�ei�̌v�Z��S�����Ă���گ�ޔԍ��o��
	       
			PART[i].PND=0;//������
			PART[i].PND2=0;
			PART[i].N=0;
			PART[i].N2=0;
			PART[i].N3=0;
			////���q�����x����
			double pnd=0;//���q�����x
			double pnd2=0;//���׼�ݗp���q�����x
			double pnd4=0;//�\�ʔ���p
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
					
				if(dis<=CON->get_re()*le && j!=i)//���z�E���U
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
				if(dis<=CON->get_re3()*le && j!=i)//�\�ʒ���re3
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
				//���̗��q���̐��ڂ��O���t�ɏo��
				if(wz==1)//1�X�e�b�v�ڂ̓t�@�C�������Z�b�g
				{
					ofstream fout("wallZ_pnd1.dat");
					fout.close();
				}
				if(wz==1)//1�X�e�b�v�ڂ̓t�@�C�������Z�b�g
				{
					ofstream fout2("wallZ_pnd2.dat");
					fout2.close();
				}
				if(wz==1)//1�X�e�b�v�ڂ̓t�@�C�������Z�b�g
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
	cout<<"�Ǐd�݂̌v�Z����"<<endl;

}