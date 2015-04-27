#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
#include"define.h"	//#define �i�[

#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
#include<vector>
#define FULL 1
#define HALF 2
#define HALFD 3
#define HALFD_shell 4

void writedata2(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM);
double get_volume(mpsconfig *CON);

//�����𕪊����邳���̍œK�ȕ������ƕ��������̎Z�o�֐� �����Ŏg������function.h�Ɉړ�
void calc_N_and_L(double dis,double le,int *N,double *L);
//�~���������v�Z�֐�
int calc_division_N_circle(double dis,double le);
//���aR�̉~�̊O��
void set_circle_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R);
//���aR�̉~����
void set_circle_in(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &type1,int *number,double le,double R,int edge_startID,int edge_lastID);
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int edge_startID,int edge_lastID);
//���aR�̋��쐬�֐�
//�����`�쐬�֐�
void set_rectangular(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height);
void set_sphere(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag);
void set_sphere2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag,int *suf_num);
//�����`�̕Ӎ쐬�֐�
void set_rectangular_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V);
//�����`�����쐬�֐�
void set_rectangular_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V,int edge_startID,int edge_lastID);
//�~���\�ʍ쐬�֐�
void set_cylinder_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int circle_start_id,int circle_end_id,int top_flag);
//�~���\�ʍ쐬�֐�
void set_circular_cone_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R2,double height,int circle_start_id,int circle_end_id,int top_flag);
//�~�������ݒu�֐�
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int flag);
//�h�[�i�c�쐬
void set_doughnut2D(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_big,double R_smal,int edge_startID,int edge_lastID);
//FSW�v���[�u�������q�Z�b�g�֐�
void set_hat_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID);
void set_hat_in_2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_mid,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID);
//���쐬�֐�
void set_box(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &surface,int *number,double le,double Width,double Height,double Depth);
//BOX���쐬�֐�
void set_box_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height,double Depth,int BO_startID,int BO_lastID);
//��ڕǓ��쐬�֐��i��ڂ̓������j
void set_crucible_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R_out,double height,double height_out,int flag,int fluid_number);

//���q���͊w�֐�
void MD_2D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN);
void MD_3D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN,double r,double region[3][2]);

//���������֐�
void make_fusion3D(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<double> &X2,vector<double> &Y2,vector<double> &Z2,vector<int> &surface2,int *number,double le);

void set_initial_placement_using_MD(mpsconfig *CON,int *particle_number)
{
	cout<<"�������q�z�u��MD�ɂ��œK����--";
	unsigned int timeA=GetTickCount();
	int number=0;	//���q���@�Ō�ɂ�particle_number�Ɋi�[
	int model=CON->get_model_number();
	int Dim=CON->get_dimention();		//��͎���

	double le=CON->get_distancebp();	//�������q�ԋ���
	double A=sqrt(3.0)*0.5;				//�悭�g���W��
	double B=sqrt(2.0/3);				//�悭�g���W�� ���������Ԋu
	

	vector<double> X;
	vector<double> Y;
	vector<double> Z;
	vector<int> type1;
	vector<int> surface;

	vector<double> X2;					//�g���񂵗p
	vector<double> Y2;
	vector<double> Z2;
	vector<int> surface2;			//ON�Ȃ�\�ʁ@OFF�Ȃ����

	vector<double> X3;					//�g���񂵗p �e�ۂQ�p
	vector<double> Y3;
	vector<double> Z3;
	vector<int> surface3;			//ON�Ȃ�\�ʁ@OFF�Ȃ����

	vector<double> X4;					//�g���񂵗p�@�e�ۂQ�p
	vector<double> Y4;
	vector<double> Z4;

	vector<double> X5;					//�g���񂵗p�@�e��1�p
	vector<double> Y5;
	vector<double> Z5;
	vector<int> surface5;			//ON�Ȃ�\�ʁ@OFF�Ȃ����
	
	vector<double> X6;					//�g���񂵗p�@�e��1�p
	vector<double> Y6;
	vector<double> Z6;
	
	ofstream fq("initial_input.dat");
	if(model==2)//////���f��2�@�������H
	{
		if(Dim==2)
		{
			double origin[3]={0,0,0};//�l�p�`�̍����̓_�̍��W
			double Width=CON->get_fluidwidth()*le;				//������������
			double Height=CON->get_fluidwidth()*le;				//������������
			set_rectangular(X,Y,Z,&number,le,Width,Height);		//���W�̌��_�͒����`�̒��S�B����Ă��Ƃňړ����邱�ƁB
			
			//������񏑂�����
			for(int i=0;i<number;i++) writedata2(fq,i,0.5*Width+origin[A_X]+X[i],Height+origin[A_Y]+Y[i],Z[i],FRFLUID,1,0,0,0,0,0,0,0,0);

		}
		else if(Dim==3)//�����̂̔�
		{
			double Width=CON->get_fluidwidth()*le;		
			double Height=CON->get_fluidwidth()*le;	
			double Depth=CON->get_fluidwidth()*le;
			set_box(X,Y,Z,surface,&number,le,Width,Height,Depth);//�Ō�3�̈����͉��A�����A���s���B�f�J���g���W�̌��_�ɑ΂��AX�������ɉ����AY�������ɉ��s���AZ�������ɍ���
			for(int i=0;i<number;i++) writedata2(fq,  i,X[i],Y[i],Z[i]+0.19,FLUID,1, OFF,0, 0,0, 0, 0, 0, 0);
		}
	}
	else if(model==3)//////���f���R�@�~�`���H
	{
		double R=CON->get_fluidwidth()*le*0.5;			//�쐬����~�̔��a
		double Zg=CON->get_height();					//���̒��S����
		if(Dim==2)
		{
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////�~����  vector�z��͎Q�Ɠn�����Ă���
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
		}
		else if(Dim==3)
		{
			//�~�쐬
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////�~����  vector�z��͎Q�Ɠn�����Ă���
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���

			//���쐬
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);
		}
		//������񏑂�����
		//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,Y[i]*40,0,0,0,0);
		for(int i=0;i<number;i++) writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,FLUID,1, OFF,0, 0,0, 0, 0, 0, ON);
		//for(int i=0;i<number;i++) writedata2(fq,  i,X[i],Y[i],Z[i],FLUID,1, OFF,0, 0,Y[i]*40, 0, 0, 0, ON);
	}
	else if(model==14)//14 �Ód����
	{
		double R=CON->get_fluidwidth()*le;
		double height=6*le*A;
		if(Dim==3)
		{
			//�~�쐬
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			int circle_end_id=number;	//�~�̗��qid���L��
			////////

			for(int i=0;i<number;i++)//X2�Ȃǂɉ~���q�̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			//�����쐬
			int flag=HALF;
			set_sphere(X,Y,Z,&number,le,R,flag);//�����쐬
			/////////////////////

			for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,0,0,0,0,0,0,1);//�����܂ł̗��q�͂��ׂė���

			int beforeN=number;
			int top_flag=ON;		//�~���̏�ʂ��쐬���邩��t���O��ON
			set_cylinder_face(X,Y,Z,&number,le,R,height,circle_start_id,circle_end_id,top_flag);//�~���\�ʍ쐬

			for(int i=beforeN;i<number;i++)//X2�Ȃǂɉ~���\�ʂ̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			for(int i=beforeN;i<number;i++) Z[i]*=-1.0;//Z���W�𔽓]
			for(int i=beforeN;i<number;i++)
			{
				if(Z[i]>-1.5*le) writedata2(fq,i,X[i],Y[i],Z[i],INWALL,1,0,0,0,0,0,0,0,1);//���q��INWALL
			}
			int num1=beforeN;//�����܂łɔ������Ă���OUTWALL�͂܂��������߂Ȃ��B�����Ő����L��
			int num2=number;

			

			beforeN=number;			//beforeN���X�V
			int number2=(int) X2.size();//X2�Ȃǂ��i�[���Ă��闱�q��
			int n=number2;		//�l��ۑ�
			set_cylinder_in(X2,Y2,Z2,&number2,le,R,height,1);//�~������ �܂���X2�Ȃǂɍ��W���i�[������
			
			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//�~�������̍��W��X,Y,Z�ɺ�߰�B������Z�͔��]����B
				Y.push_back(Y2[i]);
				Z.push_back(-1.0*Z2[i]);
				number++;
			}


			int count=beforeN;
			for(int i=beforeN;i<number;i++) 
			{
				if(Z[i]>-1.5*le) 
				{
					writedata2(fq,count,X[i],Y[i],Z[i],INWALL,1,0,0,0,0,0,0,0,1);//���q��INWALL
					count++;
				}
			}
			for(int i=beforeN;i<number;i++) 
			{
				if(Z[i]<=-1.5*le) 
				{
					writedata2(fq,count,X[i],Y[i],Z[i],OUTWALL,1,0,0,0,0,0,0,0,0);
					count++;
				}
			}
			for(int i=num1;i<num2;i++)
			{
				if(Z[i]<=-1.5*le) writedata2(fq,i,X[i],Y[i],Z[i],OUTWALL,1,0,0,0,0,0,0,0,1);//���q��OUTWALL
			}
			
		}
	}
	else if(model==16)//////���f��16�@�d�����V
	{
		double R=CON->get_fluidwidth()*le*0.5;			//�쐬����~�̔��a
		double Zg=CON->get_height();					//���̒��S����
		if(Dim==2)
		{
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////�~����  vector�z��͎Q�Ɠn�����Ă���
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
		}
		else if(Dim==3)
		{
			//�~�쐬
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////�~����  vector�z��͎Q�Ɠn�����Ă���
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���

			//���쐬
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);
		}
		//������񏑂�����
		for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zg,FLUID,1,0,0,0,0,0,0,0,1);
	}
	else if(model==20)//////��ڂƋ��`�n�Z����
	{
		double R=CON->get_fluidwidth()*0.001;//3;			//�쐬����~�̔��a
		double Zg=CON->get_height();					//���̒��S����
		double height=0.1;//10;							//��ډ~�����̍���
		double C_R=0.03;//3;			//��ڂ̔��a
		
		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//���x��֌W�̕ϐ���`
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		//double T=CON->get_roomT();
		double T=CON->get_initialT();//�������x
		double V=CON->get_particle_volume();
		double mass=density*V;	//���q����

		double h;//�G���^���s�[
		if(T>=CON->get_MP())
		{
			h=T*mass*Cp+latent_H*mass;
		}
		else h=T*mass*Cp;
									
		double wallmass=V*CON->get_wall_density();	//�Ǘ��q�̎���
		double wall_h=T*wallmass*CON->get_wall_Cp();//�Ǘ��q�̃G���^���s�[

		height=height/2;
		//C_R-=le;	//��1mm�̃R�[�e�B���O�w��z��				
		if(Dim==3)
		{
		//////////�n�Z����
			//���`
			//�~�쐬
			//R=C_R-le;		//��ڂɐڂ�����
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���

			//cout<<"�~�쐬"<<endl;

			//���쐬
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);

			//cout<<"���쐬"<<endl;

			int fluid_number=number;
			//cout<<"�����q��="<<fluid_number<<endl;

			if(CON->get_mesher()==0) for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zg+0.13125/*+0.18125-0.002*/,FLUID,materialID,1,0,0,0,0,0,h,1);
			if(CON->get_mesher()==1) for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zg,FLUID,materialID,1,0,0,0,0,0,h,1);

			//cout<<"�o�͊���="<<fluid_number<<endl;

		//////////���

			

			////����
			//��~
			int circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			int circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			//��ڑ��ʍ쐬
			int top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//�~��
			//�����k�쐬
			int flag2=HALFD_shell;//�������̋��k
			set_sphere(X,Y,Z,&number,le,C_R,flag2);

			int in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
		
		//�O��
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//��~
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			//set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			//��ڑ��ʍ쐬
			top_flag=HALF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//�~��
			//�����k�쐬
			flag2=HALFD_shell;//�������̋��k
			set_sphere(X,Y,Z,&number,le,C_R_out,flag2);
			//for(int i=in_number;i<number;i++) Z[i]-=4*le+C_R;

			set_crucible_in(X,Y,Z,&number,le,C_R,C_R_out,height,height_out,1,fluid_number);//��ړ����i�������j

			//OUTWALL�o��
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		///*////*/
		

		

	      ////////�o��
			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125+R,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			
			
		}

	}

	else if(model==26)//////IH�����f�� //�~����̊ȈՃ��f���@�{���͌��݂������̂�outwall�̗ʂɉ����𑜓x���l������K�v�����邪�A���肵�����x���������̂܂܋��E�����Ƃ��ė��p���邽�ߋC�ɂ����ɍ쐬
	{
		double R=CON->get_fluidwidth()*0.001;//3;			//
		double Zg=CON->get_height();					//
		double C_h=0.08;//10;							//�����ʕ��̍���
		double C_hout=C_h+4*le;//10;							//�����ʕ��̍���
		double C_R=0.07;//3;			//���ꕔ�̔��a
		double C_Rout=C_R+4*le;

		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//���x��֌W�̕ϐ���`
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		//double T=CON->get_roomT();
		double T=CON->get_initialT();//�������x
		double V=CON->get_particle_volume();
		double mass=density*V;	//���q����

		double h;//�G���^���s�[
		if(T>=CON->get_MP())
		{
			h=T*mass*Cp+latent_H*mass;
		}
		else h=T*mass*Cp;
									
		double wallmass=V*CON->get_wall_density();	//�Ǘ��q�̎���
		double wall_h=T*wallmass*CON->get_wall_Cp();//�Ǘ��q�̃G���^���s�[
			

		if(Dim==2)
		{
		}
		if(Dim==3)
		{
		
		//////////���Ɨ��̂𓯎��ɍ쐬����

			int fluid_number=0;

			//�~�쐬
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_Rout,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			int circle_end_id=number;	//�~�̗��qid���L��
			////////

			for(int i=0;i<number;i++)//X2�Ȃǂɉ~���q�̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			/////////////////////

			int beforeN=number;
			int top_flag=OFF;		//�~���̏�ʂ��쐬���邩��t���O��ON
			set_cylinder_face(X,Y,Z,&number,le,C_Rout,C_hout,circle_start_id,circle_end_id,top_flag);//�~���\�ʍ쐬

			for(int i=beforeN;i<number;i++)//X2�Ȃǂɉ~���\�ʂ̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////
			
			beforeN=number;			//beforeN���X�V
			int number2=(int) X2.size();//X2�Ȃǂ��i�[���Ă��闱�q��
			int n=number2;		//�l��ۑ�
			set_cylinder_in(X2,Y2,Z2,&number2,le,C_Rout,C_hout,1);//�~������ �܂���X2�Ȃǂɍ��W���i�[������

			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//�~�������̍��W��X,Y,Z�ɺ�߰�B
				Y.push_back(Y2[i]);
				Z.push_back(Z2[i]);
				number++;
			}

			for(int i=0;i<number;i++) Z[i]-=C_hout-C_h;
			
			///////////////////���q�̍ގ�����

			for(int i=0;i<number;i++) type1.push_back(FLUID);//�܂��͂Ƃ肠����FLUID����

			for(int i=0;i<number;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>C_R-0.4*le)//���̑���
				{
					if(r<C_R+0.4*le && Z[i]>0.2*le) type1[i]=INWALL;
					else type1[i]=OUTWALL;	
				}
				else//���̓�������ђ��
				{
					if(Z[i]<0.2*le && Z[i]>-0.2*le)
					{
						type1[i]=INWALL;
					}
					else if(Z[i]<=-0.2*le) type1[i]=OUTWALL;
				}
			}
			//////////////////////*/


			//���q���W�o��
			int count=0;
			int count_outf=0;//fluid�̂����AZ���W�ŏ��O����鐔
			for(int i=0;i<number;i++)
			{
				//if(type1[i]==FLUID)
				if(type1[i]==FLUID)
				{
					if(Z[i]<=0.05)
					{
						count++;
						writedata2(fq,i,X[i],Y[i],Z[i],FLUID,materialID,0,0,0,0,0,0,h,1);
					}
					else count_outf++;
				}
			}

			for(int i=0;i<number;i++)
			{
				if(type1[i]==INWALL)
				{
					count++;
					writedata2(fq,i,X[i],Y[i],Z[i],INWALL,materialID,0,0,0,0,0,0,wall_h,1);
				}
			}
					for(int i=0;i<number;i++)
			{
				if(type1[i]==OUTWALL)
				{
					count++;
					writedata2(fq,i,X[i],Y[i],Z[i],OUTWALL,materialID,0,0,0,0,0,0,wall_h,1);
				}
			}

			number-=count_outf;

	      ////////�o��
			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			//for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125+R,FLUID,materialID,0,0,0,0,0,0,0,1);
			//for(int i=fluid_number;i<in_number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,materialID,0,0,0,0,0,0,0,1);
			//for(int i=in_number;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,materialID,0,0,0,0,0,0,0,1);

			
			
		}

	}
	
	else if(model==21)//��ڂ̒��ɒe�ی^�̗��́F
	{
		double R=CON->get_fluidwidth()*0.001;//3; //���̂̔��a
		//double R=0.03-le/sqrt(2.0);//3; //���̂̔��a
		double fluid_h=CON->get_fluid_h()*0.001;//3; //���̉~�����̍���
		double height=0.1;//10;							//��ډ~�����̍���
		double C_R=0.03;//3;			//��ڂ̔��a
		//double Zf=0.13125-(0.03-R)+le/sqrt(2.0);//���̂�Z�����̈ʒu�C����
		double Zf=0.13125;

		height/=2.0;//��̂ق��̕ǂ��J�b�g

		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//���x��֌W�̕ϐ���`
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		double T=CON->get_roomT();
		double V=CON->get_particle_volume();
		double mass=density*V;									//���q����
		double h=CON->get_MP()*mass*Cp+latent_H*mass;//�G���^���s�[
		//double h;								//�G���^���s�[
		double wallmass=V*CON->get_wall_density();	//�Ǘ��q�̎���
		double wall_h=T*wallmass*CON->get_wall_Cp();//�Ǘ��q�̃G���^���s�[
		

		if(Dim==3)
		{
			//����

			//�~�쐬
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			int circle_end_id=number;	//�~�̗��qid���L��
			////////

			for(int i=0;i<number;i++)//X2�Ȃǂɉ~���q�̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////


			//�����쐬
			int beforeN=number;
			int flag=HALF;
			set_sphere(X,Y,Z,&number,le,R,flag);//�����쐬

			for(int i=beforeN;i<number;i++) Z[i]*=-1.0;//Z���W�𔽓]
			/////////////////////

			beforeN=number;
			int top_flag=ON;		//�~���̏�ʂ��쐬���邩��t���O��ON
			set_cylinder_face(X,Y,Z,&number,le,R,fluid_h,circle_start_id,circle_end_id,top_flag);//�~���\�ʍ쐬

			for(int i=beforeN;i<number;i++)//X2�Ȃǂɉ~���\�ʂ̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////
			
			beforeN=number;			//beforeN���X�V
			int number2=(int) X2.size();//X2�Ȃǂ��i�[���Ă��闱�q��
			int n=number2;		//�l��ۑ�
			set_cylinder_in(X2,Y2,Z2,&number2,le,R,fluid_h,1);//�~������ �܂���X2�Ȃǂɍ��W���i�[������

			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//�~�������̍��W��X,Y,Z�ɺ�߰�B
				Y.push_back(Y2[i]);
				Z.push_back(Z2[i]);
				number++;
			}
			
			for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+Zf,FLUID,1,0,0,0,0,0,0,h,1);//�����܂ł̗��q�͂��ׂė���
			//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]+0.13125-le-le/sqrt(2.0),FLUID,1,0,0,0,0,0,0,0,1);//�����܂ł̗��q�͂��ׂė���

			int fluid_number=number;
			//X.clear(); Y.clear(); Z.clear(); X2.clear(); Y2.clear(); Z2.clear();
			cout<<"���̗��q�o��"<<endl;

		//���
			//number=0;
		////����
			//��~
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			//��ڑ��ʍ쐬
			top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//�~��
			//�����k�쐬
			int flag2=HALFD_shell;//�������̋��k
			set_sphere(X,Y,Z,&number,le,C_R,flag2);

			int in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
		
		//�O��
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//��~
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			//set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			//��ڑ��ʍ쐬
			top_flag=HALF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//�~��
			//�����k�쐬
			flag2=HALFD_shell;//�������̋��k
			set_sphere(X,Y,Z,&number,le,C_R_out,flag2);
			//for(int i=in_number;i<number;i++) Z[i]-=4*le+C_R;

			int mostout_number=number;

			set_crucible_in(X,Y,Z,&number,le,C_R,C_R_out,height,height_out,1,fluid_number);//��ړ����i�������j

			//OUTWALL�o��
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		}
	}

	else if(model==23)//�ǋ��E��r���f��
	{
		double R=0.1;//3; //���̂̔��a
		//double R=0.03-le/sqrt(2.0);//3; //���̂̔��a
		double fluid_h=0.05;//3; //���̉~�����̍���
		double height=0.1;//10;							//��ډ~�����̍���
		double C_R=0.1+le;//3;			//��ڂ̔��a
		//double Zf=0.13125-(0.03-R)+le/sqrt(2.0);//���̂�Z�����̈ʒu�C����
		double Zf=0;//0.13125;

		

		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//���x��֌W�̕ϐ���`
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		double T=CON->get_roomT();
		double V=CON->get_particle_volume();
		double mass=density*V;									//���q����
		double h=CON->get_MP()*mass*Cp+latent_H*mass;//�G���^���s�[
		//double h;								//�G���^���s�[
		double wallmass=V*CON->get_wall_density();	//�Ǘ��q�̎���
		double wall_h=T*wallmass*CON->get_wall_Cp();//�Ǘ��q�̃G���^���s�[
		

		if(Dim==3)
		{
			//����

			//�~�쐬
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			int circle_end_id=number;	//�~�̗��qid���L��
			////////

			for(int i=0;i<number;i++)//X2�Ȃǂɉ~���q�̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////

			int beforeN=number;
			int top_flag=ON;		//�~���̏�ʂ��쐬���邩��t���O��ON
			set_cylinder_face(X,Y,Z,&number,le,R,fluid_h,circle_start_id,circle_end_id,top_flag);//�~���\�ʍ쐬

			for(int i=beforeN;i<number;i++)//X2�Ȃǂɉ~���\�ʂ̏����߰
			{
				X2.push_back(X[i]);
				Y2.push_back(Y[i]);
				Z2.push_back(Z[i]);
			}//////////////
			
			beforeN=number;			//beforeN���X�V
			int number2=(int) X2.size();//X2�Ȃǂ��i�[���Ă��闱�q��
			int n=number2;		//�l��ۑ�
			set_cylinder_in(X2,Y2,Z2,&number2,le,R,fluid_h,1);//�~������ �܂���X2�Ȃǂɍ��W���i�[������

			for(int i=n;i<number2;i++)
			{
				X.push_back(X2[i]);			//�~�������̍��W��X,Y,Z�ɺ�߰�B
				Y.push_back(Y2[i]);
				Z.push_back(Z2[i]);
				number++;
			}

			int in_number=number;
			
			for(int i=0;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i],FLUID,1,0,0,0,0,0,0,h,1);
			}
		
			int fluid_number=number;
			cout<<"���̗��q�o��"<<endl;

			////����
			//��~
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			//��ڑ��ʍ쐬
			top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//�~��

			in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]-le,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
			cout<<"����"<<endl;
		//�O��
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//��~
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			//��ڑ��ʍ쐬
			top_flag=HALF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//�~��

			int mostout_number=number;

			set_cylinder_in(X,Y,Z,&number,le,R,fluid_h,1);//�~������ �܂���X2�Ȃǂɍ��W���i�[������

			cout<<"�O�Ǔ���"<<endl;


			//OUTWALL�o��
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]-le,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		}
	}
	else if(model==25)//////�Q�d����͊ȈՃ��f��
	{
		double R=CON->get_fluidwidth()*0.001;//3;			//�쐬����~�̔��a
		double Zg=CON->get_height();					//���̒��S����
		double height=0.1;//10;							//��ډ~�����̍���
		double C_R=0.03;//3;			//��ڂ̔��a
		
		int materialID=1;
		double Cp;
		double density;
		double latent_H; 
		//���x��֌W�̕ϐ���`
		if(materialID==1)
		{
			density=CON->get_density();
			Cp=CON->get_Cp();
			latent_H=CON->get_latent_H();
		}
		else if(materialID==2)
		{
			density=CON->get_density2();
			Cp=CON->get_Cp2();
			latent_H=CON->get_latent_H2();
		}

		
		//double T=CON->get_roomT();
		double T=CON->get_initialT();//�������x
		double V=CON->get_particle_volume();
		double mass=density*V;	//���q����

		double h;//�G���^���s�[
		if(T>=CON->get_MP())
		{
			h=T*mass*Cp+latent_H*mass;
		}
		else h=T*mass*Cp;
									
		double wallmass=V*CON->get_wall_density();	//�Ǘ��q�̎���
		double wall_h=T*wallmass*CON->get_wall_Cp();//�Ǘ��q�̃G���^���s�[

		height=height/2;
		//C_R-=le;	//��1mm�̃R�[�e�B���O�w��z��				
		if(Dim==3)
		{
		//////////�n�Z����
			//���`
			//�~�쐬
			//R=C_R-le;		//��ڂɐڂ�����
			set_circle_edge(X,Y,Z,&number,le,R);//�~�O��
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���

			//cout<<"�~�쐬"<<endl;

			//���쐬
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);

			//cout<<"���쐬"<<endl;

			int fluid_number=number;
			//cout<<"�����q��="<<fluid_number<<endl;

			
			for(int i=0;i<fluid_number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FLUID,materialID,1,0,0,0,0,0,h,1);

			//cout<<"�o�͊���="<<fluid_number<<endl;

		/*/////////���

			

			////����
			//��~
			int circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			int circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			//��ڑ��ʍ쐬
			int top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R,height,circle_start_id,circle_end_id,top_flag);//�~��
			//�����k�쐬
			int flag2=HALFD_shell;//�������̋��k
			set_sphere(X,Y,Z,&number,le,C_R,flag2);

			int in_number=number;
			
			for(int i=fluid_number;i<in_number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,INWALL,1,0,0,0,0,0,0,wall_h,1);
			}
		
		//�O��
			double C_R_out=C_R+4*le;
			double height_out=height+4*le+C_R;
			//��~
			circle_start_id=number;
			set_circle_edge(X,Y,Z,&number,le,C_R_out);//�~�O��   vector�z��͎Q�Ɠn�����Ă��� ���̉~����ɉ~���A���k������
			circle_end_id=number;//�~���ɔz�u���ꂽ���q��
			//set_circle_in_using_6_pieces(X,Y,Z,&number,le,C_R_out,circle_start_id,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			//��ڑ��ʍ쐬
			top_flag=HALF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,C_R_out,height,circle_start_id,circle_end_id,top_flag);//�~��
			//�����k�쐬
			flag2=HALFD_shell;//�������̋��k
			set_sphere(X,Y,Z,&number,le,C_R_out,flag2);
			//for(int i=in_number;i<number;i++) Z[i]-=4*le+C_R;

			set_crucible_in(X,Y,Z,&number,le,C_R,C_R_out,height,height_out,1,fluid_number);//��ړ����i�������j

			//OUTWALL�o��
			for(int i=in_number;i<number;i++)
			{
				writedata2(fq,i,X[i],Y[i],Z[i]+0.13125,OUTWALL,1,0,0,0,0,0,0,wall_h,1);
			}			
		///*////*/	
		}
	}
	else if(model==19)//19 FSW
	{
		double probe_R=2.5*1e-3;	//�v���[�u���a  //CON->get_fluidwidth()*le;
		double probe_R_2=1.0*1e-3;	//�v���[�u���a(�~���`��̉���)  tool_type=1�̂Ƃ��ɗ��p;
		double shold_R=6*1e-3;	//�V�����_�[���a
		double height=4*1e-3;	//�v���[�u���� //6*le*A;
		if(CON->get_tool_type()==2) height=3*1e-3+B*le;//�\���c�[���̏ꍇ�A�c�[���̒���*2�����̗̈�̌��݂Ɠ����ɂȂ�Ȃ��Ƃ����Ȃ�

		double shold_height=10*B*le;
		double rpm=500;//�c�[����]���x
		double rps=rpm/60;
		double w=rps*2*PI;		//�p���x
		double U=CON->get_move_speed();//�v���[�u�̈ړ����x[m/sec]		
		double pich=0.7e-3;		//�v���[�u�̂˂��̃s�b�`0.7[mm]
		double T=CON->get_roomT();	//�������x
		double vol=get_volume(CON);				//���q�̑̐�
		double h;								//�G���^���s�[
		double wallmass=vol*CON->get_wall_density();	//�Ǘ��q�̎���
		double wall_h=T*wallmass*CON->get_wall_Cp();//�Ǘ��q�̃G���^���s�[
		double val=0;
		int materialID=1;
		int count2;
		int beforeN;

		//calclation type
		int plunge=0;
		int traverse=1;
		int calc_type=CON->get_process_type();		//0:plunge 1:traverse
		if(Dim==3)
		{
			///�v���[�u��ʍ쐬
			if(CON->get_tool_type()==0 ||CON->get_tool_type()==2)
			{
				set_circle_edge(X,Y,Z,&number,le,probe_R);//�~�O��
				set_circle_in_using_6_pieces(X,Y,Z,&number,le,probe_R,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			}
			if(CON->get_tool_type()==1)//�~���c�[��
			{
				set_circle_edge(X,Y,Z,&number,le,probe_R_2);//�~�O��
				set_circle_in_using_6_pieces(X,Y,Z,&number,le,probe_R_2,0,number);//�~����    vector�z��͎Q�Ɠn�����Ă���
			}
			int circle_start_id=0;
			int circle_end_id=number;
			int circle_num=number;		//�~���ɔz�u���ꂽ���q��
			////////
			int top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			//�v���[�u���ʍ쐬
			if(CON->get_tool_type()==0 ||CON->get_tool_type()==2)
			{
				top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
				set_cylinder_face(X,Y,Z,&number,le,probe_R,height,circle_start_id,circle_end_id,top_flag);//�~���\�ʍ쐬
			}
			if(CON->get_tool_type()==1)
			{
				top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
				set_circular_cone_face(X,Y,Z,&number,le,probe_R,probe_R_2,height,circle_start_id,circle_end_id,top_flag);//�~���\�ʍ쐬
			}

			//�V�����_�[��ʍ쐬
			beforeN=number;
			set_doughnut2D(X,Y,Z,&number,le,shold_R,probe_R,number-circle_num,number);//�Ō�4�̈����͑傫�����a�A���������a�A���łɍ쐬���Ă���~����startID,endID
			for(int i=beforeN;i<number;i++) Z[i]=height;//��̊֐��ō쐬�����O����Z���W�̓[���̂܂܂Ȃ̂ŁA�����ŏC��
			int shld_botom_startID=beforeN;
			int shld_botom_endID=number;

			//�V�����_�[���ʍ쐬
			beforeN=number;
			top_flag=OFF;		//�~���̏�ʂ͍쐬���Ȃ��Ă�������t���O��OFF 
			set_cylinder_face(X,Y,Z,&number,le,shold_R,shold_height,shld_botom_startID,shld_botom_endID,top_flag);//�~���\�ʍ쐬 
			for(int i=beforeN;i<number;i++) Z[i]+=height;

			//�V�����_�[���
			count2=0;//�V�����_��ʂɂ����鋫�E���q�́A���set_cylinder_face()�ōŌ�ɕt��������ꂽ���q�Q�ł���B�܂��͂��̐��𐔂���
			double H=height+shold_height;//�V�����_�[�̂Ă��؂�̍���
			for(int i=beforeN;i<number;i++) if(Z[i]<H+0.1*le && Z[i]>H-0.1*le) count2++;
			beforeN=number;
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,shold_R,number-count2,number);//�~����
			for(int i=beforeN;i<number;i++) Z[i]=H;

			for(int i=0;i<number;i++) type1.push_back(INWALL);//�����܂ł̗��q�͂��ׂ�INWALL

			//�����쐬
			beforeN=number;
			if(CON->get_tool_type()==0 ||CON->get_tool_type()==2) set_hat_in(X,Y,Z,&number,le,probe_R,shold_R,height,shold_height,0,number);
			if(CON->get_tool_type()==1) set_hat_in_2(X,Y,Z,&number,le,probe_R_2,probe_R,shold_R,height,shold_height,0,number);

			for(int i=beforeN;i<number;i++) type1.push_back(OUTWALL);//�����܂ł̗��q�͂��ׂ�OUTWALL
			
			beforeN=number;

			if(CON->get_tool_type()==2)///
			{
				for(int i=0;i<beforeN;i++)
				{
					if(Z[i]!=0)//Z=0�ʂ����\�R�s�[�̋��E��
					{
						X.push_back(X[i]);
						Y.push_back(Y[i]);
						Z.push_back(-Z[i]);
						type1.push_back(type1[i]);
						number++;
					}
					double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
					if(Z[i]==0 && r<probe_R-0.1*le)//���E�ʂ̊O�����ȊO
					{
						type1[i]=OUTWALL;
					}
				}
			}
			
			double max_dis_z=0.0;
			if(CON->get_tool_angle()>0)//�c�[����x���܂��ɉ�]������B�i�s������+y
			{
				double theta=PI/180*CON->get_tool_angle();//��]����p�x
				
				for(int i=0;i<number;i++)
				{
					X3.push_back(X[i]);//�c�[���p�x��ς����ۂ̑��x�����̐ݒ���ȒP�ɂ��邽�߁A��]�O�̍��W���L��
					Y3.push_back(Y[i]);
					Z3.push_back(Z[i]+2*1e-3);
					//////////////
					double x=X[i];
					double y=cos(theta)*Y[i]-sin(theta)*Z[i];
					double z=sin(theta)*Y[i]+cos(theta)*Z[i];
					double dis_z=z-Z[i];
					if(max_dis_z<dis_z) max_dis_z=dis_z;
					X[i]=x;
					Y[i]=y;
					Z[i]=z;

				}		
			}

			if(CON->get_tool_type()!=2) for(int i=0;i<number;i++) Z[i]+=4*1e-3;//�d�S���ړ�

			if(CON->get_tool_angle()>0) 
			{
				double theta=PI/180*CON->get_tool_angle();//��]����p�x
				for(int i=0;i<number;i++)
				{
					//Z[i]-=max_dis_z;//�d�S���ړ�
					Z[i]+=0.006*(1-cos(theta));//�d�S���ړ�
				}
			}
			
			//////////////////////////////�c�[���쐬�I��*/

			int tool_number=number;
			
			//�}�����̉�͉͂��̍s��ON�ɂ��A�c�[���̈ʒu��������
			if(calc_type==plunge) for(int i=0;i<number;i++) Z[i]+=4*1e-3;//�d�S���ړ�

			/////////////////////////////���̍쐬
			int fluid_number=0;
			double width0=18*1e-3;//25*1e-3;//18*1e-3;//���̂���߂镝
			double Width=width0+6*le;		//��18mm+�Ǘ��q�����E��4���q��
			double Height=6*1e-3+4*le*B;	//����6mm+�Ǘ��q������4���q��
			if(CON->get_tool_type()==2) Height=6*1e-3;	//�\���c�[���̏ꍇ�A�����̕ǂ���蕥���Ă���̂ŕǗ��q�̐ݒ�̊֌W�ŉ����ɒǉ�����K�v���Ȃ�

			double Depth=18*1e-3+6*le*A;//27*1e-3+6*le*A;	//���s��27mm+�Ǘ��q�����E��4���q��
			double depth0=9e-3;				//�c�[�����S�ƁA��O�̕ǂƂ̋���
			int BOX_startID=0;				//set_box()�J�n�O�̗��q��
			set_box(X2,Y2,Z2,surface2,&fluid_number,le,Width,Height,Depth);//�Ō�3�̈����͉��A�����A���s���B�f�J���g���W�̌��_�ɑ΂��AX�������ɉ����AY�������ɉ��s���AZ�������ɍ���
			int BOX_lastID=fluid_number;	//set_box()�I������̗��q��

			//double Z_mod=le*B+(1-sqrt(2.0)/2)*le;//�c�[�������̂ւ߂荞�ނ̂�h�����߂̒����ʁ@���̒l�����c�[���ȊO�̕��̂�������
			double Z_mod=le+(1-sqrt(2.0)/2)*le;
			
			for(int i=BOX_startID;i<BOX_lastID;i++)
			{
				X2[i]-=width0*0.5+3*le;		//�d�S�ړ�
				Y2[i]-=depth0+3*le*A;
				Z2[i]-=4*le*B+Z_mod;
				if(CON->get_tool_type()==2) Z2[i]+=4*le*B+Z_mod-0.5*Height;
				//Z2[i]-=4*le*B+le;//����1���q�������āA�V�����_�[�����̕�����������Ȃ��悤�ɂ���// B��������1�w�����ƐڐG�������H
			}/////////////////////////////////////////
			
			////////////////////////////�c�[���Ɣ�闬�̗��q���폜����
			beforeN=number;
			make_fusion3D(X,Y,Z,X2,Y2,Z2,surface2,&number,le);
			////////////////////////////////////

			///////////////////box���q�̍ގ�����

			for(int i=beforeN;i<number;i++) type1.push_back(FRFLUID);//�܂��͂Ƃ肠����FRFLUID����


			if(CON->get_tool_type()==0 ||CON->get_tool_type()==1)
			{
				Width-=6*le;
				Height-=4*le*B;
				Depth-=6*le*A;
				for(int i=beforeN;i<number;i++)
				{
					if(X[i]>(Width*0.5)+0.2*le)//�E�̕�
					{
						if(X[i]<(Width*0.5)+1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
					else if(X[i]<-(Width*0.5)-0.2*le)//���̕�
					{
						if(X[i]>-(Width*0.5)-1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
				
					if(Y[i]<-(depth0)-0.2*le)//��O�̕�
					{
						//if(Y[i]>-(depth0)-1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]>-(depth0)-1.2*le && Z[i]>-1.2*le-Z_mod) 	
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;
						}
						else type1[i]=OUTWALL;
					}
					else if(Y[i]>(Depth-depth0)+0.2*le)//���̕�
					{
						//if(Y[i]<(Depth-depth0)+1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]<(Depth-depth0)+1.2*le && Z[i]>-1.2*le-Z_mod)
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;// ����if���Ȃ��ƁA���ɍ��E��outwall�Ɣ��f����Ă���ꕔ��inwall�ɂȂ��Ă��܂�
						}
						else type1[i]=OUTWALL;
					}

					if(CON->get_tool_type()==0 || CON->get_tool_type()==1)
					{
						if(X[i]<=(Width*0.5)+0.2*le && X[i]>=-(Width*0.5)-0.2*le)
						{
							if(Y[i]>=-(depth0)-0.2*le && Y[i]<=(Depth-depth0)+0.2*le)
							{
								if(Z[i]<-0.2*le-Z_mod)
								{
									if(Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;//���̕�
									else type1[i]=OUTWALL;
								}
							}
						}
					}
				}
			}
			if(CON->get_tool_type()==2)//�\���c�[��
			{
				Width-=6*le;
				//Height-=4*le*B;
				Depth-=6*le*A;
				for(int i=beforeN;i<number;i++)
				{
					if(X[i]>(Width*0.5)+0.2*le)//�E�̕�
					{
						if(X[i]<(Width*0.5)+1.2*le) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
					else if(X[i]<-(Width*0.5)-0.2*le)//���̕�
					{
						if(X[i]>-(Width*0.5)-1.2*le) type1[i]=INWALL;
						else type1[i]=OUTWALL;	
					}
				
					if(Y[i]<-(depth0)-0.2*le)//��O�̕�
					{
						//if(Y[i]>-(depth0)-1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]>-(depth0)-1.2*le) 	
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;
						}
						else type1[i]=OUTWALL;
					}
					else if(Y[i]>(Depth-depth0)+0.2*le)//���̕�
					{
						//if(Y[i]<(Depth-depth0)+1.2*le && Z[i]>-1.2*le-Z_mod) type1[i]=INWALL;
						if(Y[i]<(Depth-depth0)+1.2*le)
						{
							if(X[i]>-(Width*0.5)-0.2*le && X[i]<(Width*0.5)+0.2*le) type1[i]=INWALL;// ����if���Ȃ��ƁA���ɍ��E��outwall�Ɣ��f����Ă���ꕔ��inwall�ɂȂ��Ă��܂�
						}
						else type1[i]=OUTWALL;
					}
					////�����̕ǂ͕K�v�Ȃ�
				}
			}
			//////////////////////*/


			//���q���W�o��
			count2=0;
			double Zmax=0;
			for(int i=0;i<number;i++)
			{
				if(type1[i]==FRFLUID)
				{
					//if(X[i]>0) materialID=1;
					//else       materialID=2;
					materialID=1;//�Ƃ肠�������ސڍ�
					double density,Cp,latent_H;
					if(materialID==1)
					{
						density=CON->get_density();
						Cp=CON->get_Cp();
						latent_H=CON->get_latent_H();
					}
					else if(materialID==2)
					{
						density=CON->get_density2();
						Cp=CON->get_Cp2();
						latent_H=CON->get_latent_H2();
					}
					double mass=density*vol;									//���q����
					//h=T*mass*Cp+latent_H*mass;//�G���^���s�[

					if(T>=CON->get_MP())
					{
						h=T*mass*Cp+latent_H*mass;
					}
					else h=T*mass*Cp;
					writedata2(fq, count2, X[i],Y[i],Z[i], FLUID,materialID,OFF, val,0,0,0,0,h,1);//�����܂ł̗��q�͂��ׂė���

					count2++;
					if(Z[i]>=Zmax) Zmax=Z[i];
				}
			}

			//�c�[�������̑��x������^������ŏo��
			if(CON->get_tool_angle()==0)//�c�[������]�����Ȃ��ꍇ
			{
				if(CON->get_tool_type()==0 || CON->get_tool_type()==1)
				{
					for(int i=0;i<number;i++)
					{
					
						if(type1[i]==INWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//�c�[���̒��S�����_�łȂ��Ƃ��͒���
							int AA=1;
							if(r>0.1*le && r<=shold_R+le && Z[i]>0)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								if(Z[i]<Height && r>probe_R-0.3*le) uw=-pich*rps;//INWALL�̂����A�v���[�u�̑��ʂ̂�(��ʂ͏���)�ɉ��������x�ǉ�
								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								if(calc_type==plunge) uw-=U;//Z�����Ƀc�[�����ړ�(plange phase)
							}
							if(r<=shold_R+le && Z[i]>0) AA=MOVE;
							materialID=1;		
							writedata2(fq, count2, X[i],Y[i],Z[i], INWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
							count2++;
						}
					}
					for(int i=0;i<number;i++)
					{
						if(type1[i]==OUTWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//�c�[���̒��S�����_�łȂ��Ƃ��͒���
							int AA=1;
							if(r>0.1*le && r<=shold_R+le && Z[i]>0)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								uw=-pich*rps;
								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								else if(calc_type==plunge) uw-=U;//Z�����Ƀc�[�����ړ�(plange phase)
							}
							if(r<=shold_R+le && Z[i]>0) AA=MOVE;
							materialID=1;
							writedata2(fq, count2, X[i],Y[i],Z[i], OUTWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);

							count2++;
						}
					}
				}
				else if(CON->get_tool_type()==2)
				{
					for(int i=0;i<number;i++)
					{
					
						if(type1[i]==INWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//�c�[���̒��S�����_�łȂ��Ƃ��͒���
							int AA=1;
							if(r>0.1*le && r<=shold_R+le)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								//if(Z[i]<Height && r>probe_R-0.3*le) uw=-pich*rps;//INWALL�̂����A�v���[�u�̑��ʂ̂�(��ʂ͏���)�ɉ��������x�ǉ�
								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								if(calc_type==plunge) uw-=U;//Z�����Ƀc�[�����ړ�(plange phase)
							}
							if(r<=shold_R+le) AA=MOVE;
							materialID=1;		
							writedata2(fq, count2, X[i],Y[i],Z[i], INWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
							count2++;
						}
					}
					for(int i=0;i<number;i++)
					{
						if(type1[i]==OUTWALL)
						{
							double speed=0;
							double u=0;
							double v=0;
							double uw=0;
							double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);//�c�[���̒��S�����_�łȂ��Ƃ��͒���
							int AA=1;
							if(r>0.1*le && r<=shold_R+le)
							{
								speed=r*w;
								u=speed*(-Y[i]/r);
								v=speed*(X[i]/r);
								uw=-pich*rps;
								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								else if(calc_type==plunge) uw-=U;//Z�����Ƀc�[�����ړ�(plange phase)
							}
							if(r<=shold_R+le) AA=MOVE;
							materialID=1;
							writedata2(fq, count2, X[i],Y[i],Z[i], OUTWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);

							count2++;
						}
					}
				}
			}
			else if(CON->get_tool_angle()>0)//�c�[����x���܂��ɉ�]������B�i�s������+y
			{
				double theta=PI/180*CON->get_tool_angle();//��]����p�x
				
				for(int i=0;i<number;i++)
				{
					if(type1[i]==INWALL)
					{
						double speed=0;
						double u=0;
						double v=0;
						double uw=0;
						int AA=1;

						if(i<tool_number)
						{
							double r=sqrt(X3[i]*X3[i]+Y3[i]*Y3[i]);//�c�[���̒��S�����_�łȂ��Ƃ��͒���
							if(r>0.1*le && r<=shold_R+le && Z3[i]>0)
							{
								speed=r*w;
								u=speed*(-Y3[i]/r);
								v=speed*(X3[i]/r);
								if(Z[i]<Height && r>probe_R-0.3*le) uw=-pich*rps;//INWALL�̂����A�v���[�u�̑��ʂ̂�(��ʂ͏���)�ɉ��������x�ǉ�
								//�����܂ł̑��x�����́A�c�[���̉�]�Ƌ��ɉ�]����̂ō��W�Ɠ��l�ɉ�]������

								double u_t=u;
								double v_t=cos(theta)*v-sin(theta)*uw;
								double uw_t=sin(theta)*v+cos(theta)*uw;
			
								u=u_t;
								v=v_t;
								uw=uw_t;

								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								if(calc_type==plunge) uw-=U;//Z�����Ƀc�[�����ړ�(plange phase)
							}
							if(r<=shold_R+le && Z3[i]>0) AA=MOVE;
						}

						materialID=1;		
						writedata2(fq, count2, X[i],Y[i],Z[i], INWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
						count2++;
					}
				}
				//cout<<"inwall����"<<endl;

				for(int i=0;i<number;i++)
				{
					if(type1[i]==OUTWALL)
					{
						//cout<<"out="<<i<<endl;
						double speed=0;
						double u=0;
						double v=0;
						double uw=0;
						
						int AA=1;

						if(i<tool_number)
						{
							double r=sqrt(X3[i]*X3[i]+Y3[i]*Y3[i]);//�c�[���̒��S�����_�łȂ��Ƃ��͒���
						
							if(r>0.1*le && r<=shold_R+le && Z3[i]>0)
							{
								speed=r*w;
								u=speed*(-Y3[i]/r);
								v=speed*(X3[i]/r);
								uw=-pich*rps;

								double u_t=u;
								double v_t=cos(theta)*v-sin(theta)*uw;
								double uw_t=sin(theta)*v+cos(theta)*uw;
			
								u=u_t;
								v=v_t;
								uw=uw_t;

								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								else if(calc_type==plunge) uw-=U;//Z�����Ƀc�[�����ړ�(plange phase)
							}
							if(r<=shold_R+le && Z3[i]>0) AA=MOVE;
						}
							
						materialID=1;
						writedata2(fq, count2, X[i],Y[i],Z[i], OUTWALL,materialID,OFF,0, u,v,uw,0,wall_h,AA);
						count2++;
					}
						
				}
			}

		}
		else while(1) cout<<"model�����G���["<<endl;
	}


	else while(1) cout<<"model�G���["<<endl;
	fq.close();


	//gnuplot�p�o��
	ofstream fp("plot.dat");
	for(int i=0;i<number;i++) fp<<X[i]<<" "<<Y[i]<<" "<<Z[i]<<endl;
	fp.close();
	/////////////////////
	
	*particle_number=number;
	ofstream fn("particle_number.dat");
	fn<<number<<endl;
	fn.close();
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

void writedata2(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM)
{	
	//�t�@�C���o��
	fp<<id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<P<<"\t"<<h<<"\t"<<val<<"\t"<<type<<"\t"<<materialID<<"\t"<<surface<<"\t"<<toBEM<<"\t"<<endl;
}

/*void writedata(ofstream &fp, int number, double x, double y,double z, int type,double vx,double vy,double vz,double P,double h,int toFEM)
{	
	//fprintf( fp, "%5.10f\t", 0); �Ƃ����������̕��荞�݂͂��߁H
	
	double angle=0;//��]�p
	double angle2=0;
	double angle3=0;
	double angle_s=1;
	double anglar_u=0;//�p���x
	double anglar_u2=0;//�p���x
	double anglar_u3=0;//�p���x
	
	fp<<number<<"\t";
	fp<<x<<"\t";
	fp<<y<<"\t";
	fp<<z<<"\t";
	fp<<vx<<"\t";					//���xx����
	fp<<vy<<"\t";					//���xy����
	fp<<vz<<"\t";					//���xz����
	fp<<P<<"\t";					//����
	fp<<h<<"\t";					//�G���^���s�[
	fp<<angle<<"\t";					//��]�p
	fp<<type<<"\t";
	fp<<toFEM<<endl;
}*/


//���aR�̉~�̊O��
void set_circle_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R)
{
	//�Ώ̐����l�������A�~�����L�q���闱�q�͋����łȂ���΂Ȃ�Ȃ��B

	int N=calc_division_N_circle(2*PI*R,le);//�~���̕�����
	
	double L=2*PI*R/N;				//���q�ԋ���
	double theta=L/R;				//���q��z�u����p�x

	for(int n=0;n<N;n++)
	{
		X.push_back(R*cos(theta*n));
		Y.push_back(R*sin(theta*n));
		Z.push_back(0);
	}
	*number=*number+N;
}

//�~���������v�Z�֐�
int calc_division_N_circle(double dis,double le)
{
	//�Ώ̐����l�������A�~�����L�q���闱�q�͋����łȂ���΂Ȃ�Ȃ��B�����瑼�̕ӕ������Ƃ͈�������������
	//dis:�������鋗��(�~��)
	double temp_num=dis/le;		//�~�O���ɐݒu����w���́x���q���B�������O�������܂�le�Ŋ���؂��Ƃ͌���Ȃ�

	int N1=(int)(temp_num/2);
	N1*=2;							
	int N2=N1+2;					//temp_num��N1��N2�̊Ԃɂ���B������N1,N2�͋���

	double dif1=temp_num-N1;		//�eN�Ƃ̍�
	double dif2=N2-temp_num;
	int N=N1;						//������������
	if(dif2<dif1) N=N2;				//���̏���������N�Ƃ��č̗p����B

	return N;
}

//���aR�̉~����
void set_circle_in(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &type1,int *number,double le,double R,int edge_startID,int edge_lastID)
{
	//edge_startID����(edge_lastID-1)�܂ł̗��q���A�~�̊O�����\�����闱�q�ɊY������
	//vector�^�z��͎Q�Ɠn�����Ă���Bvector<double> *X�ł͂Ȃ�vector<double> &X�ł��邱�Ƃɒ��ӁB����Ŋe�z��͒ʏ�ʂ�Ɏd�l�\�B�A���[���Z�q������Ȃ�
	//�Q�Ɠn���łȂ��ʏ�̂�肩���ł��������ǁA���̏ꍇ�A�Ⴆ��a=X[5]�Ə����Ă����p�ł��Ȃ��Ba=(*X)[5]�ȂǂƂ��Ȃ���΂Ȃ�Ȃ�

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;		//�悭�g���W��

	int half_WX=(int)(R/le)+1;  //�~���\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R/(le*A))+1;  //�~���\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-le*0.5;				//���������߂̔��a��ݒ�
	
	//�����ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii*ii+jj*jj<R2*R2)
			{
				X.push_back(ii);
				Y.push_back(jj);
				Z.push_back(0);
				newN++;
			}
		}
	}

	//���q���͊w�ɂ��ʒu���œK��
	MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);

	*number=*number+newN;
}

//���aR�̉~����
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int edge_startID,int edge_lastID)
{
	//set_circle_in()�ƈႢ�A60�x�������v�Z���A�����6�º�߰���ĉ~���\������B���ԒZ�k�Ɣz�u�̋ϓ������ړI
	//edge_startID����(edge_lastID-1)�܂ł̗��q���A�~�̊O�����\�����闱�q�ɊY������
	//vector�^�z��͎Q�Ɠn�����Ă���Bvector<double> *X�ł͂Ȃ�vector<double> &X�ł��邱�Ƃɒ��ӁB����Ŋe�z��͒ʏ�ʂ�Ɏd�l�\�B�A���[���Z�q������Ȃ�
	//�Q�Ɠn���łȂ��ʏ�̂�肩���ł��������ǁA���̏ꍇ�A�Ⴆ��a=X[5]�Ə����Ă����p�ł��Ȃ��Ba=(*X)[5]�ȂǂƂ��Ȃ���΂Ȃ�Ȃ�

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;		//�悭�g���W��

	int half_WX=(int)(R/le)+1;  //�~���\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R/(le*A))+1;  //�~���\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-le*0.5;				//���������߂̔��a��ݒ�


	//////////////////////�~��6�ɂ킯�邽�߂̒���6�𐶐�

	double temp_R_num=R/le;			//���a�����ɐݒu����w���́x���q��
	int    R_num=(int)temp_R_num;	//�^�̗��q���@�Ƃ肠�������̗��q���̐����ԂƂ���B�����ŁAtemp_R_num>R_num���������Ă���B
	double difference=temp_R_num-R_num;	//���̐��Ɛ^�̐��̍�
	
	//���̐��Ɛ^�̐��̍���0.5�܂łȂ�A�^�̐���N�Ƃ���B0.5�ȏ�Ȃ�N+1�Ƃ���
	if(difference>0.5) R_num++;
	double L=R/R_num;					//���q�ԋ���

	X.push_back(0);						//���S���q��ǉ�
	Y.push_back(0);
	Z.push_back(0);
	newN++;

	for(int k=0;k<6;k++)				//6�̒�����loop
	{
		double theta=PI/3*k;			//�����̊p�x
		for(int n=1;n<R_num;n++)		//���S���q�ƍŊO�����q�͂������邩��A�����ł�loop�͂�����J�E���g���Ȃ�
		{
			double r=L*n;				//���S����̋���
			X.push_back(r*cos(theta));
			Y.push_back(r*sin(theta));
			Z.push_back(0);
			newN++;
		}
	}
	*number=*number+newN;
	newN=0;
	beforeN=*number;
	edge_lastID=beforeN;
	/////////////////����6�𐶐�����

	//���������ʒu �������ŏ���1�s�[�X�̂�
	for(int i=0;i<=half_WX;i++)
	{
		for(int j=1;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii*ii+jj*jj<R2*R2)
			{
				if(jj<sqrt(3.0)*ii-le)		//�ŏ��̃s�[�X�̎΂ߐ����Ⴂ�̈�ɐݒu�B�������������肬��͂܂����̂ŁA�ی���-le���Ă���
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(0);
					newN++;
				}
			}
		}
	}////////////////////////

	//���q���͊w�ɂ��ʒu���œK��
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	///�������q����������6�º�߰
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//��]����p�x
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*X[i]-sin(theta)*Y[i];//��]��̍��W
			double y=sin(theta)*X[i]+cos(theta)*Y[i];

			X.push_back(x);
			Y.push_back(y);
			Z.push_back(0);
			//newN++;
		}
	}///////////////////*/

	*number=*number+newN*6;//newN�͂ЂƂ̃s�[�X���̗��q����\���Ă��邩�炱���ł�6�{
}

//���aR�̋��쐬�֐�
void set_sphere(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag)
{
	//�܂��͔��������B���̂��߂ɂ͔����\�ʂ��쐬����K�v������B
	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��
	int beforeN2=*number;		//�֐��Ăяo�����̗��q���B���̊֐��̍Ō�܂ŋL�����Ă���

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	int half_WX=(int)(R/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WZ=(int)(R/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-0.5*le;				//���������߂̔��a��ݒ�

	///////////�����\��
	int Nt;						//���\�ʂ́A�ƕ����̕�����
	double Lt;					//���\�ʂ́A�ƕ����̕�������
	calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//���~�̕������͋����E��ǂ���ł��悢
	double d_theta=Lt/R;		//�ʂ̒�����Lt�ɂȂ�p�x

	for(int k=0;k<Nt;k++)//loop��k<Nt�ŏI��点��BNt�ɊY������Ƃ���͂��łɐݒu�ς�
	{
		double THETA=k*d_theta;	//��
		double r=R*sin(THETA);	//���̍����ɂ�����~�̔��a
		double round=2*PI*r;//���̍����ɂ�����~��

		int Nf=calc_division_N_circle(round,le);//���\�ʂ́A�ƕ����̕�����
		double Lf=round/Nf;						//���\�ʂ́A�ƕ����̕�������
		double d_fai=Lf/r;						//�ʂ̒�����Lf�ɂȂ�p�x
		
		for(int i=0;i<Nf;i++)
		{
			double fai=d_fai*i;
			if(Nt%2==0)
			{
				if(k%2!=0) fai+=0.5*d_fai;//Nt�������Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͊�ԖځB����Ċ�����炷
			}
			else
			{
				if(k%2==0) fai+=0.5*d_fai;//Nt����Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͋����ԖځB����Ċ�����炷
			}
			double x=r*cos(fai);
			double y=r*sin(fai);
			double z=R*cos(THETA);
			X.push_back(x);
			Y.push_back(y);
			Z.push_back(z);
			
			newN++;
		}
	}
	if(Nt%2!=0)//Nt����̂Ƃ��́A����ɗ��q���u����Ȃ���΂Ȃ�Ȃ��B���������loop�͂��ꂪ�s�\�B����Ă����Œǉ�
	{
		X.push_back(0);
		Y.push_back(0);
		Z.push_back(R);
		newN++;
	}
	//////////////////////////////////////

	*number=*number+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��
	
	//���������ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(ii*ii+jj*jj+kk*kk<R2*R2*0.6*0.6)//���a*0.7���̗��q�͈ʒu�Œ�
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;

	newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��
	///////////////////////*/

	

	//���������ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				//if(ii*ii+jj*jj+kk*kk<R2*R2)
				if(ii*ii+jj*jj+kk*kk<R2*R2 && ii*ii+jj*jj+kk*kk>=R2*R2*0.6*0.6)
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}///////////////////////*/
	

	//���q���͊w�ɂ��ʒu���œK��
	//cout<<"�ʒu�œK��(sphere)---------";
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-0.2*R; rigion[A_Z][1]=1.2*R;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;
	}
	//cout<<"����"<<endl;
	///�㔼�����������ֺ�߰����������
	if(flag==FULL)
	{
		newN=0;					//�V�����ǉ����闱�q��
		beforeN=*number;		//���̎��ɂ����闱�q��

		for(int k=0;k<beforeN;k++)
		{
			if(Z[k]>0.4*le)
			{
				X.push_back(X[k]);
				Y.push_back(Y[k]);
				Z.push_back(-Z[k]);
				newN++;
			}
		}
		*number=*number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//���������ق����Ƃ��ɁA�������㔼�����㉺���]������
	{
		for(int k=beforeN2;k<*number;k++) Z[k]*=-1;
	}
	
}

void set_sphere2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag,int *suf_num)
{
	//�܂��͔��������B���̂��߂ɂ͔����\�ʂ��쐬����K�v������B
	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��
	int beforeN2=*number;		//�֐��Ăяo�����̗��q���B���̊֐��̍Ō�܂ŋL�����Ă���

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	int half_WX=(int)(R/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WZ=(int)(R/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-0.5*le;				//���������߂̔��a��ݒ�

	///////////�����\��
	int Nt;						//���\�ʂ́A�ƕ����̕�����
	double Lt;					//���\�ʂ́A�ƕ����̕�������
	calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//���~�̕������͋����E��ǂ���ł��悢
	double d_theta=Lt/R;		//�ʂ̒�����Lt�ɂȂ�p�x

	for(int k=0;k<Nt;k++)//loop��k<Nt�ŏI��点��BNt�ɊY������Ƃ���͂��łɐݒu�ς�
	{
		double THETA=k*d_theta;	//��
		double r=R*sin(THETA);	//���̍����ɂ�����~�̔��a
		double round=2*PI*r;//���̍����ɂ�����~��

		int Nf=calc_division_N_circle(round,le);//���\�ʂ́A�ƕ����̕�����
		double Lf=round/Nf;						//���\�ʂ́A�ƕ����̕�������
		double d_fai=Lf/r;						//�ʂ̒�����Lf�ɂȂ�p�x
		
		for(int i=0;i<Nf;i++)
		{
			double fai=d_fai*i;
			if(Nt%2==0)
			{
				if(k%2!=0) fai+=0.5*d_fai;//Nt�������Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͊�ԖځB����Ċ�����炷
			}
			else
			{
				if(k%2==0) fai+=0.5*d_fai;//Nt����Ȃ�A�쐬�ς݂̉~�Ɛڂ���Ƃ��͋����ԖځB����Ċ�����炷
			}
			double x=r*cos(fai);
			double y=r*sin(fai);
			double z=R*cos(THETA);
			X.push_back(x);
			Y.push_back(y);
			Z.push_back(z);
			
			newN++;
		}
	}
	if(Nt%2!=0)//Nt����̂Ƃ��́A����ɗ��q���u����Ȃ���΂Ȃ�Ȃ��B���������loop�͂��ꂪ�s�\�B����Ă����Œǉ�
	{
		X.push_back(0);
		Y.push_back(0);
		Z.push_back(R);
		newN++;
	}
	//////////////////////////////////////

	*number=*number+newN;
	*suf_num=*suf_num+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��
	
	//���������ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(ii*ii+jj*jj+kk*kk<R2*R2*0.6*0.6)//���a*0.7���̗��q�͈ʒu�Œ�
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;

	newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��
	///////////////////////*/

	

	//���������ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				//if(ii*ii+jj*jj+kk*kk<R2*R2)
				if(ii*ii+jj*jj+kk*kk<R2*R2 && ii*ii+jj*jj+kk*kk>=R2*R2*0.6*0.6)
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}///////////////////////*/
	

	//���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-0.2*R; rigion[A_Z][1]=1.2*R;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;
	}

	///�㔼�����������ֺ�߰����������
	if(flag==FULL)
	{
		newN=0;					//�V�����ǉ����闱�q��
		beforeN=*number;		//���̎��ɂ����闱�q��

		for(int k=0;k<beforeN;k++)
		{
			if(Z[k]>0.4*le)
			{
				X.push_back(X[k]);
				Y.push_back(Y[k]);
				Z.push_back(-Z[k]);
				newN++;
			}
		}
		*number=*number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//���������ق����Ƃ��ɁA�������㔼�����㉺���]������
	{
		for(int k=beforeN2;k<*number;k++) Z[k]*=-1;
	}
	
}


//�����𕪊����邳���̍œK�ȕ������ƕ��������̎Z�o�֐�
void calc_N_and_L(double dis,double le,int *N,double *L)
{
	double temp_N=dis/le;			//���̕������Ble�Ŋ���؂ꂽ���Ԃ������ǁA�����������Ȃ��Ƃ�������
	int Ns=(int) temp_N;				//�^�̕�����
	double difference=temp_N-Ns;		//���Ɛ^�̍�
	if(difference>0.5) Ns++;
	*L=dis/Ns;			//���q�̋���
	*N=Ns;
}

//���q���͊w�֐�
void MD_2D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN)
{
	//���q���͊w�ɂ��newN�̗��q�̈ʒu���œK���@ID��BstartID����BendID�܂ł̂͋��E���q�Ȃ̂œ������Ȃ�

	double region[2][2];	//��͗̈�

	/////////////////////��͗̈�̌���
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	for(int i=BstartID;i<BendID;i++)
	{
		if(X[i]<region[A_X][0]) region[A_X][0]=X[i];
		else if(X[i]>region[A_X][1]) region[A_X][1]=X[i];

		if(Y[i]<region[A_Y][0]) region[A_Y][0]=Y[i];
		else if(Y[i]>region[A_Y][1]) region[A_Y][1]=Y[i];
	}
	for(int D=0;D<2;D++)
	{
		region[D][0]-=5*le;	//�����̈���L�߂ɂƂ�
		region[D][1]+=5*le;
	}//////////////////////////

	//�p�����[�^
	double k0=1;
	double r=1.5;
	double dt=0.001;
	
	//�͂�ax^3+bx^2+d�̎����̗p�B����[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]���Q��
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;

	vector<double> Fx(newN);	//�e���q�ɓ���X������
	vector<double> Fy(newN);	//�e���q�ɓ���Y������
	vector<double> Ax(newN,0);	//X���������x
	vector<double> Ay(newN,0);	//Y���������x
	vector<double> U(newN,0);	//X�������x
	vector<double> V(newN,0);	//Y�������x
	vector<double> visX(newN);	//X�����S���W��
	vector<double> visY(newN);	//Y�����S���W��

	//�v�Z�̍������̂��߂Ɋi�q���`�� ��͕���r*le�Ŋ���؂��Ƃ͌���Ȃ��̂ŁA�͂ݏo�����Ƃ���͐؂�̂āB�Ȃ̂Ŋe���Ƃ����̕����ɂ͗]�T��������
	double grid_width=le*((int)(r+1));								//�i�q�̕��Br���܂ސ���*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X�����̊i�q�̌�
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//�e�������q���܂ފi�q�ԍ�
	vector<int> *MESH=new vector<int>[plane_SIZE];				//�e���b�V���Ɋi�[����闱�qID�i�[

	for(int i=BstartID;i<BendID;i++)	//�܂��͋��E���q���i�q�Ɋi�[
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);	//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);	//Y�����ɉ��ڂ̊i�q��
		int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//���ɓ������q���i�[
	{
		int i=beforeN+k;
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
		index[k]=number;
	}//////////////////////////////////////////

	
	//�v�Z�J�n
	for(int t=0;t<100;t++)
	{
		if(t%10==0 &&t>0)//MESH����蒼��
		{
			//�܂���MESH����x�j�󂷂�B
			for(int n=0;n<plane_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<BendID;i++)	//�܂��͋��E���q���i�q�Ɋi�[
			{
				int xn=(int)((X[i]-region[A_X][0])/grid_width);	//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);	//Y�����ɉ��ڂ̊i�q��
				int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//���ɓ������q���i�[
			{
				int i=beforeN+k;
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
				int number=yn*grid_sizeX+xn;					//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
				index[k]=number;
			}
		}////////////

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0;					//������
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			double kx=0;						//X�����o�l�W��
			double ky=0;
			int G_id=index[k];				//�i�[����i�q�ԍ�
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=X[j]-X[i];
						double y=Y[j]-Y[i];
						double dis=sqrt(x*x+y*y);
						if(dis<r*le && dis!=0)			//����loop�͎������g���ʉ߂��邩��Adis!=0�͕K�v
						{
							double F=a*dis*dis*dis+b*dis*dis+d;
							Fx[k]-=F*x/dis;					//F�̒l�����̂Ƃ��͐˗͂Ȃ̂ŁA-=�ɂ���
							Fy[k]-=F*y/dis;
							double K=3*a*dis*dis+2*b*dis;//�o�l�W���@�͂̎��̔����ɑ���
							K=sqrt(K*K);					//���̒l���~�����B�����畉�̂Ƃ��ɔ����Đ��ɕϊ�
							kx+=K*x*x/(dis*dis);			//k���e�����ɕ��z�B�����ŁA��ɐ��̗ʂ����z�����悤��x*x/(dis*dis)�ƂȂ��Ă���
							ky+=K*y*y/(dis*dis);
						}
					}
				}
			}
			visX[k]=1.414*sqrt(kx);//���̂悤�Ɋe�������̔S���W�������߂�B�����u�������f���ɂ�鎩�����b�V�������vP6�Q�ƁB���������ʂ�1�Ƃ��Ă���B
			visY[k]=1.414*sqrt(ky);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
		}//�e���q�̉����x�����܂����B
		
		if(t==0)	//�ŏ��̽ï�ߎ���dt������
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//�ő�����x�����܂���
			dt=sqrt(0.02*le/MaxAccel);
		}

		for(int k=0;k<newN;k++)//���x�ƈʒu�̍X�V
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			X[i]+=dt*(U[k]+u)*0.5;
			Y[i]+=dt*(V[k]+v)*0.5;
		}

		//�ċߐڋ�����le�ȉ��̏ꍇ�͂�����C��
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			int G_id=index[k];				//�i�[����i�q�ԍ�
			double mindis=le;
			int J=k;						//�ŋߐڋ����̑��藱�q
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=X[j]-X[i];
						double y=Y[j]-Y[i];
						double dis=sqrt(x*x+y*y);
						if(dis<mindis && i!=j)
						{
							mindis=dis;
							J=j;
						}
					}
				}
			}
			if(J!=i && J<beforeN)//le���ߐڂ��Ă��鑊�肪���E���q�Ȃ�
			{
				double L=le-mindis;//�J���ׂ�����
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;			
			}
			else if(J!=i && J>=beforeN)//le���ߐڂ��Ă��鑊�肪�������q�Ȃ�
			{
				double L=0.5*(le-mindis);//�J���ׂ�����
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;
				X[J]+=dX/mindis*L;
				Y[J]+=dY/mindis*L;
			}
		}//////////*/
	}/////MD�I��

	delete [] index;
	delete [] MESH;
}

//���q���͊w�֐�
void MD_3D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN,double r,double region[3][2])
{
	//���q���͊w�ɂ��newN�̗��q�̈ʒu���œK���@ID��BstartID����BendID�܂ł̂͋��E���q�Ȃ̂œ������Ȃ�
	double k0=1;
	double dt=0.001;
	
	//�͂�ax^3+bx^2+d�̎����̗p�B����[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]���Q��
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;

	//cout<<"F="<<a*le*le*le+b*le*le+d<<" "<<a*1.5*le*1.5*le*1.5*le+b*1.5*le*1.5*le+d<<endl;

	vector<double> Fx(newN);	//�e���q�ɓ���X������
	vector<double> Fy(newN);	//�e���q�ɓ���Y������
	vector<double> Fz(newN);	//�e���q�ɓ���Z������
	vector<double> Ax(newN,0);	//X���������x
	vector<double> Ay(newN,0);	//Y���������x
	vector<double> Az(newN,0);	//Z���������x
	vector<double> U(newN,0);	//X�������x
	vector<double> V(newN,0);	//Y�������x
	vector<double> W(newN,0);	//Z�������x
	vector<double> visX(newN);	//X�����S���W��
	vector<double> visY(newN);	//Y�����S���W��
	vector<double> visZ(newN);	//Y�����S���W��
	vector<double> KX(newN);	//X�����o�l�W��
	vector<double> KY(newN);	//Y�����o�l�W��
	vector<double> KZ(newN);	//Y�����o�l�W��

	//�v�Z�̍������̂��߂Ɋi�q���`�� ��͕���r*le�Ŋ���؂��Ƃ͌���Ȃ��̂ŁA�͂ݏo�����Ƃ���͐؂�̂āB�Ȃ̂Ŋe���Ƃ����̕����ɂ͗]�T��������
	double grid_width=le*((int)(r+1));								//�i�q�̕��Br���܂ސ���*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X�����̊i�q�̌�
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//�e�������q���܂ފi�q�ԍ�
	vector<int> *MESH=new vector<int>[grid_SIZE];				//�e���b�V���Ɋi�[����闱�qID�i�[

	for(int i=BstartID;i<BendID;i++)	//�܂��͋��E���q���i�q�Ɋi�[
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//���ɓ������q���i�[
	{
		int i=beforeN+k;
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
		MESH[number].push_back(i);
		index[k]=number;
	}

	//�v�Z�J�n
	for(int t=0;t<100;t++)
	{
		if(t%10==0 && t>0)
		{
			//MESH����x�j�󂷂�B
			for(int n=0;n<grid_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<BendID;i++)	//�܂��͋��E���q���i�q�Ɋi�[
			{
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
				int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//���ɓ������q���i�[
			{
				int i=beforeN+k;
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
				int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
				MESH[number].push_back(i);
				index[k]=number;
			}
		}

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0, Fz[k]=0;			//������
			KX[k]=0;KY[k]=0; KZ[k]=0;			//�o�l�W��
		}

		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			int G_id=index[k];				//�i�[����i�q�ԍ�
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							if(j>=beforeN && j>i)	//�����̈���ł���i���傫�Ȕԍ��Ȃ�
							{
								int J=j-beforeN;	//newN���ł̔ԍ�
								double x=X[j]-X[i];
								double y=Y[j]-Y[i];
								double z=Z[j]-Z[i];
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le)			//����loop�͎������g���ʉ߂��邩��Adis!=0�͕K�v
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//F�̒l�����̂Ƃ��͐˗͂Ȃ̂ŁA-=�ɂ���
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									Fx[J]+=F*x/dis;					//���藱�q�̗͂������Ōv�Z�B�����͔��]������
									Fy[J]+=F*y/dis;
									Fz[J]+=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//�o�l�W���@�͂̎��̔����ɑ���
									K=sqrt(K*K);					//���̒l���~�����B�����畉�̂Ƃ��ɔ����Đ��ɕϊ�
									KX[k]+=K*x*x/(dis*dis);			//k���e�����ɕ��z�B�����ŁA��ɐ��̗ʂ����z�����悤��x*x/(dis*dis)�ƂȂ��Ă���
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
									KX[J]+=K*x*x/(dis*dis);			//k�𑊎藱�q�ɂ����z
									KY[J]+=K*y*y/(dis*dis);
									KZ[J]+=K*z*z/(dis*dis);
								}
							}
							if(j<BendID && j>=BstartID)
							{
								double x=X[j]-X[i];
								double y=Y[j]-Y[i];
								double z=Z[j]-Z[i];
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le && dis>0)			//����loop�͎������g�͒ʉ߂��Ȃ��Adis!=0�͕s�v
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//F�̒l�����̂Ƃ��͐˗͂Ȃ̂ŁA-=�ɂ���
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//�o�l�W���@�͂̎��̔����ɑ���
									K=sqrt(K*K);					//���̒l���~�����B�����畉�̂Ƃ��ɔ����Đ��ɕϊ�
									KX[k]+=K*x*x/(dis*dis);			//k���e�����ɕ��z�B�����ŁA��ɐ��̗ʂ����z�����悤��x*x/(dis*dis)�ƂȂ��Ă���
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
								}
							}
						}
					}
				}
			}
			//visX[k]=1.414*sqrt(KX[k]);//���̂悤�Ɋe�������̔S���W�������߂�B�����u�������f���ɂ�鎩�����b�V�������vP6�Q�ƁB���������ʂ�1�Ƃ��Ă���B
			//visY[k]=1.414*sqrt(KY[k]);
			//visZ[k]=1.414*sqrt(KZ[k]);
			visX[k]=1.414*sqrt(KX[k]);//���̂悤�Ɋe�������̔S���W�������߂�B�����u�������f���ɂ�鎩�����b�V�������vP6�Q�ƁB���������ʂ�1�Ƃ��Ă���B
			visY[k]=1.414*sqrt(KY[k]);
			visZ[k]=1.414*sqrt(KZ[k]);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
			Az[k]=(Fz[k]-visZ[k]*W[k]);
		}//�e���q�̉����x�����܂����B
		
		if(t==0)	//�ŏ��̽ï�ߎ���dt������
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k]+Az[k]*Az[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//�ő�����x�����܂���
			if(MaxAccel!=0)
			{
				dt=sqrt(0.02*le/MaxAccel);
			}
		}

		for(int k=0;k<newN;k++)//���x�ƈʒu�̍X�V
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			double w=W[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			W[k]+=dt*Az[k];
			X[i]+=dt*(U[k]+u)*0.5;
			Y[i]+=dt*(V[k]+v)*0.5;
			Z[i]+=dt*(W[k]+w)*0.5;
		}

		//�ċߐڋ�����le�ȉ��̏ꍇ�͂�����C��
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//�Ή����闱�q�ԍ�
			int G_id=index[k];				//�i�[����i�q�ԍ�
			double mindis=le;
			int J=k;						//�ŋߐڋ����̑��藱�q
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							double x=X[j]-X[i];
							double y=Y[j]-Y[i];
							double z=Z[j]-Z[i];
							double dis=sqrt(x*x+y*y+z*z);
							if(dis<mindis && i!=j)
							{
								mindis=dis;
								J=j;
							}
						}
					}
				}
			}
			if(J!=i && J<beforeN)//le���ߐڂ��Ă��鑊�肪���E���q�Ȃ�
			{
				double L=le-mindis;//�J���ׂ�����
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				double dZ=Z[J]-Z[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;	
				Z[i]-=dZ/mindis*L;	
			}
			else if(J!=i && J>=beforeN)//le���ߐڂ��Ă��鑊�肪�������q�Ȃ�
			{
				double L=0.5*(le-mindis);//�J���ׂ�����
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				double dZ=Z[J]-Z[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;
				Z[i]-=dZ/mindis*L;
				X[J]+=dX/mindis*L;
				Y[J]+=dY/mindis*L;
				Z[J]+=dZ/mindis*L;
			}
		}//////////*/
	}/////MD�I��

	delete [] index;
	delete [] MESH;
}


//�����`�̕Ӎ쐬�֐�
void set_rectangular_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V)
{
	//������������Len_H,��������Len_V�̒����`�̕ӂ��쐬����
	
	int newN=0;					//���̊֐��Œǉ�����闱�q��
	
	////////////////////�܂����������쐬

	int Nh;				//�����ӂ̕�����
	double dL_H;		//�����ӂ̕�������
	calc_N_and_L(Len_H,le,&Nh,&dL_H);

	for(int n=0;n<=Nh;n++)//��� �����ł�Y,Z���W�̓[���Ƃ��Ă����B
	{
		X.push_back(-Len_H*0.5+dL_H*n);
		Y.push_back(-Len_V*0.5);
		Z.push_back(0);
		newN++;
	}
	for(int n=0;n<=Nh;n++)//���
	{
		X.push_back(-Len_H*0.5+dL_H*n);
		Y.push_back(Len_V*0.5);
		Z.push_back(0);
		newN++;
	}////////////////////////////////////////////

	////////////////////���ɐ��������쐬

	int Nv;				//�����ӂ̕�����
	double dL_V;		//�����ӂ̕�������
	calc_N_and_L(Len_V,le,&Nv,&dL_V);

	for(int n=1;n<Nv;n++)//���� �����ł�X���W�̓[���Ƃ��Ă����B�܂�n=0�ɊY������_�͂��łɐݒu�ς�
	{
		X.push_back(-Len_H*0.5);
		Y.push_back(-Len_V*0.5+n*dL_V);
		Z.push_back(0);
		newN++;
	}
	for(int n=1;n<Nv;n++)//�E�� 
	{
		X.push_back(Len_H*0.5);
		Y.push_back(-Len_V*0.5+n*dL_V);
		Z.push_back(0);
		newN++;
	}////////////////////////////////////////////

	*number=*number+newN;
}

//�����`�쐬�֐�
void set_rectangular_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V,int edge_startID,int edge_lastID)
{
	//������������Len_H,��������Len_V�̒����`�̕ӂ��쐬����
	
	int newN=0;					//���̊֐��Œǉ�����闱�q��
	int beforeN=*number;			//���̊֐��Ăяo�����ɂ��łɑ��݂��Ă��闱�q��

	double A=sqrt(3.0)/2;		//�悭�g���W��

	int half_WX=(int)(Len_H/le)+1;  //�����`���\���܂ޕ�
	int half_WY=(int)(Len_V/(le*A))+1;
	double gap=0.4*le;				//�ӂ��肬��ɓ������q��z�u���Ȃ��悤�A���Ԃ�݂���
	
	//�����ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WX;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii>-Len_H*0.5+gap && ii<Len_H*0.5-gap)
			{
				if(jj>-Len_V*0.5+gap && jj<Len_V*0.5-gap)
				{
					if(ii*ii<Len_H*0.25*Len_H*0.25 && jj*jj<Len_V*0.25*Len_V*0.25)//�\�������͔z�u�����ߑł�
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(0);
						newN++;
					}
				}
			}
		}
	}
	*number=*number+newN;
	newN=0;					
	beforeN=*number;
	edge_lastID=*number;//���EID��ύX�@�����܂łɐݒu���ꂽ���q�����E���q�ɂȂ�B

	//�����ʒu
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WX;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii>-Len_H*0.5+gap && ii<Len_H*0.5-gap)
			{
				if(jj>-Len_V*0.5+gap && jj<Len_V*0.5-gap)
				{
					if(ii*ii>=Len_H*0.25*Len_H*0.25 || jj*jj>=Len_V*0.25*Len_V*0.25)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(0);
						newN++;
					}
				}
			}
		}
	}

	//���q���͊w�ɂ��ʒu���œK��
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	*number=*number+newN;
}

//�~���\�ʍ쐬�֐�
void set_cylinder_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int circle_start_id,int circle_end_id,int top_flag)
{
	//���aR,����height�̉~���̖ʂ��쐬����B���������̊֐��Ăяo�����ɂ����āA���łɉ��ʂ̉~(Z=0)�͍쐬�ς݂Ƃ���
	//top_flag=ON�Ȃ�~����ʂ��쐬����BOFF�Ȃ炵�Ȃ����A���ʂ����͍쐬����B
	int beforeN=*number;
	int newN=0;

	int Nv;				//�����̕�����
	double dL_V;		//�����̕�������
	double A=sqrt(3.0)/2;		//�悭�g���W��
	calc_N_and_L(height,le*A,&Nv,&dL_V);

	int Nr=calc_division_N_circle(2*PI*R,le);//�~���̕�����
	double Lr=2*PI*R/Nr;				//�~����������

	double gap=0.4*le;				//�ӂ��肬��ɓ������q��z�u���Ȃ��悤�A���Ԃ�݂���

	///////////////////////////////////����
	for(int i=0;i<Nr;i++)
	{
		for(int j=1;j<Nv;j++)//j=0,j=Nv�͉��ʁA��ʂɊY������̂ł����ł͂ʂ���
		{
			double jj=j*dL_V;
			double ii=i*Lr;
			if(j%2!=0) ii+=0.5*Lr;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii<2*PI*R-gap)
			{
				if(jj<height-gap)	
				{
					double theta=2*PI*(ii/(2*PI*R));
					X.push_back(R*cos(theta));
					Y.push_back(R*sin(theta));
					Z.push_back(jj);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	////////////////////////


	//��ʍ쐬�i���ʂ̃R�s�[�j������Nv����Ȃ��ʂ͉��ʂƔ��i�q����Ȃ���΂Ȃ�Ȃ�
	beforeN=*number;
	newN=0;
	if(Nv%2==0)	//�����Ȃ炻�̂܂ܺ�߰
	{
		if(top_flag==ON)
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				X.push_back(X[i]);
				Y.push_back(Y[i]);
				Z.push_back(height);
				newN++;
			}
		}
		if(top_flag==HALF)//���Ǖ��������͂Ȃ�
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R && r<R+4*le-0.1*le)//�O���̂ݍ쐬
				{
				X.push_back(X[i]);
				Y.push_back(Y[i]);
				Z.push_back(height);
				newN++;
				}
			}
		}
		else if(top_flag==OFF)//��ʂ��K�v�Ȃ��Ȃ�
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R-0.1*le)//�O���̂ݍ쐬
				{
					X.push_back(X[i]);
					Y.push_back(Y[i]);
					Z.push_back(height);
					newN++;
				}
			}
		}
	}
	else
	{
		double d_theta=0.5*Lr/R;//���̔����p�x������]������B
		if(top_flag==ON)
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//��]��̍��W
				double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
				X.push_back(X2);
				Y.push_back(Y2);
				Z.push_back(height);
				newN++;
			}
		}
		else if(top_flag==HALF)//��ʂ��K�v�Ȃ��Ȃ�
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				//if(r>R-0.1*le)//�O���̂ݍ쐬
				if(r>R && r<R+4*le-0.1*le)//�O���̂ݍ쐬
				{
					double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//��]��̍��W
					double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
					X.push_back(X2);
					Y.push_back(Y2);
					Z.push_back(height);
					newN++;
				}
			}
		}
		else if(top_flag==OFF)//��ʂ��K�v�Ȃ��Ȃ�
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R-0.1*le)//�O���̂ݍ쐬
				{
					double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//��]��̍��W
					double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
					X.push_back(X2);
					Y.push_back(Y2);
					Z.push_back(height);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	
	/////////////////////////////

}

//�~���\�ʍ쐬�֐�
void set_circular_cone_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R2,double height,int circle_start_id,int circle_end_id,int top_flag)
{
	//�㕔���aR,�������aR2,����height�̉~���̖ʂ��쐬����B���������̊֐��Ăяo�����ɂ����āA���łɉ��ʂ̉~(Z=0)�͍쐬�ς݂Ƃ���
	//top_flag=ON�Ȃ�~����ʂ��쐬����BOFF�Ȃ炵�Ȃ����A���ʂ����͍쐬����B
	int beforeN=*number;
	int newN=0;

	int Nv;				//�����̕�����
	double dL_V;		//�����̕�������
	double A=sqrt(3.0)/2;		//�悭�g���W��
	calc_N_and_L(height,le*A,&Nv,&dL_V);

	

	double gap=0.4*le;				//�ӂ��肬��ɓ������q��z�u���Ȃ��悤�A���Ԃ�݂���

	///////////////////////////////////����
	
	for(int j=1;j<=Nv;j++)//j=0,j=Nv�͉��ʁA��ʂɊY������̂ł����ł͂ʂ���
	{
		double jj=j*dL_V;
		double r=jj*(R-R2)/height + R2;
		int Nr=calc_division_N_circle(2*PI*r,le);//�~���̕�����
		double Lr=2*PI*r/Nr;				//�~����������

		for(int i=0;i<Nr;i++)
		{
			
			double ii=i*Lr;
			if(j%2!=0) ii+=0.5*Lr;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii<2*PI*R-gap)
			{
				if(jj<height-gap)	
				{
					
					double theta=2*PI*(ii/(2*PI*r));
					X.push_back(r*cos(theta));
					Y.push_back(r*sin(theta));
					Z.push_back(jj);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	////////////////////////


	//��ʍ쐬
	beforeN=*number;
	newN=0;
	
		//if(top_flag==ON)
		//{
		//	for(int i=circle_start_id;i<circle_end_id;i++)
		//	{
		//		X.push_back(X[i]);
		//		Y.push_back(Y[i]);
		//		Z.push_back(height);
		//		newN++;
		//	}
		//}
		//if(top_flag==HALF)//���Ǖ��������͂Ȃ�
		//{
		//	for(int i=circle_start_id;i<circle_end_id;i++)
		//	{
		//		double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
		//		if(r>R && r<R+4*le-0.1*le)//�O���̂ݍ쐬
		//		{
		//		X.push_back(X[i]);
		//		Y.push_back(Y[i]);
		//		Z.push_back(height);
		//		newN++;
		//		}
		//	}
		//}
		//else if(top_flag==OFF)//��ʂ��K�v�Ȃ��Ȃ�
		//{
		//	for(int i=circle_start_id;i<circle_end_id;i++)
		//	{
		//		double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
		//		if(r>R-0.1*le)//�O���̂ݍ쐬
		//		{
		//			X.push_back(X[i]);
		//			Y.push_back(Y[i]);
		//			Z.push_back(height);
		//			newN++;
		//		}
		//	}
		//}
	
	
	*number=*number+newN;
	
	/////////////////////////////

}

//�~�������ݒu�֐�
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int flag)
{
	//���aR,����height�̉~���������쐬����B���̊֐��Ăяo������0<=i<number�̗��q�ŉ~���\�ʂ��`������Ă���Ƃ���B
	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	int WX=(int)(R/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WY=(int)(R/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WZ=(int)(height/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R-0.3*le;				//���������߂̔��a��ݒ�
	double height2=height-0.3*le;				//���������߂̍�����ݒ�

	//���������ʒu
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(ii*ii+jj*jj<R2*R2)
				{
					if(kk<height2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;

}

void set_crucible_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R_out,double height,double height_out,int flag,int fluid_number)
{
	//���aR,����height�̉~���������쐬����B���̊֐��Ăяo������0<=i<number�̗��q�ŉ~���\�ʂ��`������Ă���Ƃ���B
	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	int WX=(int)(R_out/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WY=(int)(R_out/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WZ=(int)(height/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R+0.3*le;				//���������߂̔��a��ݒ�
	double R2_out=R_out-0.3*le;				//���������߂̔��a��ݒ�
	double height2=height-0.3*le;				//���������߂̍�����ݒ�

	//���������ʒu
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=-WZ;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(kk>=0 && kk<height2)//���ǉ~����
				{
					if(ii*ii+jj*jj<R2_out*R2_out &&ii*ii+jj*jj>R2*R2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<0 && kk>-R2_out)//���ǉ~����
				{
					if(ii*ii+jj*jj+kk*kk<R2_out*R2_out &&ii*ii+jj*jj+kk*kk>R2*R2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				
			}
		}
	}///////////////////////*/

	///���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R_out; rigion[A_X][1]=1.2*R_out;
	rigion[A_Y][0]=-1.2*R_out; rigion[A_Y][1]=1.2*R_out;
	rigion[A_Z][0]=-1.2*R_out;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,fluid_number,beforeN,beforeN,newN,r,rigion);
	
	*number=*number+newN;

}

//�h�[�i�c�쐬
void set_doughnut2D(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_big,double R_smal,int edge_startID,int edge_lastID)
{
	//���łɍ쐬���ꂽ�~����p���āA������傫�Ȕ��a�̃h�[�i�c���쐬����B
	//edge_startID����(edge_lastID-1)�܂ł̗��q���A���łɍ쐬����Ă���~�̊O�����\�����闱�q�ɊY������
	//vector�^�z��͎Q�Ɠn�����Ă���Bvector<double> *X�ł͂Ȃ�vector<double> &X�ł��邱�Ƃɒ��ӁB����Ŋe�z��͒ʏ�ʂ�Ɏd�l�\�B�A���[���Z�q������Ȃ�
	//�Q�Ɠn���łȂ��ʏ�̂�肩���ł��������ǁA���̏ꍇ�A�Ⴆ��a=X[5]�Ə����Ă����p�ł��Ȃ��Ba=(*X)[5]�ȂǂƂ��Ȃ���΂Ȃ�Ȃ�

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;		//�悭�g���W��

	int half_WX=(int)(R_big/le)+1;  //�~���\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int half_WY=(int)(R_big/(le*A))+1;  //�~���\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double R2=R_big-le*0.4;				//���������߂̔��a��ݒ�
	double R1=R_smal+le*0.4;				//�����傫�߂̔��a��ݒ�


	//�����ЂƂ̉~�����쐬 ���q���͓����ő������邱�Ƃɒ���
	set_circle_edge(X,Y,Z,number,le,R_big);

	edge_lastID=*number;//MD2D()�ɂ����鋫�E���q��edge_startID<=i<edge_lastID
	beforeN=*number;

	//���������ʒu �������ŏ���1�s�[�X�̂�
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
			if(ii*ii+jj*jj<R2*R2 && ii*ii+jj*jj>R1*R1)
			{
				X.push_back(ii);
				Y.push_back(jj);
				Z.push_back(0);
				newN++;
			}
		}
	}////////////////////////

	//���q���͊w�ɂ��ʒu���œK��
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	*number=*number+newN;
}

//FSW�v���[�u�������q�Z�b�g�֐�
void set_hat_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID)
{
	//�X�q�^�̕����������쐬����B�Ⴆ��FSW�̃c�[���`��B
	
	//�X�q�̓��ɊY�����锼�a��R_smal,�΂ɊY�����锼�a��R_big,�΂̕���H_flange,���̕���H_hat
	//�����ō쐬����X�q�^�̎p���́AFSW�̃c�[���Ɠ����ŁA�������ɂ��Ă΂���B
	//�v���[�u��ӂ�Z=0�Ƃ���

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	double height=H_hat+H_flange;

	int WX=(int)(R_big/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WY=(int)(R_big/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WZ=(int)(height/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double gap=0.4*le;					//����
	double R_big2=R_big-gap;				//���������߂̔��a��ݒ�
	double R_smal2=R_smal-gap;				//���������߂̔��a��ݒ�
	double height2=height-0.3*le;				//���������߂̍�����ݒ�
	

	//���������ʒu
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(kk<=H_hat+gap)
				{
					if(ii*ii+jj*jj<R_smal2*R_smal2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<height-gap)
				{
					if(ii*ii+jj*jj<R_big2*R_big2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R_big; rigion[A_X][1]=1.2*R_big;
	rigion[A_Y][0]=-1.2*R_big; rigion[A_Y][1]=1.2*R_big;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,bound_startID,bound_endID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//FSW�v���[�u�������q�Z�b�g�֐�
void set_hat_in_2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_mid,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID)
{
	//�X�q�^�̕����������쐬����B�Ⴆ��FSW�̃c�[���`��B
	
	//�X�q�̂Ă��؂�ɊY�����锼�a��R_smal,���Ԃ镔���ɊY�����锼�a��R_mid,�΂ɊY�����锼�a��R_big,�΂̕���H_flange,���̕���H_hat
	//�����ō쐬����X�q�^�̎p���́AFSW�̃c�[���Ɠ����ŁA�������ɂ��Ă΂���B
	//�v���[�u��ӂ�Z=0�Ƃ���

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��

	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��
	double height=H_hat+H_flange;

	int WX=(int)(R_big/le)+1;		 //�����\���܂ގl�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WY=(int)(R_big/(le*A))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	int WZ=(int)(height/(le*B))+1;  //�����\���܂ސ��l�p�`��z�肷��B���̎l�p�`�̕�*0.5
	double gap=0.4*le;					//����
	double R_big2=R_big-gap;				//���������߂̔��a��ݒ�
	double R_mid2=R_mid-gap;				//���������߂̔��a��ݒ�
	double R_smal2=R_smal-gap;				//���������߂̔��a��ݒ�
	double height2=height-0.3*le;				//���������߂̍�����ݒ�
	

	//���������ʒu
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(kk<=H_hat+gap)
				{
					double r=kk*(R_mid2-R_smal2)/H_hat + R_smal2;
					if(ii*ii+jj*jj<r*r)
					//if(ii*ii+jj*jj<R_smal2*R_smal2) double r=jj*(R-R2)/height + R2;
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<height-gap)
				{
					if(ii*ii+jj*jj<R_big2*R_big2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-1.2*R_big; rigion[A_X][1]=1.2*R_big;
	rigion[A_Y][0]=-1.2*R_big; rigion[A_Y][1]=1.2*R_big;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,bound_startID,bound_endID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//���쐬�֐�
void set_box(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &surface,int *number,double le,double Width,double Height,double Depth)
{
	//��(Width)�~����(Height)�~���s��(Depth)�̔����쐬����B�f�J���g���W�̌��_�ɑ΂��AX�������ɉ����AY�������ɉ��s���AZ�������ɍ���

	int BOX_startID=*number;	//���̊֐��J�n���̗��q��
	int beforeN=*number;
	int newN=0;
	double A=sqrt(3.0)/2;		//�悭�g���W��

	vector<double> X2;					//�g���񂵗p
	vector<double> Y2;
	vector<double> Z2;
	int number2;						//�g���܂킵�悤���q��

	int Ndepth;					//���s�������̕�����
	int Nwidth;					//���s�������̕�����
	int Nheight;				//���s�������̕�����
	double dL_depth;			//�s�������̕�������
	double dL_width;			//�s�������̕�������
	double dL_height;			//�s�������̕�������

	double gap=0.4*le;
	
	calc_N_and_L(Depth,le,&Ndepth,&dL_depth);//�e�����̕������Ƃ��̒��������܂�
	calc_N_and_L(Width,le,&Nwidth,&dL_width);
	calc_N_and_L(Height,le,&Nheight,&dL_height);

	//XY����(�㉺��)�쐬//////�����`�쐬�֐����g���̂ł͂Ȃ��A�����ł�����ƍ쐬����B���R�́A�ЂƂ̒����`�̕ӂ�ʂ̒����`���g�p���邩��

	set_rectangular(X,Y,Z,number,le,Width,Depth);		//���W�̌��_�͒����`�̒��S�B����Ă��Ƃňړ����邱�ƁB
	for(int i=beforeN;i<*number;i++)					//�d�S�ړ�
	{
		X[i]+=Width*0.5;
		Y[i]+=Depth*0.5;
	}

	for(int i=beforeN;i<*number;i++)					//��ʂɺ�߰
	{
		X.push_back(X[i]);
		Y.push_back(Y[i]);
		Z.push_back(Height);
		newN++;
	}
	*number=*number+newN;
	/////////////////////////////////////

	//XZ����(���ʁE�w��)�쐬
	number2=0;
	set_rectangular(X2,Y2,Z2,&number2,le,Width,Height);		//���W�̌��_�͒����`�̒��S�B����Ă��Ƃňړ����邱�ƁB
	for(int i=0;i<number2;i++)								//�d�S�ړ�
	{
		X2[i]+=Width*0.5;
		Y2[i]+=Height*0.5;		//set_rectangular()��2D�p�Ȃ̂ŁAZ�͒l���[���ɒ���
	}
		//���q�ǉ�
	newN=0;
	beforeN=*number;
	for(int i=0;i<number2;i++)								//�����ȍ��W�Ɉړ����A�����z��ɺ�߰
	{
		if(Y2[i]>gap && Y2[i]<Height-gap)					//Y2=0,Y2=Height�̗��q�͂��łɍ쐬�ς݂Ȃ̂ŏȂ�
		{
			X.push_back(X2[i]);
			Y.push_back(0.0);			//Y=0���Ȃ킿����
			Z.push_back(Y2[i]);			//Y2����l�����炤���Ƃɒ���
			newN++;
		}
	}
	*number=*number+newN;				//���ʗ��q�z�u����

	newN=0;
	for(int i=beforeN;i<*number;i++)
	{
		X.push_back(X[i]);
		Y.push_back(Depth);			//Y=Depth���Ȃ킿�w��
		Z.push_back(Z[i]);
		newN++;
	}
	*number=*number+newN;				//�w�ʗ��q�z�u����
	////////////////////////////////////////////////////////////


	//YZ����(����)�쐬
	newN=0;
	size_t vector_size=X2.size();
	for(int k=0;k<vector_size;k++)//X2,Y2,Z2�������������
	{
		X2.pop_back();
		Y2.pop_back();
		Z2.pop_back();
	}

	number2=0;
	set_rectangular(X2,Y2,Z2,&number2,le,Depth,Height);		//���W�̌��_�͒����`�̒��S�B����Ă��Ƃňړ����邱�ƁB
	for(int i=0;i<number2;i++)								//�d�S�ړ�
	{
		X2[i]+=Depth*0.5;
		Y2[i]+=Height*0.5;									//set_rectangular()��2D�p�Ȃ̂ŁAZ�͒l���[���ɒ���
	}
		//���q�ǉ�
	newN=0;
	beforeN=*number;
	for(int i=0;i<number2;i++)								//�����ȍ��W�Ɉړ����A�����z��ɺ�߰
	{
		if(X2[i]>gap && X2[i]<Depth-gap)					//�ӗ��q�͂��ł�4�ӂƂ��쐬���݂Ȃ̂ŏȂ�
		{
			if(Y2[i]>gap && Y2[i]<Height-gap)				
			{
				X.push_back(0.0);			//Z=0���Ȃ킿����
				Y.push_back(X2[i]);			//X2����l�����炤���Ƃɒ���
				Z.push_back(Y2[i]);			//Y2����l�����炤���Ƃɒ���
				newN++;
			}
		}
	}
	*number=*number+newN;				//���ʗ��q�z�u����

	newN=0;
	for(int i=beforeN;i<*number;i++)
	{
		X.push_back(Width);
		Y.push_back(Y[i]);			//Y=Depth���Ȃ킿�w��
		Z.push_back(Z[i]);
		newN++;
	}
	*number=*number+newN;				//�w�ʗ��q�z�u����
	////////////////////////////////////////////////////////////
	
	for(int i= BOX_startID;i<*number;i++)  surface.push_back(ON);//�����܂ł͕\��
	beforeN=*number;

	set_box_in(X,Y,Z,number,le,Width,Height,Depth,BOX_startID,*number);
	
	for(int i=beforeN;i<*number;i++) surface.push_back(OFF);//����
}

//�����`�쐬�֐�
void set_rectangular(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height)
{
	int BO_startID=*number;
	set_rectangular_edge(X,Y,Z,number,le,Width,Height);	//���W�̌��_�͒����`�̒��S�B����Ă��Ƃňړ����邱�ƁB
	int BO_lastID=*number;
	set_rectangular_in(X,Y,Z,number,le,Width,Height,BO_startID,BO_lastID);
}

//BOX���쐬�֐�
void set_box_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height,double Depth,int BO_startID,int BO_lastID)
{
	//��(Width)�~����(Height)�~���s��(Depth)�̔����쐬����B�f�J���g���W�̌��_�ɑ΂��AX�������ɉ����AY�������ɉ��s���AZ�������ɍ���

	int newN=0;					//�̊֐��ŐV�����ǉ����闱�q��
	int beforeN=*number;		//���̊֐��Ăяo�����ɂ����闱�q��
	
	double A=sqrt(3.0)/2;				//�悭�g���W��
	double B=sqrt(2.0/3);						////�悭�g���W��

	int WX=(int)(Width/le)+1;		
	int WY=(int)(Depth/(le*A))+1; 
	int WZ=(int)(Height/(le*B))+1;  
	double gap=0.4*le;					//����
	
	//�����Œ菉���ʒu
	for(int i=1;i<=WX;i++)
	{
		for(int j=1;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(kk<Height*0.75 && jj<Depth*0.75 && ii<Width*0.75)
				{
					if(kk>Height*0.25 && jj>Depth*0.25 && ii>Width*0.25)
					{
						X.push_back(ii);//�\�������͌��ߑł�
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}
	*number=*number+newN;
	newN=0;					//���̊֐��ŐV�����ǉ����闱�q��
	beforeN=*number;
	///////////////////////*/

	//�������������ʒu
	for(int i=1;i<=WX;i++)
	{
		for(int j=1;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//j����Ȃ�ii��0.5�i�q�������炷
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//k����Ȃ�ii��jj�����炷
				if(kk<Height-gap && jj<Depth-gap && ii<Width-gap)
				{
					if(kk<=Height*0.25 || jj<=Depth*0.25 || ii<=Width*0.25 || kk>=Height*0.75 || jj>=Depth*0.75 || ii>=Width*0.75)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}


	//���q���͊w�ɂ��ʒu���œK��
	double r=1.5;
	double rigion[3][2];	//��͗̈�
	rigion[A_X][0]=-5*le; rigion[A_X][1]=1.5*Width;
	rigion[A_Y][0]=-5*le; rigion[A_Y][1]=1.5*Depth;
	rigion[A_Z][0]=-5*le; rigion[A_Z][1]=Height*1.5;

	MD_3D(X,Y,Z,le,BO_startID,BO_lastID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//���������֐�
void make_fusion3D(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<double> &X2,vector<double> &Y2,vector<double> &Z2,vector<int> &surface2,int *number,double le)
{
	//���݂��D�悳��闱�q�̍��W��X,Y,Z�ɁA���݂��D�悳��Ȃ����q�̍��W��X2,Y2,Z2�Ɋi�[����Ă���B
	//�܂��A���݂��D�悳��Ȃ����q���A�\�ʗ��q���������q����surface2[i]�Ɋi�[����Ă���B�\�ʗ��q�͕��q���͊w�œ������Ȃ��悤�ɂ��Ȃ��Ƃ����Ȃ�

	size_t pri_num=X.size();			//�D�敨�̂̍\�����q��
	size_t neg_num=X2.size();			//������镨�̂̍\�����q��

	//cout<<X2.size()<<" "<<surface2.size()<<endl;

	double r=1.5;
	double region[3][2];	//��͗̈�
	int *flag=new int[neg_num];			//flag��ON�Ȃ瑶�݂��������BOFF�Ȃ�������
	int newN;
	int beforeN=*number;

	////////////////////��͗̈�̌���
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	region[A_Z][0]=100; region[A_Z][1]=-100;
	for(int i=0;i<pri_num;i++)
	{
		if(X[i]<region[A_X][0]) region[A_X][0]=X[i];
		else if(X[i]>region[A_X][1]) region[A_X][1]=X[i];
		if(Y[i]<region[A_Y][0]) region[A_Y][0]=Y[i];
		else if(Y[i]>region[A_Y][1]) region[A_Y][1]=Y[i];
		if(Z[i]<region[A_Z][0]) region[A_Z][0]=Z[i];
		else if(Z[i]>region[A_Z][1]) region[A_Z][1]=Z[i];
	}
	for(int i=0;i<neg_num;i++)
	{
		if(X2[i]<region[A_X][0]) region[A_X][0]=X2[i];
		else if(X2[i]>region[A_X][1]) region[A_X][1]=X2[i];
		if(Y2[i]<region[A_Y][0]) region[A_Y][0]=Y2[i];
		else if(Y2[i]>region[A_Y][1]) region[A_Y][1]=Y2[i];
		if(Z2[i]<region[A_Z][0]) region[A_Z][0]=Z2[i];
		else if(Z2[i]>region[A_Z][1]) region[A_Z][1]=Z2[i];
	}

	for(int D=0;D<3;D++)
	{
		region[D][0]-=5*le;//�ی��̈Ӗ��ŏ����L�߂ɂƂ�
		region[D][1]+=5*le;
	}
	//////////////////////////////

	//�v�Z�̍������̂��߂Ɋi�q���`�� ��͕���r*le�Ŋ���؂��Ƃ͌���Ȃ��̂ŁA�͂ݏo�����Ƃ���͐؂�̂āB�Ȃ̂Ŋe���Ƃ����̕����ɂ͗]�T��������
	double grid_width=le*((int)(r+1));								//�i�q�̕��Br���܂ސ���*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X�����̊i�q�̌�
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[neg_num];									//������镨�̂̍\�����q���܂ފi�q�ԍ�
	vector<int> *MESH_pri=new vector<int>[grid_SIZE];				//�e���b�V���Ɋi�[�����D�旱�qID�i�[
	//vector<int> *MESH_neg=new vector<int>[grid_SIZE];				//�e���b�V���Ɋi�[������D�旱�qID�i�[

	for(int i=0;i<pri_num;i++)	//�܂��͗D�旱�q���i�q�Ɋi�[
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
		MESH_pri[number].push_back(i);
	}
	for(int i=0;i<neg_num;i++)	//���ɔ�D�旱�q���i�q�Ɋi�[
	{
		int xn=(int)((X2[i]-region[A_X][0])/grid_width);//X�����ɉ��ڂ̊i�q�� 
		int yn=(int)((Y2[i]-region[A_Y][0])/grid_width);//Y�����ɉ��ڂ̊i�q��
		int zn=(int)((Z2[i]-region[A_Z][0])/grid_width);//Z�����ɉ��ڂ̊i�q��
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//���qi���܂ފi�q�̔ԍ�
		//MESH_neg[number].push_back(i);
		index[i]=number;
	}

	//flag[i]�v�Z�J�n
	for(int i=0;i<neg_num;i++)
	{
		int G_id=index[i];				//�i�[����i�q�ԍ�
		flag[i]=ON;
		for(int II=G_id-1;II<=G_id+1;II++)
		{       
			for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
			{
				for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
				{
					int M_id=II+JJ+KK;
					for(int L=0;L<MESH_pri[M_id].size();L++)//�ߗׂ̗D�旱�q��T��
					{
						int j=MESH_pri[M_id][L];
						
						double x=X[j]-X2[i];
						double y=Y[j]-Y2[i];
						double z=Z[j]-Z2[i];
						double dis=sqrt(x*x+y*y+z*z);
						if(dis<0.7*le)
						{
							flag[i]=OFF;	
						}
					}
				}
			}
		}
	}///////////////////

	//flag[i]=ON�̗��q�̂�X,Y,Z�ɒǉ�
	newN=0;
	for(int i=0;i<neg_num;i++)
	{
		if(flag[i]==ON && surface2[i]==ON)//flag=ON���\�ʗ��q
		{
			X.push_back(X2[i]);
			Y.push_back(Y2[i]);
			Z.push_back(Z2[i]);
			newN++;
		}
	}
	*number=*number+newN;//���q���͊w�ł͂����܂ł̗��q���Œ藱�q����
	beforeN=*number;

	newN=0;
	for(int i=0;i<neg_num;i++)
	{
		if(flag[i]==ON && surface2[i]==OFF)//flag=ON���������q
		{
			X.push_back(X2[i]);
			Y.push_back(Y2[i]);
			Z.push_back(Z2[i]);
			newN++;
		}
	}

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,region);

	*number=*number+newN;

	delete [] index;
	delete [] MESH_pri;
	//delete [] MESH_neg;
	delete [] flag;
	
}

