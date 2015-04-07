#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
#include"define.h"	//#define �i�[
#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢
#include"PART.h"		//class PART��`
#include"BEMclass.h"	//BEM2D�֌W��class ��`
#include"FEM3Dclass.h"	//FEM3D�֌W��class ��`
#include<omp.h>
#include<vector>
#include"function.h"

///�\�ʒ��͌v�Z�J�n�֐�
void calc_surface_tension(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double dt,int particle_number,double n0,double **potential,int t)
{
	int d=CON->get_dimention();

	if(CON->get_surface_tension()>0)
	{
		if(CON->get_surface_tension()==1)////�d�݂��ŏ����@�����Aver4�ƈႢ�A�@�����U�ł͂Ȃ��A�\�ʋȖ�(��)�����߂��̒��ڔ����ɂ��ȗ������Ƃ߂�
		{
			surface_tension1(CON,PART,fluid_number,potential,particle_number);
		}
		else if(CON->get_surface_tension()==2)////���q�ԃ|�e���V����
		{
			surface_tension2(CON,PART,fluid_number,potential,particle_number);
		}
		else cout<<"�w�肳�ꂽ�\�ʒ��͌v�Z���@����������܂���"<<endl;
		
		//�X���[�W���O
		if(CON->get_smooth()==ON) smoothing(CON,PART,particle_number,fluid_number,potential,n0);

		/*/////�\�ʒ��͕\��
		ofstream vec("vector.dat");
		double le=CON->get_distancebp();
		double times=CON->get_times()*le*le;
		if(d==2) for(int i=0;i<fluid_number;i++) vec<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<potential[A_X][i]*times<<" "<<potential[A_Y][i]*times<<endl;
		else if(d==3) for(int i=0;i<fluid_number;i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) vec<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<potential[A_X][i]*times<<"\t"<<potential[A_Z][i]*times<<endl;
		vec.close();
		////////////////*/

		for(int D=0;D<d;D++) for(int i=0;i<fluid_number;i++) PART[i].potential[D]=potential[D][i];

		//�\�ʒ��͕\��
		plot_ST(CON,PART,fluid_number,potential,t);
		
		if(CON->get_ST_interval()>0)
		{
			//if(t==1 || t%CON->get_ST_interval()==0) plot_ST_each(CON,PART,fluid_number,potential,t);
			if(t==1 || (t-1)%CON->get_EM_interval()==0) plot_ST_each(CON,PART,fluid_number,potential,t);
		}//*/

		/*/�\�ʒ���AVS�t�@�C���o�͊֐�
		if(CON->get_avs_stension_interval()>0)
		{
			if(CON->get_current_step()==1 || CON->get_current_step()%CON->get_avs_stension_interval()==0) plot_avs_stension(CON,PART,fluid_number);
		}//*/

	}
	else for(int D=0;D<d;D++) for(int i=0;i<fluid_number;i++) potential[D][i]=0;//�v�Z���Ȃ��ꍇ���������������Ă���
}

////WLSM�ɂ��A�\�ʋȖ�(��)��@
void surface_tension1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int particle_number)
{   
	cout<<"WLSM�ɂ��\�ʒ��͌v�Zver.1---------";

	unsigned int timeA=GetTickCount();
    double le=CON->get_distancebp();		//���ϗ��q�ԋ���
    double R=CON->get_re3()*le;				//���U�p�e�����a
    double mass=CON->get_particle_mass();

	int dim=CON->get_dimention();
	int N=0;							//�W���s��̌�
	int order=2;						//�ߎ��Ȗʂ̵��ް�B 2=�� 3=3��
	double err=(0.01*le)*(0.01*le);		//�G���[�l

    double *curv=new double [fluid_number];	//�e���q�̋ȗ��i�[
	double *W=new double [fluid_number];	//�e���q�̏d�݊i�[

    double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [particle_number];//�������@���x�N�g��

	double *modify[DIMENTION];
    for(int D=0;D<DIMENTION;D++) modify[D]=new double [fluid_number];//�������@���x�N�g��

	int *particle_J=new int[particle_number];	//particle_J[i]=ON�Ȃ�A���qi�͋ȗ��v�Z�ɍl�������B�ʏ�͗��̕\�ʗ��q�̂�ON.�������ݽد�ߍl���̍ۂ͕Ǖ\�ʗ��q��ON
	int *calced=new int[fluid_number];	//calced[i]=ON�Ȃ�A���qi�͋ȗ��v�Z�����s���ꂽ�BOFF�Ȃ�ߗח��q���s���Ŏ��s����Ă��Ȃ�

	//������
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENTION;D++)
		{
			potential[D][i]=0;
			modify[D][i]=0;	
		}
		curv[i]=0;
		W[i]=0;
		if(PART[i].surface==ON) particle_J[i]=ON;//�ȗ��v�Z�Ɋ�^����
		else particle_J[i]=OFF;
		calced[i]=OFF;
	}
	if(CON->get_non_slip()==ON)
	{
		for(int i=fluid_number;i<particle_number;i++)
		{
			if(PART[i].type==INWALL && PART[i].surface==ON) particle_J[i]=ON;//�Ǖ\�ʗ��q���ȗ��v�Z�Ɋ�^����
			else particle_J[i]=OFF;
		}
	}
	else if(CON->get_non_slip()==OFF) for(int i=fluid_number;i<particle_number;i++) particle_J[i]=OFF;//�ݽد�߂łȂ��Ȃ痬�̕\�ʗ��q�����ŋȗ��v�Z
		
	//�W���s��̑傫���̌���
	if(dim==2)
	{
		if(order==2) N=3;
		else if(order==3) N=4;
	}
	else if(dim==3)
	{
		if(order==2) N=6;
		else if(order==3) cout<<"order=3??"<<endl;
	}
	////////////////////////////////

	double *matrix=new double [N*N];	//N�~N�̌W���s��
	double *Bx=new double [N];	//N�̉��s��
	double *By=new double [N];	//N�̉��s��
	double *Bz=new double [N];	//N�̉��s��
    
	//////�@���޸�ٌv�Z  
	ofstream fq("direct.dat");
    for(int i=0;i<particle_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
	
		if(PART[i].type==INWALL && PART[i].surface==ON)
		{
			if(CON->get_wall_direct()==OFF) for(int D=0;D<DIMENTION;D++) direct[D][i]=0;
			else if(CON->get_wall_direct()==2)//2DFEM�Ód�����̂Ƃ��͂���
			{
				if(CON->get_dimention()==2)
				{
					direct[A_X][i]/=sqrt(direct[A_X][i]*direct[A_X][i]);//���������̂ݎc��
					direct[A_Y][i]=0;
				}
				else if(CON->get_dimention()==3)
				{
					double r=sqrt(direct[A_X][i]*direct[A_X][i]+direct[A_Y][i]*direct[A_Y][i]);
					if(r!=0)
					{
						direct[A_X][i]/=r;
						direct[A_Y][i]/=r;
					}
					direct[A_Z][i]=0;//���������̂ݎc��
				}
			}
			else if(CON->get_wall_direct()==3)
			{
				if(CON->get_dimention()==2)
				{
					direct[A_Y][i]/=sqrt(direct[A_Y][i]*direct[A_Y][i]);//��������
					direct[A_X][i]=0;
				}
				else if(CON->get_dimention()==3)
				{
					direct[A_X][i]=0;
					direct[A_Y][i]=0;
					direct[A_Z][i]/=sqrt(direct[A_Z][i]*direct[A_Z][i]);//����
				}
			}
		}

		if(CON->get_dimention()==2) fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<direct[A_X][i]*0.0001<<" "<<direct[A_Y][i]*0.0001<<endl;
		else if(CON->get_dimention()==3)
		{
			if(PART[i].r[A_Y]>-le && PART[i].r[A_Y]<le)
			{
				fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<direct[A_X][i]*0.0001<<" "<<direct[A_Z][i]*0.0001<<endl;
			}
		}
    }
	fq.close();
    ////////////////*/

	if(dim==2 && order==2)//�񎟌�2����
	{
		//stand_d[i]��A_X��A_Y�̂ǂ��炩
		//�ߎ��Ȑ���y=ax2+bx+c �����ł�x�͎c���@y�͐^�l
		///�W���s���
		///  ����x4   ����x3  ����x2  a = ����x2*yj  
		///  ����x3   ����x2  ����x   b = ����x*yj 
		///  ����x2   ����x   ��1     c = ��yj

		//���̂Ƃ�y��2�K������y'',1�K������y'�Ƃ����Ȃ�A�ȗ�k��y''/(1+y'^2)^1.5�ƂȂ�B�����ŁAy''=2a, y'=2ax+b�ł���(���������qi�ɂ����Ă�x=0�ƂȂ邩��A���� y'=b)
		
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;//������
				for(int n=0;n<N;n++){By[n]=0;}
				int num=1;

				double Cx=PART[i].r[A_X]+direct[A_X][i]*le;//�_(Cx,Cy)�܂��ɁA���ӗ��q���ړ�������
				double Cy=PART[i].r[A_Y]+direct[A_Y][i]*le;

				double SIN=-direct[A_X][i];//���̂悤�ɒ�`����cos�Ƃ�sin�Ƃ��������ē_(Cx,Cy)�܂��ɉ�]������B
				double COS=-direct[A_Y][i];

				double Xi=COS*PART[i].r[A_X]-SIN*PART[i].r[A_Y]+Cx-Cx*COS+Cy*SIN;//�_(Cx,Cy)�܂��ɃƉ�]���������qi�̍��W
				double Yi=SIN*PART[i].r[A_X]+COS*PART[i].r[A_Y]+Cy-Cx*SIN-Cy*COS;

				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(particle_J[j]==ON)//ON�̗��q�̂݋ȗ��v�Z�Ɋ�^
					{
						double product=direct[A_X][i]*direct[A_X][j]+direct[A_Y][i]*direct[A_Y][j];//���qi�Ɨ��qj�̖@���x�N�g���̓���

						if(product>0)//���ς����̗��q�͍l�����Ȃ��B���ɗ��̂������re3>dis�ɂȂ�Ƃ܂�������
						{
							//double Xj=COS*PART[j].r[A_X]-SIN*PART[j].r[A_Y]+Cx-Cx*COS+Cy*SIN;//�_(Cx,Cy)�܂��ɃƉ�]���������qj�̍��W
							double Yj=SIN*PART[j].r[A_X]+COS*PART[j].r[A_Y]+Cy-Cx*SIN-Cy*COS;

							double dX=PART[j].r[A_X]-PART[i].r[A_X];//���ۂ̍��W�n�ł̗��qj�Ɨ��qi��X���W�̍�
							double dY=PART[j].r[A_Y]-PART[i].r[A_Y];//���ۂ̍��W�n�ł̗��qj�Ɨ��qi��Y���W�̍�
	
							double X=COS*dX-SIN*dY;//��]���������W�n�ł̗����q��X���W�̍�.�������傭�޸�ق̉�]�ł��邱�ƂɋC�Â�����.(Xj-Xi�ƌv�Z���Ă���������)
	
							double Y=Yj;
							double dis=sqrt(X*X);//�����ł̋�����x�����̋���
							
							double w=1;
							if(dis>le) w=le*le/(dis*dis);
							
							matrix[0]+=X*X*X*X*w;			
							matrix[1]+=X*X*X*w;
							matrix[2]+=X*X*w;		
								
							matrix[5]+=X*w;		
			
							matrix[8]+=w;			
		
							By[0]+=X*X*Y*w;
							By[1]+=X*Y*w;
							By[2]+=Y*w;
		
							num++;
							
						}
					}
				}
				//if(i==0) cout<<By[0]<<" "<<By[1]<<" "<<By[2]<<endl;
				if(num>2)//���ӕ\�ʗ��q���������菭�Ȃ��ƃG���[
				{
					calced[i]=ON;		//�v�Z�����Ƃ���������Ă���
					matrix[3]=matrix[1];		
					matrix[4]=matrix[2];
					matrix[6]=matrix[2];		
					matrix[7]=matrix[5];
		
					matrix[8]+=1;//�������g
					//By[2]+=PART[i].r[stand_d[i]];//�������g
					By[2]+=Yi;//�������g
					//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�
					
					//�s����K�E�X�̏����@�ŉ����@����By�Ɋi�[�����
					gauss(matrix,By,N);

					double a=By[0];
					double b=By[1];
					double c=By[2];
					
					double y2=2*a;					//y''
					double y1=b;					//y' 
					double A=pow((1+y1*y1),1.5);	//(1+y'^2)^1.5
						
					curv[i]=-y2/A;//���ϲŽ������
					//cout<<By[0]<<" "<<By[1]<<" "<<By[2]<<endl;

					//�C���ʌv�Z
					double Yc=c-Yi;//�C����
					if(fabs(Yc)>0.5*le) Yc=0;
					for(int D=0;D<3;D++) modify[D][i]+=-Yc*direct[D][i];
					W[i]+=1;
					
				}
				else
				{
					calced[i]=OFF;		//�v�Z���Ă��Ȃ��Ƃ���������Ă���
					cout<<"���q���s���H"<<endl;
				}
			}
		}
	}
	else if(dim==2 && order==3)//�񎟌�3����
	{
		//stand_d[i]��A_X��A_Y�̂ǂ��炩
		//�ߎ��Ȑ���y=ax3+bx2+cx+d �����ł�x�͎c���@y�͐^�l
		///�W���s���
		///  ����x6   ����x5  ����x4  ����x3  a = ����x3*yj  
		///  ����x5   ����x4  ����x3  ����x2  b = ����x2*yj 
		///  ����x4   ����x3  ����x2  ����x   c = ����x*yj
		///  ����x3   ����x2  ����x1  ��1     d = ����yj

		//���̂Ƃ�y��2�K������y'',1�K������y'�Ƃ����Ȃ�A�ȗ�k��y''/(1+y'^2)^1.5�ƂȂ�B�����ŁAy''=6ax+2b, y'=3ax2+2bx+c�ł���(���������qi�ɂ����Ă�x=0�ƂȂ邩��A���� y''=2b, y'=c�ƂȂ�)
		
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;//������
				for(int n=0;n<N;n++){By[n]=0;}
				int num=1;

				double Cx=PART[i].r[A_X]+direct[A_X][i]*le;//�_(Cx,Cy)�܂��ɁA���ӗ��q���ړ�������
				double Cy=PART[i].r[A_Y]+direct[A_Y][i]*le;

				double SIN=-direct[A_X][i];//���̂悤�ɒ�`����cos�Ƃ�sin�Ƃ��������ē_(Cx,Cy)�܂��ɉ�]������B
				double COS=-direct[A_Y][i];

				double Xi=COS*PART[i].r[A_X]-SIN*PART[i].r[A_Y]+Cx-Cx*COS+Cy*SIN;//�_(Cx,Cy)�܂��ɃƉ�]���������qi�̍��W
				double Yi=SIN*PART[i].r[A_X]+COS*PART[i].r[A_Y]+Cy-Cx*SIN-Cy*COS;

				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(particle_J[j]==ON)//ON�̗��q�̂݋ȗ��v�Z�Ɋ�^
					{
						//double Xj=COS*PART[j].r[A_X]-SIN*PART[j].r[A_Y]+Cx-Cx*COS+Cy*SIN;//�_(Cx,Cy)�܂��ɃƉ�]���������qj�̍��W
						double Yj=SIN*PART[j].r[A_X]+COS*PART[j].r[A_Y]+Cy-Cx*SIN-Cy*COS;

						double dX=PART[j].r[A_X]-PART[i].r[A_X];//���ۂ̍��W�n�ł̗��qj�Ɨ��qi��X���W�̍�
						double dY=PART[j].r[A_Y]-PART[i].r[A_Y];//���ۂ̍��W�n�ł̗��qj�Ɨ��qi��Y���W�̍�

						double X=COS*dX-SIN*dY;//��]���������W�n�ł̗����q��X���W�̍�.�������傭�޸�ق̉�]�ł��邱�ƂɋC�Â�����.(Xj-Xi�ƌv�Z���Ă���������)

						double Y=Yj;
						double dis=sqrt(X*X);//�����ł̋�����x�����̋���
						
						double w=1;
						if(dis>le) w=le*le/(dis*dis);
						
						matrix[0]+=X*X*X*X*X*X*w;			
						matrix[1]+=X*X*X*X*X*w;	
						matrix[2]+=X*X*X*X*w;
						matrix[3]+=X*X*X*w;
							
						matrix[7]+=X*X*w;	

						matrix[11]+=X*w;
		
						matrix[15]+=w;
								
						By[0]+=X*X*X*Y*w;
						By[1]+=X*X*Y*w;
						By[2]+=X*Y*w;
						By[3]+=Y*w;
	
						num++;
					}
				}
				if(num>4)//���ӕ\�ʗ��q���������菭�Ȃ��ƃG���[
				{
					calced[i]=ON;		//�v�Z�����Ƃ���������Ă���

					matrix[4]=matrix[1];
					matrix[5]=matrix[2];
					matrix[6]=matrix[3];

					matrix[8]=matrix[5];
					matrix[9]=matrix[6];
					matrix[10]=matrix[7];

					matrix[12]=matrix[9];
					matrix[13]=matrix[10];
					matrix[14]=matrix[11];

					matrix[15]+=1;//�������g
					By[3]+=Yi;//�������g
					//����ȊO��X=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�

					//�s����K�E�X�̏����@�ŉ����@����By�Ɋi�[�����
					gauss(matrix,By,N);

					double a=By[0];
					double b=By[1];
					double c=By[2];
					double d=By[3];

					double y2=2*b;					//y''
					double y1=c;					//y' 
					double A=pow((1+y1*y1),1.5);	//(1+y'^2)^1.5
						
					curv[i]=-y2/A;//���ϲŽ������

					//�C���ʌv�Z
					double Yc=d-Yi;//�C����
					for(int D=0;D<3;D++) modify[D][i]+=-Yc*direct[D][i];
					W[i]+=1;
				}
				else cout<<"error in surface_tension ver.6 ���q��<=4"<<endl;
			}
		}
			
	}
	else if(dim==3 && order==2)//3����1����
	{
		//�ߎ��Ȑ���z=ax2+by2+cxy+dx+ey+f �����ł�x,y�͎c���@z�͐^�l
		///�W���s���
		///  ����x4      ����x2��y2  ����x3��y  ����x3     ����x2��y  ����x2    a = ����x2*zj  
		///  ����x2��y2  ����y4      ����x��y3  ����x��y2  ����y3     ����y2    b = ����y2*zj
		///  ����x3��y   ����x��y3   ����x2��y2 ����x2��y  ����x��y2  ����x��y  c = ����x��y*zj
		///  ����x3      ����x��y2   ����x2��y  ����x2     ����x��y   ����x     d = ����x*zj
		///  ����x2��y   ����y3      ����x��y2  ����x��y   ����y2     ����y     e = ����y*zj
		///  ����x2      ����y2      ����x��y   ����x      ����y      ��1       f = ��zj

		//���̂Ƃ�z��2�K������Zxx,Zyy,1�K������Zx,Zy�Ƃ����Ȃ�A�ȗ�k��-div{��Z/(1+|��Z|^2)^1.5} = -(Zxx+Zyy+Zxx*Zy^2+Zyy*Zx^2)/(1+Zx^2+Zy^2)^1.5�ƂȂ�B�����ŁA���qi�ɂ����Ă�x=y=0�ƂȂ邩��A���� Zxx=2a, Zyy=2b, Zx=d, Zy=e�ƂȂ�)
		
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;//������
				for(int n=0;n<N;n++){Bz[n]=0;}
				int num=1;
				
				double Cx=PART[i].r[A_X]+direct[A_X][i]*le;//�_(Cx,Cy,Cz)�܂��ɁA���ӗ��q���ړ�������
				double Cy=PART[i].r[A_Y]+direct[A_Y][i]*le;
				double Cz=PART[i].r[A_Z]+direct[A_Z][i]*le;

				double normal=sqrt(direct[A_Y][i]*direct[A_Y][i]+direct[A_Z][i]*direct[A_Z][i]);

				double SIN1=-direct[A_Y][i]/normal;//���̂悤�ɒ�`����cos�Ƃ�sin�Ƃ��������ē_(Cx,Cy,Cz)��ʂ�X���܂��ɂ܂��͉�]������B
				double COS1=-direct[A_Z][i]/normal;

				double X1=PART[i].r[A_X];
				double Y1=COS1*PART[i].r[A_Y]-SIN1*PART[i].r[A_Z]+Cy-COS1*Cy+SIN1*Cz;//�_(Cx,Cy,Cz)��ʂ�X���܂��ɉ�]��������̍��W
				double Z1=SIN1*PART[i].r[A_Y]+COS1*PART[i].r[A_Z]+Cz-SIN1*Cy-COS1*Cz;

				normal=sqrt(direct[A_X][i]*direct[A_X][i]+(SIN1*direct[A_Y][i]+COS1*direct[A_Z][i])*(SIN1*direct[A_Y][i]+COS1*direct[A_Z][i]));

				double SIN2=direct[A_X][i]/normal;//���̂悤�ɒ�`����cos�Ƃ�sin�Ƃ���������,���͓_(Cx,Cy,Cz)��ʂ�Y���܂��ɉ�]������B
				double COS2=-1*(SIN1*direct[A_Y][i]+COS1*direct[A_Z][i])/normal;//����͉�]���nz�ɑ���

				double Xi=COS2*X1+SIN2*Z1+Cx-COS2*Cx-SIN2*Cz;//�_(Cx,Cy,Cz)��ʂ�Y���܂��ɉ�]��������̍��W
				double Yi=Y1;
				double Zi=-SIN2*X1+COS2*Z1+Cz+SIN2*Cx-COS2*Cz;

				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
					if(particle_J[j]==ON)//ON�̗��q�̂݋ȗ��v�Z�Ɋ�^
					{
						double product=direct[A_X][i]*direct[A_X][j]+direct[A_Y][i]*direct[A_Y][j]+direct[A_Z][i]*direct[A_Z][j];//���qi�Ɨ��qj�̖@���x�N�g���̓���
						//if(product>0.5)//���ς����̗��q�͍l�����Ȃ��B���ɗ��̂������re3>dis�ɂȂ�Ƃ܂�������
						if(product>0)//���ς����̗��q�͍l�����Ȃ��B���ɗ��̂������re3>dis�ɂȂ�Ƃ܂�������
						{
							
							double Xj1=PART[j].r[A_X];
							double Yj1=COS1*PART[j].r[A_Y]-SIN1*PART[j].r[A_Z]+Cy-COS1*Cy+SIN1*Cz;//�_(Cx,Cy,Cz)��ʂ�X���܂��ɉ�]��������̍��W
							double Zj1=SIN1*PART[j].r[A_Y]+COS1*PART[j].r[A_Z]+Cz-SIN1*Cy-COS1*Cz;

							double Xj=COS2*Xj1+SIN2*Zj1+Cx-COS2*Cx-SIN2*Cz;//�_(Cx,Cy,Cz)��ʂ�Y���܂��ɉ�]��������̍��W
							double Yj=Yj1;
							double Zj=-SIN2*Xj1+COS2*Zj1+Cz+SIN2*Cx-COS2*Cz;

							double X=Xj-Xi;
							double Y=Yj-Yi;
							double Z=Zj;
							double dis=sqrt(X*X+Y*Y);//�����ł̋�����XY���ʏ�̋���
						
							double w=1;
							if(dis>le) w=le*le/(dis*dis);
						
							matrix[0]+=X*X*X*X*w;			
							matrix[1]+=X*X*Y*Y*w;	
							matrix[2]+=X*X*X*Y*w;
							matrix[3]+=X*X*X*w;
							matrix[4]+=X*X*Y*w;
							matrix[5]+=X*X*w;
							
							matrix[7]+=Y*Y*Y*Y*w;
							matrix[8]+=X*Y*Y*Y*w;
							matrix[9]+=X*Y*Y*w;
							matrix[10]+=Y*Y*Y*w;
							matrix[11]+=Y*Y*w;

							matrix[17]+=X*Y*w;

							matrix[23]+=X*w;

							matrix[29]+=Y*w;

							matrix[35]+=w;
								
							Bz[0]+=X*X*Z*w;
							Bz[1]+=Y*Y*Z*w;
							Bz[2]+=X*Y*Z*w;
							Bz[3]+=X*Z*w;
							Bz[4]+=Y*Z*w;
							Bz[5]+=Z*w;

							num++;
						}
					}
				}
				if(num>4)//���ӕ\�ʗ��q���������菭�Ȃ��ƃG���[
				{
					calced[i]=ON;		//�v�Z�����Ƃ���������Ă���

					matrix[6]=matrix[1];

					matrix[12]=matrix[2];
					matrix[13]=matrix[8];
					matrix[14]=matrix[1];
					matrix[15]=matrix[4];
					matrix[16]=matrix[9];

					matrix[18]=matrix[3];
					matrix[19]=matrix[9];
					matrix[20]=matrix[15];
					matrix[21]=matrix[5];
					matrix[22]=matrix[17];

					matrix[24]=matrix[4];
					matrix[25]=matrix[10];
					matrix[26]=matrix[16];
					matrix[27]=matrix[22];
					matrix[28]=matrix[11];

					matrix[30]=matrix[5];
					matrix[31]=matrix[11];
					matrix[32]=matrix[17];
					matrix[33]=matrix[23];
					matrix[34]=matrix[29];

					matrix[35]+=1;//�������g
					Bz[5]+=Zi;
					//����ȊO��X=0,Y=0�ɂȂ�̂ŉ��Z����K�v�͂Ȃ�

					//�s����K�E�X�̏����@�ŉ����@����Bz�Ɋi�[�����
					gauss(matrix,Bz,N);

					double a=Bz[0];
					double b=Bz[1];
					double c=Bz[2];
					double d=Bz[3];
					double e=Bz[4];
					double f=Bz[5];
					
					double Zxx=2*a;
					double Zyy=2*b;
					double Zx=d;
					double Zy=e;
					double A=pow((1+Zx*Zx+Zy*Zy),1.5);
					double A2=Zxx*(1+Zy*Zy)+Zyy*(1+Zx*Zx);
					
					curv[i]=-A2/A;
					

					//�C���ʌv�Z
					double Zc=f-Zi;//�C����
					
					if(2*curv[i]!=curv[i]+curv[i])//�G���[�Ȃ�
					{
						curv[i]=0;//���ׂĂ��[����
						Zc=0;
					}
					//if(Zc>0.1*le || Zc<-0.1*le) Zc=0;//�C���ʂɏ���A������݂���
					if(Zc>0.1*le) Zc=0.1*le;
					if(Zc<-0.1*le) Zc=-0.1*le;

					for(int D=0;D<3;D++) modify[D][i]+=-Zc*direct[D][i];
					W[i]+=1;
				}
				else curv[i]=0;
			}
		}
	}

	//�ʒu�C��
	if(CON->get_suf_modify()==ON)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(W[i]>0) for(int D=0;D<DIMENTION;D++) PART[i].r[D]+=modify[D][i]/W[i];
		}
	}

	int count;
	count=0;
	if(CON->get_interpolate_curv()==ON)
	{
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].surface==ON)
			{
				if(calced[i]==OFF)//�\�ʗ��q�Ȃ̂ɋȗ����v�Z����Ă��Ȃ�������
				{
					double W=0;
					for(int k=0;k<PART[i].N3;k++)
					{
						int j=PART[i].NEI3[k];
						if(PART[j].type==FLUID)
						{
							if(calced[j]==ON)//�������q��calced[]��OFF�Ȃ̂ŁA�����I�ɂ����ł͕\�ʗ��q�̂ݑΏۂɍi����
							{
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
								double w=(1-dis/R)*(1-dis/R);
								W+=w;
								curv[i]+=curv[j]*w;
							}
						}
					}
					if(W>0) curv[i]/=W;
					else 
					{
						//cout<<"in interpolate_curv() ���ӗ��q���s��"<<endl;
						count++;
					}
				}
			}
		}
	}

	if(count>0) cout<<"in interpolate_curv() ���ӂ��s�݂̗��q���� n="<<count<<endl;
	////////////////////////////////////curv�̃X���[�W���O
    if(CON->get_smth_sumP()!=OFF)
    {
        double *newP=new double [fluid_number];
        for(int n=0;n<CON->get_smth_sumP();n++)
        {
            for(int i=0;i<fluid_number;i++) 
            {  
				if(PART[i].surface==OFF) newP[i]=0;
				else
				{  
					newP[i]=curv[i];
					int num=1;
					double W=1;
					for(int k=0;k<PART[i].N3;k++)
					{
						int j=PART[i].NEI3[k];
						if(PART[j].type==FLUID &&PART[j].surface==ON ) 
						{ 	
							double dif=fabs(curv[i]-curv[j]);//�ȗ��̍�
						//	if(dif<400)//�ȗ��̍������܂�ɂ����������͕̂������ɍl�����Ȃ�
							{
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
								double w=(1-dis/R)*(1-dis/R);
								W+=w;
								num++;
								newP[i]+=curv[j]*w;
							}
						}
					}
					newP[i]/=W;
				}
            } 
            for(int i=0;i<fluid_number;i++) curv[i]=newP[i];
		}
		delete [] newP;
    }
    ////////////////////*/

	///�ȗ���̧�ُo��
	ofstream fp2("curv.dat");///�ȗ��o��
	if(dim==2) for(int i=0;i<fluid_number;i++) if(PART[i].surface==ON) fp2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<curv[i]<<endl;
	if(dim==3) for(int i=0;i<fluid_number;i++) if(PART[i].surface==ON) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fp2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<curv[i]<<endl;
	fp2.close();///////////*/

	double *sigma=new double[fluid_number];
	for(int i=0;i<fluid_number;i++) sigma[i]=0;//������
	///potential[D][i]�ɑ��
	if(CON->get_SFT()==OFF)
	{
		double L=1.5*le;//�\�ʂ̌���
		double density;
		for(int i=0;i<fluid_number;i++)
		{   
			sigma[i]=CON->get_sigma();
		    if(PART[i].surface==ON)
			{
				if(PART[i].materialID==1) density=CON->get_density();
				else if(PART[i].materialID==2) density=CON->get_density2();
				for(int D=0;D<DIMENTION;D++) potential[D][i]+=CON->get_sigma()*curv[i]/(density*L)*direct[D][i];
			}
		}
	}
	else 
	{
		for(int i=0;i<fluid_number;i++)
		{   
			sigma[i]=CON->get_sigma();//�ЂƂ܂����̐ݒ�l�ŏ�����
		}
		//cout<<"���݂͕\�ʒ��͂̉��x�ˑ����͌v�Z�s��"<<endl;
		double L=1.5*le;//�\�ʂ̌���
		double *density=new double[fluid_number];

		calc_physical_property(CON,PART,fluid_number,density,particle_number,1);//���x�̉��x�ˑ�
		calc_physical_property(CON,PART,fluid_number,sigma,particle_number,3);//�\�ʒ��͌W���̉��x�ˑ�

		for(int i=0;i<fluid_number;i++)
		{   
		    if(PART[i].surface==ON)
			{
				for(int D=0;D<DIMENTION;D++) potential[D][i]+=sigma[i]*curv[i]/(density[i]*L)*direct[D][i];
			}
		}
		
		delete [] density;
		
	}
	
	
	if(CON->get_dir_for_P()==1 ||CON->get_dir_for_P()==3 )//���͂̌v�Z���ɁA�\�ʒ��͂��f�B���N���l�Ƃ��Ďg�p����
	{
		//ofstream dir("surface_tension_P.dat");
		for(int i=0;i<fluid_number;i++) PART[i].dir_Pst=sigma[i]*curv[i];			//�f�B���N���l���Z�b�g
		//for(int i=0;i<fluid_number;i++) PART[i].dir_Pst=sigma[i]*curv[i];			//�f�B���N���l���Z�b�g
		//for(int i=0;i<fluid_number;i++)  dir<<CON->get_sigma()*curv[i]<<endl;
		//dir.close();
	}

    for(int D=0;D<DIMENTION;D++) delete [] direct[D];
    delete [] curv;
	delete [] W;
	delete [] matrix;
	delete [] Bx;
	delete [] By;
	delete [] Bz;

	for(int D=0;D<DIMENTION;D++) delete [] modify[D];
	delete [] particle_J;
	delete [] calced;
	delete []sigma;

	cout<<"ok  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//���q�ԃ|�e���V����
void surface_tension2(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int particle_number)
{
	cout<<"���q�Ԉ��͂ɂ��ƂÂ��\�ʒ��͌v�Z�J�n--------";

	unsigned int timeA=GetTickCount();		//�v�Z�J�n����
    double le=CON->get_distancebp();		//���ϗ��q�ԋ���
    double r=CON->get_re3();				//�\�ʒ��͗p�e�����a
	double mass=CON->get_particle_mass();	//���q�̎���[kg]
	double Cst=CON->get_Cst();				//�ꉞ�ŏ��Ɍv�Z���ꂽ�|�e���V�����W��
    double C;								//�|�e���V�����W��
	double C2=Cst*CON->get_C_times();		//�|�e���V�����W���̕␳(��������Ȃ��Ɨ��q�ԗ͂������Ȃ肷����)
	double wall_C=CON->get_wall_C();		//�ǂƂ̐e�a�͌W��

	if(CON->get_freeon()==1)
	{
		if(CON->get_SFT()==OFF)//�\�ʒ��͂̉��x�ˑ����Ȃ�
		{
		    for(int i=0;i<fluid_number;i++)
			{   
				for(int D=0;D<DIMENTION;D++) potential[D][i]=0;
				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
				           
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
						    
					C=C2;
			  
					double p_slope;//p(r)�̌��z
					//  else   C=Cst*(1+cos(PI/2))/2; //�_���ʂ�
	
					if(PART[j].type!=FLUID) C=C2*wall_C;
					p_slope=C*(dis-le)*(dis-le*r);
					p_slope/=mass;	//���̈ʒu�E���x�̌v�Z��mass�Ŋ���̂ł����ł�[N]�ŋ��߂�
					potential[A_X][i]-=p_slope*X/dis;
					potential[A_Y][i]-=p_slope*Y/dis;	
					potential[A_Z][i]-=p_slope*Z/dis;		
				}
			}
		}
    }//potential[D][i]�����܂���*/
	
	if(CON->get_freeon()!=1)
	{
		if(CON->get_SFT()==OFF)//�\�ʒ��͂̉��x�ˑ����Ȃ�
		{
			int *check=new int[fluid_number];//0�Ȃ疢�����@1�͏����ς�
			for(int i=0;i<fluid_number;i++)
			{
				check[i]=0;//������
				for(int D=0;D<DIMENTION;D++) potential[D][i]=0;
			}
		    for(int i=0;i<fluid_number;i++)
			{   
				for(int k=0;k<PART[i].N3;k++)
				{
					int j=PART[i].NEI3[k];
				    if(j<fluid_number)//���qj�����̗��q�Ȃ�
					{
						if(check[j]==0)//�܂��������I����ĂȂ��Ȃ�
						{
							double X=PART[j].r[A_X]-PART[i].r[A_X];
							double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
							double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
							    
							C=C2;
					
							double p_slope;//p(r)�̌��z
							
							p_slope=C*(dis-le)*(dis-le*r);
							p_slope/=mass;
							potential[A_X][i]-=p_slope*X/dis;//���qj�����qi�ɂ���ڂ���
							potential[A_Y][i]-=p_slope*Y/dis;
							potential[A_Z][i]-=p_slope*Z/dis;

							potential[A_X][j]+=p_slope*X/dis;//���qi�����qj�ɂ���ڂ���
							potential[A_Y][j]+=p_slope*Y/dis;
							potential[A_Z][j]+=p_slope*Z/dis;
						}
					}
					else
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
							    
						C=C2*wall_C;
					
						double p_slope;//p(r)�̌��z
						
						p_slope=C*(dis-le)*(dis-le*r);
						p_slope/=mass;
						potential[A_X][i]-=p_slope*X/dis;//���qj�����qi�ɂ���ڂ���
						potential[A_Y][i]-=p_slope*Y/dis;
						potential[A_Z][i]-=p_slope*Z/dis;
					}
				}
				check[i]=1;
			}
			delete [] check;
	    }//potential[D][i]�����܂���
	}////*/

    /////�\�ʒ��͂̉��x�ˑ������l����ꍇ
    if(CON->get_SFT()==1)
    {    
        double hs0=mass*CON->get_Cp()*CON->get_MP();//�Z���J�n�_�̃G���^���s�[
		double hs1=hs0+CON->get_latent_H()*mass;    //�Z���I���_�̃G���^���s�[
		double *C1=new double[fluid_number];//���qi�̕\�ʒ��͌W��
		///�G���^���s�[����e���q�̉��xT[i]�����߂�
		for(int i=0;i<fluid_number;i++)
		{
		    //���q�͂��ׂĉt��
		    double T=CON->get_MP()+(PART[i].h-hs1)/mass/CON->get_Cp();
		    C1[i]=-0.0003*T*T+1.0343*T+951.65;//���̒l��mN/m
		    C1[i]/=1000;//�P�ʂ�N/m�ɂȂ����B
		    C1[i]/=CON->get_sigma();//�W���l�Ƃ̔䗦�������B
		    //cout<<C1[i]<<endl;
		}
		  
		for(int i=0;i<fluid_number;i++)
		{
		    for(int D=0;D<DIMENTION;D++) potential[D][i]=0;
		    for(int k=0;k<PART[i].N3;k++)
		    {
		        int j=PART[i].NEI3[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				
				if(dis<r*le && i!=j)
				{   
					double p_slope;//p(r)�̌��z
					if(PART[j].type==FRFLUID ||PART[j].type==BOFLUID) C=Cst*C1[i]*C1[j]/(C1[i]+C1[j])*2;
					else if(PART[j].type!=BOFLUID && PART[j].type!=FRFLUID) C=1.0*Cst*C1[i];
					p_slope=C*(dis-le)*(dis-le*r);
					p_slope/=mass;
					potential[A_X][i]-=p_slope*X/dis;
					potential[A_Y][i]-=p_slope*Y/dis;
					potential[A_Z][i]-=p_slope*Z/dis;
				}
			} 
		}//potential[D][i]�����܂��� */
		delete [] C1;
    }
	/////////////////*/

	
	cout<<"ok  time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}


//�\�ʒ��̓X���[�W���O�֐�
void smoothing(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double *potential[DIMENTION],double n0)
{
	int d=CON->get_dimention();
	double le=CON->get_distancebp();
    double r=CON->get_re()*le;  //�Ѱ��ݸނ͈̔�

	///potential[D][i]�̃X���[�W���O////////
	double *newtension[DIMENTION];
	for(int D=0;D<DIMENTION;D++) newtension[D]=new double [fluid_number];
	
	for(int n=0;n<CON->get_smn();n++)
    {
        for(int i=0;i<fluid_number;i++) 
        {  
            for(int D=0;D<3;D++) newtension[D][i]=potential[D][i];
			int num=1; //�������g���Ă��邩��1
			for(int k=0;k<PART[i].N;k++)
			{       
				int j=PART[i].NEI[k];
	        
				if(PART[j].type==FLUID)
				{
					num++;
					for(int D=0;D<3;D++) newtension[D][i]+=potential[D][j];
				}
			}
			for(int D=0;D<3;D++) newtension[D][i]/=num;
        } 
        for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) potential[D][i]=newtension[D][i];
    }
    /////////////////////////////////////*/
	
	for(int D=0;D<DIMENTION;D++) delete [] newtension[D];
}

double calc_Cst(mpsconfig *CON)
{
	double p;//�|�e���V����
	double sumP=0;//��p(r)
	double le=CON->get_distancebp();	//�������q�ԋ���
	double re=CON->get_re3()*le;		//�e�����a
	int calc_type=CON->get_model_set_way();	//���f���̃Z�b�g�^�C�v(���� or MD)

	if(CON->get_dimention()==2)
	{
		for(int k=1;k<10;k++)
		{       
			int count=0;
			for(int i=-10;i<=10;i++)
			{
				for(int j=k;j<k+10;j++)
				{
					double dis=sqrt((double)(i*i+j*j));
					dis*=le;

					if(dis<re)
					{
						p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
						sumP+=p;
						count++;
					}
				}
			}
		}
		return 2*le*CON->get_sigma()/sumP;
	}

    else if(CON->get_dimention()==3) 
	{
		//////�����i�q�z�u�̏ꍇ
		if(calc_type==0)
		{
			//////�P���q�Ƃ̑��a
			if(CON->get_C_type()==0)
			{
				for(int k=1;k<10;k++)
				{       
					int count=0;
					for(int i=-10;i<=10;i++)
					{
						for(int j=-10;j<=10;j++)
						{
							double dis=sqrt((double)(i*i+j*j+k*k));
							dis*=le;

							if(dis<re)
							{
								p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
								sumP+=p;
								count++;
							}
						}
					}
				}
			}
			//////�ߓ���̕��@
			else if(CON->get_C_type()==1)
			{
				for(int k0=(-1);k0>-10;k0--)
				{
					for(int k=1;k<10;k++)
					{       
						int count=0;
						for(int i=-10;i<=10;i++)
						{
							for(int j=-10;j<=10;j++)
							{
								double dis=sqrt((double)(i*i+j*j+(k-k0)*(k-k0)));
								dis*=le;

								if(dis<re)
								{
									p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
									sumP+=p;
									count++;
								}
							}
						}
					}
				}
			}
			return 2*le*le*CON->get_sigma()/sumP;
		}

		//////MD�ɂ��z�u�̏ꍇ
		else if(calc_type==1)	
		{
			//�K�v�ϐ��̐錾
			double xa,ya,za,xb,yb,zb;
			double dx=1;
			double dy=sqrt(3.0)/2;
			double dz=sqrt(6.0)/3;
			int size=20;
			int countA=0;
			int countB=0;

			//////�P���q�Ƃ̑��a
			if(CON->get_C_type()==0)
			{
				for(int k=1;k<size;k++)
				{       
					int count=0;
					for(int i=-size;i<=size;i++)
					{
						for(int j=-size;j<=size;j++)
						{
							xa=i*dx;
							ya=j*dy;
							za=k*dz;
							if(j%2==1)	xa=i*dx+(dx/2);
							if(k%2==1)	ya=j*dy+(dy*2/3);

							double dis=sqrt((double)(xa*xa+ya*ya+za*za));
							dis*=le;

							if(dis<re)
							{
								p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
								sumP+=p;
								count++;
							}
						}
					}
				}
			}//*/

			//////�ߓ���̕��@
			if(CON->get_C_type()==1)
			{
				for(int k0=(-1);k0>-10;k0--)
				{
					for(int k=0;k<size;k++)
					{       
						int count=0;
						for(int i=-size;i<=size;i++)
						{
							for(int j=-size;j<=size;j++)
							{
								xa=i*dx;
								ya=j*dy;
								za=k*dz;
								zb=k0*dz;
								if(j%2==1)	xa=i*dx+(dx/2);
								if(k%2==1)	ya=j*dy+(dy*2/3);

								double dis=sqrt((double)(xa*xa+ya*ya+(za-zb)*(za-zb)));
								dis*=le;

								if(dis<re)
								{
									p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3;
									sumP+=p;
									count++;
								}
							}
						}
					}
				}
			}//*/

			//////�񕪊������e�����a���̑�������
			if(CON->get_C_type()==2)
			{
				for(int ka=0;ka<=size;ka++)
				{
					for(int ja=-size;ja<=size;ja++)
					{
						for(int ia=-size;ia<=size;ia++)
						{
							xa=ia*dx;
							ya=ja*dy;
							za=((double)ka+0.5)*dz;
							if(ja%2==1)	xa=ia*dx+(dx/2);
							if(ka%2==1)	ya=ja*dy+(dy*2/3);
							double disA=sqrt((double)(xa*xa+ya*ya+za*za));
							disA*=le;

							if(disA<re)
							{
								countA++;
								for(int kb=0;kb<=size;kb++)
								{
									for(int jb=-size;jb<=size;jb++)
									{
										for(int ib=-size;ib<=size;ib++)
										{
											xb=ib*dx;
											yb=jb*dy;
											zb=((double)kb+0.5)*(-dz);
											if(jb%2==1)	xb=ib*dx+(dx/2);
											if(kb%2==1)	yb=jb*dy+(dy*2/3);
											double disB=sqrt((double)(xb*xb+yb*yb+zb*zb));
											disB*=le;
											
											if(disB<re)
											{
												countB++;
												double dis=sqrt((double)((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya)+(zb-za)*(zb-za)));
												dis*=le;
												
												if(dis<re)
												{
													p=(dis-1.5*le+0.5*re)*(dis-re)*(dis-re)/3.0;
													sumP+=p;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}//*/
			return 2*le*le*sqrt(3.0)/4.0*CON->get_sigma()/sumP;
			//return 2*le*le*CON->get_sigma()/sumP;

			cout << "countA=" << countA << endl;
			cout << "countB=" << countB << endl;
		}
	}
	else return 0;

	return 0;	//�x���u�l��Ԃ��Ȃ��R���g���[���p�X������܂��v����������
}

//�\�ʒ��͏o�́i�ŐV�X�e�b�v�j
void plot_ST(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int t)
{
	int d=CON->get_dimention();
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le;
	//double times=CON->get_times()*le*le;
	double times=CON->get_times()*CON->get_density();

	ofstream st("ST.dat");
	

	if(d==2) for(int i=0;i<fluid_number;i++) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<potential[A_X][i]*times<<" "<<potential[A_Y][i]*times<<endl;
	else if(d==3) for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<potential[A_X][i]*times<<"\t"<<potential[A_Z][i]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
	}

	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	//if(CON->get_legend_F()>0) st<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��

	st.close();
}

//�\�ʒ��̓v���b�g�֐�(���X�e�b�v�o��)
void plot_ST_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int t)
{
	int d=CON->get_dimention();
	double xmax=-100;						//�o�͗��q�̍ő剡���W
	double ymax=-100;						//�o�͗��q�̍ő�c���W
	double le=CON->get_distancebp();
	//double times=CON->get_times()/CON->get_density()/le;
	//double times=CON->get_times()*le*le;
	double times=CON->get_times()*CON->get_density();

	char filename[20];
	sprintf_s(filename,"ST%d.dat", t);
	ofstream st(filename);

	if(d==2) for(int i=0;i<fluid_number;i++) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<potential[A_X][i]*times<<" "<<potential[A_Y][i]*times<<endl;
	
	else if(d==3) for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]<le*0.5 && PART[i].r[A_Y]>-le*0.5) st<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<potential[A_X][i]*times<<"\t"<<potential[A_Z][i]*times<<endl;
		if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
		if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
	}
	
	xmax+=4*le;//�}����o���ʒu��ی��������ď����΂߂Ɉړ�
	ymax+=4*le;
	//if(CON->get_legend_F()>0) st<<xmax<<" "<<ymax<<" "<<CON->get_legend_F()*times<<" "<<0*times<<endl;//�Ō�ɖ}��o��

	st.close();
}