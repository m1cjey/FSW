#include "stdafx.h"	
#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B
#include"define.h"	//#define �i�[
#include"CONFIG.h"	//CONF.cpp�̓C���N���[�h���Ȃ��Ă��ACONFIG.h���C���N���[�h���邾���ł悢


void writedata(ofstream &fp, ofstream &gnu,int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM)
{	
	//�t�@�C���o��
	fp<<id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<P<<"\t"<<h<<"\t"<<val<<"\t"<<type<<"\t"<<materialID<<"\t"<<surface<<"\t"<<toBEM<<"\t"<<endl;
	
	//////gnuplot�p�ɏo��
	gnu<<x<<"\t"<<y<<"\t"<<z<<endl;
}


void set_initial_placement(mpsconfig *CON,int *particle_number)
{
	ofstream fp("initial_input.dat");//���q�f�[�^�i�[̧��
	ofstream gnu("plot.dat");		//gnuplot�p�����z�u�o��̧��
	
    double ii=0;
	double jj=0;
	double kk=0;
	int number=0;	//���q��
	int model=CON->get_model_number();
	double mass=CON->get_particle_mass();//���q�̎���
	

	//�Ǐd�݊֐��v�Z�p���f��
	if(model==0)
	{
		double T=293;//���x
		double h=T*mass*CON->get_Cp();//�G���^���s�[
		int materialID=1;
		double val=0;
		if(CON->get_dimention()==2)
		{
			//if(CON->get_wall_poly()==0)
			{
				//
				//////�eriw(�ǋ��E�Ɨ��q�̋���)�ł̏d�݂𒲂ׂ邽�߂̓_�B
				ii = 0;
				jj = 0;
				writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
				number++;


				//�ώ��ȕ�
				for(int i=-5;i<=5;i++)
				{
					for(int j=-5;j<=0;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
			}
			
		}
	}
	//�Ǖt�����`�t��
	if(model==1)
	{
		double T=293;//���x
		double h=T*mass*CON->get_Cp();//�G���^���s�[
		int materialID=1;
		double val=0;
		if(CON->get_dimention()==2)
		{
			//////����
			for(int i=-CON->get_fluidwidth()/2;i<=CON->get_fluidwidth()/2;i++)
			{
				
				for(int j=-CON->get_fluidwidth();j<=CON->get_fluidwidth();j++)
				{
					ii = i*CON->get_distancebp();
					jj = j*CON->get_distancebp();
					writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
					number++;
				}
			}

			if(CON->get_wall_poly()==0)
			{
				//////����
				//��
				for(int i=-CON->get_fluidwidth()/2-1;i<=-CON->get_fluidwidth()/2-1;i++)
				{
					for(int j=-CON->get_fluidwidth();j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//�E
				for(int i=CON->get_fluidwidth()/2+1;i<=CON->get_fluidwidth()/2+1;i++)
				{
					for(int j=-CON->get_fluidwidth();j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//��
				for(int i=-CON->get_fluidwidth()/2-1;i<=CON->get_fluidwidth()/2+1;i++)
				{
					for(int j=-CON->get_fluidwidth()-1;j<=-CON->get_fluidwidth()-1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//////�O��
				//��
				for(int i=-CON->get_fluidwidth()/2-5;i<=-CON->get_fluidwidth()/2-2;i++)
				{
					for(int j=-CON->get_fluidwidth()-5;j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//�E
				for(int i=CON->get_fluidwidth()/2+2;i<=CON->get_fluidwidth()/2+5;i++)
				{
					for(int j=-CON->get_fluidwidth()-5;j<=CON->get_fluidwidth()+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}

				//��
				for(int i=-CON->get_fluidwidth()/2-1;i<=CON->get_fluidwidth()/2+1;i++)
				{
					for(int j=-CON->get_fluidwidth()-5;j<=-CON->get_fluidwidth()-2;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
			}
			
		}
	}


	//�����̉t�w
	if(model==2)
	{
		double T=293;//���x
		double h=T*mass*CON->get_Cp();//�G���^���s�[
		int materialID=1;
		double val=0;
		if(CON->get_dimention()==3)
		{
			for(int i=-CON->get_fluidwidth()/2;i<=CON->get_fluidwidth()/2;i++)
			{
				for(int j=-CON->get_fluidwidth()/2;j<=CON->get_fluidwidth()/2;j++)
				{
					for(int k=-CON->get_fluidwidth()/2;k<=CON->get_fluidwidth()/2;k++)		
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						kk = k*CON->get_distancebp();	
						writedata( fp, gnu, number,  ii,  jj, kk+0.19,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
			}
		}
	}
	//////���f���R�@�~�`���H
	if(model==3)
	{
		double T=293;//���x
		double h=T*mass*CON->get_Cp();//�G���^���s�[
		int materialID=1;
		double val=0;
		//���̗��q��������
		for(int i=-CON->get_fluidwidth()/2;i<=CON->get_fluidwidth()/2;i++)
		{
			for(int j=-CON->get_fluidwidth()/2;j<=CON->get_fluidwidth()/2;j++)
			{
				if(CON->get_dimention()==2)
				{
					if(sqrt((double)(i*i+j*j))<CON->get_fluidwidth()/2)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();	
						writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
				else if(CON->get_dimention()==3)
				{
					for(int k=-CON->get_fluidwidth()/2;k<=CON->get_fluidwidth()/2;k++)
					{
			    		if(sqrt((double)(i*i+j*j+k*k))<CON->get_fluidwidth()/2)
						{
							ii = i*CON->get_distancebp();
							jj = j*CON->get_distancebp();
							kk = k*CON->get_distancebp();	
							writedata( fp, gnu, number,  ii,  jj, kk,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
							number++;
						}
					}
				}
			}
		}
	}
	//////*/
	
	/////////���f��1 �Ód����
	///�����@distancebp=0.000025 dt=0.000001 or0.0000005  fluidwidth=20
	if(model==14)
	{
		if(CON->get_dimention()==2)
		{
		    double le=CON->get_distancebp();
		    double val=0;
			int materialID=1;
		    int L=10;//��͗̈�
		    int W=CON->get_fluidwidth();
		    ///���H
		    for(int i=-W;i<=W;i++)
		    {
		        for(int j=0;j<=W;j++)
				{
					if(i*i+j*j<W*W)
					{
						ii = i*le;
			    		jj = j*le;
			    		double T=293;//���x
			    		double h=T*mass*CON->get_Cp();//�G���^���s�[
						writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
						number++;
					}
				}
		    }
		    
		    ////�~���d��
		    for(int i=-W*2;i<=W*2;i++)
		    {
		        for(int j=-1;j>=-W/3;j--)
				{
		            ii = i*le;
					jj = j*le;
					double T=293;//���x
					double h=T*mass*CON->get_Cp();//�G���^���s�[
					if(i*i<=W*W)
					{
						if(j==-1 ||i==W || i==-(W))
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
							number++;
						}
						else
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, h, 0);
							number++;
						}
					}
			  
				}
			}
		}
		/*else if(CON->get_dimention()==3)
		{
		    double le=CON->get_distancebp();
			double mass=CON->get_particle_mass();
		    int W=CON->get_fluidwidth();
		    ///���H
			 for(int i=-W;i<=W;i++)
		    {
		        for(int j=-W;j<=W;j++)
					{
					for(int k=0;k<=W;k++)
					{
						if(i*i+j*j+k*k<W*W)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
		    				double T=293;//���x
				 			double h=T*mass*CON->get_Cp();//�G���^���s�[
							writedata(fp, number, ii, jj,kk, FRFLUID,gnu,0,0,0,0,h,0);
							number++;
						}
					}
				}
			}

			///�㋅
			for(int i=-W;i<=W;i++)
		    {
		        for(int j=-W;j<=W;j++)
				{
					for(int k=0;k<=W;k++)
					{
						if(i*i+j*j+k*k<3)
						{
							ii = i*le;
			    			jj = j*le;
							kk = k*le+1.2*W*le;
			    			double T=293;//���x
			    			double h=T*mass*CON->get_Cp();//�G���^���s�[
							//writedata(fp, number, ii, jj,kk, FRFLUID,gnu,0,0,0,0,h,0);
							//number++;
						}
					}
				}
		    }//////
	    
		    ////�~���d��(INWALL)
		    for(int i=-2*W;i<=2*W;i++)
		    {
		        for(int j=-2*W;j<=2*W;j++)
				{
		            for(int k=-1;k>=-5;k--)
					{
		                ii = i*le;
						jj = j*le;
						kk = k*le;
						double T=293;//���x
						double h=T*mass*CON->get_Cp();//�G���^���s�[
						int r=i*i+j*j;
						//if(r<=(W+1)*(W+1))//1�w�����d�ɁB���H�����ɂ��ڂ��̂�h���H
						if(r<=W*W)
						{
							if(k>=-2 || r>=(W-1)*(W-1))
							{
								writedata(fp, number, ii, jj,kk, INWALL,gnu,0,0,0,0,h,0);
								number++;
							}
						}
					}	  
				}
			}
			////�~���d��(OUTWALL)
			for(int i=-2*W;i<=2*W;i++)
			{
			    for(int j=-2*W;j<=2*W;j++)
				{
		         for(int k=-1;k>=-5;k--)
					{
			            ii = i*le;
						jj = j*le;
						kk = k*le;
						double T=293;//���x
						double h=T*mass*CON->get_Cp();//�G���^���s�[
						int r=i*i+j*j;
						if(r<(W-1)*(W-1))
						{
							if(k<-2)
							{
								writedata(fp, number, ii, jj,kk, OUTWALL,gnu,0,0,0,0,h,0);
								number++;
							}
						}
					}
			  
				}
			}
		}*/
	}
	///////////////////*/

	////���f��19�@FSW
	//���a9mm�A�[��6mm�̉t���ɁA���a2.5mm�A�[��4mm����۰�ޑ}���B�V�����_�[���a��6mm
	if(model==19)
	{
		if(CON->get_dimention()==3)
		{
			double le=CON->get_distancebp();
			double mass=CON->get_particle_mass();
			int W=CON->get_fluidwidth();	//���̂̂w��
			int YW=2*W;						//���̂̂x��
			int fH=W*2/3;					//���̍���
			int WH=W*2/3;					//�ǂ̍���
			int shoR=W*2/3;					//�V�����_�[���a
			double proR=2.5e-3;				//��۰�ޔ��a
			double proL=4e-3;				//��۰�ޒ���
			double rpm=500;
			double rps=rpm/60;
			double w=rps*2*PI;				//�p���x
			
			double U=CON->get_move_speed();	//�v���[�u�̈ړ����x[m/sec]		
			double pich=0.7e-3;				//�v���[�u�̂˂��̃s�b�`0.7[mm]
			double T=CON->get_roomT();		//�������x
			double wall_h=T*(CON->get_wall_density()*le*le*le)*CON->get_wall_Cp();//�Ǘ��q�̃G���^���s�[

			int toBEM=OFF;//���̃��f���ł͂e�d�l�͎g�p���Ȃ�����ǂ���ł��悢
			double val=0;
			int materialID=1;

			//calclation type
			int plunge=0;
			int traverse=1;
			int calc_type=0;
			if(CON->get_process_type()==2)	calc_type=0;	//1:plunge 2:traverse
			else
			{
				calc_type=CON->get_process_type();
			}
			double plunge_H=4*1e-3+1*le;	//plunge��͂̍ۂ́A�c�[�����グ�鍂��

			//���� ���a9mm�A�[��6mm
			for(int i=-W;i<=W;i++)
			{
				for(int j=-W;j<=YW;j++)
				{
				
					for(int k=1;k<=fH;k++)
					{
						int flag=0;
						ii = i*le;
						jj = j*le;
						kk = k*le;
						double r=(ii*ii+jj*jj);
						if(calc_type==traverse)
						{
							if(r<proR*proR && kk>=fH*le-proL) flag=1;//��۰�ޗ̈� plunge��͂̂Ƃ��͂��ׂė��̂ŗǂ�
						}
						if(flag==0)
						{
							if(ii>0) materialID=1;
							else materialID=2;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
							writedata( fp, gnu, number,  ii,  jj, kk,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, toBEM);
							number++;
						}
					}
				}	
			}
			materialID=1;
			//��INWALL����
			for(int i=-W-2;i<=W+2;i++)
			{
				for(int j=-W-2;j<=YW+2;j++)
				{
					if(i*i>W*W || j<-W || j>YW)
					{
						for(int k=1;k<=fH;k++)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
							writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
							number++;
						}
					}
				}
			}
			//��INWALL��
			for(int i=-W-2;i<=W+2;i++)
			{
				for(int j=-W-2;j<=YW+2;j++)
				{
					for(int k=0;k>=-1;k--)
					{
						ii = i*le;
						jj = j*le;
						kk = k*le;
						double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
						writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
						number++;
					}
				}
			}
				//�v���[�uINWALL
			for(int i=-W;i<=W;i++)
			{
				for(int j=-W;j<=W;j++)
				{
					if(i*i+j*j<W*W)
					{
						for(int k=1;k<=fH;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							if(r<proR*proR && kk>=fH*le-proL) flag=1;//�v���[�u
							if(r<(proR-le)*(proR-le) && kk>=fH*le-proL+le) flag=0;
							if(flag==1)
							{
								double speed=0;
								double u=0;double v=0;double uw=0;
								if(r!=0)
								{
									speed=sqrt(r)*w;
									u=speed*(-jj/sqrt(r));
									v=speed*(ii/sqrt(r));
									if(r>=(proR-le)*(proR-le))
									{
										uw=-pich*rps;
									}
								}
								if(calc_type==plunge)
								{
									kk+=plunge_H;
									uw-=U;//Z�����ɑ}�����x����悹
								}
								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
								writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  u, v, uw, 0, wall_h, MOVE);
								number++;
							}
						}
					}
				}
			}
			//�V�����_�[INWALL
			for(int i=-shoR;i<=shoR;i++)
			{
				for(int j=-shoR;j<=shoR;j++)
				{
					if(i*i+j*j<=shoR*shoR)
					{
						for(int k=fH+1;k<=fH+2;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							
							double speed=0;
							double u=0;double v=0; double uw=0;
							if(r!=0)
							{
								speed=sqrt(r)*w;
								u=speed*(-jj/sqrt(r));
								v=speed*(ii/sqrt(r));
							}
							if(r>=(proR-le)*(proR-le) && r<proR*proR) uw=-pich*rps;//�ǂɂ����鑬�x���U�[�����̂��߁A�V�����_�[���ɂ��˂����x����
							if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
							if(calc_type==plunge)
							{
								kk+=plunge_H;
								uw-=U;//Z�����ɑ}�����x����悹
							}
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
							writedata( fp, gnu, number,  ii,  jj, kk,  INWALL,materialID, OFF, val,  u, v, uw, 0, wall_h, MOVE);
							number++;
						}
					}
				}
		    }
		    //OUTWALL����
			for(int i=-W-4;i<=W+4;i++)
		    {
				for(int j=-W-4;j<=YW+4;j++)
				{
					if(i*i>(W+2)*(W+2) || j<-W-2 || j>YW+2)
					{
						for(int k=-1;k<=fH;k++)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
							writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
							number++;
						}
					}
				}
			}
			//OUTWALL����
			for(int i=-W-4;i<=W+4;i++)
			{
				for(int j=-W-4;j<=YW+4;j++)
				{
					//if(i*i+j*j<=(W+4)*(W+4))
					{
						for(int k=-2;k>=-3;k--)
						{
							ii = i*le;
							jj = j*le;
							kk = k*le;
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
							writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, toBEM);
							number++;
						}
					}
				}
			}
			//�v���[�uOUTWALL
			for(int i=-W;i<=W;i++)
			{
				for(int j=-W;j<=W;j++)
				{	
					if(i*i+j*j<W*W)
					{
						for(int k=1;k<=fH;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							if(r<(proR-le)*(proR-le) && kk>=fH*le-proL+le) flag=1;
							if(flag==1)
							{
								double speed=0;
								double u=0;double v=0; double uw=0;
								if(r!=0)
								{
									speed=sqrt(r)*w;
									u=speed*(-jj/sqrt(r));
									v=speed*(ii/sqrt(r));
									//uw=-pich*rps;
								}
								if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
								if(calc_type==plunge)
								{
									kk+=plunge_H;
									uw-=U;//Z�����ɑ}�����x����悹
								}
								double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
								writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val, u, v, uw, 0, wall_h, MOVE);
								number++;
							}
						}
					}
				}
		    }
			//�V�����_�[OUTWALL
			for(int i=-shoR;i<=shoR;i++)
		    {
				for(int j=-shoR;j<=shoR;j++)
				{
					if(i*i+j*j<=shoR*shoR)
					{
						for(int k=fH+3;k<=fH+20;k++)
						{
							int flag=0;
							ii = i*le;
							jj = j*le;
							kk = k*le;
							
							double r=(ii*ii+jj*jj);
							
							double speed=0;
							double u=0;double v=0; double uw=0;
							if(r!=0)
							{
								speed=sqrt(r)*w;
								u=speed*(-jj/sqrt(r));
								v=speed*(ii/sqrt(r));
							}
							if(r>=(proR-le)*(proR-le) && r<proR*proR) uw=-pich*rps;//�ǂɂ����鑬�x���U�[�����̂��߁A�V�����_�[���ɂ��˂����x����

							if(calc_type==traverse) v+=U;//y�����Ƀc�[�����ړ�
							if(calc_type==plunge)
							{
								kk+=plunge_H;
								uw-=U;//Z�����ɑ}�����x����悹
							}
							double h=T*mass*CON->get_Cp()+CON->get_latent_H()*mass;//�G���^���s�[
							writedata( fp, gnu, number,  ii,  jj, kk,  OUTWALL,materialID, OFF, val, u, v, uw, 0, wall_h, MOVE);
							number++;
						}
					}
				}
			}
		}
	}

	//ih��
	if(model==26)
	{
		double le=CON->get_distancebp();
		double R=CON->get_fluidwidth()*0.001;//3;			//
		double Zg=CON->get_height();					//
		double C_h=0.1;//10;							//�����ʕ��̍���
		double C_hout=C_h+4*le;//10;							//�����ʕ��̍���
		double C_R=0.08;//3;			//���ꕔ�̔��a
		double C_Rout=C_R+4*le;
		double real_width=C_R/le;
		double real_height=C_h/le/2;

		int width=(int) real_width;
		int height=(int) real_height;

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
		
		double air_mass=V*1.205;//��C�̏d��
		double air_Cp=1006;//��C�̔�M
		double air_h=T*air_mass*air_Cp;//�ő�

		
		double val=0;

		if(CON->get_dimention()==2)
		{
			//////����
			for(int i=-width;i<=width;i++)
			{
				
				for(int j=-height;j<=height;j++)
				{
					ii = i*CON->get_distancebp();
					jj = j*CON->get_distancebp();
					if(jj<=0.0)//0.04+0.01��50mm�B���a80mm�̉~���e���1l�̂��̂���ꂽ�ۂ̂����悻�̍���
					{
					writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
					//writedata2(fp,i,ii,j,Z[i],FLUID,materialID,1,0,0,0,0,0,h,1);
					number++;
					}
				}
			}

			if(CON->get_wall_poly()==0)
			{
				//////����
				//��
				for(int i=-width-1;i<=-width-1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=-height;j<=height+1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						//if(jj<=0.0)//0.04+0.01��50mm�B���a80mm�̉~���e���1l�̂��̂���ꂽ�ۂ̂����悻�̍���
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}

				//�E
				for(int i=width+1;i<=width+1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=-height;j<=height+1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						//if(jj<=0.0)//0.04+0.01��50mm�B���a80mm�̉~���e���1l�̂��̂���ꂽ�ۂ̂����悻�̍���
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}

				//��
				for(int i=-width-1;i<=width+1;i++)
				{
					for(int j=-height-1;j<=-height-1;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}

				//////
				if(CON->get_airwall()==1)
				{
					//��//�\�ʗ�����}����p
					for(int i=-width;i<=width;i++)
					{
						for(int j=height+1;j<=height+1;j++)
						{
							ii = i*CON->get_distancebp();
							jj = j*CON->get_distancebp();
							writedata( fp, gnu, number,  ii,  jj-0.05+3*CON->get_distancebp(), 0,  INWALL,materialID, OFF, val,  0, 0, 0, 0, air_h, 0);//+3*CON->get_distancebp()
							number++;
						}
					}
					////*/
				}

				//////�O��
				//��
				for(int i=-width-5;i<=-width-2;i++)
				{
					for(int j=-height-5;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}
				/*/////
				for(int i=-width-1;i<=-width-1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=height+1;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						if(jj>0.0)//0.04+0.01��50mm�B���a80mm�̉~���e���1l�̂��̂���ꂽ�ۂ̂����悻�̍���
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}
				/////*/

				//�E
				for(int i=width+2;i<=width+5;i++)
				{
					for(int j=-height-5;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}
				/*/////
				//�E
				for(int i=width+1;i<=width+1;i++)
				{
					for(int j=-height;j<=height+5;j++)
					//for(int j=height+1;j<=height+5;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						if(jj>0.0)//0.04+0.01��50mm�B���a80mm�̉~���e���1l�̂��̂���ꂽ�ۂ̂����悻�̍���
						{
							writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
							number++;
						}
					}
				}
				////*/

				//��
				for(int i=-width-1;i<=width+1;i++)
				{
					for(int j=-height-5;j<=-height-2;j++)
					{
						ii = i*CON->get_distancebp();
						jj = j*CON->get_distancebp();
						writedata( fp, gnu, number,  ii,  jj, 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, wall_h, 0);
						number++;
					}
				}
				if(CON->get_airwall()==ON)
				{
					//////
					//��
					for(int i=-width;i<=width;i++)
					{
						for(int j=height+2;j<=height+5;j++)
						{
							ii = i*CON->get_distancebp();
							jj = j*CON->get_distancebp();
							writedata( fp, gnu, number,  ii,  jj-0.05+3*CON->get_distancebp(), 0,  OUTWALL,materialID, OFF, val,  0, 0, 0, 0, air_h, 0);//+3*CON->get_distancebp()
							number++;
						}
					}
					//////*/
				}
			}
			
		}
	}

	//�n��
	if(model==27)
	{
		double le=CON->get_distancebp();
		double R=CON->get_fluidwidth()*0.001;//3;			//

		double real_width=0.025/le/2;
		double real_height=0.01/le/2;

		int width=(int) real_width;
		int height=(int) real_height;

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
		

		
		double val=0;

		if(CON->get_dimention()==2)
		{
			//////����
			for(int i=-width;i<=width;i++)
			{
				
				for(int j=-height;j<=height;j++)
				{
					ii = i*CON->get_distancebp();
					jj = j*CON->get_distancebp();
					
					writedata( fp, gnu, number,  ii,  jj, 0,  FLUID,materialID, OFF, val,  0, 0, 0, 0, h, 0);
					//writedata2(fp,i,ii,j,Z[i],FLUID,materialID,1,0,0,0,0,0,h,1);
					number++;
					
				}
			}

		}
	}

	
	///////////////////*/

	fp.close();
	gnu.close();
	*particle_number=number;
}