////�V�������q��ǉ������Ƃ�������Ȃ��Ƃ����Ȃ��͇̂@u_laplacian_f���ݽد�߇AP_gradient�BP_gradient2�Ccalc_Temperature��k[i]
////�D freeon
#ifndef mmmconfig
#define mmmconfig
class MMM_config
{       
	///////��͏���

	
	
	int Hf_type;	//�O������@0:��l���� 1:����
	double Hf_H;			//�O������H�̋���

	
	public:
	MMM_config();
	int		get_Hf_type()	{return Hf_type;}
	double get_Hf_H()		{return Hf_H;}
	
};


#endif
