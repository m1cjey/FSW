////新しい粒子を追加したときいじらないといけないのは�@u_laplacian_fのﾉﾝｽﾘｯﾌﾟ�AP_gradient�BP_gradient2�Ccalc_Temperatureのk[i]
////�D freeon
#ifndef mmmconfig
#define mmmconfig
class MMM_config
{       
	///////解析条件

	
	
	int Hf_type;	//外部磁場　0:一様磁場 1:磁石
	double Hf_H;			//外部磁場Hの強さ

	
	public:
	MMM_config();
	int		get_Hf_type()	{return Hf_type;}
	double get_Hf_H()		{return Hf_H;}
	
};


#endif
