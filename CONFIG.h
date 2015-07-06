////�V�������q��ǉ������Ƃ�������Ȃ��Ƃ����Ȃ��͇̂@u_laplacian_f���ݽد�߇AP_gradient�BP_gradient2�Ccalc_Temperature��k[i]
////�D freeon
#ifndef config
#define config

#include"header.h"	//��v�ȃw�b�_�[�t�@�C���͂܂Ƃ߂Ă��̂Ȃ��B

class mpsconfig
{       
	///////��͏���

	double dt;//���ԍ��ݕ� 
	int    step;
	int    dimention;//����
	double maxX;//��͗̈�
	double minX;
	double maxY;
	double minY;
	double maxZ;
	double minZ;

	int material;
	
	///////����1�̕����l
	
	double density;//���x
	double nensei;//�S���W��
	double sigma;//�\�ʒ��͌W��
	double Cp;//�舳��M
	double k;//�M�`����
	double latent_H;//���M(latent_heat)
	double MP;//�Z�_(melting point)
	double CTE;//�M�c���W��[1/K] 

	int temperature_depend;	//�����l�����x�ˑ��̕ω������邩�ǂ���

	///////����2�̕����l
	
	double density2;//���x
	double nensei2;//�S���W��
	double sigma2;//�\�ʒ��͌W��
	double Cp2;//�舳��M
	double k2;//�M�`����
	double latent_H2;//���M(latent_heat)
	double MP2;//�Z�_(melting point)
	double CTE2;//�M�c���W��[1/K] 
	
	////////���q�z�u�p
	
	int fluidwidth;//���̂̑�\����
	int fluid_h;//���̂̑�\����
	double distancebp;//�������q�ԋ���
	double wlength;//���E�̕ǂ̋����B���̂̉��{��
	double height;//���̂̏d�S����
	double tool_angle;//FSW�ɂ����āA�c�[�����X����p�x
	int tool_type;
	int process_type;
	double dwelling_time;
	int airwall;
	int change_step;
	
	///////���q�@�p�p�����[�^
	
	double re;		//��ʓI�ȗ��q���are
	double re2;		//���u���V�A���p��re
	double re3;		//�\�ʒ��͗p��re
	double re4;		//freeon
	double beta;		//��
	int    dx;		//index�p�i�q���B�ő嗱�q���a���܂ސ����ɂ��邱��
	double times;		//�޸����ۯėp�̔{��
	
	////////�\�ʒ��͊֌W
	
	int    surface_tension; //�\�ʒ��̓X�C�b�`�@0=OFF 1=���ȏ� 2=���q�Ԉ���
	int    smooth;			//�X���[�W���O�X�C�b�`�@1=ON 0=OFF
	int    SFT;				//1=ON 0=OFF �\�ʒ���(surface tension)�̉��x(T)�ˑ����X�C�b�`
	int    smn;				//�Ѱ��ݸނ̉�
	int    smth_sumP;       //for ver3. curv�̽Ѱ��ݸ� 0=OFF 1=ON
    int    wall_direct;		//for ver2. BDWALL�̖@���޸�ٌv�Z���@�@0=���� 1=�ʏ� 2=�������� 3=��������
	int    non_slip;		//for ver.5 �S������ ON�Ȃ�BDWALL��ɉ��z�I�ɉt�̂����݂���ݒ�B���̍ۂ�BDWALL�̔z�u�ɗ���
	int    suf_modify;		//�ȗ��v�Z�œ���ꂽ�ȖʂɈ�v����悤�ɗ��q�ʒu���C�����邩�A���Ȃ���
	int    interpolate_curv;//�ȗ����v�Z����Ȃ��������q�̋ȗ����A���ӗ��q�̒l������
	double Cst;				//�|�e���V�����W��
	int    C_type;			//�|�e���V�����W���̌�����@
	double C_times;			//�|�e���V�����W���̔{��
	double wall_C;			//�ǂƂ̐e�a�͌W��
	
	////////���͊֌W

	int iteration_count;//���͌v�Z���s���� �ʏ��1�ɐݒ�
	int    solution;	//�A���������̉���@0=CG method 1=ICCG method
	int    B_of_P;  	//���͌v�Z�̉��s�� 0=���ȏ� 1=���x���U 2=0+1
	double w_div;		//���U�̏d�� 
	int    divU_method;	//���x���U�̌v�Z���@ 1=MPS,2=WLSM,3=���g��ʂ�Ȃ�WLSM 4=������
	int    divU_order;	//WSLM��divU���v�Z����Ƃ��̃I�[�_�[ MAX=3
	int    HL_sw;		//���׼�݂ɑ΂������̗��U�����s�� 1=ON 0=OFF
	int	   dir_for_P;	//���͌v�Z�̍ہA�\�ʂɔ�[���̃f�B���N���l�������邩�A���Ȃ���
	int    initialP;	//���͌v�Z���ɏ����l�Ƃ��Č����͂������邩�ǂ����@1=ON 0=OFF
	int    P_twice;
	double CGep;		//���͌v�Z���ɂ������������
	int	   omp_P;		//���͌v�Z����CG�@��openMP���g�p���邩 1=ON 0=OFF
	int    negativeP;	//�����@1=���� 0=�s����
	int    set_P;		//���͋�������X�C�b�` 1=ON 0=OFF
	int    Pgrad;		//���͌��z�v�Z�@ 0=���ȏ� 1=�\�ʂ̂ݖ@��
	int    minP;		////�ŏ����ͽ��� 1=ON 0=OFF
	int    ave_Pdirect;	//���͌��z�p�@���޸�ق̕��ω� 1=ON 0=OFF
	int    Pgrad_order;	//���͌��zver.4�ɂ����鐸�x 1=���` 2=2��
	int	   artP_sw;		//�l�H���͌v�Z�t���O 0=OFF 1=artP���� 2=Monaghan
	double artP;		//�l�H���͒l[Pa]
	double Pgrad_times;	//���͌��z����ۯĂ���ۂ́A�\�ʒ��͂ɑ΂���{�� �ʏ��1�ɐݒ�
	int    P_AVS;		////microAVS�p�̈��̓t�@�C�����o�͂���step�Ԋu�B0�Ȃ�OFF
	double Pgrad_limit;
	int    gridless_P;	//�O���b�h���X�@�ɂ�鈳�͌v�Z���s����(ON)�AMPS�ɂ��v�Z���s����(OFF)
	int	   wall_poly;
	int    interpolate_surface_u; //���͌��z�ɂ�鑬�x�C����̕\�ʗ��q�̑��x��������q�̑��x�ŕ⊮���邩���Ȃ����B

	////////���x��֌W
	
	int    T_field;		//���x���́@�@1=ON 0=OFF
	int    insulate;	//�ǂƂ̒f�M��� 0=�f�M�@1=��f�M
	int    T_laplacian;	//���x�����׼�݁B�@0=���ȏ��@1=��[i] 2=���U�E���z
	double wall_density;		//�ǂ̖��x[kg/m^3]
	double wall_Cp;			//�ǂ̔�M[J/(kg�EK)]
	double wall_k;				//�ǂ̔M�`����[W/(m�EK)]
	double    initialT;		//�������x
	double    roomT;		//���� int�^�ł悢
	int    air_cool;	//��C�Ƃ̔M�`�B���l�����邩���Ȃ��� 1=ON 0=OFF
	int    T_expansion;	//�M�c���v�Z�@1=ON 0=OFF
	int    T_CG;
	int    T_solution;
	double T_CGep;
	int    buoyant;		//����(���x���޼�Ƚ��ߎ�)  1=ON 0=OFF
	double TplotZ;		//3D��͂ɂ����āAXY���ʂ̉��x���o�͂���Ƃ��̂y���W
	int    T_AVS;		//microAVS�p�̉��x�t�@�C�����o�͂���step�Ԋu�B0�Ȃ�OFF
	int output_temperature_face;
	int output_equivalent_strain_rate_face;
	int output_forward;
	int output_backward;
	int output_another_face;

	////////�d����̉�@
	int		EM_method;			//�d����̉�@ 0=OFF 1=FEM 2=BEM 3=���CӰ��Ė@
	int		EM_calc_type;		//0=OFF 1=�d�� 2=���� 3=����(�Q�d��)
	int		EM_interval;    //�d����v�Z�����ï�߂Ɉ��s�����B�ʏ��1�ɐݒ�
	double  EM_distance;	//�d���͌v�Z���s���p�x���A���q�̑��ړ������Ŕ��肷�� 0:OFF ���̒l*le���ړ���������������FEM�����s
	int    region_shape;	//��͗̈�`��@0=������ 1=�~��
	double XL;				//��͗̈捶�[
	double XR;              //��͗̈�E�[
	double YU;              //��͗̈��[
	double YD;              //��͗̈扺�[
	double ZU;				//��͗̈�Z������[
	double ZD;				//��͗̈�Z�������[
	double RU;				//��͗̈悪�~���`�ƂȂ�Ƃ��̂��̔��a

	////////���b�V��
	int    mesher;			//���b�V���쐬���@ 0:�g�삳��̃f���[�j�����v���O����
	int    mesh_input;		//mesh�̍쐬���@ 0:���� 1:Magnet
	int    remesh_sw;		//ON:remesh�̈���l�����A��������remesh OFF:FULL�̈�𖈉񕪊�
	double boxalpha;		//���߰�ޯ���̑傫�� 1�`10�̎����ɂ��邱�ƁB�ʏ��2�ɐݒ�
	int    fine;			//�ߓ_�̎����ǉ����s�����ǂ��� 0:OFF 1:ON
	double co_fine;			//�ߓ_�����ǉ��ɂ�����W��.�傫���قǕ��̂��痣�ꂽ��e���Ȃ�
	int    add_points;		//�����ߓ_�ǉ���
	int    poly_memory;		//poly3D()�Ŋm�ۂ��郁�����̐�
	int	   CDT_sw;
	int	   defer_f;
	int	   mesh_smth;
	double material_judge;
	double del_length;
	int	   BD_fluid;
	int    output_wall_way;	//�Ǘ��q���o�͂��ă��b�V�����쐬���邩�ǂ����@0=���̂̂� 1=���̂ƕǂ��ׂā@2=���̂�INWALL�̂�
	double layer_thin;			//��C�w�̌��� distancebp*layer_thin�����ۂ̌���
	int    layer_num;			//��C�w�̑w��

	////////FEM�֌W
	int	   BEM_elm_type;		//�v�f���� 0:�ߓ_�v�f 1:�ӗv�f
	int	   FEM_elm_type;	//�v�f���� 0:�ߓ_�v�f 1:�ӗv�f
	int    FEM_smn;			//�d���ͽѰ��ݸމ񐔁@0�Ȃ�OFF
	int    max_DN;			//�ިظڌ^���E�������Ƃ�ő�ߓ_��
	int	   FEMCG;			//FEM�ɂ�����s���@ 0:CG 1:ICCG
	double CGaccl;			//CG,ICCG�@�ɂ���������t�@�N�^�@CGaccelerator
	double FEMCGep;			//ICCG�̎�������
	double FEMtimes;		//�d���͂���ۯĂ���ۂ́A�\�ʒ��͂ɑ΂���{�� �ʏ��1�ɐݒ�
	double legend_F;		//F.dat�̖}��ɏo�͂����[N]
	int	   node_sort;
	int	   thinout_fluid;
	
	///////�d�E�v�Z
	double V;		//�d��
	int    V_step;		//�d������J�n�X�e�b�v
	double r_perm;		//��U�d�� (relative permittivity)
	int    V_con;		//�d�������@0:��ٽ�@1:�Ʊ�@2:���萔
	double initial_V;	//�d�������Ʊ�̂Ƃ��̏����d��
	double E_times;		//̧�ُo�͂���ۂ́A�d�E�̔{��
	int    eleforce;	//�Ód�͌v�Z���@ 1:�ߓ_�͖@ 2:�ϕ���
	int	   charge;		//�d�ׂ��l�����邩���Ȃ����B0=OFF 1=ON
	
	///////����v�Z
	
	int    J_input_way;	//�d�����x������@ 0:���� 1:�\�t�g
	double J0;//�����d�����x[A/m2]
	double I0;//�����d���l[A]
	double RP;//�䓧����(Relative Permeability)
	double ele_conduc;//�d�C�`����
	double ele_conduc2;//�d�C�`����
	double  Hz;///�𗬂̎��g��
	int    div_Hz;//�P�����̕�����(��͐��x)
	int    jheat;		//�Q�d���ɂ�锭�M���l�����邩�@0=OFF 1=ON
	int	   m_force;		//�d���͌v�Z���� 0=ϸ���ق̉��� 1=�̐ϗ�
	int    NLBHI;		//�̐ϗ͂ɂ����āA�v�f�a����v�f�g�����߂�ۂɔ���`�����l�����邩�A���Ȃ���(non linier B H inverter)
	int    NLMH;		//�l�̎Z�o�ɔ���`�����l�����邩�A���Ȃ���
	double magnet_H;	//�i�v���΂̑傫��
	double magnet_B;	//�i�v���΂̋���[T]
	int    uniform_B_sw;	//��͗̈撆�Ɉ�l����𔭐������邩�ۂ� 0=OFF 1=ON
	double uniform_B;	//��l����̑傫��[T]
	double B_times;		//̧�ُo�͂���ۂ́A�������x�̔{��
	int    plot_B_type;	//�������x�o�̓^�C�v 1=�޸�� 2=�X�J���[
	int	   Je_crucible; //�Q�d���̌v�Z�Ώہ@0=�n�Z���� 1=�n�Z�����Ƃ�ځi�R�C�����j
	int	   J0eqJe;
	int    m_A;
	int    jw_interval;
	int    jw_Faverage;
	int    static_dirichlet;
	double eddy_times;
	int    A_phi;
	int    parabolic_node_element;
	
	////////�e��X�C�b�`�@�ǂ̂悤�Ȍv�Z���l�����邩�A���Ȃ���
	
	double g;		//�d�͉����x�@�ʏ��-9.8�Ƃ��邱��
	int    restart;		//restart�p�X�C�b�`
	int    autosave;    //��ľ��ފԊu�B�����ɂ������Ƃ��͑傫�Ȑ����������Ă���
	double curan;		//�N�[���������� 0�Ȃ�OFF
	int    modify_position;//���q�ԋ�����le��菬�����ꍇ������C�����邩�A���Ȃ���
	int    vis_calc_type;//�S�����v�Z��@�@0=�z��@ 1=�A��@
	int    wall_adheision;	//�ǂ̔S����ԁ@1=�ݽد�� 0=�ذ�د��
	int    laplacian;	//���׼�݁B�@0=���ȏ��@1=��[i] 2=���U�E���z
	int	   vis_solver;  //�S�������A��͂ŉ����ۂ̍s����ް 0:CG 1:ICCG
	int    initial_u;	//�S�������A�I�ɉ����ۂɁA�����l�Ƃ��Č��ݑ��x����͂��邩�A���Ȃ���
	int    temporary_r; //�z��͌�ɉ��̈ʒu���v�Z���邩���Ȃ����B 1=ON 0=OFF �ʏ��ON�ɐݒ�
	
	int    fix_center;	//1=ON 0=OFF
	int    freeon;		//���q�ˑ��֌W�֐��@1:���񉻉\ 2:����s��
	int    freeon3sw;	//freeon3���v�Z���邩���Ȃ��� 1=ON 0=OFF
	int    surface_judge2;//surface_judge2()���g�p���邩���Ȃ���  1=ON 0=OFF
	int    move_prtcl;	//�ړ����q���l�����邩���Ȃ��� 1=ON 0=OFF
	int    move_u_dirct;//�ړ����q���ړ���������@���݂́}X����=�}1,�}Y����=�}2,�}Z����=�}3
	double move_speed;	//�ړ����q�̈ړ����x[m/s]
	double move_speed2; //�ړ����q�̈ړ����x[m/s]
	int	   check_region;//�����𖞂��������q���폜���邩���Ȃ���
	double max_speed;	//chec_region��ON�ɂȂ��Ă���Ƃ��A���q�̑��x�����̒l�ȏ�Ȃ�폜����
	int    check_something;//check_something()�����s���邩���Ȃ��� 1=ON 0=OFF
	int	   AVS_HALF;
	int    adaptive_sw;	//�𑜓x���ςɂ��邩�A���Ȃ��� 1=ON 0=OFF
	double threshold;	//�A�_�v�e�B�u�ɂ���ۂ́A���͌덷��臒l
	int    fix_surface;
	int output_viscousity_face;
	
	int    model_number;
	int    model_set_way;		//model���Z�b�g������@�@0=�����i�q 1=MD
	int    model_inherit;		//�O�̉�͂Ŏg�������q���f���������p��
	int	   ea_model;

	///////���x��ۯĕϐ�
	int    speed_plot_particle;	//���x����ۯĂ��闱�q�̎�� 1=���ׂ� 2=fluid
	double speedtimes;			//���x��ۯĎ��́A���W�ɑ΂��鑬�x�̔{��
	int    speed_face;			//speed.dat�̏o�͖� 0=YZ���� 1=XZ 2=XY
	double speed_face_p;		//3D��͎���speed.dat�̏o�͖ʂ̍��W
	int    ax_sym_modify;		//3D����speed.dat�Ɋւ��āA���Ώ̂ɂ��o�͏C�����s�����ۂ��@1=ON 0=OF
	int    flat_speed_plot;		//���������̑��x����ۯĂ��邩���Ȃ���
	double flat_speed_p;		//flat_speed.dat�̏o�͖ʂ̍��W
	int	   relative_speed;		//�d�S�ɑ΂��鑊�Α��x���o�͂��邩���Ȃ��� 1=ON 0=OFF
	int    speed_AVS;			//microAVS�ɂ��3D���x���z�o�͂��邩���Ȃ��� 1=ON 0=OFF
	double legend_speed;		//speed.dat�̖}��ɏo�͂��鑬�x[m/s]
	int    set_zero_speed;//restart���ɑ��x���[���Z�b�g���邩���Ȃ���  1=ON 0=OFF

	///////GPU�֌W
	int    M_form;				//CG�@�ɂ�����W���s��̊i�[���@ CSR_scl,CSR_vec,ELL��3���߰�
	int    MAX_thread;			//�ЂƂ�SM������̍ő�گ�ސ��@�ӂ���512
	
	///////̧�ُo�͕ϐ�
	
	int    interval;//AVS�p�X�e�b�v�Ԋu
	int    AVS;	//0:���ʁ@1:���́@2:���x
	double maxT;    //AVS(2)�ɂ�����ő剷�x
	double minT;    //AVS(2)�ɂ�����ŏ����x
	int P_interval;
	int ST_interval;
	int PND_interval;
	int mesh_output_interval;
	int plot_nearbywall;
	
	public:
	mpsconfig();
	double get_dt()		{return dt;}
	int get_step()	{return step;}
	int    get_dimention()	{return dimention;}
	double get_maxX()	{return maxX;}
	double get_minX()	{return minX;}
	double get_maxY()	{return maxY;}
	double get_minY()	{return minY;}
	double get_maxZ()	{return maxZ;}
	double get_minZ()	{return minZ;}

	int get_material() {return material;}
	
	double get_density()	{return density;}
	double get_nensei()	{return nensei;}
	double get_sigma()      {return sigma;}
	double get_vis()	{return nensei/density;}//���S���W��
	double get_Cp() 	{return Cp;}
	double get_k()		{return k;}
	double get_latent_H()	{return latent_H;}
	double get_MP()		{return MP;}
	double get_CTE()	{return CTE;}

	int get_temperature_depend(){return temperature_depend;}
	double get_density2()	{return density2;}
	double get_nensei2()	{return nensei2;}
	double get_sigma2()      {return sigma2;}
	double get_vis2()	{return nensei2/density2;}//���S���W��
	double get_Cp2() 	{return Cp2;}
	double get_k2()		{return k2;}
	double get_latent_H2()	{return latent_H2;}
	double get_MP2()		{return MP2;}
	double get_CTE2()	{return CTE2;}
	
	int get_fluidwidth() {return fluidwidth;}
	int get_fluid_h() {return fluid_h;}
	double get_distancebp() {return distancebp;}
	double get_wlength()    {return wlength;}
	double get_height()	{return height;}
	double get_tool_angle() {return tool_angle;}
	int    get_tool_type() {return tool_type;}
	int	get_process_type()	{return process_type;}
	double get_dwelling_time(){return dwelling_time;}
	int    get_airwall() {return airwall;}
	int get_change_step(){return change_step;}
	
	double get_re()		{return re;}
	double get_re2()	{return re2;}
	double get_re3()	{return re3;}
	double get_re4()	{return re4;}	
	double get_beta()	{return beta;}
	int    get_dx()		{return dx;}
	double get_times()	{return times;}
	int get_X_mesh()	{return (int)((maxX-minX)/(distancebp*dx)+0.00000000000001);}//X�������̊i�q�� �ۂߌ덷��h�����߂�0.001�𑫂��Ă���
	int get_Y_mesh()	{return (int)((maxY-minY)/(distancebp*dx)+0.00000000000001);}//Y�������̊i�q�� �ۂߌ덷��h�����߂�0.001�𑫂��Ă���
	int get_Z_mesh()	{return (int)((maxZ-minZ)/(distancebp*dx)+0.00000000000001);}//Z�������̊i�q�� �ۂߌ덷��h�����߂�0.001�𑫂��Ă���
	int get_number_of_mesh(){return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);}//�i�q���FX_mesh*Y_mesh*Z_mesh
	
	int    get_surface_tension()	{return surface_tension;}
	int    get_smooth()	{return smooth;}
	int    get_SFT()	{return SFT;}
	int    get_smn()	{return smn;}
	int    get_smth_sumP()  {return smth_sumP;}
	int    get_wall_direct(){return wall_direct;}
	int    get_non_slip()	{return non_slip;}
	int    get_suf_modify(){return suf_modify;}
	int    get_interpolate_curv() {return interpolate_curv;}
	double get_Cst() {return Cst;}				//�|�e���V�����W��
	void   set_Cst(double calc_Cst) {Cst=calc_Cst;}	//�|�e���V�����W��
	int    get_C_type() {return C_type;}			//�|�e���V�����W���̌�����@
	double get_C_times() {return C_times;}			//�|�e���V�����W���̔{��
	double get_wall_C() {return wall_C;}			//�ǂƂ̐e�a�͌W��

	int    get_iteration_count(){return iteration_count;}
	int    get_Pgrad()	{return Pgrad;}
	int    get_ave_Pdirect(){return ave_Pdirect;}
	int    get_solution()	{return solution;}
	int    get_B_of_P()	{return B_of_P;}
	double get_w_div()	{return w_div;}
	int    get_divU_method() {return divU_method;}//���x���U�̌v�Z���@ 1=MPS,2=WLSM,3=���g��ʂ�Ȃ�WLSM 4=������
	int    get_divU_order()	{return divU_order;}//WSLM��divU���v�Z����Ƃ��̃I�[�_�[ MAX=3
	int    get_HL_sw()	{return HL_sw;}
	int    get_dir_for_P(){return dir_for_P;}
	int    get_initialP()	{return initialP;}
	int    get_P_twice(){return P_twice;}
	double get_CGep()	{return CGep;}
	int    get_omp_P() {return omp_P;}
	int    get_negativeP(){return negativeP;}
	int    get_minP()	{return minP;}
	int    get_set_P()	{return set_P;}
	int    get_artP_sw(){return artP_sw;}
	double get_artP()	{return artP;}
	double get_Pgrad_times(){return Pgrad_times;}
	int get_Pgrad_order(){return Pgrad_order;}
	int    get_P_AVS()	{return P_AVS;}
	double get_Pgrad_limit(){return Pgrad_limit;}
	int    get_gridless_P()	{return gridless_P;}
	int    get_wall_poly()	{return wall_poly;}
	int    get_interpolate_surface_u(){return interpolate_surface_u;}
	
	int    get_T_field()	{return T_field;}
	int    get_insulate()	{return insulate;}
	int    get_T_laplacian(){return T_laplacian;}
	double get_wall_density(){return wall_density;}
	double get_wall_Cp()    {return wall_Cp;}
	double get_wall_k()		{return wall_k;}
	double get_initialT()	{return initialT;}
	double get_roomT()	{return roomT;}
	int    get_air_cool() {return air_cool;}
	int    get_T_expansion(){return T_expansion;}
	int    get_T_solution(){return T_solution;}
    int    get_T_CG(){return T_CG;}
	double get_T_CGep(){return T_CGep;}
	int    get_buoyant()	{return buoyant;}
	double get_TplotZ()		{return TplotZ;}
	int    get_T_AVS()		{return T_AVS;}
	int get_output_temperature_face()	{return		output_temperature_face;}

	int    get_EM_method()	{return EM_method;}
	int    get_EM_calc_type() {return EM_calc_type;}
	int    get_EM_interval(){return EM_interval;}
	double get_EM_distance(){return EM_distance;}
	int    get_region_shape() {return region_shape;}
	double get_XL()		{return XL;}
	double get_XR()		{return XR;}
	double get_YU()		{return YU;}
	double get_YD()		{return YD;}
	double get_ZU()		{return ZU;}
	double get_ZD()		{return ZD;}
	double get_RU()		{return RU;}

	int    get_mesher() {return mesher;}
	int    get_mesh_input() {return mesh_input;}
	int    get_remesh_sw()	{return remesh_sw;}
	double get_boxalpha(){return boxalpha;}
	int    get_fine()	{return fine;}
	double get_co_fine(){return co_fine;}
	int    get_add_points()	{return add_points;}
	int    get_poly_memory(){return poly_memory;}
	int	   get_CDT_sw(){return CDT_sw;}
	int    get_defer_f(){return defer_f;}
	int    get_mesh_smth(){return mesh_smth;}
	double get_material_judge(){return material_judge;}
	double get_del_length(){return del_length;}
	int	   get_BD_fluid(){return BD_fluid;}
	int    get_output_wall_way()	{return output_wall_way;}
	double get_layer_thin()		{return layer_thin;}	//��C�w�̌��� distancebp*layer_thin�����ۂ̌���
	int    get_layer_num()				{return layer_num;}//��C�w�̑w��

	int	   get_BEM_elm_type(){return BEM_elm_type;}
	int    get_FEM_elm_type() {return FEM_elm_type;}
	int    get_FEM_smn()	{return FEM_smn;}
	int    get_max_DN()	{return max_DN;}
	int    get_FEMCG()	{return FEMCG;}
	double    get_CGaccl()	{return CGaccl;}
	double get_FEMCGep() {return FEMCGep;}
	double get_FEMtimes(){return FEMtimes;}
	double get_legend_F(){return legend_F;}
	int	   get_node_sort(){return node_sort;}
	int	   get_thinout_fluid(){return thinout_fluid;}
	
	double get_V()		{return V;}
	int    get_V_step()	{return V_step;}
	double get_r_perm()	{return r_perm;}
	int    get_V_con()	{return V_con;}
	double get_initial_V()	{return initial_V;}
	double get_E_times(){return E_times;}
	int	   get_eleforce(){return eleforce;}
	int    get_charge()	{return charge;}
	
	int    get_J_input_way() {return J_input_way;}
	double get_J0()		{return J0;}
	double get_I0()		{return I0;}
	double get_RP()		{return RP;}
	double get_ele_conduc() {return ele_conduc;}
	double get_ele_conduc2() {return ele_conduc2;}
	double    get_Hz()		{return Hz;}
	int    get_div_Hz()	{return div_Hz;}
	int    get_jheat()	{return jheat;}
	int    get_m_force(){return m_force;}
	int	   get_NLBHI()	{return NLBHI;}
	int    get_NLMH()	{return NLMH;}
	double get_magnet_H(){return magnet_H;}
	double get_magnet_B(){return magnet_B;}
	int    get_uniform_B_sw(){return uniform_B_sw;}
	double get_uniform_B(){return uniform_B;}
	double get_B_times(){return B_times;}
	int    get_plot_B_type(){return plot_B_type;}
	int	   get_Je_crucible(){return Je_crucible;}
	int    get_J0eqJe(){return J0eqJe;}
	int    get_m_A(){return m_A;}
	int    get_jw_interval(){return jw_interval;}
	int    get_jw_Faverage(){return jw_Faverage;}
	int    get_static_dirichlet(){return static_dirichlet;}
	double get_eddy_times()  {return eddy_times;}
	int		get_A_phi() {return A_phi;}
	int    get_parabolic_node_element() {return parabolic_node_element;}

	double get_g()		{return g;}
	int    get_restart()    {return restart;}
	int    get_autosave()	{return autosave;}
	double get_curan()    {return curan;}
	int    get_modify_position() {return modify_position;}
	int    get_vis_calc_type(){return vis_calc_type;}
	int    get_wall_adheision() {return wall_adheision;}
	int    get_laplacian()	{return laplacian;}
	int    get_vis_solver() {return vis_solver;}
	int    get_initial_u()  {return initial_u;}
	int    get_temporary_r(){return temporary_r;}
	
	int    get_fix_center() {return fix_center;}
	int    get_freeon()	{return freeon;}
	int    get_freeon3sw(){return freeon3sw;}
	int    get_surface_judge2(){return surface_judge2;}
	int    get_move_prtcl(){return move_prtcl;}
	int    get_move_u_dirct(){return move_u_dirct;}
	double get_move_speed() {return move_speed;}
	double get_move_speed2(){return move_speed2;}
	int    get_check_something() {return check_something;}
	int    get_check_region() {return check_region;}
	double get_max_speed() {return max_speed;}
	int    get_set_zero_speed(){return set_zero_speed;}
	int    get_AVS_HALF(){return AVS_HALF;}
	int    get_adaptive_sw(){return adaptive_sw;}	//�𑜓x���ςɂ��邩�A���Ȃ��� 1=ON 0=OFF
	double get_threshold(){return threshold;}		//�A�_�v�e�B�u�ɂ���ۂ́A���͌덷��臒l
	int get_fix_surface(){return fix_surface;}
	int get_output_viscousity_face(){return output_viscousity_face;}
	int get_output_equivalent_strain_rate_face(){return output_equivalent_strain_rate_face;}
	int get_output_forward(){return output_forward;}
	int get_output_backward(){return output_backward;}
	int get_output_another_face(){return output_another_face;}

	int    get_model_number(){return model_number;}
	int    get_model_set_way(){return model_set_way;}
	int    get_model_inherit(){return model_inherit;}
	int   get_ea_model(){return ea_model;}

	int    get_speed_plot_particle(){return speed_plot_particle;}
	double get_speedtimes(){return speedtimes;}
	int    get_speed_face() {return speed_face;}
	double get_speed_face_p(){return speed_face_p;}
	int    get_ax_sym_modify(){return ax_sym_modify;}
	int    get_flat_speed_plot(){return flat_speed_plot;}
	double get_flat_speed_p() {return flat_speed_p;}
	int    get_relative_speed(){return relative_speed;}
	int    get_speed_AVS(){return speed_AVS;}
	double get_legend_speed(){return legend_speed;}

	int    get_M_form()		{return M_form;}
	int    get_MAX_thread()	{return MAX_thread;}
	
	
	int    get_interval()	{return interval;}
	int    get_AVS()	{return AVS;}
	double get_maxT()	{return maxT;}
	double get_minT()	{return minT;}
	int		get_P_interval()	{return P_interval;}
	int		get_ST_interval()	{return ST_interval;}
	int		get_PND_interval()	{return PND_interval;}
	int    get_mesh_output_interval()	{return mesh_output_interval;}
	int    get_plot_nearbywall() {return plot_nearbywall;}

	double get_particle_mass()///////����				//���q�̑̐ρB�����z�u�̎d���ɉ����Čv�Z���@��ς���
	{
		double mass;
		if(model_set_way==0)	//�����i�q�̂Ƃ�
		{
			if(dimention==2){mass=density*distancebp*distancebp;}
			else			{mass=density*distancebp*distancebp*distancebp;}
		}
		else if(model_set_way==1)	//�ז��i�q�̂Ƃ�
		{
			if(dimention==2){mass=density*sqrt(3.0)/4*distancebp*distancebp;}
			else
			{
				//mass=density*sqrt(2.0)/12*distancebp*distancebp*distancebp;
				mass=density*distancebp*distancebp*distancebp/sqrt(2.0);
			}
		}
		else cout<<"���f���̐ςݕ����s��ł� ���ʂ��v�Z�ł��܂���"<<endl;
		return mass;
	}

	double get_particle_volume(){ return (get_particle_mass()/get_density());}
	
};


#endif
