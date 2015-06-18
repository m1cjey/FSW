#include"header.h"	//主要なヘッダーファイルはまとめてこのなか。

///initial for MPS
double kernel(double r,double dis);
double kernel2(double r,double dis,double d);
double kernel_in_WLSM(double dis, double R);
double initial_pnd(double r,int dimention,int calc_type);
double calclamda(mpsconfig *CON);
void input_particle_data(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int t);
void culan(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double *dt,double mindis,double Umax,double *g);
double get_volume(mpsconfig *CON);
void output_particle_density(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double n0,int particle_number,int t);
double calc_Cst(mpsconfig *CON);

//粒子数カウント関数 & 並び替え
void calc_numbers_of_particles_and_change_the_order(mpsconfig *CON,vector <mpsparticle> &PART,int particle_number,int *fluid_number,int *out,int *order_sw);
int check_position(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int *particle_number);

//粒子法におけるモデル作成
void writedata(ofstream fp, int id, double x, double y,double z, int type,int surface,double val,ofstream gnu,double vx,double vy,double vz,double P,double h,int toBEM);
void set_initial_placement(mpsconfig *CON,int *particle_number);
void set_initial_placement_using_MD(mpsconfig *CON,int *particle_number);

//格子(index)生成関係
void reload_INDEX(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int *INDEX);
void reload_INDEX2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int **MESH);

//表面判定
void calc_neighbor_relation(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int fluid_number,int out);
void freeon(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,int fluid_number,int out);
void freeon2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,int fluid_number,int out);
void surface_judge2(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);
void surface_judge2_old(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);
void surface_judge2_new(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);
void surface_judge3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number);
void freeon3(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int out);

//法線ベクトル
double pnd_for_direct(mpsconfig *CON,vector<mpsparticle> &PART,double x,double y,double z,double R,int i);
double pnd_for_direct2(mpsconfig *CON,vector<mpsparticle> &PART,double x,double y,double z,double R,int i);
void direct_f(mpsconfig *CON,vector<mpsparticle> &PART,int i,double *direct[DIMENTION]);
void direct_f2(mpsconfig *CON,vector<mpsparticle> &PART,int i,double *direct[DIMENTION]);

//ファイル出力
void plot_speed(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int fluid_number);
void post_processing(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double Umax,double mindis,int t,double TIME,unsigned int time0,int count_avs);
void particle_movie_AVS(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double T);
void post_processing3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double TIME);
void output_alldata_AVS(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double Umax,double mindis,int t,double TIME, unsigned int time0,int count_avs);

//行列解法
void CG_method(mpsconfig *CON,double *r,double *P,double *AP,double *val,int *ind,int *ptr,int pn,double *X,int *countN,double EP);
void iccg(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X,double *r,double *P,double EP,int *count2);
void gauss(double *matrix,double *B,int N);
void jacobi(double **matrix,double *B,int N);

//粘性
void visterm_negative(mpsconfig *CON,vector<mpsparticle> &PART,double *laplacian[DIMENTION],double n0,double lamda,int fluid_number,int particle_number,double dt,int t);
void u_laplacian_f(mpsconfig *CON,vector<mpsparticle> &PART,double *laplacian[DIMENTION],double n0,double lamda,int fluid_number,double dt);
void calc_viscous_term(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,double N0,double *laplacian[DIMENTION],double lamda,int t);
void calc_vis_value(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *vis,double dt,int t,int particle_number);
void output_viscousity_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number);
void output_flow_stress_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number);
void output_equivalent_strain_rate_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number);
void physical_quantity_movie_AVS(mpsconfig *CON,int t,vector<mpsparticle>&PART,int particle_number);

//表面張力
void smoothing(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double *potential[DIMENTION],double n0);
void surface_tension1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int particle_number);
void surface_tension2(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int particle_number);
void calc_surface_tension(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double dt,int particle_number,double n0,double **potential,int t);
//
void plot_ST(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int t);
void plot_ST_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *potential[DIMENTION],int t);
void plot_PND(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t);
void plot_PND_each(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *PND4,int t);
void plot_F(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *F[DIMENTION],int t);
void plot_speed_each(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int fluid_number,int t);
void plot_speed_tool(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int fluid_number,int t);
//

/////仮の速度および位置決定
void renewal_u_and_r_in_positive(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int t,double dt,double *Umax,double **potential,double **laplacian,double *g,double **previous_Un,double **F);

//圧力計算
void negative1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lamda,double N0,double *PND2,double n0,double **Un);
void negative1_twice(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,int t,double dt,double lamda,double N0,double *PND2,double n0,double **Un,double n0_4);
void pressure(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int out,double dt,int t,double lamda,double N0,double *PND2,double n0,int B_of_P,int negativeP);
void calc_P_main(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double lamda,double N0,double dt,double *PND2,double n0,int fluid_number,int out,int pn,int B_of_P);
void set_N0_and_PND2(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double *N0,int out);
void set_Dirichlet_P(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int flag,double *Dirichlet_P);
void output_dirichlet_vector_files(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int flag,double *Dirichlet_P);
double Dndt(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0);
void calc_P_main_with_gridless(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double lamda,double N0,double dt,double *PND2,double n0,int fluid_number,int out,int pn,int B_of_P,int t);
//
void plot_P(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int t,int fluid_number);
void plot_P_each(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,int t,int fluid_number);
//
void output_pressuer_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number);
void calc_Pgradient(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double **reU,double n0,double dt,int minPsw);
void set_minP(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *minP,int minPsw);
void plot_Pgradient(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **Pgrad);
void P_gradient3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *direct[DIMENTION],double dt,double **reU,double **Pgrad,int minPsw);
void cacl_WLSM_P_D3_order1(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP);
void cacl_WLSM_P_D3_order2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B,double **P_grad,int i,int N,double *minP);
void P_gradient_MPS(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0,double *P_grad);
void P_gradient4(mpsconfig *CON,vector<mpsparticle> &PART,double dt,int fluid_number,double *reU[3], double *Pgrad[3],int minPsw);
void P_gradient5(mpsconfig *CON,vector<mpsparticle> &PART,double dt,int fluid_number,double *reU[3],double *P_grad[3]);
void negative3(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double lamda,double N0,vector<point3D> &NODE,vector<element3D> &ELEM,int out);
void modify_reU(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,int t,double dt,double n0,double *udiv, double **reU);
void calc_inverse_matrix(mpsconfig *CON,vector<mpsparticle> &PART,int N, double *matrix);
void BiCGStab2_method_with_D_scale_for_sparse(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double *val,int *ind,int *ptr);
void MRTR(mpsconfig *CON,double *r,int pn,double *X,int *countN,double EP,double *val,int *ind,int *ptr);
double calc_PND_by_minmum_L(mpsconfig *CON,int dimention,double Re,double L);
void reset_surface_velocity_by_WLSM(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number);

//陰解析後の速度・座標更新
void modify_u_and_x_after_Pcalc(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **reU,double dt,double **Un);


//速度発散
double divergence(mpsconfig *CON,vector<mpsparticle> &PART,int i,double n0);
double divergence2(mpsconfig *CON,vector<mpsparticle> &PART,int i,int surface_sw);
double divergence3(mpsconfig *CON,vector<mpsparticle> &PART,int i,int surface_sw);
double divergence4(mpsconfig *CON,vector<mpsparticle> &PART,int i);
void return_X_for5N(double *matrix,int N,double *B1,double *B2,double *dudx,double *dudy);
double cacl_WLSM_divu_D3_order1(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);
double cacl_WLSM_divu_D3_order1_2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);
double cacl_WLSM_divu_D3_order2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);
double cacl_WLSM_divu_D3_order2_2(mpsconfig *CON,vector<mpsparticle> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N);

//粒子法その他
void modify_position(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double dt,int particle_number);
void delete_particle(mpsconfig &CON,vector<mpsparticle> &PART,int *particle_number,int *fluid_number,double n0_4,int t);
void calc_physical_property(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *val,int particle_number,int sw);
void calc_wallZ(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,double *mindis,int fluid_number,int out);

//BEM
void BEM2D(mpsconfig *CON,vector<point2D> &s_NODE,vector<element2D> &s_ELEM,int *s_node_num,int *s_elem_num,vector<mpsparticle> &PART,int fluid_number,int particle_number,double **F,int t);
void gauss(double **matrix,double *B,int N);
void BEM3D(mpsconfig *CON,vector<BEMpoint3D> &s_NODE,vector<BEMelement3D> &s_ELEM,int *s_node_num,int *s_elem_num,vector<mpsparticle> &PART,int fluid_number,int particle_number,double **F,int t);
void BEM3D_main_for_CONSTANT(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,vector<mpsparticle> &PART,int fluid_number,double **F,vector<REGION> &region);

//BEMのinputファイル生成
void set_BEM_static_model(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,int *s_node_num,int *s_elem_num);
void set_BEM2D_dynaic_model_for_CONSTANT(mpsconfig *CON,vector<mpsparticle> &PART,vector<point2D> &NODE,vector<element2D> &ELEM,int *d_node_num,int *d_elem_num,int particle_number,int fluid_number);
void couple_NODE_and_ELEM(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,vector<point2D> &s_NODE,vector<element2D> &s_ELEM,vector<point2D> &dy_NODE,vector<element2D> &dy_ELEM,vector<REGION> &region);
void set_gauss_parameter(int *NGS,double *ep1,double *w1,double *ep_ln,double *w_ln);
void BEM_main_for_CONSTANT(mpsconfig *CON,vector<point2D> &NODE,vector<element2D> &ELEM,vector<mpsparticle> &PART,int fluid_number,double **F,vector<REGION> &region);
void set_GHmatrix_for_constant(vector<element2D> &ELEM,double ep,double **H,double **G,double w,int type,int n,int m,int n1,int m1);
void smoothingF3D(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double *F[3],int t);
void set_BEM3D_dynaic_model_for_CONSTANT(mpsconfig *CON,vector<mpsparticle> &PART,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,int *d_node_num,int *d_elem_num,int particle_number,int fluid_number);
void couple3D_NODE_and_ELEM(mpsconfig *CON,vector<BEMpoint3D> &NODE,vector<BEMelement3D> &ELEM,vector<BEMpoint3D> &s_NODE,vector<BEMelement3D> &s_ELEM,vector<BEMpoint3D> &dy_NODE,vector<BEMelement3D> &dy_ELEM,vector<REGION> &region);



//磁気モーメント法
void Magnetic_Moment_Method(mpsconfig *CON,vector<mpsparticle> &PART,double **F,double n0,double lamda,int fluid_number,int particle_number);
void H_gradient1(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double **Hgrad,double **H);
void output_matrix(double **matrix,double *Bmatrix,int node_num);
//温度場
void plot_T(mpsconfig *CON ,vector<mpsparticle> &PART,int particle_number,double *T,double height);
void calc_Temperature(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t);
void calc_temperature_implicity(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,int particle_number,double n0,double lamda,double dt,int t);
void output_temperature_avs(mpsconfig *CON,vector<mpsparticle> &PART,int t,int particle_number,int fluid_number,double *T,double height);
void set_temperature_boundary(mpsconfig *CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, double *T,double dt,int t);
void set_Q_boundary(mpsconfig *CON, vector<mpsparticle> &PART, int fluid_number, int particle_number, double dt,int t);

void move_particle(mpsconfig *CON,vector<mpsparticle> &PART,int particle_number,int fluid_number,double dt);
void check_something(mpsconfig *CON,vector<mpsparticle> &PART,int fluid_number,double n0,int particle_number,int t);

//デローニ分割
void delaun3D(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,double rax,double ray,double raz,int *node_num,int *nelm,int FINE_sw,double rrm);
void data_avs(int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,double *val,mpsconfig *CON);
void data_avs2(mpsconfig *CON,int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,int ktj,int t);
void data_avs2flux(mpsconfig *CON,int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,double *val,int t,int flux);
void data_avs2node(mpsconfig *CON,int node_number,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,double *val,int t);
void data_avs3(int node,int nelm,vector <point3D> &NODE,vector <element3D> &ELEM,mpsconfig *CON,int t);
void delaun3D_main(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node_num,int *nelm,int FINE_sw);
void memorize_static_NODE_and_ELEM(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,vector<point3D> &static_NODE,vector<element3D> &static_ELEM,int node,int nelm);
void set_jnb3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb);
void set_nei3D(vector<point3D> &NODE, vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei);
double volume3D(vector<point3D> &NODE,int ia,int ib,int ic,int ip);
void sphere3D(vector<point3D> &NODE,vector<element3D> &ELEM,int ia,int ib,int ic,int ip,int i);
void poly3D(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void set_material(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm);
void FINE3D(vector<point3D> &NODE,vector<element3D> &ELEM,int KTJ,int KTE,int *node,int *nelm,mpsconfig *CON,double rrm,int startID);
void fill3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm);


//MPSTOFEM3D
void MPS_TO_FEM3Dmain(mpsconfig *CON,int *node_num,vector<point3D> &NODE,vector<mpsparticle> &PART, int fluid_number, int particle_number);

//FEM3D
int FEM3D_calculation(mpsconfig *CON,int *static_node,int *static_nelm,int *static_nedge,vector<point3D> &static_NODE, vector<element3D> &static_ELEM,vector<edge3D> &static_EDGE,int particle_node,double **F,int t,double TIME,vector<mpsparticle> &PART,int fluid_number,int particle_number,double dt,vector<point3D> &NODE_jw, vector<element3D> &ELEM_jw,mpsconfig& CONF);
int calc_matrix_width(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei);
int calc_matrix_width2(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,int *width);
void arrange_matrix(int pn,int *NUM,int **ROW,double **G);
void Incomplete_Cholesky_Decomposition(mpsconfig *CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2);
void ICCG3D2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void parallel_ICCG3D2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void check_matrix_symmetry(int pn,int *NUM,int **ROW,double **G);
void output_F_scalar_with_AVS_for_linear(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int t,vector<mpsparticle> &PART,int node);
void output_F_scalar_movie_with_AVS_for_linear(mpsconfig *CON,vector<point3D> &NODE,vector<element3D> &ELEM,int t,vector<mpsparticle> &PART,double TIME);
