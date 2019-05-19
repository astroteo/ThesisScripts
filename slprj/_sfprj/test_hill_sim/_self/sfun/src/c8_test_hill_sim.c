/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_hill_sim_sfun.h"
#include "c8_test_hill_sim.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "test_hill_sim_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static real_T _sfTime_;
static const char * c8_debug_family_names[23] = { "Power_Sun",
  "Power_Earth_reflected", "Power_Earth", "c", "P_S", "P_E", "P_EE", "S_s",
  "S_e", "a_p", "b_p", "Ap", "Cd_p", "Cs_p", "F_E", "F_S", "nargin", "nargout",
  "R_Sun_c", "R_c", "e_f", "m_cubesat", "a_Sun" };

/* Function Declarations */
static void initialize_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void initialize_params_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct *
  chartInstance);
static void enable_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void disable_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void c8_update_debugger_state_c8_test_hill_sim
  (SFc8_test_hill_simInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c8_test_hill_sim
  (SFc8_test_hill_simInstanceStruct *chartInstance);
static void set_sim_state_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_st);
static void finalize_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void sf_gateway_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void mdl_start_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void c8_chartstep_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void initSimStructsc8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c8_machineNumber, uint32_T
  c8_chartNumber, uint32_T c8_instanceNumber);
static const mxArray *c8_sf_marshallOut(void *chartInstanceVoid, void *c8_inData);
static void c8_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct *chartInstance,
  const mxArray *c8_b_a_Sun, const char_T *c8_identifier, real_T c8_y[3]);
static void c8_b_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId,
  real_T c8_y[3]);
static void c8_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static const mxArray *c8_b_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static const mxArray *c8_c_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static real_T c8_c_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId);
static void c8_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static void c8_info_helper(const mxArray **c8_info);
static const mxArray *c8_emlrt_marshallOut(const char * c8_u);
static const mxArray *c8_b_emlrt_marshallOut(const uint32_T c8_u);
static real_T c8_eml_xnrm2(SFc8_test_hill_simInstanceStruct *chartInstance,
  real_T c8_x[3]);
static void c8_below_threshold(SFc8_test_hill_simInstanceStruct *chartInstance);
static const mxArray *c8_d_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData);
static int32_T c8_d_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId);
static void c8_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData);
static uint8_T c8_e_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_b_is_active_c8_test_hill_sim, const char_T
  *c8_identifier);
static uint8_T c8_f_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId);
static void init_dsm_address_info(SFc8_test_hill_simInstanceStruct
  *chartInstance);
static void init_simulink_io_address(SFc8_test_hill_simInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c8_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c8_is_active_c8_test_hill_sim = 0U;
}

static void initialize_params_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct *
  chartInstance)
{
  (void)chartInstance;
}

static void enable_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c8_update_debugger_state_c8_test_hill_sim
  (SFc8_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c8_test_hill_sim
  (SFc8_test_hill_simInstanceStruct *chartInstance)
{
  const mxArray *c8_st;
  const mxArray *c8_y = NULL;
  int32_T c8_i0;
  real_T c8_u[3];
  const mxArray *c8_b_y = NULL;
  uint8_T c8_hoistedGlobal;
  uint8_T c8_b_u;
  const mxArray *c8_c_y = NULL;
  c8_st = NULL;
  c8_st = NULL;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_createcellmatrix(2, 1), false);
  for (c8_i0 = 0; c8_i0 < 3; c8_i0++) {
    c8_u[c8_i0] = (*chartInstance->c8_a_Sun)[c8_i0];
  }

  c8_b_y = NULL;
  sf_mex_assign(&c8_b_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c8_y, 0, c8_b_y);
  c8_hoistedGlobal = chartInstance->c8_is_active_c8_test_hill_sim;
  c8_b_u = c8_hoistedGlobal;
  c8_c_y = NULL;
  sf_mex_assign(&c8_c_y, sf_mex_create("y", &c8_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c8_y, 1, c8_c_y);
  sf_mex_assign(&c8_st, c8_y, false);
  return c8_st;
}

static void set_sim_state_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_st)
{
  const mxArray *c8_u;
  real_T c8_dv0[3];
  int32_T c8_i1;
  chartInstance->c8_doneDoubleBufferReInit = true;
  c8_u = sf_mex_dup(c8_st);
  c8_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c8_u, 0)),
                      "a_Sun", c8_dv0);
  for (c8_i1 = 0; c8_i1 < 3; c8_i1++) {
    (*chartInstance->c8_a_Sun)[c8_i1] = c8_dv0[c8_i1];
  }

  chartInstance->c8_is_active_c8_test_hill_sim = c8_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c8_u, 1)),
     "is_active_c8_test_hill_sim");
  sf_mex_destroy(&c8_u);
  c8_update_debugger_state_c8_test_hill_sim(chartInstance);
  sf_mex_destroy(&c8_st);
}

static void finalize_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  int32_T c8_i2;
  int32_T c8_i3;
  int32_T c8_i4;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 7U, chartInstance->c8_sfEvent);
  for (c8_i2 = 0; c8_i2 < 3; c8_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c8_R_Sun_c)[c8_i2], 0U);
  }

  for (c8_i3 = 0; c8_i3 < 3; c8_i3++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c8_R_c)[c8_i3], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c8_e_f, 2U);
  chartInstance->c8_sfEvent = CALL_EVENT;
  c8_chartstep_c8_test_hill_sim(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_test_hill_simMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c8_i4 = 0; c8_i4 < 3; c8_i4++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c8_a_Sun)[c8_i4], 3U);
  }

  _SFD_DATA_RANGE_CHECK(*chartInstance->c8_m_cubesat, 4U);
}

static void mdl_start_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void c8_chartstep_c8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  real_T c8_hoistedGlobal;
  real_T c8_b_hoistedGlobal;
  int32_T c8_i5;
  real_T c8_b_R_Sun_c[3];
  int32_T c8_i6;
  real_T c8_b_R_c[3];
  real_T c8_b_e_f;
  real_T c8_b_m_cubesat;
  uint32_T c8_debug_family_var_map[23];
  real_T c8_Power_Sun;
  real_T c8_Power_Earth_reflected;
  real_T c8_Power_Earth;
  real_T c8_c;
  real_T c8_P_S;
  real_T c8_P_E;
  real_T c8_P_EE;
  real_T c8_S_s[3];
  real_T c8_S_e[3];
  real_T c8_a_p;
  real_T c8_b_p;
  real_T c8_Ap;
  real_T c8_Cd_p;
  real_T c8_Cs_p;
  real_T c8_F_E[3];
  real_T c8_F_S[3];
  real_T c8_nargin = 4.0;
  real_T c8_nargout = 1.0;
  real_T c8_b_a_Sun[3];
  int32_T c8_i7;
  int32_T c8_i8;
  real_T c8_x[3];
  int32_T c8_i9;
  real_T c8_b_x[3];
  real_T c8_y;
  int32_T c8_i10;
  real_T c8_B;
  real_T c8_b_y;
  real_T c8_c_y;
  real_T c8_d_y;
  int32_T c8_i11;
  int32_T c8_i12;
  int32_T c8_i13;
  real_T c8_c_x[3];
  real_T c8_e_y;
  int32_T c8_i14;
  real_T c8_b_B;
  real_T c8_f_y;
  real_T c8_g_y;
  real_T c8_h_y;
  int32_T c8_i15;
  real_T c8_a;
  int32_T c8_i16;
  int32_T c8_i17;
  real_T c8_b_a;
  int32_T c8_i18;
  int32_T c8_i19;
  int32_T c8_i20;
  real_T c8_c_B;
  real_T c8_i_y;
  real_T c8_j_y;
  real_T c8_k_y;
  int32_T c8_i21;
  int32_T c8_i22;
  real_T c8_A[3];
  real_T c8_d_B;
  real_T c8_l_y;
  real_T c8_m_y;
  real_T c8_n_y;
  int32_T c8_i23;
  int32_T c8_i24;
  int32_T c8_i25;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 7U, chartInstance->c8_sfEvent);
  c8_hoistedGlobal = *chartInstance->c8_e_f;
  c8_b_hoistedGlobal = *chartInstance->c8_m_cubesat;
  for (c8_i5 = 0; c8_i5 < 3; c8_i5++) {
    c8_b_R_Sun_c[c8_i5] = (*chartInstance->c8_R_Sun_c)[c8_i5];
  }

  for (c8_i6 = 0; c8_i6 < 3; c8_i6++) {
    c8_b_R_c[c8_i6] = (*chartInstance->c8_R_c)[c8_i6];
  }

  c8_b_e_f = c8_hoistedGlobal;
  c8_b_m_cubesat = c8_b_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 23U, 23U, c8_debug_family_names,
    c8_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_Power_Sun, 0U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_Power_Earth_reflected, 1U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_Power_Earth, 2U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_c, 3U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_P_S, 4U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_P_E, 5U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_P_EE, 6U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_S_s, 7U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_S_e, 8U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_a_p, 9U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_b_p, 10U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_Ap, 11U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_Cd_p, 12U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_Cs_p, 13U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_F_E, 14U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_F_S, 15U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_nargin, 16U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c8_nargout, 17U, c8_b_sf_marshallOut,
    c8_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c8_b_R_Sun_c, 18U, c8_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c8_b_R_c, 19U, c8_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_b_e_f, 20U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c8_b_m_cubesat, 21U, c8_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c8_b_a_Sun, 22U, c8_sf_marshallOut,
    c8_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 4);
  for (c8_i7 = 0; c8_i7 < 3; c8_i7++) {
    c8_b_a_Sun[c8_i7] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 6);
  c8_Power_Sun = 1358.0;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 7);
  c8_Power_Earth_reflected = 600.0;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 8);
  c8_Power_Earth = 150.0;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 9);
  c8_c = 2.99792458E+8;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 10);
  c8_P_S = 4.5298004127909048E-6;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 11);
  c8_P_E = 2.0013845711889121E-6;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 12);
  c8_P_EE = 5.00346142797228E-7;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 15);
  for (c8_i8 = 0; c8_i8 < 3; c8_i8++) {
    c8_x[c8_i8] = c8_b_R_Sun_c[c8_i8];
  }

  for (c8_i9 = 0; c8_i9 < 3; c8_i9++) {
    c8_b_x[c8_i9] = c8_x[c8_i9];
  }

  c8_y = c8_eml_xnrm2(chartInstance, c8_b_x);
  for (c8_i10 = 0; c8_i10 < 3; c8_i10++) {
    c8_x[c8_i10] = c8_b_R_Sun_c[c8_i10];
  }

  c8_B = c8_y;
  c8_b_y = c8_B;
  c8_c_y = c8_b_y;
  c8_d_y = c8_c_y;
  for (c8_i11 = 0; c8_i11 < 3; c8_i11++) {
    c8_S_s[c8_i11] = c8_x[c8_i11] / c8_d_y;
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 16);
  for (c8_i12 = 0; c8_i12 < 3; c8_i12++) {
    c8_x[c8_i12] = c8_b_R_c[c8_i12];
  }

  for (c8_i13 = 0; c8_i13 < 3; c8_i13++) {
    c8_c_x[c8_i13] = c8_x[c8_i13];
  }

  c8_e_y = c8_eml_xnrm2(chartInstance, c8_c_x);
  for (c8_i14 = 0; c8_i14 < 3; c8_i14++) {
    c8_x[c8_i14] = c8_b_R_c[c8_i14];
  }

  c8_b_B = c8_e_y;
  c8_f_y = c8_b_B;
  c8_g_y = c8_f_y;
  c8_h_y = c8_g_y;
  for (c8_i15 = 0; c8_i15 < 3; c8_i15++) {
    c8_S_e[c8_i15] = c8_x[c8_i15] / c8_h_y;
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 19);
  c8_a_p = 0.3;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 20);
  c8_b_p = 0.3;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 22);
  c8_Ap = 0.09;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 25);
  c8_Cd_p = 0.1;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 26);
  c8_Cs_p = 0.1;
  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 28);
  c8_a = c8_b_e_f * 2.0013845711889121E-6 * 0.09 + 4.5031152851750521E-8;
  for (c8_i16 = 0; c8_i16 < 3; c8_i16++) {
    c8_x[c8_i16] = c8_S_e[c8_i16];
  }

  for (c8_i17 = 0; c8_i17 < 3; c8_i17++) {
    c8_F_E[c8_i17] = c8_a * c8_x[c8_i17];
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 29);
  c8_b_a = c8_b_e_f * 4.5298004127909048E-6 * 0.09;
  for (c8_i18 = 0; c8_i18 < 3; c8_i18++) {
    c8_x[c8_i18] = c8_S_s[c8_i18];
  }

  for (c8_i19 = 0; c8_i19 < 3; c8_i19++) {
    c8_F_S[c8_i19] = c8_b_a * c8_x[c8_i19];
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, 31);
  for (c8_i20 = 0; c8_i20 < 3; c8_i20++) {
    c8_x[c8_i20] = c8_F_S[c8_i20];
  }

  c8_c_B = c8_b_m_cubesat;
  c8_i_y = c8_c_B;
  c8_j_y = c8_i_y;
  c8_k_y = c8_j_y;
  for (c8_i21 = 0; c8_i21 < 3; c8_i21++) {
    c8_x[c8_i21] /= c8_k_y;
  }

  for (c8_i22 = 0; c8_i22 < 3; c8_i22++) {
    c8_A[c8_i22] = c8_F_E[c8_i22];
  }

  c8_d_B = c8_b_m_cubesat;
  c8_l_y = c8_d_B;
  c8_m_y = c8_l_y;
  c8_n_y = c8_m_y;
  for (c8_i23 = 0; c8_i23 < 3; c8_i23++) {
    c8_A[c8_i23] /= c8_n_y;
  }

  for (c8_i24 = 0; c8_i24 < 3; c8_i24++) {
    c8_b_a_Sun[c8_i24] = c8_x[c8_i24] + c8_A[c8_i24];
  }

  _SFD_EML_CALL(0U, chartInstance->c8_sfEvent, -31);
  _SFD_SYMBOL_SCOPE_POP();
  for (c8_i25 = 0; c8_i25 < 3; c8_i25++) {
    (*chartInstance->c8_a_Sun)[c8_i25] = c8_b_a_Sun[c8_i25];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 7U, chartInstance->c8_sfEvent);
}

static void initSimStructsc8_test_hill_sim(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c8_machineNumber, uint32_T
  c8_chartNumber, uint32_T c8_instanceNumber)
{
  (void)c8_machineNumber;
  (void)c8_chartNumber;
  (void)c8_instanceNumber;
}

static const mxArray *c8_sf_marshallOut(void *chartInstanceVoid, void *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_i26;
  real_T c8_b_inData[3];
  int32_T c8_i27;
  real_T c8_u[3];
  const mxArray *c8_y = NULL;
  SFc8_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc8_test_hill_simInstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  for (c8_i26 = 0; c8_i26 < 3; c8_i26++) {
    c8_b_inData[c8_i26] = (*(real_T (*)[3])c8_inData)[c8_i26];
  }

  for (c8_i27 = 0; c8_i27 < 3; c8_i27++) {
    c8_u[c8_i27] = c8_b_inData[c8_i27];
  }

  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static void c8_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct *chartInstance,
  const mxArray *c8_b_a_Sun, const char_T *c8_identifier, real_T c8_y[3])
{
  emlrtMsgIdentifier c8_thisId;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_b_a_Sun), &c8_thisId, c8_y);
  sf_mex_destroy(&c8_b_a_Sun);
}

static void c8_b_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId,
  real_T c8_y[3])
{
  real_T c8_dv1[3];
  int32_T c8_i28;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), c8_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c8_i28 = 0; c8_i28 < 3; c8_i28++) {
    c8_y[c8_i28] = c8_dv1[c8_i28];
  }

  sf_mex_destroy(&c8_u);
}

static void c8_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_b_a_Sun;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  real_T c8_y[3];
  int32_T c8_i29;
  SFc8_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc8_test_hill_simInstanceStruct *)chartInstanceVoid;
  c8_b_a_Sun = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_b_a_Sun), &c8_thisId, c8_y);
  sf_mex_destroy(&c8_b_a_Sun);
  for (c8_i29 = 0; c8_i29 < 3; c8_i29++) {
    (*(real_T (*)[3])c8_outData)[c8_i29] = c8_y[c8_i29];
  }

  sf_mex_destroy(&c8_mxArrayInData);
}

static const mxArray *c8_b_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  real_T c8_u;
  const mxArray *c8_y = NULL;
  SFc8_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc8_test_hill_simInstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  c8_u = *(real_T *)c8_inData;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", &c8_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static const mxArray *c8_c_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_i30;
  real_T c8_b_inData[3];
  int32_T c8_i31;
  real_T c8_u[3];
  const mxArray *c8_y = NULL;
  SFc8_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc8_test_hill_simInstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  for (c8_i30 = 0; c8_i30 < 3; c8_i30++) {
    c8_b_inData[c8_i30] = (*(real_T (*)[3])c8_inData)[c8_i30];
  }

  for (c8_i31 = 0; c8_i31 < 3; c8_i31++) {
    c8_u[c8_i31] = c8_b_inData[c8_i31];
  }

  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 0, 0U, 1U, 0U, 2, 3, 1), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static real_T c8_c_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId)
{
  real_T c8_y;
  real_T c8_d0;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), &c8_d0, 1, 0, 0U, 0, 0U, 0);
  c8_y = c8_d0;
  sf_mex_destroy(&c8_u);
  return c8_y;
}

static void c8_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_nargout;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  real_T c8_y;
  SFc8_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc8_test_hill_simInstanceStruct *)chartInstanceVoid;
  c8_nargout = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_y = c8_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_nargout), &c8_thisId);
  sf_mex_destroy(&c8_nargout);
  *(real_T *)c8_outData = c8_y;
  sf_mex_destroy(&c8_mxArrayInData);
}

const mxArray *sf_c8_test_hill_sim_get_eml_resolved_functions_info(void)
{
  const mxArray *c8_nameCaptureInfo = NULL;
  c8_nameCaptureInfo = NULL;
  sf_mex_assign(&c8_nameCaptureInfo, sf_mex_createstruct("structure", 2, 35, 1),
                false);
  c8_info_helper(&c8_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c8_nameCaptureInfo);
  return c8_nameCaptureInfo;
}

static void c8_info_helper(const mxArray **c8_info)
{
  const mxArray *c8_rhs0 = NULL;
  const mxArray *c8_lhs0 = NULL;
  const mxArray *c8_rhs1 = NULL;
  const mxArray *c8_lhs1 = NULL;
  const mxArray *c8_rhs2 = NULL;
  const mxArray *c8_lhs2 = NULL;
  const mxArray *c8_rhs3 = NULL;
  const mxArray *c8_lhs3 = NULL;
  const mxArray *c8_rhs4 = NULL;
  const mxArray *c8_lhs4 = NULL;
  const mxArray *c8_rhs5 = NULL;
  const mxArray *c8_lhs5 = NULL;
  const mxArray *c8_rhs6 = NULL;
  const mxArray *c8_lhs6 = NULL;
  const mxArray *c8_rhs7 = NULL;
  const mxArray *c8_lhs7 = NULL;
  const mxArray *c8_rhs8 = NULL;
  const mxArray *c8_lhs8 = NULL;
  const mxArray *c8_rhs9 = NULL;
  const mxArray *c8_lhs9 = NULL;
  const mxArray *c8_rhs10 = NULL;
  const mxArray *c8_lhs10 = NULL;
  const mxArray *c8_rhs11 = NULL;
  const mxArray *c8_lhs11 = NULL;
  const mxArray *c8_rhs12 = NULL;
  const mxArray *c8_lhs12 = NULL;
  const mxArray *c8_rhs13 = NULL;
  const mxArray *c8_lhs13 = NULL;
  const mxArray *c8_rhs14 = NULL;
  const mxArray *c8_lhs14 = NULL;
  const mxArray *c8_rhs15 = NULL;
  const mxArray *c8_lhs15 = NULL;
  const mxArray *c8_rhs16 = NULL;
  const mxArray *c8_lhs16 = NULL;
  const mxArray *c8_rhs17 = NULL;
  const mxArray *c8_lhs17 = NULL;
  const mxArray *c8_rhs18 = NULL;
  const mxArray *c8_lhs18 = NULL;
  const mxArray *c8_rhs19 = NULL;
  const mxArray *c8_lhs19 = NULL;
  const mxArray *c8_rhs20 = NULL;
  const mxArray *c8_lhs20 = NULL;
  const mxArray *c8_rhs21 = NULL;
  const mxArray *c8_lhs21 = NULL;
  const mxArray *c8_rhs22 = NULL;
  const mxArray *c8_lhs22 = NULL;
  const mxArray *c8_rhs23 = NULL;
  const mxArray *c8_lhs23 = NULL;
  const mxArray *c8_rhs24 = NULL;
  const mxArray *c8_lhs24 = NULL;
  const mxArray *c8_rhs25 = NULL;
  const mxArray *c8_lhs25 = NULL;
  const mxArray *c8_rhs26 = NULL;
  const mxArray *c8_lhs26 = NULL;
  const mxArray *c8_rhs27 = NULL;
  const mxArray *c8_lhs27 = NULL;
  const mxArray *c8_rhs28 = NULL;
  const mxArray *c8_lhs28 = NULL;
  const mxArray *c8_rhs29 = NULL;
  const mxArray *c8_lhs29 = NULL;
  const mxArray *c8_rhs30 = NULL;
  const mxArray *c8_lhs30 = NULL;
  const mxArray *c8_rhs31 = NULL;
  const mxArray *c8_lhs31 = NULL;
  const mxArray *c8_rhs32 = NULL;
  const mxArray *c8_lhs32 = NULL;
  const mxArray *c8_rhs33 = NULL;
  const mxArray *c8_lhs33 = NULL;
  const mxArray *c8_rhs34 = NULL;
  const mxArray *c8_lhs34 = NULL;
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("mrdivide"), "name", "name", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1410807648U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1370009886U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c8_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs0), "rhs", "rhs", 0);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs0), "lhs", "lhs", 0);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c8_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs1), "rhs", "rhs", 1);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs1), "lhs", "lhs", 1);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("rdivide"), "name", "name", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c8_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs2), "rhs", "rhs", 2);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs2), "lhs", "lhs", 2);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c8_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs3), "rhs", "rhs", 3);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs3), "lhs", "lhs", 3);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1286818796U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c8_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs4), "rhs", "rhs", 4);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs4), "lhs", "lhs", 4);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_div"), "name", "name", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1386423952U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c8_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs5), "rhs", "rhs", 5);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs5), "lhs", "lhs", 5);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c8_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs6), "rhs", "rhs", 6);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs6), "lhs", "lhs", 6);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "context", "context", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("norm"), "name", "name", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363713868U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c8_rhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs7, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs7), "rhs", "rhs", 7);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs7), "lhs", "lhs", 7);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1323170578U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c8_rhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs8, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs8), "rhs", "rhs", 8);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs8), "lhs", "lhs", 8);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c8_rhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs9, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs9), "rhs", "rhs", 9);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs9), "lhs", "lhs", 9);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!genpnorm"),
                  "context", "context", 10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_xnrm2"), "name", "name",
                  10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1375980692U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c8_rhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs10, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.blas.inline"),
                  "name", "name", 11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/inline.p"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c8_rhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs11, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"), "context",
                  "context", 12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.blas.xnrm2"),
                  "name", "name", 12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "resolved", "resolved", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c8_rhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs12, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.blas.use_refblas"), "name", "name", 13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/use_refblas.p"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c8_rhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs13, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p!below_threshold"),
                  "context", "context", 14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.blas.threshold"),
                  "name", "name", 14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "resolved", "resolved", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c8_rhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs14, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/threshold.p"),
                  "context", "context", 15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c8_rhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs15, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+blas/xnrm2.p"),
                  "context", "context", 16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.refblas.xnrm2"),
                  "name", "name", 16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1410807772U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c8_rhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs16, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("realmin"), "name", "name", 17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1307651242U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c8_rhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs17, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/realmin.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_realmin"), "name", "name",
                  18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1307651244U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c8_rhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs18, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_realmin.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_float_model"), "name",
                  "name", 19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1326727996U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c8_rhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs19, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.indexMinus"),
                  "name", "name", 20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexMinus.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c8_rhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs20, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.indexTimes"),
                  "name", "name", 21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexTimes.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c8_rhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs21, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.indexPlus"),
                  "name", "name", 22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/indexPlus.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1372583160U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c8_rhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs22, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1397257422U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c8_rhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs23, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("isfi"), "name", "name", 24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isfi.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1346510358U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c8_rhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs24, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isfi.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("isnumerictype"), "name",
                  "name", 25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/fixedpoint/isnumerictype.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1398875598U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c8_rhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs25, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("intmax"), "name", "name", 26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c8_rhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs26, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c8_rhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs27, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("intmin"), "name", "name", 28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1362261882U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c8_rhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs28, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m"), "context",
                  "context", 29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_switch_helper"), "name",
                  "name", 29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_switch_helper.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1393330858U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c8_rhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs29, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/+refblas/xnrm2.p"),
                  "context", "context", 30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("abs"), "name", "name", 30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1363713852U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c8_rhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs30, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c8_rhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs31, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1286818712U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c8_rhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs32, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "context", "context", 33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("eml_mtimes_helper"), "name",
                  "name", 33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1383877294U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c8_rhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs33, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m!common_checks"),
                  "context", "context", 34);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 34);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c8_info, c8_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c8_info, c8_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c8_rhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c8_lhs34, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c8_info, sf_mex_duplicatearraysafe(&c8_lhs34), "lhs", "lhs",
                  34);
  sf_mex_destroy(&c8_rhs0);
  sf_mex_destroy(&c8_lhs0);
  sf_mex_destroy(&c8_rhs1);
  sf_mex_destroy(&c8_lhs1);
  sf_mex_destroy(&c8_rhs2);
  sf_mex_destroy(&c8_lhs2);
  sf_mex_destroy(&c8_rhs3);
  sf_mex_destroy(&c8_lhs3);
  sf_mex_destroy(&c8_rhs4);
  sf_mex_destroy(&c8_lhs4);
  sf_mex_destroy(&c8_rhs5);
  sf_mex_destroy(&c8_lhs5);
  sf_mex_destroy(&c8_rhs6);
  sf_mex_destroy(&c8_lhs6);
  sf_mex_destroy(&c8_rhs7);
  sf_mex_destroy(&c8_lhs7);
  sf_mex_destroy(&c8_rhs8);
  sf_mex_destroy(&c8_lhs8);
  sf_mex_destroy(&c8_rhs9);
  sf_mex_destroy(&c8_lhs9);
  sf_mex_destroy(&c8_rhs10);
  sf_mex_destroy(&c8_lhs10);
  sf_mex_destroy(&c8_rhs11);
  sf_mex_destroy(&c8_lhs11);
  sf_mex_destroy(&c8_rhs12);
  sf_mex_destroy(&c8_lhs12);
  sf_mex_destroy(&c8_rhs13);
  sf_mex_destroy(&c8_lhs13);
  sf_mex_destroy(&c8_rhs14);
  sf_mex_destroy(&c8_lhs14);
  sf_mex_destroy(&c8_rhs15);
  sf_mex_destroy(&c8_lhs15);
  sf_mex_destroy(&c8_rhs16);
  sf_mex_destroy(&c8_lhs16);
  sf_mex_destroy(&c8_rhs17);
  sf_mex_destroy(&c8_lhs17);
  sf_mex_destroy(&c8_rhs18);
  sf_mex_destroy(&c8_lhs18);
  sf_mex_destroy(&c8_rhs19);
  sf_mex_destroy(&c8_lhs19);
  sf_mex_destroy(&c8_rhs20);
  sf_mex_destroy(&c8_lhs20);
  sf_mex_destroy(&c8_rhs21);
  sf_mex_destroy(&c8_lhs21);
  sf_mex_destroy(&c8_rhs22);
  sf_mex_destroy(&c8_lhs22);
  sf_mex_destroy(&c8_rhs23);
  sf_mex_destroy(&c8_lhs23);
  sf_mex_destroy(&c8_rhs24);
  sf_mex_destroy(&c8_lhs24);
  sf_mex_destroy(&c8_rhs25);
  sf_mex_destroy(&c8_lhs25);
  sf_mex_destroy(&c8_rhs26);
  sf_mex_destroy(&c8_lhs26);
  sf_mex_destroy(&c8_rhs27);
  sf_mex_destroy(&c8_lhs27);
  sf_mex_destroy(&c8_rhs28);
  sf_mex_destroy(&c8_lhs28);
  sf_mex_destroy(&c8_rhs29);
  sf_mex_destroy(&c8_lhs29);
  sf_mex_destroy(&c8_rhs30);
  sf_mex_destroy(&c8_lhs30);
  sf_mex_destroy(&c8_rhs31);
  sf_mex_destroy(&c8_lhs31);
  sf_mex_destroy(&c8_rhs32);
  sf_mex_destroy(&c8_lhs32);
  sf_mex_destroy(&c8_rhs33);
  sf_mex_destroy(&c8_lhs33);
  sf_mex_destroy(&c8_rhs34);
  sf_mex_destroy(&c8_lhs34);
}

static const mxArray *c8_emlrt_marshallOut(const char * c8_u)
{
  const mxArray *c8_y = NULL;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", c8_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c8_u)), false);
  return c8_y;
}

static const mxArray *c8_b_emlrt_marshallOut(const uint32_T c8_u)
{
  const mxArray *c8_y = NULL;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", &c8_u, 7, 0U, 0U, 0U, 0), false);
  return c8_y;
}

static real_T c8_eml_xnrm2(SFc8_test_hill_simInstanceStruct *chartInstance,
  real_T c8_x[3])
{
  real_T c8_y;
  real_T c8_scale;
  int32_T c8_k;
  int32_T c8_b_k;
  real_T c8_b_x;
  real_T c8_c_x;
  real_T c8_absxk;
  real_T c8_t;
  c8_below_threshold(chartInstance);
  c8_y = 0.0;
  c8_scale = 2.2250738585072014E-308;
  for (c8_k = 1; c8_k < 4; c8_k++) {
    c8_b_k = c8_k;
    c8_b_x = c8_x[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c8_b_k), 1, 3, 1, 0) - 1];
    c8_c_x = c8_b_x;
    c8_absxk = muDoubleScalarAbs(c8_c_x);
    if (c8_absxk > c8_scale) {
      c8_t = c8_scale / c8_absxk;
      c8_y = 1.0 + c8_y * c8_t * c8_t;
      c8_scale = c8_absxk;
    } else {
      c8_t = c8_absxk / c8_scale;
      c8_y += c8_t * c8_t;
    }
  }

  return c8_scale * muDoubleScalarSqrt(c8_y);
}

static void c8_below_threshold(SFc8_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *c8_d_sf_marshallOut(void *chartInstanceVoid, void
  *c8_inData)
{
  const mxArray *c8_mxArrayOutData = NULL;
  int32_T c8_u;
  const mxArray *c8_y = NULL;
  SFc8_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc8_test_hill_simInstanceStruct *)chartInstanceVoid;
  c8_mxArrayOutData = NULL;
  c8_u = *(int32_T *)c8_inData;
  c8_y = NULL;
  sf_mex_assign(&c8_y, sf_mex_create("y", &c8_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c8_mxArrayOutData, c8_y, false);
  return c8_mxArrayOutData;
}

static int32_T c8_d_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId)
{
  int32_T c8_y;
  int32_T c8_i32;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), &c8_i32, 1, 6, 0U, 0, 0U, 0);
  c8_y = c8_i32;
  sf_mex_destroy(&c8_u);
  return c8_y;
}

static void c8_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c8_mxArrayInData, const char_T *c8_varName, void *c8_outData)
{
  const mxArray *c8_b_sfEvent;
  const char_T *c8_identifier;
  emlrtMsgIdentifier c8_thisId;
  int32_T c8_y;
  SFc8_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc8_test_hill_simInstanceStruct *)chartInstanceVoid;
  c8_b_sfEvent = sf_mex_dup(c8_mxArrayInData);
  c8_identifier = c8_varName;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_y = c8_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c8_b_sfEvent),
    &c8_thisId);
  sf_mex_destroy(&c8_b_sfEvent);
  *(int32_T *)c8_outData = c8_y;
  sf_mex_destroy(&c8_mxArrayInData);
}

static uint8_T c8_e_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_b_is_active_c8_test_hill_sim, const char_T
  *c8_identifier)
{
  uint8_T c8_y;
  emlrtMsgIdentifier c8_thisId;
  c8_thisId.fIdentifier = c8_identifier;
  c8_thisId.fParent = NULL;
  c8_y = c8_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c8_b_is_active_c8_test_hill_sim), &c8_thisId);
  sf_mex_destroy(&c8_b_is_active_c8_test_hill_sim);
  return c8_y;
}

static uint8_T c8_f_emlrt_marshallIn(SFc8_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c8_u, const emlrtMsgIdentifier *c8_parentId)
{
  uint8_T c8_y;
  uint8_T c8_u0;
  (void)chartInstance;
  sf_mex_import(c8_parentId, sf_mex_dup(c8_u), &c8_u0, 1, 3, 0U, 0, 0U, 0);
  c8_y = c8_u0;
  sf_mex_destroy(&c8_u);
  return c8_y;
}

static void init_dsm_address_info(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc8_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c8_R_Sun_c = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c8_R_c = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c8_e_f = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 2);
  chartInstance->c8_a_Sun = (real_T (*)[3])ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c8_m_cubesat = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 3);
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c8_test_hill_sim_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(681416035U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1170651004U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(144559455U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2908331053U);
}

mxArray* sf_c8_test_hill_sim_get_post_codegen_info(void);
mxArray *sf_c8_test_hill_sim_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("Mj5HKuTYlGwOfi1Ziet7uB");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  {
    mxArray* mxPostCodegenInfo = sf_c8_test_hill_sim_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c8_test_hill_sim_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c8_test_hill_sim_jit_fallback_info(void)
{
  const char *infoFields[] = { "fallbackType", "fallbackReason",
    "incompatibleSymbol", };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 3, infoFields);
  mxArray *fallbackReason = mxCreateString("feature_off");
  mxArray *incompatibleSymbol = mxCreateString("");
  mxArray *fallbackType = mxCreateString("early");
  mxSetField(mxInfo, 0, infoFields[0], fallbackType);
  mxSetField(mxInfo, 0, infoFields[1], fallbackReason);
  mxSetField(mxInfo, 0, infoFields[2], incompatibleSymbol);
  return mxInfo;
}

mxArray *sf_c8_test_hill_sim_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c8_test_hill_sim_get_post_codegen_info(void)
{
  const char* fieldNames[] = { "exportedFunctionsUsedByThisChart",
    "exportedFunctionsChecksum" };

  mwSize dims[2] = { 1, 1 };

  mxArray* mxPostCodegenInfo = mxCreateStructArray(2, dims, sizeof(fieldNames)/
    sizeof(fieldNames[0]), fieldNames);

  {
    mxArray* mxExportedFunctionsChecksum = mxCreateString("");
    mwSize exp_dims[2] = { 0, 1 };

    mxArray* mxExportedFunctionsUsedByThisChart = mxCreateCellArray(2, exp_dims);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsUsedByThisChart",
               mxExportedFunctionsUsedByThisChart);
    mxSetField(mxPostCodegenInfo, 0, "exportedFunctionsChecksum",
               mxExportedFunctionsChecksum);
  }

  return mxPostCodegenInfo;
}

static const mxArray *sf_get_sim_state_info_c8_test_hill_sim(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"a_Sun\",},{M[8],M[0],T\"is_active_c8_test_hill_sim\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c8_test_hill_sim_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc8_test_hill_simInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc8_test_hill_simInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _test_hill_simMachineNumber_,
           8,
           1,
           1,
           0,
           5,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           (void *)S);

        /* Each instance must initialize its own list of scripts */
        init_script_number_translation(_test_hill_simMachineNumber_,
          chartInstance->chartNumber,chartInstance->instanceNumber);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_test_hill_simMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _test_hill_simMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"R_Sun_c");
          _SFD_SET_DATA_PROPS(1,1,1,0,"R_c");
          _SFD_SET_DATA_PROPS(2,1,1,0,"e_f");
          _SFD_SET_DATA_PROPS(3,2,0,1,"a_Sun");
          _SFD_SET_DATA_PROPS(4,1,1,0,"m_cubesat");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,549);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c8_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c8_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c8_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c8_sf_marshallOut,(MexInFcnForType)
            c8_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c8_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_VALUE_PTR(0U, *chartInstance->c8_R_Sun_c);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c8_R_c);
        _SFD_SET_DATA_VALUE_PTR(2U, chartInstance->c8_e_f);
        _SFD_SET_DATA_VALUE_PTR(3U, *chartInstance->c8_a_Sun);
        _SFD_SET_DATA_VALUE_PTR(4U, chartInstance->c8_m_cubesat);
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _test_hill_simMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "sVXdLiO3cmBabCb4yoSmzG";
}

static void sf_opaque_initialize_c8_test_hill_sim(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc8_test_hill_simInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
    chartInstanceVar);
  initialize_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c8_test_hill_sim(void *chartInstanceVar)
{
  enable_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c8_test_hill_sim(void *chartInstanceVar)
{
  disable_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c8_test_hill_sim(void *chartInstanceVar)
{
  sf_gateway_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c8_test_hill_sim(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c8_test_hill_sim(SimStruct* S, const mxArray
  *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c8_test_hill_sim(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc8_test_hill_simInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_hill_sim_optimization_info();
    }

    finalize_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
      chartInstanceVar);
    utFree(chartInstanceVar);
    if (crtInfo != NULL) {
      utFree(crtInfo);
    }

    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c8_test_hill_sim(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    initialize_params_c8_test_hill_sim((SFc8_test_hill_simInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c8_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_hill_sim_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,8);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,8,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,8,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,8);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,8,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,8,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 4; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,8);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(4148810857U));
  ssSetChecksum1(S,(3979782654U));
  ssSetChecksum2(S,(25924162U));
  ssSetChecksum3(S,(1504943838U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c8_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c8_test_hill_sim(SimStruct *S)
{
  SFc8_test_hill_simInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc8_test_hill_simInstanceStruct *)utMalloc(sizeof
    (SFc8_test_hill_simInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc8_test_hill_simInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c8_test_hill_sim;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c8_test_hill_sim;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c8_test_hill_sim;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c8_test_hill_sim;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c8_test_hill_sim;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c8_test_hill_sim;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c8_test_hill_sim;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c8_test_hill_sim;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c8_test_hill_sim;
  chartInstance->chartInfo.mdlStart = mdlStart_c8_test_hill_sim;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c8_test_hill_sim;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->chartInfo.callAtomicSubchartUserFcn = NULL;
  chartInstance->chartInfo.callAtomicSubchartAutoFcn = NULL;
  chartInstance->chartInfo.debugInstance = sfGlobalDebugInstanceStruct;
  chartInstance->S = S;
  crtInfo->checksum = SF_RUNTIME_INFO_CHECKSUM;
  crtInfo->instanceInfo = (&(chartInstance->chartInfo));
  crtInfo->isJITEnabled = false;
  crtInfo->compiledInfo = NULL;
  ssSetUserData(S,(void *)(crtInfo));  /* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  init_simulink_io_address(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c8_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c8_test_hill_sim(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c8_test_hill_sim(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c8_test_hill_sim(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c8_test_hill_sim_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
