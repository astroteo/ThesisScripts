/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_hill_sim_sfun.h"
#include "c13_test_hill_sim.h"
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
static const char * c13_debug_family_names[12] = { "h", "Z_km", "d_km", "M",
  "h_top", "i_top", "i", "h_test", "nargin", "nargout", "h_orb", "d" };

/* Function Declarations */
static void initialize_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void initialize_params_c13_test_hill_sim
  (SFc13_test_hill_simInstanceStruct *chartInstance);
static void enable_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void disable_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void c13_update_debugger_state_c13_test_hill_sim
  (SFc13_test_hill_simInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c13_test_hill_sim
  (SFc13_test_hill_simInstanceStruct *chartInstance);
static void set_sim_state_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_st);
static void finalize_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void sf_gateway_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void mdl_start_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void initSimStructsc13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c13_machineNumber, uint32_T
  c13_chartNumber, uint32_T c13_instanceNumber);
static const mxArray *c13_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static real_T c13_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_b_d, const char_T *c13_identifier);
static real_T c13_b_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId);
static void c13_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static const mxArray *c13_b_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static void c13_c_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId,
  real_T c13_y[44]);
static void c13_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static void c13_info_helper(const mxArray **c13_info);
static const mxArray *c13_emlrt_marshallOut(const char * c13_u);
static const mxArray *c13_b_emlrt_marshallOut(const uint32_T c13_u);
static const mxArray *c13_c_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData);
static int32_T c13_d_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId);
static void c13_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData);
static uint8_T c13_e_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_b_is_active_c13_test_hill_sim, const char_T
  *c13_identifier);
static uint8_T c13_f_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId);
static void init_dsm_address_info(SFc13_test_hill_simInstanceStruct
  *chartInstance);
static void init_simulink_io_address(SFc13_test_hill_simInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c13_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c13_is_active_c13_test_hill_sim = 0U;
}

static void initialize_params_c13_test_hill_sim
  (SFc13_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c13_update_debugger_state_c13_test_hill_sim
  (SFc13_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c13_test_hill_sim
  (SFc13_test_hill_simInstanceStruct *chartInstance)
{
  const mxArray *c13_st;
  const mxArray *c13_y = NULL;
  real_T c13_hoistedGlobal;
  real_T c13_u;
  const mxArray *c13_b_y = NULL;
  uint8_T c13_b_hoistedGlobal;
  uint8_T c13_b_u;
  const mxArray *c13_c_y = NULL;
  c13_st = NULL;
  c13_st = NULL;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_createcellmatrix(2, 1), false);
  c13_hoistedGlobal = *chartInstance->c13_d;
  c13_u = c13_hoistedGlobal;
  c13_b_y = NULL;
  sf_mex_assign(&c13_b_y, sf_mex_create("y", &c13_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c13_y, 0, c13_b_y);
  c13_b_hoistedGlobal = chartInstance->c13_is_active_c13_test_hill_sim;
  c13_b_u = c13_b_hoistedGlobal;
  c13_c_y = NULL;
  sf_mex_assign(&c13_c_y, sf_mex_create("y", &c13_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c13_y, 1, c13_c_y);
  sf_mex_assign(&c13_st, c13_y, false);
  return c13_st;
}

static void set_sim_state_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_st)
{
  const mxArray *c13_u;
  chartInstance->c13_doneDoubleBufferReInit = true;
  c13_u = sf_mex_dup(c13_st);
  *chartInstance->c13_d = c13_emlrt_marshallIn(chartInstance, sf_mex_dup
    (sf_mex_getcell(c13_u, 0)), "d");
  chartInstance->c13_is_active_c13_test_hill_sim = c13_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c13_u, 1)),
     "is_active_c13_test_hill_sim");
  sf_mex_destroy(&c13_u);
  c13_update_debugger_state_c13_test_hill_sim(chartInstance);
  sf_mex_destroy(&c13_st);
}

static void finalize_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  real_T c13_hoistedGlobal;
  real_T c13_b_h_orb;
  uint32_T c13_debug_family_var_map[12];
  real_T c13_h;
  real_T c13_Z_km[44];
  real_T c13_d_km[44];
  real_T c13_M;
  real_T c13_h_top;
  real_T c13_i_top;
  real_T c13_i;
  real_T c13_h_test;
  real_T c13_nargin = 1.0;
  real_T c13_nargout = 1.0;
  real_T c13_b_d;
  real_T c13_A;
  real_T c13_x;
  real_T c13_b_x;
  real_T c13_c_x;
  int32_T c13_i0;
  int32_T c13_i1;
  static real_T c13_dv0[44] = { 0.0948, 0.00407, 0.000331, 1.69E-5, 5.77E-7,
    1.7E-8, 2.96E-9, 9.65E-10, 3.9E-10, 1.75E-10, 8.47E-11, 4.31E-11, 2.3E-11,
    1.27E-11, 7.22E-12, 4.21E-12, 2.5E-12, 1.51E-12, 9.2E-13, 5.68E-13, 3.54E-13,
    2.23E-13, 1.42E-13, 9.2E-14, 6.03E-14, 4.03E-14, 2.75E-14, 1.93E-14,
    1.39E-14, 1.03E-14, 7.9E-15, 6.24E-15, 5.06E-15, 4.21E-15, 3.58E-15,
    3.09E-15, 2.7E-15, 2.39E-15, 2.13E-15, 1.91E-15, 1.73E-15, 1.56E-15,
    1.42E-15, 1.3E-15 };

  int32_T c13_b_i;
  real_T c13_b_A;
  real_T c13_d_x;
  real_T c13_e_x;
  real_T c13_f_x;
  real_T c13_y;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 12U, chartInstance->c13_sfEvent);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c13_h_orb, 0U);
  chartInstance->c13_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 12U, chartInstance->c13_sfEvent);
  c13_hoistedGlobal = *chartInstance->c13_h_orb;
  c13_b_h_orb = c13_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 12U, 12U, c13_debug_family_names,
    c13_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_h, 0U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c13_Z_km, 1U, c13_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c13_d_km, 2U, c13_b_sf_marshallOut,
    c13_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_M, 3U, c13_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_h_top, 4U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_i_top, 5U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_i, 6U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_h_test, 7U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_nargin, 8U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_nargout, 9U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c13_b_h_orb, 10U, c13_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c13_b_d, 11U, c13_sf_marshallOut,
    c13_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 5);
  c13_A = c13_b_h_orb;
  c13_x = c13_A;
  c13_b_x = c13_x;
  c13_c_x = c13_b_x;
  c13_h = c13_c_x / 1000.0;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 7);
  for (c13_i0 = 0; c13_i0 < 44; c13_i0++) {
    c13_Z_km[c13_i0] = 20.0 + 20.0 * (real_T)c13_i0;
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 8);
  for (c13_i1 = 0; c13_i1 < 44; c13_i1++) {
    c13_d_km[c13_i1] = c13_dv0[c13_i1];
  }

  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 9);
  c13_M = 44.0;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 9);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 11);
  c13_h_top = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 12);
  c13_i_top = 2.0;
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 13);
  c13_i = 2.0;
  c13_b_i = 0;
  while (c13_b_i < 43) {
    c13_i = 2.0 + (real_T)c13_b_i;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 14);
    c13_h_test = c13_Z_km[_SFD_EML_ARRAY_BOUNDS_CHECK("Z_km", (int32_T)
      _SFD_INTEGER_CHECK("i", c13_i), 1, 44, 1, 0) - 1];
    _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 16);
    if (CV_EML_IF(0, 1, 0, CV_RELATIONAL_EVAL(4U, 0U, 0, c13_h_test, c13_h, -1,
          4U, c13_h_test > c13_h))) {
      _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 18);
      c13_h_top = c13_Z_km[_SFD_EML_ARRAY_BOUNDS_CHECK("Z_km", (int32_T)
        _SFD_INTEGER_CHECK("i", c13_i), 1, 44, 1, 0) - 1];
      _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 19);
      c13_i_top = c13_i;
    }

    c13_b_i++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, 25);
  c13_b_A = c13_d_km[_SFD_EML_ARRAY_BOUNDS_CHECK("d_km", (int32_T)
    _SFD_INTEGER_CHECK("i_top", c13_i_top), 1, 44, 1, 0) - 1] -
    c13_d_km[_SFD_EML_ARRAY_BOUNDS_CHECK("d_km", (int32_T)_SFD_INTEGER_CHECK(
    "i_top-1", c13_i_top - 1.0), 1, 44, 1, 0) - 1];
  c13_d_x = c13_b_A;
  c13_e_x = c13_d_x;
  c13_f_x = c13_e_x;
  c13_y = c13_f_x / 20.0;
  c13_b_d = c13_d_km[_SFD_EML_ARRAY_BOUNDS_CHECK("d_km", (int32_T)
    _SFD_INTEGER_CHECK("i_top-1", c13_i_top - 1.0), 1, 44, 1, 0) - 1] + c13_y *
    (c13_h - c13_Z_km[_SFD_EML_ARRAY_BOUNDS_CHECK("Z_km", (int32_T)
      _SFD_INTEGER_CHECK("i_top-1", c13_i_top - 1.0), 1, 44, 1, 0) - 1]);
  _SFD_EML_CALL(0U, chartInstance->c13_sfEvent, -25);
  _SFD_SYMBOL_SCOPE_POP();
  *chartInstance->c13_d = c13_b_d;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 12U, chartInstance->c13_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_test_hill_simMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c13_d, 1U);
}

static void mdl_start_c13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void initSimStructsc13_test_hill_sim(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c13_machineNumber, uint32_T
  c13_chartNumber, uint32_T c13_instanceNumber)
{
  (void)c13_machineNumber;
  (void)c13_chartNumber;
  (void)c13_instanceNumber;
}

static const mxArray *c13_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  real_T c13_u;
  const mxArray *c13_y = NULL;
  SFc13_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc13_test_hill_simInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_u = *(real_T *)c13_inData;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", &c13_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, false);
  return c13_mxArrayOutData;
}

static real_T c13_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_b_d, const char_T *c13_identifier)
{
  real_T c13_y;
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_d), &c13_thisId);
  sf_mex_destroy(&c13_b_d);
  return c13_y;
}

static real_T c13_b_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId)
{
  real_T c13_y;
  real_T c13_d0;
  (void)chartInstance;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), &c13_d0, 1, 0, 0U, 0, 0U, 0);
  c13_y = c13_d0;
  sf_mex_destroy(&c13_u);
  return c13_y;
}

static void c13_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_b_d;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y;
  SFc13_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc13_test_hill_simInstanceStruct *)chartInstanceVoid;
  c13_b_d = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_d), &c13_thisId);
  sf_mex_destroy(&c13_b_d);
  *(real_T *)c13_outData = c13_y;
  sf_mex_destroy(&c13_mxArrayInData);
}

static const mxArray *c13_b_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_i2;
  real_T c13_b_inData[44];
  int32_T c13_i3;
  real_T c13_u[44];
  const mxArray *c13_y = NULL;
  SFc13_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc13_test_hill_simInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  for (c13_i2 = 0; c13_i2 < 44; c13_i2++) {
    c13_b_inData[c13_i2] = (*(real_T (*)[44])c13_inData)[c13_i2];
  }

  for (c13_i3 = 0; c13_i3 < 44; c13_i3++) {
    c13_u[c13_i3] = c13_b_inData[c13_i3];
  }

  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 0, 0U, 1U, 0U, 2, 1, 44),
                false);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, false);
  return c13_mxArrayOutData;
}

static void c13_c_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId,
  real_T c13_y[44])
{
  real_T c13_dv1[44];
  int32_T c13_i4;
  (void)chartInstance;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), c13_dv1, 1, 0, 0U, 1, 0U, 2, 1,
                44);
  for (c13_i4 = 0; c13_i4 < 44; c13_i4++) {
    c13_y[c13_i4] = c13_dv1[c13_i4];
  }

  sf_mex_destroy(&c13_u);
}

static void c13_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_d_km;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  real_T c13_y[44];
  int32_T c13_i5;
  SFc13_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc13_test_hill_simInstanceStruct *)chartInstanceVoid;
  c13_d_km = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_d_km), &c13_thisId, c13_y);
  sf_mex_destroy(&c13_d_km);
  for (c13_i5 = 0; c13_i5 < 44; c13_i5++) {
    (*(real_T (*)[44])c13_outData)[c13_i5] = c13_y[c13_i5];
  }

  sf_mex_destroy(&c13_mxArrayInData);
}

const mxArray *sf_c13_test_hill_sim_get_eml_resolved_functions_info(void)
{
  const mxArray *c13_nameCaptureInfo = NULL;
  c13_nameCaptureInfo = NULL;
  sf_mex_assign(&c13_nameCaptureInfo, sf_mex_createstruct("structure", 2, 7, 1),
                false);
  c13_info_helper(&c13_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c13_nameCaptureInfo);
  return c13_nameCaptureInfo;
}

static void c13_info_helper(const mxArray **c13_info)
{
  const mxArray *c13_rhs0 = NULL;
  const mxArray *c13_lhs0 = NULL;
  const mxArray *c13_rhs1 = NULL;
  const mxArray *c13_lhs1 = NULL;
  const mxArray *c13_rhs2 = NULL;
  const mxArray *c13_lhs2 = NULL;
  const mxArray *c13_rhs3 = NULL;
  const mxArray *c13_lhs3 = NULL;
  const mxArray *c13_rhs4 = NULL;
  const mxArray *c13_lhs4 = NULL;
  const mxArray *c13_rhs5 = NULL;
  const mxArray *c13_lhs5 = NULL;
  const mxArray *c13_rhs6 = NULL;
  const mxArray *c13_lhs6 = NULL;
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("mrdivide"), "name", "name",
                  0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1410807648U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1370009886U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c13_rhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c13_lhs0, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.assert"),
                  "name", "name", 1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/assert.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1389717774U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c13_rhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c13_lhs1, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("rdivide"), "name", "name", 2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1363713880U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c13_rhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c13_lhs2, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1395931856U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c13_rhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c13_lhs3, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1286818796U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c13_rhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c13_lhs4, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("eml_div"), "name", "name", 5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1386423952U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c13_rhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c13_lhs5, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("coder.internal.div"), "name",
                  "name", 6);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c13_info, c13_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/coder/coder/+coder/+internal/div.p"), "resolved",
                  "resolved", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(1410807770U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c13_info, c13_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c13_rhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_assign(&c13_lhs6, sf_mex_createcellmatrix(0, 1), false);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c13_info, sf_mex_duplicatearraysafe(&c13_lhs6), "lhs", "lhs",
                  6);
  sf_mex_destroy(&c13_rhs0);
  sf_mex_destroy(&c13_lhs0);
  sf_mex_destroy(&c13_rhs1);
  sf_mex_destroy(&c13_lhs1);
  sf_mex_destroy(&c13_rhs2);
  sf_mex_destroy(&c13_lhs2);
  sf_mex_destroy(&c13_rhs3);
  sf_mex_destroy(&c13_lhs3);
  sf_mex_destroy(&c13_rhs4);
  sf_mex_destroy(&c13_lhs4);
  sf_mex_destroy(&c13_rhs5);
  sf_mex_destroy(&c13_lhs5);
  sf_mex_destroy(&c13_rhs6);
  sf_mex_destroy(&c13_lhs6);
}

static const mxArray *c13_emlrt_marshallOut(const char * c13_u)
{
  const mxArray *c13_y = NULL;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", c13_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c13_u)), false);
  return c13_y;
}

static const mxArray *c13_b_emlrt_marshallOut(const uint32_T c13_u)
{
  const mxArray *c13_y = NULL;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", &c13_u, 7, 0U, 0U, 0U, 0), false);
  return c13_y;
}

static const mxArray *c13_c_sf_marshallOut(void *chartInstanceVoid, void
  *c13_inData)
{
  const mxArray *c13_mxArrayOutData = NULL;
  int32_T c13_u;
  const mxArray *c13_y = NULL;
  SFc13_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc13_test_hill_simInstanceStruct *)chartInstanceVoid;
  c13_mxArrayOutData = NULL;
  c13_u = *(int32_T *)c13_inData;
  c13_y = NULL;
  sf_mex_assign(&c13_y, sf_mex_create("y", &c13_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c13_mxArrayOutData, c13_y, false);
  return c13_mxArrayOutData;
}

static int32_T c13_d_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId)
{
  int32_T c13_y;
  int32_T c13_i6;
  (void)chartInstance;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), &c13_i6, 1, 6, 0U, 0, 0U, 0);
  c13_y = c13_i6;
  sf_mex_destroy(&c13_u);
  return c13_y;
}

static void c13_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c13_mxArrayInData, const char_T *c13_varName, void *c13_outData)
{
  const mxArray *c13_b_sfEvent;
  const char_T *c13_identifier;
  emlrtMsgIdentifier c13_thisId;
  int32_T c13_y;
  SFc13_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc13_test_hill_simInstanceStruct *)chartInstanceVoid;
  c13_b_sfEvent = sf_mex_dup(c13_mxArrayInData);
  c13_identifier = c13_varName;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c13_b_sfEvent),
    &c13_thisId);
  sf_mex_destroy(&c13_b_sfEvent);
  *(int32_T *)c13_outData = c13_y;
  sf_mex_destroy(&c13_mxArrayInData);
}

static uint8_T c13_e_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_b_is_active_c13_test_hill_sim, const char_T
  *c13_identifier)
{
  uint8_T c13_y;
  emlrtMsgIdentifier c13_thisId;
  c13_thisId.fIdentifier = c13_identifier;
  c13_thisId.fParent = NULL;
  c13_y = c13_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c13_b_is_active_c13_test_hill_sim), &c13_thisId);
  sf_mex_destroy(&c13_b_is_active_c13_test_hill_sim);
  return c13_y;
}

static uint8_T c13_f_emlrt_marshallIn(SFc13_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c13_u, const emlrtMsgIdentifier *c13_parentId)
{
  uint8_T c13_y;
  uint8_T c13_u0;
  (void)chartInstance;
  sf_mex_import(c13_parentId, sf_mex_dup(c13_u), &c13_u0, 1, 3, 0U, 0, 0U, 0);
  c13_y = c13_u0;
  sf_mex_destroy(&c13_u);
  return c13_y;
}

static void init_dsm_address_info(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc13_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c13_h_orb = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c13_d = (real_T *)ssGetOutputPortSignal_wrapper
    (chartInstance->S, 1);
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

void sf_c13_test_hill_sim_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(961822236U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3867818094U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1101675715U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(810680582U);
}

mxArray* sf_c13_test_hill_sim_get_post_codegen_info(void);
mxArray *sf_c13_test_hill_sim_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("YSdLoE9uyGdH1Q5rvyQYMB");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
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
      pr[0] = (double)(1);
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
    mxArray* mxPostCodegenInfo = sf_c13_test_hill_sim_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c13_test_hill_sim_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c13_test_hill_sim_jit_fallback_info(void)
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

mxArray *sf_c13_test_hill_sim_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c13_test_hill_sim_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c13_test_hill_sim(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[6],T\"d\",},{M[8],M[0],T\"is_active_c13_test_hill_sim\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c13_test_hill_sim_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc13_test_hill_simInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc13_test_hill_simInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _test_hill_simMachineNumber_,
           13,
           1,
           1,
           0,
           2,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"h_orb");
          _SFD_SET_DATA_PROPS(1,2,0,1,"d");
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
        _SFD_CV_INIT_EML(0,1,1,1,0,0,0,1,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",2,-1,955);
        _SFD_CV_INIT_EML_IF(0,1,0,794,807,-1,865);
        _SFD_CV_INIT_EML_FOR(0,1,0,755,765,878);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,0,797,807,-1,4);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c13_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c13_sf_marshallOut,(MexInFcnForType)c13_sf_marshallIn);
        _SFD_SET_DATA_VALUE_PTR(0U, chartInstance->c13_h_orb);
        _SFD_SET_DATA_VALUE_PTR(1U, chartInstance->c13_d);
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
  return "W5KsT1MrGCMRR4kZKqcWCG";
}

static void sf_opaque_initialize_c13_test_hill_sim(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc13_test_hill_simInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
    chartInstanceVar);
  initialize_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c13_test_hill_sim(void *chartInstanceVar)
{
  enable_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c13_test_hill_sim(void *chartInstanceVar)
{
  disable_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c13_test_hill_sim(void *chartInstanceVar)
{
  sf_gateway_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c13_test_hill_sim(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c13_test_hill_sim(SimStruct* S, const
  mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c13_test_hill_sim(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc13_test_hill_simInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_hill_sim_optimization_info();
    }

    finalize_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
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
  initSimStructsc13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c13_test_hill_sim(SimStruct *S)
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
    initialize_params_c13_test_hill_sim((SFc13_test_hill_simInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c13_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_hill_sim_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      13);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,13,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,13,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,13);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,13,1);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,13,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 1; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,13);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(4191620109U));
  ssSetChecksum1(S,(3701427927U));
  ssSetChecksum2(S,(1430838932U));
  ssSetChecksum3(S,(4248730034U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c13_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c13_test_hill_sim(SimStruct *S)
{
  SFc13_test_hill_simInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc13_test_hill_simInstanceStruct *)utMalloc(sizeof
    (SFc13_test_hill_simInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc13_test_hill_simInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c13_test_hill_sim;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c13_test_hill_sim;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c13_test_hill_sim;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c13_test_hill_sim;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c13_test_hill_sim;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c13_test_hill_sim;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c13_test_hill_sim;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c13_test_hill_sim;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c13_test_hill_sim;
  chartInstance->chartInfo.mdlStart = mdlStart_c13_test_hill_sim;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c13_test_hill_sim;
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

void c13_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c13_test_hill_sim(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c13_test_hill_sim(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c13_test_hill_sim(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c13_test_hill_sim_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
