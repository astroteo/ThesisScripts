/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_hill_sim_sfun.h"
#include "c15_test_hill_sim.h"
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
static const char * c15_debug_family_names[6] = { "j", "nargin", "nargout",
  "a_sat", "a_c", "a_c_SAT" };

/* Function Declarations */
static void initialize_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void initialize_params_c15_test_hill_sim
  (SFc15_test_hill_simInstanceStruct *chartInstance);
static void enable_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void disable_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void c15_update_debugger_state_c15_test_hill_sim
  (SFc15_test_hill_simInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c15_test_hill_sim
  (SFc15_test_hill_simInstanceStruct *chartInstance);
static void set_sim_state_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_st);
static void finalize_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void sf_gateway_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void mdl_start_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void initSimStructsc15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c15_machineNumber, uint32_T
  c15_chartNumber, uint32_T c15_instanceNumber);
static const mxArray *c15_sf_marshallOut(void *chartInstanceVoid, void
  *c15_inData);
static void c15_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_b_a_c_SAT, const char_T *c15_identifier,
  real_T c15_y[3]);
static void c15_b_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId,
  real_T c15_y[3]);
static void c15_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c15_mxArrayInData, const char_T *c15_varName, void *c15_outData);
static const mxArray *c15_b_sf_marshallOut(void *chartInstanceVoid, void
  *c15_inData);
static real_T c15_c_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId);
static void c15_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c15_mxArrayInData, const char_T *c15_varName, void *c15_outData);
static const mxArray *c15_c_sf_marshallOut(void *chartInstanceVoid, void
  *c15_inData);
static int32_T c15_d_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId);
static void c15_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c15_mxArrayInData, const char_T *c15_varName, void *c15_outData);
static uint8_T c15_e_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_b_is_active_c15_test_hill_sim, const char_T
  *c15_identifier);
static uint8_T c15_f_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId);
static void init_dsm_address_info(SFc15_test_hill_simInstanceStruct
  *chartInstance);
static void init_simulink_io_address(SFc15_test_hill_simInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c15_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c15_is_active_c15_test_hill_sim = 0U;
}

static void initialize_params_c15_test_hill_sim
  (SFc15_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c15_update_debugger_state_c15_test_hill_sim
  (SFc15_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c15_test_hill_sim
  (SFc15_test_hill_simInstanceStruct *chartInstance)
{
  const mxArray *c15_st;
  const mxArray *c15_y = NULL;
  int32_T c15_i0;
  real_T c15_u[3];
  const mxArray *c15_b_y = NULL;
  uint8_T c15_hoistedGlobal;
  uint8_T c15_b_u;
  const mxArray *c15_c_y = NULL;
  c15_st = NULL;
  c15_st = NULL;
  c15_y = NULL;
  sf_mex_assign(&c15_y, sf_mex_createcellmatrix(2, 1), false);
  for (c15_i0 = 0; c15_i0 < 3; c15_i0++) {
    c15_u[c15_i0] = (*chartInstance->c15_a_c_SAT)[c15_i0];
  }

  c15_b_y = NULL;
  sf_mex_assign(&c15_b_y, sf_mex_create("y", c15_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c15_y, 0, c15_b_y);
  c15_hoistedGlobal = chartInstance->c15_is_active_c15_test_hill_sim;
  c15_b_u = c15_hoistedGlobal;
  c15_c_y = NULL;
  sf_mex_assign(&c15_c_y, sf_mex_create("y", &c15_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c15_y, 1, c15_c_y);
  sf_mex_assign(&c15_st, c15_y, false);
  return c15_st;
}

static void set_sim_state_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_st)
{
  const mxArray *c15_u;
  real_T c15_dv0[3];
  int32_T c15_i1;
  chartInstance->c15_doneDoubleBufferReInit = true;
  c15_u = sf_mex_dup(c15_st);
  c15_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c15_u, 0)),
                       "a_c_SAT", c15_dv0);
  for (c15_i1 = 0; c15_i1 < 3; c15_i1++) {
    (*chartInstance->c15_a_c_SAT)[c15_i1] = c15_dv0[c15_i1];
  }

  chartInstance->c15_is_active_c15_test_hill_sim = c15_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c15_u, 1)),
     "is_active_c15_test_hill_sim");
  sf_mex_destroy(&c15_u);
  c15_update_debugger_state_c15_test_hill_sim(chartInstance);
  sf_mex_destroy(&c15_st);
}

static void finalize_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  int32_T c15_i2;
  real_T c15_hoistedGlobal;
  real_T c15_b_a_sat;
  int32_T c15_i3;
  real_T c15_b_a_c[3];
  uint32_T c15_debug_family_var_map[6];
  real_T c15_j;
  real_T c15_nargin = 2.0;
  real_T c15_nargout = 1.0;
  real_T c15_b_a_c_SAT[3];
  int32_T c15_i4;
  int32_T c15_b_j;
  real_T c15_d0;
  real_T c15_d1;
  int32_T c15_i5;
  int32_T c15_i6;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 14U, chartInstance->c15_sfEvent);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c15_a_sat, 0U);
  for (c15_i2 = 0; c15_i2 < 3; c15_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c15_a_c)[c15_i2], 1U);
  }

  chartInstance->c15_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 14U, chartInstance->c15_sfEvent);
  c15_hoistedGlobal = *chartInstance->c15_a_sat;
  c15_b_a_sat = c15_hoistedGlobal;
  for (c15_i3 = 0; c15_i3 < 3; c15_i3++) {
    c15_b_a_c[c15_i3] = (*chartInstance->c15_a_c)[c15_i3];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c15_debug_family_names,
    c15_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c15_j, 0U, c15_b_sf_marshallOut,
    c15_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c15_nargin, 1U, c15_b_sf_marshallOut,
    c15_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c15_nargout, 2U, c15_b_sf_marshallOut,
    c15_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c15_b_a_sat, 3U, c15_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c15_b_a_c, 4U, c15_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c15_b_a_c_SAT, 5U, c15_sf_marshallOut,
    c15_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, 3);
  for (c15_i4 = 0; c15_i4 < 3; c15_i4++) {
    c15_b_a_c_SAT[c15_i4] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, 4);
  c15_j = 1.0;
  c15_b_j = 0;
  while (c15_b_j < 3) {
    c15_j = 1.0 + (real_T)c15_b_j;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, 6);
    c15_d0 = c15_b_a_c[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c", (int32_T)
      _SFD_INTEGER_CHECK("j", c15_j), 1, 3, 1, 0) - 1];
    if (CV_EML_IF(0, 1, 0, CV_RELATIONAL_EVAL(4U, 0U, 0, c15_d0, c15_b_a_sat, -1,
          4U, c15_d0 > c15_b_a_sat))) {
      _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, 8);
      c15_b_a_c_SAT[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c_SAT", (int32_T)
        _SFD_INTEGER_CHECK("j", c15_j), 1, 3, 1, 0) - 1] = c15_b_a_sat;
    } else {
      _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, 10);
      c15_d1 = c15_b_a_c[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c", (int32_T)
        _SFD_INTEGER_CHECK("j", c15_j), 1, 3, 1, 0) - 1];
      if (CV_EML_IF(0, 1, 1, CV_RELATIONAL_EVAL(4U, 0U, 1, c15_d1, -c15_b_a_sat,
            -1, 2U, c15_d1 < -c15_b_a_sat))) {
        _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, 12);
        c15_b_a_c_SAT[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c_SAT", (int32_T)
          _SFD_INTEGER_CHECK("j", c15_j), 1, 3, 1, 0) - 1] = -c15_b_a_sat;
      } else {
        _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, 14);
        c15_b_a_c_SAT[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c_SAT", (int32_T)
          _SFD_INTEGER_CHECK("j", c15_j), 1, 3, 1, 0) - 1] =
          c15_b_a_c[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c", (int32_T)
          _SFD_INTEGER_CHECK("j", c15_j), 1, 3, 1, 0) - 1];
      }
    }

    c15_b_j++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c15_sfEvent, -14);
  _SFD_SYMBOL_SCOPE_POP();
  for (c15_i5 = 0; c15_i5 < 3; c15_i5++) {
    (*chartInstance->c15_a_c_SAT)[c15_i5] = c15_b_a_c_SAT[c15_i5];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 14U, chartInstance->c15_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_test_hill_simMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c15_i6 = 0; c15_i6 < 3; c15_i6++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c15_a_c_SAT)[c15_i6], 2U);
  }
}

static void mdl_start_c15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void initSimStructsc15_test_hill_sim(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c15_machineNumber, uint32_T
  c15_chartNumber, uint32_T c15_instanceNumber)
{
  (void)c15_machineNumber;
  (void)c15_chartNumber;
  (void)c15_instanceNumber;
}

static const mxArray *c15_sf_marshallOut(void *chartInstanceVoid, void
  *c15_inData)
{
  const mxArray *c15_mxArrayOutData = NULL;
  int32_T c15_i7;
  real_T c15_b_inData[3];
  int32_T c15_i8;
  real_T c15_u[3];
  const mxArray *c15_y = NULL;
  SFc15_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc15_test_hill_simInstanceStruct *)chartInstanceVoid;
  c15_mxArrayOutData = NULL;
  for (c15_i7 = 0; c15_i7 < 3; c15_i7++) {
    c15_b_inData[c15_i7] = (*(real_T (*)[3])c15_inData)[c15_i7];
  }

  for (c15_i8 = 0; c15_i8 < 3; c15_i8++) {
    c15_u[c15_i8] = c15_b_inData[c15_i8];
  }

  c15_y = NULL;
  sf_mex_assign(&c15_y, sf_mex_create("y", c15_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c15_mxArrayOutData, c15_y, false);
  return c15_mxArrayOutData;
}

static void c15_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_b_a_c_SAT, const char_T *c15_identifier,
  real_T c15_y[3])
{
  emlrtMsgIdentifier c15_thisId;
  c15_thisId.fIdentifier = c15_identifier;
  c15_thisId.fParent = NULL;
  c15_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c15_b_a_c_SAT), &c15_thisId,
    c15_y);
  sf_mex_destroy(&c15_b_a_c_SAT);
}

static void c15_b_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId,
  real_T c15_y[3])
{
  real_T c15_dv1[3];
  int32_T c15_i9;
  (void)chartInstance;
  sf_mex_import(c15_parentId, sf_mex_dup(c15_u), c15_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c15_i9 = 0; c15_i9 < 3; c15_i9++) {
    c15_y[c15_i9] = c15_dv1[c15_i9];
  }

  sf_mex_destroy(&c15_u);
}

static void c15_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c15_mxArrayInData, const char_T *c15_varName, void *c15_outData)
{
  const mxArray *c15_b_a_c_SAT;
  const char_T *c15_identifier;
  emlrtMsgIdentifier c15_thisId;
  real_T c15_y[3];
  int32_T c15_i10;
  SFc15_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc15_test_hill_simInstanceStruct *)chartInstanceVoid;
  c15_b_a_c_SAT = sf_mex_dup(c15_mxArrayInData);
  c15_identifier = c15_varName;
  c15_thisId.fIdentifier = c15_identifier;
  c15_thisId.fParent = NULL;
  c15_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c15_b_a_c_SAT), &c15_thisId,
    c15_y);
  sf_mex_destroy(&c15_b_a_c_SAT);
  for (c15_i10 = 0; c15_i10 < 3; c15_i10++) {
    (*(real_T (*)[3])c15_outData)[c15_i10] = c15_y[c15_i10];
  }

  sf_mex_destroy(&c15_mxArrayInData);
}

static const mxArray *c15_b_sf_marshallOut(void *chartInstanceVoid, void
  *c15_inData)
{
  const mxArray *c15_mxArrayOutData = NULL;
  real_T c15_u;
  const mxArray *c15_y = NULL;
  SFc15_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc15_test_hill_simInstanceStruct *)chartInstanceVoid;
  c15_mxArrayOutData = NULL;
  c15_u = *(real_T *)c15_inData;
  c15_y = NULL;
  sf_mex_assign(&c15_y, sf_mex_create("y", &c15_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c15_mxArrayOutData, c15_y, false);
  return c15_mxArrayOutData;
}

static real_T c15_c_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId)
{
  real_T c15_y;
  real_T c15_d2;
  (void)chartInstance;
  sf_mex_import(c15_parentId, sf_mex_dup(c15_u), &c15_d2, 1, 0, 0U, 0, 0U, 0);
  c15_y = c15_d2;
  sf_mex_destroy(&c15_u);
  return c15_y;
}

static void c15_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c15_mxArrayInData, const char_T *c15_varName, void *c15_outData)
{
  const mxArray *c15_nargout;
  const char_T *c15_identifier;
  emlrtMsgIdentifier c15_thisId;
  real_T c15_y;
  SFc15_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc15_test_hill_simInstanceStruct *)chartInstanceVoid;
  c15_nargout = sf_mex_dup(c15_mxArrayInData);
  c15_identifier = c15_varName;
  c15_thisId.fIdentifier = c15_identifier;
  c15_thisId.fParent = NULL;
  c15_y = c15_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c15_nargout),
    &c15_thisId);
  sf_mex_destroy(&c15_nargout);
  *(real_T *)c15_outData = c15_y;
  sf_mex_destroy(&c15_mxArrayInData);
}

const mxArray *sf_c15_test_hill_sim_get_eml_resolved_functions_info(void)
{
  const mxArray *c15_nameCaptureInfo = NULL;
  c15_nameCaptureInfo = NULL;
  sf_mex_assign(&c15_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), false);
  return c15_nameCaptureInfo;
}

static const mxArray *c15_c_sf_marshallOut(void *chartInstanceVoid, void
  *c15_inData)
{
  const mxArray *c15_mxArrayOutData = NULL;
  int32_T c15_u;
  const mxArray *c15_y = NULL;
  SFc15_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc15_test_hill_simInstanceStruct *)chartInstanceVoid;
  c15_mxArrayOutData = NULL;
  c15_u = *(int32_T *)c15_inData;
  c15_y = NULL;
  sf_mex_assign(&c15_y, sf_mex_create("y", &c15_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c15_mxArrayOutData, c15_y, false);
  return c15_mxArrayOutData;
}

static int32_T c15_d_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId)
{
  int32_T c15_y;
  int32_T c15_i11;
  (void)chartInstance;
  sf_mex_import(c15_parentId, sf_mex_dup(c15_u), &c15_i11, 1, 6, 0U, 0, 0U, 0);
  c15_y = c15_i11;
  sf_mex_destroy(&c15_u);
  return c15_y;
}

static void c15_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c15_mxArrayInData, const char_T *c15_varName, void *c15_outData)
{
  const mxArray *c15_b_sfEvent;
  const char_T *c15_identifier;
  emlrtMsgIdentifier c15_thisId;
  int32_T c15_y;
  SFc15_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc15_test_hill_simInstanceStruct *)chartInstanceVoid;
  c15_b_sfEvent = sf_mex_dup(c15_mxArrayInData);
  c15_identifier = c15_varName;
  c15_thisId.fIdentifier = c15_identifier;
  c15_thisId.fParent = NULL;
  c15_y = c15_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c15_b_sfEvent),
    &c15_thisId);
  sf_mex_destroy(&c15_b_sfEvent);
  *(int32_T *)c15_outData = c15_y;
  sf_mex_destroy(&c15_mxArrayInData);
}

static uint8_T c15_e_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_b_is_active_c15_test_hill_sim, const char_T
  *c15_identifier)
{
  uint8_T c15_y;
  emlrtMsgIdentifier c15_thisId;
  c15_thisId.fIdentifier = c15_identifier;
  c15_thisId.fParent = NULL;
  c15_y = c15_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c15_b_is_active_c15_test_hill_sim), &c15_thisId);
  sf_mex_destroy(&c15_b_is_active_c15_test_hill_sim);
  return c15_y;
}

static uint8_T c15_f_emlrt_marshallIn(SFc15_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c15_u, const emlrtMsgIdentifier *c15_parentId)
{
  uint8_T c15_y;
  uint8_T c15_u0;
  (void)chartInstance;
  sf_mex_import(c15_parentId, sf_mex_dup(c15_u), &c15_u0, 1, 3, 0U, 0, 0U, 0);
  c15_y = c15_u0;
  sf_mex_destroy(&c15_u);
  return c15_y;
}

static void init_dsm_address_info(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc15_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c15_a_sat = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c15_a_c = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c15_a_c_SAT = (real_T (*)[3])ssGetOutputPortSignal_wrapper
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

void sf_c15_test_hill_sim_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2545629185U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2839487140U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3530825079U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(193220667U);
}

mxArray* sf_c15_test_hill_sim_get_post_codegen_info(void);
mxArray *sf_c15_test_hill_sim_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("WG6ut1EsSACocMxs7LMPlE");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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
    mxArray* mxPostCodegenInfo = sf_c15_test_hill_sim_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c15_test_hill_sim_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c15_test_hill_sim_jit_fallback_info(void)
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

mxArray *sf_c15_test_hill_sim_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c15_test_hill_sim_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c15_test_hill_sim(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"a_c_SAT\",},{M[8],M[0],T\"is_active_c15_test_hill_sim\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c15_test_hill_sim_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc15_test_hill_simInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc15_test_hill_simInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _test_hill_simMachineNumber_,
           15,
           1,
           1,
           0,
           3,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"a_sat");
          _SFD_SET_DATA_PROPS(1,1,1,0,"a_c");
          _SFD_SET_DATA_PROPS(2,2,0,1,"a_c_SAT");
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
        _SFD_CV_INIT_EML(0,1,1,2,0,0,0,1,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,266);
        _SFD_CV_INIT_EML_IF(0,1,0,90,109,154,257);
        _SFD_CV_INIT_EML_IF(0,1,1,154,177,218,257);
        _SFD_CV_INIT_EML_FOR(0,1,0,71,81,266);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,0,94,109,-1,4);
        _SFD_CV_INIT_EML_RELATIONAL(0,1,1,162,177,-1,2);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c15_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c15_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c15_sf_marshallOut,(MexInFcnForType)
            c15_sf_marshallIn);
        }

        _SFD_SET_DATA_VALUE_PTR(0U, chartInstance->c15_a_sat);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c15_a_c);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c15_a_c_SAT);
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
  return "VK9KEuOAeSTh4h3tQ9nTIF";
}

static void sf_opaque_initialize_c15_test_hill_sim(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc15_test_hill_simInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
    chartInstanceVar);
  initialize_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c15_test_hill_sim(void *chartInstanceVar)
{
  enable_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c15_test_hill_sim(void *chartInstanceVar)
{
  disable_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c15_test_hill_sim(void *chartInstanceVar)
{
  sf_gateway_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c15_test_hill_sim(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c15_test_hill_sim(SimStruct* S, const
  mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c15_test_hill_sim(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc15_test_hill_simInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_hill_sim_optimization_info();
    }

    finalize_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
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
  initSimStructsc15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c15_test_hill_sim(SimStruct *S)
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
    initialize_params_c15_test_hill_sim((SFc15_test_hill_simInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c15_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_hill_sim_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      15);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,15,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,15,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,15);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,15,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,15,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 2; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,15);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(24831336U));
  ssSetChecksum1(S,(980655408U));
  ssSetChecksum2(S,(2512327427U));
  ssSetChecksum3(S,(751968196U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c15_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c15_test_hill_sim(SimStruct *S)
{
  SFc15_test_hill_simInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc15_test_hill_simInstanceStruct *)utMalloc(sizeof
    (SFc15_test_hill_simInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc15_test_hill_simInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c15_test_hill_sim;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c15_test_hill_sim;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c15_test_hill_sim;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c15_test_hill_sim;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c15_test_hill_sim;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c15_test_hill_sim;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c15_test_hill_sim;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c15_test_hill_sim;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c15_test_hill_sim;
  chartInstance->chartInfo.mdlStart = mdlStart_c15_test_hill_sim;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c15_test_hill_sim;
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

void c15_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c15_test_hill_sim(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c15_test_hill_sim(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c15_test_hill_sim(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c15_test_hill_sim_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
