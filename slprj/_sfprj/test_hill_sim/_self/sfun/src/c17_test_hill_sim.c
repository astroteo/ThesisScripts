/* Include files */

#include <stddef.h>
#include "blas.h"
#include "test_hill_sim_sfun.h"
#include "c17_test_hill_sim.h"
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
static const char * c17_debug_family_names[6] = { "j", "nargin", "nargout",
  "a_sat", "a_c", "a_c_SAT" };

/* Function Declarations */
static void initialize_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void initialize_params_c17_test_hill_sim
  (SFc17_test_hill_simInstanceStruct *chartInstance);
static void enable_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void disable_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void c17_update_debugger_state_c17_test_hill_sim
  (SFc17_test_hill_simInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c17_test_hill_sim
  (SFc17_test_hill_simInstanceStruct *chartInstance);
static void set_sim_state_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_st);
static void finalize_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void sf_gateway_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void mdl_start_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void initSimStructsc17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c17_machineNumber, uint32_T
  c17_chartNumber, uint32_T c17_instanceNumber);
static const mxArray *c17_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static void c17_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_b_a_c_SAT, const char_T *c17_identifier,
  real_T c17_y[3]);
static void c17_b_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId,
  real_T c17_y[3]);
static void c17_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static const mxArray *c17_b_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static const mxArray *c17_c_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static real_T c17_c_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void c17_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static const mxArray *c17_d_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData);
static int32_T c17_d_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void c17_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData);
static uint8_T c17_e_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_b_is_active_c17_test_hill_sim, const char_T
  *c17_identifier);
static uint8_T c17_f_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId);
static void init_dsm_address_info(SFc17_test_hill_simInstanceStruct
  *chartInstance);
static void init_simulink_io_address(SFc17_test_hill_simInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c17_sfEvent = CALL_EVENT;
  _sfTime_ = sf_get_time(chartInstance->S);
  chartInstance->c17_is_active_c17_test_hill_sim = 0U;
}

static void initialize_params_c17_test_hill_sim
  (SFc17_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static void enable_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void disable_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  _sfTime_ = sf_get_time(chartInstance->S);
}

static void c17_update_debugger_state_c17_test_hill_sim
  (SFc17_test_hill_simInstanceStruct *chartInstance)
{
  (void)chartInstance;
}

static const mxArray *get_sim_state_c17_test_hill_sim
  (SFc17_test_hill_simInstanceStruct *chartInstance)
{
  const mxArray *c17_st;
  const mxArray *c17_y = NULL;
  int32_T c17_i0;
  real_T c17_u[3];
  const mxArray *c17_b_y = NULL;
  uint8_T c17_hoistedGlobal;
  uint8_T c17_b_u;
  const mxArray *c17_c_y = NULL;
  c17_st = NULL;
  c17_st = NULL;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_createcellmatrix(2, 1), false);
  for (c17_i0 = 0; c17_i0 < 3; c17_i0++) {
    c17_u[c17_i0] = (*chartInstance->c17_a_c_SAT)[c17_i0];
  }

  c17_b_y = NULL;
  sf_mex_assign(&c17_b_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_setcell(c17_y, 0, c17_b_y);
  c17_hoistedGlobal = chartInstance->c17_is_active_c17_test_hill_sim;
  c17_b_u = c17_hoistedGlobal;
  c17_c_y = NULL;
  sf_mex_assign(&c17_c_y, sf_mex_create("y", &c17_b_u, 3, 0U, 0U, 0U, 0), false);
  sf_mex_setcell(c17_y, 1, c17_c_y);
  sf_mex_assign(&c17_st, c17_y, false);
  return c17_st;
}

static void set_sim_state_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_st)
{
  const mxArray *c17_u;
  real_T c17_dv0[3];
  int32_T c17_i1;
  chartInstance->c17_doneDoubleBufferReInit = true;
  c17_u = sf_mex_dup(c17_st);
  c17_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c17_u, 0)),
                       "a_c_SAT", c17_dv0);
  for (c17_i1 = 0; c17_i1 < 3; c17_i1++) {
    (*chartInstance->c17_a_c_SAT)[c17_i1] = c17_dv0[c17_i1];
  }

  chartInstance->c17_is_active_c17_test_hill_sim = c17_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c17_u, 1)),
     "is_active_c17_test_hill_sim");
  sf_mex_destroy(&c17_u);
  c17_update_debugger_state_c17_test_hill_sim(chartInstance);
  sf_mex_destroy(&c17_st);
}

static void finalize_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void sf_gateway_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  int32_T c17_i2;
  real_T c17_hoistedGlobal;
  real_T c17_b_a_sat;
  int32_T c17_i3;
  real_T c17_b_a_c[3];
  uint32_T c17_debug_family_var_map[6];
  real_T c17_j;
  real_T c17_nargin = 2.0;
  real_T c17_nargout = 1.0;
  real_T c17_b_a_c_SAT[3];
  int32_T c17_i4;
  int32_T c17_b_j;
  real_T c17_d0;
  real_T c17_d1;
  int32_T c17_i5;
  int32_T c17_i6;
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = sf_get_time(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 15U, chartInstance->c17_sfEvent);
  _SFD_DATA_RANGE_CHECK(*chartInstance->c17_a_sat, 0U);
  for (c17_i2 = 0; c17_i2 < 3; c17_i2++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c17_a_c)[c17_i2], 1U);
  }

  chartInstance->c17_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 15U, chartInstance->c17_sfEvent);
  c17_hoistedGlobal = *chartInstance->c17_a_sat;
  c17_b_a_sat = c17_hoistedGlobal;
  for (c17_i3 = 0; c17_i3 < 3; c17_i3++) {
    c17_b_a_c[c17_i3] = (*chartInstance->c17_a_c)[c17_i3];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 6U, 6U, c17_debug_family_names,
    c17_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_j, 0U, c17_c_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_nargin, 1U, c17_c_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c17_nargout, 2U, c17_c_sf_marshallOut,
    c17_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c17_b_a_sat, 3U, c17_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c17_b_a_c, 4U, c17_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c17_b_a_c_SAT, 5U, c17_sf_marshallOut,
    c17_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 3);
  for (c17_i4 = 0; c17_i4 < 3; c17_i4++) {
    c17_b_a_c_SAT[c17_i4] = 0.0;
  }

  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 4);
  c17_j = 1.0;
  c17_b_j = 0;
  while (c17_b_j < 3) {
    c17_j = 1.0 + (real_T)c17_b_j;
    CV_EML_FOR(0, 1, 0, 1);
    _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 6);
    c17_d0 = c17_b_a_c[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c", (int32_T)
      _SFD_INTEGER_CHECK("j", c17_j), 1, 3, 1, 0) - 1];
    if (CV_EML_IF(0, 1, 0, CV_RELATIONAL_EVAL(4U, 0U, 0, c17_d0, c17_b_a_sat, -1,
          4U, c17_d0 > c17_b_a_sat))) {
      _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 8);
      c17_b_a_c_SAT[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c_SAT", (int32_T)
        _SFD_INTEGER_CHECK("j", c17_j), 1, 3, 1, 0) - 1] = c17_b_a_sat;
    } else {
      _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 10);
      c17_d1 = c17_b_a_c[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c", (int32_T)
        _SFD_INTEGER_CHECK("j", c17_j), 1, 3, 1, 0) - 1];
      if (CV_EML_IF(0, 1, 1, CV_RELATIONAL_EVAL(4U, 0U, 1, c17_d1, -c17_b_a_sat,
            -1, 2U, c17_d1 < -c17_b_a_sat))) {
        _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 12);
        c17_b_a_c_SAT[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c_SAT", (int32_T)
          _SFD_INTEGER_CHECK("j", c17_j), 1, 3, 1, 0) - 1] = -c17_b_a_sat;
      } else {
        _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, 14);
        c17_b_a_c_SAT[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c_SAT", (int32_T)
          _SFD_INTEGER_CHECK("j", c17_j), 1, 3, 1, 0) - 1] =
          c17_b_a_c[_SFD_EML_ARRAY_BOUNDS_CHECK("a_c", (int32_T)
          _SFD_INTEGER_CHECK("j", c17_j), 1, 3, 1, 0) - 1];
      }
    }

    c17_b_j++;
    _SF_MEX_LISTEN_FOR_CTRL_C(chartInstance->S);
  }

  CV_EML_FOR(0, 1, 0, 0);
  _SFD_EML_CALL(0U, chartInstance->c17_sfEvent, -14);
  _SFD_SYMBOL_SCOPE_POP();
  for (c17_i5 = 0; c17_i5 < 3; c17_i5++) {
    (*chartInstance->c17_a_c_SAT)[c17_i5] = c17_b_a_c_SAT[c17_i5];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 15U, chartInstance->c17_sfEvent);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_test_hill_simMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
  for (c17_i6 = 0; c17_i6 < 3; c17_i6++) {
    _SFD_DATA_RANGE_CHECK((*chartInstance->c17_a_c_SAT)[c17_i6], 2U);
  }
}

static void mdl_start_c17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void initSimStructsc17_test_hill_sim(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_script_number_translation(uint32_T c17_machineNumber, uint32_T
  c17_chartNumber, uint32_T c17_instanceNumber)
{
  (void)c17_machineNumber;
  (void)c17_chartNumber;
  (void)c17_instanceNumber;
}

static const mxArray *c17_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i7;
  real_T c17_b_inData[3];
  int32_T c17_i8;
  real_T c17_u[3];
  const mxArray *c17_y = NULL;
  SFc17_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc17_test_hill_simInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  for (c17_i7 = 0; c17_i7 < 3; c17_i7++) {
    c17_b_inData[c17_i7] = (*(real_T (*)[3])c17_inData)[c17_i7];
  }

  for (c17_i8 = 0; c17_i8 < 3; c17_i8++) {
    c17_u[c17_i8] = c17_b_inData[c17_i8];
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 1, 3), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static void c17_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_b_a_c_SAT, const char_T *c17_identifier,
  real_T c17_y[3])
{
  emlrtMsgIdentifier c17_thisId;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_b_a_c_SAT), &c17_thisId,
    c17_y);
  sf_mex_destroy(&c17_b_a_c_SAT);
}

static void c17_b_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId,
  real_T c17_y[3])
{
  real_T c17_dv1[3];
  int32_T c17_i9;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), c17_dv1, 1, 0, 0U, 1, 0U, 1, 3);
  for (c17_i9 = 0; c17_i9 < 3; c17_i9++) {
    c17_y[c17_i9] = c17_dv1[c17_i9];
  }

  sf_mex_destroy(&c17_u);
}

static void c17_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_b_a_c_SAT;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y[3];
  int32_T c17_i10;
  SFc17_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc17_test_hill_simInstanceStruct *)chartInstanceVoid;
  c17_b_a_c_SAT = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_b_a_c_SAT), &c17_thisId,
    c17_y);
  sf_mex_destroy(&c17_b_a_c_SAT);
  for (c17_i10 = 0; c17_i10 < 3; c17_i10++) {
    (*(real_T (*)[3])c17_outData)[c17_i10] = c17_y[c17_i10];
  }

  sf_mex_destroy(&c17_mxArrayInData);
}

static const mxArray *c17_b_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_i11;
  real_T c17_b_inData[3];
  int32_T c17_i12;
  real_T c17_u[3];
  const mxArray *c17_y = NULL;
  SFc17_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc17_test_hill_simInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  for (c17_i11 = 0; c17_i11 < 3; c17_i11++) {
    c17_b_inData[c17_i11] = (*(real_T (*)[3])c17_inData)[c17_i11];
  }

  for (c17_i12 = 0; c17_i12 < 3; c17_i12++) {
    c17_u[c17_i12] = c17_b_inData[c17_i12];
  }

  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", c17_u, 0, 0U, 1U, 0U, 2, 3, 1), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static const mxArray *c17_c_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  real_T c17_u;
  const mxArray *c17_y = NULL;
  SFc17_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc17_test_hill_simInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_u = *(real_T *)c17_inData;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", &c17_u, 0, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static real_T c17_c_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  real_T c17_y;
  real_T c17_d2;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_d2, 1, 0, 0U, 0, 0U, 0);
  c17_y = c17_d2;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void c17_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_nargout;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  real_T c17_y;
  SFc17_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc17_test_hill_simInstanceStruct *)chartInstanceVoid;
  c17_nargout = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_nargout),
    &c17_thisId);
  sf_mex_destroy(&c17_nargout);
  *(real_T *)c17_outData = c17_y;
  sf_mex_destroy(&c17_mxArrayInData);
}

const mxArray *sf_c17_test_hill_sim_get_eml_resolved_functions_info(void)
{
  const mxArray *c17_nameCaptureInfo = NULL;
  c17_nameCaptureInfo = NULL;
  sf_mex_assign(&c17_nameCaptureInfo, sf_mex_create("nameCaptureInfo", NULL, 0,
    0U, 1U, 0U, 2, 0, 1), false);
  return c17_nameCaptureInfo;
}

static const mxArray *c17_d_sf_marshallOut(void *chartInstanceVoid, void
  *c17_inData)
{
  const mxArray *c17_mxArrayOutData = NULL;
  int32_T c17_u;
  const mxArray *c17_y = NULL;
  SFc17_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc17_test_hill_simInstanceStruct *)chartInstanceVoid;
  c17_mxArrayOutData = NULL;
  c17_u = *(int32_T *)c17_inData;
  c17_y = NULL;
  sf_mex_assign(&c17_y, sf_mex_create("y", &c17_u, 6, 0U, 0U, 0U, 0), false);
  sf_mex_assign(&c17_mxArrayOutData, c17_y, false);
  return c17_mxArrayOutData;
}

static int32_T c17_d_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  int32_T c17_y;
  int32_T c17_i13;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_i13, 1, 6, 0U, 0, 0U, 0);
  c17_y = c17_i13;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void c17_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c17_mxArrayInData, const char_T *c17_varName, void *c17_outData)
{
  const mxArray *c17_b_sfEvent;
  const char_T *c17_identifier;
  emlrtMsgIdentifier c17_thisId;
  int32_T c17_y;
  SFc17_test_hill_simInstanceStruct *chartInstance;
  chartInstance = (SFc17_test_hill_simInstanceStruct *)chartInstanceVoid;
  c17_b_sfEvent = sf_mex_dup(c17_mxArrayInData);
  c17_identifier = c17_varName;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c17_b_sfEvent),
    &c17_thisId);
  sf_mex_destroy(&c17_b_sfEvent);
  *(int32_T *)c17_outData = c17_y;
  sf_mex_destroy(&c17_mxArrayInData);
}

static uint8_T c17_e_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_b_is_active_c17_test_hill_sim, const char_T
  *c17_identifier)
{
  uint8_T c17_y;
  emlrtMsgIdentifier c17_thisId;
  c17_thisId.fIdentifier = c17_identifier;
  c17_thisId.fParent = NULL;
  c17_y = c17_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c17_b_is_active_c17_test_hill_sim), &c17_thisId);
  sf_mex_destroy(&c17_b_is_active_c17_test_hill_sim);
  return c17_y;
}

static uint8_T c17_f_emlrt_marshallIn(SFc17_test_hill_simInstanceStruct
  *chartInstance, const mxArray *c17_u, const emlrtMsgIdentifier *c17_parentId)
{
  uint8_T c17_y;
  uint8_T c17_u0;
  (void)chartInstance;
  sf_mex_import(c17_parentId, sf_mex_dup(c17_u), &c17_u0, 1, 3, 0U, 0, 0U, 0);
  c17_y = c17_u0;
  sf_mex_destroy(&c17_u);
  return c17_y;
}

static void init_dsm_address_info(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  (void)chartInstance;
}

static void init_simulink_io_address(SFc17_test_hill_simInstanceStruct
  *chartInstance)
{
  chartInstance->c17_a_sat = (real_T *)ssGetInputPortSignal_wrapper
    (chartInstance->S, 0);
  chartInstance->c17_a_c = (real_T (*)[3])ssGetInputPortSignal_wrapper
    (chartInstance->S, 1);
  chartInstance->c17_a_c_SAT = (real_T (*)[3])ssGetOutputPortSignal_wrapper
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

void sf_c17_test_hill_sim_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3493999906U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(865193842U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2008263417U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3747457003U);
}

mxArray* sf_c17_test_hill_sim_get_post_codegen_info(void);
mxArray *sf_c17_test_hill_sim_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals", "postCodegenInfo" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1, 1, sizeof
    (autoinheritanceFields)/sizeof(autoinheritanceFields[0]),
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("7JHkTr2Jfy800CjMIZ8MuG");
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
    mxArray* mxPostCodegenInfo = sf_c17_test_hill_sim_get_post_codegen_info();
    mxSetField(mxAutoinheritanceInfo,0,"postCodegenInfo",mxPostCodegenInfo);
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c17_test_hill_sim_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c17_test_hill_sim_jit_fallback_info(void)
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

mxArray *sf_c17_test_hill_sim_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

mxArray* sf_c17_test_hill_sim_get_post_codegen_info(void)
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

static const mxArray *sf_get_sim_state_info_c17_test_hill_sim(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[5],T\"a_c_SAT\",},{M[8],M[0],T\"is_active_c17_test_hill_sim\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c17_test_hill_sim_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc17_test_hill_simInstanceStruct *chartInstance;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
    chartInstance = (SFc17_test_hill_simInstanceStruct *)
      chartInfo->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _test_hill_simMachineNumber_,
           17,
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
          (MexFcnForType)c17_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 1;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c17_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c17_sf_marshallOut,(MexInFcnForType)
            c17_sf_marshallIn);
        }

        _SFD_SET_DATA_VALUE_PTR(0U, chartInstance->c17_a_sat);
        _SFD_SET_DATA_VALUE_PTR(1U, *chartInstance->c17_a_c);
        _SFD_SET_DATA_VALUE_PTR(2U, *chartInstance->c17_a_c_SAT);
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
  return "v5J5KGSenGCKU3rBrfLczB";
}

static void sf_opaque_initialize_c17_test_hill_sim(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc17_test_hill_simInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
    chartInstanceVar);
  initialize_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c17_test_hill_sim(void *chartInstanceVar)
{
  enable_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c17_test_hill_sim(void *chartInstanceVar)
{
  disable_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c17_test_hill_sim(void *chartInstanceVar)
{
  sf_gateway_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

static const mxArray* sf_opaque_get_sim_state_c17_test_hill_sim(SimStruct* S)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  return get_sim_state_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
    chartInfo->chartInstance);         /* raw sim ctx */
}

static void sf_opaque_set_sim_state_c17_test_hill_sim(SimStruct* S, const
  mxArray *st)
{
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
  ChartInfoStruct * chartInfo = (ChartInfoStruct *)(crtInfo->instanceInfo);
  set_sim_state_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
    chartInfo->chartInstance, st);
}

static void sf_opaque_terminate_c17_test_hill_sim(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc17_test_hill_simInstanceStruct*) chartInstanceVar)->S;
    ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)(ssGetUserData(S));
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_test_hill_sim_optimization_info();
    }

    finalize_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
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
  initSimStructsc17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c17_test_hill_sim(SimStruct *S)
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
    initialize_params_c17_test_hill_sim((SFc17_test_hill_simInstanceStruct*)
      (chartInfo->chartInstance));
  }
}

static void mdlSetWorkWidths_c17_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_test_hill_sim_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(sf_get_instance_specialization(),infoStruct,
      17);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(sf_get_instance_specialization(),
                infoStruct,17,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop
      (sf_get_instance_specialization(),infoStruct,17,
       "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(sf_get_instance_specialization(),infoStruct,17);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,17,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,17,1);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,17);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2896044748U));
  ssSetChecksum1(S,(292310822U));
  ssSetChecksum2(S,(300788411U));
  ssSetChecksum3(S,(76161944U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c17_test_hill_sim(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c17_test_hill_sim(SimStruct *S)
{
  SFc17_test_hill_simInstanceStruct *chartInstance;
  ChartRunTimeInfo * crtInfo = (ChartRunTimeInfo *)utMalloc(sizeof
    (ChartRunTimeInfo));
  chartInstance = (SFc17_test_hill_simInstanceStruct *)utMalloc(sizeof
    (SFc17_test_hill_simInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc17_test_hill_simInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c17_test_hill_sim;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c17_test_hill_sim;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c17_test_hill_sim;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c17_test_hill_sim;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c17_test_hill_sim;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c17_test_hill_sim;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c17_test_hill_sim;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c17_test_hill_sim;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c17_test_hill_sim;
  chartInstance->chartInfo.mdlStart = mdlStart_c17_test_hill_sim;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c17_test_hill_sim;
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

void c17_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c17_test_hill_sim(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c17_test_hill_sim(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c17_test_hill_sim(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c17_test_hill_sim_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
