#include "tfhe_bootstrapping.h"

/* namespaces */
using namespace std;
using namespace std::chrono;

static void this_MuxRotate(TLweSample *result, const TLweSample *accum, const TGswSample *bki, const int32_t barai,
    const TGswParams *bk_params) {
  // ACC = BKi*[(X^barai-1)*ACC]+ACC
  // temp = (X^barai-1)*ACC
  tLweMulByXaiMinusOne(result, barai, accum, bk_params->tlwe_params);
  // temp *= BKi
  tGswExternMulToTLwe(result, bki, bk_params);
  // ACC += temp
  tLweAddTo(result, accum, bk_params->tlwe_params);
}

/**
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
static void
this_blindRotate(TLweSample *accum, const TGswSample *bk, const int32_t *bara, const int32_t n, const TGswParams *bk_params) {

  //TGswSample* temp = new_TGswSample(bk_params);
  TLweSample *temp = new_TLweSample(bk_params->tlwe_params);
  TLweSample *temp2 = temp;
  TLweSample *temp3 = accum;

  for (int32_t i = 0; i < n; i++) {
    //cout << i << ":" << n << endl;
    const int32_t barai = bara[i];
    if (barai == 0) continue; //indeed, this is an easy case!

    this_MuxRotate(temp2, temp3, bk + i, barai, bk_params);
    swap(temp2, temp3);

  }
  if (temp3 != accum) {
    tLweCopy(accum, temp3, bk_params->tlwe_params);
  }

  delete_TLweSample(temp);
  //delete_TGswSample(temp);
}

/**
 * result = LWE(v_p) where p=barb-sum(bara_i.s_i) mod 2N
 * @param result the output LWE sample
 * @param v a 2N-elt anticyclic function (represented by a TorusPolynomial)
 * @param bk An array of n TGSW samples where bk_i encodes s_i
 * @param barb A coefficients between 0 and 2N-1
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
static void this_blindRotateAndExtract(LweSample *result,
    const TorusPolynomial *v,
    const TGswSample *bk,
    const int32_t barb,
    const int32_t *bara,
    const int32_t n,
    const TGswParams *bk_params) {

  const TLweParams *accum_params = bk_params->tlwe_params;
  const LweParams *extract_params = &accum_params->extracted_lweparams;
  const int32_t N = accum_params->N;
  const int32_t _2N = 2 * N;

  TorusPolynomial *testvectbis= new_TorusPolynomial(N);
  TLweSample *acc = new_TLweSample(accum_params);

  if (barb != 0) torusPolynomialMulByXai(testvectbis, _2N - barb, v);
  else torusPolynomialCopy(testvectbis, v);
  tLweNoiselessTrivial(acc, testvectbis, accum_params);
  this_blindRotate(acc, bk, bara, n, bk_params);
  tLweExtractLweSample(result, acc, extract_params, accum_params);

  delete_TLweSample(acc);
  delete_TorusPolynomial(testvectbis);
}

//Identity half-domain naive approach: to see the negaciclicity problem
void tfhe_bootstrap_woKS_naive_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod)
{
  const TGswParams *bk_params = bk->bk_params;
  const TLweParams *accum_params = bk->accum_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t N = accum_params->N;
  const int32_t Nx2 = 2 * N;
  const int32_t n = in_params->n;

  TorusPolynomial* testvect = new_TorusPolynomial(N);
  int32_t *bara = new int32_t[N];
  Torus32 offset = modSwitchToTorus32(1,(2*plain_mod));
  int32_t barb = modSwitchFromTorus32(x->b+offset, Nx2);
  for (int32_t i = 0; i < n; i++) {
    bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
  }

  assert(plain_mod < N);
  uint32_t plain_mods2 = plain_mod/2;
  uint32_t step = N/plain_mods2;
  for (int32_t i = 0; i < N; i+=step) {
    Torus32 tmp = modSwitchToTorus32((uint32_t)i/step, plain_mod);
    for(int32_t j = 0; j < step; ++j) {
      testvect->coefsT[i+j] = tmp;
    }
  }
  this_blindRotateAndExtract(result, testvect, bk->bk, barb, bara, n, bk_params);

  delete[] bara;
  delete_TorusPolynomial(testvect);
}

void tfhe_bootstrap_naive_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod)
{
  LweSample* u = new_LweSample(&bk->accum_params->extracted_lweparams);
  tfhe_bootstrap_woKS_naive_Id(u, bk, x, plain_mod);
  lweKeySwitch(result, bk->ks, u);

  delete_LweSample(u);
}

//Identity half-domain
void tfhe_bootstrap_woKS_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod)
{
  const TGswParams *bk_params = bk->bk_params;
  const TLweParams *accum_params = bk->accum_params;
  const LweParams *in_params = bk->in_out_params;
  const int32_t N = accum_params->N;
  const int32_t Nx2 = 2 * N;
  const int32_t n = in_params->n;

  TorusPolynomial* testvect = new_TorusPolynomial(N);
  int32_t *bara = new int32_t[N];
  Torus32 offset = modSwitchToTorus32(1,(2*plain_mod));
  int32_t barb = modSwitchFromTorus32(x->b+offset, Nx2);
  for (int32_t i = 0; i < n; i++) {
    bara[i] = modSwitchFromTorus32(x->a[i], Nx2);
  }

  assert(plain_mod < N);
  uint32_t plain_mods2 = plain_mod/2;
  uint32_t step = N/plain_mods2;
  for (int32_t i = 0; i < N; i+=step) {
    Torus32 tmp = modSwitchToTorus32((uint32_t)i/step, plain_mod);
    for(int32_t j = 0; j < step; ++j) {
      testvect->coefsT[i+j] = tmp + offset; //offset makes the difference
    }
  }
  this_blindRotateAndExtract(result, testvect, bk->bk, barb, bara, n, bk_params);

  delete[] bara;
  delete_TorusPolynomial(testvect);
}

void tfhe_bootstrap_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod)
{
  LweSample* u = new_LweSample(&bk->accum_params->extracted_lweparams);
  tfhe_bootstrap_woKS_Id(u, bk, x, plain_mod);
  lweKeySwitch(result, bk->ks, u);

  delete_LweSample(u);
}

/* test */
void test_bootstrapping_wParams(uint64_t plain_mod, bool print_info)
{
  cout << "**********************************" << endl;
  cout << "**** test tfhe bootstrapping  ****" << endl;
  cout << "**********************************" << endl;
  const int32_t N = 1024;
  const int32_t k = 1;
  const int32_t n = 1024;
  const int32_t bk_l = 13;
  const int32_t bk_Bgbit = 2;
  const int32_t ks_basebit = 2;
  const int32_t ks_l = 13;
  const double ks_stdev = pow(10.,-8); //standard deviation
  const double bk_stdev = pow(10.,-8); //standard deviation
  const double max_stdev = 1; //max standard deviation for a 1/4 msg space

  /*** lwe params and key ***/
  LweParams* lwe_params = new_LweParams(n, ks_stdev, max_stdev);
  LweKey* lwe_key = new_LweKey(lwe_params);

  if (print_info)
    cout << "[tfhe] input tlwe params and key generation..." << endl;

  /*** tlwe params for bs key ***/
  TLweParams* accum_params = new_TLweParams(N, k, bk_stdev, max_stdev);

  if (print_info)
    cout << "[tfhe] accumulator trlwe params generation..." << endl;

  /*** bs params and key ***/
  TGswParams* bs_params = new_TGswParams(bk_l, bk_Bgbit, accum_params);
  TGswKey* tgsw_key = new_TGswKey(bs_params);
  LweBootstrappingKey* bs_key = new_LweBootstrappingKey(ks_l, ks_basebit, lwe_params, bs_params);
  tfhe_createLweBootstrappingKey(bs_key, lwe_key, tgsw_key);

  if (print_info)
    cout << "[tfhe] bootstrapping params and key generation..." << endl;

  /*** bootstrapping example ***/
  LweSample* ct = new_LweSample(lwe_params);
  LweSample* ct_bs = new_LweSample(lwe_params);
  LweSample* ct_bs_id0 = new_LweSample(lwe_params);
  LweSample* ct_bs_id1 = new_LweSample(lwe_params);
  Torus32 res = 0;

  for (uint32_t pt = 0; pt < plain_mod; ++pt) {
    Torus32 pt_torus = modSwitchToTorus32(pt, plain_mod);
    if (print_info)
      cout << "[tfhe] plaintext: " << pt << "/" << plain_mod << endl;

    lweSymEncrypt(ct, pt_torus, lwe_params->alpha_min, lwe_key);

    /* default bootstrapping tfhe */
    tfhe_bootstrap(ct_bs, bs_key, modSwitchToTorus32(1, plain_mod), ct);
    res = lweSymDecrypt(ct_bs, lwe_key, plain_mod);
    if (print_info) {
      cout << "[default bs] decrypted_in_torus: " << res
        << " ***** decrypted_in_int: " << modSwitchFromTorus32(res, plain_mod)
        << endl;
    }


    /* bootstrapping naive id */
    tfhe_bootstrap_naive_Id(ct_bs_id0, bs_key, ct, plain_mod);
    res = lweSymDecrypt(ct_bs_id0, lwe_key, plain_mod);
    if (print_info) {
      cout << "[naive id bs]  decrypted_in_torus: " << res
        << " ***** decrypted_in_int: " << modSwitchFromTorus32(res, plain_mod)
        << endl;
    }

    /* bootstrapping id */
    tfhe_bootstrap_Id(ct_bs_id1, bs_key, ct, plain_mod);
    Torus32 offset = modSwitchToTorus32(1,(2*plain_mod));
    ct_bs_id1->b -= offset;
    res = lweSymDecrypt(ct_bs_id1, lwe_key, plain_mod);
    if (print_info) {
      cout << "[id bs: step1] decrypted_in_torus: " << res
        << " ***** decrypted_in_int: " << modSwitchFromTorus32(res, plain_mod)
        << endl;
    }

    tfhe_bootstrap_Id(ct_bs_id1, bs_key, ct_bs_id1, plain_mod);
    ct_bs_id1->b -= offset;
    res = lweSymDecrypt(ct_bs_id1, lwe_key, plain_mod);
    if (print_info) {
      cout << "[id bs: step2] decrypted_in_torus: " << res
        << " ***** decrypted_in_int: " << modSwitchFromTorus32(res, plain_mod)
        << endl << endl;
    }
  }

  /*** free the memory ***/
  delete_LweSample(ct_bs_id0);
  delete_LweSample(ct_bs_id1);
  delete_LweSample(ct_bs);
  delete_LweSample(ct);

  delete_LweBootstrappingKey(bs_key);
  delete_TGswKey(tgsw_key);
  delete_TGswParams(bs_params);

  delete_TLweParams(accum_params);

  delete_LweKey(lwe_key);
  delete_LweParams(lwe_params);
}
