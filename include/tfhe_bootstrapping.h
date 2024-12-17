#ifndef _TFHE_BOOTSTRAPPING_H_
#define _TFHE_BOOTSTRAPPING_H_

/* includes */
#include <chrono>
#include <iostream>
#include <cassert>
#include <tfhe/tfhe_io.h>
#include <tfhe/tfhe_garbage_collector.h>
#include <tfhe/tfhe.h>

/* prototypes */
void tfhe_bootstrap_woKS_naive_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod);
void tfhe_bootstrap_naive_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod);

void tfhe_bootstrap_woKS_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod);
void tfhe_bootstrap_Id(LweSample* result, const LweBootstrappingKey* bk, const LweSample* x, uint64_t plain_mod);

void test_bootstrapping_wParams(uint64_t plain_mod, bool print_info = 1);

#endif
