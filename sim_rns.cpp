
#include "sim_rns.h"
#include "nocopy_integer.h"

template <class T_I, size_t N>
class RNS {
    static RNS<T_I, NoCopyInteger, N>::ReducedIntVector RNS<T_I, NoCopyInteger, N>::naive_reduce(const vector<const NoCopyInteger*>& inputs, const ConstNumPtrArray<NoCopyInteger, N>& modulis)
    {
        size_t len_inputs = inputs.size();
        ReducedIntVector vec(len_inputs);
        for (size_t input_i = 0; input_i < len_inputs; input_i++) {
            const NoCopyInteger* ith_input = inputs[input_i];
            const PtrAllocator<NoCopyInteger>& residuals_alloc = [&](size_t res_index) {
                NoCopyInteger* p = new NoCopyInteger;
                p->EXPENSIVE_COPY(*ith_input);
                const NoCopyInteger& m = modulis[res_index];
                p->operator%=(m);
                return p;
            };
            vec[input_i] = new ReducedInt(modulis, residuals_alloc);
        }
        return vec;
    }
};
