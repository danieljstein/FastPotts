#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
bool openmp_enabled_cpp() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

// [[Rcpp::export]]
int omp_max_threads_cpp() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}
