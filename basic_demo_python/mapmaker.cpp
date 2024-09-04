// Compile as:
// g++ -shared -O3 -fPIC -fopenmp mapmaker.cpp -o mapmaker.so
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

extern "C"
void mapmaker(double *map, double *map_inv_var, double *tod, int64_t *pix, double *scan_rms, int64_t scan_len, int64_t num_scans, int64_t num_pix){
    /** Simple serial mapmaker, with the two first arguments being empty signal and inv-var maps to be filled.
     * 
     *  Args:
     *      map (OUTPUT) -- 1D array of length 'num_pix', representing the signal map, which will be populated by this function.
     *      map_inv_var (OUTPUT) -- 1D array of length 'num_pix', reprensenting inverse map variance, which will be populated with this function.
     *      tod -- 2D array, containing the TOD of 'num_scan' different scans, each of length 'scan_len'.
     *      pix -- 2D array, containing the pixel pointing index of each element in tod.
     *      scan_rms -- 1D array of length 'num_scans', containing the sigma0 estimate of each scan.
     *      scan_len -- Length of each scan as int.
     *      num_scans -- Number of scans as int.
     *      num_pix -- Number of pixels in map and map_rms.
     */

    for(int64_t iscan=0; iscan<num_scans; iscan++){
        double scan_inv_var = 1.0/(scan_rms[iscan]*scan_rms[iscan]);
        for(int64_t itod=0; itod<scan_len; itod++){
            int64_t idx = iscan*scan_len + itod;
            map[pix[idx]] += tod[idx]*scan_inv_var;
            map_inv_var[pix[idx]] += scan_inv_var;
        }
    }
    // #pragma omp parallel for
    for(int ipix=0; ipix<num_pix; ipix++){
        map[ipix] /= map_inv_var[ipix];
    }
}


extern "C"
void mapmaker_OMP(double *map, double *map_inv_var, double *tod, int64_t *pix, double *scan_rms, int64_t scan_len, int64_t num_scans, int64_t num_pix){
    /** OpenMP parallelized version of the above mapmaker, which allocates full maps per thread, and divides up the TOD.
     *  Only gives a benefit over the serial mapmaker if each pixel is hit many times (>>10).
     * 
     *  Args:
     *      map (OUTPUT) -- 1D array of length 'num_pix', representing the signal map, which will be populated by this function.
     *      map_inv_var (OUTPUT) -- 1D array of length 'num_pix', reprensenting inverse map variance, which will be populated with this function.
     *      tod -- 2D array, containing the TOD of 'num_scan' different scans, each of length 'scan_len'.
     *      pix -- 2D array, containing the pixel pointing index of each element in tod.
     *      scan_rms -- 1D array of length 'num_scans', containing the sigma0 estimate of each scan.
     *      scan_len -- Length of each scan as int.
     *      num_scans -- Number of scans as int.
     *      num_pix -- Number of pixels in map and map_rms.
     */

    #pragma omp parallel
    {
        std::vector<double> map_private(num_pix);
        std::vector<double> map_inv_var_private(num_pix);
        #pragma omp for
        for(int64_t iscan=0; iscan<num_scans; iscan++){
            double scan_inv_var = 1.0/(scan_rms[iscan]*scan_rms[iscan]);
            for(int64_t itod=0; itod<scan_len; itod++){
                int64_t idx = iscan*scan_len + itod;
                map_private[pix[idx]] += tod[idx]*scan_inv_var;
                map_inv_var_private[pix[idx]] += scan_inv_var;
            }
        }
        #pragma omp critical
        {
            for(int64_t ipix=0; ipix<num_pix; ipix++){
                map[ipix] += map_private[ipix];
                map_inv_var[ipix] += map_inv_var_private[ipix];
            }
        }
    }
    for(int ipix=0; ipix<num_pix; ipix++){
        map[ipix] /= map_inv_var[ipix];
    }
}