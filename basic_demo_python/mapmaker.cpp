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


extern "C"
void mapmaker_IQU(double *map, double *map_inv_var, double *tod, int64_t *pix, double *cos2psi, double *sin2psi, double *scan_rms, int64_t scan_len, int64_t num_scans, int64_t num_pix){
    /** Simple serial mapmaker, with the two first arguments being empty signal and inv-var maps to be filled.
     * 
     *  Args:
     *      map (OUTPUT) -- 2D array of shape '(3, num_pix)', representing the signal map, which will be populated by this function.
     *      map_inv_var (OUTPUT) -- 2D array of length '(6, num_pix)', reprensenting inverse variance matrix, which will be populated with this function.
     *      tod -- 2D array, containing the TOD of 'num_scan' different scans, each of length 'scan_len'.
     *      pix -- 2D array, containing the pixel pointing index of each element in tod.
     *      psi -- 2D array, containing the polarization angle of each element in tod.
     *      scan_rms -- 1D array of length 'num_scans', containing the sigma0 estimate of each scan.
     *      scan_len -- Length of each scan as int.
     *      num_scans -- Number of scans as int.
     *      num_pix -- Number of pixels in map and map_rms.
     */

    for(int64_t iscan=0; iscan<num_scans; iscan++){
        double scan_inv_var = 1.0/(scan_rms[iscan]*scan_rms[iscan]);
        for(int64_t itod=0; itod<scan_len; itod++){
            int64_t idx = iscan*scan_len + itod;
            map[pix[idx]]             += tod[idx]*scan_inv_var;
            map_inv_var[pix[idx]]             += scan_inv_var;
        }
        for(int64_t itod=0; itod<scan_len; itod++){
            int64_t idx = iscan*scan_len + itod;
	    map[pix[idx] +   num_pix] += tod[idx]*scan_inv_var*cos2psi[idx];
            map_inv_var[pix[idx] +   num_pix] += scan_inv_var*cos2psi[idx];
        }
        for(int64_t itod=0; itod<scan_len; itod++){
            int64_t idx = iscan*scan_len + itod;
	    map[pix[idx] +   2*num_pix] += tod[idx]*scan_inv_var*sin2psi[idx];
            map_inv_var[pix[idx] +   2*num_pix] += scan_inv_var*sin2psi[idx];
        }
        for(int64_t itod=0; itod<scan_len; itod++){
            int64_t idx = iscan*scan_len + itod;
            map_inv_var[pix[idx] +   3*num_pix] += scan_inv_var*pow(cos2psi[idx],2);
        }
        for(int64_t itod=0; itod<scan_len; itod++){
            int64_t idx = iscan*scan_len + itod;
            map_inv_var[pix[idx] +   4*num_pix] += scan_inv_var*cos2psi[idx]*sin2psi[idx];
        }
        for(int64_t itod=0; itod<scan_len; itod++){
            int64_t idx = iscan*scan_len + itod;
            map_inv_var[pix[idx] +   5*num_pix] += scan_inv_var*pow(sin2psi[idx],2);
        }
    }

    std::vector<double> A(num_pix);
    std::vector<double> B(num_pix);
    std::vector<double> C(num_pix);
    std::vector<double> E(num_pix);
    std::vector<double> F(num_pix);
    std::vector<double> I(num_pix);
    std::vector<double> det(num_pix);



    for(int ipix=0; ipix<num_pix; ipix++){
	int64_t i_I  = ipix;
	int64_t i_Q  = ipix +   num_pix;
	int64_t i_U  = ipix + 2*num_pix;
	int64_t i_IQ = ipix + 3*num_pix;
	int64_t i_IU = ipix + 4*num_pix;
	int64_t i_QU = ipix + 5*num_pix;
	A[ipix] = map_inv_var[i_IQ]*map_inv_var[i_QU] - pow(map_inv_var[i_U],2);
	B[ipix] = -(map_inv_var[i_Q]*map_inv_var[i_QU] - map_inv_var[i_U]*map_inv_var[i_IU]);
	C[ipix] = map_inv_var[i_U]*map_inv_var[i_IU] - map_inv_var[i_U]*map_inv_var[i_IQ];
	E[ipix] = map_inv_var[i_I]*map_inv_var[i_QU] - pow(map_inv_var[i_U],2);
	F[ipix] = -(map_inv_var[i_I]*map_inv_var[i_IU] - map_inv_var[i_Q]*map_inv_var[i_U]);
	I[ipix] = map_inv_var[i_I]*map_inv_var[i_IQ] - pow(map_inv_var[i_Q],2);
	det[ipix] = map_inv_var[i_I]*A[ipix] + map_inv_var[i_Q]*B[ipix] + map_inv_var[i_U]*C[ipix];
    }

    std::vector<double> T(num_pix);
    std::vector<double> Q(num_pix);
    std::vector<double> U(num_pix);

    for(int ipix=0; ipix<num_pix; ipix++){
	    T[ipix] = (A[ipix]*map[ipix] + B[ipix]*map[ipix + num_pix] + C[ipix]*map[ipix + 2*num_pix])/det[ipix];
	    Q[ipix] = (B[ipix]*map[ipix] + E[ipix]*map[ipix + num_pix] + F[ipix]*map[ipix + 2*num_pix])/det[ipix];
	    U[ipix] = (C[ipix]*map[ipix] + F[ipix]*map[ipix + num_pix] + I[ipix]*map[ipix + 2*num_pix])/det[ipix];
    }

    for(int ipix=0; ipix<num_pix; ipix++){
	    map[ipix] = T[ipix];
    }
    for(int ipix=0; ipix<num_pix; ipix++){
	    map[ipix + num_pix] = Q[ipix];
    }
    for(int ipix=0; ipix<num_pix; ipix++){
	    map[ipix + 2*num_pix] = U[ipix];
    }


    for(int ipix=0; ipix<num_pix; ipix++){
	    map_inv_var[ipix] = A[ipix]/det[ipix];
    }
    for(int ipix=0; ipix<num_pix; ipix++){
	    map_inv_var[ipix + num_pix] = E[ipix]/det[ipix];
    }
    for(int ipix=0; ipix<num_pix; ipix++){
	    map_inv_var[ipix + 2*num_pix] = I[ipix]/det[ipix];
    }
    for(int ipix=0; ipix<num_pix; ipix++){
	    map_inv_var[ipix + 3*num_pix] = B[ipix]/det[ipix];
    }
    for(int ipix=0; ipix<num_pix; ipix++){
	    map_inv_var[ipix + 4*num_pix] = C[ipix]/det[ipix];
    }
    for(int ipix=0; ipix<num_pix; ipix++){
	    map_inv_var[ipix + 5*num_pix] = F[ipix]/det[ipix];
    }

}
