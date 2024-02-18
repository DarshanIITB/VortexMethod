#include <vector>
#include <cmath>
#include <corecrt_math_defines.h>
#include <numeric>
#include <tuple>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "Utils.hpp"
#include "json.hpp"
//#include "./coordinates/Coordinates.h"

using std::vector, std::pow, std::sqrt, std::cout;

const float rho = 1.225f;

const int span_sections = 10; // number of points used per blade
const int N = 24; // number if turns in one rotation
const int NR = 7; // number of rotations
const int NR_F = 2;
const int FAR_FIELD_CUTOFF = (NR - NR_F) * N;
const float dpsi_d = (360 / N);
const float mu = 0.0f;

const float alpha_d = 0;
const float alpha = alpha_d * (M_PI / 180);
const float alpha_p = -alpha;

const float RPM = 383.0f; // rpm at cruise
const float OMEGA = 2 * M_PI * RPM / 60;

const float a = 5.75; // lift curve slope per radian from literature
const float Cd0 = 0.022;  // co - efficient of drag from literature

const float Vv = 0.0f; // tip path plane in radians

const float R_tip = 5.5f;
const float u_inf = mu * (OMEGA * R_tip) / cos(alpha); // this is Y direction
const float R_hub = 0.1 * R_tip;
vector<float> radii = linspace(R_hub, R_tip, span_sections);
const float chord_at_hub = 0.325;
const float taper_ratio = 1.0f;
const float chord_at_tip = chord_at_hub * taper_ratio;
const float taper_slope = (chord_at_tip - chord_at_hub) / (R_tip - R_hub);
const float chord_at_root = chord_at_hub - taper_slope * R_hub;
vector<float> chord = scalarAdd(chord_at_root, scalarMul(taper_slope, radii));

const int N_BLADES = 2;
const int N_ROTORS = 1;
const int BLADES_PER_ROTOR = N_BLADES / N_ROTORS;

const float x_displacement = 1.5 * R_tip;
const float y_displacement = 0 * R_tip;
const float z_displacement = 0.2 * R_tip;

const float d = sqrt(pow(x_displacement, 2) + pow(y_displacement, 2));

const float OVERLAP_AREA = 2 * pow(R_tip, 2) * std::acos(d / (2 * R_tip)) - 0.5 * sqrt(4 * pow(R_tip, 2) - pow(d, 2));
const float DISC_AREA = (N_ROTORS > 1) ? M_PI * pow(R_tip, 2) * N_ROTORS - OVERLAP_AREA : M_PI * pow(R_tip, 2) * N_ROTORS;

const float dpsi = dpsi_d * (M_PI / 180);
const float NT = (NR * (360)) / dpsi_d;
const int JM = int(NT);

const float alp = 1.2;
const float K_visc = 0.000015;
const float del_s = 1000;

const float THETA_TIP = 2.0f;
const float THETA_HUB = 12.0f;
const float F = (THETA_TIP - THETA_HUB) / (R_tip - R_hub);
const float E = THETA_HUB - F * R_hub;
const float lc = Vv / (OMEGA * R_tip);

const float blade_twist = (THETA_TIP - THETA_HUB) * (M_PI / 180);

const float h_by_r = 0.0f; // for ground effect

vector<float> theta = scalarAdd(E * M_PI / 180, (scalarMul(F * M_PI / 180, radii)));
vector<float> solidity = scalarMul(N_BLADES/(M_PI * R_tip), chord);

vector<vector<float>> bsl_2(vector<vector<vector<vector<float>>>>& positions, vector<vector<vector<float>>>& trail_circulation, int blade_index, int span_index, vector <vector< float >> &v_ind)
{
    for (int j = 0; j < FAR_FIELD_CUTOFF; j++)
    {
        const int xc = positions[j][span_index][0][blade_index];
        const int yc = positions[j][span_index][1][blade_index];
        const int zc = positions[j][span_index][2][blade_index];
        for (int blade_index = 0; blade_index < N_BLADES; blade_index++)
        {
            for (int i = 0; i < span_sections; i++)
            {
                for (int k = 0; k < JM - 1; k++)
                {
                    vector<float> v_ind_c = vector<float>(3, 0.0f);
                    if (k != j)
                    {
                        const float d_gamma = trail_circulation[k][i][blade_index] - trail_circulation[k+1][i][blade_index];

                        const float xa = positions[k][i][0][blade_index];
                        const float xb = positions[k+1][i][0][blade_index];
                        const float ya = positions[k][i][1][blade_index];
                        const float yb = positions[k+1][i][1][blade_index];
                        const float za = positions[k][i][2][blade_index];
                        const float zb = positions[k+1][i][2][blade_index];

                        // calculating differences for optimization
                        const float xb_minus_xa = xb - xa;
                        const float yb_minus_ya = yb - ya;
                        const float zb_minus_za = zb - za;
                        
                        const float xc_minus_xb = xc - xb;
                        const float yc_minus_yb = yc - yb;
                        const float zc_minus_zb = zc - zb;

                        // distance b/w A and B
                        const float L = sqrt(pow(xb_minus_xa, 2) + pow(yb_minus_ya, 2) + pow(zb_minus_za, 2));

                        const float cutoff = sqrt(pow(0.25f * L, 4) + 4 * alp * del_s * K_visc * ((k + 1) * dpsi / OMEGA));
                        // distance b/w A and C
                        const float r1 = sqrt(pow(xc - xa, 2) + pow(yc - ya, 2) + pow(zc - za, 2)) + cutoff;

                        // distance b/w B and C
                        const float r2 = sqrt(pow(xc_minus_xb, 2) + pow(yc_minus_yb, 2) + pow(zc_minus_zb, 2)) + cutoff;

                        const float multiplier = d_gamma * (r1 + r2) / (2 * M_PI * r1 * r2 * (pow(r1 + r2, 2) - pow(L, 2)));

                        float u_bs = multiplier * (zc_minus_zb * yb_minus_ya - yc_minus_yb * zc_minus_zb);
                        float v_bs = multiplier * (xc_minus_xb * zb_minus_za - zc_minus_zb * xb_minus_xa);
                        float w_bs = multiplier * (yc_minus_yb * xb_minus_xa - xc_minus_xb * yb_minus_ya);

                        if (blade_index + 1 > BLADES_PER_ROTOR)
                        {
                            u_bs *= -1;
                            v_bs *= -1;
                            w_bs *= -1;
                        }
                        v_ind_c = { u_bs, v_bs, w_bs };
                    }
                    for (int t = 0; t < 3; t++)
                    {
                        v_ind[j][t] += v_ind_c[t];
                    }
                }
            }
        }
    }
    for (int rot = 0; rot < NR_F; rot++)
    {
        for (int i = 0; i <= N; i++)
        {
            v_ind[FAR_FIELD_CUTOFF + N * rot + i] = v_ind[FAR_FIELD_CUTOFF - N + i];
        }
    }
    return v_ind;
}

vector<float> bsl_perf(vector<vector<vector<vector<float>>>>& positions, vector<vector<vector<float>>>& trail_circulation, int blade_index, int span_index, vector<float>& v_ind)
// not using shed circulation for now
{
    const int xc = positions[0][span_index][0][blade_index];
    const int yc = positions[0][span_index][1][blade_index];
    const int zc = positions[0][span_index][2][blade_index];
    for (int blade_index = 0; blade_index < N_BLADES; blade_index++)
    {
        for (int i = 0; i < span_sections; i++)
        {
            for (int k = 0; k < JM - 1; k++)
            {
                // d_gamma should be the difference in the trail circulation
                const float d_gamma = trail_circulation[k][i][blade_index] - trail_circulation[k + 1][i][blade_index];

                const float xa = positions[k][i][0][blade_index];
                const float xb = positions[k + 1][i][0][blade_index];
                const float ya = positions[k][i][1][blade_index];
                const float yb = positions[k + 1][i][1][blade_index];
                const float za = positions[k][i][2][blade_index];
                const float zb = positions[k + 1][i][2][blade_index];

                vector<float> v_ind_c = vector<float>(3, 0.0f);

                // calculating differences for optimization
                const float xb_minus_xa = xb - xa;
                const float yb_minus_ya = yb - ya;
                const float zb_minus_za = zb - za;

                const float xc_minus_xb = xc - xb;
                const float yc_minus_yb = yc - yb;
                const float zc_minus_zb = zc - zb;

                // distance b/w A and B
                const float L = sqrt(pow(xb_minus_xa, 2) + pow(yb_minus_ya, 2) + pow(zb_minus_za, 2));

                const float cutoff = sqrt(pow(0.25f * L, 4) + 4 * alp * del_s * K_visc * ((k + 1) * dpsi / OMEGA));
                // distance b/w A and C
                const float r1 = sqrt(pow(xc - xa, 2) + pow(yc - ya, 2) + pow(zc - za, 2)) + cutoff;

                // distance b/w B and C
                const float r2 = sqrt(pow(xc_minus_xb, 2) + pow(yc_minus_yb, 2) + pow(zc_minus_zb, 2)) + cutoff;

                const float multiplier = d_gamma * (r1 + r2) / (2 * M_PI * r1 * r2 * (pow(r1 + r2, 2) - pow(L, 2)));

                float u_bs = multiplier * (zc_minus_zb * yb_minus_ya - yc_minus_yb * zc_minus_zb);
                float v_bs = multiplier * (xc_minus_xb * zb_minus_za - zc_minus_zb * xb_minus_xa);
                float w_bs = multiplier * (yc_minus_yb * xb_minus_xa - xc_minus_xb * yb_minus_ya);

                if (blade_index + 1 > BLADES_PER_ROTOR)
                {
                    u_bs *= -1;
                    v_bs *= -1;
                    w_bs *= -1;
                }
                v_ind_c = { u_bs, v_bs, w_bs };
                for (int t = 0; t < 3; t++)
                {
                    v_ind[t] += v_ind_c[t];
                }
            }
        }
    }
    return v_ind;
}

//void save4DArray(const vector<vector<vector<vector<float>>>>& array, const char* filename) {
//    std::ofstream file(filename, std::ios::binary);
//    if (file.is_open()) {
//        // Write dimensions
//        size_t dim1 = array.size();
//        size_t dim2 = array[0].size();
//        size_t dim3 = array[0][0].size();
//        size_t dim4 = array[0][0][0].size();
//        file.write(reinterpret_cast<const char*>(&dim1), sizeof(size_t));
//        file.write(reinterpret_cast<const char*>(&dim2), sizeof(size_t));
//        file.write(reinterpret_cast<const char*>(&dim3), sizeof(size_t));
//        file.write(reinterpret_cast<const char*>(&dim4), sizeof(size_t));
//
//        // Write data
//        for (const auto& dim1_slice : array) {
//            for (const auto& dim2_slice : dim1_slice) {
//                for (const auto& dim3_slice : dim2_slice) {
//                    for (const auto& row : dim3_slice) {
//                        file.write(reinterpret_cast<const char*>(&row), dim4 * sizeof(float));
//                    }
//                }
//            }
//        }
//        file.close();
//        std::cout << "4D array saved successfully.\n";
//    }
//    else {
//        std::cerr << "Unable to open file.\n";
//    }
//}

//void save4DArray(const std::vector<std::vector<std::vector<std::vector<float>>>>& array, const char* filename) {
//    std::ofstream file(filename, std::ios::binary);
//    if (file.is_open()) {
//        // Write magic number for npy format
//        const char magic[] = "\x93NUMPY";
//        file.write(magic, sizeof(magic) - 1); // -1 to exclude null terminator
//
//        // Write version number
//        const char version[] = "\x01\x00";
//        file.write(version, sizeof(version) - 1); // -1 to exclude null terminator
//
//        // Write header size (int32)
//        int32_t header_size = 0;
//        file.write(reinterpret_cast<const char*>(&header_size), sizeof(int32_t));
//
//        // Write header content
//        std::string header = "{\"descr\": \"<f4\", \"fortran_order\": False, \"shape\": (";
//        header += std::to_string(array.size()) + ", " + std::to_string(array[0].size()) + ", " + std::to_string(array[0][0].size()) + ", " + std::to_string(array[0][0][0].size());
//        header += "), }\n";
//        file << header;
//
//        // Write data
//        for (const auto& dim1_slice : array) {
//            for (const auto& dim2_slice : dim1_slice) {
//                for (const auto& dim3_slice : dim2_slice) {
//                    for (const auto& row : dim3_slice) {
//                        file.write(reinterpret_cast<const char*>(&row), sizeof(float));
//                    }
//                }
//            }
//        }
//        file.close();
//        std::cout << "4D array saved successfully.\n";
//    }
//    else {
//        std::cerr << "Unable to open file.\n";
//    }
//}

void saveData(vector<vector<vector<vector<float>>>>& array4D)
{
    std::ofstream outfile("array4D.bin", std::ios::binary);

    int dim1 = array4D.size();
    int dim2 = array4D.empty() ? 0 : array4D[0].size();
    int dim3 = array4D.empty() || array4D[0].empty() ? 0 : array4D[0][0].size();
    int dim4 = array4D.empty() || array4D[0].empty() || array4D[0][0].empty() ? 0 : array4D[0][0][0].size();

    outfile.write(reinterpret_cast<const char*>(&dim1), sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&dim2), sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&dim3), sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&dim4), sizeof(int));

    for (const auto& threeD : array4D) {
        for (const auto& twoD : threeD) {
            for (const auto& oneD : twoD) {
                for (const auto& element : oneD) {
                    outfile.write(reinterpret_cast<const char*>(&element), sizeof(float));
                }
            }
        }
    }
    outfile.close();
}

vector< vector<vector<vector<float>>>> HoverWake(float CT_input)
{
    const float CT = CT_input / N_ROTORS;
    const float A = 0.78;
    const float lam = 0.145 + 27 * CT;
    vector<vector<float>> Zv_array(JM, vector<float>(span_sections, 0));
    vector<vector<vector<float>>> arrayX1(JM, vector<vector<float>>(span_sections, vector<float>(N_BLADES, 0)));
    vector<vector<vector<float>>> arrayY1(JM, vector<vector<float>>(span_sections, vector<float>(N_BLADES, 0)));
    vector<vector<vector<float>>> arrayZ1(JM, vector<vector<float>>(span_sections, vector<float>(N_BLADES, 0)));
    for (int i = 0; i < span_sections; i++)
    {
        if (i < span_sections - 1)
        {
            for (int j = 0; j < JM; j++) // time loop
            {
                const float psi = j * dpsi;
                const float R_vortex = radii[i] * (A + (1 - A) * std::exp(-lam * psi));
                const float chord_ = chord[i];
                const float solidity_ = (BLADES_PER_ROTOR * chord_) / (M_PI * R_tip);
                const float K1 = -2.2 * sqrt(0.5 * CT);
                const float K2 = -2.7 * sqrt(0.5 * CT);
                if (psi <= (2 * M_PI / BLADES_PER_ROTOR))
                {
                    Zv_array[j][i] = radii[i] * (K1 * psi);
                }
                else
                {
                    Zv_array[j][i] = radii[i] * (K1 * (2 * M_PI / BLADES_PER_ROTOR) + K2 * (psi - (2 * M_PI / BLADES_PER_ROTOR)));
                }
                for (int bn = 0; bn < N_BLADES; bn++)
                {
                    if (bn < BLADES_PER_ROTOR)
                    {
                        std::pair<float, float> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
                        arrayX1[j][i][bn] = new_coords.first;
                        arrayY1[j][i][bn] = new_coords.second;
                        for (int a = 0; a < JM; a++)
                        {
                            for (int b = 0; b < span_sections; b++)
                            {
                                arrayZ1[a][b][bn] = Zv_array[a][b] + h_by_r * R_tip;
                            }
                        }
                    }
                    else
                    {
                        std::pair<float, float> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
                        arrayX1[j][i][bn] = new_coords.first;
                        arrayY1[j][i][bn] = new_coords.second;
                        if (h_by_r == 0)
                        {
                            for (int a = 0; a < JM; a++)
                            {
                                for (int b = 0; b < span_sections; b++)
                                {
                                    arrayZ1[a][b][bn] = Zv_array[a][b] + z_displacement;
                                }
                            }
                        }
                        else
                        {
                            for (int a = 0; a < JM; a++)
                            {
                                for (int b = 0; b < span_sections; b++)
                                {
                                    arrayZ1[a][b][bn] = scalarMul(-1, Zv_array)[a][b] - h_by_r * R_tip;
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (i == span_sections-1)
        {
            for (int j = 0; j < JM; j++) // time loop
            {
                const float psi = j * dpsi;
                const float R_vortex = radii[i] * (A + (1 - A) * std::exp(-lam * psi));
                const float chord_ = chord[i];
                const float solidity_ = (BLADES_PER_ROTOR * chord_) / (M_PI * R_tip);
                const float K1 = -0.25 * ((CT / solidity_) + 0.001 * blade_twist * (180 / M_PI));
                const float K2 = -(1.41 + 0.0141 * blade_twist * (180 / M_PI)) * sqrt(0.5 * CT);
                if (psi <= (2 * M_PI / BLADES_PER_ROTOR))
                {
                    Zv_array[j][i] = radii[i] * (K1 * psi);
                }
                else
                {
                    Zv_array[j][i] = radii[i] * (K1 * (2 * M_PI / BLADES_PER_ROTOR) + K2 * (psi - (2 * M_PI / BLADES_PER_ROTOR)));
                }
                for (int bn = 0; bn < N_BLADES; bn++)
                {
                    if (bn < BLADES_PER_ROTOR)
                    {
                        std::pair<float, float> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
                        arrayX1[j][i][bn] = new_coords.first;
                        arrayY1[j][i][bn] = new_coords.second;
                        for (int a = 0; a < JM; a++)
                        {
                            for (int b = 0; b < span_sections; b++)
                            {
                                arrayZ1[a][b][bn] = Zv_array[a][b] + h_by_r * R_tip;
                            }
                        }
                    }
                    else
                    {
                        std::pair<float, float> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
                        arrayX1[j][i][bn] = new_coords.first;
                        arrayY1[j][i][bn] = new_coords.second;
                        if (h_by_r == 0)
                        {
                            for (int a = 0; a < JM; a++)
                            {
                                for (int b = 0; b < span_sections; b++)
                                {
                                    arrayZ1[a][b][bn] = Zv_array[a][b] + z_displacement;
                                }
                            }
                        }
                        else
                        {
                            for (int a = 0; a < JM; a++)
                            {
                                for (int b = 0; b < span_sections; b++)
                                {
                                    arrayZ1[a][b][bn] = scalarMul(-1, Zv_array)[a][b] - h_by_r * R_tip;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return { arrayX1, arrayY1, arrayZ1 };
}

int main()
{
    vector<float> psi_array = linspace(0, 2 * M_PI, N);
    const float location = R_hub;
    vector<float> PTC_inflow_arr(span_sections);
    vector<float> location_arr;

    for (int i = 0; i < span_sections; i++)
    {
        float diff = 1.0f;
        float lmd_1 = 1.0f;
        while (diff > 0.001)
        {
            float f = 0.5 * N_BLADES * (1 - radii[i] / R_tip) / lmd_1;
            float FL;
            if (f < 0.00001)
            {
                FL = (2 / M_PI) * std::acos(std::exp(-0.00001));
            }
            else
            {
                FL = (2 / M_PI) * std::acos(std::exp(-f));
            }
            double lmd_2 = std::sqrt(pow((a * solidity[i] / (16 * FL)) - (0.5 * lc), 2) +
                ((solidity[i] * a * theta[i] * radii[i]) / (FL * 8 * R_tip))) -
                ((a * solidity[i] / (16 * FL)) - 0.5 * lc);
            //cout << "lmd_2 term 1 part 1: " << pow((a * solidity[i] / (16 * FL)) - 0.5 * lc, 2)
            //    << ", lmd_2 term 1 part 2: " << ((solidity[i] * a * theta[i] * radii[i]) / (FL * 8 * R_tip))
            //    << ", lmd_2 term 1 original: " << pow((a * solidity[i] / (16 * FL)) - (0.5 * lc), 2) +
            //    ((solidity[i] * a * theta[i] * radii[i]) / (FL * 8 * R_tip))
            //    << ", lmd_2 term 1 root: " << std::sqrt(pow((a * solidity[i] / (16 * FL)) - (0.5 * lc), 2) +
            //        ((solidity[i] * a * theta[i] * radii[i]) / (FL * 8 * R_tip)))
            //    << std:: endl;
            //cout << "solidity[i]: " << solidity[i] << ", lc: " << lc << ", lmd_2: " << lmd_2 << std::endl;
            diff = abs(lmd_1 - lmd_2);
            lmd_1 = lmd_2;
        }
        PTC_inflow_arr[i] = lmd_1;
    }
    vector<float> phi_array = scalarMul(R_tip, vectorDiv(PTC_inflow_arr, radii));
    std::transform(phi_array.begin(), phi_array.end(), phi_array.begin(), arctan);
    vector<float> alpha_arr = vectorSub(theta, phi_array);
    vector<float> Cd_arr(alpha_arr);
    std::transform(Cd_arr.begin(), Cd_arr.end(), Cd_arr.begin(), square);
    Cd_arr = scalarAdd(Cd0, scalarMul(1.25, Cd_arr));
    vector<float> Cl_arr = scalarMul(a, alpha_arr);
    vector<float> circulation_ = scalarMul(0.5f, elementWiseMul(scalarAdd(u_inf * std::cos(alpha) * std::sin(M_PI / 2), scalarMul(OMEGA, radii)), elementWiseMul(chord, Cl_arr)));
    vector<float> circulation_old = circulation_;
    
    const float gamma_max = *std::max_element(circulation_.begin(), circulation_.end());
    vector<float> Thrust2_arr = vector<float>(span_sections, 0.0f);
    vector<float> Power2_arr = vector<float>(span_sections, 0.0f);
    vector<float> UP_arr = vector<float>(span_sections, 0.0f);
    const float T_a = 0.5 * rho * BLADES_PER_ROTOR; // constant part in thrust equation
    vector<float> Ut_arr = scalarMul(OMEGA, radii);
    for (int i = 0; i < span_sections; i++)
    {
        const float Ut = OMEGA * radii[i];  // tangential velocity
        const float Up = PTC_inflow_arr[i] * OMEGA * R_tip;  // perpendicular velocity
        const float Cd = Cd0 + 1.25 * pow(theta[i] - atan(Up / Ut), 2);  // drag equation in terms of r
        const float Cl = a * (theta[i] - atan(Up / Ut));  // lift equation in terms of r

        const float dT = 0.5 * rho * N_BLADES * (pow(Up, 2) + pow(Ut, 2)) * chord[i] * (Cl * cos(atan(Up / Ut)) - Cd * sin(atan(Up / Ut)));
        const float dP = 0.5 * rho * N_BLADES * OMEGA * radii[i] * (pow(Up, 2) + pow(Ut, 2)) * chord[i] * (Cd * cos(atan(Up / Ut)) + Cl * sin(atan(Up / Ut)));
        Thrust2_arr[i] = dT;
        Power2_arr[i] = dP;
        UP_arr[i] = Up;
    }
    const float Thrust = trapz(Thrust2_arr, radii) / cos(alpha);
    const float Power = trapz(Power2_arr, radii) / cos(alpha);
    const float CT_A1 = Thrust / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));
    float CP = Power / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));
    vector<float> circulation_A1 = scalarMul(0.5, elementWiseMul(Cl_arr, elementWiseMul(chord, Ut_arr)));

    const float LHS = CT_A1 / N_ROTORS;
    float il = 0.01f; // initial guess for Newton Raphson method

    for (int i = 0; i < 1000; i++)
    {
        const float h = - h1_ig(il, LHS, mu, alpha) / h2_ig(il, LHS, mu, alpha);
        il += h;
    }
    const float L_ig_calculated = il;
    cout << "L_ig_calculated: " << L_ig_calculated << std::endl;
    //const float LG = L_ig_calculated + L_i;
    vector<float> dOPM_array(N, 0.0f);
    vector<float> dORM_array(N, 0.0f);
    vector<float> dT_blade_array(N, 0.0f);
    vector<float> dOQi_array(N, 0.0f);
    vector<float> dOQo_array(N, 0.0f);
    vector<float> dT_array(N, 0.0f);
    vector<float> dOPi_array(N, 0.0f);
    vector<float> dOPo_array(N, 0.0f);

    const float dbeta_dt = 0;  // TPP as a reference frame
    const float beta_o_d = 0;  // TPP frame of reference # beta   in degrees
    const float beta_1s_d = 0;  // beta_1s kept 0 or as pilot input in degrees
    const float beta_1c_d = 0;  // beta_1c kept 0 pilot input  degrees
    const float theta_1s = (0) * (M_PI / 180);  // theta_1s  in degrees at trim conditio
    const float theta_1c = (0) * (M_PI / 180);  // theta_1c  in degrees at trim conditio
    vector<vector<float>> Thrust_plot_array(N + 1, vector<float>(span_sections, 0.0f));
    vector<vector<float>> circulation(JM, vector<float>(span_sections, 0.0f));
    vector<vector<float>> Up_BEM_array(N, vector<float>(span_sections, 0.0f));

    for (int j = 0; j < N; j++)
    {
        const float beta = (beta_o_d + beta_1c_d * cos(j * dpsi) + beta_1s_d * sin(j * dpsi)) * (M_PI / 180);  // beta
        // arrays to store sectional Lift, induced & profile Torque, induced & profile Power at each radial location
        vector<float> dL_array(span_sections, 0.0f);
        vector<float> dP_array(span_sections, 0.0f);
        vector<float> dR_array(span_sections, 0.0f);
        vector<float> dQi_array(span_sections, 0.0f);
        vector<float> dQo_array(span_sections, 0.0f);
        vector<float> dPi_array(span_sections, 0.0f);
        vector<float> dPo_array(span_sections, 0.0f);

        for (int i = 0; i < span_sections; i++) // loop to solve the dL, dQi, dQo, dPi, dPo, dPM, dRM
        {
            float inr = 0.0f;
            float difference = 1.0f;
            while (difference > 0.001f)
            {
                const float h = -h1_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip) / h2_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip);
                //cout << h1_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip) << " " << h2_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip) << ", ";
                difference = abs(h);
                inr = inr + h;
            }
            //cout << std::endl;
            const float L_i_calculated = inr; // calculated value of L_i (inflow ratio) at given radial location
            
            const float theta_equation = (theta[i] + theta_1c * cos(j * dpsi) + theta_1s * sin(j * dpsi));  // theta in tersm of r i.e.theta_twd = E + Fr
            /*cout << "u_inf: " << u_inf << ", alpha: " << alpha << ", L_i_calculated: " << L_i_calculated << ", Vv: " << Vv << ", radii[i]" << radii[i] << 
                ", dbeta_dt:" << dbeta_dt << ", beta: " << beta << ", j*dpsi" << j*dpsi << std::endl;*/
            const float Up = u_inf * sin((alpha)) + (OMEGA * R_tip * (L_i_calculated)+Vv) + radii[i] * dbeta_dt + u_inf * sin(beta) * cos(j * dpsi);
            const float Ut = OMEGA * radii[i] + u_inf * cos(alpha) * sin(j * dpsi);

            Up_BEM_array[j][i] = Up;

            //cout << "rho: " << rho << ", chord[i]: " << chord[i] << ", a" << a << ", theta_equation: " << theta_equation << ", Ut: " << Ut << ", Up: " << Up << std::endl;

            dL_array[i] = (0.5 * rho * chord[i] * a * (theta_equation * pow(Ut, 2) - Ut * Up)); // filling the  corresponding array
            dP_array[i] = dL_array[i] * radii[i] * cos(j * dpsi);  // filling the  corresponding array
            dR_array[i] = dL_array[i] * radii[i] * sin(j * dpsi);  // filling the  corresponding array

            const float dDi = 0.5 * rho * chord[i] * a * (theta_equation - atan(Up / Ut)) * Up * Ut; // Di term in induced torque equation
            dQi_array[i] = dDi * radii[i]; // filling the  corresponding array

            const float dDo = 0.5 * rho * Ut * Ut * chord[i] * (Cd0); // Do term in  profile torque equation
            dQo_array[i] = dDo * radii[i]; // filling the  corresponding array

            dPi_array[i] = dDi * Ut; // filling the  corresponding array
            dPo_array[i] = dDo * Ut; // filling the  corresponding array

            Thrust_plot_array[j][i] = dL_array[i];
            circulation[j][i] = (dL_array[i]) / (rho * Ut);
        }
        dT_array[j] = trapz(dL_array, radii); // thust produced at given azimuthal location

        dOPM_array[j] = trapz(dP_array, radii); // pitching moment at given azimuthal 
        dORM_array[j] = trapz(dR_array, radii); // rolling moment at given azimuthal 
        dOQi_array[j] = trapz(dQi_array, radii); // induced torque at given azimuthal location                                              #filling in  corresponding array
        dOQo_array[j] = trapz(dQo_array, radii); // profile torque at given azimuthal 
        dOPi_array[j] = trapz(dPi_array, radii); // induced power at given azimuthal location
        dOPo_array[j] = trapz(dPo_array, radii); // profile power at given azimuthal location
    }

    const float THRUST = trapz(dT_array, psi_array) * (N_BLADES / (2 * M_PI)); // thust produced at given azimuthal location

    const float OVERALL_PITCHING_MOMENT = trapz(dOPM_array, psi_array) * (N_BLADES / (2 * M_PI));   // Overall Pitching Moment produced by rotor
    const float OVERALL_ROLLING_MOMENT = trapz(dORM_array, psi_array) * (N_BLADES / (2 * M_PI));  // Overall Rolling Moment produced by rotor
    const float INDUCED_TORQUE = trapz(dOQi_array, psi_array) * (N_BLADES / (2 * M_PI)); // Induced Torque produced by rotor
    const float PROFILE_TORQUE = trapz(dOQo_array, psi_array) * (N_BLADES / (2 * M_PI)); // Profile Torque produced by rotor

    const float TOTAL_TORQUE = INDUCED_TORQUE + PROFILE_TORQUE; // Total Torque produced by rotor

    const float INDUCED_POWER = trapz(dOPi_array, psi_array) * (N_BLADES / (2 * M_PI)); // Induced Power produced by rotor
    const float PROFILE_POWER = trapz(dOPo_array, psi_array) * (N_BLADES / (2 * M_PI)); // Profile Power produced by rotor
    const float TOTAL_POWER = (INDUCED_POWER)+PROFILE_POWER; // Total Power produced by rotor

    cout << ">>>>>>>>>>>>  RESULTS  <<<<<<<<<<<<" << std::endl;
    cout << "Thrust at Trim condition (N):" << (round(THRUST * 1000.0) / 1000.0) << std::endl;
    cout << "Thrust at Trim condition (N):" << THRUST << std::endl;
    cout << "Overall Pitching Moment (Nm):" << (round(OVERALL_PITCHING_MOMENT * 1000.0) / 1000.0) << std::endl;
    cout << "Overall Rolling Moment (Nm):" << (round(OVERALL_ROLLING_MOMENT * 1000.0) / 1000.0) << std::endl;
    cout << "Induced Torque (Nm):" << (round(INDUCED_TORQUE * 1000.0) / 1000.0) << std::endl;
    cout << "Profile Torque (Nm):" << (round(PROFILE_TORQUE * 1000.0) / 1000.0) << std::endl;
    cout << "Total Torque (Nm):" << (round(TOTAL_TORQUE * 1000.0) / 1000.0) << std::endl;
    cout << "Induced Power (kW):" << abs(INDUCED_POWER * 0.001) << std::endl;
    cout << "Profile Power (kW):" << (round(PROFILE_POWER * 0.001 * 1000.0) / 1000.0) << std::endl;
    cout << "Total Power (kW):" << (round(TOTAL_POWER * 0.001 * 1000.0) / 1000.0) << std::endl;

    const float CT = THRUST / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));  // Thrust co efficient 
    CP = TOTAL_POWER / (rho * DISC_AREA * pow(OMEGA * R_tip, 3));  // Thrust co efficient 
    cout << "ThrustA2:" << round(Thrust * 1000.0)/1000.0 << "N and CT:" << CT << std::endl;
    cout << "PowerA2:" << round(Power)/1000.0 << "kW and CP:" << CP << std::endl;
    cout << "Thrust:" << round(THRUST * 1000.0)/1000.0 << "N and CT:" << CT << std::endl;
    cout << "Power:" << round(TOTAL_POWER) / 1000.0 << "kW and CP:" << CP << std::endl;

    //const float THRUST_ = 8786.417;
    //const float CT_ = THRUST / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));
    //vector< vector<vector<vector<float>>>> hoverWakerSoln = HoverWake(CT);
    ////save4DArray(hoverWakerSoln, "hover_wake_sol.npy");
    ////saveData(hoverWakerSoln);
    //std::cout << "JM: " << JM;
    //std::cout << "CT: " << CT;
    //std::cout << "dpsi: " << dpsi;

    return 0;
}


void EX2BMethod(const float& CT, vector<vector<float>>& circulation) {
    vector<vector<vector<vector<float>>>> positions_old(JM, vector<vector<vector<float>>>(span_sections, vector<vector<float>>(3, vector<float>(N_BLADES, 0))));
    vector<vector<vector<vector<float>>>> positions_new(JM, vector<vector<vector<float>>>(span_sections, vector<vector<float>>(3, vector<float>(N_BLADES, 0))));
    vector<vector<vector<vector<float>>>> positions_pc(JM, vector<vector<vector<float>>>(span_sections, vector<vector<float>>(3, vector<float>(N_BLADES, 0))));
    vector<vector<vector<vector<float>>>> blade_positions(JM, vector<vector<vector<float>>>(span_sections, vector<vector<float>>(3, vector<float>(N_BLADES, 0))));
    vector<vector<vector<vector<float>>>> positions_compare(JM, vector<vector<vector<float>>>(span_sections, vector<vector<float>>(3, vector<float>(N_BLADES, 0))));

    vector<vector<vector<vector<float>>>> result = HoverWake(CT);
    vector<vector<vector<float>>> array_Xprs = result[0];
    vector<vector<vector<float>>> array_Yprs = result[1];
    vector<vector<vector<float>>> array_Zprs = result[2];

    for (int j = 0; j < JM; j++)
    {
        for (int i = 0; i < span_sections; i++)
        {
            for (int k = 0; k < N_BLADES; k++)
            {
                positions_old[j][i][0][k] = array_Xprs[j][i][k];
                positions_old[j][i][1][k] = array_Yprs[j][i][k];
                positions_old[j][i][2][k] = array_Zprs[j][i][k];


                positions_pc[j][i][0][k] = array_Xprs[j][i][k];
                positions_pc[j][i][0][k] = array_Yprs[j][i][k];
                positions_pc[j][i][0][k] = array_Zprs[j][i][k];
                
                positions_compare = positions_old;
                
                const float C_psi = 0;
                vector<float> RMS_arr;
                vector<float> Power_arr_iter;
                vector<float> Thrust_arr_iter;
                const float rms_change = 1;
                const int Rotation_count = 0;
                const float rotation = 0;
                const int count = 0;
                const int az_index = 0;
                const int shed_index = 1;
                const float dbeta_dt = 0;
                const float beta = 0;

                vector<vector<float>> Azimuthal_dT_array(N, vector<float>(N_BLADES, 0.0f));
                vector<vector<float>> Azimuthal_dOPi_array(N, vector<float>(N_BLADES, 0.0f));
                vector<vector<float>> Azimuthal_dOPo_array(N, vector<float>(N_BLADES, 0.0f));
                vector<vector<float>> Azimuthal_PM_array(N, vector<float>(N_BLADES, 0.0f));
                vector<vector<float>> Azimuthal_RM_array(N, vector<float>(N_BLADES, 0.0f));
                vector<vector<float>> Azimuthal_dOQo_array(N, vector<float>(N_BLADES, 0.0f));
                vector<vector<float>> Azimuthal_dOQi_array(N, vector<float>(N_BLADES, 0.0f));

                vector<vector<vector<vector<float>>>> V_IND(JM, vector<vector<vector<float>>>(3, vector<vector<float>>(3, vector<float>(N_BLADES, 0.0f))));
                vector<vector<vector<vector<float>>>> CirVM_array(JM+1, vector<vector<vector<float>>>(span_sections+1, vector<vector<float>>(N_BLADES, vector<float>(2, 0.0f))));
                vector<vector<vector<vector<float>>>> Thrust_contours(JM, vector<vector<vector<float>>>(span_sections, vector<vector<float>>(N_BLADES, vector<float>(2, 0.0f))));

                vector<vector<vector<vector<float>>>> Shed_circulation(JM, vector<vector<vector<float>>>(span_sections+1, vector<vector<float>>(N_BLADES, vector<float>(2, 0.0f))));
                vector<vector<vector<vector<float>>>> Trail_circulation(JM, vector<vector<vector<float>>>(span_sections+1, vector<vector<float>>(N_BLADES, vector<float>(2, 0.0f))));

                for (int bc = 0; bc < N_BLADES; bc++)
                {
                    for (int ic = 0; ic < span_sections; ic++)
                    {
                        for (int j = 0; j < JM; j++)
                        {
                            CirVM_array[j][ic][bc][0] = circulation[j][ic];
                            CirVM_array[j][ic][bc][1] = circulation[j][ic];
                        }
                    }
                }

            }
        }
    }

}

