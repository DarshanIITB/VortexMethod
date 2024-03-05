#include <vector>
#include <cmath>
#include <corecrt_math_defines.h>
#include <numeric>
#include <tuple>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include "Utils.hpp"
#include "json.hpp"
//#include "./coordinates/Coordinates.h"

using std::vector, std::pow, std::sqrt, std::cout;

const double rho = 1.225f;

const int span_sections = 10; // number of points used per blade
const int N = 24; // number if turns in one rotation
const int NR = 7; // number of rotations
const int NR_F = 2;
const int FAR_FIELD_CUTOFF = (NR - NR_F) * N;
const double dpsi_d = (360 / N);
const double mu = 0.0f;

const double alpha_d = 0;
const double alpha = alpha_d * (M_PI / 180);
double alpha_p = -alpha;

const double RPM = 383.0f; // rpm at cruise
const double OMEGA = 2 * M_PI * RPM / 60;

const double a = 5.75; // lift curve slope per radian from literature
const double Cd0 = 0.022;  // co - efficient of drag from literature

const double Vv = 0.0f; // tip path plane in radians

const double R_tip = 5.5f;
const double u_inf = mu * (OMEGA * R_tip) / cos(alpha); // this is Y direction
const double v_inf = 0; // this is X direction
const double w_inf = 0 - u_inf * sin(alpha) - Vv;
vector<double>V_vec = { v_inf, u_inf, w_inf };

const double R_hub = 0.1 * R_tip;
vector<double> radii = linspace(R_hub, R_tip, span_sections);
const double chord_at_hub = 0.325;
const double taper_ratio = 1.0f;
const double chord_at_tip = chord_at_hub * taper_ratio;
const double taper_slope = (chord_at_tip - chord_at_hub) / (R_tip - R_hub);
const double chord_at_root = chord_at_hub - taper_slope * R_hub;
vector<double> chord = scalarAdd(chord_at_root, scalarMul(taper_slope, radii));

const int N_BLADES = 2;
const int N_ROTORS = 1;
const int BLADES_PER_ROTOR = N_BLADES / N_ROTORS;

const double x_displacement = 1.5 * R_tip;
const double y_displacement = 0 * R_tip;
const double z_displacement = 0.2 * R_tip;

const double d = sqrt(pow(x_displacement, 2) + pow(y_displacement, 2));

const double OVERLAP_AREA = 2 * pow(R_tip, 2) * std::acos(d / (2 * R_tip)) - 0.5 * sqrt(4 * pow(R_tip, 2) - pow(d, 2));
const double DISC_AREA = (N_ROTORS > 1) ? M_PI * pow(R_tip, 2) * N_ROTORS - OVERLAP_AREA : M_PI * pow(R_tip, 2) * N_ROTORS;

const double dpsi = dpsi_d * (M_PI / 180);
const double NT = (NR * (360)) / dpsi_d;
const int JM = int(NT);

const double alp = 1.2;
const double K_visc = 0.000015;
const double del_s = 1000;

const double THETA_TIP = 2.0f;
const double THETA_HUB = 12.0f;
const double F = (THETA_TIP - THETA_HUB) / (R_tip - R_hub);
const double E = THETA_HUB - F * R_hub;
const double lc = Vv / (OMEGA * R_tip);

const double blade_twist = (THETA_TIP - THETA_HUB) * (M_PI / 180);

const double h_by_r = 0.0f; // for ground effect

vector<double> theta = scalarAdd(E * M_PI / 180, (scalarMul(F * M_PI / 180, radii)));
vector<double> solidity = scalarMul(N_BLADES/(M_PI * R_tip), chord);

const int MAX_ITERS_ALLOWED = 50;

vector<vector<double>> bsl_2(vector<vector<vector<vector<double>>>>& positions, vector<vector<vector<vector<double>>>>& trail_circulation, int tc_index, int blade_index, int span_index, vector <vector< double >>& v_ind)
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
                    vector<double> v_ind_c = vector<double>(3, 0.0f);
                    if (k != j)
                    {
                        const double d_gamma = trail_circulation[k][i][blade_index][tc_index] - trail_circulation[k + 1][i][blade_index][tc_index];

                        const double xa = positions[k][i][0][blade_index];
                        const double xb = positions[k+1][i][0][blade_index];
                        const double ya = positions[k][i][1][blade_index];
                        const double yb = positions[k+1][i][1][blade_index];
                        const double za = positions[k][i][2][blade_index];
                        const double zb = positions[k+1][i][2][blade_index];

                        // calculating differences for optimization
                        const double xb_minus_xa = xb - xa;
                        const double yb_minus_ya = yb - ya;
                        const double zb_minus_za = zb - za;
                        
                        const double xc_minus_xb = xc - xb;
                        const double yc_minus_yb = yc - yb;
                        const double zc_minus_zb = zc - zb;

                        // distance b/w A and B
                        const double L = sqrt(pow(xb_minus_xa, 2) + pow(yb_minus_ya, 2) + pow(zb_minus_za, 2));

                        const double cutoff = sqrt(pow(0.25f * L, 4) + 4 * alp * del_s * K_visc * ((k + 1) * dpsi / OMEGA));
                        // distance b/w A and C
                        const double r1 = sqrt(pow(xc - xa, 2) + pow(yc - ya, 2) + pow(zc - za, 2)) + cutoff;

                        // distance b/w B and C
                        const double r2 = sqrt(pow(xc_minus_xb, 2) + pow(yc_minus_yb, 2) + pow(zc_minus_zb, 2)) + cutoff;

                        const double multiplier = d_gamma * (r1 + r2) / (2 * M_PI * r1 * r2 * (pow(r1 + r2, 2) - pow(L, 2)));

                        double u_bs = multiplier * (zc_minus_zb * yb_minus_ya - yc_minus_yb * zc_minus_zb);
                        double v_bs = multiplier * (xc_minus_xb * zb_minus_za - zc_minus_zb * xb_minus_xa);
                        double w_bs = multiplier * (yc_minus_yb * xb_minus_xa - xc_minus_xb * yb_minus_ya);

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
        for (int i = 0; i < N; i++)
        {
            v_ind[FAR_FIELD_CUTOFF + N * rot + i] = v_ind[FAR_FIELD_CUTOFF - N + i];
        }
    }
    return v_ind;
}

vector<double> bsl_perf(vector<vector<vector<vector<double>>>>& positions, vector<vector<vector<vector<double>>>>& trail_circulation, int tc_index, int blade_index, int span_index, vector<double>& v_ind)
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
                const double d_gamma = trail_circulation[k][i][blade_index][tc_index] - trail_circulation[k + 1][i][blade_index][tc_index];

                const double xa = positions[k][i][0][blade_index];
                const double xb = positions[k + 1][i][0][blade_index];
                const double ya = positions[k][i][1][blade_index];
                const double yb = positions[k + 1][i][1][blade_index];
                const double za = positions[k][i][2][blade_index];
                const double zb = positions[k + 1][i][2][blade_index];

                vector<double> v_ind_c = vector<double>(3, 0.0f);

                // calculating differences for optimization
                const double xb_minus_xa = xb - xa;
                const double yb_minus_ya = yb - ya;
                const double zb_minus_za = zb - za;

                const double xc_minus_xb = xc - xb;
                const double yc_minus_yb = yc - yb;
                const double zc_minus_zb = zc - zb;

                // distance b/w A and B
                const double L = sqrt(pow(xb_minus_xa, 2) + pow(yb_minus_ya, 2) + pow(zb_minus_za, 2));

                const double cutoff = sqrt(pow(0.25f * L, 4) + 4 * alp * del_s * K_visc * ((k + 1) * dpsi / OMEGA));
                // distance b/w A and C
                const double r1 = sqrt(pow(xc - xa, 2) + pow(yc - ya, 2) + pow(zc - za, 2)) + cutoff;

                // distance b/w B and C
                const double r2 = sqrt(pow(xc_minus_xb, 2) + pow(yc_minus_yb, 2) + pow(zc_minus_zb, 2)) + cutoff;

                const double multiplier = d_gamma * (r1 + r2) / (2 * M_PI * r1 * r2 * (pow(r1 + r2, 2) - pow(L, 2)));

                double u_bs = multiplier * (zc_minus_zb * yb_minus_ya - yc_minus_yb * zc_minus_zb);
                double v_bs = multiplier * (xc_minus_xb * zb_minus_za - zc_minus_zb * xb_minus_xa);
                double w_bs = multiplier * (yc_minus_yb * xb_minus_xa - xc_minus_xb * yb_minus_ya);

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

//void save4DArray(const vector<vector<vector<vector<double>>>>& array, const char* filename) {
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
//                        file.write(reinterpret_cast<const char*>(&row), dim4 * sizeof(double));
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

//void save4DArray(const std::vector<std::vector<std::vector<std::vector<double>>>>& array, const char* filename) {
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
//                        file.write(reinterpret_cast<const char*>(&row), sizeof(double));
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

void saveData(vector<vector<vector<vector<double>>>>& array4D)
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
                    outfile.write(reinterpret_cast<const char*>(&element), sizeof(double));
                }
            }
        }
    }
    outfile.close();
}

vector< vector<vector<vector<double>>>> HoverWake(double CT_input)
{
    const double CT = CT_input / N_ROTORS;
    const double A = 0.78;
    const double lam = 0.145 + 27 * CT;
    vector<vector<double>> Zv_array(JM, vector<double>(span_sections, 0));
    vector<vector<vector<double>>> arrayX1(JM, vector<vector<double>>(span_sections, vector<double>(N_BLADES, 0)));
    vector<vector<vector<double>>> arrayY1(JM, vector<vector<double>>(span_sections, vector<double>(N_BLADES, 0)));
    vector<vector<vector<double>>> arrayZ1(JM, vector<vector<double>>(span_sections, vector<double>(N_BLADES, 0)));
    for (int i = 0; i < span_sections; i++)
    {
        if (i < span_sections - 1)
        {
            for (int j = 0; j < JM; j++) // time loop
            {
                const double psi = j * dpsi;
                const double R_vortex = radii[i] * (A + (1 - A) * std::exp(-lam * psi));
                const double chord_ = chord[i];
                const double solidity_ = (BLADES_PER_ROTOR * chord_) / (M_PI * R_tip);
                const double K1 = -2.2 * sqrt(0.5 * CT);
                const double K2 = -2.7 * sqrt(0.5 * CT);
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
                        std::pair<double, double> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
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
                        std::pair<double, double> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
                        arrayX1[j][i][bn] = new_coords.first + x_displacement;
                        arrayY1[j][i][bn] = new_coords.second + y_displacement;
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
                const double psi = j * dpsi;
                const double R_vortex = radii[i] * (A + (1 - A) * std::exp(-lam * psi));
                const double chord_ = chord[i];
                const double solidity_ = (BLADES_PER_ROTOR * chord_) / (M_PI * R_tip);
                const double K1 = -0.25 * ((CT / solidity_) + 0.001 * blade_twist * (180 / M_PI));
                const double K2 = -(1.41 + 0.0141 * blade_twist * (180 / M_PI)) * sqrt(0.5 * CT);
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
                        std::pair<double, double> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
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
                        std::pair<double, double> new_coords = rotate(R_vortex, 0.0f, 1.0f * psi + (bn * 2 * M_PI / BLADES_PER_ROTOR));
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

void EX2BMethod(const double& CT, vector<vector<double>>& circulation)
{
    vector<double> psi_array = linspace(0, 2 * M_PI, N);
    vector<vector<double>>V_INF(JM, vector<double>(3, 0.0f));
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < JM; j++)
        {
            V_INF[j][i] = V_vec[i];
        }
    }


    const double theta_1s = (0) * (M_PI / 180);  // theta_1s  in degrees at trim condition
    const double theta_1c = (0) * (M_PI / 180);  // theta_1c  in degrees at trim condition

    vector<vector<vector<vector<double>>>> positions_old(JM, vector<vector<vector<double>>>(span_sections, vector<vector<double>>(3, vector<double>(N_BLADES, 0))));
    //vector<vector<vector<vector<double>>>> positions_old(N_BLADES, vector<vector<vector<double>>>(JM, vector<vector<double>>(span_sections, vector<double>(3, 0))));
    vector<vector<vector<vector<double>>>> positions_new(JM, vector<vector<vector<double>>>(span_sections, vector<vector<double>>(3, vector<double>(N_BLADES, 0))));
    vector<vector<vector<vector<double>>>> positions_pc(JM, vector<vector<vector<double>>>(span_sections, vector<vector<double>>(3, vector<double>(N_BLADES, 0))));
    vector<vector<vector<vector<double>>>> blade_positions(JM, vector<vector<vector<double>>>(span_sections, vector<vector<double>>(3, vector<double>(N_BLADES, 0))));
    vector<vector<vector<vector<double>>>> positions_compare(JM, vector<vector<vector<double>>>(span_sections, vector<vector<double>>(3, vector<double>(N_BLADES, 0))));

    vector<double>Thrust_new(N_ROTORS, 0.0f);
    vector<double>Induced_Power(N_ROTORS, 0.0f);
    vector<double>Profile_Power(N_ROTORS, 0.0f);
    vector<double>Power_new(N_ROTORS, 0.0f);
    vector<double>Torque_new(N_ROTORS, 0.0f);
    vector<double>OPM_new(N_ROTORS, 0.0f);
    vector<double>ORM_new(N_ROTORS, 0.0f);

    vector<vector<vector<vector<double>>>> result = HoverWake(CT);
    vector<vector<vector<double>>> array_Xprs = result[0];
    vector<vector<vector<double>>> array_Yprs = result[1];
    vector<vector<vector<double>>> array_Zprs = result[2];

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
            }
        }
    }
    int C_psi = 0;
    vector<double> RMS_arr;
    vector<double> Power_arr_iter;
    vector<double> Thrust_arr_iter;
    double rms_change = 1;
    int Rotation_count = 0;
    double rotation = 0;
    int count = 0;
    int az_index = 0;
    const int shed_index = 1;
    const double dbeta_dt = 0;
    const double beta = 0;

    vector<vector<double>> Azimuthal_dT_array(N, vector<double>(N_BLADES, 0.0f));
    vector<vector<double>> Azimuthal_dOPi_array(N, vector<double>(N_BLADES, 0.0f));
    vector<vector<double>> Azimuthal_dOPo_array(N, vector<double>(N_BLADES, 0.0f));
    vector<vector<double>> Azimuthal_PM_array(N, vector<double>(N_BLADES, 0.0f));
    vector<vector<double>> Azimuthal_RM_array(N, vector<double>(N_BLADES, 0.0f));
    vector<vector<double>> Azimuthal_dOQo_array(N, vector<double>(N_BLADES, 0.0f));
    vector<vector<double>> Azimuthal_dOQi_array(N, vector<double>(N_BLADES, 0.0f));

    vector<vector<vector<vector<double>>>> V_IND(JM, vector<vector<vector<double>>>(3, vector<vector<double>>(3, vector<double>(N_BLADES, 0.0f))));
    vector<vector<vector<vector<double>>>> CirVM_array(JM+1, vector<vector<vector<double>>>(span_sections+1, vector<vector<double>>(N_BLADES, vector<double>(2, 0.0f))));
    vector<vector<vector<vector<double>>>> Thrust_contours(JM, vector<vector<vector<double>>>(span_sections, vector<vector<double>>(N_BLADES, vector<double>(2, 0.0f))));

    vector<vector<vector<vector<double>>>> Shed_circulation(JM, vector<vector<vector<double>>>(span_sections+1, vector<vector<double>>(N_BLADES, vector<double>(2, 0.0f))));
    vector<vector<vector<vector<double>>>> Trail_circulation(JM, vector<vector<vector<double>>>(span_sections+1, vector<vector<double>>(N_BLADES, vector<double>(2, 0.0f))));

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
    for (int j = 0; j < JM; j++)
    {
        for (int i = 0; i < span_sections; i++)
        {
            Trail_circulation[j][i] = vectorSub(CirVM_array[j][i], CirVM_array[j][i + 1]);
        }
    }
    vector<vector<vector<double>>> Blade_velocities(span_sections, vector<vector<double>>(3, vector<double>(N_BLADES, 0.0f)));
    vector<vector<double>> Up_VM_arr(N, vector<double>(span_sections, 0.0f));
    while (rms_change > 0.01)
    {
        //cout << "Iteration " << count;
        count++;
        rotation = rotation + dpsi * (180 / M_PI);
        //C_psi = C_psi + 1;
        C_psi++;
        vector<vector<double>>Radial_dT_array(span_sections, vector<double>(N_BLADES, 0.0f));
        vector<vector<double>>Radial_dOPi_array(span_sections, vector<double>(N_BLADES, 0.0f));
        vector<vector<double>>Radial_dOPo_array(span_sections, vector<double>(N_BLADES, 0.0f));
        vector<vector<double>>Radial_PM_array(span_sections, vector<double>(N_BLADES, 0.0f));
        vector<vector<double>>Radial_RM_array(span_sections, vector<double>(N_BLADES, 0.0f));
        vector<vector<double>>Radial_dOQo_array(span_sections, vector<double>(N_BLADES, 0.0f));
        vector<vector<double>>Radial_dOQi_array(span_sections, vector<double>(N_BLADES, 0.0f));
        for (int bl = 0; bl < N_BLADES; bl++)
        {
            for (int ir = 0; ir < span_sections; ir++)
            {
                double psi, psi_p;
                if (bl <= BLADES_PER_ROTOR)
                {

                    //Here was the error
                    //alpha_p *= -alpha;
                    alpha_p = -alpha;
                    psi = (1 * count * dpsi + ((bl) / BLADES_PER_ROTOR) * 2 * M_PI);
                    psi_p = (-1 * count * dpsi + ((bl) / BLADES_PER_ROTOR) * 2 * M_PI);
                    positions_new[0][ir][1][bl] = radii[ir] * ((cos(psi_p) * cos(beta)) + sin(beta) * sin(alpha_p));
                    positions_new[0][ir][0][bl] = radii[ir] * (sin(psi_p) * cos(beta));
                    positions_new[0][ir][2][bl] = radii[ir] * (sin(beta) * cos(alpha_p) - cos(psi_p) * cos(beta) * sin(alpha_p)) + h_by_r * R_tip;
                }
                else
                {
                    psi = (1 * count * dpsi + ((bl) / BLADES_PER_ROTOR) * 2 * M_PI);
                    psi_p = (1 * count * dpsi + ((bl) / BLADES_PER_ROTOR) * 2 * M_PI);
                    positions_new[0][ir][1][bl] = 1 * x_displacement + radii[ir] * ((cos(psi_p) * cos(beta)) + sin(beta) * sin(alpha_p));
                    positions_new[0][ir][0][bl] = 1 * y_displacement + radii[ir] * (sin(psi_p) * cos(beta));
                    positions_new[0][ir][2][bl] = array_Zprs[0][ir][bl] + radii[ir] * (sin(beta) * cos(alpha_p) - cos(psi_p) * cos(beta) * sin(alpha_p));
                }
                vector<vector<double>>V_ind_geometry(JM, vector<double>(3, 0.0f));
                vector<vector<double>> res = bsl_2(positions_old, Trail_circulation, 0, bl, ir, V_ind_geometry);
                //cout << "res dimensions: (" << res.size() << ", " << res[0].size() << ")" << std::endl;
                //cout << "V_ind dimensions: (" << V_IND.size() << ", " << V_IND[0].size() << ", " << V_IND[0][0].size() << ", " << V_IND[0][0][0].size() << ")" << std::endl;
                for (int j = 0; j < JM; j++)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        //cout << "j: " << j << ", i: " << i << std::endl;
                        V_IND[j][i][0][bl] = res[j][i];
                    }
                }

                for (int j = 1; j < JM; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        positions_new[j][ir][k][bl] = positions_old[j - 1][ir][k][bl] + (dpsi / OMEGA) * (V_INF[j][k] + (1.0f / 3) * (V_IND[j][k][0][bl] + V_IND[j][k][1][bl] + V_IND[j][k][2][bl]));
                    }
                }
                vector<double> V_ind_perf(3, 0.0f);
                vector<double> res_perf = bsl_perf(positions_old, Trail_circulation, 0, bl, ir, V_ind_perf);
                for (int k = 0; k < 3; k++)
                {
                    Blade_velocities[ir][k][bl] = res_perf[k];
                }
                for(int j = 0; j < JM; j++)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        V_IND[j][i][2][bl] = V_IND[j][i][1][bl];
                        V_IND[j][i][1][bl] = V_IND[j][i][0][bl];
                    }
                }
                const double theta_equation = (theta[ir] + theta_1c * cos(psi) + theta_1s * sin(psi));

                const double Up = abs(Blade_velocities[ir][2][bl]) + u_inf * sin(alpha) + u_inf * sin(beta) * cos(psi); // - abs(Vv * (1 + sin(10 * psi)))
                    // Up = PTC_inflow_array[ir] * omega * R_tip + u_inf * sin(alpha) + u_inf * sin(beta) * cos(psi)  #  - abs(Vv * (1 + sin(10 * psi)))

                const double Ut = OMEGA * radii[ir] + u_inf * cos(alpha) * sin(psi);
                Up_VM_arr[az_index][ir] = Up;
                const double L_bar = 0.5 * rho * chord[ir] * a * (theta_equation - atan(abs(Up) / Ut)) * (pow(Up, 2) + pow(Ut, 2)); // filling the  corresponding array

                const double dDi = 0.5 * rho * chord[ir] * a * (theta_equation - atan(Up / Ut)) * Up * Ut; // Di term in induced torque equation
                const double dDo = 0.5 * rho * Ut * Ut * chord[ir] * (Cd0); // Do term in  profile torque equation
                const double dQi = dDi * radii[ir];
                const double dQo = dDo * radii[ir];
                const double dRM = L_bar * radii[ir] * sin(abs(psi));
                const double dPM = 1 * L_bar * radii[ir] * cos(psi);

                for (int j = 1; j < JM; j++)
                {
                    CirVM_array[j][ir][bl][1] = CirVM_array[j-1][ir][bl][1];
                    Thrust_contours[j][ir][bl][1] = Thrust_contours[j-1][ir][bl][1];
                }

                CirVM_array[0][ir][bl][1] = CirVM_array[0][ir][bl][0] + 0.1 * (abs(L_bar / (rho * Ut)) - CirVM_array[0][ir][bl][0]);

                Thrust_contours[0][ir][bl][1] = L_bar;

                for (int j = 0; j < JM; j++)
                {
                    Shed_circulation[j][ir][bl][1] = CirVM_array[j][ir][bl][1] - CirVM_array[0][ir][bl][0];
                }
                for (int j = 2; j < JM; j++)
                {
                    Trail_circulation[j][ir][bl][1] = Trail_circulation[j - 1][ir][bl][1];
                }
                Trail_circulation[1][ir][bl][1] = CirVM_array[0][ir][bl][1] - CirVM_array[0][ir + 1][bl][1];
                Trail_circulation[0][ir][bl][1] = CirVM_array[0][ir][bl][1];

                const double dPi = dDi * Ut; // filling the  corresponding array
                const double dPo = dDo * Ut;

                Radial_dT_array[ir][bl] = L_bar;
                Radial_dOPi_array[ir][bl] = dPi;
                Radial_dOPo_array[ir][bl] = dPo;

                Radial_dOQi_array[ir][bl] = dQi;
                Radial_dOQo_array[ir][bl] = dQo;
                Radial_PM_array[ir][bl] = dPM;
                Radial_RM_array[ir][bl] = dRM;
            }
        }
        positions_old = positions_new;

        for (int bx = 0; bx < N_BLADES; bx++)
        {
            Azimuthal_dT_array[az_index][bx] = trapz(Radial_dT_array, bx, radii) * (1 / (2 * M_PI));
            Azimuthal_dOPi_array[az_index][bx] = trapz(Radial_dOPi_array, bx, radii) * (1 / (2 * M_PI));
            Azimuthal_dOPo_array[az_index][bx] = trapz(Radial_dOPo_array, bx, radii) * (1 / (2 * M_PI));
            Azimuthal_dOQi_array[az_index][bx] = trapz(Radial_dOQi_array, bx, radii) * (1 / (2 * M_PI));
            Azimuthal_dOQo_array[az_index][bx] = trapz(Radial_dOQo_array, bx, radii) * (1 / (2 * M_PI));
            Azimuthal_PM_array[az_index][bx] = trapz(Radial_PM_array, bx, radii) * (1 / (2 * M_PI));
            Azimuthal_RM_array[az_index][bx] = trapz(Radial_RM_array, bx, radii) * (1 / (2 * M_PI));
        }

        for (int j = 0; j < JM; j++)
        {
            for (int i = 0; i < span_sections; i++)
            {
                for (int k = 0; k < N_BLADES; k++)
                {
                    CirVM_array[j][i][k][0] = CirVM_array[j][i][k][1];
                    Shed_circulation[j][i][k][0] = Shed_circulation[j][i][k][1];
                    Trail_circulation[j][i][k][0] = Trail_circulation[j][i][k][1];
                    Thrust_contours[j][i][k][0] = Thrust_contours[j][i][k][1];
                }
            }
        }
        az_index++;
        //cout << count << " ";
        if (C_psi % N == 0)
        {
            //cout << std::endl;
            Rotation_count += 1;
            vector<vector<double>> Change_arr(JM, vector<double>(3));
            for (int j = 0; j < JM; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    Change_arr[j][k] = positions_new[j][span_sections - 1][k][N_BLADES - 1] - positions_compare[j][span_sections - 1][k][N_BLADES - 1];
                }
            }
            double shift_sum = 0;
            /*dx2_dy2_dz2 = np.square(Change_array);
            dx2_dy2_dz2_sum = dx2_dy2_dz2[:, 0] + dx2_dy2_dz2[:, 1] + dx2_dy2_dz2[:, 2];*/
            for (int s = 0; s < JM; s++)
            {
                for (int t = 0; t < 3; t++)
                {
                    shift_sum += pow(Change_arr[s][t], 2);
                }
            }
            // Eucladian_distance = np.sqrt(dx2_dy2_dz2[:, 0] + dx2_dy2_dz2[:, 1] + dx2_dy2_dz2[:, 2]);
            //shift_sum = np.sum(dx2_dy2_dz2_sum);
            const double err1 = sqrt(shift_sum / JM);
            rms_change = err1 / R_tip;

            RMS_arr.push_back(rms_change);

            for (int by = 0; by < N_BLADES; by++)
            {
                int nx = by <= BLADES_PER_ROTOR ? 0 : 1;
                Thrust_new[nx] = Thrust_new[nx] + trapz(Azimuthal_dT_array, by, psi_array);
                Induced_Power[nx] = Induced_Power[nx] + (trapz(Azimuthal_dOPi_array, by, psi_array));
                Profile_Power[nx] = Profile_Power[nx] + (trapz(Azimuthal_dOPo_array, by, psi_array));
                Power_new[nx] = Power_new[nx] + trapz(Azimuthal_dOPi_array, by, psi_array) + trapz(Azimuthal_dOPo_array, by, psi_array);
                Torque_new[nx] = Torque_new[nx] + trapz(Azimuthal_dOQi_array, by, psi_array) + trapz(Azimuthal_dOQo_array, by, psi_array);
                OPM_new[nx] = OPM_new[nx] + trapz(Azimuthal_PM_array, by, psi_array);
                ORM_new[nx] = ORM_new[nx] + trapz(Azimuthal_RM_array, by, psi_array);
            }

            vector<vector<double>>Azimuthal_dT_array(N, vector<double>(N_BLADES, 0.0f));
            vector<vector<double>>Azimuthal_dOPi_array(N, vector<double>(N_BLADES, 0.0f));
            vector<vector<double>>Azimuthal_dOPo_array(N, vector<double>(N_BLADES, 0.0f));
            vector<vector<double>>Azimuthal_PM_array(N, vector<double>(N_BLADES, 0.0f));
            vector<vector<double>>Azimuthal_RM_array(N, vector<double>(N_BLADES, 0.0f));
            vector<vector<double>>Azimuthal_dOQo_array(N, vector<double>(N_BLADES, 0.0f));
            vector<vector<double>>Azimuthal_dOQi_array(N, vector<double>(N_BLADES, 0.0f));


            cout << Rotation_count << " " << round_(rms_change, 5) << " " << round_(sumVecElem(Thrust_new), 3) << " " << round_(sumVecElem(scalarMul(0.001, Induced_Power)), 3) << " " << round_(sumVecElem(scalarMul(0.001, Profile_Power)), 3) << " " << round_(sumVecElem(scalarMul(0.001, Power_new)), 3) << " " << round_(sumVecElem(Torque_new), 2) << " " << round_(sumVecElem(OPM_new), 2) << " " << round_(sumVecElem(ORM_new), 2) << std::endl;

            positions_compare = positions_new;

            az_index = 0;

        }
        //Rotation_count++;
        if (Rotation_count == MAX_ITERS_ALLOWED) break;
    }
    cout << " Rotor Thrust  Induced Power   Profile Power Total Power" << std::endl;
    for (int nx = 0; nx < N_ROTORS; nx++)
    {
        cout << nx + 1 << " " << round_(Thrust_new[nx], 3) << " " << round_(Induced_Power[nx], 3) << " " << round_(Profile_Power[nx], 3) << " " << round_(Induced_Power[nx] + Profile_Power[nx], 3);
    }
    //return positions_new, RMS_arr, CirVM_array, Thrust_new, Power_new, Blade_velocities, Thrust_contours, Shed_circulation, Trail_circulation, Up_VM_arr;
}


int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    vector<double> wake_age_array = linspace(0, JM * dpsi * 180 / M_PI, JM);

    vector<double> psi_array = linspace(0, 2 * M_PI, N);
    const double location = R_hub;
    vector<double> PTC_inflow_arr(span_sections);
    vector<double> location_arr;

    //cout << "JM: " << JM << ", FAR_FIELD_CUTOFF: " << FAR_FIELD_CUTOFF << std::endl;

    for (int i = 0; i < span_sections; i++)
    {
        double diff = 1.0f;
        double lmd_1 = 1.0f;
        while (diff > 0.001)
        {
            double f = 0.5 * N_BLADES * (1 - radii[i] / R_tip) / lmd_1;
            double FL;
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
    vector<double> phi_array = scalarMul(R_tip, vectorDiv(PTC_inflow_arr, radii));
    std::transform(phi_array.begin(), phi_array.end(), phi_array.begin(), arctan);
    vector<double> alpha_arr = vectorSub(theta, phi_array);
    vector<double> Cd_arr(alpha_arr);
    std::transform(Cd_arr.begin(), Cd_arr.end(), Cd_arr.begin(), square);
    Cd_arr = scalarAdd(Cd0, scalarMul(1.25, Cd_arr));
    vector<double> Cl_arr = scalarMul(a, alpha_arr);
    vector<double> circulation_ = scalarMul(0.5f, elementWiseMul(scalarAdd(u_inf * std::cos(alpha) * std::sin(M_PI / 2), scalarMul(OMEGA, radii)), elementWiseMul(chord, Cl_arr)));
    vector<double> circulation_old = circulation_;

    const double gamma_max = *std::max_element(circulation_.begin(), circulation_.end());
    vector<double> Thrust2_arr = vector<double>(span_sections, 0.0f);
    vector<double> Power2_arr = vector<double>(span_sections, 0.0f);
    vector<double> UP_arr = vector<double>(span_sections, 0.0f);
    const double T_a = 0.5 * rho * BLADES_PER_ROTOR; // constant part in thrust equation
    vector<double> Ut_arr = scalarMul(OMEGA, radii);
    for (int i = 0; i < span_sections; i++)
    {
        const double Ut = OMEGA * radii[i];  // tangential velocity
        const double Up = PTC_inflow_arr[i] * OMEGA * R_tip;  // perpendicular velocity
        const double Cd = Cd0 + 1.25 * pow(theta[i] - atan(Up / Ut), 2);  // drag equation in terms of r
        const double Cl = a * (theta[i] - atan(Up / Ut));  // lift equation in terms of r

        const double dT = 0.5 * rho * N_BLADES * (pow(Up, 2) + pow(Ut, 2)) * chord[i] * (Cl * cos(atan(Up / Ut)) - Cd * sin(atan(Up / Ut)));
        const double dP = 0.5 * rho * N_BLADES * OMEGA * radii[i] * (pow(Up, 2) + pow(Ut, 2)) * chord[i] * (Cd * cos(atan(Up / Ut)) + Cl * sin(atan(Up / Ut)));
        Thrust2_arr[i] = dT;
        Power2_arr[i] = dP;
        UP_arr[i] = Up;
    }
    const double Thrust = trapz(Thrust2_arr, radii) / cos(alpha);
    const double Power = trapz(Power2_arr, radii) / cos(alpha);
    const double CT_A1 = Thrust / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));
    double CP = Power / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));
    vector<double> circulation_A1 = scalarMul(0.5, elementWiseMul(Cl_arr, elementWiseMul(chord, Ut_arr)));

    const double LHS = CT_A1 / N_ROTORS;
    double il = 0.01f; // initial guess for Newton Raphson method

    for (int i = 0; i < 1000; i++)
    {
        const double h = -h1_ig(il, LHS, mu, alpha) / h2_ig(il, LHS, mu, alpha);
        il += h;
    }
    const double L_ig_calculated = il;
    //cout << "L_ig_calculated: " << L_ig_calculated << std::endl;
    //const double LG = L_ig_calculated + L_i;
    vector<double> dOPM_array(N, 0.0f);
    vector<double> dORM_array(N, 0.0f);
    vector<double> dT_blade_array(N, 0.0f);
    vector<double> dOQi_array(N, 0.0f);
    vector<double> dOQo_array(N, 0.0f);
    vector<double> dT_array(N, 0.0f);
    vector<double> dOPi_array(N, 0.0f);
    vector<double> dOPo_array(N, 0.0f);

    const double dbeta_dt = 0;  // TPP as a reference frame
    const double beta_o_d = 0;  // TPP frame of reference # beta   in degrees
    const double beta_1s_d = 0;  // beta_1s kept 0 or as pilot input in degrees
    const double beta_1c_d = 0;  // beta_1c kept 0 pilot input  degrees
    const double theta_1s = (0) * (M_PI / 180);  // theta_1s  in degrees at trim conditio
    const double theta_1c = (0) * (M_PI / 180);  // theta_1c  in degrees at trim conditio
    vector<vector<double>> Thrust_plot_array(N + 1, vector<double>(span_sections, 0.0f));
    vector<vector<double>> circulation(JM, vector<double>(span_sections, 0.0f));
    vector<vector<double>> Up_BEM_array(N, vector<double>(span_sections, 0.0f));

    for (int j = 0; j < N; j++)
    {
        const double beta = (beta_o_d + beta_1c_d * cos(j * dpsi) + beta_1s_d * sin(j * dpsi)) * (M_PI / 180);  // beta
        // arrays to store sectional Lift, induced & profile Torque, induced & profile Power at each radial location
        vector<double> dL_array(span_sections, 0.0f);
        vector<double> dP_array(span_sections, 0.0f);
        vector<double> dR_array(span_sections, 0.0f);
        vector<double> dQi_array(span_sections, 0.0f);
        vector<double> dQo_array(span_sections, 0.0f);
        vector<double> dPi_array(span_sections, 0.0f);
        vector<double> dPo_array(span_sections, 0.0f);

        for (int i = 0; i < span_sections; i++) // loop to solve the dL, dQi, dQo, dPi, dPo, dPM, dRM
        {
            double inr = 0.0f;
            double difference = 1.0f;
            while (difference > 0.001f)
            {
                const double h = -h1_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip) / h2_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip);
                //cout << h1_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip) << " " << h2_i(inr, L_ig_calculated, mu, j * dpsi, radii[i], R_tip) << ", ";
                difference = abs(h);
                inr = inr + h;
            }
            //cout << std::endl;
            const double L_i_calculated = inr; // calculated value of L_i (inflow ratio) at given radial location

            const double theta_equation = (theta[i] + theta_1c * cos(j * dpsi) + theta_1s * sin(j * dpsi));  // theta in tersm of r i.e.theta_twd = E + Fr
            /*cout << "u_inf: " << u_inf << ", alpha: " << alpha << ", L_i_calculated: " << L_i_calculated << ", Vv: " << Vv << ", radii[i]" << radii[i] <<
                ", dbeta_dt:" << dbeta_dt << ", beta: " << beta << ", j*dpsi" << j*dpsi << std::endl;*/
            const double Up = u_inf * sin((alpha)) + (OMEGA * R_tip * (L_i_calculated)+Vv) + radii[i] * dbeta_dt + u_inf * sin(beta) * cos(j * dpsi);
            const double Ut = OMEGA * radii[i] + u_inf * cos(alpha) * sin(j * dpsi);

            Up_BEM_array[j][i] = Up;

            //cout << "rho: " << rho << ", chord[i]: " << chord[i] << ", a" << a << ", theta_equation: " << theta_equation << ", Ut: " << Ut << ", Up: " << Up << std::endl;

            dL_array[i] = (0.5 * rho * chord[i] * a * (theta_equation * pow(Ut, 2) - Ut * Up)); // filling the  corresponding array
            dP_array[i] = dL_array[i] * radii[i] * cos(j * dpsi);  // filling the  corresponding array
            dR_array[i] = dL_array[i] * radii[i] * sin(j * dpsi);  // filling the  corresponding array

            const double dDi = 0.5 * rho * chord[i] * a * (theta_equation - atan(Up / Ut)) * Up * Ut; // Di term in induced torque equation
            dQi_array[i] = dDi * radii[i]; // filling the  corresponding array

            const double dDo = 0.5 * rho * Ut * Ut * chord[i] * (Cd0); // Do term in  profile torque equation
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

    const double THRUST = trapz(dT_array, psi_array) * (N_BLADES / (2 * M_PI)); // thust produced at given azimuthal location

    const double OVERALL_PITCHING_MOMENT = trapz(dOPM_array, psi_array) * (N_BLADES / (2 * M_PI));   // Overall Pitching Moment produced by rotor
    const double OVERALL_ROLLING_MOMENT = trapz(dORM_array, psi_array) * (N_BLADES / (2 * M_PI));  // Overall Rolling Moment produced by rotor
    const double INDUCED_TORQUE = trapz(dOQi_array, psi_array) * (N_BLADES / (2 * M_PI)); // Induced Torque produced by rotor
    const double PROFILE_TORQUE = trapz(dOQo_array, psi_array) * (N_BLADES / (2 * M_PI)); // Profile Torque produced by rotor

    const double TOTAL_TORQUE = INDUCED_TORQUE + PROFILE_TORQUE; // Total Torque produced by rotor

    const double INDUCED_POWER = trapz(dOPi_array, psi_array) * (N_BLADES / (2 * M_PI)); // Induced Power produced by rotor
    const double PROFILE_POWER = trapz(dOPo_array, psi_array) * (N_BLADES / (2 * M_PI)); // Profile Power produced by rotor
    const double TOTAL_POWER = (INDUCED_POWER)+PROFILE_POWER; // Total Power produced by rotor

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

    const double CT = THRUST / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));  // Thrust co efficient 
    CP = TOTAL_POWER / (rho * DISC_AREA * pow(OMEGA * R_tip, 3));  // Thrust co efficient 
    cout << "ThrustA2:" << round(Thrust * 1000.0) / 1000.0 << "N and CT:" << CT << std::endl;
    cout << "PowerA2:" << round(Power) / 1000.0 << "kW and CP:" << CP << std::endl;
    cout << "Thrust:" << round(THRUST * 1000.0) / 1000.0 << "N and CT:" << CT << std::endl;
    cout << "Power:" << round(TOTAL_POWER) / 1000.0 << "kW and CP:" << CP << std::endl;

    //const double THRUST_ = 8786.417;
    //const double CT_ = THRUST / (rho * DISC_AREA * pow(OMEGA * R_tip, 2));
    vector< vector<vector<vector<double>>>> hoverWakerSoln = HoverWake(CT);
    //save4DArray(hoverWakerSoln, "hover_wake_sol.npy");
    saveData(hoverWakerSoln);
    //std::cout << "JM: " << JM;
    //std::cout << "CT: " << CT;
    //std::cout << "dpsi: " << dpsi;

    cout << "   RMS_CHANGE         Thrust(N)       Power(kW)" << std::endl;

    EX2BMethod(CT, circulation);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Total time taken: " << duration.count() << " seconds.";

    return 0;
}
