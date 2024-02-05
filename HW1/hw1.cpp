//
// Created by Vitto R on 2/1/24.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

// Q1: Read Data
tuple<int, vector<vector<double>>, vector<double>> read_txt_file(string& file_name){
    // Read in input file
    ifstream inputFile("/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW1/HW1/"+file_name);

    // Throw error if file was not opened
    if (!inputFile) {
        cerr << "Error opening file." << endl;
    }

    // Set total number of atoms
    int n;
    inputFile >> n;

    // Read in atom identity and xyz coordinates
    vector<vector<double>> xyz_list;
    vector<double> atom_list;
    for (int i = 0; i < n; ++i) {
        double atom, x, y, z;
        inputFile >> atom >> x >> y >> z ;
        if (atom != 79) {
            cerr << "Atom No." << i+1 << ": This atom is not a gold atom!" << endl;
        }
        atom_list.push_back(atom);
        xyz_list.push_back({x, y, z});
    }
    inputFile.close();

    // Echo Input
    for (int i = 0; i < n; ++i) {
        cout << atom_list[i]<<  '(' << xyz_list[i][0] << ',' << xyz_list[i][1] << ',' << xyz_list[i][2] << ')' << endl;
    }
    return make_tuple(n, xyz_list, atom_list);
}

// Q1: Driver function for selective epsilon and sigma values based on atom identity
tuple<double, double, double, double> e_s_select(const vector<double>& atom_list, int i, int j) {
    // Assign constants here
    vector<int>    atom_library    = {79   ,47   };
    vector<double> epsilon_library = {5.29 ,4.56 };
    vector<double> sigma_library   = {2.951,2.955};

    double epsilon_i = epsilon_library[find(atom_library.begin(), atom_library.end(), atom_list[i]) -
                                       atom_library.begin()];
    double epsilon_j = epsilon_library[find(atom_library.begin(), atom_library.end(), atom_list[j]) -
                                       atom_library.begin()];
    double sigma_i = sigma_library[find(atom_library.begin(), atom_library.end(), atom_list[i]) -
                                   atom_library.begin()];
    double sigma_j = sigma_library[find(atom_library.begin(), atom_library.end(), atom_list[j]) -
                                   atom_library.begin()];

    return make_tuple(epsilon_i, epsilon_j, sigma_i, sigma_j);
}

// Q1: Function to select coordinates from xyz_list
tuple<double, double, double, double, double, double> xyz_select(const vector<vector<double>>& xyz_list, int i, int j) {
    double xi = xyz_list[i][0];
    double yi = xyz_list[i][1];
    double zi = xyz_list[i][2];

    double xj = xyz_list[j][0];
    double yj = xyz_list[j][1];
    double zj = xyz_list[j][2];

    return make_tuple(xi, xj, yi, yj, zi, zj);
}

// Q1: Lennard-Jones Potential function
double LJP(double e, double s, double R) {
    return e * (pow(s / R, 12) - 2 * pow(s / R, 6));
}

// Q1: Vector norm calculator function
double Rij(double xi, double xj, double yi, double yj, double zi, double zj) {
    return sqrt(pow(xi - xj, 2) + pow(yi - yj, 2) + pow(zi - zj, 2));
}
// Q1: Calculate LJP Energy
double CalculateLJPEnergy(int n , vector<vector<double>> xyz_list, vector<double> atom_list){
    // Initialize LJP
    double E_LJ = 0.0;

    // Iterate through every possible 2-body interaction
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            // Collect the atomic epsilon and sigma values
            auto [epsilon_i, epsilon_j, sigma_i, sigma_j] = e_s_select(atom_list, i, j);

            // Calculate the weighted two-body epsilon and sigma values
            double epsilon_ij = sqrt(epsilon_i * epsilon_j);
            double sigma_ij = sqrt(sigma_i * sigma_j);

            // Collect the atomic coordinates
            auto [xi, xj, yi, yj, zi, zj] = xyz_select(xyz_list, i, j);

            // Calculate inter-atomic distance
            double R_ij = Rij(xi, xj, yi, yj, zi, zj);
            // Calculate LJ potential of this specific 2-body interaction
            double E_ij = LJP(epsilon_ij, sigma_ij, R_ij);

            // Add all of the LJPs
            E_LJ += E_ij;
        }
    }
    return E_LJ;
}

// Q2: Analytical LJP Force
double AnalytLJP_Force(double e,double s,double R,double xk,double xi){
    return -e*(12*(pow(s,12)/pow(R,13))-12*(pow(s,6)/pow(R,7)))*(xk-xi)/(R);
}

mat CalculateAnalyticalForce(int n , vector<vector<double>> xyz_list, vector<double> atom_list){
    // Initialize Analytical Force
    mat F(3,n,fill::zeros);

    // Iterate through every possible 2-body interaction
    for (int i = 0; i < n; ++i) {
        double F_x_i = 0,F_y_i = 0,F_z_i = 0;
        for (int k = 0; k < n; ++k) {
            if (i != k){
                // Collect the atomic epsilon and sigma values
                auto [epsilon_i, epsilon_k, sigma_i, sigma_k] = e_s_select(atom_list, i, k);

                // Calculate the weighted two-body epsilon and sigma values
                double epsilon_ik = sqrt(epsilon_i * epsilon_k);
                double sigma_ik = sqrt(sigma_i * sigma_k);

                // Collect the atomic coordinates
                auto [xi, xk, yi, yk, zi, zk] = xyz_select(xyz_list, i, k);

                // Calculate inter-atomic distance
                double R_ik = Rij(xi, xk, yi, yk, zi, zk);

                // Calculate LJ potential of this specific 2-body interaction
                double F_x_ik = AnalytLJP_Force(epsilon_ik,sigma_ik,R_ik,xk,xi);
                double F_y_ik = AnalytLJP_Force(epsilon_ik,sigma_ik,R_ik,yk,yi);
                double F_z_ik = AnalytLJP_Force(epsilon_ik,sigma_ik,R_ik,zk,zi);

                F_x_i += F_x_ik;
                F_y_i += F_y_ik;
                F_z_i += F_z_ik;
            }
        }
        F(0,i) = F_x_i;
        F(1,i) = F_y_i;
        F(2,i) = F_z_i;
    }
    return F;
}

// Q2: Lennard-Jones Potential Forward-Difference Force
double ApproxLJP_Force(const string& force_type,const string& dir,
                        double e,double s,
                        double xi,double xj,double yi,
                        double yj,double zi,double zj,double h){
    double Rij1, Rij2, approxF;
    if (force_type == "for"){
        if (dir == "x"){
            Rij1 = Rij(xi+h,xj,yi,yj,zi,zj);
            Rij2 = Rij(xi  ,xj,yi,yj,zi,zj);
        } else if (dir == "y"){
            Rij1 = Rij(xi,xj,yi+h,yj,zi,zj);
            Rij2 = Rij(xi,xj,yi  ,yj,zi,zj);
        } else if (dir == "z"){
            Rij1 = Rij(xi,xj,yi,yj,zi+h,zj);
            Rij2 = Rij(xi,xj,yi,yj,zi  ,zj);
        }
        approxF = -pow(  h,-1)*(LJP(e,s,Rij1)-LJP(e,s,Rij2));
    } else if (force_type == "cen"){
        if (dir == "x"){
            Rij1 = Rij(xi+h,xj,yi,yj,zi,zj);
            Rij2 = Rij(xi-h,xj,yi,yj,zi,zj);
        } else if (dir == "y"){
            Rij1 = Rij(xi,xj,yi+h,yj,zi,zj);
            Rij2 = Rij(xi,xj,yi-h,yj,zi,zj);
        } else if (dir == "z"){
            Rij1 = Rij(xi,xj,yi,yj,zi+h,zj);
            Rij2 = Rij(xi,xj,yi,yj,zi-h,zj);
        }
        approxF = -pow(2*h,-1)*(LJP(e,s,Rij1)-LJP(e,s,Rij2));
    } else {
        cerr << "Check approx force inputs" << endl;
        return 1;
    }
    return approxF;
}

// Q2: Calculate approximate LJ force (forward diff & central diff
tuple<mat,mat> CalculateApproxForce(int n , vector<vector<double>> xyz_list, vector<double> atom_list,double h){
    mat For_F(3,n,fill::zeros);
    mat Cen_F(3,n,fill::zeros);
    for (int i = 0; i < n; ++i) {
        double For_F_x_i = 0, For_F_y_i = 0, For_F_z_i = 0;
        double Cen_F_x_i = 0, Cen_F_y_i = 0, Cen_F_z_i = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                // Collect the atomic epsilon and sigma values
                auto [epsilon_i, epsilon_j, sigma_i, sigma_j] = e_s_select(atom_list, i, j);

                // Calculate the weighted two-body epsilon and sigma values
                double epsilon_ij = sqrt(epsilon_i * epsilon_j);
                double sigma_ij = sqrt(sigma_i * sigma_j);

                // Collect the atomic coordinates
                auto [xi, xj, yi, yj, zi, zj] = xyz_select(xyz_list, i, j);

                // calculate forward and central difference force of this specific 2-body interaction
                double For_F_x_ij = ApproxLJP_Force("for","x",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                double For_F_y_ij = ApproxLJP_Force("for","y",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                double For_F_z_ij = ApproxLJP_Force("for","z",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                double Cen_F_x_ij = ApproxLJP_Force("cen","x",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                double Cen_F_y_ij = ApproxLJP_Force("cen","y",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                double Cen_F_z_ij = ApproxLJP_Force("cen","z",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                // Add contribution of two-body interactions to atom-specific force calculations
                For_F_x_i += For_F_x_ij;
                For_F_y_i += For_F_y_ij;
                For_F_z_i += For_F_z_ij;
                Cen_F_x_i += Cen_F_x_ij;
                Cen_F_y_i += Cen_F_y_ij;
                Cen_F_z_i += Cen_F_z_ij;
            }
        }
        // Add force calculations to matrices
        For_F(0,i) = For_F_x_i;
        For_F(1,i) = For_F_y_i;
        For_F(2,i) = For_F_z_i;
        Cen_F(0,i) = Cen_F_x_i;
        Cen_F(1,i) = Cen_F_y_i;
        Cen_F(2,i) = Cen_F_z_i;
    }
    return make_tuple(For_F,Cen_F);
}

void TruncationErrorCalculations(vector<double> For_error,vector<double> Cen_error,vector<double> h_values){
    // Calculate slope of truncation error vs h
    vector<double> F_slopes;
    vector<double> C_slopes;
    for (int i = 0; i < h_values.size()-1;++i){
        double F_slope = (For_error[i+1]-For_error[i])/(h_values[i+1]-h_values[i]);
        F_slopes.push_back(F_slope);
        double C_slope =(Cen_error[i+1]-Cen_error[i])/(h_values[i+1]-h_values[i]);
        C_slopes.push_back(C_slope);
    }
    double F_avg_slope = reduce(F_slopes.begin(),F_slopes.end())/F_slopes.size();
    double C_avg_slope = reduce(C_slopes.begin(),C_slopes.end())/C_slopes.size();

    // Print error of the forward/central difference approximations for the diff h step sizes
    cout << "Forward difference error for the diff h step sizes: ";
    for (int i = 0; i < For_error.size(); i++)
    {
        cout << For_error[i];
    }
    cout << "with slope " << F_avg_slope;
    cout << endl << "Central difference error for the diff h step sizes: ";
    for (int i = 0; i < Cen_error.size(); i++)
    {
        cout << Cen_error[i];
    }
    cout << "with slope " << C_avg_slope;
}
//Question 3
mat VecVec2Mat(int n, vector<vector<double>> xyz_list){
    mat new_matrix(3,n,fill::zeros);
    for (int i = 0; i < n; ++i){
        for (int j=0; j < 3; ++j){
            new_matrix(j,i)=xyz_list[i][j];
        }
    }
    return new_matrix;
}

vector<vector<double>> Mat2VecVec(int n, mat coord_mat){
    vector<vector<double>> new_vec;
    for (int i = 0; i < n; ++i){
        vector<double> coord = {coord_mat(0,i),coord_mat(1,i),coord_mat(2,i)};
        new_vec.push_back(coord);
    }
    return new_vec;
}

double VecNorm(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
}

mat UnitGrad(int n,mat Gradient){
    for (int i = 0; i < n; ++i){
        Gradient(0,i) *= 1/VecNorm(Gradient(0,i),Gradient(1,i),Gradient(2,i));
        Gradient(1,i) *= 1/VecNorm(Gradient(0,i),Gradient(1,i),Gradient(2,i));
        Gradient(2,i) *= 1/VecNorm(Gradient(0,i),Gradient(1,i),Gradient(2,i));
    }
    return Gradient;
}

double GradNorm(int n,mat Gradient){
    double x_comp = 0;
    double y_comp = 0;
    double z_comp = 0;
    for (int i = 0; i < n; ++i){
        x_comp += abs(Gradient(0,i));
        y_comp += abs(Gradient(1,i));
        z_comp += abs(Gradient(2,i));
    }
    double norm_val = VecNorm(x_comp,y_comp,z_comp);
    return norm_val;
}

int main() {
    // Read in file
    string file_name = "2.txt";
    auto [n, xyz_list, atom_list] = read_txt_file(file_name);

    //Question 1: Calculate + print out LJP energy
    double E_LJ = CalculateLJPEnergy(n,xyz_list,atom_list);
    cout << "E_LJ = " << E_LJ << endl;

    // Question 2: Calculate + print out Analytical LJP Force
    mat F = CalculateAnalyticalForce(n,xyz_list,atom_list);
    F.print("F_LJ_analytical");

    // Question 2: Calculate + print out Analytical LJP Force
    // h, Step Sizes
    vector<double> h_values = {0.1,0.01,0.001,0.0001};

    // Initialize Error
    vector<double> For_error;
    vector<double> Cen_error;

    // Iterate through diff h values & calculate approximate force
    bool run = true;
    if (run){
        for (double h : h_values){
            cout << "Stepsize for finite difference:" << h << endl;
            // Calculate Approximate Force (Forward Difference & Central Difference)
            auto [For_F,Cen_F] = CalculateApproxForce(n , xyz_list, atom_list,h);
            // Print force approximations
            For_F.print("F_LJ_forward_difference");
            Cen_F.print("F_LJ_central_difference");
            For_error.push_back(accu(abs(For_F-F)));
            Cen_error.push_back(accu(abs(Cen_F-F)));
        }
    }


    // Calculate slope of truncation error vs h
    TruncationErrorCalculations(For_error,Cen_error,h_values);

    // Question 3
    cout << "start steepest descent with golden section line search" << endl;
    double E = CalculateLJPEnergy(n,xyz_list,atom_list);
    cout << "Initial energy: " << E << endl;
    double h = 0.0001;
    double s2 = 0.3; //standard step
    double golden = (3-sqrt(5))/2;
    double s1 = 10; // initial step
    double ForceConvergenceThreshold = 0.01;
    double EnergyConvergenceThreshold = 1e-8;
    cout << "Stepsize for central difference is:"<< h <<";Initial stepsize for line search is:"<< s2 <<";Threshold for convergence in force is:" << ForceConvergenceThreshold << endl;
    auto [For_F,Cen_F] = CalculateApproxForce(n , xyz_list, atom_list,h);
    Cen_F.print("Central Difference Force");

    //UnitGrad(n,F).print("Hopefully unit vector");
    cout << "Start steepest descent with golden section line search using central difference force" << endl;
    vector<vector<double>> A = xyz_list;
    int numberIterations;
    while (GradNorm(n,F)> ForceConvergenceThreshold){
        numberIterations += 1;
        cout << "Start golden section search. Interation " << numberIterations << endl;
        mat C = VecVec2Mat(n,A)+UnitGrad(n,F)*s1;
        mat B = VecVec2Mat(n,A)+UnitGrad(n,Cen_F)*s1*golden;
        double s1 = 0;
        while (CalculateLJPEnergy(n,Mat2VecVec(n,B),atom_list) > CalculateLJPEnergy(n,Mat2VecVec(n,C),atom_list) or  CalculateLJPEnergy(n,Mat2VecVec(n,B),atom_list) > CalculateLJPEnergy(n,A,atom_list)){
            s1 += s2;
            mat C = VecVec2Mat(n,A)+UnitGrad(n,F)*s1;
            mat B = VecVec2Mat(n,A)+UnitGrad(n,Cen_F)*s1*golden;
        }
        mat X1 = VecVec2Mat(n,A) + UnitGrad(n,F)*s1*golden*golden;
        mat X2 = VecVec2Mat(n,A) + UnitGrad(n,F)*s1*golden*(1+golden);
        mat D;
        if (CalculateLJPEnergy(n,Mat2VecVec(n,X1),atom_list) <= CalculateLJPEnergy(n,Mat2VecVec(n,X2),atom_list)){
            mat D = X1;
        } else {
            mat D = X2;
        }
        D.print("New point");
        // Exit loop if convergence reached
        if (abs(CalculateLJPEnergy(n,Mat2VecVec(n,D),atom_list)-CalculateLJPEnergy(n,A,atom_list)) <= EnergyConvergenceThreshold){
            break;
        }
        vector<vector<double>> A = Mat2VecVec(n,D);
        F = CalculateAnalyticalForce(n,A,atom_list);
        auto [For_F,Cen_F] = CalculateApproxForce(n , xyz_list, atom_list,h);
        Cen_F.print("Central Difference Force");

    }
    double FinalE = CalculateLJPEnergy(n,A,atom_list);
    cout << "Final energy: " << FinalE << endl;
    cout << "Optimized structure: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << atom_list[i]<<  '(' << A[i][0] << ',' << A[i][1] << ',' << A[i][2] << ')' << endl;
    }
}