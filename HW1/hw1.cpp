//
// Created by Vitto Resnick on 2/1/24.
//
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

// Q1: Read In Data from txt file
tuple<int, vector<vector<double>>, vector<double>> read_txt_file(string& file_name){
    // Read in input file
    ifstream inputFile(file_name);

    // Throw error if file was not opened correctly
    if (!inputFile) {
        cerr << "Error opening file." << endl;
    }

    int n;          // Initialize total number of atoms = n
    inputFile >> n; // Set total number of atoms = n


    vector<vector<double>> xyz_list;       // Initialize list for atoms' xyz coordinates
    vector<double> atom_list;              // Initialize list for atoms' atomic numbers
    // Read in atom identity and xyz coordinates
    for (int i = 0; i < n; ++i) {          // Iterate through every atom
        double atom, x, y, z;              // Initialize atomic number/atom identity and xyz coordinates
        inputFile >> atom >> x >> y >> z ; // Set atomic number/atom identity and xyz coordinates
        if (atom != 79) {                  // If a given atom is not gold, throw an error
            cerr << "Atom No." << i+1 << ": This atom is not a gold atom!" << endl;
        }
        atom_list.push_back(atom);         // Append this atom's atomic number/atom identity to list
        xyz_list.push_back({x, y, z});     // Append this atom's xyz coordinates to list
    }
    inputFile.close();                     // Close the txt file

    // Echo Input
    for (int i = 0; i < n; ++i) { // Iterate through every atom, and output atomic number and xyz coordinates
        cout << atom_list[i]<<  '(' << xyz_list[i][0] << ',' << xyz_list[i][1] << ',' << xyz_list[i][2] << ')' << endl;
    }
    return make_tuple(n, xyz_list, atom_list); // Output # of atoms, coordinates, and identities
}

// Q1: Driver function for selective epsilon and sigma values based on atom identity
tuple<double, double, double, double> e_s_select(const vector<double>& atom_list, int i, int j) {
    // Assign constants here
    vector<int>    atom_library    = {79   ,47   };
    vector<double> epsilon_library = {5.29 ,4.56 };
    vector<double> sigma_library   = {2.951,2.955};

    // Select epsilon i and epsilon j
    double epsilon_i = epsilon_library[find(atom_library.begin(), atom_library.end(), atom_list[i]) -
                                       atom_library.begin()];
    double epsilon_j = epsilon_library[find(atom_library.begin(), atom_library.end(), atom_list[j]) -
                                       atom_library.begin()];
    // Select sigma i and sigma j
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

// Q2: 1D Analytical LJP Force between 2 atoms at xi and xk, distance R away from each other
double AnalytLJP_Force(double e,double s,double R,double xk,double xi){
    return -e*(12*(pow(s,12)/pow(R,13))-12*(pow(s,6)/pow(R,7)))*(xk-xi)/(R);
}

// Q2: Calculate Analytical LJP Force for n atoms at positions xyz_list with identity atom_list
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

// Q2: Lennard-Jones Potential 3D Forward-Difference or Central-Difference Force
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

// Q2: Calculate approximate 3D LJ force (forward diff & central diff)
// for n atoms at positions xyz_list with identity atom_list
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

// Q2: Calcuate slope of truncation error vs h values for forward and central differences
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
// Q3: Convert 2D vector vector<vector<double>> to matrix
mat VecVec2Mat(int n, vector<vector<double>> xyz_list){
    mat new_matrix(3,n,fill::zeros);
    for (int i = 0; i < n; ++i){
        for (int j=0; j < 3; ++j){
            new_matrix(j,i)=xyz_list[i][j];
        }
    }
    return new_matrix;
}

// Q3: Convert matrix to 2D vector vector<vector<double>>
vector<vector<double>> Mat2VecVec(int n, mat coord_mat){
    vector<vector<double>> new_vec;
    for (int i = 0; i < n; ++i){
        vector<double> coord = {coord_mat(0,i),coord_mat(1,i),coord_mat(2,i)};
        new_vec.push_back(coord);
    }
    return new_vec;
}

// Q3: Vectorize Gradient -> Calculate Norm value of Gradient
double GradNorm(int n,mat Gradient){
    return norm(vectorise(Gradient));
}

// Q3: Calculate Unit Vector (matrix) of Gradient or Normalize Gradient
mat UnitGrad(int n,mat Gradient){
    return Gradient/GradNorm(n,Gradient);
}

// Q3: LJP Calculations For Matrices inputs of xyz_list coords
double LJ(int n , mat M, vector<double> ats){
    // Initialize LJP
    double E_LJ = 0.0;

    // Iterate through every possible 2-body interaction
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            // Collect the atomic epsilon and sigma values
            auto [epsilon_i, epsilon_j, sigma_i, sigma_j] = e_s_select(ats, i, j);

            // Calculate the weighted two-body epsilon and sigma values
            double epsilon_ij = sqrt(epsilon_i * epsilon_j);
            double sigma_ij = sqrt(sigma_i * sigma_j);

            // Collect the atomic coordinates
            vector<vector<double>> xyz_list = Mat2VecVec(n,M);
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
mat    aF(int n , mat M, vector<double> ats){
    // Initialize Analytical Force
    mat F(3,n,fill::zeros);

    // Iterate through every possible 2-body interaction
    for (int i = 0; i < n; ++i) {
        double F_x_i = 0,F_y_i = 0,F_z_i = 0;
        for (int k = 0; k < n; ++k) {
            if (i != k){
                // Collect the atomic epsilon and sigma values
                auto [epsilon_i, epsilon_k, sigma_i, sigma_k] = e_s_select(ats, i, k);

                // Calculate the weighted two-body epsilon and sigma values
                double epsilon_ik = sqrt(epsilon_i * epsilon_k);
                double sigma_ik = sqrt(sigma_i * sigma_k);

                // Collect the atomic coordinates
                vector<vector<double>> xyz_list = Mat2VecVec(n,M);
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
mat    cF(int n , mat M, vector<double> ats,double h){
    mat Cen_F(3,n,fill::zeros);
    for (int i = 0; i < n; ++i) {
        double Cen_F_x_i = 0, Cen_F_y_i = 0, Cen_F_z_i = 0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                // Collect the atomic epsilon and sigma values
                auto [epsilon_i, epsilon_j, sigma_i, sigma_j] = e_s_select(ats, i, j);

                // Calculate the weighted two-body epsilon and sigma values
                double epsilon_ij = sqrt(epsilon_i * epsilon_j);
                double sigma_ij = sqrt(sigma_i * sigma_j);

                // Collect the atomic coordinates
                vector<vector<double>> xyz_list = Mat2VecVec(n,M);
                auto [xi, xj, yi, yj, zi, zj] = xyz_select(xyz_list, i, j);

                // calculate forward and central difference force of this specific 2-body interaction
                double Cen_F_x_ij = ApproxLJP_Force("cen","x",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                double Cen_F_y_ij = ApproxLJP_Force("cen","y",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                double Cen_F_z_ij = ApproxLJP_Force("cen","z",epsilon_ij,sigma_ij,xi,xj,yi,yj,zi,zj,h);
                // Add contribution of two-body interactions to atom-specific force calculations
                Cen_F_x_i += Cen_F_x_ij;
                Cen_F_y_i += Cen_F_y_ij;
                Cen_F_z_i += Cen_F_z_ij;
            }
        }
        // Add force calculations to matrices
        Cen_F(0,i) = Cen_F_x_i;
        Cen_F(1,i) = Cen_F_y_i;
        Cen_F(2,i) = Cen_F_z_i;
    }
    return Cen_F;
}

// Q3: Golden Section Line Search
double GoldenSectionSearch(double EnergyConvergenceThreshold,
                           int n,mat F, mat A,
                           double a, double b, double c,vector<double> ats){
    double golden = (3-sqrt(5))/2;
    mat B = A + UnitGrad(n,F)*b;
    mat C = A + UnitGrad(n,F)*c;
    double LJ_A = LJ(n,A,ats);
    double LJ_B = LJ(n,B,ats);
    double LJ_C = LJ(n,C,ats);

    if (abs(LJ_A-LJ_B) <= EnergyConvergenceThreshold or abs(LJ_B-LJ_C) <= EnergyConvergenceThreshold){
        return b;
    } else if (abs(LJ_A-LJ_C) <= EnergyConvergenceThreshold) {
        return a;
    } else {
        double x = a + golden*(c-a);
        mat X = A + UnitGrad(n,F)*x;
        double LJ_X = LJ(n,X,ats);
                     // if f(x) < f(b)
        if (LJ_X < LJ_B){
                     // if f(x) < f(b) and a < b < x, new (a,b,c)=(b,x,c)
            if (x > b){
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,b,x,c,ats);
            } else { // if f(x) < f(b) and a < x < b, new (a,b,c)=(a,x,b)
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,a,x,b,ats);
            }
        } else {     // if f(x) > f(b)
                     // if f(x) > f(b) and a < b < x, new (a,b,c)=(a,b,x)
            if (x > b){
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,a,b,x,ats);
            } else { // if f(x) > f(b) and a < x < b, new (a,b,c)=(x,b,c)
                return GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,x,b,c,ats);
            }
        }
    }
}

// Q3:　1D Steepest Descent
mat steepest_descent(mat A, int& i, double ForceConvergenceThreshold,double EnergyConvergenceThreshold,double h,int n, double s, vector<double> ats){
    i += 1; // This is the i-th iteration

    // 1. Input initial point guess (Point A)
    double a = 0;
    double LJ_A = LJ(n,A,ats);

    // 2.1. Compute the gradient, take a step in the opposite direction, check for convergence

    // Compute the gradient
    bool useAnalytical = false;
    mat F;
    if (useAnalytical){
        F = aF(n,A,ats);
    } else{
        F = cF(n,A,ats,h);
    }

    // Print New point & grad
    if (!(i==0)){
        cout << "Start golden section search" << endl;
        A.print("new_point");
        cout << "current energy: " << LJ(n,A,ats) << endl;
        mat Cen_F = cF(n,A,ats,h);
        Cen_F.print("Central Difference Force");
    }


    // check for convergence
    if (GradNorm(n,F) <= ForceConvergenceThreshold){
        return A;
    } else {
        // take a step in the opposite direction
        double b = s;
        mat B = A + UnitGrad(n,F)*b;
        double LJ_B = LJ(n,B,ats);

        // 2.2. Search for Point B in the opposite direction of the gradient (E(B) < E(A))
        while (!(LJ_B<LJ_A)){
            b /= 2;
            B = A + UnitGrad(n,F)*b;
            LJ_B = LJ(n,B,ats);
        }

        // 2.3. Search for Point C in the opposite direction of the gradient (E(C) > E(B))
        double c = 10;
        mat C = A + UnitGrad(n,F)*c;
        double LJ_C = LJ(n,C,ats);
        while (!(LJ_C > LJ_B)){
            c *= 2;
            C = A + UnitGrad(n,F)*c;
            LJ_C = LJ(n,C,ats);
        }

        // 2.5. If can find both B and C, use golden section search
        // for the minimum in the opposite direction of the gradient,
        double step_size = GoldenSectionSearch(EnergyConvergenceThreshold,n,F,A,a,b,c,ats);
        A = A + UnitGrad(n,F)*step_size;

        // return to step 2.1
        return steepest_descent(A,i,ForceConvergenceThreshold,EnergyConvergenceThreshold,h,n,s,ats);
    }
}

int main(int argc, char* argv[]) {
    // Read in file
    string file_name = "/Users/vittor/Documents/CLASSES/SPRING 2024/CHEM_179_HW1/HW1/2.txt";
    auto [n, xyz_list, atom_list] = read_txt_file(file_name);

    // Question 1: Calculate + print out LJP energy
    cout << "Question 1:" << endl;
    double E_LJ = CalculateLJPEnergy(n,xyz_list,atom_list);
    cout << "E_LJ = " << E_LJ << endl;

    // Question 2: Calculate + print out Analytical LJP Force
    cout << "\n Question 2:" << endl;
    mat F = CalculateAnalyticalForce(n,xyz_list,atom_list);
    F.print("F_LJ_analytical");

    // Question 2: Calculate + print out Analytical LJP Force
    vector<double> h_values = {0.1,0.01,0.001,0.0001}; // h, Step Sizes

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
        // Calculate slope of truncation error vs h
        TruncationErrorCalculations(For_error,Cen_error,h_values);
    }

    // Question 3
    cout << endl;
    cout << endl;
    cout << "Question 3:" << endl;
    cout << "start steepest descent with golden section line search" << endl;
    double E = CalculateLJPEnergy(n,xyz_list,atom_list);
    cout << "Initial energy: " << E << endl;

    double h = 0.0001;
    double s = 0.3; // initial step
    double ForceConvergenceThreshold = 0.01;
    double EnergyConvergenceThreshold = 1e-8;
    cout << "Stepsize for central difference is:"<< h <<";Initial stepsize for line search is:"<< s <<";Threshold for convergence in force is:" << ForceConvergenceThreshold << endl;

    // 1. Input initial point guess (Point A)
    mat A = VecVec2Mat(n,xyz_list);
    vector<double> ats = atom_list;
    mat Cen_F = cF(n,A,ats,h);
    Cen_F.print("Central Difference Force");

    cout << "Start steepest descent with golden section line search using central difference force" << endl;
    int it = 0;
    A = steepest_descent(A,it,ForceConvergenceThreshold,EnergyConvergenceThreshold,h,n, s, ats);

    cout << "Total iterations: " << it << endl;
    double FinalE = LJ(n,A,ats);
    A.print("This is A");
    cout << "Final energy: " << FinalE << endl;
    cout << "Optimized structure: " << endl;
    for (int i = 0; i < n; ++i) {
        cout << atom_list[i] <<  '(' << A(0,i) << ',' << A(1,i) << ',' << A(2,i) << ')' << endl;
    }

}