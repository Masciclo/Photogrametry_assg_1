#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Calibration/matrix.h"
#include "Calibration/matrix_algo.h"
#include "Calibration/vector.h"

using namespace easy3d;



// Return double to handle precise coordinates
std::vector<double> read_file(std::string path) {
    std::ifstream myfile(path);
    std::vector<double> all_data;
    double val; // Changed to double

    if (myfile.is_open()) {
        while (myfile >> val) {
            all_data.push_back(val);
        }
        myfile.close();
    } else {
        std::cerr << "Error: Could not open file at: " << path << std::endl;
    }
    return all_data;
}

// Define coord structure to output two vectors in flat_to matrix function
struct coord {
    std::vector<std::vector<double>> pw; // Changed to double for precision
    std::vector<std::vector<double>> pc;
};

coord flat_to_matrix(std::vector<double> flat_data) {
    // 1. Catch the data as doubles

    coord result; // This specify the result format of the data.

    // 2. Loop through the flat vector 5 elements at a time
    for (size_t i = 0; i + 4 < flat_data.size(); i += 5) {
        // Create 3D point {X, Y, Z}
        std::vector<double> point3D = { flat_data[i], flat_data[i+1], flat_data[i+2], 1 };
        result.pw.push_back(point3D);

        // Create 2D point {u, v}
        std::vector<double> point2D = { flat_data[i+3], flat_data[i+4] };
        result.pc.push_back(point2D);
    }

    return result;
}

Matrix vector_to_matrix(std::vector<std::vector<double>> pw, std::vector<std::vector<double>> pc) {
    int n_points = pw.size();

    // The DLT Design Matrix 'L' usually has 2*n rows and 11 columns
    Matrix L(2 * n_points, 12, 0.0);

    for (int i = 0; i < n_points; i++) {
        // 1. Extract our point values for readability
        double X = pw[i][0];
        double Y = pw[i][1];
        double Z = pw[i][2];
        double u = pc[i][0];
        double v = pc[i][1];

        // 2. Build the first row (for the u-coordinate)
        // Structure: [ X, Y, Z, 1, 0, 0, 0, 0, -u*X, -u*Y, -u*Z ]
        L.set_row(2 * i, {X, Y, Z, 1.0, 0.0, 0.0, 0.0, 0.0, -u*X, -u*Y, -u*Z, -u});

        // 3. Build the second row (for the v-coordinate)
        // Structure: [ 0, 0, 0, 0, X, Y, Z, 1, -v*X, -v*Y, -v*Z ]
        L.set_row(2 * i + 1, {0.0, 0.0, 0.0, 0.0, X, Y, Z, 1.0, -v*X, -v*Y, -v*Z, -u});
    }

    return L;
}




int main()
{
    std::string myPath = R"(D:\TU Delft\GEO1016 Photogrametry\Photogrametry_assg_1\A1_Calibration_Code\resources\data\test_data_1(6_points)-normal.txt)";
    std::vector<double> flat_data = read_file(myPath);

    coord coord_result  = flat_to_matrix(flat_data);

    std::vector<std::vector<double>> point3D = coord_result.pw;
    std::cout << "Number of points " << point3D.size() << std::endl;

    std::vector<std::vector<double>> point2D = coord_result.pc;

    Matrix m = vector_to_matrix(coord_result.pw, coord_result.pc);
    std::cout << m << std::endl;

    // P is (2n x 12)
    const int r = m.rows();
    const int c = m.cols(); // should be 12

    Matrix U(r, r, 0.0);
    Matrix S(r, c, 0.0);
    Matrix V(c, c, 0.0);

    svd_decompose(m, U, S, V);

}
