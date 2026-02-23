#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

struct coord {
    std::vector<std::vector<double>> pw; // Changed to double for precision
    std::vector<std::vector<double>> pc;
};

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



coord flat_to_matrix(std::vector<double> flat_data) {
    // 1. Catch the data as doubles

    coord result; // This specify the result format of the data.

    // 2. Loop through the flat vector 5 elements at a time
    for (size_t i = 0; i + 4 < flat_data.size(); i += 5) {
        // Create 3D point {X, Y, Z}
        std::vector<double> point3D = { flat_data[i], flat_data[i+1], flat_data[i+2] };
        result.pw.push_back(point3D);

        // Create 2D point {u, v}
        std::vector<double> point2D = { flat_data[i+3], flat_data[i+4] };
        result.pc.push_back(point2D);
    }

    return result;
}

    // P is (2n x 12)
    const int r = P.rows();
    const int c = P.cols(); // should be 12

    Matrix U(r, r, 0.0);
    Matrix S(r, c, 0.0);
    Matrix V(c, c, 0.0);

    svd_decompose(P, U, S, V);

    // Solution is last column of V  (smallest singular value)
    Vector m = V.get_column(V.cols() - 1);   // length 12

    // Reshape m -> 3x4 projection matrix M
    M(0,0) = m[0];   M(0,1) = m[1];   M(0,2) = m[2];   M(0,3) = m[3];
    M(1,0) = m[4];   M(1,1) = m[5];   M(1,2) = m[6];   M(1,3) = m[7];
    M(2,0) = m[8];   M(2,1) = m[9];   M(2,2) = m[10];  M(2,3) = m[11];

    //for loop on i: multiplying M times world coordinate coord.pw should give image coordinate coord.pc
    // -> coord-pw - coord.pc ≈≈ 0±0.1
    //preliminary code
    double tol = 0.1;     // ±0.1 requirement
    int ok = 0;

    for (size_t i = 0; i < coords.size(); ++i) {
        double X = coords[i].pw.x;
        double Y = coords[i].pw.y;
        double Z = coords[i].pw.z;

        double u_obs = coords[i].pc.x;
        double v_obs = coords[i].pc.y;

        // Multiply M * [X Y Z 1]^T -> (x, y, w)
        double x = M(0,0)*X + M(0,1)*Y + M(0,2)*Z + M(0,3)*1.0;
        double y = M(1,0)*X + M(1,1)*Y + M(1,2)*Z + M(1,3)*1.0;
        double w = M(2,0)*X + M(2,1)*Y + M(2,2)*Z + M(2,3)*1.0;

        // avoid divide-by-zero
        if (std::abs(w) < 1e-12) {
            std::cout << "Point " << i << ": w ~ 0, cannot project\n";
            continue;
        }

        // Normalize to get pixel coords
        double u_pred = x / w;
        double v_pred = y / w;

        // Compare
        double du = u_pred - u_obs;
        double dv = v_pred - v_obs;

        bool good = (std::abs(du) <= tol) && (std::abs(dv) <= tol);
        if (good) ok++;

        std::cout << "Point " << i
                  << "  obs=(" << u_obs << "," << v_obs << ")"
                  << "  pred=(" << u_pred << "," << v_pred << ")"
                  << "  err=(" << du << "," << dv << ")"
                  << (good ? "  OK\n" : "  BAD\n");
    }

    std::cout << "OK: " << ok << " / " << coords.size() << "\n";

int main()
{
    std::string myPath = R"(D:\TU Delft\GEO1016 Photogrametry\Assignment_n_1\Photogrametry_assg_1\A1_Calibration_Code\resources\data\test_data_1(6_points)-normal.txt)";
    std::vector<double> flat_data = read_file(myPath);

    coord coord_result  = flat_to_matrix(flat_data);

    std::vector<std::vector<double>> point3D = coord_result.pw;
    std::cout << "Number of points " << point3D.size();

    std::vector<std::vector<double>> point2D = coord_result.pc;

}
