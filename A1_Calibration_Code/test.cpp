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


int main()
{
    std::string myPath = R"(D:\TU Delft\GEO1016 Photogrametry\Assignment_n_1\Photogrametry_assg_1\A1_Calibration_Code\resources\data\test_data_1(6_points)-normal.txt)";
    std::vector<double> flat_data = read_file(myPath);

    coord coord_result  = flat_to_matrix(flat_data);

    std::vector<std::vector<double>> point3D = coord_result.pw;
    std::cout << "Number of points " << point3D.size();

    std::vector<std::vector<double>> point2D = coord_result.pc;

}
