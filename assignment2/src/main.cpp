#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/bounding_box.h>
#include <igl/octree.h>
#include <igl/knn.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Normals evaluated via PCA method, #P x3
Eigen::MatrixXd NP;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
int resolution = 20;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

/*** Additional variables for spatial indexing ***/

// Number of buckets per dimension for spatial indexing through a uniform grid
const int buckets_per_dim = 5;

// Spatial Index using a uniform grid: 4D Vector, 3 Dimensions for indexing, last 4th Dimension for creating a bucket of points.
std::vector<int> spatialIndex[buckets_per_dim][buckets_per_dim][buckets_per_dim];

// BB to be set during preprocessSpatialIndex
Eigen::ArrayXXd minP, maxP;

// Grid Padding multiplier
double padding_mult = 0.01;

// "Bool" for alignment option
int alignP = 0;

// Neighbourhood search param
int k_neighbourhood = 5;

// Data path
string data_path;

// Functions
void createGrid();
void evaluateImplicitFunc();
void evaluateImplicitFunc_PolygonSoup();
void getLines();
void pcaNormal();
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid()
{
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;

    bb_min = constrained_points.colwise().minCoeff().array();
    bb_max = constrained_points.colwise().maxCoeff().array();

    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Enlarge Grid slightly

    // Add padding based on BB dim
    bb_min -= dim*padding_mult;
    bb_max += dim*padding_mult;

    // New bounding box dimensions
    dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines()
{
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1)
                {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1)
                {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1)
                {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

// Find the closest points index to p using a naive loop
int closestPointLoop(Eigen::RowVectorXd p){
    // Calculate distances to p for points in P
    Eigen::VectorXd distances (P.rows());
    for(int i = 0; i < P.rows(); i++){
        distances(i) = (P.row(i) - p).norm();
    }
    // Find min index
    double min = std::numeric_limits<double>::max();
    int index = -1;
    for(int i = 0; i < distances.size(); i++){
        if(distances(i) < min){
            min = distances(i);
            index = i;
        }
    }
    return index;
}

// Find the point indices within distance h to p using a naive loop
std::vector<int> withinDistanceLoop(Eigen::RowVectorXd p, double h){
    // Calculate distances to p for points in P
    Eigen::VectorXd distances (P.rows());
    for(int i = 0; i < P.rows(); i++){
        distances(i) = (P.row(i) - p).norm();
    }
    std::vector<int> indices;
    // Save distances(i) <= h
    for(int i = 0; i < distances.size(); i++){
        if(distances(i) <= h){
            indices.push_back(i);
        }
    }
    return indices;
}

// Find which bucket point p belongs to (within the limits of the bounding box)
Eigen::ArrayXi getSpatialIndex(Eigen:: RowVectorXd p){
    return floor((p.array() - minP) / ((maxP - minP) / buckets_per_dim)).cast<int>();
}

// Get bounding box and split it into buckets filling spatialIndex, padding to ensure imprecision errors do not occur
void preprocessSpatialIndex(){
    // Clear for possible repeat calls
    for(int i = 0; i < buckets_per_dim; i++){
        for(int j = 0; j < buckets_per_dim; j++){
            for(int k = 0; k < buckets_per_dim; k++){
                spatialIndex[i][j][k].clear();
            }
        }
    }

    double epsilon = 1e-7;
    minP = P.colwise().minCoeff().array()-epsilon;
    maxP = P.colwise().maxCoeff().array()+epsilon;

    for(int i = 0; i < P.rows(); i++){
        Eigen::ArrayXi si = getSpatialIndex(P.row(i));
        spatialIndex[si[0]][si[1]][si[2]].push_back(i);
    }    
}

// if in range [0, buckets_per_dim)
bool inValidBucketInterval(Eigen::ArrayXi si){
    return (0 <= si[0]) && (si[0] < buckets_per_dim) && (0 <= si[1]) && (si[1] < buckets_per_dim) && (0 <= si[2]) && (si[2] < buckets_per_dim);
}

// Find the closest points index to p using the spatial index (slow for invalid search points that are outside of BB due to searchRadius having to increase over buckets_per_dim)
int closestPointSI(Eigen::RowVectorXd p){
    int searchRadius = 0;
    Eigen::ArrayXi si = getSpatialIndex(p);
    std::vector<int> searchspace;
    while(true){
        for(int xi = -searchRadius; xi <= searchRadius; xi++){
            for(int yi = -searchRadius; yi <= searchRadius; yi++){
                for(int zi = -searchRadius; zi <= searchRadius; zi++){
                    // Bucket in valid interval [0, buckets_per_dim]^3, iterate over the newly added buckets only
                    if((abs(xi) == searchRadius || abs(yi) == searchRadius || abs(zi) == searchRadius)
                        && inValidBucketInterval(Eigen::Array3i(si[0]+xi, si[1]+yi, si[2]+zi))){
                        std::vector<int> bucket = spatialIndex[si[0]+xi][si[1]+yi][si[2]+zi];
                        if(!bucket.empty()){
                            // add bucket to search space
                            searchspace.insert(searchspace.end(), bucket.begin(), bucket.end());
                        }
                    }
                }
            }
        }
        // If nothing has been found, increase search radius
        if(searchspace.empty()){
            searchRadius += 1;
        }
        else{
            break;
        }
    }
    // Calculate distances to p for points in searchspace
    Eigen::VectorXd distances (searchspace.size());
    for(int i = 0; i < searchspace.size(); i++){
        distances(i) = (P.row(searchspace[i]) - p).norm();
    }
    // Find min index
    double min = std::numeric_limits<double>::max();
    int index = -1;
    for(int i = 0; i < distances.size(); i++){
        if(distances(i) < min){
            min = distances(i);
            index = i;
        }
    }
    return searchspace[index];
}

// Find the point indices within distance h to p using the spatial index
std::vector<int> withinDistanceSI(Eigen::RowVectorXd p, double h){
    int searchRadiusX = (int)ceil(h/((maxP(0) - minP(0))/buckets_per_dim));
    int searchRadiusY = (int)ceil(h/((maxP(1) - minP(1))/buckets_per_dim));
    int searchRadiusZ = (int)ceil(h/((maxP(2) - minP(2))/buckets_per_dim));
    // cout << searchRadiusX << " " << searchRadiusY << " " << searchRadiusZ << endl;
    Eigen::ArrayXi si = getSpatialIndex(p);
    std::vector<int> searchspace;
    for(int xi = -searchRadiusX; xi <= searchRadiusX; xi++){
        for(int yi = -searchRadiusY; yi <= searchRadiusY; yi++){
            for(int zi = -searchRadiusZ; zi <= searchRadiusZ; zi++){
                if(inValidBucketInterval(Eigen::Array3i(si[0]+xi, si[1]+yi, si[2]+zi))){
                    // add bucket to search space
                    std::vector<int> bucket = spatialIndex[si[0]+xi][si[1]+yi][si[2]+zi];
                    if(!bucket.empty()){
                        //cout << si[0]+xi << " " << si[1]+yi << " " << si[2]+zi << endl;
                        searchspace.insert(searchspace.end(), bucket.begin(), bucket.end());
                    }
                }
            }
        }
    }
    // Calculate distances to p for points in searchspace
    Eigen::VectorXd distances (searchspace.size());
    for(int i = 0; i < searchspace.size(); i++){
        distances(i) = (P.row(searchspace[i]) - p).norm();
    }
    std::vector<int> indices;
    // Save distances(i) <= h
    for(int i = 0; i < distances.size(); i++){
        if(distances(i) <= h){
            indices.push_back(searchspace[i]);
        }
    }
    return indices;
}

void calculateConstraints(Eigen::MatrixXd &constrained_points_e, Eigen::VectorXd &constrained_values_e, double e_input){
        // Calculate constraint points & values for epsilon
        for(int i = 0; i < P.rows(); i++){
            Eigen::VectorXd p_e;
            double e = e_input * 2;
            do{
                e = e/2;
                p_e = P.row(i) + e * N.row(i); 
            }
            while(closestPointSI(p_e) != i);
            constrained_points_e.row(i) = p_e;
            constrained_values_e(i) = e;
        }
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc_Sphere()
{
    // Sphere center
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    Eigen::RowVector3d center = 0.5 * (bb_min + bb_max);

    double radius = 0.5 * (bb_max - bb_min).minCoeff();

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; x < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);

                // Value at (x,y,z) = implicit function for the sphere
                grid_values[index] = (grid_points.row(index) - center).norm() - radius;
            }
        }
    }
}

// returns the polynomial basis function for point p_i based on polyDegree
Eigen::RowVectorXd polynomialBasis(Eigen::RowVectorXd p_i){
    double x = p_i(0);
    double y = p_i(1);
    double z = p_i(2);
    if(polyDegree == 2){
        Eigen::RowVectorXd bf(10);
        bf << 1, x, y, z, x*y, x*z, y*z, x*x, y*y, z*z;
        return bf;
    }
    // For invalid polyDegree return as if polyDegree == 1.
    else if(polyDegree == 1){
        Eigen::RowVectorXd bf(4);
        bf << 1, x, y, z;
        return bf;
    }
    else{
        Eigen::RowVectorXd bf(1);
        bf << 1;
        return bf;
    }
}

int polynomialBasisLength(){
    if(polyDegree == 2){
        return 10;
    }
    else if(polyDegree == 1){
        return 4;
    }
    // For invalid polyDegree return as if polyDegree == 0.
    else{
        return 1;
    }
}

// Well defined within interval [0, wendlandRadius] for given distance d = (x - p_i).norm()
double wendlandFunction(Eigen::RowVectorXd x, Eigen::RowVectorXd p_i){
    double d = (x-p_i).norm();
    return pow((1-d/wendlandRadius), 4) * (4*d/wendlandRadius+1);
}

// Called upon when "key == '3'"
void evaluateImplicitFunc()
{
    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate signed distance function at each gridpoint via MLS


    cout << "Evaluating with polyDegree=" << polyDegree << ", wendlandRadius=" << wendlandRadius << ", resolution=" << resolution << endl; 


    // Iterate over grid points computing weights for each point and performing one step of MLS
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                int grid_index = x + resolution * (y + resolution * z);
                Eigen::RowVector3d gridPointCoords = grid_points.row(grid_index);

                // Get neighbours in P
                std::vector<int> neighbours = withinDistanceSI(gridPointCoords, wendlandRadius);

                //cout << neighbours.size();

                // If within wendlandRadius empty, then set to large positive number
                if(neighbours.empty()){
                    grid_values[grid_index] = DBL_MAX;
                }
                // Otherwise find weights within neighbourhood and perform MLS step
                else{
                    // Add additional neighbours in contrained_points
                    int num_p_neighbours = neighbours.size();
                    // Add +e neighbours
                    for(int i = 0; i < num_p_neighbours; i++){
                        if((constrained_points.row(neighbours[i] + P.rows()) - gridPointCoords).norm() <= wendlandRadius){
                            neighbours.push_back(neighbours[i] + P.rows());
                        }
                    }
                    // Add -e neighbours
                    for(int i = 0; i < num_p_neighbours; i++){
                        if((constrained_points.row(neighbours[i] + 2*P.rows()) - gridPointCoords).norm() <= wendlandRadius){
                            neighbours.push_back(neighbours[i] + 2*P.rows());
                        }
                    }

                    // Create the polynomial basis functions in matrix form
                    Eigen::MatrixXd B(neighbours.size(), polynomialBasisLength());
                    for(int i = 0; i < B.rows(); i++){
                        B.row(i) = polynomialBasis(constrained_points.row(neighbours[i]));
                    }

                    // Set W
                    Eigen::VectorXd weights = Eigen::VectorXd::Zero(neighbours.size());
                    for(int i = 0; i < neighbours.size(); i++){
                        weights(i) = wendlandFunction(gridPointCoords, constrained_points.row(neighbours[i]));
                    }
                    Eigen::MatrixXd W = weights.asDiagonal();

                    // Set fb
                    Eigen::VectorXd fb(neighbours.size());
                    for(int i = 0; i < neighbours.size(); i++){
                        fb(i) = constrained_values(neighbours[i]);
                    }

                    //cout << "B dims: "<< B.rows() << " " << B.cols() << endl;
                    //cout << "W dims: "<< W.rows() << " " << W.cols() << endl;
                    //cout << "fb dims: "<< fb.rows() << " " << fb.cols() << endl;

                    Eigen::VectorXd c = (W*B).householderQr().solve(W*fb);

                    //cout << "c:" << c.transpose() << endl;

                    grid_values[grid_index] = polynomialBasis(gridPointCoords).dot(c);
                }
            }
        }
    }
}

void alignPoints(Viewer &viewer){
    // Center
    P = P.rowwise() - P.colwise().mean();
    // Get centered covariance matrix
    Eigen::MatrixXd centered_cov = P.adjoint() * P;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(centered_cov);

    Eigen::RowVectorXd ev0 = eig.eigenvectors().col(0).transpose();
    Eigen::RowVectorXd ev1 = eig.eigenvectors().col(1).transpose();
    Eigen::RowVectorXd ev2 = eig.eigenvectors().col(2).transpose();

    // get BB for calcs
    Eigen::RowVector3d bb_min = P.colwise().minCoeff().array();
    double DBB = igl::bounding_box_diagonal(P);
    
    // Display ev in 3D
    // Eigen::MatrixXd ev_vertices(4,3);
    // ev_vertices << bb_min, bb_min+ev0*DBB, bb_min+ev1*DBB, bb_min+ev2*DBB;
    // Eigen::MatrixXd ev_line_start(3,3);
    // ev_line_start << bb_min, bb_min, bb_min; 
    // Eigen::MatrixXd ev_line_end(3,3);
    // ev_line_end << bb_min+ev0*DBB, bb_min+ev1*DBB, bb_min+ev2*DBB;           
    // viewer.data().point_size = 5;            
    // viewer.data().add_points(ev_vertices, Eigen::RowVector3d(0.5, 0.5, 0.5));
    // viewer.data().add_edges(ev_line_start, ev_line_end, Eigen::RowVector3d(0.5, 0.5, 0.5));

    // Align data with xyz axes
    P *= eig.eigenvectors().transpose().inverse();
    N *= eig.eigenvectors().transpose().inverse();
}

double S_k(Eigen::RowVectorXd x, int p_i){
    return (x - P.row(p_i)).dot(N.row(p_i));
}

// Called upon when "key == '5'"
void evaluateImplicitFunc_PolygonSoup()
{
    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate signed distance function at each gridpoint via MLS


    cout << "Evaluating with polyDegree=" << polyDegree << ", wendlandRadius=" << wendlandRadius << ", resolution=" << resolution << endl; 


    // Iterate over grid points computing weights for each point and performing one step of MLS
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                int grid_index = x + resolution * (y + resolution * z);
                Eigen::RowVector3d gridPointCoords = grid_points.row(grid_index);

                // Get neighbours in P
                std::vector<int> neighbours = withinDistanceSI(gridPointCoords, wendlandRadius);

                //cout << neighbours.size();

                // If within wendlandRadius empty, then set to large positive number
                if(neighbours.empty()){
                    grid_values[grid_index] = DBL_MAX;
                }
                // Otherwise find weights within neighbourhood and perform MLS step
                else{
                    // Create the polynomial basis functions in matrix form
                    Eigen::MatrixXd B(neighbours.size(), polynomialBasisLength());
                    for(int i = 0; i < B.rows(); i++){
                        B.row(i) = polynomialBasis(constrained_points.row(neighbours[i]));
                    }

                    // Set W
                    Eigen::VectorXd weights = Eigen::VectorXd::Zero(neighbours.size());
                    for(int i = 0; i < neighbours.size(); i++){
                        weights(i) = wendlandFunction(gridPointCoords, constrained_points.row(neighbours[i]));
                    }
                    Eigen::MatrixXd W = weights.asDiagonal();

                    // Set fb
                    Eigen::VectorXd fb(neighbours.size());
                    for(int i = 0; i < neighbours.size(); i++){
                        fb(i) = S_k(gridPointCoords, neighbours[i]);
                    }

                    //cout << "B dims: "<< B.rows() << " " << B.cols() << endl;
                    //cout << "W dims: "<< W.rows() << " " << W.cols() << endl;
                    //cout << "fb dims: "<< fb.rows() << " " << fb.cols() << endl;

                    Eigen::VectorXd c = (W*B).householderQr().solve(W*fb);

                    //cout << "c:" << c.transpose() << endl;

                    grid_values[grid_index] = polynomialBasis(gridPointCoords).dot(c);
                }
            }
        }
    }
}

// Estimation of the normals via PCA.
void pcaNormal()
{
    NP.resize(P.rows(), 3);

    // Build octree
    std::vector<std::vector<int > > O_PI;
    Eigen::MatrixXi O_CH;
    Eigen::MatrixXd O_CN;
    Eigen::VectorXd O_W;
    igl::octree(P,O_PI,O_CH,O_CN,O_W);

    Eigen::MatrixXi I;
    igl::knn(P,k_neighbourhood,O_PI,O_CH,O_CN,O_W,I);

    for(int i=0; i<P.rows(); i++){
        Eigen::MatrixXd C(k_neighbourhood, 3);
        for(int j = 0; j<k_neighbourhood; j++){
            C.row(j) = P.row(I(i,j));
        }
        // Get centered covariance matrix
        C = C.rowwise() - C.colwise().mean();
        Eigen::MatrixXd centered_cov_C = C.adjoint() * C;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_C(centered_cov_C);
        Eigen::VectorXd unsigned_normal = eig_C.eigenvectors().col(0);

        // Negate if pointing towards wrong direction
        if(N.row(i).dot(unsigned_normal) < 0){
            NP.row(i) = -unsigned_normal;
        }
        else{
            NP.row(i) = unsigned_normal;
        }
    }
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '0')
    {
        // Added key==0 for debugging / testing purposes
        auto start_time = std::chrono::high_resolution_clock::now();
        // Code to test time on
        auto end_time = std::chrono::high_resolution_clock::now();
        auto time = end_time - start_time;
        cout << "Time:" << time/std::chrono::nanoseconds(1) << endl;
    }
    if (key == '1')
    {
        viewer.data().clear();

        // Display axis in 3D
        double DBB = igl::bounding_box_diagonal(P);
        Eigen::MatrixXd o_vertices(4,3);
        o_vertices << 0,0,0,DBB,0,0,0,DBB,0,0,0,DBB;
        Eigen::MatrixXd o_line_start = Eigen::MatrixXd::Zero(3,3);
        Eigen::MatrixXd o_line_end(3,3);
        o_line_end << DBB,0,0,0,DBB,0,0,0,DBB;
        viewer.data().point_size = 5;            
        viewer.data().add_points(o_vertices, Eigen::RowVector3d(1, 1, 1));
        viewer.data().add_edges(o_line_start, o_line_end, Eigen::RowVector3d(1, 1, 1));

        // Align P with xyz axes using PCA
        if(alignP){
            alignPoints(viewer);
        }

        // Show imported points
        viewer.core().align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }

    if (key == '2')
    {

        // Preprocess the spatial index
        auto start_time = std::chrono::high_resolution_clock::now();
        preprocessSpatialIndex();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto time = end_time - start_time;
        cout << "Spatial Index Preprocessing took " << time/std::chrono::milliseconds(1) << "ms." << endl;

        // Start timer
        start_time = std::chrono::high_resolution_clock::now();
        
        // Show all constraints
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // get DBB
        double DBB = igl::bounding_box_diagonal(P);
        double e = 0.01;

        // Calculate constraint points & values for positive epsilon
        Eigen::MatrixXd constrained_points_e_pos(P.rows(), P.cols());
        Eigen::VectorXd constrained_values_e_pos(P.rows());

        calculateConstraints(constrained_points_e_pos, constrained_values_e_pos, e*DBB);

        // Calculate constraint points & values for positive epsilon
        Eigen::MatrixXd constrained_points_e_neg(P.rows(), P.cols());
        Eigen::VectorXd constrained_values_e_neg(P.rows());

        calculateConstraints(constrained_points_e_neg, constrained_values_e_neg, -e*DBB);

        constrained_points.resize(3*P.rows(), P.cols());
        constrained_values.resize(3*P.rows());

        constrained_points << P, constrained_points_e_pos, constrained_points_e_neg;
        constrained_values << Eigen::VectorXd::Zero(P.rows()), constrained_values_e_pos, constrained_values_e_neg;

        // Code for displaying all points, as above
        viewer.data().point_size = 11;
        viewer.data().add_points(constrained_points_e_pos, Eigen::RowVector3d(1, 0, 0));        
        viewer.data().add_points(constrained_points_e_neg, Eigen::RowVector3d(0, 1, 0));
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 1));

        end_time = std::chrono::high_resolution_clock::now();
        time = end_time - start_time;
        cout << "Constraints Calculation took " << time/std::chrono::milliseconds(1) << "ms." << endl;

        // cout << "=============================" << endl;
        // cout << "  constrained_points_e_pos   "   << endl;
        // cout << "=============================" << endl;
        // cout << constrained_points_e_pos << endl;

        // cout << "=============================" << endl;
        // cout << "  constrained_points_e_neg   "   << endl;
        // cout << "=============================" << endl;
        // cout << constrained_points_e_neg << endl;

        // cout << "=============================" << endl;
        // cout << "              P              " << endl;
        // cout << "=============================" << endl;
        // cout << P << endl;


    }
    if (key == '3')
    {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        // Do evaluateImplicitFunc:
        // Preprocess the spatial index
        auto start_time = std::chrono::high_resolution_clock::now();
        evaluateImplicitFunc();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto time = end_time - start_time;
        cout << "MLS took " << time/std::chrono::milliseconds(1) << "ms." << endl;

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
    }

    if (key == '4')
    {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
    }

    if (key == '5')
    {
        // Use the structure for key=='3' but replace the function evaluateImplicitFunc();
        // with a function performing the approximation of the implicit surface from polygon soup
        // Ref: Chen Shen, James F. O’Brien, and Jonathan Richard Shewchuk. Interpolating and approximating implicit surfaces from polygon soup.

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function --> Function to be modified here
        evaluateImplicitFunc_PolygonSoup();

        // get grid lines
        getLines();

        // Display the reconstruction
        // callback_key_down(viewer, '4', modifiers);

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));

    }

    if (key == '6' || key == '7' || key == '8')
    {
        // Implement PCA Normal Estimation --> Function to be modified here
        pcaNormal();

        // To use the normals estimated via PCA instead of the input normals and then restaurate the input normals
        Eigen::MatrixXd N_tmp = N;
        N = NP;

        switch (key)
        {
        case '6':
            callback_key_down(viewer, '2', modifiers);
            break;
        case '7':
            callback_key_down(viewer, '3', modifiers);
            break;
        case '8':
            callback_key_down(viewer, '4', modifiers);
            break;
        default:
            break;
        }

        // Restore input normals
        N = N_tmp;
    }
    // Save reconstructed mesh to /res
    if (key == '9'){
        string data_name = data_path.substr(8, data_path.substr(8).find("."));
        string write_path = "../res/" + data_name + "_reconstruction.off";
        igl::writeOFF(write_path, V,F);
    }

    return true;
}

bool callback_load_mesh(Viewer &viewer, string filename)
{
    igl::readOFF(filename, P, F, N);
    callback_key_down(viewer, '1', 0);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/sphere.off", P, F, N);
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1], P, F, N);
        data_path = argv[1];
    }

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputInt("Resolution", &resolution, 0, 0);
            ImGui::InputInt("Align P?", &alignP);
            ImGui::InputDouble("Padding Mult.", &padding_mult, 0.01);
            ImGui::InputDouble("Wendland Radius", &wendlandRadius, 0.01);
            ImGui::InputInt("Polynomial Degree", &polyDegree);
            ImGui::InputInt("k-Neighbourhood", &k_neighbourhood, 1);
            if (ImGui::Button("Reset Grid", ImVec2(-1, 0)))
            {
                cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer, '3', 0);
            }
        }
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
