#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>

/*** additional headers ***/
#include <igl/dijkstra.h>
#include <igl/colormap.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

bool showingUV = false;
bool freeBoundary = true;
bool showDistortion = false;
static int distortionType = 0;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

// Moved here from computeParametrisation
VectorXi fixed_UV_indices;
MatrixXd fixed_UV_positions;

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);
	// Show distortion
	if(showDistortion){
		// set distortion values
		// Get Dx and Dy
		SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);
		Eigen::VectorXd distortion_values(F.rows());
		//Conformal ie angle preserving
		if(distortionType == 0){
			for(int i = 0; i < F.rows(); i++){
				Eigen::Matrix2d J;
				J << Dx.row(i) * UV.col(0), Dy.row(i) * UV.col(0), Dx.row(i) * UV.col(1), Dy.row(i) * UV.col(1);
				distortion_values(i) = pow(((J + J.transpose() - (J.trace() * Eigen::Matrix2d::Identity())).norm()), 2);
			}
		}
		// Isometric ie length preserving
		else if(distortionType == 1){
			for(int i = 0; i < F.rows(); i++){
				Eigen::Matrix2d J, U, S, V, R;
				J << Dx.row(i) * UV.col(0), Dy.row(i) * UV.col(0), Dx.row(i) * UV.col(1), Dy.row(i) * UV.col(1);
				SSVD2x2(J,U,S,V);
				double det = (U*V.transpose()).determinant();
				Eigen::Matrix2d det_term = Eigen::Matrix2d::Identity();
				det_term(1,1) = det;
				R = U * det_term * V.transpose();
				distortion_values(i) = pow(((J-R).norm()), 2);
			}
		}
		// Authalic ie area preserving
		else if(distortionType == 2){
			for(int i = 0; i < F.rows(); i++){
				Eigen::Matrix2d J, U, S, V, R;
				J << Dx.row(i) * UV.col(0), Dy.row(i) * UV.col(0), Dx.row(i) * UV.col(1), Dy.row(i) * UV.col(1);
				distortion_values(i) = pow((J.determinant()-1), 2);
			}
		}

		// compute color values from distortion values
		distortion_values.normalize();
		Eigen::MatrixXd colors(F.rows(), 3);
		igl::colormap(igl::COLOR_MAP_TYPE_TURBO, distortion_values, false, colors); // true to use the whole range
		viewer.data().set_mesh(V, F);
		viewer.data().show_texture = false;
		viewer.data().set_colors(colors);
	}

    else if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}
}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	C.resize(2*indices.rows(), 2*V.rows());
	d.resize(2*indices.rows(), 1);

	// Watch out for uv Coordinates being stacked vertically
	for(int i = 0; i < indices.rows(); i++){
		// u coord
		C.insert(i, indices[i]) = 1;
		// v coord
		C.insert(indices.rows() + i, V.rows() + indices[i]) = 1;
	}

	for(int i = 0; i < positions.rows(); i++){
		// u coord
		d[i] = positions.row(i)[0];
		// v coord
		d[positions.rows() + i] = positions.row(i)[1];
	}
}

void computeParameterization(int type)
{
	SparseMatrix<double> A;
	VectorXd b;

	Eigen::SparseMatrix<double> C;
	VectorXd d;

	// Condition to skip new boundary computation
	if(!(type == '4')){
		// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
		if (!freeBoundary)
		{
			// The boundary vertices should be fixed to positions on the unit disc. Find these position and
			// save them in the #V x 2 matrix fixed_UV_position.

			// Get boundary loop, save it in fixed_UV_indices
			igl::boundary_loop(F, fixed_UV_indices);
			igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
		}
		else
		{
			// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.

			// set size
			fixed_UV_indices.resize(2);

			// Dijkstra setup:

			// declare variables and init adjacency list
			Eigen::VectorXd min_distance(V.rows());
			Eigen::VectorXi previous(V.rows());
			std::vector<std::vector<int>> adj_list;
			igl::adjacency_list(F, adj_list);

			// Perform Dijkstra n_samples times:

			int source;
			int target;
			double maxDist = 0;

			// Heuristically set params
			int n_runs = 10;
			int d_per_run = 3;

			// Run algo for n_runs
			for(int i = 0; i < n_runs; i++){
				// Variable setup for dijkstra
				double maxDist_i;
				int source_i = std::rand() % (V.rows()); // randomly sample first source between all vertices
				int prev_source_i;
				int target_i;

				// Within each run, perform dijkstra d_per_run times.
				for(int j = 0; j < d_per_run; j++){
					std::set<int> targets = {};
					igl::dijkstra(source_i, targets, adj_list, min_distance, previous);
					maxDist_i = min_distance.maxCoeff(&target_i);

					prev_source_i = source_i;
					source_i = target_i;
				}
				if(maxDist < maxDist_i){
						//cout << "UPDATE" << endl;
						maxDist = maxDist_i;
						source = prev_source_i;
						target = source_i;
				}
			}

			cout << "------------------------" << endl;
			cout << "First index chosen:  " 	<< source << endl;
			cout << "Second index chosen: " 	<< target << endl;
			cout << "------------------------" << endl;

			// fill indices and positions
			fixed_UV_indices << source, target;
			igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
		}
	}

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh

		// init b
		b = VectorXd::Zero(2*V.rows(), 1);

		// init uniform laplacian
		Eigen::SparseMatrix<double> adjacency;
		igl::adjacency_matrix(F, adjacency);

		Eigen::VectorXd cumsum;
		igl::sum(adjacency, 1, cumsum);

		for(int i = 0; i < adjacency.rows(); i++){
			adjacency.coeffRef(i,i) = -cumsum[i];
		}


		// init A by stacking (dim = 1 vertical, dim = 2 horizontal)
		// C 0
		// 0 C

		Eigen::SparseMatrix<double> A_col1, A_col2;
		igl::cat(1, adjacency, Eigen::SparseMatrix<double>(adjacency.rows(), adjacency.cols()), A_col1);
		igl::cat(1, Eigen::SparseMatrix<double>(adjacency.rows(), adjacency.cols()), adjacency, A_col2);
		igl::cat(2, A_col1, A_col2, A);

	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~

		// init b
		b = VectorXd::Zero(2*V.rows(), 1);

		// get cotangent matrix
		Eigen::SparseMatrix<double> cotangent;
		igl::cotmatrix(V, F, cotangent);

		// init A by stacking
		// Cot 0
		// 0 Cot
		Eigen::SparseMatrix<double> A_col1, A_col2;
		igl::cat(1, cotangent, Eigen::SparseMatrix<double>(cotangent.rows(), cotangent.cols()), A_col1);
		igl::cat(1, Eigen::SparseMatrix<double>(cotangent.rows(), cotangent.cols()), cotangent, A_col2);
		igl::cat(2, A_col1, A_col2, A);

	}

	if (type == '3') {

		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!

		// init b
		b = VectorXd::Zero(2*V.rows(), 1);

		// Get Dx and Dy
		SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);

		// Get doubled area
		VectorXd dblArea;
		igl::doublearea(V, F, dblArea);

		SparseMatrix<double> Area;
		Area.resize(dblArea.rows(), dblArea.rows());
		for(int i = 0; i < dblArea.rows(); i++){
			Area.coeffRef(i,i) = sqrt(dblArea[i]/2);
		}

		// Mult A^(1/2) with Dx or Dy
		SparseMatrix<double> ADx(Area * Dx);
		SparseMatrix<double> ADy(Area * Dy);

		// very pretty looking concat
		SparseMatrix<double> A_row1, A_row2, A_row3, A_row4, A_rows12, A_rows34;
		igl::cat(2, ADx, SparseMatrix<double>(-ADy), A_row1);
		igl::cat(2, ADy, ADx, A_row2);
		igl::cat(2, ADy, ADx, A_row3);
		igl::cat(2, SparseMatrix<double>(-ADx), ADy, A_row4);
		igl::cat(1, A_row1, A_row2, A_rows12); // rows 1 and 2 into 1
		igl::cat(1, A_row3, A_row4, A_rows34); // rows 3 and 4 into 3
		igl::cat(1, A_rows12, A_rows34, A); // rows 1 and 3 into A

		// Calc A.T * A
		A = SparseMatrix<double>(A.transpose() * A);

	}

	if (type == '4') {

		// // Add your code for computing ARAP system and right-hand side
		// // Implement a function that computes the local step first
		// // Then construct the matrix with the given rotation matrices

		// -------EXPLICIT MATRIX CALC-------

		// // Get Dx and Dy
		// SparseMatrix<double> Dx, Dy;
		// computeSurfaceGradientMatrix(Dx, Dy);

		// // Get doubled area
		// VectorXd dblArea;
		// igl::doublearea(V, F, dblArea);

		// SparseMatrix<double> Area;
		// Area.resize(dblArea.rows(), dblArea.rows());
		// for(int i = 0; i < dblArea.rows(); i++){
		// 	Area.coeffRef(i,i) = sqrt(dblArea[i]/2);
		// }

		// // Mult A^(1/2) with Dx or Dy
		// SparseMatrix<double> ADx(Area * Dx);
		// SparseMatrix<double> ADy(Area * Dy);

		// // very pretty looking concat
		// SparseMatrix<double> zeros(Dx.rows(), Dx.cols());

		// SparseMatrix<double> A_row1, A_row2, A_row3, A_row4, A_rows12, A_rows34;

		// igl::cat(2, ADx, zeros, A_row1);
		// igl::cat(2, ADy, zeros, A_row2);
		// igl::cat(2, zeros, ADx, A_row3);
		// igl::cat(2, zeros, ADy, A_row4);
		// igl::cat(1, A_row1, A_row2, A_rows12); // rows 1 and 2 into 1
		// igl::cat(1, A_row3, A_row4, A_rows34); // rows 3 and 4 into 3
		// igl::cat(1, A_rows12, A_rows34, A); // rows 1 and 3 into A

		// Eigen::VectorXd Rj_vectorized(4*Dx.rows(), 1);
		// for(int i = 0; i < Dx.rows(); i++){

		// 	// J calculation
		// 	Eigen::Matrix2d J, U, S, V;
		// 	J << Dx.row(i).dot(UV.col(0)), Dy.row(i).dot(UV.col(0)), Dx.row(i).dot(UV.col(1)), Dy.row(i).dot(UV.col(1));

		// 	// Do SVD
		// 	SSVD2x2(J,U,S,V);

		// 	// Get det_term
		// 	Eigen::Matrix2d det_term;
		// 	det_term << 1, 0, 0, (U * V.transpose()).determinant();

		// 	// Get Rj
		// 	Matrix2d Rj = U * det_term * V.transpose();

		// 	// Fill Rj
		// 	double tra = sqrt(dblArea[i]/2);
		// 	Rj_vectorized[i] 			   = Rj(0,0) * tra;
		// 	Rj_vectorized[i +   Dx.rows()] = Rj(0,1) * tra;
		// 	Rj_vectorized[i + 2*Dx.rows()] = Rj(1,0) * tra;
		// 	Rj_vectorized[i + 3*Dx.rows()] = Rj(1,1) * tra;
		// }

		// // Set b for solver
		// b = A.transpose() * Rj_vectorized;

		// // Set A for solver
		// A = SparseMatrix<double>(A.transpose() * A);




		// -------IMPLICIT MATRIX CALC-------

		// Get Dx and Dy
		SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);

		// Calc Rj
		Eigen::VectorXd Rj_00 (Dx.rows());
		Eigen::VectorXd Rj_01 (Dx.rows());
		Eigen::VectorXd Rj_10 (Dx.rows());
		Eigen::VectorXd Rj_11 (Dx.rows());
		for(int i = 0; i < Dx.rows(); i++){
			// J calculation
			Eigen::Matrix2d J, U, S, V, Rj;
			J << Dx.row(i) * UV.col(0), Dy.row(i) * UV.col(0), Dx.row(i) * UV.col(1), Dy.row(i) * UV.col(1);
			// Do SVD
			SSVD2x2(J,U,S,V);
			// Get det_term
			double det = (U*V.transpose()).determinant();
			Eigen::Matrix2d det_term = Eigen::Matrix2d::Identity();
			det_term(1,1) = det;
			Rj = U * det_term * V.transpose();

			// Fill Rj
			Rj_00[i] = Rj(0,0);
			Rj_01[i] = Rj(0,1);
			Rj_10[i] = Rj(1,0);
			Rj_11[i] = Rj(1,1);
		}

		// Calc Area
		VectorXd dblArea;
		igl::doublearea(V, F, dblArea);
		SparseMatrix<double> Area;
		igl::diag(dblArea, Area);
		Area = Area/2;

		// Calc A
		SparseMatrix<double> zeros(V.rows(), V.rows());
		SparseMatrix<double> A_row1, A_row2;
		SparseMatrix<double> L(V.rows(), V.rows());
		L = Dx.transpose() * Area * Dx + Dy.transpose() * Area * Dy;
		igl::cat(2, L, zeros, A_row1);
		igl::cat(2, zeros, L, A_row2);
		igl::cat(1, A_row1, A_row2, A);


		// Calc b
		VectorXd b_row1 = Dx.transpose() * Area * Rj_00 + Dy.transpose() * Area * Rj_01;
		VectorXd b_row2 = Dx.transpose() * Area * Rj_10 + Dy.transpose() * Area * Rj_11;
		b.resize(b_row1.size() + b_row2.size());
		b << b_row1, b_row2;


	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail

	// The solver will output a vector

	// Initialize LHS SparseMatrix
	// A C.T
	// C 0
	
	Eigen::SparseMatrix<double> LHS;
	Eigen::SparseMatrix<double> LHS_col1, LHS_col2;
	igl::cat(1, A, C, LHS_col1);
    igl::cat(1, Eigen::SparseMatrix<double>(C.transpose()), Eigen::SparseMatrix<double>(C.rows(), C.rows()), LHS_col2);
    igl::cat(2, LHS_col1, LHS_col2, LHS);

	// Initialize RHS SparseMatrix
	Eigen::VectorXd RHS;
    igl::cat(1, b, d, RHS);

	Eigen::SparseLU<Eigen::SparseMatrix<double>> SparseLUSolver;
    SparseLUSolver.compute(LHS);

    Eigen::VectorXd stacked_result = SparseLUSolver.solve(RHS);

	UV.resize(V.rows(), 2);
	UV.col(0) = stacked_result.block(0, 0, V.rows(), 1);
	UV.col(1) = stacked_result.block(V.rows(), 0, V.rows(), 1);
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		break;
	case '5':
		// Update drawing:
		Redraw();
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core();
      viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core();
        viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core().align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/cathead.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

  igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parametrization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);

			// TODO: Add more parameters to tweak here...
			ImGui::Checkbox("Show Distortion", &showDistortion);
			static const char* distortions[] = {"Angle", "Area", "Length"};
			ImGui::Combo("Distortion", &distortionType, distortions, 3);
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}
