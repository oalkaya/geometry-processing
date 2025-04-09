#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
/*** insert any libigl headers here ***/
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/edge_flaps.h>
#include <igl/is_border_vertex.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int>> VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

void subdivide_sqrt3(const Eigen::MatrixXd &V,
                     const Eigen::MatrixXi &F,
                     Eigen::MatrixXd &Vout,
                     Eigen::MatrixXi &Fout)
// TO DO:
// 1)
// Add vertices to midpoint of faces, M
// Append M to V
// Add three new faces using middle point and each edge, replaces F  with F"
// 2)
// Update vertex locations of original matrix resulting in P instead of V
// now V' = [P, M]
// 3)
// Replace F" with F' as shown in figure
{
    // Compute midpoints M
    Eigen::MatrixXd M(F.rows(), F.cols());
    igl::barycenter(V, F, M);

    // for(int i = 0;i<F.rows();i++){
    //    M.row(i) = (V.row(F(i,0)) + V.row(F(i,1)) + V.row(F(i,2)))/3;
    // }

    // Append midpoints M to V
    Eigen::MatrixXd V2(V.rows() + M.rows(), V.cols());
    V2 << V, M;

    // cout << "V rows:" << V.rows() << "cols:" <<V.cols() << endl;
    // cout << "M rows:" << M.rows() << "cols:" <<M.cols() << endl;
    // cout << "V2 rows:" << V2.rows() << "cols:" <<V2.cols() << endl;

    // Replace faces F with new faces F"
    Eigen::MatrixXi F2(F.rows() * 3, F.cols());
    for (int i = 0; i < F.rows(); i++)
    {
        // Use index F.rows() + i since V2 = [V, M]
        // Triangle i*3
        F2(i * 3, 0) = V.rows() + i;
        F2(i * 3, 1) = F(i, 1);
        F2(i * 3, 2) = F(i, 2);
        // Triangle i*3 + 1
        F2(i * 3 + 1, 0) = F(i, 0);
        F2(i * 3 + 1, 1) = V.rows() + i;
        F2(i * 3 + 1, 2) = F(i, 2);
        // Triangle i*3 + 2
        F2(i * 3 + 2, 0) = F(i, 0);
        F2(i * 3 + 2, 1) = F(i, 1);
        F2(i * 3 + 2, 2) = V.rows() + i;
    }

    // New vertex positions V1
    Eigen::MatrixXd P; //(V.rows(), V.cols());
    P.setZero(V.rows(), V.cols());

    std::vector<std::vector<int>> VV2;
    igl::adjacency_list(F, VV2);

    std::vector<bool> is_border = igl::is_border_vertex(F);

    // Iterate through each vertex
    for (int i = 0; i < V.rows(); i++)
    {
        // If boundary, no update
        if (is_border[i])
        {
            P.row(i) = V.row(i);
        }
        else
        {
            // Compute coefficients
            int n = VV2[i].size();
            float an = (4 - 2 * cosf(2 * M_PI / n)) / 9;
            // Iterate through incident vertices
            Eigen::MatrixXd NS;
            NS.setZero(1, 3);

            for (int j = 0; j < VV2[i].size(); j++)
            {
                NS += V.row(VV2[i][j]);
            }
            // Update rule
            P.row(i) = (1 - an) * V.row(i) + (an / n) * NS;
        }
    }
    Eigen::MatrixXd V1(V.rows() + M.rows(), V.cols());
    V1 << P, M;

    // Replace F2 with F1, flipping the original edges of mesh connecting points in P

    //// unique_edge_map
    Eigen::MatrixXi E, uE, EMAP;
    std::vector<std::vector<int>> uE2E;
    igl::unique_edge_map(F2, E, uE, EMAP, uE2E);

    // edge_flaps
    Eigen::MatrixXi EF, EI;
    igl::edge_flaps(F2, uE, EMAP, EF, EI);

    //// We use igl::edge_topology to construct the edge topology
    // Eigen::MatrixXi EV, FE, EF;
    // igl::edge_topology(V1, F2, EV, FE, EF);

    //// Check for boundary edges with igl::boundary_edge
    // Eigen::MatrixXi B;
    // igl::is_boundary_edge(EV, F2, B);

    Eigen::MatrixXi F1(F2.rows(), F2.cols());
    int nfi = 0; // new face (F1) index, incremented when we add new faces in iteration
    for (int i = 0; i < EF.rows(); i++)
    {
        // get edge flap faces
        // If boundary edge, one of the faces will be indexed -1
        int fI1 = EF(i, 0);
        int fI2 = EF(i, 1);

        // handle non-boundary edge first
        if ((fI1 != -1) && (fI2 != -1))
        {
            // edge flap corner indices (0-2)
            int cI1 = EI(i, 0);
            int cI2 = EI(i, 1);

            // get actual corner vertices
            int cV1 = F2(fI1, cI1);
            int cV2 = F2(fI2, cI2);

            // find the two vertices on the edge flap
            // EMAP maps E(i,0) and E(i,1) to the same unique edge
            // Indexing can be done over fI1 or fI2, same result

            int index = fI2 + (F2.rows() * cI2);
            int uEV1 = uE(EMAP(index), 0);
            int uEV2 = uE(EMAP(index), 1);

            // ONLY INCLUDE IF EDGE PRESENT IN ORIGINAL
            if (uEV1 < V.rows() && uEV2 < V.rows())
            {
                //// printout test
                // cout << "Face "<< nfi << ": "<< cV1 << ", "<< uEV1 << ", "<< cV2 << endl;
                // cout << "Face "<< nfi+1 << ": "<< cV2 << ", "<< uEV2 << ", "<< cV1 << endl;

                // construct the new triangles
                F1(nfi, 0) = cV1;
                F1(nfi, 1) = uEV1;
                F1(nfi, 2) = cV2;

                F1(nfi + 1, 0) = cV2;
                F1(nfi + 1, 1) = uEV2;
                F1(nfi + 1, 2) = cV1;
                nfi = nfi + 2;
            }
        }
        // fI1 is the only valid face
        else if (fI1 > fI2)
        {
            F1.row(nfi) = F2.row(fI1);
            nfi = nfi + 1;
        }
        // fI2 is the only valid face
        else if (fI2 > fI1)
        {
            F1.row(nfi) = F2.row(fI2);
            nfi = nfi + 1;
        }
        // handle boundary edge here
    }

    // Update outputs
    Vout = V1;
    Fout = F1;
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        igl::vertex_triangle_adjacency(V, F, VF, VFi);

        // Print VF
        cout << "Vertex: Face(s)" << endl;
        for (int i = 0; i < VF.size(); i++)
        {
            cout << i << ": ";
            for (int j = 0; j < VF[i].size(); j++)
            {
                cout << VF[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;

        // // Print VFi
        // cout << "VFi:" << endl;
        // for (int i = 0; i < VFi.size(); i++)
        // {
        //     for (int j = 0; j < VFi[i].size(); j++)
        //     {
        //     cout << VFi[i][j] << " ";
        //     }
        //     cout << endl;
        // }
    }

    if (key == '2')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
        igl::adjacency_list(F, VV);

        // Print VV
        cout << "Vertex: Adjacent Vertices" << endl;
        for (int i = 0; i < VV.size(); i++)
        {
            cout << i << ": ";
            for (int j = 0; j < VV[i].size(); j++)
            {
                cout << VV[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    if (key == '3')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(), 3);
        // Add your code for computing per-face normals here: store in FN.

        igl::per_face_normals(V, F, FN);

        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.

        igl::per_vertex_normals(V, F, VN);

        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
        igl::per_corner_normals(V, F, CN);
        // Set the viewer normals
        viewer.data().set_normals(CN);
    }

    if (key == '6')
    {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(), 3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        igl::facet_components(F, cid);

        std::map<int, int> cMap;
        for (int i = 0; i < cid.size(); i++)
        {
            if (cMap.find(cid(i)) != cMap.end())
            {
                cMap[cid(i)] = cMap[cid(i)] + 1;
            }
            else
            {
                cMap[cid(i)] = 1; // new cid element
            }
        }

        // Print statements needed for report
        cout << "Number of connected components: " << cMap.size() << endl;
        map<int, int>::iterator it;
        for (it = cMap.begin(); it != cMap.end(); it++)
        {
            it->first;
            cout << it->second << " "; // count per component
        }
        cout << endl << endl;

        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        igl::jet(cid, true, component_colors_per_face);

        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7')
    {
        Eigen::MatrixXd Vout;
        Eigen::MatrixXi Fout;
        // Fill the subdivide_sqrt3() function with your code for sqrt(3) subdivision.
        subdivide_sqrt3(V, F, Vout, Fout);
        // Set up the viewer to display the new mesh
        V = Vout;
        F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

bool load_mesh(Viewer &viewer, string filename, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    igl::readOFF(filename, V, F);
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(V, F);
    return true;
}

int main(int argc, char *argv[])
{
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;

    std::string filename;
    if (argc == 2)
    {
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else
    {
        filename = std::string("../data/bunny.off"); // Default mesh
    }

    load_mesh(viewer, filename, V, F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    viewer.launch();
}