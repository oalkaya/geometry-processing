#include "Lasso.h"

#include <iostream>
#include <fstream>

#include <igl/project.h>
#include <igl/unproject.h>
#include <igl/facet_components.h>
#include <igl/barycenter.h>

#include <igl/unproject_onto_mesh.h>

using namespace igl;
using namespace std;

Lasso::Lasso(const Eigen::MatrixXd &V_,
             const Eigen::MatrixXi &F_,
             const Viewer &v):
V(V_),
F(F_),
viewer(v)
{
}

Lasso::~Lasso()
{
}

int Lasso::pickVertex(int mouse_x, int mouse_y)
{
  int vi = -1;

  int fid;
  Eigen::Vector3f bc;
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core().viewport(3) - viewer.current_mouse_y;
  if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core().view /* viewer.data().model*/,
                              viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
  {
    // paint hit red
    bc.maxCoeff(&vi);
    vi = F(fid,vi);
  }
  return vi;

}
int Lasso::strokeAdd(int mouse_x,
                    int mouse_y)
{
  // Cast a ray in the view direction starting from the mouse position
  double x = mouse_x;
  double y = viewer.core().viewport(3) - mouse_y;

  std::vector<unsigned> pt2D; pt2D.push_back(x); pt2D.push_back(y);
  stroke2DPoints.push_back(pt2D);

  Eigen::RowVector3d pt(0,0,0);
  int fi = -1;

  Eigen::Matrix4f modelview = viewer.core().view; // * viewer.data().model;

  if (d<0)//first time
  {
    Eigen::Vector3f bc;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x,y),
                           modelview,
                           viewer.core().proj,
                           viewer.core().viewport,
                           V,
                           F,
                           fi,
                           bc))
    {
      pt = V.row(F(fi,0))*bc(0) + V.row(F(fi,1))*bc(1) + V.row(F(fi,2))*bc(2);
      Eigen::Vector3f proj = igl::project(pt.transpose().cast<float>().eval(), modelview, viewer.core().proj,viewer.core().viewport);
      d = proj[2];
    }

  }

  // This is lazy, it will find more than just the first hit
  pt = igl::unproject(Eigen::Vector3f(x,y,0.95*d), modelview, viewer.core().proj, viewer.core().viewport).transpose().cast<double>();

  strokePoints.push_back(pt);

  return fi;

}

namespace {
  bool IGL_INLINE point_in_poly( const std::vector<std::vector<unsigned int > >&poly, 
            const unsigned int xt, 
            const unsigned int yt)
{
  int npoints= poly.size();
  unsigned int xnew,ynew;
  unsigned int xold,yold;
  unsigned int x1,y1;
  unsigned int x2,y2;
  int i;
  int inside=0;
  
  if (npoints < 3) {
    return(0);
  }
  xold=poly[npoints-1][0];
  yold=poly[npoints-1][1];
  for (i=0 ; i < npoints ; i++) {
    xnew=poly[i][0];
    ynew=poly[i][1];
    if (xnew > xold) {
      x1=xold;
      x2=xnew;
      y1=yold;
      y2=ynew;
    }
    else {
      x1=xnew;
      x2=xold;
      y1=ynew;
      y2=yold;
    }
    if ((xnew < xt) == (xt <= xold)          /* edge "open" at one end */
        && ((long)yt-(long)y1)*(long)(x2-x1)
        < ((long)y2-(long)y1)*(long)(xt-x1)) {
      inside=!inside;
    }
    xold=xnew;
    yold=ynew;
  }
  return (inside != 0);
}
}

void Lasso::strokeFinish(Eigen::VectorXi &selected_vertices)
{

  Eigen::Matrix4f modelview = viewer.core().view; // * viewer.data().model;

  //marker for selected vertices
  Eigen::VectorXi is_selected; is_selected.setZero(V.rows(),1);

  //project all vertices, check which ones land inside the polyline
  for (int vi =0; vi<V.rows(); ++vi)
  {
    Eigen::Vector3f vertex = V.row(vi).transpose().cast<float>();
    Eigen::Vector3f proj = igl::project(vertex, modelview, viewer.core().proj,viewer.core().viewport);
    if (point_in_poly(stroke2DPoints, proj[0], proj[1]))
      is_selected[vi] = 1;

  }

  //the selection might consist of front facing and back facing facets.
  //we will only select the connected component that is frontmost

  //first, determine faces that have at least one selected vertex
  int nf = 0;
  Eigen::MatrixXi Fsel(F.rows(), 3);
  for (int fi = 0; fi<F.rows(); ++fi) {
    for (int i = 0; i < 3; ++i) {
      if (is_selected[F(fi, i)]) {
        Fsel.row(nf++) = F.row(fi);
        break;
      }
    }
  }
  Fsel.conservativeResize(nf, Eigen::NoChange);

  // Determine the frontmost component, if there is any selection
  if (Fsel.rows() > 0) {
    //compute their barycenters
    Eigen::MatrixXd MFsel;
    igl::barycenter(V, Fsel, MFsel);

    //now, find all connected components of selected faces
    Eigen::VectorXi cid;
    igl::facet_components(Fsel, cid);

    //compute centroids of connected components
    int ncomp = cid.maxCoeff()+1;
    Eigen::MatrixXd region_centroids;
    region_centroids.setZero(ncomp,3);
    Eigen::VectorXi total; total.setZero(ncomp,1);
    for (long fi = 0; fi<Fsel.rows(); ++fi)
    {
      int r = cid[fi];
      region_centroids.row(r) += MFsel.row(fi);
      total[r]++;
    }
    for (long i = 0; i<ncomp; ++i)
      region_centroids.row(i) = region_centroids.row(i).array()/total[i];

    //project all centroids and isolate only the most frontal one
    float mind = 1e10;
    int r = -1;
    for (long i = 0; i<ncomp; ++i)
    {
      Eigen::Vector3f t = region_centroids.row(i).transpose().cast<float>();
      Eigen::Vector3f proj = igl::project(t,
                                          modelview,
                                          viewer.core().proj,
                                          viewer.core().viewport);
      float depth = proj[2];
      if (mind > depth)
      {
        r = i;
        mind = depth;
      }
    }

    //all vertices belonging to other components are unmarked
    for (long fi = 0; fi<Fsel.rows(); ++fi)
    {
      if (cid[fi] != r)
      for (int i = 0; i<3; ++i)
        is_selected[Fsel(fi,i)] = 0;
    }
  }

  //return the selected vertices
  selected_vertices.resize(is_selected.sum(), 1);
  int num = 0;
  for (int vi =0; vi < V.rows(); ++vi)
    if (is_selected[vi])
      selected_vertices[num++] = vi;

  strokePoints.clear();
  stroke2DPoints.clear();
  d = -1;
}

void Lasso::strokeReset()
{
  strokePoints.clear();
  stroke2DPoints.clear();
  d = -1;
}
