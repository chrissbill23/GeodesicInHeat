
#include "TriangleMeshHeat.h"
#include "Utils.h"
#include <igl/per_face_normals.h>
#include <igl/doublearea.h>
#include <igl/cotmatrix.h>
#include <igl/grad.h>
#include <igl/avg_edge_length.h>

void TriangleMeshHeat::init_attr(){
  laplacianTriMesh(*Vertices, *Faces, L);
  //L = -0.5*L;
  cout << " FACTORIZATION OF LAPLACIAN OPERATOR"<<endl;
  laplacianFactor = LLT<MatrixXd> (-0.5*L);
   massMatrix(*Vertices, *Faces, M, M_inv);
    findBoundaries();
}
void TriangleMeshHeat::computeTangentNormals(){
  // compute the gradient operators and normals, using igl::grad and igl::per_face_normals functions
  igl::per_face_normals(*Vertices, *Faces,normals);
  SparseMatrix<double> G;
  igl::grad(*Vertices, *Faces,G);
  MatrixXd& V = *Vertices; MatrixXi& F = *Faces;
  tangents = MatrixXd::Zero(F.rows()*3, 3);//igl::grad return 3 gradients (one for each edge ) per triangle face
  for(size_t i = 0; i < F.rows(); ++i){
     tangents(i,0) = G.coeffRef(i, F(i,0));
     tangents(i,1) = G.coeffRef(i, F(i,1));
     tangents(i,2) = G.coeffRef(i, F(i,2));
     
     tangents(i+F.rows(),0) = G.coeffRef(i+F.rows(), F(i,0));
     tangents(i+F.rows(),1) = G.coeffRef(i+F.rows(), F(i,1));
     tangents(i+F.rows(),2) = G.coeffRef(i+F.rows(), F(i,2));

     tangents(i+F.rows()*2,0) = G.coeffRef(i+F.rows()*2, F(i,0));
     tangents(i+F.rows()*2,1) = G.coeffRef(i+F.rows()*2, F(i,1));
     tangents(i+F.rows()*2,2) = G.coeffRef(i+F.rows()*2, F(i,2));
  }
}
void TriangleMeshHeat::computeTimeStep(){
     double l = igl::avg_edge_length(*Vertices,*Faces);
     time = l * l * smooth;
}
void TriangleMeshHeat::findBoundaries(){
  MatrixXd& V = *Vertices;
    for (int i = 0; i < V.rows(); i++) {
        int edge1 = he.getEdge(i);
        int face1 = he.getFace(edge1);
        int edge2 = he.getOpposite(edge1);
        int face2 = he.getFace(edge2);
        if (face1 == -1 || face2==-1) {
            boundaryIndexes.push_back(he.getTarget(edge1));
            boundaryIndexes.push_back(he.getTarget(edge2));
        }
    }
}
void TriangleMeshHeat::solveHeatNeuman(bool boundary){
     heatNeuman = heatFlowOperatorFactor.solve(sources);
}
void TriangleMeshHeat::solveHeatDirichlet(bool boundary){
     MatrixXd& V = *Vertices; MatrixXi& F = *Faces;
     if(boundary && boundaryIndexes.size() > 0){
        heatDirichlet = heatFlowOperatorFactor.solve(sources);
        std::list<int>::iterator it;
        int count = boundaryIndexes.size();
        while (count > 0) {
            int index = boundaryIndexes.front();
            boundaryIndexes.pop_front();
            //V.row(index)[0] = 0; V.row(index)[1] = 0; V.row(index)[2] = 0;
            heatDirichlet[index] = 0;
            --count;
        }
     } else{
         heatDirichlet = heatFlowOperatorFactor.solve(sources);
     }
}
void TriangleMeshHeat::solveHeat(){
     switch (boundarycondition) {
          case PURE_NEUMAN : solveHeatNeuman(); heat = heatNeuman; break;
          case PURE_DIRICHLET : solveHeatDirichlet(); heat = heatDirichlet; break;
          case BOUNDARY_COND_NEUMAN: solveHeatNeuman(true); heat = heatNeuman; break;
          case BOUNDARY_COND_DIRICHLET: solveHeatDirichlet(true); heat = heatDirichlet; break;
          default : solveHeatNeuman(true); solveHeatDirichlet(true); heat = 0.5* heatNeuman + 0.5*heatDirichlet;
     }
}
void TriangleMeshHeat::solveVectorField(){
     
     MatrixXd& V = *Vertices; MatrixXi& F = *Faces;
     poisson.setZero(V.rows());
     //X.setZero(V.rows(), 3);
     VectorXd dblA;
     igl::doublearea(V,F,dblA);
     for(int f = 0; f < F.rows(); ++f){

         int v0 = F(f,0);
         int v1 = F(f,1);
         int v2 = F(f,2);
         
         double u0 = heat(v0), u1 = heat(v1), u2 = heat(v2);
         VectorXd gradU(3);
         gradU(0) = u0 * tangents(f,0) + u1 * tangents(f,1) + u2 * tangents(f,2);
         gradU(1) = u0 * tangents(f+F.rows(),0) + u1 * tangents(f+F.rows(),1) + u2 * tangents(f+F.rows(),2);
         gradU(2) = u0 * tangents(f+F.rows()*2,0) + u1 * tangents(f+F.rows()*2,1) + u2 * tangents(f+F.rows()*2,2);
         gradU.normalize();

         poisson(v0) += 0.5* dblA(f) * ( gradU(0)*tangents(f,0) + gradU(1)*tangents(f+F.rows(),0) + gradU(2)*tangents(f+F.rows()*2,0) );
         poisson(v1) += 0.5* dblA(f) * ( gradU(0)*tangents(f,1) + gradU(1)*tangents(f+F.rows(),1) + gradU(2)*tangents(f+F.rows()*2,1) );
         poisson(v2) += 0.5* dblA(f) * ( gradU(0)*tangents(f,2) + gradU(1)*tangents(f+F.rows(),2) + gradU(2)*tangents(f+F.rows()*2,2) );

     }
     /*for(int f = 0; f < F.rows(); ++f){

         int v0 = F(f,0);
         int v1 = F(f,1);
         int v2 = F(f,2);
         
         double u0 = heat(v0), u1 = heat(v1), u2 = heat(v2);
         double molt = max( max( u0, u1 ), u2 );
         if(molt > 0){
            molt = 1./molt;
            u0*= molt; u1*= molt; u2*= molt;

            RowVector3d e0 = V.row(v2) - V.row(v1);
            RowVector3d e1 = V.row(v2) - V.row(v0);
            RowVector3d e2 = V.row(v1) - V.row(v0);

            RowVector3d N = normals.row(f).normalized();
            
            e0 *= (cotangent(e1,e2)*0.5);//cot.coeffRef(0,1);
            e1 *= (cotangent(e0,e2)*0.5);//cot.coeffRef(2,1);
            e2 *= (cotangent(e1,e0)*0.5);//cot.coeffRef(0,2);

         
            RowVector3d ne0 = N.cross(e0), ne1 = N.cross(e1), ne2 = N.cross(e2);
         
            X.row(f) += u0 * ne0;
            X.row(f) += u1 * ne1;
            X.row(f) += u2 * ne2;

            double a = dblA(f);
            a = pow(a, -1.);
            X.row(f) *= a;
            X.row(f).normalize();
            //X.row(f) *= - 1.;

            poisson(v0) -= e1.dot(X.row(f)) - e2.dot(X.row(f));
            poisson(v1) -= e2.dot(X.row(f)) - e0.dot(X.row(f));
            poisson(v2) -= e0.dot(X.row(f)) - e1.dot(X.row(f));
          }  
     }*/
}
void TriangleMeshHeat::solvePoisson(){
     poisson = laplacianFactor.solve(poisson);
     double mean = 0.;
     for(size_t i = 0; i < sources.rows(); ++i){
        if(sources(i) == 1)
           mean += poisson(i);
     }
     mean /= sources.sum();
     poisson = poisson - VectorXd::Constant(sources.rows(),mean);
     if(poisson.mean() < 0)
        poisson *= -1;
     /*poisson.setZero(V.rows(),1);
     for(int f = 0; f < F.rows(); ++f){
         int edge = he.getEdgeInFace(f);
         int opp = he.getOpposite(edge);
         int count = 3;
         while(count > 0){
            int v1 = he.getTarget(edge);
            //if(sources(v1) == 0){
            int v2 = he.getTarget(opp);
            int v3 = he.getTarget(he.getNext(edge));
            MatrixXd vtmp(3,3); SparseMatrix<double> cot;
            RowVector3d e1 = (V.row(v2) - V.row(v1)).normalized();
            RowVector3d e2 = (V.row(v3) - V.row(v1)).normalized();
            RowVector3d e3 = (V.row(v2) - V.row(v3)).normalized();
            vtmp.row(0) = e1; vtmp.row(1) = e2; vtmp.row(2) = e3;
            igl::cotmatrix(vtmp,RowVector3d(0,1,2),cot);
            double cot2 = cot.coeffRef(0,2);
            double cot1 = cot.coeffRef(1,2);

            
            int f1 = he.getFace(opp);
            int f2 = he.getFace(he.getOpposite(he.getNext(edge)));
            
            poisson(v1) += cot1*(e1.dot(X.row(f))) + cot2*(e2.dot(X.row(f)));
            if(f1 != -1)
            poisson(v1) += (cot1*(e1.dot(X.row(f1))) + cot2*(e2.dot(X.row(f1))));
            if(f2 != -1)
            poisson(v1) += (cot1*(e1.dot(X.row(f2))) + cot2*(e2.dot(X.row(f2))));

            //cout << f<<"("<<v1<<","<<cot1*e1.dot(X.row(f))<<","<<cot1*e1.dot(X.row(f1))<<","<<cot1*e1.dot(X.row(f2))<<")" << "     "
                 //<<"("<<cot2*e2.dot(X.row(f))<<","<<cot2*e2.dot(X.row(f1))<<","<<cot2*e2.dot(X.row(f2))<<")"<<endl;
            //}
            //else{ poisson(v1) = -1;}
            edge = he.getNext(edge);
            opp = he.getOpposite(edge);
            --count;
         }
         cout<<endl<<endl<<endl;
     }
     poisson*= 0.5 ;*/
}



