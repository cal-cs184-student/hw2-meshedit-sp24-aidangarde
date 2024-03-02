#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    std::vector<Vector2D> newPoints;

    for(int i = 0; i < points.size() - 1; i++) {
      Vector2D newPoint = (1-t)*points[i] + t*points[i+1];
      newPoints.push_back(newPoint);
    }

    return newPoints;
    // return std::vector<Vector2D>();
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.

    std::vector<Vector3D> newPoints;


    for (int i = 0; i < points.size() - 1; i++){
      Vector3D newPoint = (1-t) * points[i] + t*points[i+1];
      newPoints.push_back(newPoint);
    }

    return newPoints;

    // return std::vector<Vector3D>();
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    if(points.size() == 1){
      return points[0];
    }

    std::vector<Vector3D> newPoints = evaluateStep(points, t);
    return evaluate1D(newPoints, t);
    // return Vector3D();
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
    std::vector<Vector3D> newPoints;

    for (int i = 0; i < controlPoints.size(); i++) {
      Vector3D pU = evaluate1D(controlPoints[i], u);
      newPoints.push_back(pU);
    }

    return evaluate1D(newPoints, v);
    // return Vector3D();
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.

    struct area_normal_pair {
      Vector3D normal;
      double area;
    };

    HalfedgeCIter h = this->halfedge();

    std::list<area_normal_pair> pairs;
    do {
      FaceCIter f = h->face();
      HalfedgeCIter h_iter = f->halfedge();
      VertexCIter vertices[3];
      int index = 0;
      do {
        VertexCIter v_iter = h_iter->vertex();
        vertices[index] = v_iter;
        h_iter = h_iter->next();
        index++;
      } while (h_iter != f->halfedge());
    
      Vector3D edge1;
      edge1.x = vertices[1]->position.x - vertices[0]->position.x;
      edge1.y = vertices[1]->position.y - vertices[0]->position.y;
      edge1.z = vertices[1]->position.z - vertices[0]->position.z;

      // Initialize edge2
      Vector3D edge2;
      edge2.x = vertices[2]->position.x - vertices[0]->position.x;
      edge2.y = vertices[2]->position.y - vertices[0]->position.y;
      edge2.z = vertices[2]->position.z - vertices[0]->position.z;

      // Initialize normal
      Vector3D normal;
      normal.x = edge1.y * edge2.z - edge2.y * edge1.z;
      normal.y = edge2.x * edge1.x - edge1.x * edge2.z;
      normal.z = edge1.x * edge2.y - edge2.x * edge1.y;

      double area = 0.5 * std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

      area_normal_pair pair = {normal, area};
      
      pairs.push_back(pair);

      h = h->twin()->next();

    } while (h != this->halfedge());

    Vector3D weightedSum = Vector3D();
    double total = 0.0;

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      weightedSum.x += it->normal.x * it->area;
      weightedSum.y += it->normal.y * it->area;
      weightedSum.z += it->normal.z * it->area;

      total += it->area;
    }

    weightedSum.x /= total;
    weightedSum.y /= total;
    weightedSum.z /= total;

    double length = std::sqrt(weightedSum.x * weightedSum.x + weightedSum.y * weightedSum.y + weightedSum.z * weightedSum.z);
    if (length != 0.0) {
        return Vector3D(weightedSum.x / length, weightedSum.y / length, weightedSum.z / length);
    } else {
        return Vector3D();
    }
    // return weightedSum;
}






  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    if (e0 -> isBoundary()) {
      return e0;
    }

    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h5->twin();
    HalfedgeIter h7 = h4->twin();
    HalfedgeIter h8 = h2->twin();
    HalfedgeIter h9 = h1->twin();

    VertexIter v0 = h0->vertex();
    VertexIter v1 = h1->vertex();
    VertexIter v2 = h2->vertex();
    VertexIter v3 = h5->vertex();

    // EdgeIter e0 = e0;
    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h4->edge();
    EdgeIter e4 = h5->edge();

    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    

    v0->halfedge() = h4;
    v1->halfedge() = h1;
    v2->halfedge() = h2;
    v3->halfedge() = h5;

    f0->halfedge() = h0;
    f1->halfedge() = h3;

    e0->halfedge() = h0;
    e1->halfedge() = h1;
    e2->halfedge() = h2;
    e3->halfedge() = h4;
    e4->halfedge() = h5;

    h0->setNeighbors(h2, h3, v3, e0, f0);
    h1->setNeighbors(h3, h9, v1, e1, f1);
    h2->setNeighbors(h4, h8, v2, e2, f0);
    h3->setNeighbors(h5, h0, v2, e0, f1);
    h4->setNeighbors(h0, h7, v0, e3 ,f0);
    h5->setNeighbors(h1, h6, v3, e4, f1);
    // h6->setNeighbors(h6->next(), h5, v1, e4, h6->face());
    // h7->setNeighbors(h7->next(), h4, v3, e3, h8->face());
    // h8->setNeighbors(h8->next(), h2, v0, e2, h8->face());
    // h9->setNeighbors(h9->next(), h1, v2, e1, h9->face());
    return e0;

    








    // HalfedgeIter h0 = e0->halfedge();
    // HalfedgeIter h1 = h0->twin();

    // VertexIter vb = h0->vertex();
    // VertexIter vc = h1->vertex();
    // VertexIter va = h0->next()->next()->vertex();
    // VertexIter vd = h1->next()->next()->vertex();

   

    // FaceIter f0 = h0->face();
    // FaceIter f1 = h1->face();

    // va->halfedge() = h0->next()->next();
    // check_for(va);
    // vb->halfedge() = h1->next();
    // vc->halfedge() = h0->next();
    // vd->halfedge() = h1->next()->next();

    // e0->halfedge() = h0;

    // f0->halfedge() = h0;
    // check_for(f0);
    // f1->halfedge() = h1;
    // check_for(f1);
    
    // h0->setNeighbors(h1->next()->next(), h1, va, e0, f0);

    // h1->setNeighbors(va->halfedge(), h0, vd, e0, f1);
    // check_for(h0);
    // check_for(h1);

    // h0->next()->next()->next() = h0;
    // h1->next()->next()->next() = h1;

    

    
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    if (e0 -> isBoundary()) {
      return newVertex();
    }

    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h5 = h0->twin();
    HalfedgeIter h3 = h5->next();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h6 = h2->twin();
    HalfedgeIter h7 = h1->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h3->twin();

    VertexIter v0 = h2->vertex();
    VertexIter v1 = h0->vertex();
    VertexIter v2 = h5->vertex();
    VertexIter v3 = h4->vertex();

    e0;
    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h3->edge();
    EdgeIter e4 = h4->edge();

    FaceIter f0 = h0->face();
    FaceIter f1 = h5->face();

    HalfedgeIter h10 = this->newHalfedge();
    HalfedgeIter h11 = this->newHalfedge();
    HalfedgeIter h12 = this->newHalfedge();
    HalfedgeIter h13 = this->newHalfedge();
    HalfedgeIter h14 = this->newHalfedge();
    HalfedgeIter h15 = this->newHalfedge();

    VertexIter v4 = this->newVertex();
    //compute position
    Vector3D a = v1->position;  
    Vector3D b = v2->position;
    Vector3D newPos = Vector3D();
    newPos.x = (a.x + b.x) / 2;
    newPos.y = (a.y + b.y) / 2;
    newPos.z = (a.z + b.z) / 2;

    v4->position = newPos;

    EdgeIter e5 = newEdge();
    EdgeIter e6 = newEdge();
    EdgeIter e7 = newEdge();

    FaceIter f2 = newFace();
    FaceIter f3 = newFace();



    //vfeh
    
    v0->halfedge() = h2;
    v1->halfedge() = h3;
    v2->halfedge() = h1;
    v3->halfedge() = h4;
    v4->halfedge() = h12;

    f0->halfedge() = h1;
    f1->halfedge() = h4;
    f2->halfedge() = h2;
    f3->halfedge() = h3;

    e0->halfedge() = h0;
     e0->isNew = true;
    e1->halfedge() = h1;
    e2->halfedge() = h2;
    e3->halfedge() = h3;
    e4->halfedge() = h4;
    e5->halfedge() = h13;
    e5->isNew = true;
    e5->isBlue = true;
    e6->halfedge() = h11;
    e6->isNew = true;
    e7->halfedge() = h10;
    e7->isBlue = true;
    e7->isNew = true;

    h0->setNeighbors(h1, h5, v4, e0, f0);
    h1->setNeighbors(h10, h7, v2, e1, f0);
    h2->setNeighbors(h11, h6, v0, e2, f2);
    h3->setNeighbors(h13, h9, v1, e3, f3);
    h4->setNeighbors(h5, h8, v3, e4, f1);
    h5->setNeighbors(h15, h0, v2, e0, f1);
    // h6->setNeighbors(h6->next(), h2, v1, e2, h6->face());
    // h7->setNeighbors(h7->next(), h1, v0, e1, h7->face());
    // h8->setNeighbors(h8->next(), h4, v2, e4, h8->face());
    // h9->setNeighbors(h9->next(), h3, v3, e3, h9->face());
    h10->setNeighbors(h0, h12, v0, e7, f0);
    h11->setNeighbors(h12, h14, v1, e6, f2);
    h12->setNeighbors(h2, h10, v4, e7, f2);
    h13->setNeighbors(h14, h15, v3, e5, f3);
    h14->setNeighbors(h3, h11, v4, e6, f3);
    h15->setNeighbors(h4, h13, v4, e5, f1);







    return v4;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.

    for(VertexIter v = mesh.verticesBegin(); v!= mesh.verticesEnd(); v++){
      if (!(v->isNew)){
        int n = v->degree();
        Vector3D sumPositions = Vector3D();
        HalfedgeCIter h = v->halfedge();
        do {
          sumPositions += h->twin()->vertex()->position;
          h = h->twin()->next();
          
        } while (h != v->halfedge());
        double u;
        if(n == 3){
          u = 0.1875;
        } else {
          u = 3.0/(8*n);
        }

        v->newPosition = ((1 - n * u) * v->position + u * sumPositions);

      }
    }
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      


      // e->newPosition = (e->halfedge()->vertex()->position + e->halfedge()->twin()->vertex()->position) / 2;
      Vector3D A = e->halfedge()->vertex()->position;
      Vector3D B = e->halfedge()->twin()->vertex()->position;
      Vector3D D = e->halfedge()->next()->next()->vertex()->position;
      Vector3D C = e->halfedge()->twin()->next()->next()->vertex()->position;

      Vector3D newPosition = Vector3D();
      newPosition.x = 0.375 * (A.x + B.x) + 0.125 * (C.x + D.x);
      newPosition.y = 0.375 * (A.y + B.y) + 0.125 * (C.y + D.y);
      newPosition.z = 0.375 * (A.z + B.z) + 0.125 * (C.z + D.z);
      e->newPosition = newPosition;
    }
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    int i = 0;
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        if (!(e->isNew)) {
            VertexIter v = mesh.splitEdge(e);
            v->isNew = true;
            // v->position = e->newPosition;
            v->newPosition = e->newPosition;
            i++;
            if (i % 1000 == 0){
              cout << i;
            }
        } 
    }
    //return;

    
    
    // 4. Flip any new edge that connects an old and new vertex.

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        // if (e->isNew && (e->halfedge()->vertex()->isNew ^ e->halfedge()->twin()->vertex()->isNew)) {
        //     mesh.flipEdge(e);
        // }
        if (e->isNew && (e->halfedge()->vertex()->isNew ^ e->halfedge()->twin()->vertex()->isNew) && e->isBlue) {
            mesh.flipEdge(e);
        }
        e->isNew = false;
        e->isBlue = false;
    }
    // return;
    // 5. Copy the new vertex positions into final Vertex::position.
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
        // if (v->isNew){
        v->position = v->newPosition;
        v->isNew = false;
        // }
    }



    // return;
  }
}
