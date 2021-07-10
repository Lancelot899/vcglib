/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef __VCG_TRI_UPDATE_TOPOLOGY
#define __VCG_TRI_UPDATE_TOPOLOGY
#include <algorithm>
#include <vector>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/simplex/edge/topology.h>

namespace vcg {
namespace tri {
/// \ingroup trimesh

/// \headerfile topology.h vcg/complex/algorithms/update/topology.h

/// \brief Generation of per-vertex and per-face topological information.

template <class UpdateMeshType>
class UpdateTopology
{

public:
typedef UpdateMeshType MeshType;
typedef typename MeshType::ScalarType     ScalarType;
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::EdgeType       EdgeType;
typedef typename MeshType::EdgePointer    EdgePointer;
typedef typename MeshType::EdgeIterator   EdgeIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;


/// \headerfile topology.h vcg/complex/algorithms/update/topology.h

/// \brief Auxiliairy data structure for computing face face adjacency information.
/**
It identifies and edge storing two vertex pointer and a face pointer where it belong.
*/

class PEdge
{
public:

  VertexPointer  v[2];  // the two Vertex pointer are ordered!
  FacePointer    f;     // the face where this edge belong
  int            z;     // index in [0..2] of the edge of the face
  bool isBorder;

  PEdge() {}
  PEdge(FacePointer  pf, const int nz) { this->Set(pf,nz); }
  void Set( FacePointer  pf, const int nz )
  {
    assert(pf!=0);
    assert(nz>=0);
    assert(nz<pf->VN());

    v[0] = pf->V(nz);
    v[1] = pf->V(pf->Next(nz)); // pf->Next得到的是下一个指标的点
    assert(v[0] != v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

    if( v[0] > v[1] ) std::swap(v[0],v[1]);
    f    = pf;
    z    = nz;
  }

  inline bool operator <  ( const PEdge & pe ) const
  {
    if( v[0]<pe.v[0] ) return true;
    else if( v[0]>pe.v[0] ) return false;
    else return v[1] < pe.v[1];
  }

  inline bool operator == ( const PEdge & pe ) const
  {
    return v[0]==pe.v[0] && v[1]==pe.v[1];
  }
  /// Convert from edge barycentric coord to the face baricentric coord a point on the current edge.
  /// Face barycentric coordinates are relative to the edge face.
  inline Point3<ScalarType> EdgeBarycentricToFaceBarycentric(ScalarType u) const
  {
    Point3<ScalarType> interp(0,0,0);
    interp[ this->z     ] = u;
    interp[(this->z+1)%3] = 1.0f-u;
    return interp;
  }
};

/// Fill a vector with all the edges of the mesh.
/// each edge is stored in the vector the number of times that it appears in the mesh, with the referring face.
/// optionally it can skip the faux edges (to retrieve only the real edges of a triangulated polygonal mesh)
static void FillEdgeVector(MeshType &m, std::vector<PEdge> &edgeVec, bool includeFauxEdge=true) {
  edgeVec.reserve(m.fn*3);
  for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
          for (int j = 0; j < (*fi).VN(); ++j) {
              if (includeFauxEdge || !(*fi).IsF(j)) {
                  edgeVec.push_back(PEdge(&*fi, j));
              }
          }
      }
  }
}

static void FillUniqueEdgeVector(MeshType &m, std::vector<PEdge> &edgeVec, bool includeFauxEdge=true, bool computeBorderFlag=false)
{
    FillEdgeVector(m,edgeVec,includeFauxEdge);
    sort(edgeVec.begin(), edgeVec.end()); // oredering by vertex

    if (computeBorderFlag) {
        for (size_t i=0; i<edgeVec.size(); i++)
            edgeVec[ i ].isBorder = true;
        for (size_t i=1; i<edgeVec.size(); i++) {
            if (edgeVec[i]==edgeVec[i-1])
                edgeVec[i-1].isBorder = edgeVec[i-1].isBorder = false;
        }
    }

    typename std::vector< PEdge>::iterator newEnd = std::unique(edgeVec.begin(), edgeVec.end());

    edgeVec.resize(newEnd-edgeVec.begin()); // redundant! remove?
}


/*! \brief Initialize the edge vector all the edges that can be inferred from current face vector, setting up all the current adjacency relations
 *
 *
 */

static void AllocateEdge(MeshType &m)
{
  // Delete all the edges (if any)
  for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
        tri::Allocator<MeshType>::DeleteEdge(m,*ei);
  tri::Allocator<MeshType>::CompactEdgeVector(m);

  // Compute and add edges
  std::vector<PEdge> Edges;
  FillUniqueEdgeVector(m,Edges,true,tri::HasPerEdgeFlags(m) );
  assert(m.edge.empty());
  tri::Allocator<MeshType>::AddEdges(m,Edges.size());
  assert(m.edge.size()==Edges.size());

  // Setup adjacency relations
  if(tri::HasEVAdjacency(m))
  {
    for(size_t i=0; i< Edges.size(); ++i)
    {
      m.edge[i].V(0) = Edges[i].v[0];
      m.edge[i].V(1) = Edges[i].v[1];
    }
  }

  if (tri::HasPerEdgeFlags(m)){
    for(size_t i=0; i< Edges.size(); ++i) {
        if (Edges[i].isBorder) m.edge[i].SetB(); else m.edge[i].ClearB();
    }
  }

  if(tri::HasEFAdjacency(m)) // Note it is an unordered relation.
  {
    for(size_t i=0; i< Edges.size(); ++i)
    {
      std::vector<FacePointer> fpVec;
      std::vector<int> eiVec;
      face::EFStarFF(Edges[i].f,Edges[i].z,fpVec,eiVec);
      m.edge[i].EFp() = Edges[i].f;
      m.edge[i].EFi() = Edges[i].z;
    }
  }

  if(tri::HasFEAdjacency(m))
  {
    for(size_t i=0; i< Edges.size(); ++i)
    {
      std::vector<FacePointer> fpVec;
      std::vector<int> eiVec;
      face::EFStarFF(Edges[i].f,Edges[i].z,fpVec,eiVec);
      for(size_t j=0;j<fpVec.size();++j)
        fpVec[j]->FEp(eiVec[j])=&(m.edge[i]);

//      Edges[i].f->FE(Edges[i].z) = &(m.edge[i]);
//      Connect in loop the non manifold
//      FaceType* fpit=fp;
//      int eit=ei;

//      do
//      {
//        faceVec.push_back(fpit);
//        indVed.push_back(eit);
//        FaceType *new_fpit = fpit->FFp(eit);
//        int       new_eit  = fpit->FFi(eit);
//        fpit=new_fpit;
//        eit=new_eit;
//      } while(fpit != fp);


//      m.edge[i].EFp() = Edges[i].f;
//      m.edge[i].EFi() = ;
    }
  }

}

/// \brief Clear the Face-Face topological relation setting each involved pointer to null.
/// useful when you passed a mesh with ff adjacency to an algorithm that does not use it and could have messed it.
static void ClearFaceFace(MeshType &m)
{
  RequireFFAdjacency(m);
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
  {
    if( ! (*fi).IsD() )
    {
      for(int j=0;j<fi->VN();++j)
      {
        fi->FFp(j)=0;
        fi->FFi(j)=-1;
      }
    }
  }
}

/// \brief Update the Face-Face topological relation by allowing to retrieve for each face what other faces shares their edges.
static void FaceFace(MeshType &m)
{
  RequireFFAdjacency(m);
  if( m.fn == 0 ) return;

  std::vector<PEdge> e;
  FillEdgeVector(m,e); // 这一步会得到所有的边，每个边里面记录了face，同一个边(点相同)记录的face一定不一样
  sort(e.begin(), e.end());							// Lo ordino per vertici
  ///   /_\ 
  /// /_\/_\
  /// 三角形： f2
  ///      f3 f0 f1
  /// f0(513)  f1(312) f2(453) f3(501)
  /// 获得的PEdge: e(f3, 0, 1, 1) e(f3, 0, 5, 0) e(f3, 1, 2, 1) e(f0, 1, 3, 1) e(f3, 1, 3, 0) e(f0, 1,5, 0)
  ///             e(f3, 1, 5, 2) e(f1, 2, 3, 2) e(f2, 3, 4, 2) e(f0, 3, 5, 2) e(f2, 3, 5, 1) e(f2, 4, 5, 0)
  /// 第一次循环，会操作e(f3, 0, 1, 1)，由于没有和它相同的pE，因此设置f3.FFp(1) = f3
  /// 第二次循环，会操作e(f3, 0, 5, 0)，由于没有和它相同的pE，因此设置f3.FFp(0) = f3
  /// 第三次循环，会操作e(f0, 1, 3, 1)，由于没有和它和e(f3, 1, 3, 0)相同，因此设置f0.FFp(1) = f3, f3.FFp(0) = f0
  /// 最终会得到表格
  /// f0: f3, f1, f2
  /// f1: f0, f1, f1
  /// f2: f2. f0, f2
  /// f3: f3, f3, f0
  /// 正好如果有共边的，就设置为对面的面id，没有共边的就设置为自己的面id
  int ne = 0;											// Numero di edge reali

  typename std::vector<PEdge>::iterator pe,ps;
  ps = e.begin();pe=e.begin();
  //for(ps = e.begin(),pe=e.begin();pe<=e.end();++pe)	// Scansione vettore ausiliario
  do
  {
    if( pe==e.end() || !(*pe == *ps) )					// Trovo blocco di edge uguali
    {
        //这里生成vertexFace类似的链表
      typename std::vector<PEdge>::iterator q,q_next;
      for (q=ps;q<pe-1;++q)						// Scansione facce associate
      {
        assert((*q).z>=0);
        //assert((*q).z< 3);
        q_next = q;
        ++q_next;
        assert((*q_next).z>=0);
        assert((*q_next).z< (*q_next).f->VN());
        (*q).f->FFp(q->z) = (*q_next).f;				// Collegamento in lista delle facce
        (*q).f->FFi(q->z) = (*q_next).z;
      }
      assert((*q).z>=0);
      assert((*q).z< (*q).f->VN());
      (*q).f->FFp((*q).z) = ps->f;
      (*q).f->FFi((*q).z) = ps->z;
      ps = pe;
      ++ne;										// Aggiorno il numero di edge
    }
    if(pe==e.end()) break;
    ++pe;
  } while(true);
}

/// \brief Update the Vertex-Face topological relation.
/**
The function allows to retrieve for each vertex the list of faces sharing this vertex.
After this call all the VF component are initialized. Isolated vertices have a null list of faces.
\sa vcg::vertex::VFAdj
\sa vcg::face::VFAdj
*/

static void VertexFace(MeshType &m)
{
  RequireVFAdjacency(m);

  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    (*vi).VFp() = 0;
    (*vi).VFi() = 0; // note that (0,-1) means uninitiazlied while 0,0 is the valid initialized values for isolated vertices.
  }

  for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD())
      {
          for (int j = 0; j < (*fi).VN(); ++j)
          {
              (*fi).VFp(j) = (*fi).V(j)->VFp(); //首先会遍历面的所有点，然后记录点指向的面
              (*fi).VFi(j) = (*fi).V(j)->VFi(); // 这里记录那个面的对应点是第几个点
              (*fi).V(j)->VFp() = &(*fi); //把这个点所只想的面记录当前面，用于下一次迭代使用
              (*fi).V(j)->VFi() = j;
          }
      }
  }

  /// 这个过程其实是通过点存储的VFp和VFi将邻域信息记录到面的VFp和VFi中，具体流程是
  /// __ __ __  如图所示，面从左到右为0,1,2,3,4
  /// \/_\/_\/  点为012 213 324 435 456
  /// 第一次迭代时候，0号面记录的点节点为xxx(x代表空指针)，然后012三个点的VFp变成了000
  /// 第二次迭代，1号面指针变为了00x，然后213三个点的VFP变成了111
  /// 第三次迭代，2号面指针变为了11x, 然后324三个点的VFP变成了222
  /// 第四次迭代，3号面指针变为了22x, 然后435三个点的VFp变为333
  /// 第五次迭代，4号面指针变为了33x，然后456三个点VFp变为444
  /// 最后统计一下点，(0, 0) (1,1) (2,2) (3,3) (4, 4) (5, 4) (6, 4)
  /// 现在考虑迭代2号点周围的面,首先会拿出2号点对应的面，即2号面，然后2号面记录的2号点对应的面VFp指针为1，
  /// 这时候会迭代1号面，1号面对应的2号点对应的指针为0，这时会迭代0号面，最后0号面记录了x，终止迭代
}


/// \headerfile topology.h vcg/complex/algorithms/update/topology.h

/// \brief Auxiliairy data structure for computing face face adjacency information.
/**
It identifies and edge storing two vertex pointer and a face pointer where it belong.
*/

class PEdgeTex
{
public:

  typename FaceType::TexCoordType  v[2];		// the two TexCoord are ordered!
  FacePointer    f;                       // the face where this edge belong
  int      z;				      // index in [0..2] of the edge of the face

  PEdgeTex() {}

  void Set( FacePointer  pf, const int nz )
  {
    assert(pf!=0);
    assert(nz>=0);
    assert(nz<3);

    v[0] = pf->WT(nz);
    v[1] = pf->WT(pf->Next(nz));
    assert(v[0] != v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

    if( v[1] < v[0] ) std::swap(v[0],v[1]);
    f    = pf;
    z    = nz;
  }

  inline bool operator <  ( const PEdgeTex & pe ) const
  {
    if( v[0]<pe.v[0] ) return true;
    else if( pe.v[0]<v[0] ) return false;
    else return v[1] < pe.v[1];
  }
  inline bool operator == ( const PEdgeTex & pe ) const
  {
    return (v[0]==pe.v[0]) && (v[1]==pe.v[1]);
  }
  inline bool operator != ( const PEdgeTex & pe ) const
  {
    return (v[0]!=pe.v[0]) || (v[1]!=pe.v[1]);
  }

};


/// \brief Update the Face-Face topological relation so that it reflects the per-wedge texture connectivity

/**
Using this function two faces are adjacent along the FF relation IFF the two faces have matching texture coords along the involved edge.
In other words F1->FFp(i) == F2 iff F1 and F2 have the same tex coords along edge i
*/

static void FaceFaceFromTexCoord(MeshType &m)
{
  RequireFFAdjacency(m);
  RequirePerFaceWedgeTexCoord(m);
  vcg::tri::UpdateTopology<MeshType>::FaceFace(m);
  for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
  {
    if (!(*fi).IsD())
    {
      for (int i = 0; i < (*fi).VN(); i++)
      {
        if (!vcg::face::IsBorder((*fi), i))
        {
          typename MeshType::FacePointer nextFace = (*fi).FFp(i);
          int nextEdgeIndex = (*fi).FFi(i);
          bool border = false;
          if ((*fi).cV(i) == nextFace->cV(nextEdgeIndex))
          {
            if ((*fi).WT(i) != nextFace->WT(nextEdgeIndex) || (*fi).WT((*fi).Next(i)) != nextFace->WT(nextFace->Next(nextEdgeIndex)))
              border = true;
          }
          else
          {
            if ((*fi).WT(i) != nextFace->WT(nextFace->Next(nextEdgeIndex)) || (*fi).WT((*fi).Next(i)) != nextFace->WT(nextEdgeIndex))
              border = true;
          }
          if (border)
            vcg::face::FFDetach((*fi), i);

        }
      }
    }
  }
}

/// \brief Test correctness of VEtopology
static void TestVertexEdge(MeshType &m)
{
  std::vector<int> numVertex(m.vert.size(),0);
  
  tri::RequireVEAdjacency(m);
  
  for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
  {
      if (!(*ei).IsD())
      {
        assert(tri::IsValidPointer(m,ei->V(0)));
        assert(tri::IsValidPointer(m,ei->V(1)));
        if(ei->VEp(0)) assert(tri::IsValidPointer(m,ei->VEp(0)));
        if(ei->VEp(1)) assert(tri::IsValidPointer(m,ei->VEp(1)));
        numVertex[tri::Index(m,(*ei).V(0))]++;
        numVertex[tri::Index(m,(*ei).V(1))]++;
      }
  }
  
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
      if (!vi->IsD())
      {
        int cnt =0;
        int ind = tri::Index(m,*vi);
        int incidentNum = numVertex[ind];
        for(edge::VEIterator<EdgeType> vei(&*vi);!vei.End();++vei)
          cnt++;
        EdgeType *vep = vi->VEp();
        assert((incidentNum==0) == (vi->VEp()==0) );
        assert(cnt==incidentNum);        
      }
  }  
}


/// \brief Test correctness of VFtopology
static void TestVertexFace(MeshType &m)
{
    SimpleTempData<typename MeshType::VertContainer, int > numVertex(m.vert,0);

  assert(tri::HasPerVertexVFAdjacency(m));

    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
        if (!(*fi).IsD())
        {
            numVertex[(*fi).V0(0)]++;
            numVertex[(*fi).V1(0)]++;
            numVertex[(*fi).V2(0)]++;
        }
    }

    vcg::face::VFIterator<FaceType> VFi;

    for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
    {
        if (!vi->IsD())
        if(vi->VFp()!=0) // unreferenced vertices MUST have VF == 0;
        {
            int num=0;
            assert(tri::IsValidPointer(m, vi->VFp()));
            VFi.f=vi->VFp();
            VFi.z=vi->VFi();
            while (!VFi.End())
            {
                num++;
                assert(!VFi.F()->IsD());
                assert((VFi.F()->V(VFi.I()))==&(*vi));
                ++VFi;
            }
            assert(num==numVertex[&(*vi)]);
        }
    }
}

/// \brief Test correctness of FFtopology (only for 2Manifold Meshes!)
static void TestFaceFace(MeshType &m)
{
  assert(HasFFAdjacency(m));

  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
    if (!fi->IsD())
        {
      for (int i=0;i<(*fi).VN();i++)
            {
        FaceType *ffpi=fi->FFp(i);
        int e=fi->FFi(i);
        //invariant property of FF topology for two manifold meshes
        assert(ffpi->FFp(e) == &(*fi));
        assert(ffpi->FFi(e) == i);

        // Test that the two faces shares the same edge
        // Vertices of the i-th edges of the first face
        VertexPointer v0i= fi->V0(i);
        VertexPointer v1i= fi->V1(i);
        // Vertices of the corresponding edge on the other face
        VertexPointer ffv0i= ffpi->V0(e);
        VertexPointer ffv1i= ffpi->V1(e);

        assert( (ffv0i==v0i) || (ffv0i==v1i) );
        assert( (ffv1i==v0i) || (ffv1i==v1i) );
            }

        }
    }
}

/// Auxiliairy data structure for computing edge edge adjacency information.
/// It identifies an edge storing a vertex pointer and a edge pointer where it belong.
class PVertexEdge
{
public:

  VertexPointer  v;		// the two Vertex pointer are ordered!
  EdgePointer    e;		  // the edge where this vertex belong
  int      z;				      // index in [0..1] of the vertex of the edge

  PVertexEdge(  ) {}
  PVertexEdge( EdgePointer  pe, const int nz )
{
  assert(pe!=0);
  assert(nz>=0);
  assert(nz<2);

  v= pe->V(nz);
  e    = pe;
  z    = nz;
}
inline bool operator  <  ( const PVertexEdge & pe ) const { return ( v<pe.v ); }
inline bool operator ==  ( const PVertexEdge & pe ) const { return ( v==pe.v ); }
inline bool operator !=  ( const PVertexEdge & pe ) const { return ( v!=pe.v ); }
};



static void EdgeEdge(MeshType &m)
{
  RequireEEAdjacency(m);
  std::vector<PVertexEdge> v;
  if( m.en == 0 ) return;

//  printf("Inserting Edges\n");
  for(EdgeIterator pf=m.edge.begin(); pf!=m.edge.end(); ++pf)			// Lo riempio con i dati delle facce
    if( ! (*pf).IsD() )
      for(int j=0;j<2;++j)
      {
//        printf("egde %i ind %i (%i %i)\n",tri::Index(m,&*pf),j,tri::Index(m,pf->V(0)),tri::Index(m,pf->V(1)));
        v.push_back(PVertexEdge(&*pf,j));
      }

//  printf("en = %i (%i)\n",m.en,m.edge.size());
  sort(v.begin(), v.end());							// Lo ordino per vertici

  int ne = 0;											// Numero di edge reali

  typename std::vector<PVertexEdge>::iterator pe,ps;
  // for(ps = v.begin(),pe=v.begin();pe<=v.end();++pe)	// Scansione vettore ausiliario
  ps = v.begin();pe=v.begin();
  do
  {
//    printf("v %i -> e %i\n",tri::Index(m,(*ps).v),tri::Index(m,(*ps).e));
    if( pe==v.end() || !(*pe == *ps) )					// Trovo blocco di edge uguali
    {
      typename std::vector<PVertexEdge>::iterator q,q_next;
      for (q=ps;q<pe-1;++q)						// Scansione edge associati
      {
        assert((*q).z>=0);
        assert((*q).z< 2);
        q_next = q;
        ++q_next;
        assert((*q_next).z>=0);
        assert((*q_next).z< 2);
        (*q).e->EEp(q->z) = (*q_next).e;				// Collegamento in lista delle facce
        (*q).e->EEi(q->z) = (*q_next).z;
      }
      assert((*q).z>=0);
      assert((*q).z< 2);
      (*q).e->EEp((*q).z) = ps->e;
      (*q).e->EEi((*q).z) = ps->z;
      ps = pe;
      ++ne;										// Aggiorno il numero di edge
    }
    if(pe==v.end()) break;
    ++pe;
   } while(true);
}

static void VertexEdge(MeshType &m)
{
  RequireVEAdjacency(m);

  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
  {
    (*vi).VEp() = 0;
    (*vi).VEi() = 0;
  }

  for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
  if( ! (*ei).IsD() )
  {
    for(int j=0;j<2;++j)
    { assert(tri::IsValidPointer(m,ei->V(j)));
      (*ei).VEp(j) = (*ei).V(j)->VEp();
      (*ei).VEi(j) = (*ei).V(j)->VEi();
      (*ei).V(j)->VEp() = &(*ei);
      (*ei).V(j)->VEi() = j;
    }
  }
}

}; // end class

}	// End namespace
}	// End namespace


#endif
