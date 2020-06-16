"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                for (let i = 0; i < mesh.vertices.length; i++) {
                  mesh.vertices[i].index = i;
                }
                for (let i = 0; i < mesh.edges.length; i++) {
                  mesh.edges[i].index = i;
                }
                for (let i = 0; i < mesh.faces.length; i++) {
                  mesh.faces[i].index = i;
                }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                let m = new Triplet(mesh.edges.length, mesh.vertices.length);
                for (let edge of mesh.edges) {
                  m.addEntry(1, edge.index, edge.halfedge.vertex.index);
                  m.addEntry(1, edge.index, edge.halfedge.next.vertex.index);
                }
                return SparseMatrix.fromTriplet(m);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let m = new Triplet(mesh.faces.length, mesh.edges.length);
                for (let face of mesh.faces) {
                  for (let edge of face.adjacentEdges()) {
                    m.addEntry(1, face.index, edge.index);
                  }
                }
                return SparseMatrix.fromTriplet(m);
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let vector = DenseMatrix.zeros(this.mesh.vertices.length,);
                for (let vertex of subset.vertices) {
                  vector.set(1, vertex);
                }
                return vector;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let vector = DenseMatrix.zeros(this.mesh.edges.length);
                for (let edge of subset.edges) {
                  vector.set(1, vertex);
                }
                return vector;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let vector = DenseMatrix.zeros(this.mesh.faces.length);
                for (let face of subset.faces) {
                  vector.set(1, face);
                }
                return vector;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                let s = MeshSubset.deepCopy(subset);
                let star = new MeshSubset(s.vertices, s.edges, s.faces);
                for (let vertex of s.vertices) {
                  for (let edge of this.mesh.vertices[vertex].adjacentEdges()) {
                    star.addEdge(edge.index);
                  }
                  for (let face of this.mesh.vertices[vertex].adjacentFaces()) {
                    star.addFace(face.index);
                  }
                }
                for (let edge of s.edges) {
                  star.addFace(this.mesh.edges[edge].halfedge.face.index);
                  star.addFace(this.mesh.edges[edge].halfedge.twin.face.index);
                }

                return MeshSubset.deepCopy(star); // placeholder
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                let s = MeshSubset.deepCopy(subset);
                let closure = new MeshSubset(s.vertices, s.edges, s.faces);
                for (let face of s.faces) {
                  for (let edge of this.mesh.faces[face].adjacentEdges()) {
                    closure.addEdge(edge.index);
                  }
                  for (let vertex of this.mesh.faces[face].adjacentVertices()) {
                    closure.addVertex(vertex.index);
                  }
                }
                for (let edge of s.edges) {
                  closure.addVertex(this.mesh.edges[edge].halfedge.vertex.index);
                  closure.addVertex(this.mesh.edges[edge].halfedge.next.vertex.index);
                }

                return MeshSubset.deepCopy(closure);
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                let d1 = MeshSubset.deepCopy(subset);
                let d2 = MeshSubset.deepCopy(subset);
                let s1 = this.closure(this.star(d1));
                let s2 = this.star(this.closure(d2));
                let link = new MeshSubset();
                for (let vertex of s1.vertices) {
                  if (!s2.vertices.has(vertex)) {
                    link.addVertex(vertex);
                  }
                }
                for (let edge of s1.edges) {
                  if (!s2.edges.has(edge)) {
                    link.addEdge(edge);
                  }
                }
                for (let face of s1.faces) {
                  if (!s2.faces.has(face)) {
                    link.addFace(face);
                  }
                }
                return MeshSubset.deepCopy(link);
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                let closure = this.closure(subset);
                return closure.equals(subset);
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                let d = MeshSubset.deepCopy(subset);
                if (!this.isComplex(d)) {
                  return -1;
                }
                let satisfy = false;
                if (d.faces.size != 0) {
                  for (let face of d.faces) {
                    for (let edge of this.mesh.faces[face].adjacentEdges()) {
                      if (!d.edges.has(edge.index)) {
                        return -1;
                      }
                    }
                  }
                  for (let edge of d.edges) {
                    if (!(d.vertices.has(this.mesh.edges[edge].halfedge.vertex.index) && d.vertices.has(this.mesh.edges[edge].halfedge.next.vertex.index))) {
                      return -1;
                    }
                    if (!(d.faces.has(this.mesh.edges[edge].halfedge.face.index) || d.faces.has(this.mesh.edges[edge].halfedge.twin.face.index))) {
                      return -1;
                    }
                  }
                  for (let vertex of d.vertices) {
                    satisfy = false;
                    for (let edge of this.mesh.vertices[vertex].adjacentEdges()) {
                      if (d.edges.has(edge.index)) {
                        satisfy = true;
                      }
                    }
                    if (!satisfy) {
                      return -1;
                    }
                  }
                  return 2;
                }
                if (d.edges.size != 0) {
                  for (let edge of d.edges) {
                    if (!(d.vertices.has(this.mesh.edges[edge].halfedge.vertex.index) && d.vertices.has(this.mesh.edges[edge].halfedge.next.vertex.index))) {
                      return -1;
                    }
                  }
                  for (let vertex of d.vertices) {
                    satisfy = false;
                    for (let edge of this.mesh.vertices[vertex].adjacentEdges()) {
                      if (d.edges.has(edge.index)) {
                        satisfy = true;
                      }
                    }
                    if (!satisfy) {
                      return -1;
                    }
                  }
                  return 1;
                }
                if (d.vertices.size == 1) {
                  return 0;
                }
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let d = MeshSubset.deepCopy(subset);
                let k = this.isPureComplex(d);
                let boundary = new MeshSubset();
                let count = 0;
                if (k == 2) {
                  for (let edge of d.edges) {
                    count = 0;
                    if (d.faces.has(this.mesh.edges[edge].halfedge.face.index)) {
                      count += 1;
                    }
                    if (d.faces.has(this.mesh.edges[edge].halfedge.twin.face.index)) {
                      count += 1;
                    }
                    if (count == 1) {
                      boundary.addEdge(edge);
                      boundary.addVertex(this.mesh.edges[edge].halfedge.vertex.index);
                      boundary.addVertex(this.mesh.edges[edge].halfedge.twin.vertex.index);
                    }
                  }
                }
                if (k == 1) {
                  for (let vertex of d.vertices) {
                    count = 0;
                    for (let edge of this.mesh.vertices[vertex].adjacentEdges()) {
                      if (d.edges.has(edge.index)) {
                        count += 1;
                      }
                    }
                    if (count == 1) {
                      boundary.addVertex(vertex);
                    }
                  }
                }
                return boundary;
        }
}
