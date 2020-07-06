"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		let m = DenseMatrix.zeros(geometry.mesh.vertices.length);
		for (let v of geometry.mesh.vertices) {
			m.set(geometry.barycentricDualArea(v), vertexIndex[v]);
		}
		return SparseMatrix.diag(m);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		let v = DenseMatrix.zeros(geometry.mesh.edges.length);
		for (let e of geometry.mesh.edges) {
			v.set(0.5*(geometry.cotan(e.halfedge)+geometry.cotan(e.halfedge.twin)), edgeIndex[e]);
		}
		return SparseMatrix.diag(v);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		let v = DenseMatrix.zeros(geometry.mesh.faces.length);
		for (let f of geometry.mesh.faces) {
			v.set(1/geometry.area(f), faceIndex[f]);
		}
		return SparseMatrix.diag(v);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		let m = new Triplet(geometry.mesh.edges.length, geometry.mesh.vertices.length);
		for (let e of geometry.mesh.edges) {
			m.addEntry(-1, edgeIndex[e], vertexIndex[e.halfedge.vertex]);
			m.addEntry(1, edgeIndex[e], vertexIndex[e.halfedge.twin.vertex]);
		}
		return SparseMatrix.fromTriplet(m);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		let m = new Triplet(geometry.mesh.faces.length, geometry.mesh.edges.length);
		for (let e of geometry.mesh.edges) {
			m.addEntry(1, faceIndex[e.halfedge.face], edgeIndex[e]);
			m.addEntry(-1, faceIndex[e.halfedge.twin.face], edgeIndex[e]);
		}
		return SparseMatrix.fromTriplet(m);
	}
}
