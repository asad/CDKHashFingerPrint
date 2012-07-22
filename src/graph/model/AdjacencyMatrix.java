/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph.model;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author Asad
 */
public class AdjacencyMatrix implements Serializable {

    private static final long serialVersionUID = 7688786252424151L;
    private Map<Integer, AtomVertex> lookupMap;
    private Map<AtomVertex, Integer> atomIndexMap;
    private int currentIndex;
    private Integer[][] matrix;

    public AdjacencyMatrix() {
    }

    public AdjacencyMatrix(int vertices) {
        super();
        this.matrix = new Integer[vertices][vertices];
        this.atomIndexMap = new HashMap<AtomVertex, Integer>();
        this.lookupMap = new HashMap<Integer, AtomVertex>();
        initializeMatrix();
    }

    /*
     * initialize the matrix without element as null
     */
    private void initializeMatrix() {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                matrix[i][j] = null;
            }
        }
    }
    
    /**
     * Add an unweighted edge
     * @param source
     * @param sink
     */
    public void addEdge(AtomVertex source, AtomVertex sink) {
        addEdge(source, sink, 0);
    }
    
    /**
     * Add a weighted edge
     * @param source
     * @param sink
     * @param weight
     */
    public void addEdge(AtomVertex source, AtomVertex sink, int weight) {
        Integer i = atomIndexMap.get(source);
        Integer j = atomIndexMap.get(sink);

        if (i == null || j == null || i > matrix.length - 1
                || j > matrix.length - 1) {
            return;
        }

        this.matrix[i][j] = Integer.valueOf(weight);
        this.lookupMap.put(j, sink);
    }

    /*
     * Get the adjacent vertices of the query vertex
     *
     *
     * @param queryVertex atom vertex @return
     */
    public Map<AtomVertex, Integer> getAdjacentVertices(AtomVertex queryVertex) {
        Integer i = atomIndexMap.get(queryVertex);
        return i == null ? null : getAdjacentVertices(i);
    }

    private Map<AtomVertex, Integer> getAdjacentVertices(Integer i) {
        Map<AtomVertex, Integer> adjacentAtomVertexMap = new HashMap<AtomVertex, Integer>();
        for (int j = 0; j < matrix.length; j++) {
            if (matrix[i][j] != null) {
                adjacentAtomVertexMap.put(lookupMap.get(j), matrix[i][j]);
            }
        }
        return adjacentAtomVertexMap;
    }

    /**
     * Add a vertex to the adjacency matrix
     *
     * @param atomVertex
     * @return previous value associated with the key 
     */
    public Integer addVertex(AtomVertex atomVertex) {
        Integer value = this.atomIndexMap.put(atomVertex, currentIndex);
        currentIndex += 1;
        return value;
    }
}
