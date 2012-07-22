/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph.model;

import java.util.Collections;
import java.util.Map;

/**
 *
 * @author Asad
 */
public class ShortestPathContainer {

    private final AtomContainerGraph atomContainerGraph;
    private final AtomVertex start;
    private final Map<AtomVertex, Integer> path;
    private final Map<AtomVertex, AtomVertex> ancestors;

    /**
     *
     * @param path
     * @param ancestors
     * @param atomContainerGraph
     * @param start
     */
    public ShortestPathContainer(
            Map<AtomVertex, Integer> path,
            Map<AtomVertex, AtomVertex> ancestors,
            AtomContainerGraph atomContainerGraph,
            AtomVertex start) {
        super();
        this.path = path;
        this.ancestors = ancestors;
        this.atomContainerGraph = atomContainerGraph;
        this.start = start;
    }

    /**
     * @return the atomContainerGraph
     */
    public AtomContainerGraph getG() {
        return atomContainerGraph;
    }

    /**
     * @return the start
     */
    public AtomVertex getStart() {
        return start;
    }

    /**
     * @return the path
     */
    public Map<AtomVertex, Integer> getPath() {
        return Collections.unmodifiableMap(path);
    }

    /**
     * @return the ancestors
     */
    public Map<AtomVertex, AtomVertex> getAncestors() {
        return Collections.unmodifiableMap(ancestors);
    }
}
