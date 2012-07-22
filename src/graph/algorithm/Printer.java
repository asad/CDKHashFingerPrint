/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph.algorithm;

import graph.model.AtomVertex;
import graph.model.ShortestPathContainer;
import java.util.Map;
import java.util.Stack;

/**
 *
 * @author Asad
 */
class Printer {

    protected void printShortestPaths(ShortestPathContainer e, boolean withCosts) {
        for (Map.Entry<AtomVertex, AtomVertex> entry : e.getAncestors().entrySet()) {
            if (withCosts) {
                System.out.println("sink " + entry.getKey().toString() + " cost: "
                        + e.getPath().get(entry.getKey()));
            } else {
                System.out.println("\t" + entry.getKey());
            }
        }
    }

    protected void printShortestPath(Stack<AtomVertex> pathStack) {
        while (!pathStack.isEmpty()) {
            System.out.println("\t" + pathStack.pop());
        }
    }
}
