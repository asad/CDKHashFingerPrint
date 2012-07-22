/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph.model;

import java.io.Serializable;
import org.openscience.cdk.interfaces.IAtom;

/**
 *
 * @author Asad
 */
class Vertex implements
        Cloneable, Serializable, Comparable<AtomVertex> {

    private static final long serialVersionUID = 786786786110723560L;
    private int vertexId;

    public Vertex() {
        super();
    }

    public void setVertexId(int vertexId) {
        this.vertexId = vertexId;
    }

    public int getVertexId() {
        return vertexId;
    }

    @Override
    public int compareTo(AtomVertex t) {
        if (this.getVertexId() > t.getVertexId()) {
            return 1;
        } else if (this.getVertexId() < t.getVertexId()) {
            return -1;
        } else {
            return 0;
        }
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Vertex other = (Vertex) obj;
        if (this.vertexId != other.getVertexId()) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 37 * hash + this.getVertexId();
        return hash;
    }
}

public class AtomVertex extends Vertex {

    private static final long serialVersionUID = 136767675678688282L;
    private IAtom atom;

    public AtomVertex() {
    }

    @Override
    public int compareTo(AtomVertex o) {
        return this.getAtom().getSymbol().compareTo(o.getAtom().getSymbol());
    }

    public AtomVertex(IAtom name) {
        super();
        this.setAtom(name);
    }

    public AtomVertex(IAtom name, int id) {
        super();
        this.setAtom(name);
        this.setVertexId(id);
    }

    @Override
    public AtomVertex clone() {
        return new AtomVertex(this.getAtom(), getVertexId());
    }

    @Override
    public String toString() {
        return "Vertex {" + getAtom().getSymbol() + ", " + getVertexId() + "}";
    }

    public void setAtom(IAtom name) {
        this.atom = name;
    }

    public IAtom getAtom() {
        return atom;
    }
}
