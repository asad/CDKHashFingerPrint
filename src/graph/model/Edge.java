package graph.model;

public class Edge {

    private AtomVertex toVertex;
    private int weight;

    public Edge() {
    }

    public Edge(AtomVertex sinkVertex) {
        super();
        this.toVertex = sinkVertex;
    }

    public Edge(AtomVertex to, int weight) {
        super();
        this.toVertex = to;
        this.weight = weight;
    }

    @Override
    public String toString() {
        return "Edge [sink=" + getSinkVertex() + ", weight=" + weight + "]";
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result
                + ((getSinkVertex() == null) ? 0 : getSinkVertex().hashCode());
        result = prime * result + weight;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        Edge other = (Edge) obj;
        if (getSinkVertex() == null) {
            if (other.getSinkVertex() != null) {
                return false;
            }
        } else if (!toVertex.equals(other.getSinkVertex())) {
            return false;
        }
        if (weight != other.getWeight()) {
            return false;
        }
        return true;
    }

    public int getWeight() {
        return weight;
    }

    /**
     * @return the toVertex
     */
    public AtomVertex getSinkVertex() {
        return toVertex;
    }
}
