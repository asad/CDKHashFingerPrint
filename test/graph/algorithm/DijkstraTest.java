/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph.algorithm;

import graph.model.AtomContainerGraph;
import org.junit.*;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Asad
 */
public class DijkstraTest {

    public DijkstraTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of main method, of class Dijkstra.
     */
    @Test
    public void testMain() throws Exception {
        System.out.println("main");
        IAtomContainer atomContainer = new AtomContainer();
        IAtom atom1 = new Atom("C");
        IAtom atom2 = new Atom("N");
        IAtom atom3 = new Atom("O");
        IAtom atom4 = new Atom("S");
        IAtom atom5 = new Atom("C");

        IBond bond1 = new Bond(atom1, atom2, IBond.Order.SINGLE);
        IBond bond2 = new Bond(atom1, atom3, IBond.Order.SINGLE);
        IBond bond3 = new Bond(atom1, atom4, IBond.Order.SINGLE);
        IBond bond4 = new Bond(atom2, atom5, IBond.Order.SINGLE);

        atomContainer.addAtom(atom1);
        atomContainer.addAtom(atom2);
        atomContainer.addAtom(atom3);
        atomContainer.addAtom(atom4);
        atomContainer.addAtom(atom5);


        atomContainer.addBond(bond1);
        atomContainer.addBond(bond2);
        atomContainer.addBond(bond3);
        atomContainer.addBond(bond4);

        AtomContainerGraph g = new AtomContainerGraph(atomContainer, true);
        Dijkstra dijkstra = new Dijkstra(g, g.getVertexLookupMap().get(atom1));
        dijkstra.printShortestPath(dijkstra.getSinkShorestPath(g.getVertexLookupMap().get(atom5)));

    }
}
