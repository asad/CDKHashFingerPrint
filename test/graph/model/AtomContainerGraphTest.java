/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package graph.model;

import org.junit.*;

/**
 *
 * @author Asad
 */
public class AtomContainerGraphTest {

    public AtomContainerGraphTest() {
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
     * Test of main method, of class AtomContainerGraph.
     */
    @Test
    public void testMain() throws Exception {
        System.out.println("main");
        AtomContainerGraph g = new AtomContainerGraph();
        g.generateExampleGraphContainer();
    }
}
