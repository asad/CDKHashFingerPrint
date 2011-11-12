/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package isomorphism;

import fingerprints.Fingerprinter;
import fingerprints.IFingerprinter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.BitSet;
import junit.framework.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.algorithm.vflib.substructure.VF2;

/**
 *
 * @author Asad
 */
public class SubgraphTest {

    public static void main(String[] args) throws InvalidSmilesException, Exception {
        testIsSubgraph1();
        testIsSubgraph2();
    }

    public SubgraphTest() {
    }

    /**
     * Test of Fingerprinter method
     * @throws InvalidSmilesException
     * @throws CDKException  
     */
    @Test
    public static void testIsSubgraph1() throws InvalidSmilesException, CDKException {

        String smiles = "CCCCC1C(=O)N(N(C1=O)C1=CC=CC=C1)C1=CC=CC=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new Fingerprinter(1024);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getFingerprint(molecule);
        System.out.println("fp " + fingerprint1.toString());
    }

    /**
     * Test of Fingerprinter method
     * @throws InvalidSmilesException
     * @throws CDKException
     * @throws FileNotFoundException  
     */
    @Test
    public static void testIsSubgraph2() throws InvalidSmilesException, CDKException, FileNotFoundException {


        FileReader smilesQ =
                new FileReader(System.getProperty("user.home") + File.separator + "Software/GITROOT/Fingerprint/test/data/mol/C00107.mol");
        FileReader smilesT = new FileReader(System.getProperty("user.home") + File.separator + "Software/GITROOT/Fingerprint/test/data/mol/C00196.mol");
        MDLV2000Reader readerQ = new MDLV2000Reader(smilesQ);
        MDLV2000Reader readerT = new MDLV2000Reader(smilesT);
        IAtomContainer moleculeQ = (IAtomContainer) readerQ.read(new Molecule());
        IAtomContainer moleculeT = (IAtomContainer) readerT.read(new Molecule());
        moleculeQ.setID("C00107");
        moleculeT.setID("C00196");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        moleculeQ = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        moleculeT = AtomContainerManipulator.removeHydrogens(moleculeT);

        boolean uit = UniversalIsomorphismTester.isSubgraph(moleculeT, moleculeQ);
        boolean vf2 = new VF2(moleculeQ, moleculeT, true, false).isSubgraph();

        System.out.println("isSubset: " + vf2);

        Assert.assertFalse(vf2);
        Assert.assertTrue(uit);
    }
}
