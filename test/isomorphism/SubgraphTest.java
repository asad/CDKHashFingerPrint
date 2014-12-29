package isomorphism;

import fingerprints.hashed.HashedFingerprinter;
import fingerprints.interfaces.IFingerprinter;
import java.io.FileNotFoundException;
import java.util.BitSet;
import junit.framework.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.algorithm.vflib.substructure.VF2;

/**
 *
 * @author Asad
 */
public class SubgraphTest {

    final static SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

    public static void main(String[] args) throws InvalidSmilesException, Exception {
        testIsSubgraph1();
        testIsSubgraph2();
    }

    public SubgraphTest() {
    }

    /**
     * Test of Fingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public static void testIsSubgraph1() throws InvalidSmilesException, CDKException {

        String smiles = "CCCCC1C(=O)N(N(C1=O)C1=CC=CC=C1)C1=CC=CC=C1";
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new HashedFingerprinter(1024);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    /**
     * Test of Fingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     * @throws FileNotFoundException
     */
    @Test
    public static void testIsSubgraph2() throws InvalidSmilesException, CDKException, FileNotFoundException {
        String smilesq = "NC([*])C(=O)NC([*])C(O)=O";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        String smilest = "OC(=O)C1=CC=CC(O)=C1O";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);

        moleculeQ.setID("C00107");
        moleculeT.setID("C00196");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        moleculeQ = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        moleculeT = AtomContainerManipulator.removeHydrogens(moleculeT);

        boolean uit = new UniversalIsomorphismTester().isSubgraph(moleculeT, moleculeQ);
        boolean vf2 = new VF2(moleculeQ, moleculeT, true, true, false).isSubgraph();

        System.out.println("isSubset: " + vf2);

        Assert.assertFalse(vf2);
        Assert.assertFalse(uit);
    }
}
