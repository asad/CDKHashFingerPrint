/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fingerprints;

import fingerprints.interfaces.IFingerprinter;
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
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class HashedBloomFingerprinterTest {

    public static void main(String[] args) throws InvalidSmilesException, Exception {
        testGenerateFingerprint();
        testGenerateFingerprintIsSubset();
//        testGenerateFingerprintIsNotASubset1();
//        testGenerateFingerprintIsNotASubset2();
//        testGenerateFingerprintIsNotASubset3();
        testGenerateFingerprintIsNotASubset4();
    }

    public HashedBloomFingerprinterTest() {
    }

    /**
     * Test of HashedBloomFingerprinter method
     * @throws InvalidSmilesException
     * @throws CDKException  
     */
    @Test
    public static void testGenerateFingerprint() throws InvalidSmilesException, CDKException {

        String smiles = "CCCCC1C(=O)N(N(C1=O)C1=CC=CC=C1)C1=CC=CC=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new HashedBloomFingerprinter(1024);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getFingerprint(molecule);
        System.out.println("fp " + fingerprint1.toString());
    }

    /**
     * Test of HashedBloomFingerprinter method
     * @throws InvalidSmilesException
     * @throws CDKException  
     */
    @Test
    public static void testGenerateFingerprintIsSubset() throws InvalidSmilesException, CDKException {

        String smilesT =
                "NC(=O)C1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        String smilesQ = "CC1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new HashedBloomFingerprinter(1024);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getFingerprint(moleculeQ);
        fingerprintT = fingerprint.getFingerprint(moleculeT);

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintT, fingerprintQ));

        Assert.assertTrue(FingerprinterTool.isSubset(fingerprintT, fingerprintQ));
    }

    /**
     * Test of HashedBloomFingerprinter method
     * @throws InvalidSmilesException
     * @throws CDKException
     * @throws FileNotFoundException  
     */
    @Test
    public static void testGenerateFingerprintIsNotASubset1() throws InvalidSmilesException, CDKException, FileNotFoundException, FileNotFoundException {

        String smilesT =
                "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O";
        String smilesQ = "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(O)=O";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        smilesParser.setPreservingAromaticity(true);
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);

        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new HashedBloomFingerprinter(1024);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getFingerprint(moleculeQ);
        fingerprintT = fingerprint.getFingerprint(moleculeT);

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintT, fingerprintQ));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintT, fingerprintQ));
    }

    @Test
    public static void testGenerateFingerprintIsNotASubset2() throws InvalidSmilesException, Exception {

        FileReader smilesQ =
                new FileReader(System.getProperty("user.home") + File.separator + "Software/GITROOT/Fingerprint/test/data/mol/C00137.mol");
        FileReader smilesT = new FileReader(System.getProperty("user.home") + File.separator + "Software/GITROOT/Fingerprint/test/data/mol/C00257.mol");
        MDLV2000Reader readerQ = new MDLV2000Reader(smilesQ);
        MDLV2000Reader readerT = new MDLV2000Reader(smilesT);
        IAtomContainer moleculeQ = (IAtomContainer) readerQ.read(new Molecule());
        IAtomContainer moleculeT = (IAtomContainer) readerT.read(new Molecule());
        moleculeQ.setID((smilesQ.toString().split(".mol"))[0]);
        moleculeT.setID((smilesT.toString().split(".mol"))[0]);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new HashedBloomFingerprinter(1024);
        fingerprint.setRespectRingMatches(true);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getFingerprint(moleculeQ);
        fingerprintT = fingerprint.getFingerprint(moleculeT);

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintT, fingerprintQ));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintT, fingerprintQ));
    }

    @Test
    public static void testGenerateFingerprintIsNotASubset3() throws InvalidSmilesException, Exception {

        FileReader smilesQ =
                new FileReader(System.getProperty("user.home") + File.separator + "Software/GITROOT/Fingerprint/test/data/mol/C00186.mol");
        FileReader smilesT = new FileReader(System.getProperty("user.home") + File.separator + "Software/GITROOT/Fingerprint/test/data/mol/C00021.mol");
        MDLV2000Reader readerQ = new MDLV2000Reader(smilesQ);
        MDLV2000Reader readerT = new MDLV2000Reader(smilesT);
        IAtomContainer moleculeQ = (IAtomContainer) readerQ.read(new Molecule());
        IAtomContainer moleculeT = (IAtomContainer) readerT.read(new Molecule());
        moleculeQ.setID("C00186");
        moleculeT.setID("C00021");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new HashedBloomFingerprinter(1024);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        IAtomContainer removeHydrogens = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        IAtomContainer removeHydrogens1 = AtomContainerManipulator.removeHydrogens(moleculeT);

        fingerprintQ = fingerprint.getFingerprint(removeHydrogens);
        fingerprintT = fingerprint.getFingerprint(removeHydrogens1);

        System.out.println(moleculeQ.getID() + " fpQ " + fingerprintQ.toString());
        System.out.println(moleculeT.getID() + " fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintQ, fingerprintT));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintQ, fingerprintT));
    }

    @Test
    public static void testGenerateFingerprintIsNotASubset4() throws InvalidSmilesException, Exception {

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
        IFingerprinter fingerprint = new HashedBloomFingerprinter(1024);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        IAtomContainer removeHydrogens = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        IAtomContainer removeHydrogens1 = AtomContainerManipulator.removeHydrogens(moleculeT);

        fingerprintQ = fingerprint.getFingerprint(removeHydrogens);
        fingerprintT = fingerprint.getFingerprint(removeHydrogens1);

        System.out.println(moleculeQ.getID() + " fpQ " + fingerprintQ.toString());
        System.out.println(moleculeT.getID() + " fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintQ, fingerprintT));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintQ, fingerprintT));
    }
}
