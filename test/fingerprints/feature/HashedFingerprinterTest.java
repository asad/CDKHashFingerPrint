/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fingerprints.feature;

import java.io.FileNotFoundException;
import java.util.BitSet;
import junit.framework.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class HashedFingerprinterTest {

    /**
     * Test of HashedFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGenerateFingerprint() throws InvalidSmilesException, CDKException {

        String smiles = "CCCCC1C(=O)N(N(C1=O)C1=CC=CC=C1)C1=CC=CC=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new FeatureFingerprinter();
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    /**
     * Test of HashedFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGenerateFingerprintIsASubsetCase1() throws InvalidSmilesException, CDKException {

        String smilesT = "NC(=O)C1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        String smilesQ = "CC1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());

        IFingerprinter fingerprint = new FeatureFingerprinter();
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getBitFingerprint(moleculeQ).asBitSet();
        fingerprintT = fingerprint.getBitFingerprint(moleculeT).asBitSet();

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintT, fingerprintQ));

        Assert.assertTrue(FingerprinterTool.isSubset(fingerprintT, fingerprintQ));
    }

    /**
     * Test of HashedFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     * @throws FileNotFoundException
     */
    @Test
    public void testGenerateFingerprintIsNotASubsetCase2() throws InvalidSmilesException, CDKException, FileNotFoundException, FileNotFoundException {

        String smilesT
                = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O";
        String smilesQ = "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(O)=O";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new FeatureFingerprinter();
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getBitFingerprint(moleculeQ).asBitSet();
        fingerprintT = fingerprint.getBitFingerprint(moleculeT).asBitSet();

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintT, fingerprintQ));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintT, fingerprintQ));
    }

    @Test
    public void testGenerateFingerprintAnthracene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC3=CC=CC=C3C=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new FeatureFingerprinter();
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintNaphthalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=CC=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new FeatureFingerprinter();
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintMultiphtalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=C3C4=CC5=CC6=CC=CC=C6C=C5C=C4C=CC3=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new FeatureFingerprinter();
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }
}
