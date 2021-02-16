/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.bioinception.chem.fp.fingerprints.circle;

import com.bioinception.chem.fp.fingerprints.helper.FingerprinterTool;
import java.io.FileNotFoundException;
import java.util.BitSet;
import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class CircularFingerprinterTest {

    final static SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

    /**
     * Test of HashedFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGenerateFingerprint() throws InvalidSmilesException, CDKException {

        String smiles = "CCCCC1C(=O)N(N(C1=O)C1=CC=CC=C1)C1=CC=CC=C1";
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
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
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());

        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getBitFingerprint(moleculeQ).asBitSet();
        fingerprintT = fingerprint.getBitFingerprint(moleculeT).asBitSet();

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintT, fingerprintQ));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintT, fingerprintQ));
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
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
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
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintNaphthalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=CC=C2C=C1";
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintMultiphtalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=C3C4=CC5=CC6=CC=CC=C6C=C5C=C4C=CC3=C2C=C1";
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        System.out.println("Atom count " + molecule.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    /**
     * Test of HashedBloomFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     * @throws FileNotFoundException
     */
    @Test
    public void testGenerateFingerprintIsNotASubset1() throws InvalidSmilesException, CDKException, FileNotFoundException, FileNotFoundException {

        String smilesT
                = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O";
        String smilesQ = "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(O)=O";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);

        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
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
    public void testGenerateFingerprintIsNotASubset2() throws InvalidSmilesException, Exception {

        String smilesq = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        String smilest = "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(O)=O";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);
        moleculeQ.setID("C00137");
        moleculeT.setID("C00257");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
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
    public void testGenerateFingerprintIsNotASubset3() throws InvalidSmilesException, Exception {

        String smilesq = "C[C@H](O)C(O)=O";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        String smilest = "N[C@@H](CCSC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=NC2=C1N=CN=C2N)C(O)=O";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);
        moleculeQ.setID("C00186");
        moleculeT.setID("C00021");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        IAtomContainer removeHydrogens = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        IAtomContainer removeHydrogens1 = AtomContainerManipulator.removeHydrogens(moleculeT);

        fingerprintQ = fingerprint.getBitFingerprint(removeHydrogens).asBitSet();
        fingerprintT = fingerprint.getBitFingerprint(removeHydrogens1).asBitSet();

        System.out.println(moleculeQ.getID() + " fpQ " + fingerprintQ.toString());
        System.out.println(moleculeT.getID() + " fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintQ, fingerprintT));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintQ, fingerprintT));
    }

    @Test
    public void testGenerateFingerprintIsNotASubset4() throws InvalidSmilesException, Exception {

        String smilesq = "NC([*])C(=O)NC([*])C(O)=O";
        String smilest = "OC(=O)C1=CC=CC(O)=C1O";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);

        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);

        moleculeQ.setID("C00107");
        moleculeT.setID("C00196");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        IAtomContainer removeHydrogens = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        IAtomContainer removeHydrogens1 = AtomContainerManipulator.removeHydrogens(moleculeT);

        fingerprintQ = fingerprint.getBitFingerprint(removeHydrogens).asBitSet();
        fingerprintT = fingerprint.getBitFingerprint(removeHydrogens1).asBitSet();

        System.out.println(moleculeQ.getID() + " fpQ " + fingerprintQ.toString());
        System.out.println(moleculeT.getID() + " fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintQ, fingerprintT));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintQ, fingerprintT));
    }

    @Test
    public void testSim() throws InvalidSmilesException, Exception {

        String smilesq = "N[C@@H](CC1=CC=C(O)C=C1)C(O)=O";
////        String s1 = "N[C@@H](CC1=CC=C(O)C=C1)C(*)=O";
        String smilest = "OC(=O)C(=O)CC1=CC=C(O)C=C1";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);

        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        CircularFingerprinter fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        IAtomContainer removeHydrogens = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        IAtomContainer removeHydrogens1 = AtomContainerManipulator.removeHydrogens(moleculeT);

        fingerprintQ = fingerprint.getBitFingerprint(removeHydrogens).asBitSet();
        fingerprintT = fingerprint.getBitFingerprint(removeHydrogens1).asBitSet();

        System.out.println(moleculeQ.getID() + " fpQ " + fingerprintQ.toString());
        System.out.println(moleculeT.getID() + " fpT " + fingerprintT.toString());
        System.out.println("Sim " +Tanimoto.calculate(fingerprintQ, fingerprintT));
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintQ, fingerprintT));

        Assert.assertFalse(FingerprinterTool.isSubset(fingerprintQ, fingerprintT));
    }
}
