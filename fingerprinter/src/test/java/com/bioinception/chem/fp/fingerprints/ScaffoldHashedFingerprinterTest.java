/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.bioinception.chem.fp.fingerprints;

import com.bioinception.chem.fp.fingerprints.bi.ScaffoldHashedFingerprinter;
import com.bioinception.chem.fp.fingerprints.helper.FingerprinterTool;
import static com.bioinception.chem.fp.fingerprints.helper.FingerprinterTool.isSubset;
import java.io.FileNotFoundException;
import java.util.BitSet;
import org.freehep.graphicsbase.util.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class ScaffoldHashedFingerprinterTest {

    final static SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

    /**
     * Test of HashedBloomFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGenerateFingerprint() throws InvalidSmilesException, CDKException {

        String smiles = "CCCCC1C(=O)N(N(C1=O)C1=CC=CC=C1)C1=CC=CC=C1";
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getBitFingerprint(molecule).asBitSet();
        System.out.println("fp " + fingerprint1.toString());
    }

    /**
     * Test of HashedBloomFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGenerateFingerprintIsSubset() throws InvalidSmilesException, CDKException {

        String smilesT = "NC(=O)C1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        String smilesQ = "CC1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
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
     * Test of HashedBloomFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     * @throws FileNotFoundException
     */
    @Test
    public void testGenerateFingerprintIsNotASubset1() throws InvalidSmilesException, CDKException, FileNotFoundException, FileNotFoundException {

        String smilesT = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O";
        String smilesQ = "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(O)=O";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesQ);
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesT);
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getBitFingerprint(moleculeQ).asBitSet();
        fingerprintT = fingerprint.getBitFingerprint(moleculeT).asBitSet();

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + isSubset(fingerprintT, fingerprintQ));

        Assert.assertFalse(isSubset(fingerprintT, fingerprintQ));
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
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
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
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
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
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        String smilest = "OC(=O)C1=CC=CC(O)=C1O";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);

        moleculeQ.setID("C00107");
        moleculeT.setID("C00196");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
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
    public void testGenerateFingerprintIsNotASubset5() throws InvalidSmilesException, Exception {

        String smilesq = "OC[C@H]1OC(=O)C[C@@H]1O";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        String smilest = "CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);
        moleculeQ.setID("CHEBI:17281");
        moleculeT.setID("CHEBI:35020");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
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
    public void testGenerateFingerprintIsNotASubset6() throws InvalidSmilesException, Exception {

        String smilesq = "CC1=C2CC[C@]3(C)CC[C@H](O)C(=C)[C@H]3C[C@H](CC1)C2(C)C";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        String smilest = "CC(C)=CCC[C@]1(C)[C@H]2CC=C(C)[C@@H]1C2";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);
        moleculeQ.setID("CHEBI:30038");
        moleculeT.setID("CHEBI:62756");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        IFingerprinter fingerprint = new ScaffoldHashedFingerprinter(1024);
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
}
