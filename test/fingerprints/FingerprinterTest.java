/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fingerprints;

import fingerprints.Fingerprinter;
import java.util.BitSet;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.FingerprinterTool;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.fingerprint.IFingerprinter;

/**
 *
 * @author Asad
 */
public class FingerprinterTest {

    public static void main(String[] args) throws InvalidSmilesException, CDKException {
//        testGenerateFingerprint();
        testGenerateFingerprintIsSubset();
    }

    public FingerprinterTest() {
    }

    /**
     * Test of Fingerprinter method
     * @throws InvalidSmilesException
     * @throws CDKException  
     */
    @Test
    public static void testGenerateFingerprint() throws InvalidSmilesException, CDKException {

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
        IFingerprinter fingerprint = new Fingerprinter(1024);
        BitSet fingerprintQ;
        BitSet fingerprintT;
        fingerprintQ = fingerprint.getFingerprint(moleculeQ);
        fingerprintT = fingerprint.getFingerprint(moleculeT);

        System.out.println("fpQ " + fingerprintQ.toString());
        System.out.println("fpT " + fingerprintT.toString());
        System.out.println("isSubset: " + FingerprinterTool.isSubset(fingerprintT, fingerprintQ));
    }
}
