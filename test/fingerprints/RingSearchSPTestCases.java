/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fingerprints;

import java.util.BitSet;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class RingSearchSPTestCases {

    @Test
    public void testGenerateFingerprintNaphthalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=CC=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        HashedSPFingerprinter fingerprint = new HashedSPFingerprinter(1024);
        fingerprint.setRespectRingMatches(true);
        fingerprint.setRespectFormalCharges(true);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getFingerprint(molecule);
        System.out.println("Naphthalene fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintAnthracene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC3=CC=CC=C3C=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        System.out.println("Atom count " + molecule.getAtomCount());
        HashedSPFingerprinter fingerprint = new HashedSPFingerprinter(1024);
        fingerprint.setRespectRingMatches(true);
        fingerprint.setRespectFormalCharges(true);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getFingerprint(molecule);
        System.out.println("Anthracene fp " + fingerprint1.toString());
    }

    @Test
    public void testGenerateFingerprintMultiphtalene() throws InvalidSmilesException, Exception {

        String smiles = "C1=CC2=CC=C3C4=CC5=CC6=CC=CC=C6C=C5C=C4C=CC3=C2C=C1";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        System.out.println("Atom count " + molecule.getAtomCount());
        HashedSPFingerprinter fingerprint = new HashedSPFingerprinter(1024);
        fingerprint.setRespectRingMatches(true);
        fingerprint.setRespectFormalCharges(true);
        BitSet fingerprint1;
        fingerprint1 = fingerprint.getFingerprint(molecule);
        System.out.println(" Multiphtalene fp " + fingerprint1.toString());
    }
}
