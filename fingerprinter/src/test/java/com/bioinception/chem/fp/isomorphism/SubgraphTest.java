package com.bioinception.chem.fp.isomorphism;

import com.bioinception.chem.fp.fingerprints.hashed.HashedFingerprinter;
import com.bioinception.chem.fp.fingerprints.interfaces.IFingerprinter;
import com.google.common.collect.FluentIterable;
import java.io.FileNotFoundException;
import java.util.BitSet;
import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 *
 * @author Asad
 */
public class SubgraphTest {

    final static SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());

    /**
     * Test of Fingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testIsSubgraph1() throws InvalidSmilesException, CDKException {

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
    public void testIsSubgraph2() throws InvalidSmilesException, CDKException, FileNotFoundException {
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

        BondMatcher bondmatcher = BondMatcher.forOrder();
        AtomMatcher atommatcher = AtomMatcher.forAny();

        int count_bond_match = FluentIterable.from(
                VentoFoggia.findSubstructure(moleculeT,
                        atommatcher, bondmatcher)
                        .matchAll(moleculeQ)).size();
        int count_match = FluentIterable.from(
                VentoFoggia.findSubstructure(moleculeT)
                        .matchAll(moleculeQ)).size();
        boolean trueMatch = count_bond_match > 0;
        boolean match = count_match > 0;

        System.out.println("isSubset match bonds too: " + trueMatch);
        System.out.println("isSubset: " + match);

        Assert.assertEquals(false, trueMatch);
        Assert.assertEquals(false, match);
    }

    /**
     * Test of HashedFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGenerateFingerprintIsSubset() throws InvalidSmilesException, CDKException {

        String smilesT
                = "NC(=O)C1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesT);
        String smilesQ = "CC1=C2C=CC(Br)=CC2=C(Cl)C=C1";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilesQ);

        moleculeQ.setID("Q");
        moleculeT.setID("T");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        moleculeQ = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        moleculeT = AtomContainerManipulator.removeHydrogens(moleculeT);

        BondMatcher bondmatcher = BondMatcher.forOrder();
        AtomMatcher atommatcher = AtomMatcher.forElement();

        int count_bond_match = FluentIterable.from(
                VentoFoggia.findSubstructure(moleculeT,
                        atommatcher, bondmatcher)
                        .matchAll(moleculeQ)).size();
        int count_match = FluentIterable.from(
                VentoFoggia.findSubstructure(moleculeT)
                        .matchAll(moleculeQ)).size();
        boolean trueMatch = count_bond_match > 0;
        boolean match = count_match > 0;

        System.out.println("isSubset match bonds too: " + trueMatch);
        System.out.println("isSubset: " + match);

        Assert.assertEquals(true, trueMatch);
        Assert.assertEquals(true, match);
    }

    /**
     * Test of HashedFingerprinter method
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGenerateFingerprintIsSubset6() throws InvalidSmilesException, CDKException {

        String smilesq = "CC1=C2CC[C@]3(C)CC[C@H](O)C(=C)[C@H]3C[C@H](CC1)C2(C)C";
        IAtomContainer moleculeQ = smilesParser.parseSmiles(smilesq);
        String smilest = "CC(C)=CCC[C@]1(C)[C@H]2CC=C(C)[C@@H]1C2";
        IAtomContainer moleculeT = smilesParser.parseSmiles(smilest);
        moleculeQ.setID("CHEBI:30038");
        moleculeT.setID("CHEBI:62756");
        System.out.println("Atom count Q:" + moleculeQ.getAtomCount());
        System.out.println("Atom count T:" + moleculeT.getAtomCount());
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeQ);
        moleculeQ = AtomContainerManipulator.removeHydrogens(moleculeQ);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(moleculeT);
        moleculeT = AtomContainerManipulator.removeHydrogens(moleculeT);

        BondMatcher bondmatcher = BondMatcher.forOrder();
        AtomMatcher atommatcher = AtomMatcher.forElement();

        int count_bond_match = FluentIterable.from(
                VentoFoggia.findSubstructure(moleculeT,
                        atommatcher, bondmatcher)
                        .matchAll(moleculeQ)).size();
        int count_match = FluentIterable.from(
                VentoFoggia.findSubstructure(moleculeT)
                        .matchAll(moleculeQ)).size();
        boolean trueMatch = count_bond_match > 0;
        boolean match = count_match > 0;

        System.out.println("isSubset match bonds too: " + trueMatch);
        System.out.println("isSubset: " + match);

        Assert.assertEquals(false, trueMatch);
        Assert.assertEquals(false, match);
    }
}
