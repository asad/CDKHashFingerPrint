/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cdkhashedfingerprint;

import helper.RandomNumber;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sf.jniinchi.INCHI_RET;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.signature.MoleculeFromSignatureBuilder;
import org.openscience.cdk.signature.MoleculeSignature;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.Molecule;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import signature.AbstractVertexSignature;
import signature.ColoredTree;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Base {

    public Base() {
    }

    /**
     *
     * @param ac
     * @return
     * @throws CDKException
     */
    public static String generateInchiKey(IAtomContainer ac) throws CDKException {
        InChIGenerator gen = InChIGeneratorFactory.getInstance().getInChIGenerator(ac);
        if (gen.getReturnStatus() != INCHI_RET.OKAY) {
            //System.err.println("inchi failed: " + gen.getMessage());
            return null;
        }
        return gen.getInchiKey();
    }

    /**
     *
     * @param ac
     * @return
     */
    public static Map<String, IAtomContainer> getSignatureFragments(IAtomContainer ac) {
        Map<String, IAtomContainer> fragments = new HashMap<String, IAtomContainer>();
        for (int i = 1; i < ac.getAtomCount(); i++) {
            int height = (int) RandomNumber.generateMersenneTwisterRandomNumber(ac.getAtomCount(), (long) i);
            MoleculeSignature sig = new MoleculeSignature(new Molecule(ac), height);
//            System.out.println("mol size " + ac.getAtomCount() + " height " + height);
            IAtomContainer fragment = makeMoleculeFromSignature(sig.toCanonicalString());
            String inchiKey;
            try {
                inchiKey = generateInchiKey(fragment);
                if (inchiKey == null || fragments.containsKey(inchiKey)) {
                    continue;
                }
                fragments.put(inchiKey, fragment);
            } catch (CDKException ex) {
                Logger.getLogger(Base.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return fragments;
    }

    /**
     * Convert a signature string into a molecule.
     *
     * @param signatureString
     * @return
     */
    public static IAtomContainer makeMoleculeFromSignature(String signatureString) {
        MoleculeFromSignatureBuilder builder = new MoleculeFromSignatureBuilder(new AtomContainer().getBuilder());
        ColoredTree tree = AbstractVertexSignature.parse(signatureString);
        builder.makeFromColoredTree(tree);
        IAtomContainer atomContainer = builder.getAtomContainer();
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
            CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
        } catch (CDKException ex) {
            Logger.getLogger(Base.class.getName()).log(Level.SEVERE, null, ex);
        }
        return atomContainer;
    }

    /**
     *
     * @param dir
     * @param cutoff 
     * @return
     * @throws FileNotFoundException
     * @throws CDKException
     * @throws IOException  
     */
    public static Map<String, IAtomContainer> readMDLMolecules(File dir, int cutoff) throws FileNotFoundException, CDKException, IOException {
        Map<String, IAtomContainer> inchiMolMap = new HashMap<String, IAtomContainer>();

        if (dir.isDirectory()) {
            File[] listFiles = dir.listFiles();
            for (File fileIndex : listFiles) {
                if (fileIndex.isFile() && fileIndex.getName().contains(".mol")) {
                    MDLV2000Reader reader = new MDLV2000Reader(new FileReader(fileIndex));
                    IAtomContainer ac = (IAtomContainer) reader.read(new Molecule());
                    String inchiKey;
                    try {
                        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
                        AtomContainerManipulator.removeHydrogens(ac);
                        CDKHueckelAromaticityDetector.detectAromaticity(ac);
                        inchiKey = generateInchiKey(ac);
                        if (inchiKey == null || inchiMolMap.containsKey(inchiKey) || ac.getAtomCount() < 3) {
                            continue;
                        }
                        inchiMolMap.put(inchiKey, ac);
                    } catch (Exception ex) {
                        Logger.getLogger(Base.class.getName()).log(Level.SEVERE, null, ex);
                    }

                    if (cutoff == inchiMolMap.size()) {
                        break;
                    }
                    reader.close();
                }
            }
        }
        return inchiMolMap;

    }
}
