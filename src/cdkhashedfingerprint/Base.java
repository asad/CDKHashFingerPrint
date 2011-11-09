/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cdkhashedfingerprint;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sf.jniinchi.INCHI_RET;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.Molecule;

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
                    ac.setID((fileIndex.getName().split(".mol"))[0]);
                    try {
                        String inchiKey = generateInchiKey(ac);
                        if (inchiKey == null || ac.getAtomCount() < 3) {
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
